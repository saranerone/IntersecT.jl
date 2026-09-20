"""IntersecT"""
module IntersecT

using DataFrames
using CSV
using Statistics

export IntersecTResult, IntersecTLog, run_intersect, generate_log, calc_obs_err, load_model_csv, load_measurements_csv, list_measured_phases, select_phases

# ============================================================
# Constants
# ============================================================

const MIN_COEFF_ERR       = 1
const MAX_COEFF_ERR       = 6
const MIN_ERROR_THRESHOLD = 0.01

"""Empirical parameters for uncertainty estimation by analytical technique."""
const ANALYSIS_PARAMS = Dict(
    "EDS"      => (0.0703, 0.3574, 0.01,  0.1),
    "WDS map"  => (0.0434, 0.3451, 0.005, 0.05),
    "WDS spot" => (0.023,  0.2772, 0.005, 0.05),
)

# ============================================================
# Result struct
# ============================================================

"""IntersecTResult"""
struct IntersecTResult
    x                     :: Vector{Float64}
    y                     :: Vector{Float64}
    x_label               :: String
    y_label               :: String
    phase_names           :: Vector{String}
    element_names         :: Vector{String}
    Qcmp_elem             :: Matrix{Float64}
    Qcmp_phase            :: Matrix{Float64}
    redchi2_phase         :: Matrix{Float64}
    redchi2_tot           :: Vector{Float64}
    Qcmp_unweighted       :: Vector{Float64}
    Qcmp_weighted         :: Vector{Float64}
    min_redchi2           :: Vector{Float64}   # length n_phases
    n_elements_per_phase  :: Vector{Int}       # length n_phases
    phase_bases           :: Vector{String}    # phase names without domain label
    analysis_type         :: String
end

# ============================================================
# Unit conversion
# ============================================================

"""convert_coordinate(values, colname) -> (converted_values, new_label)"""
function convert_coordinate(values::Vector{Float64}, colname::String)
    n = lowercase(colname)
    if occursin("t(k)", n) || occursin("t[k]", n) || occursin("t_k", n)
        return values .- 273.0, "T(C)"
    elseif occursin("t(c)", n) || occursin("t[c]", n) || occursin("t_c", n)
        return values, "T(C)"
    elseif occursin("kbar", n)
        return values ./ 10.0, "P(GPa)"
    elseif occursin("p(bar)", n) || occursin("p[bar]", n) || occursin("p_bar", n)
        return values ./ 10_000.0, "P(GPa)"
    else
        return values, colname
    end
end

# ============================================================
# Uncertainty calculation
# ============================================================

"""calc_obs_err(apfu_obs, analysis_type) -> Vector{Float64}"""
function calc_obs_err(apfu_obs::Vector{Float64}, analysis_type::String)::Vector{Float64}
    params = get(ANALYSIS_PARAMS, analysis_type, nothing)
    if isnothing(params)
        valid = join(keys(ANALYSIS_PARAMS), ", ")
        error("Unknown analysis type: " * analysis_type * ". Valid options: " * valid)
    end
    coeff, exp_, min_err, max_err = params
    return clamp.(coeff .* apfu_obs .^ exp_, min_err, max_err)
end

# ============================================================
# Core statistical functions
# Each operates on vectors of observed/modelled values for one grid point.
# Parallelisation happens at the grid-point loop level in run_intersect.
# ============================================================

"""q_elem(apfu_obs, model_val, obs_err) -> Float64"""
@inline function q_elem(apfu_obs::Float64, model_val::Float64, obs_err::Float64)::Float64
    sigma = max(obs_err, MIN_ERROR_THRESHOLD)
    diff  = abs(apfu_obs - model_val)
    num   = clamp(diff - sigma / MIN_COEFF_ERR, 0.0, MAX_COEFF_ERR * sigma)
    return 100.0 * abs(1.0 - num / (MAX_COEFF_ERR * sigma))^(model_val + 1.0)
end

"""q_phase(apfu_obs, model_vals, obs_errs) -> Float64"""
@inline function q_phase(
    apfu_obs   :: AbstractVector{Float64},
    model_vals :: AbstractVector{Float64},
    obs_errs   :: AbstractVector{Float64}
)::Float64
    n     = length(apfu_obs)
    total = 0.0
    @inbounds for k in 1:n
        sigma = max(obs_errs[k], MIN_ERROR_THRESHOLD)
        diff  = abs(apfu_obs[k] - model_vals[k])
        num   = clamp(diff - sigma / MIN_COEFF_ERR, 0.0, MAX_COEFF_ERR * sigma)
        total += abs(1.0 - num / (MAX_COEFF_ERR * sigma))^(model_vals[k] + 1.0)
    end
    return total / n * 100.0
end

"""chi2_stat(apfu_obs, model_vals, obs_errs) -> Float64"""
@inline function chi2_stat(
    apfu_obs   :: AbstractVector{Float64},
    model_vals :: AbstractVector{Float64},
    obs_errs   :: AbstractVector{Float64}
)::Float64
    total = 0.0
    @inbounds for k in eachindex(apfu_obs)
        total += (apfu_obs[k] - model_vals[k])^2 / obs_errs[k]^2
    end
    return total
end

"""red_chi2_stat(apfu_obs, model_vals, obs_errs, dof) -> Float64"""
@inline function red_chi2_stat(
    apfu_obs   :: AbstractVector{Float64},
    model_vals :: AbstractVector{Float64},
    obs_errs   :: AbstractVector{Float64},
    dof        :: Int
)::Float64
    return chi2_stat(apfu_obs, model_vals, obs_errs) / (dof - 1)
end

# ============================================================
# Input parsing
# ============================================================

_isblank(x) = ismissing(x) || (x isa AbstractString && isempty(strip(x)))
_isauto(x)  = x isa AbstractString && lowercase(strip(x)) == "auto"
_isnumeric(x) = !_isblank(x) && !isnothing(tryparse(Float64, strip(string(x))))

# Removes the suffix CSV.jl adds to duplicated column names (Grt_Mg_1 -> Grt_Mg).
function _strip_dedup(names_raw::Vector{String})
    out = similar(names_raw)
    for j in eachindex(names_raw)
        m = match(r"^(.*_.*)_(\d+)$", names_raw[j])
        out[j] = isnothing(m) ? names_raw[j] : String(m.captures[1])
    end
    return out
end

# The row below the header is a domain row if it holds at least one label.
function _has_domain_row(row)
    all(_isblank, row) &&
        error("The row below the header is empty. Remove it, or fill in the domain labels.")
    return any(x -> !_isblank(x) && !_isnumeric(x), row)
end

_qualify(base::String, domain::String) = isempty(domain) ? base : base * " " * domain

function _tofloat(x, col::AbstractString)
    ismissing(x) && error("Empty cell in column '$col': expected a number.")
    x isa AbstractString || return Float64(x)
    v = tryparse(Float64, strip(x))
    isnothing(v) && error("Cannot read '$x' as a number in column '$col'.")
    return v
end

"""parse_measurements(df) -> (element_names, element_labels, apfu_obs, obs_err,
                              phase_names, phase_bases, phase_ids)"""
function parse_measurements(df::DataFrame)
    raw_names     = String.(names(df))
    element_names = _strip_dedup(raw_names)          # used to query the model
    row1          = Vector(df[1, :])

    has_domain = _has_domain_row(row1)
    nrow(df) < (has_domain ? 3 : 2) &&
        error("The measurements file has too few rows.")

    if has_domain
        domains = [_isblank(x) ? "" : String(strip(string(x))) for x in row1]
        any(d -> lowercase(d) == "auto", domains) &&
            error("'auto' is not a valid domain label.")
        row_obs = Vector(df[2, :])
        row_err = Vector(df[3, :])
    else
        domains = fill("", length(element_names))
        row_obs = row1
        row_err = Vector(df[2, :])
    end

    apfu_obs = [_tofloat(row_obs[j], element_names[j]) for j in eachindex(row_obs)]

    all_nan = all(x -> ismissing(x) || (x isa Number && isnan(Float64(x))), row_err)

    obs_err = if _isauto(row_err[1]) || all_nan
        if _isauto(row_err[1]) && !all(_isblank, row_err[2:end])
            error("Uncertainty row: when 'auto' is used, all remaining cells must be empty.")
        end
        Float64[]
    else
        if any(_isblank, row_err) || any(_isauto, row_err)
            error("The uncertainty row must be either fully numeric, or 'auto' in the first cell only.")
        end
        [_tofloat(row_err[j], element_names[j]) for j in eachindex(row_err)]
    end

    element_labels = String[]
    phase_names    = String[]    # qualified, for display
    phase_bases    = String[]    # base names, for the model
    phase_ids      = Int[]

    for (j, name) in enumerate(element_names)
        parts   = split(name, "_"; limit=2)
        base    = length(parts) >= 2 ? String(parts[1]) : name
        element = length(parts) >= 2 ? String(parts[2]) : ""
        key     = _qualify(base, domains[j])

        push!(element_labels, isempty(element) ? key : key * "_" * element)

        idx = findfirst(==(key), phase_names)
        if isnothing(idx)
            push!(phase_names, key)
            push!(phase_bases, base)
            push!(phase_ids, length(phase_names))
        else
            push!(phase_ids, idx)
        end
    end

    return element_names, element_labels, apfu_obs, obs_err,
           phase_names, phase_bases, phase_ids
end

"""list_measured_phases(df) -> (labels, bases)"""
function list_measured_phases(df::DataFrame)
    _, _, _, _, phase_names, phase_bases, _ = parse_measurements(df)
    return phase_names, phase_bases
end

"""select_phases(df, selected) -> DataFrame"""
function select_phases(df::DataFrame, selected::Vector{String})
    _, _, _, _, phase_names, phase_bases, phase_ids = parse_measurements(df)
    idx_sel = [findfirst(==(s), phase_names) for s in selected]
    any(isnothing, idx_sel) &&
        error("Unknown phase selection. Available: " * join(phase_names, ", "))
    idx_sel = Int.(idx_sel)
    length(unique(phase_bases[idx_sel])) == length(idx_sel) ||
        error("Only one domain per phase can be selected in a single run.")
    return df[:, findall(p -> p in idx_sel, phase_ids)]
end

"""extract_model_data(model_df, element_names, x_col, y_col)"""
function extract_model_data(
    model_df      :: DataFrame,
    element_names :: Vector{String},
    x_col         :: String,
    y_col         :: String
)
    cols = names(model_df)

    function find_col(target)
        idx = findfirst(c -> lowercase(c) == lowercase(target), cols)
        if isnothing(idx)
            error("Column not found: " * target * ". Available: " * join(cols, ", "))
        end
        return cols[idx]
    end

    xc    = find_col(x_col)
    yc    = find_col(y_col)
    x_raw = Float64.(model_df[!, xc])
    y_raw = Float64.(model_df[!, yc])
    x, x_label = convert_coordinate(x_raw, xc)
    y, y_label = convert_coordinate(y_raw, yc)

    n_points      = nrow(model_df)
    n_elements    = length(element_names)
    model_matrix  = fill(NaN, n_points, n_elements)

    for (j, elem_name) in enumerate(element_names)
        cidx = findfirst(c -> lowercase(c) == lowercase(elem_name), cols)
        if !isnothing(cidx)
            col_vals = model_df[!, cols[cidx]]
            for i in 1:n_points
                v = col_vals[i]
                model_matrix[i, j] = (ismissing(v) || (v isa Number && isnan(Float64(v)))) ? NaN : Float64(v)
            end
        else
            @warn "Column '$elem_name' not found in the model output. " *
                  "The corresponding map will be empty. The model is queried " *
                  "with the base phase name, not with the domain label."
        end
    end

    return x, y, x_label, y_label, model_matrix
end

# ============================================================
# I/O helpers for standalone use
# ============================================================

"""load_model_csv(path; x_col="T(K)", y_col="P(bar)") -> DataFrame"""
function load_model_csv(path::String)::DataFrame
    return CSV.read(path, DataFrame)
end

"""load_measurements_csv(path) -> DataFrame"""
function load_measurements_csv(path::String)::DataFrame
    return CSV.read(path, DataFrame; header=true)
end

# ============================================================
# Main calculation pipeline
# ============================================================

"""run_intersect(model_df, measurements_df;"""
function run_intersect(
    model_df        :: DataFrame,
    measurements_df :: DataFrame;
    x_col           :: String = "T(K)",
    y_col           :: String = "P(bar)",
    analysis_type   :: String = "WDS spot"
)::IntersecTResult

    # --- parse measurements ---
    element_names, element_labels, apfu_obs, obs_err, phase_names, phase_bases, phase_ids =
        parse_measurements(measurements_df)

    n_elements = length(element_names)
    n_phases   = length(phase_names)

    if isempty(obs_err)
        obs_err = calc_obs_err(apfu_obs, analysis_type)
    end
    obs_err = max.(obs_err, MIN_ERROR_THRESHOLD)

    # --- extract model data ---
    x, y, x_label, y_label, model_matrix =
        extract_model_data(model_df, element_names, x_col, y_col)

    n_points = length(x)

    # --- precompute per-phase element indices ---
    phase_elem_idx = [findall(==(p), phase_ids) for p in 1:n_phases]

    # --- warn if any phase has too few elements for reliable chi-squared ---
    for p in 1:n_phases
        f = length(phase_elem_idx[p])
        if f == 1
            @warn "Phase $(phase_names[p]) has only 1 element. χ2 is unreliable."
        elseif f == 2
            @warn "Phase $(phase_names[p]) has only 2 elements. χ2 (not reduced) will be used."
        end
    end

    # --- allocate output arrays (NaN = absent, not yet computed) ---
    Qcmp_elem_mat     = fill(NaN, n_points, n_elements)
    Qcmp_phase_mat    = fill(NaN, n_points, n_phases)
    redchi2_phase_mat = fill(NaN, n_points, n_phases)
    redchi2_tot_vec   = fill(NaN, n_points)
    Qcmp_unweighted   = fill(NaN, n_points)
    Qcmp_weighted_vec = fill(NaN, n_points)

    # --- main loop: parallelised over grid points ---
    # Each grid point is independent -> trivially parallel with Threads.@threads.
    # Start Julia with `julia --threads auto` to use all CPU cores.
    Threads.@threads for i in 1:n_points

        model_row = @view model_matrix[i, :]

        # --- per-phase calculations ---
        # Each phase is evaluated independently: if a phase is absent at this
        # point (any of its elements is NaN), that phase's outputs stay NaN
        # but other phases are still computed normally.
        for p in 1:n_phases
            eidx = phase_elem_idx[p]

            # Check if this phase is present (no NaN among its elements)
            phase_present = true
            for j in eidx
                if isnan(model_row[j])
                    phase_present = false
                    break
                end
            end
            phase_present || continue   # leave NaN in outputs for this phase

            mv = model_row[eidx]
            ao = apfu_obs[eidx]
            oe = obs_err[eidx]

            # element-level quality factors for elements of this phase
            for j in eidx
                Qcmp_elem_mat[i, j] = q_elem(apfu_obs[j], model_row[j], obs_err[j])
            end

            # phase-level quality factor
            Qcmp_phase_mat[i, p] = q_phase(ao, mv, oe)

            # phase-level chi-squared (reduced if >2 elements, plain chi2 otherwise)
            f = length(eidx)
            redchi2_phase_mat[i, p] = f > 2 ?
                red_chi2_stat(ao, mv, oe, f) :
                chi2_stat(ao, mv, oe)
        end

        # --- totals: only where ALL phases are present ---
        # If any phase is absent at this point, unweighted, weighted, and
        # redchi2_tot remain NaN. This delimits the full-assemblage stability field.
        all_present = !any(isnan, @view Qcmp_phase_mat[i, :])
        all_present || continue

        redchi2_tot_vec[i] = red_chi2_stat(apfu_obs, model_row, obs_err, n_elements)

    end  # end Threads.@threads

    # --- minimum reduced chi-squared per phase (used for weighting) ---
    # Computed only over points where the phase is present (finite values).
    # For phases with <=2 elements, chi2 is used instead of reduced chi2,
    # and 1 is added so the minimum floor is consistent with the >2 case
    # (where reduced chi2 = 1 means perfect fit within uncertainties).
    min_redchi2 = Vector{Float64}(undef, n_phases)
    for p in 1:n_phases
        finite_vals = filter(isfinite, @view redchi2_phase_mat[:, p])
        min_val     = isempty(finite_vals) ? 1.0 : minimum(finite_vals)
        f           = length(phase_elem_idx[p])
        min_redchi2[p] = f <= 2 ? 1.0 + min_val : min_val
    end

    # --- unweighted total Qcmp ---
    # Equal weight for all phases. NaN where any phase is absent (propagates automatically).
    w_equal = fill(1.0 / n_phases, n_phases)
    for i in 1:n_points
        row = @view Qcmp_phase_mat[i, :]
        any(isnan, row) && continue
        Qcmp_unweighted[i] = dot_safe(row, w_equal)
    end

    # --- weighted total Qcmp ---
    # Weight = 1 / min_redchi2, floored at 1 before inversion.
    min_rc2_clamped = max.(min_redchi2, 1.0)
    weight_raw      = 1.0 ./ min_rc2_clamped
    w_weighted      = weight_raw ./ sum(weight_raw)
    for i in 1:n_points
        row = @view Qcmp_phase_mat[i, :]
        any(isnan, row) && continue
        Qcmp_weighted_vec[i] = dot_safe(row, w_weighted)
    end

    n_elements_per_phase = [length(phase_elem_idx[p]) for p in 1:n_phases]

    return IntersecTResult(
        x, y,
        x_label, y_label,
        phase_names,
        element_labels,
        Qcmp_elem_mat,
        Qcmp_phase_mat,
        redchi2_phase_mat,
        redchi2_tot_vec,
        Qcmp_unweighted,
        Qcmp_weighted_vec,
        min_redchi2,
        n_elements_per_phase,
        phase_bases,
        analysis_type,
    )
end

# Simple dot product without importing LinearAlgebra
@inline function dot_safe(a::AbstractVector{Float64}, b::AbstractVector{Float64})::Float64
    total = 0.0
    @inbounds for k in eachindex(a)
        total += a[k] * b[k]
    end
    return total
end

# ============================================================
# Log
# ============================================================

"""IntersecTLog"""
struct IntersecTLog
    phase_names                :: Vector{String}
    element_names              :: Vector{String}
    Qcmp_elem_max              :: Vector{Float64}    # max Qcmp per element (n_elements,)
    Qcmp_elem_max_x            :: Vector{Float64}
    Qcmp_elem_max_y            :: Vector{Float64}
    Qcmp_phase_max             :: Vector{Float64}    # max Qcmp per phase (n_phases,)
    Qcmp_phase_max_x           :: Vector{Float64}
    Qcmp_phase_max_y           :: Vector{Float64}
    redchi2_phase_min          :: Vector{Float64}    # min chi2 per phase (n_phases,)
    redchi2_phase_label        :: Vector{String}     # "chi2" or "reduced chi2"
    redchi2_warnings           :: Vector{String}
    redchi2_tot_min            :: Float64
    redchi2_tot_min_x          :: Float64
    redchi2_tot_min_y          :: Float64
    phase_weights              :: Vector{Float64}    # normalised weights (n_phases,)
    Qcmp_weighted_max          :: Float64
    Qcmp_weighted_max_x        :: Float64
    Qcmp_weighted_max_y        :: Float64
    Qcmp_phase_at_weighted_max :: Vector{Float64}    # per-phase Qcmp at weighted max
    Qcmp_unweighted_max        :: Float64
    Qcmp_unweighted_max_x      :: Float64
    Qcmp_unweighted_max_y      :: Float64
end

# Returns (mean_x, mean_y, extremal_value) ignoring NaN.
# If multiple points share the extremum, returns their mean position.
function _extremum_position(values, x, y; find=:max)
    finite_vals = filter(isfinite, values)
    isempty(finite_vals) && return NaN, NaN, NaN
    ext = find == :max ? maximum(finite_vals) : minimum(finite_vals)
    idx = findall(v -> isfinite(v) && v == ext, values)
    return sum(x[i] for i in idx) / length(idx), sum(y[i] for i in idx) / length(idx), ext
end

"""generate_log(result) -> IntersecTLog"""
function generate_log(result::IntersecTResult)::IntersecTLog
    n_phases   = length(result.phase_names)
    n_elements = length(result.element_names)
    x, y       = result.x, result.y

    # element maxima
    Qcmp_elem_max   = Vector{Float64}(undef, n_elements)
    Qcmp_elem_max_x = Vector{Float64}(undef, n_elements)
    Qcmp_elem_max_y = Vector{Float64}(undef, n_elements)
    for j in 1:n_elements
        Qcmp_elem_max_x[j], Qcmp_elem_max_y[j], Qcmp_elem_max[j] =
            _extremum_position(result.Qcmp_elem[:, j], x, y)
    end

    # phase maxima
    Qcmp_phase_max   = Vector{Float64}(undef, n_phases)
    Qcmp_phase_max_x = Vector{Float64}(undef, n_phases)
    Qcmp_phase_max_y = Vector{Float64}(undef, n_phases)
    for p in 1:n_phases
        Qcmp_phase_max_x[p], Qcmp_phase_max_y[p], Qcmp_phase_max[p] =
            _extremum_position(result.Qcmp_phase[:, p], x, y)
    end

    # per-phase chi-squared minima and warnings
    redchi2_phase_min   = Vector{Float64}(undef, n_phases)
    redchi2_phase_label = Vector{String}(undef, n_phases)
    redchi2_warnings    = String[]
    for p in 1:n_phases
        finite_vals = filter(isfinite, result.redchi2_phase[:, p])
        redchi2_phase_min[p] = isempty(finite_vals) ? NaN : minimum(finite_vals)
        f = result.n_elements_per_phase[p]
        if f == 1
            redchi2_phase_label[p] = "χ2"
            push!(redchi2_warnings,
                "WARNING: Phase $(result.phase_names[p]) has only 1 element. " *
                "χ2 is unreliable.")
        elseif f == 2
            redchi2_phase_label[p] = "χ2"
        else
            redchi2_phase_label[p] = "reduced χ2"
        end
    end

    # total redchi2 minimum
    rc2_min_x, rc2_min_y, rc2_min = _extremum_position(result.redchi2_tot, x, y; find=:min)

    # weights
    min_rc2_clamped = max.(result.min_redchi2, 1.0)
    weight_raw      = 1.0 ./ min_rc2_clamped
    phase_weights   = weight_raw ./ sum(weight_raw)

    # weighted max
    wx, wy, wmax = _extremum_position(result.Qcmp_weighted, x, y)
    best_idx     = findfirst(v -> isfinite(v) && v == wmax, result.Qcmp_weighted)
    Qcmp_at_best = isnothing(best_idx) ? fill(NaN, n_phases) :
                   Vector{Float64}(result.Qcmp_phase[best_idx, :])

    # unweighted max
    ux, uy, umax = _extremum_position(result.Qcmp_unweighted, x, y)

    return IntersecTLog(
        result.phase_names, result.element_names,
        Qcmp_elem_max, Qcmp_elem_max_x, Qcmp_elem_max_y,
        Qcmp_phase_max, Qcmp_phase_max_x, Qcmp_phase_max_y,
        redchi2_phase_min, redchi2_phase_label, redchi2_warnings,
        rc2_min, rc2_min_x, rc2_min_y,
        phase_weights,
        wmax, wx, wy, Qcmp_at_best,
        umax, ux, uy,
    )
end

"""format_log(log, x_label, y_label) -> String"""
function format_log(log::IntersecTLog, x_label::String, y_label::String)::String
    L = String[]
    sep = "=" ^ 60

    push!(L, sep, "IntersecT — Analysis Report", sep)
    push!(L, "\nPhases:   " * join(log.phase_names, ", "))
    push!(L, "Elements: " * join(log.element_names, ", "))

    push!(L, "\n--- Maximum Qcmp per element ---")
    for j in eachindex(log.element_names)
        push!(L, "  $(log.element_names[j]) : $(round(log.Qcmp_elem_max[j], digits=2))  " *
              "at $x_label = $(round(log.Qcmp_elem_max_x[j], digits=2)), " *
              "$y_label = $(round(log.Qcmp_elem_max_y[j], digits=4))")
    end

    push!(L, "\n--- Maximum Qcmp per phase ---")
    for p in eachindex(log.phase_names)
        push!(L, "  $(log.phase_names[p]) : $(round(log.Qcmp_phase_max[p], digits=2))  " *
              "at $x_label = $(round(log.Qcmp_phase_max_x[p], digits=2)), " *
              "$y_label = $(round(log.Qcmp_phase_max_y[p], digits=4))")
    end

    push!(L, "\n--- Minimum χ2 statistic per phase ---")
    for p in eachindex(log.phase_names)
        push!(L, "  $(log.phase_names[p]) : min $(log.redchi2_phase_label[p]) = " *
              "$(round(log.redchi2_phase_min[p], digits=4))  " *
              "(weight = $(round(log.phase_weights[p], digits=4)))")
    end
    isempty(log.redchi2_warnings) || append!(L, ["", log.redchi2_warnings...])

    push!(L, "\n--- Minimum total reduced χ2 ---")
    push!(L, "  $(round(log.redchi2_tot_min, digits=4))  " *
          "at $x_label = $(round(log.redchi2_tot_min_x, digits=2)), " *
          "$y_label = $(round(log.redchi2_tot_min_y, digits=4))")

    push!(L, "\n--- Phase weights ---")
    for p in eachindex(log.phase_names)
        push!(L, "  $(log.phase_names[p]) : $(round(log.phase_weights[p], digits=4))")
    end

    push!(L, "\n--- Maximum Q*cmp weighted ---")
    push!(L, "  $(round(log.Qcmp_weighted_max, digits=2))  " *
          "at $x_label = $(round(log.Qcmp_weighted_max_x, digits=2)), " *
          "$y_label = $(round(log.Qcmp_weighted_max_y, digits=4))")
    for p in eachindex(log.phase_names)
        push!(L, "  $(log.phase_names[p]) : $(round(log.Qcmp_phase_at_weighted_max[p], digits=2))")
    end

    push!(L, "\n--- Maximum Q*cmp unweighted ---")
    push!(L, "  $(round(log.Qcmp_unweighted_max, digits=2))  " *
          "at $x_label = $(round(log.Qcmp_unweighted_max_x, digits=2)), " *
          "$y_label = $(round(log.Qcmp_unweighted_max_y, digits=4))")

    push!(L, "\n" * sep)
    return join(L, "\n")
end

end  # module IntersecT