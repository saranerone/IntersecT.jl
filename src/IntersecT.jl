"""
    IntersecT

Quantitative Isopleth Thermobarometry from Thermodynamic Models.
Nerone et al. (2025), https://doi.org/10.1016/j.cageo.2025.105949

Calculates quality factors (Qcmp) and reduced chi-squared statistics
to assess the agreement between measured and modelled mineral compositions.

Based on the quality factor approach of Duesterhoeft & Lanari (2020),
doi:10.1111/jmg.12538.

# Main entry point
    run_intersect(model_df, measurements_df; analysis_type="WDS spot") -> IntersecTResult

# I/O helpers (for standalone use; replaced by MAGEMin interface in production)
    load_model_csv(path; x_col, y_col)      -> DataFrame
    load_measurements_csv(path)             -> DataFrame

# Output struct
    IntersecTResult with fields:
    - x, y              : coordinate arrays (n_points,)
    - x_label, y_label  : axis labels after unit conversion
    - phase_names       : names of the phases considered
    - element_names     : "Phase_Element" labels
    - Qcmp_elem         : quality factor per element       (n_points × n_elements)
    - Qcmp_phase        : quality factor per phase         (n_points × n_phases)
    - redchi2_phase     : reduced chi-squared per phase    (n_points × n_phases)
    - redchi2_tot       : total reduced chi-squared        (n_points,)
    - Qcmp_unweighted   : unweighted total quality factor  (n_points,)
    - Qcmp_weighted     : weighted total quality factor    (n_points,)
    - min_redchi2       : minimum reduced chi-squared per phase (n_phases,)

NaN policy
----------
Where a phase is absent at a grid point, all per-phase and total outputs
(Qcmp_phase, redchi2_phase, Qcmp_unweighted, Qcmp_weighted, redchi2_tot)
are set to NaN for that point. Only grid points where ALL requested phases
are present receive a numeric result. Per-element outputs (Qcmp_elem) follow
the same rule: NaN where the element's phase is absent.
"""
module IntersecT

using DataFrames
using Statistics

export IntersecTResult, run_intersect, calc_obs_err, load_model_csv, load_measurements_csv

# ============================================================
# Constants
# ============================================================

const MIN_COEFF_ERR       = 1
const MAX_COEFF_ERR       = 6
const MIN_ERROR_THRESHOLD = 0.01

"""
Empirical parameters for uncertainty estimation by analytical technique.
Tuple layout: (coefficient, exponent, min_error, max_error)
Formula: calc_err = clamp(coeff * apfu^exponent, min_error, max_error)
"""
const ANALYSIS_PARAMS = Dict(
    "EDS"      => (0.0703, 0.3574, 0.01,  0.1),
    "WDS map"  => (0.0434, 0.3451, 0.005, 0.05),
    "WDS spot" => (0.023,  0.2772, 0.005, 0.05),
)

# ============================================================
# Result struct
# ============================================================

"""
    IntersecTResult

Container for all outputs of a quality factor analysis.
All matrices have dimensions (n_points × n_phases) or (n_points × n_elements).
Vectors have length n_points unless noted otherwise.
See module docstring for the NaN policy at absent-phase grid points.
"""
struct IntersecTResult
    x                :: Vector{Float64}
    y                :: Vector{Float64}
    x_label          :: String
    y_label          :: String
    phase_names      :: Vector{String}
    element_names    :: Vector{String}
    Qcmp_elem        :: Matrix{Float64}
    Qcmp_phase       :: Matrix{Float64}
    redchi2_phase    :: Matrix{Float64}
    redchi2_tot      :: Vector{Float64}
    Qcmp_unweighted  :: Vector{Float64}
    Qcmp_weighted    :: Vector{Float64}
    min_redchi2      :: Vector{Float64}   # length n_phases
end

# ============================================================
# Unit conversion
# ============================================================

"""
    convert_coordinate(values, colname) -> (converted_values, new_label)

Converts temperature (K -> degrees C) and pressure (bar or kbar -> GPa).
Returns values and colname unchanged if no known pattern matches.
"""
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

"""
    calc_obs_err(apfu_obs, analysis_type) -> Vector{Float64}

Estimate analytical uncertainties from a.p.f.u. values using empirical
power-law relationships calibrated for each technique:

  EDS:      sigma = clamp(0.0703 * x^0.3574,  0.010, 0.10)
  WDS map:  sigma = clamp(0.0434 * x^0.3451,  0.005, 0.05)
  WDS spot: sigma = clamp(0.023  * x^0.2772,  0.005, 0.05)
"""
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

"""
    q_elem(apfu_obs, model_val, obs_err) -> Float64

Quality factor for a single element at a single grid point. Range: [0, 100].

Formula (Duesterhoeft & Lanari 2020):
    Q = 100 * |1 - clamp(|obs - model| - sigma, 0, 6*sigma) / (6*sigma)| ^ (model + 1)
"""
@inline function q_elem(apfu_obs::Float64, model_val::Float64, obs_err::Float64)::Float64
    sigma = max(obs_err, MIN_ERROR_THRESHOLD)
    diff  = abs(apfu_obs - model_val)
    num   = clamp(diff - sigma / MIN_COEFF_ERR, 0.0, MAX_COEFF_ERR * sigma)
    return 100.0 * abs(1.0 - num / (MAX_COEFF_ERR * sigma))^(model_val + 1.0)
end

"""
    q_phase(apfu_obs, model_vals, obs_errs) -> Float64

Quality factor for a phase: mean of q_elem over all elements. Range: [0, 100].
"""
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

"""
    chi2_stat(apfu_obs, model_vals, obs_errs) -> Float64

Chi-squared statistic: sum of [(obs - model)^2 / sigma^2] over all elements.
"""
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

"""
    red_chi2_stat(apfu_obs, model_vals, obs_errs, dof) -> Float64

Reduced chi-squared: chi2 / (dof - 1).
dof is the number of elements; used only when dof > 2.
"""
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

"""
    parse_measurements(df) -> (element_names, apfu_obs, obs_err, phase_names, phase_ids)

Parse a measurements DataFrame into internal arrays.

Expected format:
  - Column names  : "Phase_Element" strings, e.g. "Grt_Mg", "Grt_Ca", "Bt_Fe"
  - Row 1         : observed a.p.f.u. values
  - Row 2         : uncertainties, or NaN/missing to trigger automatic calculation

Returns:
  element_names  Vector{String}  e.g. ["Grt_Mg", "Grt_Ca", "Bt_Fe"]
  apfu_obs       Vector{Float64} observed values
  obs_err        Vector{Float64} uncertainties (empty if all NaN -> auto)
  phase_names    Vector{String}  unique phases in order of appearance
  phase_ids      Vector{Int}     index into phase_names for each element
"""
function parse_measurements(df::DataFrame)
    element_names = String.(names(df))
    apfu_obs      = Float64.(Vector(df[1, :]))
    err_row       = Vector(df[2, :])

    obs_err = if all(x -> ismissing(x) || (x isa Number && isnan(Float64(x))), err_row)
        Float64[]
    else
        Float64.(coalesce.(err_row, NaN))
    end

    phase_names = String[]
    phase_ids   = Int[]
    for name in element_names
        parts = split(name, "_"; limit=2)
        phase = length(parts) >= 2 ? String(parts[1]) : name
        idx   = findfirst(==(phase), phase_names)
        if isnothing(idx)
            push!(phase_names, phase)
            push!(phase_ids, length(phase_names))
        else
            push!(phase_ids, idx)
        end
    end

    return element_names, apfu_obs, obs_err, phase_names, phase_ids
end

"""
    extract_model_data(model_df, element_names, x_col, y_col)
        -> (x, y, x_label, y_label, model_matrix)

Build coordinate vectors and a (n_points x n_elements) model composition matrix
from the MAGEMin/Perple_X output DataFrame.

Column matching is case-insensitive. Missing phases at a grid point are NaN.

# NOTE FOR MAGEMIN INTEGRATION
# -----------------------------
# This function currently expects a DataFrame loaded from a CSV file.
# When IntersecT is called directly from MAGEMin, replace this function
# (or add an overload) that accepts MAGEMin's internal data structure
# directly, bypassing CSV I/O entirely. The model_matrix output format
# must remain unchanged: Float64 matrix, NaN for absent phases.
"""
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
        end
        # column not found -> stays NaN for all points (phase never present)
    end

    return x, y, x_label, y_label, model_matrix
end

# ============================================================
# I/O helpers for standalone use
# ============================================================

"""
    load_model_csv(path; x_col="T(K)", y_col="P(bar)") -> DataFrame

Load a MAGEMin or Perple_X model output CSV file into a DataFrame.
The file must have a header row with column names.

# NOTE FOR MAGEMIN INTEGRATION
# -----------------------------
# In production, the model output is provided directly by MAGEMin as an
# internal data structure. This function is only needed for standalone use
# (testing, benchmarking, the Python-compatible workflow).
# Replace the call to load_model_csv with whatever MAGEMin exposes.
"""
function load_model_csv(path::String)::DataFrame
    return CSV.read(path, DataFrame)
end

"""
    load_measurements_csv(path) -> DataFrame

Load a measurements CSV file into the format expected by parse_measurements.

Expected CSV format (matches the MAGEMin bulk composition CSV style):
  - Header row  : column names = "Phase_Element" labels (e.g. Grt_Mg, Grt_Ca)
  - Row 1       : observed a.p.f.u. values
  - Row 2       : uncertainties, or leave blank/NaN for automatic calculation

Example:
    Grt_Mg, Grt_Ca, Grt_Fe, Bt_Fe
    1.20,   0.80,   0.95,   2.10
    0.05,   0.04,   0.04,   0.10

# NOTE FOR MAGEMIN INTEGRATION
# -----------------------------
# The GUI will collect measurements interactively. When that is in place,
# this function can be removed or kept as a fallback for batch processing.
"""
function load_measurements_csv(path::String)::DataFrame
    return CSV.read(path, DataFrame; header=true)
end

# ============================================================
# Main calculation pipeline
# ============================================================

"""
    run_intersect(model_df, measurements_df;
                  x_col="T(K)", y_col="P(bar)",
                  analysis_type="WDS spot") -> IntersecTResult

Run the full IntersecT quality factor analysis.

# Arguments
- model_df        DataFrame from MAGEMin/Perple_X. Must have coordinate columns
                  and one column per "Phase_Element" (e.g. "Grt_Mg", "Bt_Fe").
                  Use NaN or missing where a phase is absent at a grid point.
- measurements_df DataFrame: column names = "Phase_Element", row 1 = observed
                  a.p.f.u., row 2 = uncertainties (NaN triggers auto-calculation).
- x_col           Name of the x-coordinate column in model_df.
- y_col           Name of the y-coordinate column in model_df.
- analysis_type   Analytical technique for auto-uncertainty if not provided:
                  "EDS", "WDS map", or "WDS spot".

# NaN policy
- Qcmp_phase[i,p] and redchi2_phase[i,p]: NaN only where phase p is absent at point i.
  Other phases at the same point are computed normally.
- Qcmp_elem[i,j]: NaN where the element's phase is absent at point i.
- Qcmp_unweighted[i], Qcmp_weighted[i], redchi2_tot[i]: NaN if ANY phase is absent
  at point i. These outputs are only meaningful where the full assemblage is stable.

# Parallelisation
The main loop over grid points uses Threads.@threads.
Launch Julia with `julia --threads auto` to use all available CPU cores.
"""
function run_intersect(
    model_df        :: DataFrame,
    measurements_df :: DataFrame;
    x_col           :: String = "T(K)",
    y_col           :: String = "P(bar)",
    analysis_type   :: String = "WDS spot"
)::IntersecTResult

    # --- parse measurements ---
    element_names, apfu_obs, obs_err, phase_names, phase_ids =
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

    return IntersecTResult(
        x, y,
        x_label, y_label,
        phase_names,
        element_names,
        Qcmp_elem_mat,
        Qcmp_phase_mat,
        redchi2_phase_mat,
        redchi2_tot_vec,
        Qcmp_unweighted,
        Qcmp_weighted_vec,
        min_redchi2,
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

end  # module IntersecT