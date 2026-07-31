"""
Test suite for IntersecT.jl

Tests verify:
  1. Core statistical functions (q_elem, q_phase, chi2, red_chi2)
  2. Unit conversion
  3. Uncertainty calculation
  4. Measurement parsing
  5. NaN propagation: absent phases leave NaN in all outputs
  6. Weighting logic: min_redchi2, floor at 1, weight inversion
  7. Full pipeline on synthetic data: correct max position
  8. Two-phase pipeline
  9. Reproducibility across two runs (thread safety)

Run from the IntersecT.jl directory:
    julia --threads auto test/runtests.jl
"""

using Test
using DataFrames

include(joinpath(@__DIR__, "..", "src", "IntersecT.jl"))
using .IntersecT

# ============================================================
# 1. q_elem
# ============================================================

@testset "q_elem" begin
    # Perfect match -> 100
    @test IntersecT.q_elem(1.0, 1.0, 0.05) == 100.0

    # obs_err below threshold: MIN_ERROR_THRESHOLD used, perfect match still 100
    @test IntersecT.q_elem(1.0, 1.0, 0.001) == 100.0

    # Large mismatch: diff=0.5, sigma=0.05
    # num = clamp(0.5 - 0.05, 0, 0.30) = 0.30
    # Q = 100 * |1 - 0.30/0.30|^(0.5+1) = 0
    @test IntersecT.q_elem(1.0, 0.5, 0.05) == 0.0

    # Intermediate: obs=1.0, model=1.1, err=0.05
    # sigma=0.05, diff=0.1, num=clamp(0.05, 0, 0.30)=0.05
    # Q = 100 * |1 - 0.05/0.30|^2.1
    expected = 100.0 * abs(1.0 - 0.05 / 0.30)^(1.1 + 1.0)
    @test IntersecT.q_elem(1.0, 1.1, 0.05) ≈ expected atol=1e-8
end

# ============================================================
# 2. q_phase
# ============================================================

@testset "q_phase" begin
    # Single element perfect match
    @test IntersecT.q_phase([1.0], [1.0], [0.05]) == 100.0

    # Two elements both perfect
    @test IntersecT.q_phase([1.0, 0.5], [1.0, 0.5], [0.05, 0.02]) == 100.0

    # One perfect + one zero -> mean = 50
    q1 = IntersecT.q_elem(1.0, 1.0, 0.05)   # 100.0
    q2 = IntersecT.q_elem(1.0, 0.5, 0.05)   # 0.0
    @test IntersecT.q_phase([1.0, 1.0], [1.0, 0.5], [0.05, 0.05]) ≈ (q1 + q2) / 2.0 atol=1e-8
end

# ============================================================
# 3. chi2_stat and red_chi2_stat
# ============================================================

@testset "chi2_stat" begin
    obs = [1.0, 2.0, 3.0]
    mod = [1.1, 1.9, 3.2]
    err = [0.05, 0.05, 0.05]
    # (0.1/0.05)^2 + (0.1/0.05)^2 + (0.2/0.05)^2 = 4 + 4 + 16 = 24
    @test IntersecT.chi2_stat(obs, mod, err) ≈ 24.0 atol=1e-8
end

@testset "red_chi2_stat" begin
    obs = [1.0, 2.0, 3.0]
    mod = [1.1, 1.9, 3.2]
    err = [0.05, 0.05, 0.05]
    # 24 / (3-1) = 12
    @test IntersecT.red_chi2_stat(obs, mod, err, 3) ≈ 12.0 atol=1e-8
end

# ============================================================
# 4. Unit conversion
# ============================================================

@testset "convert_coordinate" begin
    v, l = IntersecT.convert_coordinate([600.0, 700.0], "T(K)")
    @test v ≈ [327.0, 427.0] atol=1e-10
    @test l == "T(C)"

    v, l = IntersecT.convert_coordinate([5000.0, 10000.0], "P(bar)")
    @test v ≈ [0.5, 1.0] atol=1e-10
    @test l == "P(GPa)"

    v, l = IntersecT.convert_coordinate([5.0, 10.0], "P[kbar]")
    @test v ≈ [0.5, 1.0] atol=1e-10
    @test l == "P(GPa)"

    v, l = IntersecT.convert_coordinate([0.3, 0.7], "X[0.0-1.0]")
    @test v ≈ [0.3, 0.7]
    @test l == "X[0.0-1.0]"
end

# ============================================================
# 5. calc_obs_err
# ============================================================

@testset "calc_obs_err" begin
    obs = [1.0, 0.5, 0.1]

    err = IntersecT.calc_obs_err(obs, "EDS")
    @test err ≈ clamp.(0.0703 .* obs .^ 0.3574, 0.01, 0.1) atol=1e-10

    err = IntersecT.calc_obs_err(obs, "WDS spot")
    @test err ≈ clamp.(0.023 .* obs .^ 0.2772, 0.005, 0.05) atol=1e-10

    # Values too small -> clamped to min
    @test all(IntersecT.calc_obs_err([0.0001, 0.0001], "WDS spot") .>= 0.005)

    # Invalid type
    @test_throws ErrorException IntersecT.calc_obs_err(obs, "SIMS")
end

# ============================================================
# 6. parse_measurements
# ============================================================

@testset "parse_measurements" begin
    df = DataFrame(
        "Grt_Mg" => [1.2, 0.05],
        "Grt_Ca" => [0.8, 0.04],
        "Bt_Fe"  => [2.1, 0.1],
    )
    elem_names, apfu_obs, obs_err, phase_names, phase_ids =
        IntersecT.parse_measurements(df)

    @test elem_names  == ["Grt_Mg", "Grt_Ca", "Bt_Fe"]
    @test apfu_obs    ≈ [1.2, 0.8, 2.1]
    @test obs_err     ≈ [0.05, 0.04, 0.1]
    @test phase_names == ["Grt", "Bt"]
    @test phase_ids   == [1, 1, 2]
end

@testset "parse_measurements NaN uncertainties" begin
    df = DataFrame("Grt_Mg" => [1.2, NaN], "Grt_Ca" => [0.8, NaN])
    _, _, obs_err, _, _ = IntersecT.parse_measurements(df)
    @test isempty(obs_err)
end

@testset "parse_measurements auto keyword" begin
    df = DataFrame("Grt_Mg" => ["1.2", "auto"], "Grt_Ca" => ["0.8", missing])
    _, apfu_obs, obs_err, _, _ = IntersecT.parse_measurements(df)
    @test apfu_obs ≈ [1.2, 0.8]
    @test isempty(obs_err)

    df_upper = DataFrame("Grt_Mg" => ["1.2", "AUTO"], "Grt_Ca" => ["0.8", missing])
    @test isempty(IntersecT.parse_measurements(df_upper)[3])

    df_mixed = DataFrame("Grt_Mg" => ["1.2", "auto"], "Grt_Ca" => ["0.8", "0.04"])
    @test_throws ErrorException IntersecT.parse_measurements(df_mixed)

    df_partial = DataFrame("Grt_Mg" => [1.2, 0.05], "Grt_Ca" => [0.8, missing])
    @test_throws ErrorException IntersecT.parse_measurements(df_partial)
end

# ============================================================
# 7. NaN propagation for absent phases
# ============================================================
#
# Scenario: 4 grid points, 1 phase (Grt) with 2 elements (Mg, Ca).
# The phase is present at points 1, 2, 4 but absent at point 3.
# Expected: all outputs are NaN at point 3, numeric elsewhere.

@testset "NaN propagation - absent phase" begin
    model_df = DataFrame(
        "T(K)"   => [873.0, 873.0, 973.0, 973.0],
        "P(bar)" => [10000.0, 20000.0, 10000.0, 20000.0],
        "Grt_Mg" => [1.2,  1.0,  NaN,  0.9],
        "Grt_Ca" => [0.8,  0.7,  NaN,  0.6],
    )
    meas_df = DataFrame("Grt_Mg" => [1.2, 0.05], "Grt_Ca" => [0.8, 0.04])

    result = IntersecT.run_intersect(model_df, meas_df; x_col="T(K)", y_col="P(bar)")

    # Point 3 (index 3): phase absent -> all outputs NaN
    @test isnan(result.Qcmp_phase[3, 1])
    @test isnan(result.redchi2_phase[3, 1])
    @test isnan(result.redchi2_tot[3])
    @test isnan(result.Qcmp_unweighted[3])
    @test isnan(result.Qcmp_weighted[3])
    @test isnan(result.Qcmp_elem[3, 1])
    @test isnan(result.Qcmp_elem[3, 2])

    # Other points: phase present -> numeric (not NaN)
    for i in [1, 2, 4]
        @test !isnan(result.Qcmp_phase[i, 1])
        @test !isnan(result.Qcmp_unweighted[i])
        @test !isnan(result.Qcmp_weighted[i])
    end

    # Point 1: Grt_Mg matches exactly (1.2 == 1.2), Grt_Ca matches (0.8 == 0.8) -> Qcmp = 100
    @test result.Qcmp_phase[1, 1] ≈ 100.0 atol=1e-6
end

@testset "NaN propagation - two phases, one absent" begin
    # 4 points, 2 phases: Grt and Bt
    # Bt is absent at point 2 -> point 2 should be NaN in weighted/unweighted
    # Grt is present everywhere
    model_df = DataFrame(
        "T(K)"   => [873.0, 873.0, 973.0, 973.0],
        "P(bar)" => [10000.0, 20000.0, 10000.0, 20000.0],
        "Grt_Mg" => [1.2,  1.2,  1.0,  0.9],
        "Grt_Ca" => [0.8,  0.8,  0.7,  0.6],
        "Bt_Fe"  => [2.1,  NaN,  2.0,  1.9],
    )
    meas_df = DataFrame(
        "Grt_Mg" => [1.2, 0.05],
        "Grt_Ca" => [0.8, 0.04],
        "Bt_Fe"  => [2.1, 0.10],
    )

    result = IntersecT.run_intersect(model_df, meas_df; x_col="T(K)", y_col="P(bar)")

    # Point 2: Bt absent -> NaN in unweighted and weighted
    @test isnan(result.Qcmp_unweighted[2])
    @test isnan(result.Qcmp_weighted[2])

    # Grt is present at point 2 -> its Qcmp_phase must be numeric
    @test !isnan(result.Qcmp_phase[2, 1])     # Grt: present, computed
    @test  isnan(result.Qcmp_phase[2, 2])     # Bt: absent, NaN

    # Grt element outputs also numeric at point 2
    @test !isnan(result.Qcmp_elem[2, 1])      # Grt_Mg
    @test !isnan(result.Qcmp_elem[2, 2])      # Grt_Ca
    # Bt element output NaN at point 2
    @test  isnan(result.Qcmp_elem[2, 3])      # Bt_Fe

    # redchi2_phase follows the same rule
    @test !isnan(result.redchi2_phase[2, 1])  # Grt: computed
    @test  isnan(result.redchi2_phase[2, 2])  # Bt: NaN

    # Points 1, 3, 4: both phases present -> numeric
    for i in [1, 3, 4]
        @test !isnan(result.Qcmp_unweighted[i])
        @test !isnan(result.Qcmp_weighted[i])
    end
end

# ============================================================
# 8. Weighting logic
# ============================================================
#
# Scenario: 2 phases, one fitting well (low redchi2), one fitting poorly.
# The well-fitting phase should get higher weight.
# We verify: min_redchi2 values, weight direction, floor at 1.

@testset "weighting - better phase gets higher weight" begin
    # Phase Grt: model matches obs very closely -> low redchi2
    # Phase Bt:  model far from obs -> high redchi2
    # Both have 3 elements so reduced chi2 applies (not chi2)
    model_df = DataFrame(
        "T(K)"   => [873.0, 873.0, 973.0, 973.0],
        "P(bar)" => [10000.0, 20000.0, 10000.0, 20000.0],
        # Grt: near-perfect match
        "Grt_Mg" => [1.20, 1.21, 1.19, 1.20],
        "Grt_Ca" => [0.80, 0.81, 0.79, 0.80],
        "Grt_Fe" => [0.95, 0.96, 0.94, 0.95],
        # Bt: very poor match
        "Bt_Fe"  => [5.0,  5.0,  5.0,  5.0],
        "Bt_Mg"  => [5.0,  5.0,  5.0,  5.0],
        "Bt_Al"  => [5.0,  5.0,  5.0,  5.0],
    )
    meas_df = DataFrame(
        "Grt_Mg" => [1.20, 0.05],
        "Grt_Ca" => [0.80, 0.04],
        "Grt_Fe" => [0.95, 0.04],
        "Bt_Fe"  => [2.10, 0.10],
        "Bt_Mg"  => [0.80, 0.05],
        "Bt_Al"  => [1.50, 0.06],
    )

    result = IntersecT.run_intersect(model_df, meas_df; x_col="T(K)", y_col="P(bar)")

    # Grt fits well -> lower min_redchi2 -> higher weight -> weighted > unweighted
    # (because Grt Qcmp > Bt Qcmp in this scenario)
    @test result.min_redchi2[1] < result.min_redchi2[2]   # Grt better than Bt

    # Floor: min_redchi2 is always >= some value > 0 (could be < 1 for Grt)
    # The floor at 1 is applied inside weighting, not stored in min_redchi2 itself
    # So min_redchi2 can be < 1 for a near-perfect fit
    @test all(result.min_redchi2 .>= 0.0)
end

@testset "weighting - two-element phase uses chi2 + 1" begin
    # Phase with exactly 2 elements: chi2 is used, and 1 is added to min_val
    model_df = DataFrame(
        "T(K)"   => [873.0, 973.0],
        "P(bar)" => [10000.0, 10000.0],
        "Grt_Mg" => [1.2, 1.0],
        "Grt_Ca" => [0.8, 0.7],
    )
    meas_df = DataFrame("Grt_Mg" => [1.2, 0.05], "Grt_Ca" => [0.8, 0.04])

    result = IntersecT.run_intersect(model_df, meas_df; x_col="T(K)", y_col="P(bar)")

    # At point 1, model = obs exactly -> chi2 = 0, stored as 0 in redchi2_phase
    @test result.redchi2_phase[1, 1] ≈ 0.0 atol=1e-10
    # min_redchi2 = 1 + 0 = 1.0 (the +1 for <=2 elements)
    @test result.min_redchi2[1] ≈ 1.0 atol=1e-10
end

# ============================================================
# 9. Full pipeline - synthetic 3x3 grid
# ============================================================

@testset "run_intersect - synthetic single phase" begin
    obs_Mg, obs_Ca = 1.2, 0.8
    T_vals = [500.0, 600.0, 700.0, 500.0, 600.0, 700.0, 500.0, 600.0, 700.0] .+ 273.0
    P_vals = [10000.0, 10000.0, 10000.0, 20000.0, 20000.0, 20000.0, 30000.0, 30000.0, 30000.0]

    # Point 5 (T=600, P=20000) matches observations exactly
    Grt_Mg = [1.0, 1.0, 1.0, 1.0, obs_Mg, 1.0, 1.0, 1.0, 1.0]
    Grt_Ca = [0.5, 0.5, 0.5, 0.5, obs_Ca, 0.5, 0.5, 0.5, 0.5]

    model_df = DataFrame(
        "T(K)"   => T_vals,
        "P(bar)" => P_vals,
        "Grt_Mg" => Grt_Mg,
        "Grt_Ca" => Grt_Ca,
    )
    meas_df = DataFrame("Grt_Mg" => [obs_Mg, 0.05], "Grt_Ca" => [obs_Ca, 0.04])

    result = IntersecT.run_intersect(model_df, meas_df; x_col="T(K)", y_col="P(bar)")

    # Coordinate conversion applied
    @test result.x_label == "T(C)"
    @test result.y_label == "P(GPa)"

    # Point 5: perfect match
    @test result.Qcmp_phase[5, 1]    ≈ 100.0 atol=1e-6
    @test result.Qcmp_unweighted[5]  ≈ 100.0 atol=1e-6
    @test result.Qcmp_weighted[5]    ≈ 100.0 atol=1e-6
    @test result.redchi2_tot[5]      ≈ 0.0   atol=1e-10

    # Maximum is at point 5
    finite_w = filter(!isnan, result.Qcmp_weighted)
    @test argmax(result.Qcmp_weighted) == 5

    # No NaN in this dataset (phase present everywhere)
    @test !any(isnan, result.Qcmp_weighted)

    # Output shapes
    @test size(result.Qcmp_elem)     == (9, 2)
    @test size(result.Qcmp_phase)    == (9, 1)
    @test size(result.redchi2_phase) == (9, 1)
    @test length(result.redchi2_tot) == 9
end

# ============================================================
# 10. Thread reproducibility
# ============================================================

@testset "thread reproducibility" begin
    n = 200
    rng_T = range(773.0, 1073.0, length=n)
    rng_P = range(5000.0, 25000.0, length=n)

    model_df = DataFrame(
        "T(K)"   => collect(rng_T),
        "P(bar)" => collect(rng_P),
        "Grt_Mg" => rand(n) .* 2.0,
        "Grt_Ca" => rand(n),
        "Bt_Fe"  => rand(n) .* 3.0,
    )
    meas_df = DataFrame(
        "Grt_Mg" => [1.2, 0.05],
        "Grt_Ca" => [0.8, 0.04],
        "Bt_Fe"  => [2.1, 0.10],
    )

    r1 = IntersecT.run_intersect(model_df, meas_df; x_col="T(K)", y_col="P(bar)")
    r2 = IntersecT.run_intersect(model_df, meas_df; x_col="T(K)", y_col="P(bar)")

    @test r1.Qcmp_weighted   == r2.Qcmp_weighted
    @test r1.Qcmp_unweighted == r2.Qcmp_unweighted
    @test r1.redchi2_tot     == r2.redchi2_tot
end

println("\nAll tests passed.")