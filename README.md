# IntersecT.jl

Julia core for **IntersecT** — Quantitative Isopleth Thermobarometry from Thermodynamic Models.
Reference: Nerone et al. (2025), doi:[10.1016/j.cageo.2025.105949](https://doi.org/10.1016/j.cageo.2025.105949)

Computes quality factors (Qcmp) and reduced chi-squared statistics to assess the agreement between measured and modelled mineral compositions, following the approach of Duesterhoeft & Lanari (2020), doi:[10.1111/jmg.12538](https://doi.org/10.1111/jmg.12538).

This package implements the computational core only. Visualisation and user interface are handled externally (MAGEMin GUI).

---

## Installation

This package is not yet registered. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/YOUR-USERNAME/IntersecT.jl")
```

Or clone the repository and activate the local environment:

```julia
using Pkg
Pkg.activate("path/to/IntersecT.jl")
Pkg.instantiate()
```

---

## Usage

```julia
using IntersecT

result = IntersecT.run_intersect(
    model_df,       # DataFrame: one row per grid point, columns "Phase_Element" in a.p.f.u.
    measurements_df;# DataFrame: row 1 = observed values, row 2 = uncertainties or "auto"
    x_col = "T(K)", # name of the x-coordinate column in model_df
    y_col = "P(bar)",
    analysis_type = "WDS map"  # required if row 2 of measurements_df is "auto"
)
```

### Input format — model DataFrame

One row per grid point. Coordinate columns plus one column per phase-element combination:

| T(K)  | P(bar)  | Grt_Mg | Grt_Ca | Grt_Fe | Ms_Si | Ms_Al | ... |
|-------|---------|--------|--------|--------|-------|-------|-----|
| 873.0 | 10000.0 | 1.20   | 0.80   | 0.95   | 3.57  | 1.96  | ... |
| 873.0 | 20000.0 | NaN    | NaN    | NaN    | 3.61  | 1.94  | ... |

Use `NaN` or `missing` where a phase is absent at a grid point.

### Input format — measurements DataFrame

Column names are `Phase_Element` labels. Row 1 = observed a.p.f.u., row 2 = uncertainties. Write `auto` in the first cell of row 2 (leaving all other cells empty) to trigger automatic uncertainty estimation from `analysis_type`. Row 2 must be either fully numeric or `auto`; mixed values are not accepted.

| Grt_Mg | Grt_Ca | Grt_Fe | Ms_Si | Ms_Al |
|--------|--------|--------|-------|-------|
| 1.20   | 0.80   | 0.95   | 3.57  | 1.96  |
| 0.05   | 0.04   | 0.04   | 0.21  | 0.21  |

Or, using automatic uncertainty estimation:

| Grt_Mg | Grt_Ca | Grt_Fe | Ms_Si | Ms_Al |
|--------|--------|--------|-------|-------|
| 1.20   | 0.80   | 0.95   | 3.57  | 1.96  |
| auto   |        |        |       |       |

### Output

`run_intersect` returns an `IntersecTResult` struct:

| Field | Type | Description |
|-------|------|-------------|
| `x`, `y` | `Vector{Float64}` | Grid coordinates after unit conversion |
| `x_label`, `y_label` | `String` | Axis labels (e.g. `"T(C)"`, `"P(GPa)"`) |
| `phase_names` | `Vector{String}` | Phase names in order |
| `element_names` | `Vector{String}` | `"Phase_Element"` labels in order |
| `Qcmp_elem` | `Matrix{Float64}` | Quality factor per element `(n_points × n_elements)` |
| `Qcmp_phase` | `Matrix{Float64}` | Quality factor per phase `(n_points × n_phases)` |
| `redchi2_phase` | `Matrix{Float64}` | Reduced χ² per phase `(n_points × n_phases)` |
| `redchi2_tot` | `Vector{Float64}` | Total reduced χ² `(n_points,)` |
| `Qcmp_unweighted` | `Vector{Float64}` | Unweighted total Qcmp `(n_points,)` |
| `Qcmp_weighted` | `Vector{Float64}` | Weighted total Qcmp `(n_points,)` |
| `min_redchi2` | `Vector{Float64}` | Minimum reduced χ² per phase, used for weighting `(n_phases,)` |

### NaN policy

- `Qcmp_phase[i, p]` and `redchi2_phase[i, p]` are `NaN` only where phase `p` is absent at point `i`. Other phases at the same point are computed normally.
- `Qcmp_weighted[i]`, `Qcmp_unweighted[i]`, and `redchi2_tot[i]` are `NaN` where **any** phase is absent. These outputs are only defined within the full assemblage stability field.

---

## Coordinate conversion

The following column name patterns are recognised and converted automatically:

| Pattern | Conversion | Output label |
|---------|-----------|--------------|
| `T(K)`, `T[K]`, `T_K` | − 273 | `T(C)` |
| `T(C)`, `T[C]`, `T_C` | none | `T(C)` |
| `P(bar)`, `P[bar]`, `P_bar` | ÷ 10000 | `P(GPa)` |
| `P(kbar)`, `P[kbar]`, `P_kbar`, `kbar` | ÷ 10 | `P(GPa)` |
| anything else | none | original name |

---

## Parallelisation

The main loop over grid points uses `Threads.@threads`. Launch Julia with `--threads auto` to use all available CPU cores:

```
julia --threads auto
```

---

## Running the tests

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("test/runtests.jl")
```

---

## Dependencies

- [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl)
- [CSV.jl](https://github.com/JuliaData/CSV.jl) — only for standalone I/O helpers

---

## References

Duesterhoeft, E. & Lanari, P. (2020). Iterative thermodynamic modelling — Part 1: A new forward approach. *Journal of Metamorphic Geology*, 38, 733–752. doi:[10.1111/jmg.12538](https://doi.org/10.1111/jmg.12538)

Nerone, S, Lanari, P., Dominguez, H., Forshaw, J. B., Groppo, C., Rolfo, F. (2025). IntersecT: a Python script for quantitative isopleth thermobarometry of equilibrium and disequilibrium systems. *Computers & Geosciences*, 202, 105949. doi:[10.1016/j.cageo.2025.105949](https://doi.org/10.1016/j.cageo.2025.105949)

Warr, L. N. (2021). IMA-CNMNC approved mineral symbols. *Mineralogical Magazine*, 85, 291–320. doi:[10.1180/mgm.2021.43](https://doi.org/10.1180/mgm.2021.43)
