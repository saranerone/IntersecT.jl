# IntersecT.jl

Julia core for **IntersecT** — Quantitative Isopleth Thermobarometry from Thermodynamic Models.
Reference: Nerone et al. (2025), doi: [10.1016/j.cageo.2025.105949](https://doi.org/10.1016/j.cageo.2025.105949)

Computes quality factors (Qcmp) and reduced chi-squared statistics to assess the agreement between measured and modelled mineral compositions, following the approach of Duesterhoeft & Lanari (2020), doi: [10.1111/jmg.12538](https://doi.org/10.1111/jmg.12538).

This package implements the computational core only. Visualisation and user interface are handled externally (MAGEMin GUI).

---

## Installation

Users of the MAGEMin GUI do not need to install this package separately: it is
installed as a dependency of MAGEMinApp. Install it directly only for standalone
use.

This package is not yet registered. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/saranerone/IntersecT.jl")
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
    measurements_df;# DataFrame: observed values, uncertainties or "auto", optional domain row
    x_col = "T(K)", # name of the x-coordinate column in model_df
    y_col = "P(bar)",
    analysis_type = "WDS map"  # required if the uncertainty row is "auto"
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

Column names are `Phase_Element` labels. Row 1 = observed a.p.f.u., row 2 = uncertainties. Write `auto` in the first cell of the uncertainty row (leaving all other cells empty) to trigger automatic uncertainty estimation from `analysis_type`. The uncertainty row must be either fully numeric or `auto`; mixed values are not accepted.

| Grt_Mg | Grt_Ca | Grt_Fe | Ms_Si | Ms_Al |
|--------|--------|--------|-------|-------|
| 1.20   | 0.80   | 0.95   | 3.57  | 1.96  |
| 0.05   | 0.04   | 0.04   | 0.21  | 0.21  |

Or, using automatic uncertainty estimation:

| Grt_Mg | Grt_Ca | Grt_Fe | Ms_Si | Ms_Al |
|--------|--------|--------|-------|-------|
| 1.20   | 0.80   | 0.95   | 3.57  | 1.96  |
| auto   |        |        |       |       |

### Input format — compositional domains

Compositional domains (e.g. garnet core and rim, biotite syn- or post- main foliation) are declared in an optional row placed directly below the header, above the observed values:

| Grt_Mg | Grt_Ca | Grt_Mg | Grt_Ca | Bt_Mg | Bt_Fe |
|--------|--------|--------|--------|-------|-------|
| core   | core   | rim    | rim    |       |       |
| 1.20   | 0.80   | 0.60   | 0.95   | 1.10  | 1.40  |
| 0.05   | 0.04   | 0.05   | 0.04   | 0.05  | 0.05  |

Domain labels are free text. Empty cells mean that the column carries no domain. The row is detected automatically: a first row that is fully numeric is read as the observed values, and existing files without domains are therefore unaffected. `auto` is not accepted as a domain label.

Domains are kept separate at every stage of the calculation. A phase with a domain is displayed with its qualified name (`Grt core`), while the model is always queried with the base name (`Grt`). Duplicated column names are expected in files with domains, and the suffix added by CSV.jl (`Grt_Mg_1`) is removed internally.

Two domains of the same phase cannot enter the same run: domains grown at different conditions do not constitute an equilibrium assemblage, which is what the quality factor assumes. Each domain is instead evaluated in a separate run.

```julia
labels, bases = IntersecT.list_measured_phases(measurements_df)
# labels = ["Grt core", "Grt rim", "Bt"]
# bases  = ["Grt", "Grt", "Bt"]

core_df = IntersecT.select_phases(measurements_df, ["Grt core", "Bt"])
rim_df  = IntersecT.select_phases(measurements_df, ["Grt rim", "Bt"])

res_core = IntersecT.run_intersect(model_df, core_df)
res_rim  = IntersecT.run_intersect(model_df, rim_df)
```

`list_measured_phases` returns the qualified names, to be displayed in the interface, and the base names, to be intersected with the stable phases of the diagram. `select_phases` extracts the columns of the selected phases and throws an error if two domains of the same phase are selected together.

Absolute Q*cmp values are not comparable between runs that include different phases or different elements, because the weights are normalised within each run. Comparable is the position of the maximum, not its value.

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
| `n_elements_per_phase` | `Vector{Int}` | Number of measured elements per phase `(n_phases,)` |
| `phase_bases` | `Vector{String}` | Phase names without domain label `(n_phases,)` |
| `analysis_type` | `String` | Analysis type used for the run |

### NaN policy

- `Qcmp_phase[i, p]` and `redchi2_phase[i, p]` are `NaN` only where phase `p` is absent at point `i`. Other phases at the same point are computed normally.
- `Qcmp_weighted[i]`, `Qcmp_unweighted[i]`, and `redchi2_tot[i]` are `NaN` where **any** phase is absent. These outputs are only defined within the full assemblage stability field.
- A measured column with no matching column in the model output produces a warning and an empty map.

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
