# NLVE_SCHAPERY – nonlinear viscoelastic Schapery model for Abaqus (UMAT & VUMAT)

A user material for Abaqus/Standard (UMAT) and Abaqus/Explicit (VUMAT) that describes the
nonlinear viscoelastic, stress-dependent creep behavior of thermoplastics after Schapery.
It was developed and calibrated for polyamide 6 (PA6) at different climatic conditions.

<br/><br/>

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/Extraweich/Schapery_visco_elasticity/blob/main/Figures/Kelvin-Voigt3_light.svg?raw=true">
    <img alt="Schematic of implemented homogenization procedures" src="https://github.com/Extraweich/Schapery_visco_elasticity/blob/main/Figures/Kelvin-Voigt3.svg?raw=true", width="50%">
  </picture>
</p>

<br/><br/>

The model and its implementation were developed during my PhD project at the
Karlsruhe Institute of Technology (KIT). Theory, derivation, experimental calibration and
verification are documented in the dissertation:

> Nicolas Christ: On the influence of temperature and humidity on interfaces in carbon fiber reinforced polyamide 6. Doctoral thesis, Karlsruhe Institute of
> Technology (KIT), 2026. DOI: [10.5445/IR/1000194396](https://publikationen.bibliothek.kit.edu/1000194396)

**If you use this code in published work, please cite the dissertation** (see
[Citation](#citation)).

## Features

- Generalized Kelvin-Voigt model (spring + 3 Kelvin-Voigt elements) extended by the
  Schapery functions g0, g1·g2 and a_σ of the von Mises stress
- Formulation for large stress increments (e.g. suddenly applied creep loads), implicit
  stress update with a local Newton-Raphson iteration
- **One source file for both solvers**: the same `src/nlve_schapery.f` is passed to
  Abaqus/Standard and Abaqus/Explicit
- 3D continuum, plane strain and axisymmetric elements
- Identical state-variable layout in UMAT and VUMAT: an implicit analysis can be continued
  in Abaqus/Explicit with `*IMPORT` (example 03)
- Consistent tangent (DDSDDE) for quadratic convergence in Abaqus/Standard
- Large rotations (NLGEOM): tensor-valued state variables are rotated
- No external dependencies (no LAPACK, no solver-specific Abaqus utilities)
- Material-point test suite that runs without Abaqus

## Repository structure

```
nlve-schapery/
├── src/
│   └── nlve_schapery.f              UMAT + VUMAT (+ SDVINI), the only file Abaqus needs
├── examples/
│   ├── 01_creep_3D/                 uniaxial creep test, one C3D8 element, 160 h
│   ├── 02_tube_axisymmetric/        thick-walled tube under internal pressure (CAX4)
│   ├── 03_implicit_to_explicit/     creep in Standard, continued in Explicit via *IMPORT
│   └── reference/
│       └── analytical_creep.py      closed-form solution (Eq. 52) for examples 01 and 03
├── tests/
│   └── nlve_tests.f90               material-point tests (gfortran, no Abaqus needed)
├── .github/workflows/tests.yml      runs the tests on every push
├── CITATION.cff
├── LICENSE
└── README.md
```

## Usage

```
abaqus job=<job> user=src/nlve_schapery.f                  # Abaqus/Standard
abaqus job=<job> user=src/nlve_schapery.f double=both      # Abaqus/Explicit
```

Requirements: Abaqus with a configured Intel Fortran compiler (the file starts with the
Intel directive `!DEC$ FREEFORM`). Abaqus/Explicit must run in double precision
(the default `double=explicit` or `double=both`).

### Material definition

```
*Depvar
 29,
*User Material, constants=23, unsymm      ** UNSYMM only in Abaqus/Standard
```

| PROPS | Meaning |
|---|---|
| 1 | E0 – instantaneous modulus |
| 2 | ν – Poisson's ratio (identical for all springs) |
| 3, 4 | E1, τ1 – Kelvin-Voigt element 1 |
| 5, 6 | E2, τ2 – Kelvin-Voigt element 2 |
| 7, 8 | E3, τ3 – Kelvin-Voigt element 3 |
| 9, 10–13 | type and parameters p1–p4 of g0 |
| 14, 15–18 | type and parameters p1–p4 of g1·g2 |
| 19, 20–23 | type and parameters p1–p4 of a_σ |

Function types of f(σ̄), σ̄ = von Mises stress:

| Type | f(σ̄) |
|---|---|
| 0 | 1 (linear viscoelasticity) |
| 1 | p4 + (p1 − p4) / (1 + (σ̄/p3)^p2) (sigmoid with a, b, c, d = p1…p4) |
| 2 | p1 + p2 · exp(−(σ̄ − p3)² / (2 p4²)) (Gaussian) |

Units are arbitrary but must be consistent; the functions take σ̄ in the stress unit of
the moduli. The examples use mm, N, MPa, s and the parameters of PA6 at 23 °C and 10 % r.H.
from the thesis.

With all three types set to 0 the model reduces to linear viscoelasticity, which can be
verified against the built-in `*VISCOELASTIC` model of Abaqus (see the thesis).

### State variables

Tensors in Abaqus/Standard notation (order 11, 22, 33, 12, 13, 23, engineering shear
strains), identical for UMAT and VUMAT. In 2D the components 13 and 23 are zero.

| SDV | Meaning |
|---|---|
| 1–6, 7–12, 13–18 | inherited strains of Kelvin-Voigt elements 1, 2, 3 |
| 19, 20, 21 | g0, g1·g2, a_σ at the end of the increment |
| 22 | relative residual of the local Newton iteration |
| 23 | von Mises (indicator) stress |
| 24–29 | total strain (used by the VUMAT and for `*IMPORT`) |

`*INITIAL CONDITIONS, TYPE=SOLUTION, USER` (SDVINI) is optional.

### Limitations

- Plane stress, shell and beam elements are not supported. Abaqus/Explicit passes the same
  dimensions for plane stress as for plane strain, so this cannot be detected in the VUMAT.
- Small-strain formulation; with NLGEOM the model is applied to the corotational strain
  measures of Abaqus.
- Temperature and humidity enter only through the parameter set (one set per climate).

## Examples

| Example | Shows | Check |
|---|---|---|
| 01_creep_3D | uniaxial creep, 30 MPa for 160 h | creep strain vs. Eq. 52: `python examples/reference/analytical_creep.py --stress 30` |
| 02_tube_axisymmetric | multiaxial, inhomogeneous stress state | stresses redistribute over time; with function types 0 they stay constant |
| 03_implicit_to_explicit | Standard → Explicit with `*IMPORT` | both jobs form one continuous creep curve |

Each input file states in its header how to run it and what to expect.

## Tests

The material routine can be tested without Abaqus:

```
gfortran -O2 -ffree-form src/nlve_schapery.f tests/nlve_tests.f90 -o nlve_tests
./nlve_tests
```

The tests cover creep against the closed-form solution (uniaxial and shear, 3D and
axisymmetric), identity of UMAT and VUMAT on random multiaxial strain paths, the
Standard → Explicit hand-over, rigid-body rotations, zero time and strain increments, the
VUMAT data-check call and the consistent tangent.

## Citation

If you use this code, please cite the dissertation:

```bibtex
@phdthesis{Christ_Thesis_2026,
  author       = {Christ, Nicolas},
  year         = {2026},
  title        = {On the influence of temperature and humidity on interfaces in carbon fiber reinforced polyamide 6},
  doi          = {10.5445/IR/1000194396},
  publisher    = {{Karlsruher Institut für Technologie (KIT)}},
  keywords     = {nonlinear viscoelasticity, Schapery, single fiber pull-out test, climbing drum peel test},
  pagetotal    = {261},
  school       = {Karlsruher Institut für Technologie (KIT)},
  language     = {english}
}
```

GitHub shows the same information via "Cite this repository" (`CITATION.cff`).

## Development history

The model and the original Abaqus implementation were developed by the author during the PhD project, without AI assistance. That version is printed in the appendix of the dissertation. For the public release, the code was revised with the help of an AI assistant (Claude, Anthropic): the separate UMAT/VUMAT and 3D/axisymmetric versions (only results from the 3D UMAT were published in the thesis) were merged into one file, the parameter functions were moved to the material constants, errors in two of the original versions were corrected, the numerical robustness was improved, and the tests and examples were added. The model equations are unchanged. All changes were reviewed by the author, and results obtained with the original version are reproduced within the solver tolerance.

## License

MIT License, see [LICENSE](LICENSE). You may use, modify and redistribute the code freely,
also commercially. The citation of the dissertation is a request, not a license condition.
