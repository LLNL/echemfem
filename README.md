# EchemFEM

[![DOI](https://zenodo.org/badge/513600791.svg)](https://zenodo.org/badge/latestdoi/513600791)

This code provides Finite Element solvers for electrochemical transport.
Both continuous Galerkin (CG) and discontinuous Galerkin (DG) schemes are provided. The DG scheme for electroneutral Nernst-Planck is described in [Roy et al., 2023](https://doi.org/10.1016/j.jcp.2022.111859). The CG scheme using SUPG for the advection-migration term.

The following transport mechanisms are available: diffusion, advection, electromigration. EchemFEM supports both non-porous and porous cases. The ionic potential can either be described using an electroneutrality constraint or a Poisson equation. In the porous case, the electronic potential can be described by a Poisson equation.

LLNL-CODE-837342

## Installation

Please install the open-source finite element library [Firedrake](https://www.firedrakeproject.org/download.html).

To install EchemFEM, simply run the following in the parent echemfem folder:
```
pip install -e .
```
