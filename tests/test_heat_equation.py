import pytest
from firedrake import *
from echemfem import TransientEchemSolver
import numpy as np


class HeatEquationSolver(TransientEchemSolver):
    def __init__(self, N, family="DG"):
        mesh = UnitSquareMesh(N, N, quadrilateral=True)
        conc_params = []

        self.time = Constant(0)  # for now, need to define this here
        t = self.time
        x, y = SpatialCoordinate(mesh)
        C1ex = (sin(x) + cos(y)) * exp(t)
        C2ex = (cos(x) + sin(y)) * exp(-t)
        self.C1ex = C1ex
        self.C2ex = C2ex

        def f(C):
            f1 = C1ex - div(grad(C1ex))
            f2 = -C2ex - div(grad(C2ex))
            return [f1, f2]

        conc_params.append({"name": "C1",
                            "diffusion coefficient": 1.0,
                            "bulk": C1ex,
                            })

        conc_params.append({"name": "C2",
                            "diffusion coefficient": 1.0,
                            "bulk": C2ex,
                            })
        physical_params = {"flow": ["diffusion"],
                           "bulk reaction": f,
                           }

        super().__init__(conc_params, physical_params, mesh, family=family)

    def set_boundary_markers(self):
        self.boundary_markers = {"bulk dirichlet": (1, 2, 3, 4,),
                                 }


def test_convergence():
    err_old = 1e6
    for i in range(4):
        solver = HeatEquationSolver(2**(i + 1))
        solver.setup_solver()
        times = np.linspace(0, 1, 1+2**(2*(i+1)))
        solver.solve(times)
        c1, c2 = solver.u.subfunctions
        err = errornorm(solver.C1ex, c1) + errornorm(solver.C2ex, c2)
        assert err < 0.27 * err_old
        err_old = err


def test_convergence_CG():
    err_old = 1e6
    for i in range(5):
        solver = HeatEquationSolver(2**(i + 1), family="CG")
        solver.setup_solver()
        times = np.linspace(0, 1, 1+2**(2*i))
        solver.solve(times)
        c1, c2 = solver.u.subfunctions
        err = errornorm(solver.C1ex, c1) + errornorm(solver.C2ex, c2)
        assert err < 0.26 * err_old
        err_old = err


def test_convergence_BE():
    err_old = 1e6
    solver = HeatEquationSolver(16, family="CG")
    for i in range(5):
        solver.time.assign(0)
        solver.setup_solver()
        solver.u_old.assign(solver.u)
        times = np.linspace(0, 1, 1+2**(i+1))
        solver.solve(times)
        c1, c2 = solver.u.subfunctions
        err = errornorm(solver.C1ex, c1) + errornorm(solver.C2ex, c2)
        assert err < 0.52 * err_old
        err_old = err
