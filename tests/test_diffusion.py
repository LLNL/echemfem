import pytest
from firedrake import *
from echemfem import EchemSolver


class DiffusionSolver(EchemSolver):
    def __init__(self, N, family="DG"):
        mesh = UnitSquareMesh(N, N, quadrilateral=True)
        conc_params = []

        x, y = SpatialCoordinate(mesh)
        C1ex = sin(x) + cos(y)
        C2ex = cos(x) + sin(y)
        self.C1ex = C1ex
        self.C2ex = C2ex

        def f(C):
            f1 = - div(grad(C1ex))
            f2 = - div(grad(C2ex))
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
        self.boundary_markers = {"bulk dirichlet": (1, 2, 3, 4),
                                 }

def test_convergence():
    err_old = 1e6
    for i in range(5):
        solver = DiffusionSolver(2**(i + 1))
        solver.setup_solver()
        solver.solve()
        c1, c2 = solver.u.split()
        err = errornorm(solver.C1ex, c1) + errornorm(solver.C2ex, c2)
        assert err < 0.29 * err_old
        err_old = err

def test_convergence_CG():
    err_old = 1e6
    for i in range(5):
        solver = DiffusionSolver(2**(i + 1), family="CG")
        solver.setup_solver()
        solver.solve()
        c1, c2 = solver.u.split()
        err = errornorm(solver.C1ex, c1) + errornorm(solver.C2ex, c2)
        assert err < 0.26 * err_old
        err_old = err
