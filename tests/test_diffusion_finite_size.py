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
            num = grad(C1ex) + grad(C2ex)
            den = 1.0 - C1ex - C2ex
            f1 = - div(grad(C1ex) + C1ex * num/den)
            f2 = - div(grad(C2ex) + C2ex * num/den)
            return [f1, f2]

        conc_params.append({"name": "C1",
                            "diffusion coefficient": 1.0,
                            "bulk": C1ex,
                            "solvated diameter": 1.0  # m
                            })

        conc_params.append({"name": "C2",
                            "diffusion coefficient": 1.0,
                            "bulk": C2ex,
                            "solvated diameter": 1.0  # m
                            })
        physical_params = {"flow": ["diffusion", "finite size"],
                           "vacuum permittivity": 1.0,  # F/m
                           "relative permittivity": 1.0,
                           "Avogadro constant": 1.0,  # 1/mol
                           "bulk reaction": f,
                           }

        super().__init__(conc_params, physical_params, mesh, family=family)

    def set_boundary_markers(self):
        self.boundary_markers = {"bulk dirichlet": (1, 2, 3, 4),
                                 }


def test_DG_error():
    with pytest.raises(NotImplementedError):
        solver = DiffusionSolver(2, family="DG")


def test_convergence_CG():
    err_old = 1e6
    for i in range(5):
        solver = DiffusionSolver(2**(i + 1), family="CG")
        solver.setup_solver()
        solver.solve()
        c1, c2 = solver.u.subfunctions
        err = errornorm(solver.C1ex, c1) + errornorm(solver.C2ex, c2)
        assert err < 0.26 * err_old
        err_old = err
