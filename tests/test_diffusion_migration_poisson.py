import pytest
from firedrake import *
from echemfem import EchemSolver


class DiffusionMigrationSolver(EchemSolver):
    def __init__(self, N, family="DG"):
        mesh = UnitSquareMesh(N, N, quadrilateral=True)
        conc_params = []

        x, y = SpatialCoordinate(mesh)
        C1ex = cos(x) + sin(y) + 3
        C2ex = cos(x) - sin(y) + 3
        Uex = sin(x) + cos(y) + 3
        self.C1ex = C1ex
        self.C2ex = C2ex
        self.Uex = Uex

        def f(C):
            D1 = 0.5
            D2 = 1.0
            z1 = 2.0
            z2 = -2.0
            K = 2e-1
            f1 = div((- D1 * z1 * grad(Uex)) * C1ex) - div(D1 * grad(C1ex))
            f2 = div((- D2 * z2 * grad(Uex)) * C2ex) - div(D2 * grad(C2ex))
            # f3 = div(( - D1 * z1**2 * grad(Uex)) * C1ex) + div(( - D2 * z2**2 * grad(Uex)) * C2ex)
            f3 = div((-K * grad(Uex))) - z1 * C1ex - z2 * C2ex
            return [f1, f2, f3]

        conc_params.append({"name": "C1",
                            "diffusion coefficient": 0.5,
                            "z": 2.,
                            "bulk": C1ex,
                            })

        conc_params.append({"name": "C2",
                            "diffusion coefficient": 1.0,
                            "z": -2.,
                            "bulk": C2ex,
                            })
        physical_params = {"flow": ["diffusion", "migration", "poisson"],
                           "F": 1.0,
                           "R": 1.0,
                           "T": 1.0,
                           "vacuum permittivity": 2.0,
                           "relative permittivity": 1e-1,
                           "U_app": Uex,
                           "bulk reaction": f,
                           }

        super().__init__(conc_params, physical_params, mesh, family=family)

    def set_boundary_markers(self):
        self.boundary_markers = {"bulk dirichlet": (1, 2, 3, 4,),
                                 "applied": (1, 2, 3, 4,),
                                 }


def test_convergence():
    errC1_old = 1e6
    errC2_old = 1e6
    errU_old = 1e6
    for i in range(5):
        solver = DiffusionMigrationSolver(2**(i + 1))
        solver.setup_solver()
        solver.solve()
        c1, c2, U = solver.u.subfunctions
        errC1 = errornorm(solver.C1ex, c1)
        errC2 = errornorm(solver.C2ex, c2)
        errU = errornorm(solver.Uex, U)
        assert errC1 < 0.26 * errC1_old
        assert errC2 < 0.26 * errC2_old
        assert errU < 0.26 * errU_old
        errC1_old = errC1
        errC2_old = errC2
        errU_old = errU


def test_convergence_CG():
    errC1_old = 1e6
    errC2_old = 1e6
    errU_old = 1e6
    for i in range(5):
        solver = DiffusionMigrationSolver(2**(i + 1), family="CG")
        solver.setup_solver()
        solver.solve()
        c1, c2, U = solver.u.subfunctions
        errC1 = errornorm(solver.C1ex, c1)
        errC2 = errornorm(solver.C2ex, c2)
        errU = errornorm(solver.Uex, U)
        assert errC1 < 0.27 * errC1_old
        assert errC2 < 0.27 * errC2_old
        assert errU < 0.27 * errU_old
        errC1_old = errC1
        errC2_old = errC2
        errU_old = errU


test_convergence_CG()
