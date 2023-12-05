import pytest
from firedrake import *
from echemfem import EchemSolver


class DiffusionMigrationSolver(EchemSolver):
    def __init__(self, N, family="DG"):
        mesh = UnitSquareMesh(N, N, quadrilateral=True)
        conc_params = []

        x, y = SpatialCoordinate(mesh)
        C1ex = cos(x) + sin(y) + 3
        C2ex = C1ex
        Uex = sin(x) + cos(y) + 3
        self.C1ex = C1ex
        self.C2ex = C2ex
        self.Uex = Uex

        D1 = 0.5
        D2 = 1.0
        z1 = 2.0
        z2 = -2.0

        def f(C):
            f1 = div((- D1 * z1 * grad(Uex)) * C1ex) - div(D1 * grad(C1ex))
            f2 = div((- D2 * z2 * grad(Uex)) * C2ex) - div(D2 * grad(C2ex))
            return [f1, f2]

        conc_params.append({"name": "C1",
                            "diffusion coefficient": D1,
                            "z": z1,
                            "bulk": C1ex,
                            })

        conc_params.append({"name": "C2",
                            "diffusion coefficient": D2,
                            "z": z2,
                            "bulk": C2ex,
                            })
        physical_params = {
            "flow": [
                "diffusion",
                "migration",
                "electroneutrality"],
            "F": 1.0,
            "R": 1.0,
            "T": 1.0,
            "U_app": Uex,
            "bulk reaction": f,
        }

        super().__init__(conc_params, physical_params, mesh, p=1, family=family)

    def set_boundary_markers(self):
        self.boundary_markers = {"bulk dirichlet": (1, 2, 3, 4,),
                                 "applied": (1, 2, 3, 4,),
                                 }


def test_convergence():
    errC_old = 1e6
    errU_old = 1e6
    for i in range(5):
        solver = DiffusionMigrationSolver(2**(i + 1))
        solver.setup_solver()
        solver.solve()
        c1, U = solver.u.subfunctions
        errC = errornorm(solver.C1ex, c1)
        errU = errornorm(solver.Uex, U)
        assert errC < 0.26 * errC_old
        assert errU < 0.26 * errU_old
        errC_old = errC
        errU_old = errU


def test_convergence_CG():
    errC_old = 1e6
    errU_old = 1e6
    for i in range(5):
        solver = DiffusionMigrationSolver(2**(i + 1), family="CG")
        solver.setup_solver()
        solver.solve()
        c1, U = solver.u.subfunctions
        errC = errornorm(solver.C1ex, c1)
        errU = errornorm(solver.Uex, U)
        assert errC < 0.26 * errC_old
        assert errU < 0.26 * errU_old
        errC_old = errC
        errU_old = errU


test_convergence_CG()
