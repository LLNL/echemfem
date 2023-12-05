import pytest
from firedrake import *
from echemfem import EchemSolver


class AdvectionDiffusionMigrationSolver(EchemSolver):
    def __init__(self, N, vavg, family="DG"):
        mesh = UnitSquareMesh(N, N, quadrilateral=True)
        conc_params = []

        x, y = SpatialCoordinate(mesh)
        C1ex = cos(x) + sin(y) + 3
        C2ex = C1ex
        Uex = sin(x) + cos(y) + 3
        self.C1ex = C1ex
        self.C2ex = C2ex
        self.Uex = Uex

        def f(C):
            D1 = 0.5
            D2 = 1.0
            z1 = 2.0
            z2 = -2.0
            f1 = div((self.vel - D1 * z1 * grad(Uex))
                     * C1ex) - div(D1 * grad(C1ex))
            f2 = div((self.vel - D2 * z2 * grad(Uex))
                     * C2ex) - div(D2 * grad(C2ex))
            return [f1, f2]

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
        physical_params = {"flow": ["diffusion", "advection", "migration", "electroneutrality"],
                           "F": 1.0,  # C/mol
                           "R": 1.0,  # J/K/mol
                           "T": 1.0,  # K
                           "U_app": Uex,  # V
                           "bulk reaction": f,
                           "v_avg": Constant(vavg),
                           }

        super().__init__(conc_params, physical_params, mesh, family=family)

    def set_boundary_markers(self):
        self.boundary_markers = {  # "inlet": (1,2,3,4,),
            "applied": (1, 2, 3, 4,),
            # "outlet": (1,2,3,4,),
            "bulk dirichlet": (1, 2, 3, 4,),
        }

    def set_velocity(self):
        _, y = SpatialCoordinate(self.mesh)
        h = 1.0
        self.vel = as_vector(
            (6. * self.physical_params["v_avg"] / h**2 * y * (h - y), Constant(0.)))


def test_convergence_low_peclet():
    errC_old = 1e6
    errU_old = 1e6
    for i in range(6):
        solver = AdvectionDiffusionMigrationSolver(2**(i + 1), 1.0)
        solver.setup_solver()
        solver.solve()
        c1, U = solver.u.subfunctions
        errC = errornorm(solver.C1ex, c1)
        errU = errornorm(solver.Uex, U)
        assert errC < 0.26 * errC_old
        assert errU < 0.26 * errU_old
        errC_old = errC
        errU_old = errU


def test_convergence_high_peclet():
    errC_old = 1e6
    errU_old = 1e6
    for i in range(6):
        solver = AdvectionDiffusionMigrationSolver(2**(i + 1), 1e5)
        solver.setup_solver()
        solver.solve()
        c1, U = solver.u.subfunctions
        errC = errornorm(solver.C1ex, c1)
        errU = errornorm(solver.Uex, U)
        assert errC < 0.36 * errC_old  # p+1/2 convergence
        assert errU < 0.26 * errU_old
        errC_old = errC
        errU_old = errU


def test_convergence_low_peclet_CG():
    errC_old = 1e6
    errU_old = 1e6
    for i in range(6):
        solver = AdvectionDiffusionMigrationSolver(2**(i + 1), 1.0, family="CG")
        solver.setup_solver()
        solver.solve()
        c1, U = solver.u.subfunctions
        errC = errornorm(solver.C1ex, c1)
        errU = errornorm(solver.Uex, U)
        assert errC < 0.26 * errC_old
        assert errU < 0.26 * errU_old
        errC_old = errC
        errU_old = errU


def test_convergence_high_peclet_CG():
    errC_old = 1e6
    errU_old = 1e6
    for i in range(6):
        solver = AdvectionDiffusionMigrationSolver(2**(i + 2), 1e5, family="CG")
        solver.setup_solver()
        solver.solve()
        c1, U = solver.u.subfunctions
        errC = errornorm(solver.C1ex, c1)
        errU = errornorm(solver.Uex, U)
        assert errC < 0.66 * errC_old
        assert errU < 0.51 * errU_old
        errC_old = errC
        errU_old = errU
