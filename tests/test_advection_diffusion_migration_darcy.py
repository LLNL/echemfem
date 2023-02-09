import pytest
from firedrake import *
from echemfem import EchemSolver


class AdvectionDiffusionMigrationSolver(EchemSolver):
    def __init__(self, N, vavg):
        mesh = UnitSquareMesh(N, N, quadrilateral=True)
        conc_params = []

        x, y = SpatialCoordinate(mesh)
        C1ex = cos(x) + sin(y) + 3
        C2ex = C1ex
        U1ex = sin(x) + cos(y) + 3
        U2ex = sin(x) + cos(y) + 3
        pex = (0.5 * cos(x) + sin(y) + 2) * vavg
        vel = -grad(pex)
        self.C1ex = C1ex
        self.C2ex = C2ex
        self.U1ex = U1ex
        self.U2ex = U2ex
        self.pex = pex
        D1 = 0.5
        D2 = 1.0
        z1 = 2.0
        z2 = -2.0
        K = 1.0

        def f(C):
            f1 = div((vel - self.effective_diffusion(D1) * z1 * grad(U1ex)) * C1ex) - \
                div(self.effective_diffusion(D1) * grad(C1ex))
            f2 = div((vel - self.effective_diffusion(D2) * z2 * grad(U1ex)) * C2ex) - \
                div(self.effective_diffusion(D2) * grad(C2ex))
            f3 = div(- self.effective_diffusion(K, phase="solid") * grad(U2ex))
            return [f1, f2, f3]

        def fp():
            return -div(grad(pex) + self.effective_diffusion(D1) * z1 * grad(U1ex) *
                        C1ex + self.effective_diffusion(D2) * z2 * grad(U1ex) * C2ex)

        conc_params.append({"name": "C1",
                            "diffusion coefficient": 0.5,
                            "z": 2.,
                            "bulk": C1ex,
                            "molar mass": 1.0,
                            })

        conc_params.append({"name": "C2",
                            "diffusion coefficient": 1.0,
                            "z": -2.,
                            "bulk": C2ex,
                            "molar mass": 1.0,
                            })
        physical_params = {"flow": ["diffusion", "advection", "migration", "electroneutrality", "porous", "darcy"],
                           "F": 1.0,  # C/mol
                           "R": 1.0,  # J/K/mol
                           "T": 1.0,  # K
                           "U_app": U1ex,  # V
                           "bulk reaction": f,
                           "v_avg": vavg,
                           "porosity": 0.5,
                           "solid conductivity": 1.,
                           "specific surface area": 1.,
                           "standard potential": 0.,
                           "p_gas": pex,
                           "permeability": 1.0,
                           "liquid viscosity": 1.0,
                           "liquid density": 1.0,
                           "pressure source": fp,  # only for testing
                           }

        super().__init__(conc_params, physical_params, mesh)

    def set_boundary_markers(self):
        self.boundary_markers = {"inlet": (1, 2, 3, 4,),
                                 "bulk dirichlet": (1, 2, 3, 4,),
                                 "applied": (1, 2, 3, 4,),
                                 "outlet": (1, 2, 3, 4,),
                                 "liquid applied": (1, 2, 3, 4),
                                 "gas": (1, 2, 3, 4,),
                                 }


def test_convergence_low_peclet():
    err_old = 1e6
    errp_old = 1e6
    errC_old = 1e6
    errU1_old = 1e6
    errU2_old = 1e6
    for i in range(6):
        solver = AdvectionDiffusionMigrationSolver(2**(i + 1), 1.0)
        solver.setup_solver()
        solver.solve()
        c1, U1, U2, p = solver.u.subfunctions
        err = errornorm(solver.C1ex,
                        c1) + errornorm(solver.U1ex,
                                        U1) + errornorm(solver.U2ex,
                                                        U2) + errornorm(solver.pex,
                                                                        p)
        assert err < 0.29 * err_old
        err_old = err


def test_convergence_high_peclet():
    err_old = 1e6
    errp_old = 1e6
    errC_old = 1e6
    errU1_old = 1e6
    errU2_old = 1e6
    for i in range(3, 6):
        solver = AdvectionDiffusionMigrationSolver(2**(i + 1), 1e3)
        solver.setup_solver()
        solver.solve()
        c1, U1, U2, p = solver.u.subfunctions
        err = errornorm(solver.C1ex,
                        c1) + errornorm(solver.U1ex,
                                        U1) + errornorm(solver.U2ex,
                                                        U2) + errornorm(solver.pex,
                                                                        p)
        assert err < 0.35 * err_old
        err_old = err
