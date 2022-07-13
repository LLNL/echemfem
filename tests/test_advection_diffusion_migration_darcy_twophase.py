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
        X1ex = 0.2 * (sin(x) + cos(y)) + 0.5
        X2ex = 1 - X1ex
        U1ex = sin(x) + cos(y) + 3
        U2ex = sin(x) + cos(y) + 3
        pex = (0.5 * cos(x) + sin(y) + 2) * vavg
        pgex = pex
        #pgex = cos(x) + 0.5 * sin(y) + 3
        Sex = self.saturation(pex, pgex)
        krl = self.relative_permeability(Sex)
        krg = self.relative_permeability(1 - Sex)
        vel = -krl * grad(pex)
        velg = -krg * grad(pgex)
        self.C1ex = C1ex
        self.C2ex = C2ex
        self.X1ex = X1ex
        self.X2ex = X2ex
        self.U1ex = U1ex
        self.U2ex = U2ex
        self.pex = pex
        self.pgex = pgex
        D1 = 0.5
        D2 = 1.0
        z1 = 2.0
        z2 = -2.0
        K = 1.0
        Mn = 1 / (X1ex + X2ex / 2.)
        rhog = Mn * pgex
        sumyD = X2ex / 2.0
        D1M = (1. - X1ex) / sumyD / Mn
        sumyD = X1ex / 1.0
        D2M = (1. - X2ex) / sumyD / Mn
        D1K = 2. * 0.1 / 3. * sqrt(8. / pi / 1.)
        D2K = 2. * 0.1 / 3. * sqrt(8. / pi / 2.)
        DX1 = 1 / (1 / D1M + 1 / D1K)
        DX2 = 1 / (1 / D2M + 1 / D2K)

        def f(C):
            f1 = div((vel - D1 * z1 * grad(U1ex)) * C1ex) - \
                div(self.effective_diffusion(D1) * grad(C1ex))
            f2 = div((vel - D2 * z2 * grad(U1ex)) * C2ex) - \
                div(self.effective_diffusion(D2) * grad(C2ex))
            f3 = div(- self.effective_diffusion(K, phase="solid") * grad(U2ex))
            return [f1, f2, f3]

        def fg(X):
            f1 = div(velg * X1ex * rhog) - div(rhog * self.effective_diffusion(DX1,
                                                                               phase="gas") * (grad(X1ex) + X1ex * grad(Mn) / Mn))
            f2 = div(velg * X2ex * rhog) - div(rhog * self.effective_diffusion(DX2,
                                                                               phase="gas") * (grad(X2ex) + X2ex * grad(Mn) / Mn))
            return [f1, f2]

        def fp():
            return -div(krl * grad(pex) + D1 * z1 * grad(U1ex)
                        * C1ex + D2 * z2 * grad(U1ex) * C2ex)

        # def fpg():
        #    return -div(krg * rhog * grad(pgex))

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
                           "porosity": 0.5,
                           "solid conductivity": 1.,
                           "specific surface area": 1.,
                           "standard potential": 0.,
                           "p_gas": pex,
                           "permeability": 1.0,
                           "liquid viscosity": 1.0,
                           "liquid density": 1.0,
                           "pressure source": fp,  # only for testing
                           # "gas pressure source": fpg, #only for testing
                           "gas source": fg,  # only for testing
                           "gas viscosity": 1.0,
                           "average pore radius": 0.1,
                           }
        gas_params = []
        gas_params.append({"name": "X1",
                           "diffusion coefficient": {"X2": 1.0},
                           "gas": X1ex,
                           "molar mass": 1.0,
                           })

        gas_params.append({"name": "X2",
                           "diffusion coefficient": {"X1": 1.0},
                           "gas": X2ex,
                           "n": 2,
                           "molar mass": 2.0,
                           })

        super().__init__(conc_params, physical_params, mesh, gas_params=gas_params)

    def set_boundary_markers(self):
        self.boundary_markers = {"inlet": (1, 2, 3, 4,),
                                 "applied": (1, 2, 3, 4,),
                                 "outlet": (1, 2, 3, 4,),
                                 "bulk dirichlet": (1, 2, 3, 4),
                                 "liquid applied": (1, 2, 3, 4),
                                 "gas": (1, 2, 3, 4),
                                 "gas inlet": (1, 2, 3, 4),
                                 "gas outlet": (1, 2, 3, 4),
                                 }

    def saturation(self, pl, pg):
        # eye-ball approximation of the saturation curve from Weng Weber
        S_wr = 0.25  # Residual saturation of water
        S_gr = 0.01  # Residual saturation of gas
        c = 5e-2
        p0 = 0
        pc = pl - pg  # capillary pressure
        return ((1. - S_gr) - S_wr) * tanh(c * (pc - p0)) / 2. + 0.64
        # return 0.64

    def set_gas_density(self, X, pg, gas_params):
        R = self.physical_params["R"]
        T = self.physical_params["T"]

        # average molar mass of the mixture
        Mn = 0.
        Xj = 1.0
        for i in range(self.num_g):
            if gas_params[i].get("eliminated") is None:
                Mn += X[i] / gas_params[i]["molar mass"]
                Xj -= X[i]
            else:
                j = i
        Mn += Xj / gas_params[j]["molar mass"]
        Mn = 1. / Mn
        self.Mn = Mn
        self.Xj = Xj

        # ideal gas law
        self.rhog = Mn * pg / R / T
        if self.flow["diffusion"] and (
                self.boundary_markers.get("gas") is not None):
            Mn_gas = 0.0
            for i in range(self.num_g):
                Mn_gas += gas_params[i]["gas"] / gas_params[i]["molar mass"]
            Mn_gas = 1. / Mn_gas
            self.Mn_gas = Mn_gas


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
        c1, X1, U1, U2, p, pg = solver.u.split()
        err = errornorm(solver.C1ex,
                        c1) + errornorm(solver.U1ex,
                                        U1) + errornorm(solver.U2ex,
                                                        U2) + errornorm(solver.pex,
                                                                        p) + errornorm(solver.X1ex,
                                                                                       X1) + errornorm(solver.pgex,
                                                                                                       pg)
        assert err < 0.29 * err_old
        err_old = err


def test_convergence_high_peclet():
    err_old = 1e6
    errp_old = 1e6
    errC_old = 1e6
    errU1_old = 1e6
    errU2_old = 1e6
    for i in range(4, 6):
        solver = AdvectionDiffusionMigrationSolver(2**(i + 1), 1e2)
        solver.setup_solver()
        solver.solve()
        c1, X1, U1, U2, p, pg = solver.u.split()
        err = errornorm(solver.C1ex,
                        c1) + errornorm(solver.U1ex,
                                        U1) + errornorm(solver.U2ex,
                                                        U2) + errornorm(solver.pex,
                                                                        p) + errornorm(solver.X1ex,
                                                                                       X1) + errornorm(solver.pgex,
                                                                                                       pg)
        assert err < 0.29 * err_old
        err_old = err
