import pytest
from firedrake import *
from echemfem import EchemSolver


class AdvectionDiffusionPoissonSolver(EchemSolver):
    def __init__(self, N):
        mesh = UnitSquareMesh(N, N, quadrilateral=True)
        conc_params = []

        x, y = SpatialCoordinate(mesh)
        C1ex = cos(x) + sin(y) + 3
        C2ex = cos(x) - sin(y) + 3
        Uex = sin(x) + cos(y) + 3
        U2ex = sin(x) + cos(y) + 3
        self.C1ex = C1ex
        self.C2ex = C2ex
        self.Uex = Uex
        self.U2ex = U2ex

        def f(C):
            D1 = 0.5
            D2 = 1.0
            z1 = 2.0
            z2 = -2.0
            K = 1.0
            f1 = div((self.vel) * C1ex) - \
                div(self.effective_diffusion(D1) * grad(C1ex))
            f2 = div((self.vel) * C2ex) - \
                div(self.effective_diffusion(D2) * grad(C2ex))
            f3 = div((- self.effective_diffusion(D1) * z1**2 * grad(Uex)) * C1ex) + \
                div((- self.effective_diffusion(D2) * z2**2 * grad(Uex)) * C2ex)
            f4 = div(- self.effective_diffusion(K, phase="solid") * grad(U2ex))
            return [f1, f2, f3, f4]

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
        physical_params = {
            "flow": [
                "diffusion",
                "advection",
                "poisson",
                "porous"],
            "F": 1.0,
            "R": 1.0,
            "T": 1.0,
            "U_app": Uex,
            "bulk reaction": f,
            "v_avg": 1.,
            "porosity": 0.5,
            "solid conductivity": 1.,
            "specific surface area": 1.,
            "standard potential": 1.,
        }

        super().__init__(conc_params, physical_params, mesh)

    def set_boundary_markers(self):
        self.boundary_markers = {"applied": (1, 2, 3, 4,),
                                 "liquid applied": (1, 2, 3, 4),
                                 "bulk dirichlet": (1, 2, 3, 4,),
                                 }

    def set_velocity(self):
        _, y = SpatialCoordinate(self.mesh)
        h = 1.0
        self.vel = as_vector(
            (6. * self.physical_params["v_avg"] / h**2 * y * (h - y), Constant(0.)))


def test_convergence():
    err_old = 1e6
    for i in range(5):
        solver = AdvectionDiffusionPoissonSolver(2**(i + 1))
        solver.setup_solver()
        solver.solve()
        c1, c2, U1, U2 = solver.u.subfunctions
        err = errornorm(solver.C1ex, c1) + errornorm(solver.C2ex, c2)\
            + errornorm(solver.Uex, U1) + errornorm(solver.U2ex, U2)
        assert err < 0.29 * err_old
        err_old = err
