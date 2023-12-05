import pytest
from firedrake import *
from echemfem import EchemSolver


class AdvectionSolver(EchemSolver):
    def __init__(self, N, extruded=False):
        if extruded:
            plane_mesh = UnitSquareMesh(N, N, quadrilateral=True)
            mesh = ExtrudedMesh(plane_mesh, N, layer_height=1.0 / N)
            x, y, _ = SpatialCoordinate(mesh)
        else:
            mesh = UnitSquareMesh(N, N, quadrilateral=True)
            x, y = SpatialCoordinate(mesh)

        conc_params = []

        C1ex = sin(x) + cos(y)
        C2ex = cos(x) + sin(y)
        self.C1ex = C1ex
        self.C2ex = C2ex

        def f(C):
            f1 = div(self.vel * C1ex)
            f2 = div(self.vel * C2ex)
            return [f1, f2]

        conc_params.append({"name": "C1",
                            "bulk": C1ex,
                            })

        conc_params.append({"name": "C2",
                            "bulk": C2ex,
                            })
        physical_params = {"flow": ["advection"],
                           "bulk reaction": f,
                           "v_avg": 1.,
                           }

        super().__init__(conc_params, physical_params, mesh)

    def set_boundary_markers(self):
        self.boundary_markers = {"inlet": (1, 2, 3, 4),
                                 "outlet": (1, 2, 3, 4),
                                 }

    def set_velocity(self):
        if self.mesh.layers is not None:
            _, y, _ = SpatialCoordinate(self.mesh)
        else:
            _, y = SpatialCoordinate(self.mesh)

        h = 1.0
        x_vel = 6. * self.physical_params["v_avg"] / h**2 * y * (h - y)

        if self.mesh.layers is not None:
            self.vel = as_vector(
                (x_vel, Constant(0.), Constant(0.)))
        else:
            self.vel = as_vector(
                (x_vel, Constant(0.)))


def test_convergence(extruded=False):
    err_old = 1e6
    if extruded:
        n = 3
    else:
        n = 5
    for i in range(n):
        solver = AdvectionSolver(2**(i + 1), extruded)
        solver.setup_solver()
        solver.solve()
        c1, c2 = solver.u.subfunctions
        err = errornorm(solver.C1ex, c1) + errornorm(solver.C2ex, c2)
        assert err < 0.26 * err_old
        err_old = err


def test_convergence_extruded():
    test_convergence(extruded=True)
