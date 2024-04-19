from firedrake import *
from echemfem import EchemSolver, RectangleBoundaryLayerMesh

"""
A 2D flow past the electrode toy model with two species and
advection-diffusion. Taken from
Lin, T.Y., Baker, S.E., Duoss, E.B. and Beck, V.A., 2021. Analysis of
the Reactive CO2 Surface Flux in Electrocatalytic Aqueous Flow
Reactors. Industrial & Engineering Chemistry Research, 60(31),
pp.11824-11833.
"""

class CarbonateSolver(EchemSolver):
    def __init__(self):

        Ly = 0.1
        Lx = 1.

        mesh = RectangleBoundaryLayerMesh(50, 50, Lx, Ly, 50, 1e-1, Ly_bdlayer=5e-3, boundary=(3, 1,))

        C_1_inf = 1.
        C_2_inf = Constant(0)

        def bulk_reaction(y):
            yC1 = y[0]
            yC2 = y[1]
            dC1 = -(1.)*(1e3)*yC1*yC2
            dC2 = -(2.)*(1e3)*yC1*yC2
            return [dC1, dC2]

        conc_params = []
        conc_params.append({"name": "C1",
                            "diffusion coefficient": 1.,
                            "bulk": C_1_inf,
                            })

        conc_params.append({"name": "C2",
                            "diffusion coefficient": 1.,
                            "bulk": C_2_inf,
                            })

        physical_params = {"flow": ["advection", "diffusion"],
                           "bulk reaction": bulk_reaction,
                           }

        super().__init__(conc_params, physical_params, mesh, family="DG")

    def neumann(self, C, conc_params, u):
        name = conc_params["name"]

        if name == "C1":
            return -(1.e6)*u[0]
        if name == "C2":
            return 2.*(1.e6)*u[0]

    def set_boundary_markers(self):
        self.boundary_markers = {"inlet": (1),
                                 "bulk dirichlet": (4),
                                 "outlet": (2,),
                                 "neumann": (3,),
                                 }

    def set_velocity(self):
        _, y = SpatialCoordinate(self.mesh)
        self.vel = as_vector([(1.e5)*y, Constant(0)])


solver = CarbonateSolver()
solver.setup_solver()
solver.solve()

n = FacetNormal(solver.mesh)
cC1, _, = solver.u.subfunctions
flux = assemble(dot(grad(cC1), n) * ds(3))
print("Sh = %f" % flux)
