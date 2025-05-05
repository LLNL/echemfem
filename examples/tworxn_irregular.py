from firedrake import *
from echemfem import EchemSolver, NavierStokesFlowSolver

"""
A 2D flow past an irregular electrode toy model with two species and
advection-diffusion, and Navier-Stokes for the flow. The electrode surface
consists of square blocks. The mesh can be create in GMSH using the following
command:

    gmsh -2 squares_small.geo

The model is adapated from
Lin, T.Y., Baker, S.E., Duoss, E.B. and Beck, V.A., 2021. Analysis of
the Reactive CO2 Surface Flux in Electrocatalytic Aqueous Flow
Reactors. Industrial & Engineering Chemistry Research, 60(31),
pp.11824-11833.
"""

peclet=10
damkohler=10
Ly = 0.1
Lx = 1.
diffusion_coefficient=Lx/peclet
mass_transfert_coefficient=damkohler/Lx*diffusion_coefficient
mesh=Mesh('squares_small.msh')

class CarbonateSolver(EchemSolver):
    def __init__(self):
        
        C_1_inf = 1.
        C_2_inf = Constant(0)

       
        def bulk_reaction(y):
            yC1=y[0];
            yC2=y[1];
            dC1 = -(1.)*(1e3)*yC1*yC2 
            dC2 = -(2.)*(1e3)*yC1*yC2 
            return [0,0]

        conc_params = []
        conc_params.append({"name": "C1",
                            "diffusion coefficient": diffusion_coefficient,
                            "bulk": C_1_inf,
                            })

        conc_params.append({"name": "C2",
                            "diffusion coefficient": diffusion_coefficient,
                            "bulk": C_2_inf,
                            })


        physical_params = {"flow": ["advection","diffusion"],
                           "bulk reaction": bulk_reaction,
                           }

        super().__init__(conc_params, physical_params, mesh, family="DG")

    def neumann(self, C, conc_params, u):
        name = conc_params["name"]

        if name == "C1":
            return -(mass_transfert_coefficient)*u[0]
        if name == "C2":
            return 2.*(1.e6)*u[0]

    def set_boundary_markers(self):
        self.boundary_markers = {"inlet": (12),
                                 "bulk dirichlet": (13),
                                 "outlet": (14,),
                                 "neumann": (11,), 
                                 }

    def set_velocity(self):
        boundary_markers = {"no slip": (11,10,15),
                            "inlet velocity": (12,13,),
                            "outlet velocity": (14,)
                            }

        x, y = SpatialCoordinate(mesh)
        vel = as_vector([y, Constant(0)])
        flow_params = {"inlet velocity": vel,
                       "outlet velocity": vel,
                       "Reynolds number": 100
                       }
        NS_solver = NavierStokesFlowSolver(mesh, flow_params, boundary_markers)
        NS_solver.setup_solver()
        NS_solver.solve()
        self.vel = NS_solver.vel
        

solver = CarbonateSolver()
solver.setup_solver()
solver.solve()

n = FacetNormal(solver.mesh)
cC1, _, = solver.u.subfunctions
flux = assemble(dot(grad(cC1), n) * ds(11))
flux1 = assemble(dot(grad(cC1), n)/dot(grad(cC1), n)*ds(11))


print("Flux = %f" % flux)
print("surface : ",flux1)
fichier=open('results_square.dat','a')
fichier.write(str(peclet)+' '+str(damkohler)+' '+str(flux)+' '+str(flux1)+' \n')
fichier.close()

VTKFile("SquareWave_"+str(peclet)+str(damkohler)+".pvd").write(cC1)
