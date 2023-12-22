from firedrake import *
from echemfem import EchemSolver
import numpy as np
import sys
from Navier_Stokes_irregular import navier_stokes_no_slip


diffusion_coefficient=1.
peclet=10
damkohler=10
Ly = 0.1
Lx = 1.
diffusion_coefficient=Lx/peclet
mass_transfert_coefficient=damkohler/Lx*diffusion_coefficient


class CarbonateSolver(EchemSolver):
    def __init__(self):
        """
        Two reactions example reproduced from:
        Lin, T.Y., Baker, S.E., Duoss, E.B. and Beck, V.A., 2021. Analysis of
        the Reactive CO2 Surface Flux in Electrocatalytic Aqueous Flow
        Reactors. Industrial & Engineering Chemistry Research, 60(31),
        pp.11824-11833.
        """
        mesh=Mesh('squares.msh')
        
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
        self.vel=navier_stokes_no_slip(self.mesh)
        




diffusion_coefficient=Lx/peclet
mass_transfert_coefficient=damkohler/Lx*diffusion_coefficient
solver = CarbonateSolver()
solver.setup_solver()
solver.solve()

n = FacetNormal(solver.mesh)
cC1, _, = solver.u.split()
flux = assemble(dot(grad(cC1), n) * ds(11))
flux1 = assemble(dot(grad(cC1), n)/dot(grad(cC1), n)*ds(11))

Peclet=Lx/diffusion_coefficient
Damkohler=mass_transfert_coefficient*Lx/diffusion_coefficient


print("Flux = %f" % flux)
print("surface : ",flux1)
fichier=open('results_square.dat','a')
fichier.write(str(Peclet)+' '+str(damkohler)+' '+str(flux)+' '+str(flux1)+' \n')
fichier.close()

File("SquareWave_"+str(Peclet)+str(damkohler)+".pvd").write(cC1)
