"""
Created on Friday Jan 27 2023 

@author: Nitish Govindarajan

Simple example for the BMCSL model for finite size ion effects

Code to reproduce Figure 1 from doi:10.1016/j.jcis.2007.08.006
Biesheuvel, P.M. and Van Soestbergen, M., 2007. Counterion volume effects in
mixed electrical double layers. Journal of Colloid and Interface Science,
316(2), pp.490-499.
"""

# Import required libraries

from firedrake import *
from echemfem import EchemSolver, IntervalBoundaryLayerMesh
from math import log10
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import math

# Initialize bulk concentrations (1 mM CsCl + 9 mM LiCl)

C_Cs_bulk = 1

C_Li_bulk = 9

C_Cl_bulk = 10


class BMCSLSolver(EchemSolver):

    def __init__(self):

        delta = 0.00008  # Boundary layer thickness (m)

        mesh = IntervalBoundaryLayerMesh(200, delta, 800, 1e-8, boundary=(2,))

        x = SpatialCoordinate(mesh)[0]

        conc_params = []

        conc_params.append({"name": "Cl",
                            "diffusion coefficient": 20.6e-10,  # m^2/s
                            "bulk": C_Cl_bulk,  # mol/m3
                            "z": -1,
                            "solvated diameter": 0.0  # co-ion size not relevant at sufficiently high negative surface charge density
                            })

        conc_params.append({"name": "Li",
                            "diffusion coefficient": 10.6e-10,  # m^2/s
                            "bulk": C_Li_bulk,  # mol/m3
                            "z": 1,
                            "solvated diameter": 7.6e-10  # m   (a_Li = 3.8 A)
                            })

        conc_params.append({"name": "Cs",
                            "diffusion coefficient": 20.6e-10,  # m^2/s
                            "bulk": C_Cs_bulk,  # mol/m3
                            "z": 1,
                            "solvated diameter": 6.6e-10  # m   (a_Cs = 3.3 A)
                            })

        self.U_solid = Constant(0)

        physical_params = {"flow": ["migration", "poisson", "diffusion finite size_BMCSL"],  # Poisson with BMCSL finite size correction
                           "F": 96485.,  # C/mol
                           "R": 8.3144598,  # J/K/mol
                           "T": 273.15 + 25.,  # K
                           "U_app": conditional(gt(x, 1e-6), 0.0, 0.0),
                           "vacuum permittivity": 8.8541878128e-12,  # F/m
                           "relative permittivity": 78.4,
                           "Avogadro constant": 6.02214076e23,  # 1/mol
                           "surface charge density": Constant(0),  # C/m^2
                           }

        super().__init__(conc_params, physical_params, mesh, family="CG", p=2)

    def set_boundary_markers(self):
        self.boundary_markers = {"bulk dirichlet": (1,),  # C = C_0
                                 "applied": (1,),  # U_liquid = 0
                                 "poisson neumann": (2,),
                                 }


solver = BMCSLSolver()


x = SpatialCoordinate(solver.mesh)[0]

solver.setup_solver()

solver.solve()


xmesh = Function(solver.V).interpolate(solver.mesh.coordinates[0])

np.savetxt('xmesh.tsv', xmesh.dat.data)


# Surface charge densities in C/m2

sigmas = np.arange(0, -0.625, -0.025)


store_sigmas = [-0.3, -0.4, -0.5, -0.6]  # surface charge densities reported in Figure 1

for sigma in sigmas:

    sigma = round(sigma, 2)

    solver.physical_params["surface charge density"].assign(sigma)

    solver.solve()

    C_Cl, C_Li, C_Cs, Phi = solver.u.subfunctions

    if (sigma in store_sigmas):

        print(sigma)

        np.savetxt('Cs_{0:.1g}.tsv'.format(sigma), C_Cs.dat.data)

        np.savetxt('Li_{0:.1g}.tsv'.format(sigma), C_Li.dat.data)

        np.savetxt('phi_{0:.1g}.tsv'.format(sigma), Phi.dat.data)
