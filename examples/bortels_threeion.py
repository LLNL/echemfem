from firedrake import *
from echemfem import EchemSolver
import argparse
import petsc4py
petsc4py.PETSc.Sys.popErrorHandler()
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--family", type=str, default='DG')
args, _ = parser.parse_known_args()

"""
2D Flow-plate reactor with electroneutral Nernst-Planck. Using a custom gmsh
mesh with refinement close to the electrodes. GMSH is used to get the mesh file
from the .geo file as follows:

    gmsh -2 bortels_structuredquad.geo

Three ion case taken from: 
Bortels, L., Deconinck, J. and Van Den Bossche, B., 1996. The multi-dimensional
upwinding method as a new simulation tool for the analysis of multi-ion
electrolytes controlled by diffusion, convection and migration. Part 1. Steady
state analysis of a parallel plane flow channel. Journal of Electroanalytical
Chemistry, 404(1), pp.15-26.
"""

class BortelsSolver(EchemSolver):
    def __init__(self):
        conc_params = []

        conc_params.append({"name": "Cu",
                            "diffusion coefficient": 7.2e-10,  # m^2/s
                            "z": 2.,
                            "bulk": 10.,  # mol/m^3
                            })
        conc_params.append({"name": "SO4",
                            "diffusion coefficient": 10.65e-10,  # m^2/s
                            "z": -2.,
                            "bulk": 1010.,  # mol/m^3
                            })
        conc_params.append({"name": "H",
                            "diffusion coefficient": 93.12e-10,  # m^2/s
                            "z": 1.,
                            "bulk": 2000.,  # mol/m^3
                            })

        physical_params = {"flow": ["advection", "diffusion", "migration", "electroneutrality"],
                           "F": 96485.3329,  # C/mol
                           "R": 8.3144598,  # J/K/mol
                           "T": 273.15 + 25.,  # K
                           "U_app": 0.006,  # V
                           "v_avg": 0.03  # m/s
                           }

        def reaction(u, V):
            C = u[0]
            U = u[2]
            C_b = conc_params[0]["bulk"]
            J0 = 30.  # A/m^2
            F = physical_params["F"]
            R = physical_params["R"]
            T = physical_params["T"]
            eta = V - U
            return J0 * (exp(F / R / T * (eta))
                         - (C / C_b) * exp(-F / R / T * (eta)))

        def reaction_cathode(u):
            return reaction(u, Constant(0))

        def reaction_anode(u):
            return reaction(u, Constant(physical_params["U_app"]))
        # will need to make U_app a class member if it needs updating later

        echem_params = []

        echem_params.append({"reaction": reaction_cathode,
                             "electrons": 2,
                             "stoichiometry": {"Cu": 1},
                             "boundary": "cathode",
                             })

        echem_params.append({"reaction": reaction_anode,
                             "electrons": 2,
                             "stoichiometry": {"Cu": 1},
                             "boundary": "anode",
                             })
        super().__init__(
            conc_params,
            physical_params,
            Mesh('bortels_structuredquad.msh'),
            echem_params=echem_params,
            family=args.family)

    def set_boundary_markers(self):
        # Inlet: 10, Outlet 11, Anode: 12, Cathode: 13, Wall: 14
        self.boundary_markers = {"inlet": (10,),
                                 "outlet": (11,),
                                 "anode": (12,),
                                 "cathode": (13,),
                                 }

    def set_velocity(self):
        _, y = SpatialCoordinate(self.mesh)
        h = 0.01  # m
        self.vel = as_vector(
            (6. * self.physical_params["v_avg"] / h**2 * y * (h - y), Constant(0.)))  # m/s


solver = BortelsSolver()
solver.setup_solver()
solver.solve()
