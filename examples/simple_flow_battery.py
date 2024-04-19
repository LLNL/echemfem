from firedrake import *
from echemfem import EchemSolver, NavierStokesBrinkmanFlowSolver
import argparse


"""
A simple Vanadium flow battery using advection-diffusion-reaction, Poisson for
the ionic potential with a predefinied conductivity and Navier-Stokes-Brinkman
for the flow.

Model taken from
Lin, T.Y., Baker, S.E., Duoss, E.B. and Beck, V.A., 2022. Topology optimization
of 3D flow fields for flow batteries. Journal of The Electrochemical Society,
169(5), p.050540.
"""

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--family", type=str, default='CG')
parser.add_argument("--vel_file", type=str, default=None)
args, _ = parser.parse_known_args()
if args.family == "CG":
    SUPG = True
else:
    SUPG = False

current_density = 50 * 10 # A/m2 -> 50 mA/cm2

electrode_thickness = 500e-6 # m
electrode_length = 500e-6 # m
mesh = RectangleMesh(50, 50, electrode_length, electrode_thickness, quadrilateral=True)

class FlowBatterySolver(EchemSolver):
    def __init__(self):
        conc_params = []

        conc_params.append({"name": "V2",
                            "diffusion coefficient": 2.4e-10,  # m^2/s
                            "bulk": 1000.,  # mol/m^3
                            })

        conc_params.append({"name": "V3",
                            "diffusion coefficient": 2.4e-10,  # m^2/s
                            "bulk": 1000.,  # mol/m^3
                            })

        physical_params = {"flow": ["advection", "diffusion", "poisson", "porous"],
                           "F": 96485.3329,  # C/mol
                           "R": 8.3144598,  # J/K/mol
                           "T": 273.15 + 25.,  # K
                           "solid conductivity": 1e4, # S/m
                           "liquid conductivity": 40, # S/m
                           "specific surface area": 8e4, # 1/m 
                           "porosity": 0.68,
                           "U_app": 0., # ground U_solid = 0
                           "applied current density": Constant(current_density), 
                           }

        def reaction(u):
            # Butler-Volmer
            V2 = u[0]
            V3 = u[1]
            Phi2 = u[2]
            Phi1 = u[3]
            #Phi2 = -0.2499
            #Phi1 = 0
            Cref = 1. # mol/m3
            J0 = 0.016  # A/m^2
            U0 = -0.25
            F = physical_params["F"]
            R = physical_params["R"]
            T = physical_params["T"]
            beta = 0.5 * F / R / T
            eta = Phi1 - Phi2 - U0
            return J0 / Cref * (V2 * exp(beta * eta)
                                - V3 * exp(-beta * eta))

        echem_params = []
        echem_params.append({"reaction": reaction,
                             "electrons": 1,
                             "stoichiometry": {"V2": -1,
                                               "V3": 1},
                             })
        super().__init__(
            conc_params,
            physical_params,
            mesh,
            echem_params=echem_params,
            family=args.family,
            SUPG=SUPG)

    def set_boundary_markers(self):
        self.boundary_markers = {"inlet": (1,),
                                 "outlet": (2,),
                                 "applied": (3,), # U_solid = 0
                                 "applied liquid current": (4,),
                                 }

    def set_velocity(self):
        boundary_markers = {"no slip": (3,4,),
                            "inlet pressure": (1,),
                            "outlet pressure": (2,),
                            }

        flow_params = {"outlet pressure": 0.,
                       "inlet pressure": 1e-1,
                       "density": 1e3, # kg/m3
                       "dynamic viscosity": 8.9e-4, # Pa s
                       "permeability": 5.53e-11 # m2
                       }

        NSB_solver = NavierStokesBrinkmanFlowSolver(mesh, flow_params, boundary_markers)
        NSB_solver.setup_solver()
        NSB_solver.solve()
        self.vel = NSB_solver.vel


solver = FlowBatterySolver()
solver.setup_solver(initial_solve=False)
solver.solve()
