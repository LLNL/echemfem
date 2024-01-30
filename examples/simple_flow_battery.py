from firedrake import *
from echemfem import EchemSolver, NavierStokesBrinkmanFlowSolver, NavierStokesFlowSolver
import argparse
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--family", type=str, default='CG')
parser.add_argument("--vel_file", type=str, default=None)
args, _ = parser.parse_known_args()
if args.family == "CG":
    SUPG = True
else:
    SUPG = False

#flow_rate = 2.4 * 5/3 * 1e-8 # m3/s -> x ml/min
#area = 4e-4 # m2 -> 4 cm2
#v_avg = flow_rate/area
v_avg = 1e-10#
print(v_avg)

current_density = 50 * 10 # A/m2 -> 50 mA/cm2
current_density = 1e-1 # A/m2 -> 50 mA/cm2

electrode_thickness = 500e-6 # m
electrode_length = 500e-6#2e-3 # m
inlet_length = electrode_length/4
outlet_length = electrode_length/4
mesh = RectangleMesh(50, 50, electrode_length, electrode_thickness, quadrilateral=True)
if False:
    x, y = SpatialCoordinate(mesh)
    V1 = FunctionSpace(mesh, "HDiv Trace", 0)
    f3 = Function(V1).interpolate(conditional(And(And(y < 1e-16,
                                                    x <= electrode_length - outlet_length),
                                                    x >= inlet_length), 1., 0.))
    f5 = Function(V1).interpolate(conditional(And(y < 1e-16, x < inlet_length), 1., 0.))
    f6 = Function(V1).interpolate(conditional(And(y < 1e-16, x > electrode_length - outlet_length), 1., 0.))
    mesh = RelabeledMesh(mesh,
                         [f3, f5, f6],
                         [3, 5, 6])

"""
Schematic of the domain:
     _____4_____
    |           |
    1           2
    |           |
    |_5_._3_._6_|

"""

class BortelsSolver(EchemSolver):
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

        #physical_params = {"flow": ["advection", "diffusion",  "porous"],
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
        #def reaction(u):
        #    return 1e-1



        echem_params = []
        echem_params.append({"reaction": reaction,
                             "electrons": 1,
                             "stoichiometry": {"V2": 1,
                                               "V3": -1},
                             })
        #echem_params = []
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
        boundary_markers = {"no slip": (3,4,1,2,),
                            #"inlet velocity": (5,),
                            #"outlet velocity": (6,),
                            "outlet pressure": (6,),
                            "inlet pressure": (5,),
                            }
        boundary_markers = {"no slip": (3,4,),
                            #"inlet velocity": (5,),
                            #"outlet velocity": (6,),
                            "outlet pressure": (2,),
                            "inlet pressure": (1,),
                            }

        #vel = as_vector([Constant(0), Constant(v_avg)])
        #vel_out = as_vector([Constant(0), Constant(-v_avg)])
        x, y = SpatialCoordinate(self.mesh)
        vel = as_vector([Constant(0), 4 * Constant(v_avg) * x * (Constant(inlet_length) - x)/Constant(inlet_length**2)])
        vel_out = as_vector([Constant(0), 4 * Constant(-v_avg) * (x - Constant(electrode_length-outlet_length)) * (Constant(electrode_length) - x)/Constant(outlet_length**2)])
        flow_params = {"inlet velocity": vel,
                       "outlet pressure": 0.,
                       "inlet pressure": 4e-3,
                       "outlet velocity": vel_out,
                       "density": 1e3, # kg/m3
                       "dynamic viscosity": 8.9e-4, # Pa s
                       "permeability": 5.53e-11 # m2
                       }
        NS_solver = NavierStokesBrinkmanFlowSolver(mesh, flow_params, boundary_markers)
        #NS_solver = NavierStokesFlowSolver(mesh, flow_params, boundary_markers)
        NS_solver.setup_solver()
        NS_solver.solve()
        self.vel = NS_solver.vel


solver = BortelsSolver()
solver.setup_solver(initial_solve=False)
solver.solve()
us = solver.u.subfunctions
i_n = Function(solver.V, name="Current Density").interpolate(solver.echem_params[0]["reaction"](solver.u))
File("results/current_density.pvd").write(i_n)
Vec = VectorFunctionSpace(solver.mesh, "CG", 1)
kappa_0 = solver.physical_params["liquid conductivity"]
kappa = solver.physical_params["porosity"]**1.5 * kappa_0
i_l = Function(Vec, name="ionic current").interpolate(kappa * grad(us[solver.i_Ul]))
File("results/current.pvd").write(i_l)
av = solver.physical_params["specific surface area"]
print(assemble(av * i_n * dx)/assemble(Constant(1) * dx(domain=solver.mesh)))
print(current_density/electrode_length)
