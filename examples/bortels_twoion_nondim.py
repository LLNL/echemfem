from firedrake import *
from echemfem import EchemSolver

"""
2D Flow-plate reactor with electroneutral Nernst-Planck. Using a custom gmsh
mesh with refinement close to the electrodes. GMSH is used to get the mesh file
from the .geo file as follows:

    gmsh -2 bortels_structuredquad_nondim.geo

Two ion case taken from: 
Bortels, L., Deconinck, J. and Van Den Bossche, B., 1996. The multi-dimensional
upwinding method as a new simulation tool for the analysis of multi-ion
electrolytes controlled by diffusion, convection and migration. Part 1. Steady
state analysis of a parallel plane flow channel. Journal of Electroanalytical
Chemistry, 404(1), pp.15-26.

Nondimensionalization from:
Roy, T., Andrej, J. and Beck, V.A., 2023. A scalable DG solver for the
electroneutral Nernst-Planck equations. Journal of Computational Physics, 475,
p.111859.
"""

class BortelsSolver(EchemSolver):
    def __init__(self):
        conc_params = []

        D_Cu = 7.2e-10  # m^2/s
        J0 = 30.  # A/m^2
        C_Cu = 10.  # mol/m^3
        D_SO4 = 10.65e-10  # m^2/s
        C_SO4 = 10.  # mol/m^3
        F = 96485.3329  # C/mol
        R = 8.3144598  # J/K/mol
        T = 273.15 + 25.  # K
        U_app = 0.1  # V
        v_avg = 0.01  # m/s

        L = 0.01  # m

        D_Cu_hat = D_Cu / L / v_avg
        D_SO4_hat = D_SO4 / L / v_avg
        J0_hat = J0 / v_avg / C_Cu / F
        U_app_hat = F / R / T * U_app

        conc_params = []
        conc_params.append({"name": "Cu",
                            "diffusion coefficient": D_Cu_hat,
                            "z": 2.,
                            "bulk": 1.,
                            })

        conc_params.append({"name": "SO4",
                            "diffusion coefficient": D_SO4_hat,
                            "z": -2.,
                            "bulk": 1.,
                            })

        physical_params = {
            "flow": [
                "advection",
                "diffusion",
                "migration",
                "electroneutrality"],
            "F": 1.0,
            "R": 1.0,
            "T": 1.0,
            "U_app": U_app_hat,
            "v_avg": 1.0}

        def reaction(u, V):
            C = u[0]
            U = u[1]
            C_b = conc_params[0]["bulk"]
            J0 = J0_hat
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

        mesh = Mesh('bortels_structuredquad_nondim.msh')
        super().__init__(conc_params, physical_params, mesh, echem_params=echem_params)

    def set_boundary_markers(self):
        # Inlet: 10, Outlet 11, Anode: 12, Cathode: 13, Wall: 14
        self.boundary_markers = {"inlet": (10,),
                                 "outlet": (11,),
                                 "anode": (12,),
                                 "cathode": (13,),
                                 }

    def set_velocity(self):
        _, y = SpatialCoordinate(self.mesh)
        h = 1.0
        self.vel = as_vector(
            (6. * self.physical_params["v_avg"] / h**2 * y * (h - y), Constant(0.)))  # m/s


solver = BortelsSolver()
solver.setup_solver()
solver.solve()
