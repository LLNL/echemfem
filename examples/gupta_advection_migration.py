from firedrake import *
from echemfem import EchemSolver, RectangleBoundaryLayerMesh


class GuptaSolver(EchemSolver):
    """
    Flow past a flat plate electrode for CO2 reduction using electroneutral
    Nernst-Planck. The homogenous bulk reactions and the constant-rate
    charge-transfer kinetics are taken from:

    Gupta, N., Gattrell, M. and MacDougall, B., 2006. Calculation for the
    cathode surface concentrations in the electrochemical reduction of CO2 in
    KHCO3 solutions. Journal of applied electrochemistry, 36(2), pp.161-172.
    """

    def __init__(self):

        Ly = 1e-3
        Lx = 5e-3

        mesh = RectangleBoundaryLayerMesh(100, 50, Lx, Ly, 50, 1e-6, boundary=(3,))
        x, y = SpatialCoordinate(mesh)
        active = conditional(And(x >= 1e-3, x < Lx-1e-3), 1., 0.)

        conc_params = []

        conc_params.append({"name": "CO2",
                            "diffusion coefficient": 19.1e-10,  # m^2/s
                            "bulk": 34.2,  # mol/m3
                            "z": 0,
                            })

        conc_params.append({"name": "HCO3",
                            "diffusion coefficient": 9.23e-10,  # m^2/s
                            "bulk": 499.,  # mol/m3
                            "z": -1,
                            })

        conc_params.append({"name": "CO3",
                            "diffusion coefficient": 11.9e-10,  # m^2/s
                            "bulk": 7.6e-1,  # mol/m3
                            "z": -2,
                            })

        conc_params.append({"name": "OH",
                            "diffusion coefficient": 52.7e-10,  # m^2/s
                            "bulk": 3.3e-4,  # mol/m3
                            "z": -1,
                            })

        conc_params.append({"name": "K",
                            "diffusion coefficient": 19.6E-10,  # m^2/s
                            "bulk": 499. + 7.6e-1 + 3.3e-4,  # mol/m3
                            "z": 1,
                            })

        homog_params = []

        homog_params.append({"stoichiometry": {"CO2": -1,
                                               "OH": -1,
                                               "HCO3": 1,
                                               },
                             "forward rate constant": 5.93,
                             "backward rate constant": 1.34e-4
                             })

        homog_params.append({"stoichiometry": {"HCO3": -1,
                                               "OH": -1,
                                               "CO3": 1,
                                               },
                             "forward rate constant": 1e5,
                             "backward rate constant": 2.15e4
                             })

        physical_params = {"flow": ["diffusion", "electroneutrality", "migration", "advection"],
                           "F": 96485.3329,  # C/mol
                           "R": 8.3144598,  # J/K/mol
                           "T": 273.15 + 25.,  # K
                           }

        def current(cef):
            j = 50.
            def curr(u):
                return cef * j * active
            return curr 

        echem_params = []

        echem_params.append({"reaction": current(0.1), # HCOO
                             "stoichiometry": {"CO2": -1,
                                               "OH": 1
                                               },
                             "electrons": 2,
                             "boundary": "electrode",
                             })

        echem_params.append({"reaction": current(0.05), # CO
                             "stoichiometry": {"CO2": -1,
                                               "OH": 2
                                               },
                             "electrons": 2,
                             "boundary": "electrode",
                             })

        echem_params.append({"reaction": current(0.25), # CH4
                             "stoichiometry": {"CO2": -1,
                                               "OH": 8
                                               },
                             "electrons": 8,
                             "boundary": "electrode",
                             })

        echem_params.append({"reaction": current(0.2), # C2H4
                             "stoichiometry": {"CO2": -2,
                                               "OH": 12
                                               },
                             "electrons": 12,
                             "boundary": "electrode",
                             })

        echem_params.append({"reaction": current(0.4), # H2
                             "stoichiometry": {"OH": 2
                                               },
                             "electrons": 2,
                             "boundary": "electrode",
                             })

        super().__init__(conc_params, physical_params, mesh, family="CG", echem_params=echem_params, homog_params=homog_params)

    def set_boundary_markers(self):
        self.boundary_markers = {"inlet": (1,),
                                 "outlet": (2,),
                                 "bulk": (4,),
                                 "electrode": (3,),
                                 }
    def set_velocity(self):
        _, y = SpatialCoordinate(self.mesh)
        self.vel = as_vector([1.91*y, Constant(0)])  # m/s

solver = GuptaSolver()
solver.setup_solver()
solver.solve()
