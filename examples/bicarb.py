from firedrake import *
from echemfem import EchemSolver


class CarbonateSolver(EchemSolver):
    def __init__(self):
        """
        Flow past an electrode with bicarbonate bulk reactions and linear CO2 electrolysis
        Example reproduced from:
        Lin, T.Y., Baker, S.E., Duoss, E.B. and Beck, V.A., 2021. Analysis of
        the Reactive CO2 Surface Flux in Electrocatalytic Aqueous Flow
        Reactors. Industrial & Engineering Chemistry Research, 60(31),
        pp.11824-11833.
        """

        Ly = 1e-3  # .6e-3 # m
        Lx = 1e-2  # m

        mesh = RectangleMesh(200, 100, Lx, Ly, quadrilateral=True)

        k1f = (8.42E3)/1000
        k1r = 1.97E-4
        k2f = 3.71E-2
        k2r = (8.697E4)/1000
        k3f = 1.254E6
        k3r = (6E9)/1000
        k4f = (1.242E12)/1000
        k4r = 59.44
        k5f = (2.3E10)/1000
        k5r = (2.3E-4)*1000

        C_1_inf = 0.034*1000
        C_K = 1.*1000
        Keq = k1f*k3f/(k1r*k3r)

        C_3_inf = -C_1_inf*Keq/4 + C_1_inf*Keq/4*sqrt(1+8*C_K/(C_1_inf*Keq))
        C_4_inf = (C_K - C_3_inf)/2
        C_2_inf = C_3_inf/C_1_inf*k1r/k1f
        C_5_inf = (k5r/k5f)/C_2_inf

        def bulk_reaction(y):
            yCO2 = y[0]
            yOH = y[1]
            yHCO3 = y[2]
            yCO3 = y[3]
            yH = y[4]

            dCO2 = -(k1f) * yCO2*yOH \
                + (k1r) * yHCO3 \
                - (k2f) * yCO2 \
                + (k2r) * yHCO3*yH

            dOH = -(k1f) * yCO2*yOH \
                + (k1r) * yHCO3 \
                + (k3f) * yCO3 \
                - (k3r) * yOH*yHCO3\
                - (k5f) * yOH*yH\
                + (k5r)

            dHCO3 = (k1f) * yCO2*yOH\
                - (k1r) * yHCO3\
                + (k3f) * yCO3\
                - (k3r) * yOH*yHCO3\
                + (k2f) * yCO2 \
                - (k2r) * yHCO3*yH \
                + (k4f) * yCO3*yH\
                - (k4r) * yHCO3

            dCO3 = -(k3f) * yCO3 \
                + (k3r) * yOH*yHCO3\
                - (k4f)*yCO3*yH\
                + (k4r) * yHCO3

            dH = (k2f) * yCO2 \
                - (k2r) * yHCO3*yH \
                - (k4f) * yCO3*yH\
                + (k4r) * yHCO3\
                - (k5f) * yOH*yH\
                + (k5r)

            return [dCO2, dOH, dHCO3, dCO3, dH]  # , 0.]

        C_CO2_bulk = C_1_inf
        C_OH_bulk = C_2_inf
        C_HCO3_bulk = C_3_inf
        C_CO32_bulk = C_4_inf
        C_H_bulk = C_5_inf
        C_K_bulk = C_K

        conc_params = []

        conc_params.append({"name": "CO2",
                            "diffusion coefficient": 1.91E-9,  # m^2/s
                            "bulk": C_CO2_bulk,  # mol/m3
                            })

        conc_params.append({"name": "OH",
                            "diffusion coefficient": 5.29E-9,  # m^2/s
                            "bulk": C_OH_bulk,  # mol/m3
                            })

        conc_params.append({"name": "HCO3",
                            "diffusion coefficient": 1.185E-9,  # m^2/s
                            "bulk": C_HCO3_bulk,  # mol/m3
                            })

        conc_params.append({"name": "CO3",
                            "diffusion coefficient": .92E-9,  # m^2/s
                            "bulk": C_CO32_bulk,  # mol/m3
                            })

        conc_params.append({"name": "H",
                            "diffusion coefficient": 9.311E-9,  # m^2/s
                            "bulk": C_H_bulk,  # mol/m3
                            })

        # conc_params.append({"name": "K",
        #                    "diffusion coefficient": 1.96E-9,  # m^2/s
        #                    "bulk": C_K_bulk,  # mol/m3
        #                    })

        physical_params = {"flow": ["advection", "diffusion"],
                           "bulk reaction": bulk_reaction,
                           }

        super().__init__(conc_params, physical_params, mesh)

    def neumann(self, C, conc_params, u):
        name = conc_params["name"]
        if name in ["HCO3", "CO3", "H"]:
            return Constant(0)

        if name == "CO2":
            return -(1.91E-1)*u[0]
        if name == "OH":
            return 2*(1.91E-1)*u[0]

    def set_boundary_markers(self):
        self.boundary_markers = {"inlet": (1),
                                 "bulk dirichlet": (4),
                                 "outlet": (2,),
                                 "neumann": (3,),
                                 }

    def set_velocity(self):
        _, y = SpatialCoordinate(self.mesh)
        self.vel = as_vector([1.91*y, Constant(0)])  # m/s


solver = CarbonateSolver()
solver.setup_solver()
solver.solve()

n = FacetNormal(solver.mesh)
cCO2, _, _, _, _ = solver.u.subfunctions
flux = assemble(dot(grad(cCO2), n) * ds(3))
print("flux = %f" % flux)
