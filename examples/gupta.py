from firedrake import *
from echemfem import EchemSolver, IntervalBoundaryLayerMesh


class GuptaSolver(EchemSolver):
    """
    A 1D example of diffusion-reaction for CO2 electrolysis with simplified
    bicarbonate bulk reactions.

    Steady-state version of example from
    Gupta, N., Gattrell, M. and MacDougall, B., 2006. Calculation for the
    cathode surface concentrations in the electrochemical reduction of CO2 in
    KHCO3 solutions. Journal of applied electrochemistry, 36(2), pp.161-172.
    """

    def __init__(self):

        delta = 0.0001
        mesh = IntervalBoundaryLayerMesh(100, delta, 100, 1e-6)

        def bulk_reaction(C):
            C_CO2 = C[0]
            C_HCO3 = C[1]
            C_CO3 = C[2]
            C_OH = C[3]

            k1f = 5.93e3  # 1/M/s
            k2f = 1e8  # 1/M/s
            k1r = 1.34e-4  # 1/s
            k2r = 2.15e4  # 1/s
            dCO2 = k1r * C_HCO3 - k1f * C_OH * C_CO2
            dHCO3 = k1f * C_OH * C_CO2 - k1r * C_HCO3 \
                - k2f * C_HCO3 * C_OH + k2r * C_CO3
            dCO3 = k2f * C_HCO3 * C_OH - k2r * C_CO3
            dOH = - k1f * C_CO2 * C_OH + k1r * C_HCO3 \
                - k2f * C_HCO3 * C_OH + k2r * C_CO3
            return [dCO2, dHCO3, dCO3, dOH]

        conc_params = []

        conc_params.append({"name": "CO2",
                            "diffusion coefficient": 19.1e-10,  # m^2/s
                            "bulk": 0.0342,  # M
                            })

        conc_params.append({"name": "HCO3",
                            "diffusion coefficient": 9.23e-10,  # m^2/s
                            "bulk": 0.499,  # M
                            })

        conc_params.append({"name": "CO3",
                            "diffusion coefficient": 11.9e-10,  # m^2/s
                            "bulk": 7.6e-4,  # M
                            })

        conc_params.append({"name": "OH",
                            "diffusion coefficient": 52.7e-10,  # m^2/s
                            "bulk": 3.3e-7,  # M
                            })

        physical_params = {"flow": ["diffusion"],
                           "bulk reaction": bulk_reaction,
                           }

        super().__init__(conc_params, physical_params, mesh, family="CG")

    def neumann(self, C, conc_params, u):
        name = conc_params["name"]
        if name in ["HCO3", "CO3"]:
            return Constant(0)
        j = 0.05  # 50 A/m^2. Adjustment for concentration units
        F = 96485.3329  # C/mol
        zeffCH4 = 8.
        zeffC2H4 = 12.
        zeffCO = 2.
        zeffHCOO = 2.
        zeffH2 = 2.
        cefCH4 = 0.25
        cefC2H4 = 0.20
        cefCO = 0.05
        cefHCOO = 0.10
        cefH2 = 0.40
        if name == "CO2":
            return -(j / F) * (cefHCOO / zeffHCOO + cefCO / zeffCO
                               + cefCH4 / zeffCH4 + 2 * cefC2H4 / zeffC2H4)
        if name == "OH":
            return (j / F) * (cefHCOO / zeffHCOO + 2 * cefCO / zeffCO
                              + 8 * cefCH4 / zeffCH4 + 12 * cefC2H4 / zeffC2H4
                              + 2 * cefH2 / zeffH2)  # note error in Gupta et al.

    def set_boundary_markers(self):
        self.boundary_markers = {"bulk dirichlet": (1,),
                                 "neumann": (2,),
                                 }


solver = GuptaSolver()
solver.setup_solver()
solver.solve()

## Plotting

C_CO2, C_HCO3, C_CO3, C_OH = solver.u.subfunctions
# OH boundary layer
x = solver.mesh.coordinates
C_OH_bl = Function(solver.V).assign(C_OH).dat.data[100:]
x_bl = x.dat.data[100:]

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
filename = "gupta.png"
fig = plt.figure(constrained_layout=True, figsize=(16, 8))
spec = gridspec.GridSpec(ncols=3, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[1, 0])
ax4 = fig.add_subplot(spec[1, 1])
ax5 = fig.add_subplot(spec[1, 2])

plot(C_CO2, axes=ax1)
ax1.set(xlabel='distance (m)',
        ylabel='CO$_2$ concentration (M)')
plot(C_HCO3, axes=ax2)
ax2.set(xlabel='distance (m)',
        ylabel='HCO$_3$ concentration (M)')
plot(C_CO3, axes=ax3)
ax3.set(xlabel='distance (m)',
        ylabel='CO$_3$ concentration (M)')
plot(C_OH, axes=ax4)
ax4.set(xlabel='distance (m)',
        ylabel='OH concentration (M)')
plt.plot(x_bl, C_OH_bl, color='k', linewidth=2)
ax5.set(xlabel='distance (m)',
        ylabel='OH concentration (M)')

plt.savefig(filename)
