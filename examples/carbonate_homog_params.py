from firedrake import *
from echemfem import EchemSolver, IntervalBoundaryLayerMesh

"""
A 1D example of diffusion-reaction for CO2 electrolysis with bicarbonate bulk
reactions. The charge-transfer reactions are implemented by hand but the bulk
reactions are implemented using the homog_params interface.

Steady-state version of example from
Gupta, N., Gattrell, M. and MacDougall, B., 2006. Calculation for the
cathode surface concentrations in the electrochemical reduction of CO2 in
KHCO3 solutions. Journal of applied electrochemistry, 36(2), pp.161-172.

Using the bicarbonate bulk reactions from
Schulz, K.G., Riebesell, U., Rost, B., Thoms, S. and Zeebe, R.E., 2006.
Determination of the rate constants for the carbon dioxide to bicarbonate
inter-conversion in pH-buffered seawater systems. Marine chemistry,
100(1-2), pp.53-65.
"""

class CarbonateSolver(EchemSolver):
    def __init__(self):

        delta = 0.0001
        mesh = IntervalBoundaryLayerMesh(100, delta, 100, 1e-7)

        # Reaction Rates in SI units
        k1f = 2.23  # m3/mol.s
        k1r = 9.71e-5  # 1/s
        k2f = 3.71e-2  # 1/s
        k2r = 2.67e1  # m3/mol.s
        k3f = 3.06e5  # 1/s
        k3r = 6.0e6  # m3/mol.s
        k4f = 5.0e7  # m3/mol.s
        k4r = 59.44  # 1/s
        k5f = 2.31e7  # m3/mols. (Note - typo in Schulz)
        k5r = 1.4  # mol/m3.s
        # Bulk Conditions
        C_HCO3_bulk = 500.  # mol/m3
        C_CO32_bulk = 6.5  # mol/m3
        C_OH_bulk = k3f / k3r * C_CO32_bulk / C_HCO3_bulk  # mol/m3
        C_H_bulk = (k5r / k5f) / C_OH_bulk  # mol/m3
        C_CO2_bulk = (k1r / k1f) * C_HCO3_bulk / C_OH_bulk  # mol/m3

        homog_params = []

        homog_params.append({"stoichiometry": {"CO2": -1,
                                               "OH": -1,
                                               "HCO3": 1},
                             "forward rate constant": k1f,
                             "backward rate constant": k1r
                             })

        homog_params.append({"stoichiometry": {"CO2": -1,
                                               "H": 1,
                                               "HCO3": 1},
                             "forward rate constant": k2f,
                             "backward rate constant": k2r
                             })

        homog_params.append({"stoichiometry": {"CO3": -1,
                                               "OH": 1,
                                               "HCO3": 1},
                             "forward rate constant": k3f,
                             "backward rate constant": k3r
                             })

        homog_params.append({"stoichiometry": {"CO3": -1,
                                               "H": -1,
                                               "HCO3": 1},
                             "forward rate constant": k4f,
                             "backward rate constant": k4r
                             })

        homog_params.append({"stoichiometry": {"OH": -1,
                                               "H": -1},
                             "forward rate constant": k5f,
                             "backward rate constant": k5r
                             })

        conc_params = []

        conc_params.append({"name": "CO2",
                            "diffusion coefficient": 19.1e-10,  # m^2/s
                            "bulk": C_CO2_bulk,  # mol/m3
                            })

        conc_params.append({"name": "OH",
                            "diffusion coefficient": 54.0e-10,  # m^2/s
                            "bulk": C_OH_bulk,  # mol/m3
                            })

        conc_params.append({"name": "H",
                            "diffusion coefficient": 96.0e-10,  # m^2/s
                            "bulk": C_H_bulk,  # mol/m3
                            })

        conc_params.append({"name": "CO3",
                            "diffusion coefficient": 7.0e-10,  # m^2/s
                            "bulk": C_CO32_bulk,  # mol/m3
                            })

        conc_params.append({"name": "HCO3",
                            "diffusion coefficient": 9.4e-10,  # m^2/s
                            "bulk": C_HCO3_bulk,  # mol/m3
                            })

        physical_params = {"flow": ["diffusion"],
                           "F": 96485.,  # C/mol
                           }

        super().__init__(conc_params, physical_params, mesh, family="CG", homog_params=homog_params)

    def neumann(self, C, conc_params, u):
        name = conc_params["name"]
        if name in ["HCO3", "CO3", "H"]:
            return Constant(0)
        j = 50.  # A/m^2
        F = self.physical_params["F"]
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
            print(-(j / F) * (cefHCOO / zeffHCOO + cefCO / zeffCO
                              + cefCH4 / zeffCH4 + 2 * cefC2H4 / zeffC2H4))
            return -(j / F) * (cefHCOO / zeffHCOO + cefCO / zeffCO
                               + cefCH4 / zeffCH4 + 2 * cefC2H4 / zeffC2H4)
        if name == "OH":
            return (j / F) * (cefHCOO / zeffHCOO + 2 * cefCO / zeffCO
                              + 8 * cefCH4 / zeffCH4 + 12 * cefC2H4 / zeffC2H4
                              + 2 * cefH2 / zeffH2)  # note error in Gupta et al.

    def set_boundary_markers(self):
        self.boundary_markers = {"bulk dirichlet": (1,),  # U_solid = U_app, C = C_0
                                 "neumann": (2,),
                                 }


solver = CarbonateSolver()
solver.setup_solver()
solver.solve()

### Plotting

C_CO2, C_OH, C_H, C_CO3, C_HCO3 = solver.u.subfunctions
# OH boundary layer
x = solver.mesh.coordinates
C_OH_bl = Function(solver.V).assign(C_OH).dat.data[100:]
x_bl = x.dat.data[100:]

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
filename = "carbonate.png"
fig = plt.figure(constrained_layout=True, figsize=(16, 8))
spec = gridspec.GridSpec(ncols=3, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1])
ax3 = fig.add_subplot(spec[1, 0])
ax4 = fig.add_subplot(spec[1, 1])
ax5 = fig.add_subplot(spec[0, 2])
ax6 = fig.add_subplot(spec[1, 2])

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
plot(C_H, axes=ax5)
ax5.set(xlabel='distance (m)',
        ylabel='H concentration (M)')
plt.plot(x_bl, C_OH_bl, axes=ax6, color='k', linewidth=2)
ax6.set(xlabel='distance (m)',
        ylabel='OH concentration (M)')

plt.savefig(filename)
