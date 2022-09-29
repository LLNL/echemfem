import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from firedrake import *
from echemfem import EchemSolver
import numpy as np
import csv
"""
Model fromn
Hawks, S.A., Ehlinger, V.M., Moore, T., Duoss, E.B., Beck, V.A., Weber, A.Z.
and Baker, S.E., 2022. Analyzing Production Rate and Carbon Utilization
Trade-offs in CO2RR Electrolyzers. ACS Energy Letters, 7, pp.2685-2693.
"""
# operating conditions
T = 293.15             # temperature (K)
P = 1                  # pressure (atm)
Vcell = Constant(-0.7) # potential (V)
S = 0.64               # saturation

# catalyst layer properties
eps_CL = 0.3        # catalyst layer porosity
sigma_CL = 100      # solid conductivity (S/m)
L_CL = 2.75e-7      # catalyst layer thickness (m)
RF = 11.37          # roughness
av0 = RF/L_CL       # specific surface area (1/m)

# physical constants
R = 8.3144598       # ideal gas constant (J/mol/K)
F = 96485.33289     # Faraday's constant (C/mol)
H_CO2 = 34.         # Henry's constant for CO2 (mol/m3/atm)

# CO2RR
i0CO2RR = 6.11e-6   # exchange current density (A/m2)
alphacCO2RR = 0.7    # cathode coefficient
U0CO2RR = 0.023257   # standard potential (V)

# HER
i0HER = 6.16e-7  # exchange current density (A/m2)
alphacHER = 0.6     # cathode coefficient 
U0HER = 0.0         # standard potential (V)

# electrolyte properties
cref = 1e3          # reference concentration (mol/m3)
Lelec = 7e-3        # electrode thickness (m)
Aelec = 1.44e-4     # electrode cross-sectional area (m2)
ql = 8.333e-8       # electrolyte flow rate (m3/s)
mul = 8.9e-4        # electrolyte viscosity (Pa*s)
rhol = 1000         # electrolyte density (kg/m3)
nul = ql / Aelec    # electrolyte flow velocity (m/s)

# diffusion coefficients (m2/s)
# [CO2 OH- H+ CO32- HCO3- K+]	
D = [1.910e-9, 1.185e-9, 9.311e-9, 0.923e-9, 1.185e-9, 1.957e-9]	

# mass transfer coefficients
kMT = []
for Di in D:
    kMT.append(0.664 * (Di / Lelec) * pow((rhol * nul * Lelec) / mul, 0.5) *
                 pow(mul / (rhol * Di), 1 / 3))

# charge of each species
z = [0, -1, 1, -2, -1, 1]

# equilibrium coefficients
Kw = 1e-8           # mol^2/m^6
K1 = pow(10,-3.37)  # mol/m3
K2 = pow(10,-7.32)  # mol/m3
K3 = K1/Kw
K4 = K2/Kw

# forward rate constants
k1f = 3.71e-2       # s^-1
k2f = 59.44         # s^-1
k3f = 2.23          # m3/(mol*s)
k4f = 6e6           # m3/(mol*s)
kwf = 1.4           # mol/m3/s

# backward rate constants
k1b = k1f/K1
k2b = k2f/K2
k3b = k3f/K3
k4b = k4f/K4
kwb = kwf/Kw

# Bulk Electrolyte concentrations
pH0 = 7.9
cKb = 1000                # mol/m3
c0Hb = pow(10,-pH0+3)     # mol/m3
c0OHb = Kw/c0Hb
c0CO3b = cKb/(1+c0Hb/K2+c0Hb**2/(K1*K2))
c0HCO3b = c0Hb*c0CO3b/K2
c0CO2b = H_CO2
cb = np.array([c0CO2b, c0OHb, c0Hb, c0CO3b, c0HCO3b, cKb])

class PorousSolver(EchemSolver):
    def __init__(self):

        def bulk_reaction(y):
            yCO2 = y[0]
            yOH = y[1]
            yH = y[2]
            yCO3 = y[3]
            yHCO3 = y[4]

            dCO2 = (k1b * yH + k3b) * yHCO3 - (k1f + k3f * yOH) * yCO2

            dHCO3 = (k1f + k3f * yOH) * yCO2 + (k4b + k2b * yH) * yCO3 \
                - (k1b * yH + k2f + k3b + k4f * yOH) * yHCO3

            dCO3 = (k2f + k4f * yOH) * yHCO3 - (k2b * yH + k4b) * yCO3

            dH = k1f * yCO2 + k2f * yHCO3 + kwf \
                - (k1b * yHCO3 + k2b * yCO3 + kwb * yOH) * yH

            dOH = k3b * yHCO3 + k4b * yCO3 + kwf \
                - (k3f * yCO2 + k4f * yHCO3 + kwb * yH) * yOH

            return [
                eps_CL * S * dCO2,
                eps_CL * S * dOH,
                eps_CL * S * dH,
                eps_CL * S * dCO3,
                eps_CL * S * dHCO3,
                0.]

        def reaction_CO2RR(u):
            CCO2 = u[0]
            CH = u[2]
            Phi1 = u[7]
            Phi2 = u[6]

            UCO2RR = U0CO2RR + ((R * T) / F) * ln(CH / cref)
            etaCO2RR = Phi1 - Phi2 - UCO2RR  # reaction overpotential (V)
            iCO2RR = i0CO2RR * exp(-((alphacCO2RR * F) / (R * T)) * etaCO2RR)
            return iCO2RR

        def reaction_HER(u):
            CH = u[2]
            Phi1 = u[7]
            Phi2 = u[6]
            UHER = U0HER + ((R * T) / F) * ln(CH / cref)
            etaHER = Phi1 - Phi2 - UHER  # reaction overpotential (V)
            iHER = i0HER * exp(-((alphacHER * F) / (R * T)) * etaHER)
            return iHER

        mesh = IntervalMesh(500, L_CL)
        conc_params = []

        conc_params.append({"name": "CO2",
                            "diffusion coefficient": D[0],
                            "z": z[0],
                            "bulk": cb[0],
                            "mass transfer coefficient": kMT[0],
                            "gas": H_CO2,
                            })

        conc_params.append({"name": "OH",
                            "diffusion coefficient": D[1],
                            "z": z[1],
                            "bulk": cb[1],
                            "mass transfer coefficient": kMT[1],
                            })

        conc_params.append({"name": "H",
                            "diffusion coefficient": D[2],
                            "z": z[2],
                            "bulk": cb[2],
                            "mass transfer coefficient": kMT[2],
                            })

        conc_params.append({"name": "CO3",
                            "diffusion coefficient": D[3],
                            "z": z[3],
                            "bulk": cb[3],
                            "mass transfer coefficient": kMT[3],
                            })

        conc_params.append({"name": "HCO3",
                            "diffusion coefficient": D[4],
                            "z": z[4],
                            "bulk": cb[4],
                            "mass transfer coefficient": kMT[4],
                            })

        conc_params.append({"name": "K",
                            "diffusion coefficient": D[5],
                            "z": z[5],
                            "bulk": cb[5],
                            "mass transfer coefficient": kMT[5],
                            })

        physical_params = {"flow": ["diffusion", "migration", "electroneutrality full", "porous"],
                           "F": F,  # C/mol
                           "R": R,  # J/K/mol
                           "T": T,  # K
                           "U_app": Vcell,  # V
                           "porosity": eps_CL,
                           "saturation": S,
                           "specific surface area": av0,  # 1/m
                           "solid conductivity": sigma_CL,  # S/m
                           "bulk reaction": bulk_reaction,
                           }
        echem_params = []

        echem_params.append({"reaction": reaction_CO2RR,
                             "electrons": 6,
                             "stoichiometry": {"CO2": -1, # reactant
                                               "OH": 6}, # product
                             })

        echem_params.append({"reaction": reaction_HER,
                             "electrons": 2,
                             "stoichiometry": {"OH": 2}, # product
                             })

        super().__init__(conc_params, physical_params, mesh, echem_params=echem_params, family="CG", p=2)

    def set_boundary_markers(self):
        self.boundary_markers = {"applied": (2,),  # U_solid = U_app
                                 "gas": (2,),  # C_CO2 = C_CO2_gas
                                 "bulk": (1,),  # U_liquid = 0, NC = k_x (C_0 - C)
                                 }


solver = PorousSolver()
solver.setup_solver(initial_solve=False)
#from IPython import embed; embed()
Vlist = np.linspace(-0.70, -1.5, num=41)
#Vlist = [-0.7]
sol = []
icell = []
V = solver.V
V1 = FunctionSpace(solver.mesh, "CG", 1)
for Vs in Vlist:
    solver.U_app.assign(Vs)
    print("V = %d mV" % (Vs * 1000))
    solver.solve()
    cCO2, cOH, cH, cCO3, cHCO3, cK, phi2, phi1 = solver.u.split()

    x_ = Function(V).interpolate(solver.mesh.coordinates[0]).vector().dat.data

    i1 = Function(V).interpolate(-solver.effective_diffusion(sigma_CL, phase="solid") * grad(phi1)[0])

    NCO2 = Function(V).interpolate( solver.effective_diffusion( D[0]) *
            grad(cCO2)[0])
    NOH = Function(V).interpolate((solver.effective_diffusion( D[1]) *
            grad(cOH) + F * z[1] * D[1] / R / T * cOH * grad(phi2))[0])
    NH = Function(V).interpolate((solver.effective_diffusion(D[2]) *
            grad(cH) + F * z[2] * solver.effective_diffusion(D[2]) / R / T * cH
            * grad(phi2))[0])
    NCO3 = Function(V).interpolate((solver.effective_diffusion(D[3]) *
            grad(cCO3) + F * z[3] * solver.effective_diffusion(D[3]) / R / T *
            cCO3 * grad(phi2))[0])
    NHCO3 = Function(V).interpolate((solver.effective_diffusion(D[4]) *
            grad(cHCO3) + F * z[4] * solver.effective_diffusion(D[4]) / R / T *
            cHCO3 * grad(phi2))[0])
    NK = Function(V).interpolate((solver.effective_diffusion(D[5]) *
            grad(cK) + F * z[5] * solver.effective_diffusion(D[5])/ R / T * cK
            * grad(phi2))[0])

    i2 = Function(V).assign(- F * (z[5] * NK + z[2] * NH + z[4] * NHCO3 +
              z[1] * NOH + z[3] * NCO3 + z[0] * NCO2))

    fig = plt.figure(constrained_layout=True, figsize=(8, 6))
    spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[0, 1])
    ax3 = fig.add_subplot(spec[1, 0])
    ax4 = fig.add_subplot(spec[1, 1])

    plot(i1, label="$i_s$", axes = ax1)
    plot(i2, label="$i_l$", axes = ax1, color='b')
    ax1.legend()
    ax1.set(xlabel='catalyst layer distance (m)',
            ylabel='current density (A/m$^2$)')

    plot(phi1, label="$\phi_s$", axes = ax2)
    plot(phi2, label="$\phi_l$", axes = ax2, color = 'b')
    ax2.legend()
    ax2.set(xlabel='catalyst layer distance (m)', ylabel='potential (V)')

    plot(cCO2, axes = ax3)
    ax3.set(xlabel='catalyst layer distance (m)',
            ylabel='CO$_2$ concentration (mol/m$^3$)')

    plot(NCO2, axes = ax4)
    ax4.set(xlabel='catalyst layer distance (m)',
            ylabel='CO$_2$ flux (mol/m$^2$/s)')
    filename = "results/current%dmV.png" % (-Vs * 1000)
    plt.savefig(filename)

    fig = plt.figure(constrained_layout=True, figsize=(12, 6))
    spec = gridspec.GridSpec(ncols=3, nrows=2, figure=fig)
    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[0, 1])
    ax3 = fig.add_subplot(spec[0, 2])
    ax4 = fig.add_subplot(spec[1, 0])
    ax5 = fig.add_subplot(spec[1, 1])
    ax6 = fig.add_subplot(spec[1, 2])

    plot(cK, axes = ax1)
    ax1.set(xlabel='catalyst layer distance (m)',
            ylabel='K$^+$ concentration (mol/m$^3$)')

    plot(cH, axes = ax2)
    ax2.set(xlabel='catalyst layer distance (m)',
            ylabel='H$^+$ concentration (mol/m$^3$)')

    plot(cHCO3, axes = ax3)
    ax3.set(xlabel='catalyst layer distance (m)',
            ylabel='HCO$_3^-$ concentration (mol/m$^3$)')

    plot(cOH, axes = ax4)
    ax4.set(xlabel='catalyst layer distance (m)',
            ylabel='OH$^-$ concentration (mol/m$^3$)')

    plot(cCO3, axes = ax5)
    ax5.set(xlabel='catalyst layer distance (m)',
            ylabel='CO$_3^{2-}$ concentration (mol/m$^3$)')

    plot(cCO2, axes = ax6)
    ax6.set(xlabel='catalyst layer distance (m)',
            ylabel='CO$_2$ concentration (mol/m$^3$)')
    filename = "results/electrolyte%dmV.png" % (-Vs * 1000)
    plt.savefig(filename)



    phi1_ = phi1.vector().dat.data
    phi2_ = phi2.vector().dat.data
    cK_ = cK.vector().dat.data
    cH_ = cH.vector().dat.data
    cHCO3_ = cHCO3.vector().dat.data
    cOH_ = cOH.vector().dat.data
    cCO3_ = cCO3.vector().dat.data
    cCO2_ = cCO2.vector().dat.data

    i1_ = i1.dat.data
    i2_ = i2.dat.data

    Vname = np.rint(-Vs*1000)
    filename = "full_gde_data%dmV.csv" % Vname
    with open(filename, mode = 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['x','i1','V1','i2','V2','cK','cH','cHCO3','cOH','cCO3','cCO2'])            
        for i in range(0,len(x_)):
            j = i
            if j > len(x_)-1:
                j = len(x_)-1
            csvwriter.writerow([x_[i],i1_[j],phi1_[i],i2_[j],phi2_[i],cK_[i],cH_[i],cHCO3_[i],cOH_[i],cCO3_[i],cCO2_[i]])
    icell.append(i1_[-1])
