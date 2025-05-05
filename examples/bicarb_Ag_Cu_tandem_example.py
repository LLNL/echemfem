from firedrake import *
from echemfem import EchemSolver
import numpy as np
import matplotlib.pyplot as plt
from math import log10
import os
import shutil

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

"""
This example case runs a tandem electrode system, with an Ag catalyst placed
in front of a Cu catalyst, at a Peclet number of 1E3. 
The advection/diffusion/migration/electroneutrality equations are solved 
for the bulk electrolyte, and a Tafel description is applied for the surface 
catalyst reactions.

Assumed products of Cu catalyst are C2H4, C2H6O, and H2.
Assumed products of Ag catalysis are CO and H2.

--------------------------------------------------------------------------
list of references

Bicarbonate flow reactor setup given in:
Lin, T.Y., Baker, S.E., Duoss, E.B. and Beck, V.A., 2021. Analysis of
the Reactive CO2 Surface Flux in Electrocatalytic Aqueous Flow
Reactors. Industrial & Engineering Chemistry Research, 60(31),
pp.11824-11833.

Ag Tafel parameter values taken from: 
Corpus, K. R. M., Bui, J. C., Limaye, A. M., Pant, L. M., Manthiram, K., 
Weber, A. Z., and Bell, A. T., 2023. Coupling covariance matrix adaptation 
with continuum modeling for determination of kinetic parameters associated 
with electrochemical CO2 reduction. Joule, 7, pp.1289-1307.

Cu Tafel parameter values taken from:
Li, J., Chang, X., Zhang, H., Malkani, A. S., Cheng, M. J., Xu, B., 
and Lu, Q., 2021. Electrokinetic and in situ spectroscopic investigations 
of CO electrochemical reduction on copper. Nature Communications, 12(1), pp.3264.

Diffusivities values taken from:
Weng, L. C., Bell, A. T., & Weber, A. Z., 2018. Modeling gas-diffusion 
electrodes for CO2 reduction. Physical Chemistry Chemical Physics, 20(25), 
pp. 16973-16984.
and from: 
https://en.wikipedia.org/wiki/Mass_diffusivity#cite_note-Cussler-3

The mesh has been taken from that used in (and then slightly modified for this smaller domain size):
Govindarajan, N., Lin, T. Y., Roy, T., Hahn, C., & Varley, J. B., 2023. Coupling 
Microkinetics with Continuum Transport Models to Understand Electrochemical 
CO2 Reduction in Flow Reactors. PRX Energy, 2(3), pp.033010.

--------------------------------------------------------------------------
end of list of references
"""


#---------------------------------------------------------------------------
#constants, operating conditions, extraneous values used
# operating conditions
T = 298 # 293.15             # temperature (K)
Vcell = Constant(-0.5) # potential (V)

# physical constants
R = 8.3144598       # ideal gas constant (J/mol/K)
F = 96485.33289     # Faraday's constant (C/mol)

#other constants
cref = 1 # reference concentration (mol/m^3)
cref_Ag = 1e3 #reference concentration for Ag Tafel (mol/m^3)
#---------------------------------------------------------------------------
#end ofconstants, operating conditions, extraneous values used



##--------------------------------------------------------------------------
#these parameters are for Ag Tafel expressions

i0_H2_Ag = 1.698E-6 # exchange current density (A/m^2)
alpha_c_CO = 0.544 # transfer coefficient

alpha_c_H2_Ag = 0.312 # transfer coefficient
gamma_OH_H2 = alpha_c_H2_Ag # exponent used in Tafel expression

U0_CO = -0.11 # standard potential (V)
U0_H2_Ag = 0.0 # standard potential (V)

alpha_c_CO = 0.544 # transfer coefficient
i0_CO = 1.905E-6 # exchange current density (A/m^2)
gamma_CO2_CO = -(2-alpha_c_CO)/2 # exponent used in Tafel expression
gamma_OH_CO = alpha_c_CO # exponent used in Tafel expression

##Sechenov coefficients for CO2 activity
h_s_OH = 8.39E-5 # (m^3/mol)
h_s_HCO3 = 9.67E-5 # (m^3/mol) ##!!!note that this is slightly different from the reported value of 
#5.49E-5 as found in https://jpldataeval.jpl.nasa.gov/pdf/Jpl15_Sectn5_HeterogenChem.pdf
h_s_CO3 = 14.23E-5 # (m^3/mol)
h_s_H = 0. # (m^3/mol)
h_s_K = 9.2E-5 # (m^3/mol) #reported as 9.22E-5 (different # of sig figs) in https://doi.org/10.1021/acscatal.9b05272
h_g_CO2 = -1.7159E-5 #(m^3/mol)

##--------------------------------------------------------------------------
#end of Ag Tafel expression parameters



##--------------------------------------------------------------------------
#these parameters here are for Cu Tafel expressions
#note: parameter extracted using fitted (rather than reported by Li et al Nat Comm) Tafel slopes
#uses the Tafel values reported for pH = 7.2

pH = 7.2 # pH of data used in fit

# C2H4
i0_C2H4 = 1.3835E-11 # exchange current density (A/m^2)
alpha_c_C2H4 = 0.4990 # transfer coefficient
UC2H4 = 0.17 # standard potential (V)

# C2H6O
i0_C2H6O = 8.9302E-11 # exchange current density (A/m^2)
alpha_c_C2H6O = 0.4452 # 
UC2H6O = 0.19 # standard potential (V)

# H2
i0_H2_Cu = 4.6270E-5 # exchange current density (A/m^2)
alpha_c_H2_Cu = 0.2547 #(-)
UH2 = 0. # standard potential (V)
##--------------------------------------------------------------------------
#end of Cu Tafel expression parameters



#---------------------------------------------------
#mesh/domain parameters start
La = Constant(0.00025) # offset for inlet developing region (m)
Lx = Constant(0.001) # total catalyst length (m)
total_domain_length = Constant(0.0015) # total domain length (m)
#---------------------------------------------------
#mesh/domain parameters end



class CarbonateSolver(EchemSolver):
    def __init__(self):
        
        
        mesh = Mesh('mesh_bicarb_tandem_example.msh')
        #Note: will need to use "gmsh -2 mesh_bicarb_tandem_example.geo" command to convert to .msh format


        #---------------------------------------------------
        #begin create tanh functions to locally specify Ag and Cu locations
        x, y = SpatialCoordinate(mesh)##to get coordinates of mesh
        sharpening_factor = 80e3##larger value means a steeper slope; the value needs to be of value O(10 * 1/L_x) in units of (1/m) to sufficiently demarcate between electrode and defect
        
        ##initialize Ag BC with tanh
        BC_local_Ag = 0.5 * ( tanh(sharpening_factor*(x-La))  -  tanh(sharpening_factor*(x-Lx/2-La)) )

        ##initialize Cu BC with tanh
        BC_local_Cu = 0.5 * ( tanh(sharpening_factor*(x-Lx/2-La))  -  tanh(sharpening_factor*(x-Lx-La)) )

        #---------------------------------------------------
        #end create tanh functions to differentiate Ag and Cu



        k1f = (8.42E3)/1000 # (L/mol/s)
        k1r = 1.97E-4 # (1/s)
        k2f = 3.71E-2 # (1/s)
        k2r = (8.697E4)/1000 # (L/mol/s)
        k3f = 1.254E6 # (1/s)
        k3r = (6E9)/1000 # (L/mol/s)
        k4f = (1.242E12)/1000 # (L/mol/s)
        k4r = 59.44 # (1/s)
        k5f = (2.3E10)/1000 # (L/mol/s)
        k5r = (2.3E-4)*1000 # (mol/L/s)

        C_1_inf = 0.034*1000 # (mol/m^3), obtained from Henry's law constant at p = 1atm 
        C_K =  0.5*1000  # (mol/m^3)
        Keq = k1f*k3f/(k1r*k3r)

        C_3_inf = -C_1_inf*Keq/4 + C_1_inf*Keq/4*sqrt(1+8*C_K/(C_1_inf*Keq)) # (mol/m^3)
        C_4_inf = (C_K - C_3_inf)/2 # (mol/m^3)
        C_2_inf = C_3_inf/C_1_inf*k1r/k1f # (mol/m^3)
        C_5_inf = (k5r/k5f)/C_2_inf # (mol/m^3)

        C_CO2_bulk = C_1_inf # (mol/m^3)
        C_OH_bulk = C_2_inf # (mol/m^3)
        C_HCO3_bulk = C_3_inf # (mol/m^3)
        C_CO32_bulk = C_4_inf # (mol/m^3)
        C_H_bulk = C_5_inf # (mol/m^3)
        C_K_bulk = C_K # (mol/m^3)

        #Make sure that electroneutrality holds by manually adjusting K concentration
        netcharge = C_K_bulk+ C_H_bulk -2.*C_CO32_bulk - C_HCO3_bulk - C_OH_bulk
        print("Total charge is:",netcharge )
        C_K_bulk = C_K_bulk - netcharge
        netcharge = C_K_bulk+ C_H_bulk -2.*C_CO32_bulk - C_HCO3_bulk - C_OH_bulk
        print("Total charge is:",netcharge )


        def bulk_reaction(y):
            yCO2=y[0];
            yOH=y[1];
            yHCO3=y[2];
            yCO3=y[3];
            yH=y[4];
            
            
            dCO2 = -(k1f)   *yCO2*yOH \
                    +(k1r)    *yHCO3 \
                    -(k2f)  *yCO2 \
                    +(k2r)   *yHCO3*yH

            dOH = -(k1f)       *yCO2*yOH \
                           +(k1r)        *yHCO3 \
                           +(k3f)        *yCO3 \
                           -(k3r)       *yOH*yHCO3\
                           -(k5f)      *yOH*yH\
                           +(k5r)

            dHCO3 = (k1f)        *yCO2*yOH\
                           -(k1r)       *yHCO3\
                           +(k3f)        *yCO3\
                           -(k3r)       *yOH*yHCO3\
                           +(k2f)       *yCO2 \
                           -(k2r)      *yHCO3*yH \
                           +(k4f)       *yCO3*yH\
                           -(k4r)      *yHCO3

            dCO3 = -(k3f) *yCO3 \
                           +(k3r)  *yOH*yHCO3\
                           -(k4f)*yCO3*yH\
                           +(k4r) *yHCO3

            dH = (k2f)   *yCO2 \
                           -(k2r)  *yHCO3*yH \
                           -(k4f)  *yCO3*yH\
                           +(k4r)   *yHCO3\
                           -(k5f)  *yOH*yH\
                           +(k5r)
               
            return [dCO2, dOH, dHCO3, dCO3, dH, 0., 0., 0., 0., 0.]

        


        conc_params = []

        #z is the charge for each species

        conc_params.append({"name": "CO2",
                            "diffusion coefficient": 1.91E-9,  # m^2/s
                            "bulk": C_CO2_bulk,  # (mol/m^3)
                            "z": 0,
                            })

        conc_params.append({"name": "OH",
                            "diffusion coefficient": 5.29E-9,  # m^2/s
                            "bulk": C_OH_bulk,  # (mol/m^3)
                            "z": -1,
                            })

        conc_params.append({"name": "HCO3",
                            "diffusion coefficient": 1.185E-9,  # m^2/s
                            "bulk": C_HCO3_bulk,  # (mol/m^3)
                            "z": -1,
                            })

        conc_params.append({"name": "CO3",
                            "diffusion coefficient": .92E-9,  # m^2/s
                            "bulk": C_CO32_bulk,  # (mol/m^3)
                            "z": -2,
                            })

        conc_params.append({"name": "H",
                            "diffusion coefficient": 9.311E-9,  # m^2/s
                            "bulk": C_H_bulk,  # (mol/m^3)
                            "z": 1,
                            })

        conc_params.append({"name": "H2",
                            "diffusion coefficient": 4.5E-9, #1.96E-9,  # m^2/s
                            "bulk": 0.0,  # (mol/m^3)
                            "z": 0,
                            })

        conc_params.append({"name": "CO",
                            "diffusion coefficient": 2.03E-9, #1.96E-9,  # m^2/s
                            "bulk": 0.0,  # (mol/m^3)
                            "z": 0,
                            })

        conc_params.append({"name": "C2H4",
                            "diffusion coefficient": 1.87E-9,  # m^2/s
                            "bulk": 0.,  # (mol/m^3)
                            "z": 0,
                            })

        conc_params.append({"name": "C2H6O",
                            "diffusion coefficient": 0.84E-9,  # m^2/s
                            "bulk": 0.,  # (mol/m^3)
                            "z": 0, 
                            })

        conc_params.append({"name": "K",
                    "diffusion coefficient": 1.957E-9, #1.96E-9,  # m^2/s
                    "bulk": C_K_bulk,  # (mol/m^3)
                    "z": 1,
                    })

        physical_params = {"flow": ["advection", "diffusion", "migration", "electroneutrality"],
                   "F": F,  # (C/mol)
                   "R": R,  # (J/K/mol)
                   "T": T,  # (K)
                   "U_app": Vcell, # applied voltage, (V)
                   "bulk reaction": bulk_reaction,
                   }

        def reaction_CO(u):

            CCO2 = u[0]
            COH = u[1]
            CHCO3 = u[2]
            CCO3 = u[3]
            CH = u[4]
            CH2 = u[5]
            CCO = u[6]
            CK = -u[4] + 2*u[3] + u[2] + u[1]

            Phi2 = u[9]
            Phi1 = physical_params["U_app"]
            UCO = U0_CO

            ##---------------------------calculate activities
            #charges of species
            z_CO2 = 0.
            z_OH = -1.
            z_HCO3 = -1.
            z_CO3 = -2.
            z_H = 1.
            z_H2 = 0.
            z_CO = 0.
            z_K = 1.

            #local ionic strength to calculate a_OH and a_H
            I =  0.5 * ( CCO2*pow(z_CO2,2) + COH*pow(z_OH,2) + CHCO3*pow(z_HCO3,2) + CCO3*pow(z_CO3,2) + CH*pow(z_H,2) + CK*pow(z_K,2)) #local ionic strength (mol/m^3)
            I = I/1000. #convert from mol/m^3 to M
            f_OH = pow(10.,-0.51 * pow(z_OH,2) * (pow(I,0.5)/(1. + pow(I,0.5)) - 0.3 * I ) )
            f_H = pow(10.,-0.51 * pow(z_H,2) * (pow(I,0.5)/(1. + pow(I,0.5)) - 0.3 * I ) )

            a_OH = f_OH * COH / cref_Ag
            a_OH_bulk = C_OH_bulk / cref_Ag

            a_H = f_H * CH / cref_Ag
            pH = -ln(a_H)/ln(10)

            f_CO2 = exp(COH*(h_s_OH+h_g_CO2) + CHCO3*(h_s_HCO3+h_g_CO2) + CCO3*(h_s_CO3+h_g_CO2) + CH*(h_s_H+h_g_CO2) + CK*(h_s_K+h_g_CO2))
            a_CO2 = f_CO2 * CCO2 / cref_Ag
            a_CO2_bulk = C_CO2_bulk / cref_Ag

            ##-------------------------end of calculate activities            

            eta_CO = Phi1 - Phi2 - (UCO - ((2.303*R*T)/F)*(pH) + (R*T)/(2*F)*(ln(a_CO2))) # reaction overpotential (V)
            iCO = i0_CO * (a_CO2 / a_CO2_bulk)**(-gamma_CO2_CO) * (a_OH/a_OH_bulk)**(gamma_OH_CO) * exp(-((alpha_c_CO * F) / (R * T)) * eta_CO)
            return iCO * BC_local_Ag


        def reaction_H2_Ag(u):

            CCO2 = u[0]
            COH = u[1]
            CHCO3 = u[2]
            CCO3 = u[3]
            CH = u[4]
            CH2 = u[5]
            CCO = u[6]
            CK = -u[4] + 2*u[3] + u[2] + u[1]

            Phi2 = u[9]
            Phi1 = physical_params["U_app"]
            UCO = U0_CO

            ##---------------------------calculate activities

            #charges of species
            z_CO2 = 0.
            z_OH = -1.
            z_HCO3 = -1.
            z_CO3 = -2.
            z_H = 1.
            z_H2 = 0.
            z_CO = 0.
            z_K = 1.

            #local ionic strength to calculate a_OH and a_H
            I =  0.5 * ( CCO2*pow(z_CO2,2) + COH*pow(z_OH,2) + CHCO3*pow(z_HCO3,2) + CCO3*pow(z_CO3,2) + CH*pow(z_H,2) + CK*pow(z_K,2)) #local ionic strength (mol/m^3)
            I = I/1000. #convert from mol/m^3 to M
            f_OH = pow(10.,-0.51 * pow(z_OH,2) * (pow(I,0.5)/(1. + pow(I,0.5)) - 0.3 * I ) )
            f_H = pow(10.,-0.51 * pow(z_H,2) * (pow(I,0.5)/(1. + pow(I,0.5)) - 0.3 * I ) )

            a_OH = f_OH * COH / cref_Ag
            a_OH_bulk = C_OH_bulk / cref_Ag

            a_H = f_H * CH / cref_Ag
            pH = -ln(a_H)/ln(10)

            ##-------------------------end of calculate activities 

            eta_H2_Ag = Phi1 - Phi2 - (U0_H2_Ag - ((2.303*R*T)/F)*(pH))  # reaction overpotential (V)
            iH2_Ag = i0_H2_Ag * (a_OH/a_OH_bulk)**(gamma_OH_H2) * exp(-((alpha_c_H2_Ag * F) / (R * T)) * eta_H2_Ag)
            return iH2_Ag * BC_local_Ag

        def reaction_C2H4(u):
            C_CO = u[6]
            Phi2 = u[9]
            Phi1 = physical_params["U_app"]

            eta_C2H4 = Phi1 - Phi2 # reaction overpotential (V)
            iC2H4 = i0_C2H4 * (C_CO/cref) * exp(-((alpha_c_C2H4 * F) / (R * T)) * eta_C2H4)
            return iC2H4 * BC_local_Cu

        def reaction_C2H6O(u):
            C_CO = u[6]
            C_H2 = u[5]
            Phi2 = u[9]
            Phi1 = physical_params["U_app"]

            eta_C2H6O = Phi1 - Phi2 # reaction overpotential (V) 
            iC2H6O = i0_C2H6O * (C_CO/cref) * exp(-((alpha_c_C2H6O * F) / (R * T)) * eta_C2H6O)
            return iC2H6O * BC_local_Cu

        def reaction_H2_Cu(u):
            Phi2 = u[9]
            Phi1 = physical_params["U_app"]
    
            eta_H2 = Phi1 - Phi2 # reaction overpotential (V)
            iH2_Cu = i0_H2_Cu * exp(-((alpha_c_H2_Cu * F) / (R * T)) * eta_H2)
            return iH2_Cu * BC_local_Cu


        echem_params = []

        echem_params.append({"reaction": reaction_CO,
                             "electrons": 2,
                             "stoichiometry": {"CO2": -1,
                                               "OH": 2,
                                               "CO": 1},
                             "boundary": "catalyst",
                             })

        echem_params.append({"reaction": reaction_H2_Ag,
                             "electrons": 2,
                             "stoichiometry": {"OH": 2,
                                               "H2": 1},
                             "boundary": "catalyst",
                             })

        echem_params.append({"reaction": reaction_C2H4,
                             "electrons": 8,
                             "stoichiometry": {"CO": -2,
                                               "C2H4": 1,
                                               "OH": 8}, 
                             "boundary": "catalyst",
                             })

        echem_params.append({"reaction": reaction_C2H6O,
                             "electrons": 8,
                             "stoichiometry": {"CO": -2,
                                               "C2H6O": 1,
                                               "OH": 8}, 
                             "boundary": "catalyst",
                             })
        echem_params.append({"reaction": reaction_H2_Cu,
                             "electrons": 2,
                             "stoichiometry": {"OH": 2,
                                               "H2": 1}, 
                             "boundary": "catalyst",
                             })


        super().__init__(conc_params, physical_params, mesh, echem_params=echem_params, family="DG")
        #super().__init__(conc_params, physical_params, mesh, echem_params=echem_params, family="CG", SUPG = True) #this line would be used for CG instead of DG

    def set_boundary_markers(self):
        self.boundary_markers = {"inlet": (1),
                                 "bulk dirichlet": (3),
                                 "outlet": (2,),
                                 "catalyst": (4,),
                                 "bulk": (3,), 
                                 }


    def set_velocity(self):
        _, y = SpatialCoordinate(self.mesh)
        self.vel = as_vector([1.91*y,Constant(0)]) # m/s


solver = CarbonateSolver()
solver.U_app.assign(Vcell)
solver.setup_solver()

#solver.setup_solver(initial_solve=False)
n = FacetNormal(solver.mesh)
solver.solve()
cCO2, cOH, cHCO3, cCO3, cH, cH2, cCO, cC2H4, cC2H6O, phi2 = solver.u.subfunctions


Vlist = np.linspace(-0.5,-1.3,num=17)
jCO = []
jH2_Ag = []
jC2H4 = []
jC2H6O = []
jH2_Cu = []

for Vs in Vlist:
    solver.U_app.assign(Vs)
    print("V = %d mV" % np.rint(Vs * 1000))
    solver.solve()
    cCO2, cOH, cHCO3, cCO3, cH, cH2, cCO, cC2H4, cC2H6O, phi2 = solver.u.subfunctions


    iCO_avg = assemble(solver.echem_params[0]["reaction"](solver.u.subfunctions)*ds(4))/assemble(Constant(1) * ds(4, domain=solver.mesh))
    jCO.append(iCO_avg)
    print(iCO_avg)

    iH2_Ag_avg = assemble(solver.echem_params[1]["reaction"](solver.u.subfunctions)*ds(4))/assemble(Constant(1) * ds(4, domain=solver.mesh))
    jH2_Ag.append(iH2_Ag_avg)
    print(iH2_Ag_avg)

    iC2H4_avg = assemble(solver.echem_params[2]["reaction"](solver.u.subfunctions)*ds(4))/assemble(Constant(1) * ds(4, domain=solver.mesh))
    jC2H4.append(iC2H4_avg)
    print(iC2H4_avg)

    iC2H6O_avg = assemble(solver.echem_params[3]["reaction"](solver.u.subfunctions)*ds(4))/assemble(Constant(1) * ds(4, domain=solver.mesh))
    jC2H6O.append(iC2H6O_avg)
    print(iC2H6O_avg)

    iH2_Cu_avg = assemble(solver.echem_params[4]["reaction"](solver.u.subfunctions)*ds(4))/assemble(Constant(1) * ds(4, domain=solver.mesh))
    jH2_Cu.append(iH2_Cu_avg)
    print(iH2_Cu_avg)

    ##rename results folder so that I can save each of them
    base_dir = os.getcwd()
    os.chdir(base_dir)
    name_new_results_dir = "results_V_s_"+str(Vs)+"/"
    if rank == 0:
        dest = shutil.move("results",name_new_results_dir)


np.savetxt('U_SHE.tsv', Vlist)
np.savetxt('jCO.tsv',jCO)
np.savetxt('jH2_Ag.tsv',jH2_Ag)
np.savetxt('jC2H4.tsv',jC2H4)
np.savetxt('jC2H6O.tsv',jC2H6O)
np.savetxt('jH2_Cu.tsv',jH2_Cu)
