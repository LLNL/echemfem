from echemfem import EchemSolver
import argparse
from firedrake import *
from mpi4py import MPI
from firedrake.petsc import PETSc
PETSc.Sys.popErrorHandler()

set_log_level(DEBUG)


parser = argparse.ArgumentParser()
parser.add_argument('-ref_levels', type=int, default=0)
parser.add_argument('-stats_file', type=str, default='')
parser.add_argument('-overwrite_stats_file', type=bool, default=False)
parser.add_argument('-csolver', type=str, default='asm')
parser.add_argument('-degree', type=int, default=1)
args, unknown = parser.parse_known_args()

GLOBAL_N_RANKS = MPI.Comm.Get_size(MPI.COMM_WORLD)
if GLOBAL_N_RANKS % 8 == 0:
    REDFACTOR = int(GLOBAL_N_RANKS / 8)
else:
    REDFACTOR = 1


class BortelsSolver(EchemSolver):
    def __init__(self):
        # Create an initial coarse mesh
        plane_mesh = Mesh('bortels_structuredquad_nondim_coarse1664.msh')
        plane_mesh_hierarchy = MeshHierarchy(
            plane_mesh, refinement_levels=args.ref_levels)
        hz = 6.0  # non-dim z length
        extruded_hierarchy = ExtrudedMeshHierarchy(
            plane_mesh_hierarchy, hz, 8)
        mesh = extruded_hierarchy[-1]

        D_Cu = 7.2e-10  # m^2/s
        J0 = 30.  # A/m^2
        C_Cu = 10.  # mol/m^3
        D_SO4 = 10.65e-10  # m^2/s
        C_SO4 = 1010.  # mol/m^3
        D_H = 93.12e-10  # m^2/s
        C_H = 2000.  # mol/m^3
        F = 96485.3329  # C/mol
        R = 8.3144598  # J/K/mol
        T = 273.15 + 25.  # K
        U_app = 0.03  # V
        v_avg = 0.03  # m/s

        L = 0.01  # m

        D_Cu_hat = D_Cu / L / v_avg
        D_SO4_hat = D_SO4 / L / v_avg
        D_H_hat = D_H / L / v_avg
        J0_hat = J0 / v_avg / C_Cu / F
        U_app_hat = F / R / T * U_app

        _, _, z = SpatialCoordinate(mesh)
        J_parab = J0_hat * 3 / 5 * (2 - (z - 0.5 * hz)**2 / (0.5 * hz)**2)

        conc_params = []
        conc_params.append({"name": "Cu",
                            "diffusion coefficient": D_Cu_hat,
                            "z": 2.,
                            "bulk": 1.,
                            "C_ND": 1.,
                            })

        conc_params.append({"name": "SO4",
                            "diffusion coefficient": D_SO4_hat,
                            "z": -2.,
                            "bulk": 1.,
                            "C_ND": C_SO4 / C_Cu,
                            })

        conc_params.append({"name": "H",
                            "diffusion coefficient": D_H_hat,
                            "z": 1.,
                            "bulk": 1.,
                            "C_ND": C_H / C_Cu,
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
            U = u[2]
            C_b = conc_params[0]["bulk"]
            J0 = J_parab
            F = physical_params["F"]
            R = physical_params["R"]
            T = physical_params["T"]
            eta = V - U
            return J0 * (exp(F / R / T * (eta)) -
                            (C / C_b) * exp(-F / R / T * (eta)))

        def reaction_cathode(u):
            return reaction(u, Constant(0))

        def reaction_anode(u):
            return reaction(u, Constant(physical_params["U_app"]))

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
            mesh,
            echem_params=echem_params,
            stats_file=args.stats_file,
            overwrite_stats_file=args.overwrite_stats_file,
            p=args.degree)

        U_is = self.num_mass
        is_list = [str(i) for i in range(self.num_mass)]
        C_is = ",".join(is_list)

        asm = {"pc_type": "asm",
               "pc_asm_overlap": 1,
               "sub": {
                   "pc_type": "ilu",
                   "pc_factor_levels": 0,
               }
               }
        gmg = {"pc_type": "mg",
               "mg_levels_ksp_type": "richardson",
               "mg_levels": asm,
               "mg_coarse_ksp_type": "preonly",
               "mg_coarse": {
                   "pc_type": "python",
                   "pc_python_type": "firedrake.AssembledPC",
                   "assembled": {
                       "mat_type": "aij",
                       "pc_type": "telescope",
                       "pc_telescope_reduction_factor": REDFACTOR,
                       "pc_telescope_subcomm_type": "contiguous",
                       "telescope_pc_type": "lu",
                       "telescope_pc_factor_mat_solver_type": "mumps",
                   }
               }
               }
        if args.csolver == "asm":
            csolver = asm
        else:
            csolver = gmg

        self.init_solver_parameters(
            custom_solver={
                "mat_type": "aij",
                "snes_view": None,
                "snes_converged_reason": None,
                "snes_monitor": None,
                "snes_rtol": 1e-6,
                "ksp_monitor": None,
                "ksp_converged_reason": None,
                "ksp_type": "fgmres",
                "ksp_rtol": 1e-3,
                "pc_type": "fieldsplit",
                "pc_fieldsplit_0_fields": U_is,
                "pc_fieldsplit_1_fields": 0,
                "pc_fieldsplit_2_fields": 1,
                "fieldsplit_0": {
                    "ksp_converged_reason": None,
                    "ksp_rtol": 1e-1,
                    "ksp_type": "cg",
                    "pc_type": "hypre",
                    "pc_hypre_boomeramg": {
                        "strong_threshold": 0.7,
                        "coarsen_type": "HMIS",
                        "agg_nl": 3,
                        "interp_type": "ext+i",
                        "agg_num_paths": 5,
                    },
                },
                "fieldsplit_1": {**{
                    "ksp_converged_reason": None,
                    "ksp_rtol": 1e-1,
                    "ksp_type": "gmres",
                }, **csolver},
                "fieldsplit_2": {**{
                    "ksp_converged_reason": None,
                    "ksp_rtol": 1e-1,
                    "ksp_type": "gmres",
                }, **csolver},
            },
            custom_potential_solver={
                "mat_type": "aij",
                "snes_monitor": None,
                "snes_rtol": 1e-6,
                "ksp_converged_reason": None,
                "ksp_type": "cg",
                "ksp_rtol": 1e-2,
                "pc_type": "hypre",
                "pc_hypre_boomeramg": {
                    "strong_threshold": 0.7,
                    "coarsen_type": "HMIS",
                    "agg_nl": 3,
                    "interp_type": "ext+i",
                    "agg_num_paths": 5,
                },
            })

    def set_boundary_markers(self):
        self.boundary_markers = {"inlet": (10,),
                                 "outlet": (11,),
                                 "anode": (12,),
                                 "cathode": (13,),
                                 }

    def set_velocity(self):
        _, y, z = SpatialCoordinate(self.mesh)
        h = 1.0
        self.vel = as_vector((6. *
                              self.physical_params["v_avg"] /
                              h**2 *
                              y *
                              (h -
                               y), Constant(0.), Constant(0.)))  # m/s


solver = BortelsSolver()

solver.setup_solver(initial_solve=True)
solver.solve()
