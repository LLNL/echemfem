from abc import ABC
from firedrake import *
from petsc4py import PETSc
pprint = PETSc.Sys.Print

class FlowSolver(ABC):
    """Base class for a flow solver.

    This class is used to create solvers for fluid flow decoupled from the chemistry
    """

    def __init__(self, mesh, fluid_params, boundary_markers):
        self.mesh = mesh
        self.dim = mesh.topological_dimension()
        self.fluid_params = fluid_params
        self.boundary_markers = boundary_markers
        self.setup_functions()
        self.setup_problem()
        self.setup_dirichlet_bcs()

    def setup_functions(self):
        # Taylor Hood
        mesh = self.mesh
        self.V = VectorFunctionSpace(mesh, "CG", 2)
        self.W = FunctionSpace(mesh, "CG", 1)
        self.Z = self.V * self.W
        self.soln = Function(self.Z, name="Flow solution")
        self.u, self.p = split(self.soln)
        self.v, self.q = TestFunctions(self.Z)

    def setup_problem(self):
        # set up PDE system - problem specific
        pass
    
    def setup_dirichlet_bcs(self):
        # set up Dirichlet boundary conditions - should be the same regardless of model

        Z = self.Z
        params = self.fluid_params
        bcs = []
        if self.boundary_markers.get("no slip"):
            if self.dim == 2:
                zero = Constant((0, 0))
            elif self.dim == 3:
                zero = Constant((0, 0, 0))
            bcs.append(DirichletBC(Z.sub(0), zero, self.boundary_markers["no slip"]))
        if self.boundary_markers.get("inlet velocity"):
           bcs.append(DirichletBC(Z.sub(0), params["inlet velocity"],
                                  self.boundary_markers["inlet velocity"]))
        if self.boundary_markers.get("outlet velocity"):
           bcs.append(DirichletBC(Z.sub(0), params["outlet velocity"],
                                  self.boundary_markers["outlet velocity"]))
        if self.boundary_markers.get("outlet pressure"):
           bcs.append(DirichletBC(Z.sub(1), params["outlet pressure"],
                                  self.boundary_markers["outlet pressure"]))
        if self.boundary_markers.get("inlet pressure"):
           bcs.append(DirichletBC(Z.sub(1), params["inlet pressure"],
                                  self.boundary_markers["inlet pressure"]))

        if self.boundary_markers.get("inlet pressure") or self.boundary_markers.get("outlet pressure"):
            self.nullspace = None
        else:
            self.nullspace = MixedVectorSpaceBasis( Z, [Z.sub(0),
                                                        VectorSpaceBasis(constant=True)])
        self.bcs = bcs

    def setup_solver(self):
        # set up nonlinear solver
        parameters = {"snes_monitor": None,
                      "ksp_type": "preonly",
                      "pc_type": "lu",
                      "pc_factor_mat_solver_type": "mumps",
                      }
        self.problem = NonlinearVariationalProblem(self.Form, self.soln,
                                                   bcs=self.bcs)
        self.solver = NonlinearVariationalSolver(self.problem,
                                                 solver_parameters=parameters,
                                                 nullspace=self.nullspace)

    def solve(self):
        self.solver.solve()
        u, p = self.soln.subfunctions
        u.rename("Velocity")
        p.rename("Pressure")
        if True:
            file = File("results/flow.pvd")
            file.write(u, p)

        with CheckpointFile('Velocity_field'+'.h5', 'w') as afile:
            afile.save_mesh(self.mesh)
            afile.save_function(u)
            afile.save_function(p)

        self.vel = u # velocity field for EchemSolver

class NavierStokesFlowSolver(FlowSolver):

    """ Incompressible Navier-Stokes solver

        For nondimensional Navier-Stokes, pass Reynolds number to fluid_params.
        For dimensional Navier-Stokes, pass density and kinematic viscosity.
    """
    def setup_problem(self):
        u = self.u
        p = self.p
        v = self.v
        q = self.q

        params = self.fluid_params

        # nondimensional Navier-Stokes
        if params.get("Reynolds number"):
            pprint("Using nondimensional Navier-Stokes")
            self.Re = Constant(params["Reynolds number"])
            F = 1.0/self.Re * inner(grad(u), grad(v)) * dx \
                + inner(dot(grad(u), u), v) * dx \
                - p * div(v) * dx \
                + div(u) * q * dx
                #- inner(u, grad(q)) * dx
                #+ dot(grad(p), v) * dx \
            if self.boundary_markers.get("inlet pressure"):
                n = FacetNormal(self.mesh)
                in_id = self.boundary_markers["inlet pressure"]
                p_in = params["inlet pressure"]
                F -= inner(1/self.Re * dot(grad(u), n) - p_in * n, v) * ds(in_id)
            if self.boundary_markers.get("outlet pressure"):
                n = FacetNormal(self.mesh)
                out_id = self.boundary_markers["outlet pressure"]
                p_out = params["outlet pressure"]
                F -= inner(1/self.Re * dot(grad(u), n) - p_out * n, v) * ds(out_id)
        # dimensional Navier-Stokes
        else:
            pprint("Using dimensional Navier-Stokes")
            rho = params["density"]
            if params.get("kinematic viscosity"):
                nu = params["kinematic viscosity"]
            else:
                nu = params["dynamic viscosity"] / rho
            F = nu * inner(grad(u), grad(v)) * dx \
                + inner(dot(grad(u), u), v) * dx \
                - 1.0/rho * p * div(v) * dx \
                + div(u) * q * dx
            if self.boundary_markers.get("inlet pressure"):
                n = FacetNormal(self.mesh)
                in_id = self.boundary_markers["inlet pressure"]
                p_in = params["inlet pressure"]
                F -= inner(nu * dot(grad(u), n) - 1.0/rho * p_in * n, v) * ds(in_id)
            if self.boundary_markers.get("outlet pressure"):
                n = FacetNormal(self.mesh)
                out_id = self.boundary_markers["outlet pressure"]
                p_out = params["outlet pressure"]
                F -= inner(nu * dot(grad(u), n) - 1.0/rho * p_out * n, v) * ds(out_id)

        self.Form = F

    def setup_solver(self, ksp_solver="lu"):
        """ PCD preconditioner for nondimensional Navier-Stokes
        """

        if ksp_solver == "lu":
            appctx = None
            parameters = {"snes_monitor": None,
                          "ksp_type": "preonly",
                          "pc_type": "lu",
                          "pc_factor_mat_solver_type": "mumps",
                          }

        elif ksp_solver == "pcd":
            appctx = {"Re": self.Re, "velocity_space": 0}
            parameters = {"mat_type": "matfree",
                      "snes_monitor": None,
                     "ksp_type": "fgmres",
                     "ksp_gmres_modifiedgramschmidt": None,
                     "ksp_monitor_true_residual": None,
                     "pc_type": "fieldsplit",
                     "pc_fieldsplit_type": "schur",
                     "pc_fieldsplit_schur_fact_type": "lower",
                     "fieldsplit_0_ksp_type": "preonly",
                     "fieldsplit_0_pc_type": "python",
                     "fieldsplit_0_pc_python_type": "firedrake.AssembledPC",
                     "fieldsplit_0_assembled_pc_type": "hypre",
                     "fieldsplit_1_ksp_type": "gmres",
                     "fieldsplit_1_ksp_rtol": 1e-4,
                     "fieldsplit_1_pc_type": "python",
                     "fieldsplit_1_pc_python_type": "firedrake.PCDPC",
                     "fieldsplit_1_pcd_Mp_ksp_type": "preonly",
                     "fieldsplit_1_pcd_Mp_pc_type": "bjacobi",
                     "fieldsplit_1_pcd_Mp_sub_pc_type": "ilu",
                     "fieldsplit_1_pcd_Kp_ksp_type": "preonly",
                     "fieldsplit_1_pcd_Kp_pc_type": "hypre",
                     "fieldsplit_1_pcd_Fp_mat_type": "matfree"}

        self.problem = NonlinearVariationalProblem(self.Form, self.soln,
                                                   bcs=self.bcs)
        self.solver = NonlinearVariationalSolver(self.problem, appctx=appctx,
                                                 solver_parameters=parameters,
                                                 nullspace=self.nullspace)

class NavierStokesBrinkmanFlowSolver(FlowSolver):

    """ Incompressible Navier-Stokes-Brinkman solver

        For dimensionless form, pass Reynolds number to fluid_params.
        For dimensional form, pass density and kinematic viscosity.
    """
    def setup_problem(self):
        u = self.u
        p = self.p
        v = self.v
        q = self.q
        Z = self.Z

        params = self.fluid_params

        # nondimensional Navier-Stokes-Brinkman
        if params.get("Reynolds number"):
            pprint("Using nondimensional Navier-Stokes-Brinkman")
            self.Re = params["Reynolds number"]
            Da = params["Darcy number"]
            F = 1.0/ self.Re * inner(grad(u), grad(v)) * dx \
                + inner(dot(grad(u), u), v) * dx \
                - p * div(v) * dx \
                + 1.0 / self.Re / Da * inner(u, v) * dx \
                + div(u) * q * dx
            if self.boundary_markers.get("inlet pressure"):
                n = FacetNormal(self.mesh)
                in_id = self.boundary_markers["inlet pressure"]
                p_in = params["inlet pressure"]
                F -= inner(1/self.Re * dot(grad(u), n) - p_in * n, v) * ds(in_id)
            if self.boundary_markers.get("outlet pressure"):
                n = FacetNormal(self.mesh)
                out_id = self.boundary_markers["outlet pressure"]
                p_out = params["outlet pressure"]
                F -= inner(1/self.Re * dot(grad(u), n) - p_out * n, v) * ds(out_id)
        # dimensional Navier-Stokes-Brinkman
        else:
            pprint("Using dimensional Navier-Stokes-Brinkman")
            rho = params["density"]
            if params.get("kinematic viscosity"):
                nu = params["kinematic viscosity"]
            else:
                nu = params["dynamic viscosity"] / rho
            # inverse permeability: scalar field only for now
            if params.get("permeability"):
                inv_K = 1 / params["permeability"]
            else:
                inv_K = params["inverse permeability"]
            if params.get("effective kinematic viscosity"):
                nu_eff = params.get("effective kinematic viscosity")
            else:
                nu_eff = nu
            F = nu_eff * inner(grad(u), grad(v)) * dx \
                + inner(dot(grad(u), u), v) * dx \
                - 1.0/rho * p * div(v) * dx \
                + nu * inv_K * inner(u, v) * dx \
                + div(u) * q * dx
                # - inner(u, grad(q)) * dx # does this sets u.n = 0 on Neumann BCs?
            if self.boundary_markers.get("inlet pressure"):
                n = FacetNormal(self.mesh)
                in_id = self.boundary_markers["inlet pressure"]
                p_in = params["inlet pressure"]
                F -= inner(nu_eff * dot(grad(u), n) - 1.0/rho * p_in * n, v) * ds(in_id)
            if self.boundary_markers.get("outlet pressure"):
                n = FacetNormal(self.mesh)
                out_id = self.boundary_markers["outlet pressure"]
                p_out = params["outlet pressure"]
                F -= inner(nu_eff * dot(grad(u), n) - 1.0/rho * p_out * n, v) * ds(out_id)

        self.Form = F
