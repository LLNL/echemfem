from abc import ABC
from firedrake import *

class FlowSolver(ABC):
    """Base class for a flow solver.

    This class is used to create solvers for fluid flow decoupled from the chemistry
    """

    def __init__(self, mesh, fluid_params, boundary_markers):
        self.mesh = mesh
        self.fluid_params = fluid_params
        self.boundary_markers = boundary_markers
        self.setup_functions()
        self.setup_problem()

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

    def setup_solver(self):
        # set up nonlinear solver - preconditioner is problem specific
        pass

    def solve(self):
        self.solver.solve()
        u, p = self.soln.subfunctions
        u.rename("Velocity")
        p.rename("Pressure")
        if True:
            file = File("flow.pvd")
            file.write(u, p)

        with CheckpointFile('Velocity_field'+'.h5', 'w') as afile:
            afile.save_mesh(self.mesh)
            afile.save_function(u)
            afile.save_function(p)

        self.vel = u # velocity field for EchemSolver

class NavierStokesFlowSolver(FlowSolver):

    """
        Inputs, which are useful for nondim
        Length L: default takes x direction from mesh
        Flow velocity U: reference velocity. if there is an inflow velocity, can use that. if it's pressure BC, no idea
        Density. 
        Kinematic viscosity nu: default room temp water? 1e-6 m2/s
    """
    def setup_problem(self):
        u = self.u
        p = self.p
        v = self.v
        q = self.q
        Z = self.Z

        Re = Constant(100.0)

        F = (
        1.0/Re * inner(grad(u), grad(v)) * dx +
        inner(dot(grad(u), u), v) * dx -
        p * div(v) * dx +
        div(u) * q * dx
        )
        x, y = SpatialCoordinate(self.mesh)

        # setup bcs
        bcs = [
           DirichletBC(Z.sub(0), Constant((0, 0)), (11,10)), # wall
           DirichletBC(Z.sub(0), as_vector([y, Constant(0)]), (12,13,14)), # shear flow
           DirichletBC(Z.sub(1),Constant(0),(13)) # outlet     
           ]

        self.nullspace = MixedVectorSpaceBasis(
        Z, [Z.sub(0), VectorSpaceBasis(constant=True)])

        self.Form = F
        self.bcs = bcs
        self.Re = Re

    def setup_solver(self):
        """ PCD preconditioner
        """

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



mesh=Mesh('../examples/squares_small.msh')
boundary_markers = None
solver = NavierStokesFlowSolver(mesh, None, boundary_markers)
solver.setup_solver()
solver.solve()
