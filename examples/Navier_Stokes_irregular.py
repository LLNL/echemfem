# Navier-Stokes equations
# =======================

# TODO: These are nondimensionalized equations! We need the dimensional velocity

from firedrake import *
import sys
from firedrake.petsc import PETSc


def navier_stokes_no_slip(mesh):
    M=mesh
    V = VectorFunctionSpace(M, "CG", 2)
    W = FunctionSpace(M, "CG", 1)
    Z = V * W
    n = FacetNormal(M)
    up = Function(Z,name="up")
    u, p = split(up)
    v, q = TestFunctions(Z)


    Re = Constant(100.0)

    F = (
    1.0/Re * inner(grad(u), grad(v)) * dx +
    inner(dot(grad(u), u), v) * dx -
    p * div(v) * dx +
    div(u) * q * dx
    )
    x, y = SpatialCoordinate(M)

    bcs = [
       DirichletBC(Z.sub(0), Constant((0, 0)), (11,10)),
       DirichletBC(Z.sub(0), as_vector([y, Constant(0)]), (12,13,14)),
       DirichletBC(Z.sub(1),Constant(0),(13))      
       ]

    nullspace = MixedVectorSpaceBasis(
    Z, [Z.sub(0), VectorSpaceBasis(constant=True)])

    appctx = {"Re": Re, "velocity_space": 0}
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
             "fieldsplit_0_assembled_pc_type": "lu",
             "fieldsplit_1_ksp_type": "gmres",
             "fieldsplit_1_ksp_rtol": 1e-4,
             "fieldsplit_1_pc_type": "python",
             "fieldsplit_1_pc_python_type": "firedrake.PCDPC",
             "fieldsplit_1_pcd_Mp_ksp_type": "preonly",
             "fieldsplit_1_pcd_Mp_pc_type": "lu",
             "fieldsplit_1_pcd_Kp_ksp_type": "preonly",
             "fieldsplit_1_pcd_Kp_pc_type": "lu",
             "fieldsplit_1_pcd_Fp_mat_type": "matfree"}

    up.assign(0)

    solve(F == 0, up, bcs=bcs, nullspace=nullspace, solver_parameters=parameters,
      appctx=appctx)

    u, p = up.subfunctions
    u.rename("Velocity")
    p.rename("Pressure")

    with CheckpointFile('Velocity_field'+'.h5', 'w') as afile:
        afile.save_mesh(M)  # optional
        afile.save_function(u)
        afile.save_function(p)
    return u

    

#navier_stokes_no_slip()
# A runnable python script implementing this demo file is available
# :demo:`here <navier_stokes.py>`.
