# Navier-Stokes equations
# =======================
#

from firedrake import *
import sys
'''

N = 20
L=5
#M = RectangleMesh(N*L, N,L,0.2)

Ly = 0.1
Lx = 1.
#print(" Hauteur domaine=",pow(peclet,-1./3))
##mesh = RectangleBoundaryLayerMesh(50,50+int(10*pow(peclet,-1./3)),Lx,2*pow(peclet,-1./3)+Ly,50,1e-1, Ly_bdlayer = 5e-3, boundary=(3,1,))
##mesh = RectangleBoundaryLayerMesh(50,50,Lx,Ly,50,1e-1, Ly_bdlayer = 5e-3, boundary=(3,1,))
#M=RectangleMesh(50,50,Lx,Ly,quadrilateral=True)'''




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

# Having set up the problem, we now move on to solving it.  Some
# preconditioners, for example pressure convection-diffusion (PCD), require
# information about the the problem that is not easily accessible from
# the bilinear form.  In the case of PCD, we need the Reynolds number
# and additionally which part of the mixed velocity-pressure space the
# velocity corresponds to.  We provide this information to
# preconditioners by passing in a dictionary context to the solver.
# This is propagated down through the matrix-free operators and is
# therefore accessible to custom preconditioners. ::

    appctx = {"Re": Re, "velocity_space": 0}

# Now we'll solve the problem.  First, using a direct solver.  Again, if
# MUMPS is not installed, this solve will not work, so we wrap the solve
# in a ``try/except`` block. ::

    from firedrake.petsc import PETSc


# Now we'll show an example using the :class:`~.PCDPC` preconditioner
# that implements the pressure convection-diffusion approximation to the
# pressure Schur complement.  We'll need more solver parameters this
# time, so again we'll set those up in a dictionary. ::

    parameters = {"mat_type": "matfree",
              "snes_monitor": None,

# We'll use a non-stationary Krylov solve for the Schur complement, so
# we need to use a flexible Krylov method on the outside. ::

             "ksp_type": "fgmres",
             "ksp_gmres_modifiedgramschmidt": None,
             "ksp_monitor_true_residual": None,

# Now to configure the preconditioner::

             "pc_type": "fieldsplit",
             "pc_fieldsplit_type": "schur",
             "pc_fieldsplit_schur_fact_type": "lower",

# we invert the velocity block with LU::

             "fieldsplit_0_ksp_type": "preonly",
             "fieldsplit_0_pc_type": "python",
             "fieldsplit_0_pc_python_type": "firedrake.AssembledPC",
             "fieldsplit_0_assembled_pc_type": "lu",

# and invert the schur complement inexactly using GMRES, preconditioned
# with PCD. ::

             "fieldsplit_1_ksp_type": "gmres",
             "fieldsplit_1_ksp_rtol": 1e-4,
             "fieldsplit_1_pc_type": "python",
             "fieldsplit_1_pc_python_type": "firedrake.PCDPC",

# We now need to configure the mass and stiffness solvers in the PCD
# preconditioner.  For this example, we will just invert them with LU,
# although of course we can use a scalable method if we wish. First the
# mass solve::

             "fieldsplit_1_pcd_Mp_ksp_type": "preonly",
             "fieldsplit_1_pcd_Mp_pc_type": "lu",

# and the stiffness solve.::

             "fieldsplit_1_pcd_Kp_ksp_type": "preonly",
             "fieldsplit_1_pcd_Kp_pc_type": "lu",

# Finally, we just need to decide whether to apply the action of the
# pressure-space convection-diffusion operator with an assembled matrix
# or matrix free.  Here we will use matrix-free::

             "fieldsplit_1_pcd_Fp_mat_type": "matfree"}

# With the parameters set up, we can solve the problem, remembering to
# pass in the application context so that the PCD preconditioner can
# find the Reynolds number. ::

    up.assign(0)

    solve(F == 0, up, bcs=bcs, nullspace=nullspace, solver_parameters=parameters,
      appctx=appctx)

# And finally we write the results to a file for visualisation. ::
#print(up.name)
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
