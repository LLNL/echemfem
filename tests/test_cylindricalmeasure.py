import pytest
from firedrake import *
from echemfem import CylindricalMeasure

mesh = UnitSquareMesh(2, 2)
x, y = SpatialCoordinate(mesh)
V = FunctionSpace(mesh, "DP", 0)
f0 = Function(V).interpolate(conditional(Or(And(x < .4,
                                                And(y > .4, y < .6)),
                                            And(And(x > .6, x < .9),
                                                And(y > .6, y < .9))), 1., 0.))
f1 = Function(V).interpolate(conditional(Or(And(x < .4,
                                                And(y > .4, y < .6)),
                                            And(And(x > .6, x < .9),
                                                And(y > .6, y < .9))), 0., 1.))
VS = FunctionSpace(mesh, "HDiv Trace", 0)
f0_S = Function(VS).interpolate(conditional(Or(And(x < .4,
                                               And(y > .4, y < .6)),
                                               And(And(x > .6, x < .9),
                                               And(y > .6, y < .9))), 1., 0.))
f1_S = Function(VS).interpolate(conditional(Or(And(x < .4,
                                               And(y > .4, y < .6)),
                                               And(And(x > .6, x < .9),
                                               And(y > .6, y < .9))), 0., 1.))
mesh = RelabeledMesh(mesh,
                     [f0, f1, f0_S, f1_S],
                     [1, 2, 1, 2])
Vc = FunctionSpace(mesh, "CG", 1)
u = Function(Vc).assign(Constant(1))
r, z = SpatialCoordinate(mesh)
_dx = CylindricalMeasure(r, 'cell')
_ds = CylindricalMeasure(r, 'exterior_facet')
_dS = CylindricalMeasure(r, 'interior_facet')


def test_dx():
    a = assemble(u * r * dx)
    ar = assemble(u * _dx)
    assert (abs(a-ar) < 1e-15)


def test_dx_subs():
    a = assemble(u * r * dx((1, 2,)))
    ar = assemble(u * _dx((1, 2,)))
    assert (abs(a-ar) < 1e-15)


def test_ds():
    a = assemble(u * r * ds)
    ar = assemble(u * _ds)
    assert (abs(a-ar) < 1e-15)


def test_ds_subs():
    a = assemble(u * r * ds((2, 3, 4,)))
    ar = assemble(u * _ds((2, 3, 4,)))
    assert (abs(a-ar) < 1e-15)


def test_dS():
    a = assemble(u * r * dS)
    ar = assemble(u * _dS)
    assert (abs(a-ar) < 1e-15)


def test_dS_subs():
    a = assemble(u * r * dS((1, 2,)))
    ar = assemble(u * _dS((1, 2,)))
    assert (abs(a-ar) < 1e-15)
