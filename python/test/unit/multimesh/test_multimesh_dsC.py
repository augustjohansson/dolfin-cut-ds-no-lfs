import pytest
from dolfin import *
import ufl_legacy as ufl
import numpy as np

from dolfin_utils.test import fixture, skip_in_parallel


@fixture
def multimesh():
    mesh_0 = RectangleMesh(Point(0, 0), Point(1, 1), 10, 10)
    x0 = 0.1 * np.pi
    mesh_1 = RectangleMesh(Point(x0, x0), Point(x0 + 1.0, x0 + 1.0), 10, 10)
    multimesh = MultiMesh()
    multimesh.add(mesh_0)
    multimesh.add(mesh_1)
    multimesh.build()
    return multimesh


@fixture
def V(multimesh):
    element = FiniteElement("Lagrange", triangle, 1)
    return MultiMeshFunctionSpace(multimesh, element)


@fixture
def V_high(multimesh):
    element = FiniteElement("Lagrange", triangle, 3)
    return MultiMeshFunctionSpace(multimesh, element)


@fixture
def V_vector(multimesh):
    return MultiMeshVectorFunctionSpace(multimesh, "CG", 1)


@fixture
def e_x():
    return Expression(("1.0", "0.0"), degree=1)


@fixture
def f_vec2():
    return Expression(("sin(x[1])", "1.0"), degree=3)


@fixture
def v_vec(V_vector):
    return MultiMeshFunction(V_vector)


@fixture
def f_vec():
    return Expression(("x[0]", "x[1]"), degree=2)


@fixture
def f():
    return Expression("sin(pi*x[0])*cos(2*pi*x[1])", degree=4)


@fixture
def f_2():
    return Expression("x[0]*x[1]", degree=2)


@fixture
def v(V):
    return MultiMeshFunction(V)


@fixture
def v_high(V_high):
    return MultiMeshFunction(V_high)


@skip_in_parallel
def test_measure_mul(v, multimesh):
    assert isinstance(v * dX, ufl.form.Form)


@skip_in_parallel
def test_assemble_zero(v, multimesh):
    assert np.isclose(assemble_multimesh(v * dX), 0)


@skip_in_parallel
def test_assemble_area(v, multimesh):
    mesh_0 = multimesh.part(0)
    mesh_1 = multimesh.part(1)
    xmin_0 = mesh_0.coordinates()[:, 0].min()
    ymin_0 = mesh_0.coordinates()[:, 1].min()
    xmax_0 = mesh_0.coordinates()[:, 0].max()
    ymax_0 = mesh_0.coordinates()[:, 1].max()
    xmin_1 = mesh_1.coordinates()[:, 0].min()
    ymin_1 = mesh_1.coordinates()[:, 1].min()
    xmax_1 = mesh_1.coordinates()[:, 0].max()
    ymax_1 = mesh_1.coordinates()[:, 1].max()
    A = (xmax_0 - xmin_0) * (ymax_0 - ymin_0) + (xmax_1 - xmin_1) * (ymax_1 - ymin_1)
    overlap = (xmax_0 - xmin_1) * (ymax_0 - ymin_1)
    v.vector()[:] = 1
    assert np.isclose(assemble_multimesh(v * dX), A - overlap)


# @skip_in_parallel
# def test_assemble_exterior_facet(v, multimesh):
#     mesh_0 = multimesh.part(0)
#     mesh_1 = multimesh.part(1)
#     circumference = 0.0
#     for mesh in (mesh_0, mesh_1):
#         xmin = mesh.coordinates()[:, 0].min()
#         xmax = mesh.coordinates()[:, 0].max()
#         ymin = mesh.coordinates()[:, 1].min()
#         ymax = mesh.coordinates()[:, 1].max()
#         circumference += 2 * (xmax - xmin) + 2 * (ymax - ymin)
#     v.vector()[:] = 1
#     assert np.isclose(assemble_multimesh(v * ds), circumference)


@skip_in_parallel
def test_assemble_dI(v, multimesh):
    mesh_0 = multimesh.part(0)
    mesh_1 = multimesh.part(1)
    xmin = mesh_1.coordinates()[:, 0].min()
    ymin = mesh_1.coordinates()[:, 1].min()
    xmax = mesh_0.coordinates()[:, 0].max()
    ymax = mesh_0.coordinates()[:, 1].max()
    interface_length = (xmax - xmin) + (ymax - ymin)
    v.vector()[:] = 1
    assert np.isclose(assemble_multimesh(v * dI), interface_length)


@skip_in_parallel
def test_assemble_dsC(v, multimesh):
    mesh_0 = multimesh.part(0)
    mesh_1 = multimesh.part(1)
    xmin_0 = mesh_0.coordinates()[:, 0].min()
    ymin_0 = mesh_0.coordinates()[:, 1].min()
    xmax_0 = mesh_0.coordinates()[:, 0].max()
    ymax_0 = mesh_0.coordinates()[:, 1].max()
    xmin_1 = mesh_1.coordinates()[:, 0].min()
    ymin_1 = mesh_1.coordinates()[:, 1].min()
    xmax_1 = mesh_1.coordinates()[:, 0].max()
    ymax_1 = mesh_1.coordinates()[:, 1].max()
    dsC_val = 2 * (xmax_0 - xmin_0) + 2 * (ymax_0 - ymin_0)
    dsC_val += 2 * (xmax_1 - xmax_0) + 2 * (ymax_1 - ymax_0)
    v.vector()[:] = 1
    assert np.isclose(assemble_multimesh(v * dsC), dsC_val)
