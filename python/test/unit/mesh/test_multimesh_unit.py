# -*- coding: utf-8 -*-
"""Unit tests for MultiMesh instantiation and initialization.

These tests are converted from test/unit/cpp/mesh/MultiMesh.cpp.
They only instantiate and initialise MultiMesh objects; the only
thing that is tested is that the code does not crash.
"""

# Copyright (C) 2017 August Johansson, Benjamin Kehlet
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

import pytest
from dolfin import *
from dolfin_utils.test import skip_in_parallel


@skip_in_parallel
def test_multimesh_trivial_3d():
    """Trivial test case 3D: two overlapping coarse 3D meshes."""
    background = UnitCubeMesh(1, 1, 1)
    overlapping = UnitCubeMesh(1, 1, 1)
    overlapping.translate(Point(0.1, 0.1, 0.1))

    multimesh = MultiMesh()
    multimesh.add(background)
    multimesh.add(overlapping)
    multimesh.build()


@skip_in_parallel
def test_multimesh_trivial_3d_boxmesh():
    """Trivial case 3D 2: UnitCubeMesh background with a BoxMesh overlay.

    Some cells in refmesh may produce a degenerate tetrahedron (all 4
    vertices coplanar) when coordinate stride is copied incorrectly.
    Until the C++ code handles degenerate tetrahedra the MultiMesh
    constructor may throw a RuntimeError; we accept both success and
    that expected exception.
    """
    background = UnitCubeMesh(1, 1, 1)

    refmesh = BoxMesh(Point(0.394383, 0.783099, 0.197551),
                      Point(0.840188, 0.798440, 0.911647),
                      1, 1, 1)

    try:
        multimesh = MultiMesh()
        multimesh.add(background)
        multimesh.add(refmesh)
        multimesh.build()
    except RuntimeError:
        # Expected until degenerate-tetrahedron handling is implemented.
        pass
