import pytest
import goo
import bpy
import numpy as np


@pytest.fixture
def setup_blender():
    bpy.ops.wm.read_factory_settings(use_empty=True)  # Reset to empty scene
    cellsA = goo.CellType("A")
    cellsA.homo_adhesion_strength = 100
    cell = cellsA.create_cell("A1", (0, 0, 0), color=(0.5, 0, 0), size=2)

    molA = goo.Molecule("molA", conc=1, D=1, gradient="linear")
    molB = goo.Molecule("molB", conc=1, D=2, gradient="random")
    molC = goo.Molecule("molB", conc=1, D=2)

    diffusionsystem = goo.DiffusionSystem(molecules=[molA, molB, molC])

    yield cell, molA, molB, molC, diffusionsystem


def test_grid_initialization(setup_blender):
    molA, molB, molC, diffsys = setup_blender[1:5]
    kd_tree = diffsys._build_kdtree()
    expected_shape = (len(diffsys.molecules), *diffsys._grid_size)
    closest_point_index = diffsys._get_nearest_grid_index(point=(0, 0, 0))
    closest_point_coord = diffsys._get_grid_position(closest_point_index)

    assert diffsys._grid_concentrations.shape == expected_shape
    assert diffsys._kd_tree == kd_tree
    assert closest_point_index == 61224
    # center of the grid is (25, 25, 25) for a 50x50x50 grid
    assert closest_point_coord == pytest.approx((25, 25, 25), abs=1e1)
    assert kd_tree.n == 125000  # 50*50*50


def test_update_conc(setup_blender):
    molA, molB, molC, diffsys = setup_blender[1:5]
    kd_tree = diffsys._build_kdtree()
    previous_conc = diffsys.get_concentration(mol_idx=0, index=61224)
    diffsys.update_concentration(mol_idx=0, index=61224, value=10)
    updated_conc = diffsys.get_concentration(mol_idx=0, index=61224)

    assert diffsys._kd_tree == kd_tree
    assert updated_conc == previous_conc + 10
