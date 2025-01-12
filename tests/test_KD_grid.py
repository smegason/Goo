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

    diffusionsystem = goo.DiffusionSystem(
        molecules=[molA, molB, molC],
        grid_size=(50, 50, 50),
        grid_center=(0, 0, 0),
    )

    yield cell, molA, molB, molC, diffusionsystem


def test_grid_initialization(setup_blender):
    molA, molB, molC, diffsys = setup_blender[1:5]
    kd_tree = diffsys.build_kdtree()
    expected_shape = diffsys.grid_size
    closest_point_index = diffsys._nearest_idx(point=(1, 1, 1))
    closest_point_coord = kd_tree.data[closest_point_index]

    for mol in diffsys.molecules:
        assert diffsys._grid_concentrations[mol].shape == expected_shape
    assert diffsys._kd_tree == kd_tree

    # check point in grid
    assert closest_point_coord == pytest.approx((1, 1, 1), abs=1e1)
    assert kd_tree.n == 125000  # 50*50*50


def test_update_conc(setup_blender):
    molA, molB, molC, diffsys = setup_blender[1:5]
    kd_tree = diffsys.build_kdtree()
    previous_conc = diffsys.get_concentration(molA, (1, 1, 1))
    diffsys.update_concentration(molA, (1, 1, 1), value=10)
    updated_conc = diffsys.get_concentration(molA, (1, 1, 1))

    assert diffsys._kd_tree == kd_tree
    assert updated_conc == previous_conc + 10
