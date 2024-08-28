import pytest
import bpy
import goo
from goo import Cell, CellType
from goo.force import * 
from goo.handler import AdhesionLocationHandler, GrowthPIDHandler
from goo.simulator import Simulator
import bmesh
from mathutils import Vector


@pytest.fixture
def setup_blender():
    # Create a Blender object to pass to the Cell constructor
    bpy.ops.wm.read_factory_settings(use_empty=True)  # Reset to empty scene
    cellsA = CellType("A")
    cellsA.homo_adhesion_strength = 500
    cell = cellsA.create_cell("A1", (0, 0, 0), color=(0.5, 0, 0), size=1)

    sim = Simulator([cellsA], time=5, physics_dt=1)
    # cannot use setup_world in testing because requires Blender in non-headless mode
    # sim.setup_world()
    sim.add_handlers(
        [
            GrowthPIDHandler(target_volume=50),
            AdhesionLocationHandler(),
        ]
    )
    yield cell


def test_force_follows_cell(setup_blender: Cell):
    cell = setup_blender
    assert cell.loc == cell.adhesion_forces[0].loc


def test_force_follows_cell_t1(setup_blender):
    cell = setup_blender
    bpy.context.scene.frame_set(bpy.context.scene.frame_current + 1)

    assert cell.loc == cell.adhesion_forces[0].loc


def test_force_follows_cell_t2(setup_blender):
    cell = setup_blender
    bpy.context.scene.frame_set(bpy.context.scene.frame_current + 1)

    assert cell.loc == cell.adhesion_forces[0].loc


def test_force_follows_cell_t3(setup_blender):
    cell = setup_blender
    bpy.context.scene.frame_set(bpy.context.scene.frame_current + 1)

    assert cell.loc == cell.adhesion_forces[0].loc
