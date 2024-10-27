import pytest
import bpy
import goo
from goo import Cell, CellType
from goo.force import *
from goo.handler import RecenterHandler, GrowthPIDHandler, RandomMotionHandler
from goo.simulator import Simulator
import bmesh
from mathutils import Vector


@pytest.fixture
def setup_blender():
    # Create a Blender object to pass to the Cell constructor
    bpy.ops.wm.read_factory_settings(use_empty=True)  # Reset to empty scene
    cellsA = CellType("A")
    cellsA.homo_adhesion_strength = 1000
    cell = cellsA.create_cell("A1", (0, 0, 0), color=(0.5, 0, 0), size=2)

    sim = Simulator([cellsA], time=5, physics_dt=1)
    sim.setup_world()
    sim.add_handlers(
        [
            GrowthPIDHandler(target_volume=25),
            RecenterHandler(),
        ]
    )
    yield cell
