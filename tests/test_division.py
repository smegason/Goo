import pytest
import bpy
from goo import Cell, CellType
from goo.force import * 
from goo.division import BisectDivisionLogic, BooleanDivisionLogic
import bmesh
from mathutils import Vector


@pytest.fixture
def setup_blender():
    # Create a Blender object to pass to the Cell constructor
    bpy.ops.wm.read_factory_settings(use_empty=True)  # Reset to empty scene
    cellsA = CellType("A")
    cellsA.homo_adhesion_strength = 1000
    cell_bisect_large = cellsA.create_cell("A1", (10, 0, 0), color=(0.5, 0, 0))
    cell_bisect_small = cellsA.create_cell("A2", (10, 0, 0), color=(0.5, 0.5, 0))
    oriented_cell = cellsA.create_cell("A3", (0, 0, 0), color=(0.5, 0.5, 0), scale=(2, 1.7, 1.7))
    
    yield cell_bisect_large, cell_bisect_small, oriented_cell


def test_divide_bisect_cell_volume(setup_blender: tuple[Cell, Cell]):
    cell_bisect_large, cell_bisect_small, _ = setup_blender
    mother_volume_large = cell_bisect_large.volume()
    mother_volume_small = cell_bisect_small.volume()

    logic_large = BisectDivisionLogic(margin=0.025)
    logic_small = BisectDivisionLogic(margin=0.25)
    mother_large, daughter_large = cell_bisect_large.divide(logic_large)
    mother_small, daughter_small = cell_bisect_small.divide(logic_small)
    logic_large.flush()
    logic_small.flush()

    # bpy.context.scene.frame_set(bpy.context.scene.frame_current + 1)
    
    assert isinstance(mother_large, Cell)
    assert isinstance(daughter_large, Cell)
    assert isinstance(mother_small, Cell)
    assert isinstance(daughter_small, Cell)

    assert mother_large.name == "A1.0"
    assert daughter_large.name == "A1.1"
    assert mother_small.name == "A2.0"
    assert daughter_small.name == "A2.1"

    assert mother_large.volume() > mother_small.volume()
    assert daughter_large.volume() > daughter_small.volume()

    assert mother_volume_large / 2 > daughter_large.volume()
    assert mother_volume_small / 2 > daughter_small.volume()

    assert (mother_large.volume() + daughter_large.volume()
            == pytest.approx(mother_volume_large, rel=2e-1))

    assert (mother_small.volume() + daughter_small.volume()
            == pytest.approx(mother_volume_small, rel=3e-1))


def test_divide_oriented_cell(setup_blender: Cell): 
    _, _, oriented_cell = setup_blender
    # mother_volume = oriented_cell.volume()

    com = oriented_cell.COM(local_coords=True)
    axis = oriented_cell.major_axis().axis(local_coords=True)

    assert pytest.approx(com, abs=1e-3) == (0, 0, 0)
    assert pytest.approx(axis, abs=1e-2) == (-1, 0, 0)

    logic = BisectDivisionLogic(margin=0.025)
    daughter1, daughter2 = oriented_cell.divide(logic)
    logic.flush()

    assert (pytest.approx(daughter1.volume(), abs=1e-3) 
            == pytest.approx(daughter2.volume(), abs=1e-3))
