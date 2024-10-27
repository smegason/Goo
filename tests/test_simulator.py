import pytest
import goo
import bpy


@pytest.fixture
def setup_blender():
    bpy.ops.wm.read_factory_settings(use_empty=True)  # Reset to empty scene
    cellsA = goo.CellType("A")
    cellsA.homo_adhesion_strength = 100
    cell = cellsA.create_cell(
        "A1", (0, 0, 0), color=(0.5, 0, 0), target_volume=50, size=2
    )

    molA = goo.Molecule("molA", conc=5, D=0.5, gradient="constant")
    molB = goo.Molecule("molB", conc=10, D=1, gradient="random")
    molC = goo.Molecule("molC", conc=1, D=0, gradient="linear")

    diffusionsystem = goo.DiffusionSystem(molecules=[molA, molB, molC])

    sim = goo.Simulator(
        celltypes=[cellsA],
        diffsystems=[diffusionsystem],
        time=50,
        physics_dt=1,
        molecular_dt=0.1,  # default is physics_dt / 10
    )
    sim.set_seed(2024)
    sim.add_handlers(
        [
            goo.GrowthPIDHandler(),
            goo.RecenterHandler(),
            goo.DiffusionHandler(),
        ]
    )

    yield cell, sim, molA, molB, molC, diffusionsystem


def test_toggle_gravity(setup_blender):
    cell, sim, molA, molB, molC, diffsys = setup_blender
    sim.toggle_gravity(False)

    assert bpy.context.scene.use_gravity is False


def test_add_celltype(setup_blender):
    cell, sim, molA, molB, molC, diffsys = setup_blender

    assert len(sim.celltypes) == 1

    sim.celltypes.append(goo.CellType("B"))
    assert len(sim.celltypes) == 2


def test_get_cells(setup_blender):
    cell, sim, molA, molB, molC, diffsys = setup_blender

    assert sim.get_cells() == [cell]


def test_add_handler(setup_blender):
    cell, sim, molA, molB, molC, diffsys = setup_blender

    handler = goo.SizeDivisionHandler(goo.BisectDivisionLogic, mu=60, sigma=2)
    sim.add_handler(handler)
    bpy.context.scene.frame_set(bpy.context.scene.frame_current + 1)

    assert handler.run in bpy.app.handlers.frame_change_post


def test_run_simulation_background(setup_blender):
    cell, sim, molA, molB, molC, diffsys = setup_blender

    sim.run(end=10)
    assert bpy.context.scene.frame_current == 10
