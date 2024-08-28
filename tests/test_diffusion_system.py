import pytest
import goo
import bpy


@pytest.fixture
def setup_blender():
    bpy.ops.wm.read_factory_settings(use_empty=True)  # Reset to empty scene
    cellsA = goo.CellType("A")
    cellsA.homo_adhesion_strength = 100
    cell = cellsA.create_cell("A1", (0, 0, 0), color=(0.5, 0, 0), size=2)

    molA = goo.Molecule("molA", conc=5, D=0.5, gradient="constant")
    molB = goo.Molecule("molB", conc=10, D=1, gradient="random")
    molC = goo.Molecule("molC", conc=1, D=0, gradient="linear")

    diffusionsystem = goo.DiffusionSystem(molecules=[molA, molB, molC])
    
    sim = goo.Simulator(
        celltypes=[cellsA], 
        diffsystems=[diffusionsystem], 
        time=50, 
        physics_dt=1, 
        molecular_dt=0.1  # default is physics_dt / 10 
    )
    sim.set_seed(2024)
    sim.add_handlers(
        
        [
            goo.GrowthPIDHandler(target_volume=50),
            goo.AdhesionLocationHandler(),
            goo.DiffusionHandler(diffusionSystem=diffusionsystem),
        ]
    )

    yield cell, sim, molA, molB, molC, diffusionsystem


def test_diffsys_molecules(setup_blender):
    molA, molB, molC, diffsys = setup_blender[2:6]
    assert diffsys.molecules == [molA, molB, molC]


def test_diffsys_gridsize(setup_blender):
    diffsys = setup_blender[5]
    assert diffsys.grid_size == (50, 50, 50)  # default value


def test_diffsys_timestep(setup_blender):
    cell, sim, molA, molB, molC, diffsys = setup_blender
    assert diffsys.time_step == sim.physics_dt / 10


def test_diffsys_totaltime(setup_blender):
    cell, sim, molA, molB, molC, diffsys = setup_blender
    assert diffsys.total_time == sim.physics_dt


def test_cell_intial_concentrations(setup_blender):
    cell, sim, molA, molB, molC, diffsys = setup_blender

    bpy.context.scene.frame_set(bpy.context.scene.frame_current + 1)

    assert cell.molecules_conc.get(molA.name) == 5.0
    assert cell.molecules_conc.get(molB.name) == pytest.approx(8.287, abs=1e-3)
    assert cell.molecules_conc.get(molC.name) == pytest.approx(0.347, abs=1e-3)


def test_cell_diffused_concentrations(setup_blender):
    cell, sim, molA, molB, molC, diffsys = setup_blender

    for i in range(5):
        bpy.context.scene.frame_set(bpy.context.scene.frame_current + 1)

    # constant gradient leads to no diffusion
    assert cell.molecules_conc.get(molA.name) == 5.0
    assert cell.molecules_conc.get(molB.name) == pytest.approx(10.747, abs=1e-3)
    assert cell.molecules_conc.get(molC.name) == pytest.approx(0.408, abs=1e-3)