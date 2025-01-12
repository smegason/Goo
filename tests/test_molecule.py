import pytest
import goo
import bpy


@pytest.fixture
def setup_blender():
    bpy.ops.wm.read_factory_settings(use_empty=True)  # Reset to empty scene
    cellsA = goo.CellType("A")
    cellsA.homo_adhesion_strength = 100
    cell = cellsA.create_cell(
        "A1", (0, 0, 0), color=(0.5, 0, 0), size=2, target_volume=50
    )

    molA = goo.Molecule("molA", conc=1, D=1, gradient="linear")
    molB = goo.Molecule("molB", conc=0.5, D=2, gradient="random")
    molC = goo.Molecule("molB", conc=0.5, D=2)

    diffusionsystem = goo.DiffusionSystem(molecules=[molA, molB, molC])

    sim = goo.Simulator(celltypes=[cellsA], diffsystems=[diffusionsystem], time=20)
    sim.set_seed(2024)
    sim.add_handlers(
        [
            goo.GrowthPIDHandler(),
            goo.RecenterHandler(),
            goo.DiffusionHandler(),
        ]
    )

    yield cell, molA, molB, molC


def test_molecule_name(setup_blender):
    molA, molB = setup_blender[1:3]

    assert molA.name == "molA"
    assert molB.name == "molB"

    molA.name = "molA_updated"
    molB.name = "molB_updated"

    assert molA.name == "molA_updated"
    assert molB.name == "molB_updated"


def test_molecule_initial_conc(setup_blender):
    molA, molB = setup_blender[1:3]

    assert molA.conc == 1
    assert molB.conc == 0.5

    molA.conc = 2
    molB.conc = 1

    assert molA.conc == 2
    assert molB.conc == 1


def test_molecule_diff_coeff(setup_blender):
    molA, molB = setup_blender[1:3]

    assert molA.D == 1
    assert molB.D == 2

    molA.D = 2
    molB.D = 1

    assert molA.D == 2
    assert molB.D == 1


def test_molecule_gradient(setup_blender):
    molA, molB, molC = setup_blender[1:4]

    assert molA.gradient == "linear"
    assert molB.gradient == "random"
    if molC.gradient:
        assert False

    molA.D = "random"
    molB.D = "linear"
    molC.D = "random"

    assert molA.D == "random"
    assert molB.D == "linear"
    assert molC.D == "random"
