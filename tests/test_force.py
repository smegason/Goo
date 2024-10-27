import pytest
import bpy
from mathutils import Vector, Euler
from goo.force import (
    Force,
    AdhesionForce,
    MotionForce,
    create_force,
    create_adhesion,
    ForceCollection,
)


@pytest.fixture
def setup_blender():
    bpy.ops.wm.read_factory_settings(use_empty=True)  # Reset to empty scene
    adhesion_force = create_adhesion(
        strength=100,
        name="AdhesionForce",
        loc=(0, 0, 0),
        shape="POINT",
    )

    yield adhesion_force


def test_force_creation(setup_blender):
    force = setup_blender
    assert force.name == "AdhesionForce"
    assert force.strength == 100
    assert force.loc == Vector((0, 0, 0))
    assert force.shape == "POINT"


def test_enable_disable_force(setup_blender):
    force = setup_blender
    force.disable()
    assert not force.enabled()
    force.enable()
    assert force.enabled()


def test_force_shape(setup_blender):
    force = setup_blender
    force.shape = "POINT"
    assert force.shape == "POINT"


def test_force_strength(setup_blender):
    force = setup_blender
    assert force.strength == 100


def test_motion_force_set_loc(setup_blender):
    obj = bpy.data.objects.new("MotionForce", None)
    motion_force = MotionForce(obj)
    target_loc = Vector((2, 2, 2))
    motion_force.point_towards(target_loc)
    print(motion_force.obj.rotation_euler)
    assert (
        motion_force.obj.rotation_euler == target_loc.to_track_quat("Z", "X").to_euler()
    )


def test_force_collection(setup_blender):
    force_collection = ForceCollection("TestCollection")
    assert force_collection.name == "TestCollection"
    assert len(force_collection.forces) == 0


def test_add_remove_force_collection(setup_blender):
    force_collection = ForceCollection("TestCollection")
    obj = bpy.data.objects.new("TestForce", None)
    force = Force(obj, "FORCE")
    force_collection.add(force)
    assert len(force_collection.forces) == 1
    force_collection.remove(force)
    assert len(force_collection.forces) == 0


def test_global_force_collection(setup_blender):
    global_forces = ForceCollection.global_forces()
    assert global_forces.name == "globals"
