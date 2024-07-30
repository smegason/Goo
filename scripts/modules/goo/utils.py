from functools import reduce
from typing import Union

import bpy
import bmesh
from bpy.types import Modifier
from mathutils import *


class BlenderObject:
    def __init__(self, obj: bpy.types.Object):
        self.obj = obj

    @property
    def name(self):
        return self.obj.name

    @name.setter
    def name(self, name):
        self.obj.name = name

    @property
    def loc(self):
        return self.obj.location

    @loc.setter
    def loc(self, loc):
        self.obj.location = loc

    def hide(self):
        self.obj.hide_set(True)

    # ----- CUSTOM PROPERTIES -----
    def __setitem__(self, k: str, v: Union[float, list[float], int, list[int], str]):
        self.obj[k] = v

    def __contains__(self, k: str):
        return k in self.obj.keys()

    def __getitem__(self, k):
        return self.obj[k]


class Axis:
    def __init__(self, axis, start, end, world_matrix):
        """
        :param axis: Vector of axis
        :param start: Vector of start endpoint in object space
        :param end: Vector of end endpoint in object space
        :param world_matrix: 4x4 matrix of object to world transformation
        :param local_coords: if coordinates given are in local space
        """
        self._axis = axis
        self._start = start
        self._end = end
        self._matrix_world = world_matrix

    def axis(self, local_coords=False):
        if local_coords:
            axis = self._axis.copy()
            axis.rotate(self._matrix_world.to_quaternion().inverted())
            return axis
        return self._axis

    def endpoints(self, local_coords=False):
        mat = self._matrix_world.inverted() if local_coords else Matrix.Identity(4)
        return [mat @ self._start, mat @ self._end]

    def length(self, local_coords=False):
        start, end = self.endpoints(local_coords)
        return (end - start).length


# ----- BLENDER FUNCTIONS -----
def create_mesh(
    name,
    loc,
    mesh="icosphere",
    size=1.5,
    rotation=(0, 0, 0),
    scale=(1, 1, 1),
    subdivisions=2,
    **kwargs,
):
    bm = bmesh.new()

    if isinstance(rotation, tuple):
        rotation = Euler(rotation)
    elif isinstance(rotation, Quaternion):
        rotation = rotation.to_euler()

    match mesh:
        case "icosphere":
            bmesh.ops.create_icosphere(
                bm, subdivisions=subdivisions, radius=size, **kwargs
            )
            # bmesh.ops.beautify_fill(bm, faces=bm.faces, edges=bm.edges)
            bmesh.ops.triangulate(bm, faces=bm.faces[:])
        case "plane":
            bmesh.ops.create_grid(
                bm, x_segments=1, y_segments=1, size=size / 2, **kwargs
            )
        case "monkey":
            bmesh.ops.create_monkey(bm, **kwargs)
        case "cube":
            bmesh.ops.create_cube(bm, size=size, **kwargs)
            # bmesh.ops.beautify_fill(bm, faces=bm.faces, edges=bm.edges)
            # bmesh.ops.triangulate(bm, faces=bm.faces[:])
        case _:
            raise ValueError(
                """mesh must be one of "icosphere", "plane", "monkey" or "cube"."""
            )

    me = bpy.data.meshes.new(f"{name}_mesh")
    bm.to_mesh(me)
    bm.free()

    obj = bpy.data.objects.new(name, me)
    obj.location = loc
    obj.rotation_euler = rotation
    obj.scale = scale

    return obj


def create_material(name, color):
    mat = bpy.data.materials.new(name=name)
    r, g, b = color
    mat.diffuse_color = (r, g, b, 0.8)  # viewport color
    mat.use_nodes = True
    mat.blend_method = "BLEND"

    # get the material nodes
    nodes = mat.node_tree.nodes
    nodes.clear()

    # create principled node for main color
    node_main = nodes.new(type="ShaderNodeBsdfPrincipled")
    node_main.location = -200, 100
    node_main.inputs["Base Color"].default_value = (r, g, b, 0.8)
    node_main.inputs["Metallic"].default_value = 0.036
    node_main.inputs["Roughness"].default_value = 0.318
    node_main.inputs["IOR"].default_value = 1.450

    # specular
    node_main.inputs["Anisotropic"].default_value = 0.041
    node_main.inputs["Anisotropic Rotation"].default_value = 0.048
    node_main.inputs["Alpha"].default_value = 0.414

    # create noise texture source
    node_noise = nodes.new(type="ShaderNodeTexNoise")
    node_noise.inputs["Scale"].default_value = 0.600
    node_noise.inputs["Detail"].default_value = 15.0
    node_noise.inputs["Roughness"].default_value = 0.500
    node_noise.inputs["Distortion"].default_value = 3.0

    # create HSV
    node_HSV = nodes.new(type="ShaderNodeHueSaturation")
    node_HSV.inputs["Hue"].default_value = 0.800
    node_HSV.inputs["Saturation"].default_value = 2.00
    node_HSV.inputs["Value"].default_value = 2.00
    node_HSV.inputs["Fac"].default_value = 1.00

    # create second principled node for random color variation
    node_random = nodes.new(type="ShaderNodeBsdfPrincipled")
    node_random.location = -200, -100
    node_random.inputs["Base Color"].default_value = (r, g, b, 1)
    node_random.inputs["Metallic"].default_value = 0.0

    node_random.inputs["Roughness"].default_value = 0.482
    node_random.inputs["Anisotropic"].default_value = 0.0
    node_random.inputs["Anisotropic Rotation"].default_value = 0.0
    node_random.inputs["IOR"].default_value = 1.450
    node_random.inputs["Alpha"].default_value = 0.555

    # create mix shader node
    node_mix = nodes.new(type="ShaderNodeMixShader")
    node_mix.location = 0, 0
    node_mix.inputs["Fac"].default_value = 0.079

    # create output node
    node_output = nodes.new(type="ShaderNodeOutputMaterial")
    node_output.location = 200, 0

    # link nodes
    links = mat.node_tree.links
    links.new(node_noise.outputs[1], node_HSV.inputs[4])  # link_noise_HSV
    links.new(node_HSV.outputs[0], node_random.inputs[0])  # link_HSV_random
    links.new(node_main.outputs[0], node_mix.inputs[1])  # link_main_mix
    links.new(node_random.outputs[0], node_mix.inputs[2])  # link_random_mix
    links.new(node_mix.outputs[0], node_output.inputs[0])  # link_mix_out

    return mat


# --- PHYSICS MODIFIER CONSTRUCTORS ---
class PhysicsConstructor:
    def __init__(self, *mod_contructors: "ModConstructor"):
        self.mod_constructors = [*mod_contructors]

    def __call__(self, bobj: BlenderObject):
        for mod_constructor in self.mod_constructors:
            mod_constructor().construct(bobj.obj)


class ModConstructor:
    name = ""
    type = ""

    def construct(self, obj):
        mod = obj.modifiers.new(name=self.name, type=self.type)
        self.setup_mod(mod)

    def setup_mod(self, mod: bpy.types.Modifier):
        pass


class ClothConstructor(ModConstructor):
    name = "Cloth"
    type = "CLOTH"

    def setup_mod(self, mod: bpy.types.ClothModifier):
        stiffness = 1
        pressure = 0.01

        mod.settings.quality = 10
        mod.settings.air_damping = 10
        mod.settings.bending_model = "ANGULAR"
        mod.settings.mass = 1
        mod.settings.time_scale = 1
        mod.point_cache.frame_start = bpy.context.scene.frame_start
        mod.point_cache.frame_end = bpy.context.scene.frame_end
        # Cloth > Stiffness
        mod.settings.tension_stiffness = stiffness
        mod.settings.compression_stiffness = stiffness
        mod.settings.shear_stiffness = stiffness
        mod.settings.bending_stiffness = 1
        # Cloth > Damping
        mod.settings.tension_damping = 50
        mod.settings.compression_damping = 50
        mod.settings.shear_damping = 50
        mod.settings.bending_damping = 0.5
        # Cloth > Pressure
        mod.settings.use_pressure = True
        mod.settings.uniform_pressure_force = pressure
        mod.settings.use_pressure_volume = True
        mod.settings.target_volume = 1
        mod.settings.pressure_factor = 2
        mod.settings.fluid_density = 1.05
        # Cloth > Collisions
        mod.collision_settings.collision_quality = 4
        mod.collision_settings.use_collision = True
        mod.collision_settings.use_self_collision = True
        mod.collision_settings.self_friction = 0
        mod.collision_settings.friction = 0
        mod.collision_settings.self_distance_min = 0.01
        mod.collision_settings.distance_min = 0.01
        mod.collision_settings.self_impulse_clamp = 100
        mod.collision_settings.impulse_clamp = 100

        # Cloth > Field Weights
        mod.settings.effector_weights.gravity = 0


class YolkClothConstructor(ClothConstructor):
    def setup_mod(self, mod: bpy.types.ClothModifier):
        stiffness = 5
        pressure = 5.2

        mod.settings.quality = 10
        mod.settings.air_damping = 10
        mod.settings.bending_model = "ANGULAR"
        mod.settings.mass = 2
        mod.settings.time_scale = 1
        # Cloth > Stiffness
        mod.settings.tension_stiffness = stiffness
        mod.settings.compression_stiffness = stiffness
        mod.settings.shear_stiffness = stiffness
        mod.settings.bending_stiffness = stiffness
        # Cloth > Damping
        mod.settings.tension_damping = 50
        mod.settings.compression_damping = 50
        mod.settings.shear_damping = 50
        mod.settings.bending_damping = 0.5
        # Cloth > Internal Springs
        mod.settings.use_internal_springs = True
        mod.settings.internal_spring_max_length = 1
        mod.settings.internal_spring_max_diversion = 0.785398
        mod.settings.internal_spring_normal_check = False
        mod.settings.internal_tension_stiffness = 10000
        mod.settings.internal_compression_stiffness = 10000
        mod.settings.internal_tension_stiffness_max = 10000
        mod.settings.internal_compression_stiffness_max = 10000
        # Cloth > Pressure
        mod.settings.use_pressure = True
        mod.settings.uniform_pressure_force = pressure
        mod.settings.use_pressure_volume = True
        mod.settings.target_volume = 1
        mod.settings.pressure_factor = 2
        mod.settings.fluid_density = 1.15
        # Cloth > Collisions
        mod.collision_settings.collision_quality = 10
        mod.collision_settings.use_collision = True
        mod.collision_settings.use_self_collision = True
        mod.collision_settings.self_friction = 0
        mod.collision_settings.friction = 0
        mod.collision_settings.self_distance_min = 0.005
        mod.collision_settings.distance_min = 0.005
        mod.collision_settings.self_impulse_clamp = 0


class CollisionConstructor(ModConstructor):
    name = "Collision"
    type = "COLLISION"

    def setup_mod(self, mod: bpy.types.CollisionModifier):
        mod.settings.use_culling = True
        mod.settings.damping = 1
        mod.settings.thickness_outer = 0.025
        mod.settings.thickness_inner = 0.25
        mod.settings.cloth_friction = 0
        mod.settings.use_normal = False


class BoundaryCollisionConstructor(CollisionConstructor):
    def setup_mod(self, mod: bpy.types.CollisionModifier):
        super(BoundaryCollisionConstructor, self).setup_mod(mod)
        mod.settings.use_culling = False


class SubsurfConstructor(ModConstructor):
    name = "Subdivision"
    type = "SUBSURF"

    def setup_mod(self, mod: bpy.types.SubsurfModifier):
        mod.subdivision_type = "CATMULL_CLARK"
        mod.levels = 1
        mod.render_levels = 1


# not stable
class RemeshConstructor(ModConstructor):
    name = "Remesh"
    type = "REMESH"

    def setup_mod(self, mod: bpy.types.RemeshModifier):
        mod.mode = "VOXEL"
        mod.voxel_size = 0.5
        mod.adaptivity = 0
        mod.use_remove_disconnected = True
        mod.use_smooth_shade = True
        mod.show_in_editmode = True
