import bpy
import mathutils
import numpy as np
import sys


def make_cell(name: str, location: tuple, radius: float, size: float, scale: tuple):
    """Core function: creates a Blender mesh corresponding to a Goo :class:`Cell` object. 

    :param `Cell` cell: The Goo :class:`Cell` object. 
    :returns: None
    """

    matg = bpy.data.materials.new("Green")
    matg.diffuse_color = (0,0.1,0,0.8)

    #invariant parameters
    subdiv = 2
    vertex_mass = 0.3
    air_damping = 10

    # Making cell
    print('Making cell')


    bpy.ops.mesh.primitive_ico_sphere_add(enter_editmode=False,
                                            align='WORLD',
                                            location=location,
                                            scale=scale,
                                            radius=radius)


    # Give the Blender object the cell's name
    bpy.ops.object.select = True
    obj = bpy.context.active_object
    obj.name = name
    bpy.context.view_layer.objects.active = bpy.data.objects[name]

    # Smooth the mesh
    bpy.ops.object.select = True
    bpy.ops.object.shade_smooth()

    # Add subsurface modifier to make smoother
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = subdiv

    # Add cloth settings for physics
    bpy.ops.object.modifier_add(type='CLOTH')
    bpy.context.object.modifiers["Cloth"].settings.bending_model = 'LINEAR'
    bpy.context.object.modifiers["Cloth"].settings.quality = 5
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1
    bpy.context.object.modifiers["Cloth"].settings.mass = vertex_mass
    bpy.context.object.modifiers["Cloth"].settings.air_damping = air_damping
    bpy.context.object.modifiers["Cloth"].settings.tension_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.compression_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.shear_stiffness = 5
    bpy.context.object.modifiers["Cloth"].settings.bending_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.tension_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.compression_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.shear_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.bending_damping = 0.5
    bpy.context.object.modifiers["Cloth"].settings.use_pressure = True
    bpy.context.object.modifiers["Cloth"].settings.uniform_pressure_force = 2.8
    bpy.context.object.modifiers["Cloth"].settings.pressure_factor = 1
    bpy.context.object.modifiers["Cloth"].settings.fluid_density = 1
    bpy.context.object.modifiers["Cloth"].collision_settings.use_collision = True
    bpy.context.object.modifiers["Cloth"].collision_settings.distance_min = 0.015
    bpy.context.object.modifiers["Cloth"].collision_settings.impulse_clamp = 0

    # add Collision modifier for physics
    bpy.ops.object.modifier_add(type='COLLISION')
    bpy.context.object.collision.use_culling = False
    # bpy.context.object.collision.damping = 0.579821
    # bpy.context.object.collision.thickness_outer = 0.02
    # bpy.context.object.collision.thickness_inner = 0.2
    # bpy.context.object.collision.cloth_friction = 5
    # bpy.ops.object.forcefield_toggle()
    # bpy.context.object.field.type = 'FORCE'
    # bpy.context.object.field.strength = -600
    # bpy.context.object.field.strength = 0
    # bpy.context.object.field.shape = 'POINT'
    # bpy.context.object.name = cell.name

    # add material to cell based on name of material
    obj = bpy.context.active_object
    obj.active_material = matg


make_cell(name='Cell1', location=(0,0,0), radius=1, size=1, scale=(1, 1, 1))