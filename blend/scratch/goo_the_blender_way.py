<<<<<<< HEAD
# create cells in the blender way without using Python classes
import bpy
import sys

# data initialization
# TODO: add functions to add data to dict

#immutable
cell_rendering_data = {
    'ID': 0,
    'enter_editmode': False,
    'align': 'WORLD', # used
    'size': (1, 1, 1),
    'scale': (1, 1, 1), # used
    'arcdiv': 8,
    'subdiv': 2, # used
    'vertex_mass': 0.3, # used
    'density': 1.0,
    'edges_pull': 0.0,
    'edges_bend': 0.2,
    'air_damping': 10, # used
    'self_collision': True,
    'self_collision_stiffness': 0.01,
    'ball_size': 10,
    'softbody_goal': False,
    'sotbody_friction': 0.2,
    'phi': 0,
    'theta': 0,
    'div_axis': (0, 0, 0),
    'mother': 'none',
    'daughters': ['none', 'none'],
    }


# collection of available materials and their node diffusion values
materials = {
    'Green' : (0,0.1,0,0.8), 
    'red' : (0.1,0,0,0.8), 
    'blue' : (0.0,0,1,0.8)
    }

def create_collection(collection_name):
    # create collection
    collection = bpy.data.collections.new(collection_name)
    # link the collection to the scene for visualization 
    bpy.context.scene.collection.children.link(collection)

    return collection

def link_obj2collection(cell, collection): 
    # set the cell to link as active
    bpy.context.view_layer.objects.active = cell
    # Remove object from all collections not used in a scenel; otherwise object will link to object inside the collection
    bpy.ops.collection.objects_remove_all()
    # link object to collection
    collection.objects.link(cell)

def get_vertices(obj): 
    local_vertices = obj.data.vertices
    global_vertices = [(obj.matrix_world @ v.co) for v in obj.data.vertices]
    return local_vertices, global_vertices 

def set_vertices(obj): 
    local_v, global_v = get_vertices(obj)
    obj['vertices_local'] = local_v
    obj['vertices_global'] = global_v


# create cell-like objects
def make_cell_blender(name, radius, location, material):
    # add mesh type as an cell argument? 
    """Core function: creates a Blender mesh corresponding to a Goo :class:`Cell` object. 

    :param str mesh_type: The type of the mesh used to create the cell in Blender. 
    :param str force_name: The name of the force to link to the mesh object. 
    :returns: None
    """
    
    mesh = bpy.ops.mesh.primitive_ico_sphere_add(enter_editmode=cell_rendering_data.get('enter_editmode'),
                                                    align=cell_rendering_data.get('align'),
                                                    location=location,
                                                    scale=cell_rendering_data.get('scale'),
                                                    radius=radius)

    # Give the Blender object the cell's name
    obj = bpy.context.object
    bpy.context.active_object.name = name

    # Add radius as a custom properties of the mesh object
    obj['radius'] = float(radius)


    # Add the coord of each vertex of the mesh as a custom properties
    '''
    # local vertices
    local_vertices = [obj.data.vertices]

    global_vertices = [(obj.matrix_world @ v.co) for v in obj.data.vertices]
    vertices = [obj.data.vertices[vertex].co for vertex in global_vertices]
    '''
    obj['vertices_local'] = obj.data.vertices
    obj['vertices_global'] = [(obj.matrix_world @ v.co) for v in obj.data.vertices]
    obj['vertices_coord_test'] = obj['vertices_global'][0][0]

=======
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
>>>>>>> 784994e237308e3e1c5832067ba898aaf8a42ec5

    # Smooth the mesh
    bpy.ops.object.select = True
    bpy.ops.object.shade_smooth()

    # Add subsurface modifier to make smoother
    bpy.ops.object.modifier_add(type='SUBSURF')
<<<<<<< HEAD
    bpy.context.object.modifiers["Subdivision"].levels = cell_rendering_data.get('subdiv')
=======
    bpy.context.object.modifiers["Subdivision"].levels = subdiv
>>>>>>> 784994e237308e3e1c5832067ba898aaf8a42ec5

    # Add cloth settings for physics
    bpy.ops.object.modifier_add(type='CLOTH')
    bpy.context.object.modifiers["Cloth"].settings.bending_model = 'LINEAR'
    bpy.context.object.modifiers["Cloth"].settings.quality = 5
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1
<<<<<<< HEAD
    bpy.context.object.modifiers["Cloth"].settings.mass = cell_rendering_data.get('vertex_mass')
    obj.modifiers["Cloth"].settings.air_damping = cell_rendering_data.get('air_damping')
=======
    bpy.context.object.modifiers["Cloth"].settings.mass = vertex_mass
    bpy.context.object.modifiers["Cloth"].settings.air_damping = air_damping
>>>>>>> 784994e237308e3e1c5832067ba898aaf8a42ec5
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
<<<<<<< HEAD
    bpy.context.object.collision.use_culling = True
    bpy.context.object.collision.damping = 0.579821
    bpy.context.object.collision.thickness_outer = 0.02
    bpy.context.object.collision.thickness_inner = 0.2
    bpy.context.object.collision.cloth_friction = 5

    # add material to cell based on name of material
    mat = bpy.data.materials.new(material)
    mat.diffuse_color = materials.get(material)
    # the created cell is the active object
    obj = bpy.context.active_object
    obj.active_material = mat

    return obj
    

'''
def data_handler(self, scence, depsgraph): 
    for collection in bpy.data.collections:
                if 'cells' in collection.name_full:
                    for cell in collection.objects:
                        cell['vertices_global'] = [(cell.matrix_world @ v.co) for v in cell.data.vertices]
                        cell['vertices_coord_test'] = cell['vertices_global'][0][0]
'''
=======
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
>>>>>>> 784994e237308e3e1c5832067ba898aaf8a42ec5
