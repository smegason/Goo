from collections import defaultdict
import bpy
import mathutils
from mathutils import Vector, Matrix
import numpy as np
import sys
import random
import os 
import math
from datetime import datetime
import json
import subprocess
import bmesh


"""
Refactored by Antoine Ruzette, November 2022.

Note on comments: 
Core functions are required for ``Goo`` to correctly run. 
Auxilliary function are not required for ``Goo`` to run 
but are used for supporting tasks such as retrieving 
metrics - typically they can be removed. 

Sphynx docstring was enhanced by adding types, and written in a more sphynx-ic way. 
"""


def calculate_volume(obj):
    """Core function: calculates the volume of the Blender mesh. 

    In order to retrieve the mesh as it is currrently evaluated - including 
    the effect of all modifiers - in the simulation, its corresponding evaluated 
    ID is obtained from the dependency graph for the current context. 

    .. seealso:: [Blender API Documentation > 
                ``evaluated_get(depsgraph)``]
                (https://docs.blender.org/api/current/bpy.types.ID.html?highlight=evaluated_get#bpy.types.ID.evaluated_get)

    :param bpy.data.objects['name'] obj: The Blender mesh.
    :returns: The volume of the mesh. 
    :rtype: float

    .. note:: The function may return a negative volume. 
                See ``calc_volume(signed=True)``. 

    """

    # We need to get the cell as it is evaluated in the simulation.
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    # Use BMesh to calculate volume of a mesh 
    mesh_from_eval = obj_eval.to_mesh()
    bm = bmesh.new()
    bm.from_mesh(mesh_from_eval)
    # Apply the object's scale and dimensions to the bmesh
    bm.transform(obj_eval.matrix_world)
    # Calculate volume
    volume = bm.calc_volume()
    # Output the result
    print(f"Volume of {obj.name}: {abs(volume)}")
    # Free the bmesh
    bm.free()

    return volume


def get_centerofmass(obj): 
    """Core function: calculates the center of mass of a mesh. 

    This function fetch the evaluated object's dependency graph, 
    retrieves the coordinates of each vertex then computes the center of mass 
    as the mean position among the set of vertices.  

    :param bpy.data.objects['name'] obj: The Blender mesh.
    :returns: The coordinates of the center of mass of the mesh as a tuple(x, y, z). 
    :rtype: tuple
    """

    bpy.context.view_layer.objects.active = obj
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    vertices = obj_eval.data.vertices
    vert_coords = np.asarray([(obj_eval.matrix_world @ v.co) for v in vertices])
    COM = np.mean(vert_coords, axis=0)
 
    return tuple(COM)


def get_long_axis_np(obj):
    """Calculates the long axis of a mesh.

    This function calculates the first eigenvector of the vertices in the mesh,
    which corresponds to the long axis.

    :param bpy.data.objects['name'] obj: The Blender mesh.
    :returns: The coordinates of the long axis of the mesh as a tuple(x, y, z)
              which gives direction from the origin (0, 0, 0).
    :rtype: tuple
    """
    # Get the evaluated object and its vertices
    dg = bpy.context.evaluated_depsgraph_get()
    evaluated_object = obj.evaluated_get(dg)
    vertices = evaluated_object.data.vertices
    vertex_coords = np.array([evaluated_object.matrix_world @ v.co for v in vertices])

    # Subtract the mean from each dimension
    vertex_coords -= np.mean(vertex_coords, axis=0)

    # Calculate the covariance matrix and its eigenvectors
    cov_matrix = np.cov(vertex_coords, rowvar=False)
    eigenvals, eigenvecs = np.linalg.eigh(cov_matrix)
    sort_indices = np.argsort(eigenvals)

    # Get the eigenvector corresponding to the largest eigenvalue
    long_axis = eigenvecs[:, sort_indices[-1]]

    return tuple(long_axis)[:3]


def get_long_axis_global(obj):

    # Get the evaluated object and its vertices
    dg = bpy.context.evaluated_depsgraph_get()
    evaluated_object = obj.evaluated_get(dg)

    vertices = evaluated_object.data.vertices
    vertex_coords = np.array([evaluated_object.matrix_world @ v.co for v in vertices])
    # Calculate the covariance matrix of the vertices
    covariance_matrix = np.cov(vertex_coords, rowvar=False)
    # Calculate the eigenvectors and eigenvalues of the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
    # Sort the eigenvectors by descending eigenvalues
    eigenvectors = eigenvectors[:, eigenvalues.argsort()[::-1]]
    # The major axis is the first eigenvector
    major_axis = eigenvectors[:, 0]

    # Convert the major axis to global coordinates, no need for this line
    major_axis_global = np.array(evaluated_object.matrix_world) @ np.append(major_axis, 
                                                                            1)

    return tuple(major_axis_global)[:3]


def get_long_axis(obj):
    """Core function: calculates the long axis of a mesh. 

    This function calculates the first eigenvector of the vertices in the mesh, 
    which corresponds to the long axis.

    :param bpy.data.objects['name'] obj: The Blender mesh.
    :returns: The coordinates of the long axis of the mesh as a tuple(x, y, z) 
                which gives direction from the origin (0, 0, 0). 
    :rtype: tuple
    """

    # We need to get the cell as it is evaluated in the simulation.
    # To do this, we fetch its dependency graph and obtain the
    # evaluated cell (denoted as obj_eval here)
    dg = bpy.context.evaluated_depsgraph_get()
    evaluated_object = obj.evaluated_get(dg)
    # We obtain the (x, y, z) coordinates of the vertices in the
    # evaluated cell
    vertices = evaluated_object.data.vertices
    vertex_coords = [(evaluated_object.matrix_world @ v.co) for v in vertices]
    vertex_coords = np.asarray(vertex_coords)

    # We separate the x, y, and z coordinates into their own arrays.
    # We also subtract the mean of each dimension from the corresponding
    # array values (normalization). 
    # Mesh is centered on the origin. This is part of the PCA algorithm.
    x = vertex_coords[:, 0]
    x = x - np.mean(x)
    y = vertex_coords[:, 1]
    y = y - np.mean(y)
    z = vertex_coords[:, 2]
    z = z - np.mean(z)

    # We stack the three arrays together to make the "new" coordinates
    new_coords = np.vstack([x, y, z])
    # This is then used to find the covariance matrix of the coordinates
    cov_matrix = np.cov(new_coords)
    # Per the PCA algorithm, we find the eigenalues and eigenvectors
    # of the covariance matrix
    eigenvals, eigenvecs = np.linalg.eig(cov_matrix)

    # The eigenvalues are sorted, and the primary eigenvector
    # is the major axis.
    sort_indices = np.argsort(eigenvals)
    major_x, major_y, major_z = eigenvecs[:, sort_indices[-1]]
    long_axis = (major_x, major_y, major_z)

    # We project each vertex onto the major axis and calculate
    # the distance between the vertex and the closest and farthest
    # projections. The difference between these distances is the
    # length of the axis contained within the mesh.
    # Calculate the principal components of the mesh
    covariance_matrix = np.cov(vertex_coords, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
    principal_axis = eigenvectors[:, np.argmax(eigenvalues)]
    principal_axis = np.abs(principal_axis)
    long_axis_normalized = tuple(principal_axis / np.linalg.norm(principal_axis))

    # Find the two vertices farthest apart along the principal axis
    vertex_distances = vertex_coords.dot(principal_axis)
    vertex_indices = np.argsort(vertex_distances)
    first_vertex = vertices[vertex_indices[0]]
    last_vertex = vertices[vertex_indices[-1]]

    # Calculate the length of the long axis contained in the mesh
    long_axis_length = np.linalg.norm(
        (evaluated_object.matrix_world @ last_vertex.co) -
        (evaluated_object.matrix_world @ first_vertex.co)
    )

    # Compute end points of the long axis
    endpoints = [first_vertex.co, last_vertex.co]

    return long_axis, long_axis_normalized, long_axis_length, endpoints


def get_longaxis_angles(axis):
    """Auxilliary function: retrieves the 2 angles associated 
        with the long axis of a cell. 

    :param tuple axis: Coordinates of the long axis to the division plane 
                        as retrived by ``get_major_axis(obj)``.
    :returns: 
        - phi (:py:class:`tuple`) - Phi is the angle between the division axis 
            and the z-axis in spherical coordinates. 
        - theta (:py:class:`tuple`) - Theta is the angle projected on the xy plane 
            in cartesian 3D coordinates. 
    """

    # We define the unit vector of z-axis
    z_axis = np.array((0, 0, 1))
    # We find the unit vector of the division axis
    division_axis = axis/np.linalg.norm(axis)
    # We calculate the dot product of the normalized
    # z and division axes
    dot_product = np.dot(z_axis, division_axis)
    # The inverse cosine of the dot product is phi
    phi = np.arccos(dot_product)

    # The first step to find theta is to find the projection
    # of the division axis on the z axis
    proj_axis_on_z = dot_product*z_axis
    # This projection is subtracted from the division axis
    # To find the projection on the xy-plane
    proj_xy_axis = division_axis - proj_axis_on_z
    # We normalize this projection
    proj_xy_axis = proj_xy_axis/np.linalg.norm(proj_xy_axis)
    # We take the dot product of the x-axis and the normalized projection
    dot_product = np.dot((1, 0, 0), proj_xy_axis)
    # The inverse cosin of this dot product is theta.
    theta = np.arccos(dot_product)

    print(f'Division angles:(phi={phi}; theta={theta})')
    return phi, theta


def repair_hole(obj):

    """Core function: repair the holes created by the bissection of the Blender mesh. 
    
    This function adds a new face to the cell 
    after division (bisection in two equal parts) 
    of the mesh and regularizes the mesh so that the mesh is evenly covered with faces. 

    :param bpy.data.objects['name'] obj: The Blender mesh. 
    :returns: None
    """

    # Go into Blender Edit Mode
    bpy.ops.object.mode_set(mode='EDIT')
    # Select all the edges in the mesh
    bpy.ops.mesh.select_mode(type="EDGE")
    bpy.ops.mesh.select_all(action='SELECT')
    # Add a new face (based on open eduges)
    bpy.ops.mesh.edge_face_add()
    # Deselct the edges
    bpy.ops.mesh.select_all(action='DESELECT')
    # Return to Blender Object Mode
    bpy.ops.object.mode_set(mode='OBJECT')
    # Remesh the cell with voxels of size 0.2
    bpy.ops.object.modifier_add(type='REMESH')
    bpy.context.object.modifiers["Remesh"].voxel_size = 0.2
    bpy.ops.object.modifier_apply(modifier="Remesh")


def print_info(obj):
    """Auxilliary function: Retrieves the number of vertices, 
        edges and faces of a Blender object

    :param bpy.data.objects['name'] obj: The Blender mesh. 
    :returns: The number of vertices, faces and edges of ``obj``. 
    :rtype: list of float
    """
    vertices = obj.data.vertices
    edges = obj.data.edges
    faces = obj.data.polygons

    print(f'Vertices: {len(vertices)}; Edges: {len(edges)}; Faces: {len(faces)}')
    
    return [vertices, edges, faces]


# not used
def get_curvature_deviation(obj):
    # Get the mesh data
    mesh = obj.data
    bm = bmesh.new()
    bm.from_mesh(mesh)
    
    # Compute the deviation from a perfect sphere for each face
    deviation_sum = 0
    total_faces = 0
    for face in bm.faces:
        # Compute the center of the face
        center = np.mean([v.co for v in face.verts], axis=0)
        
        # Compute the deviation from a perfect sphere
        deviation = np.linalg.norm(center) - 1
        
        deviation_sum += deviation
        total_faces += 1
    
    # Divide the sum of deviations by the total number of
    # faces to get the average deviation
    avg_deviation = deviation_sum / total_faces
    
    return avg_deviation


def select_translate_cell(obj, COM, major_axis):
    """Core function: selects one of the daughter cells to translate.
    
    After dividing the mother cell, this function compares the
    center of mass of the original mother cell to that of the
    corresponding daughter cells. To do so, the signed distance
    between the daughter cell's center of mass and the original
    plane of division is calculated.
    If this distance is positive, the function returns True
    and this selected daughter cell will be translated later.
    If this distance is negative, the function returns False
    and the other daughter cell will be translated later

    :param bpy.data.objects['name'] obj: The Blender mesh of the mother cell.
    :param tuple COM: The XYZ coordinates of the center of mass of the mother cell.
    :param tuple major_axis: The major axis of the mother mesh
        specified by XYZ coordinates from the origin.
    :returns: A boolean flag indicating which cell has to be translated
        in ``translate_cell(obj, axis)``.
    :rtype: Boolean
    """

    # Get the mother cell's center of mass
    com = np.copy(COM)
    # Set the daughter cell's location to its center of mass,
    # and set this value as the new center of mass (denoted new_COM)
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
    new_COM = obj.location
    # Calculate the signed distance between the new center of mass
    # and the plane of division
    distance = mathutils.geometry.distance_point_to_plane(new_COM, com, major_axis)
    # If the signed distance is greater than 0, return True. Else, return False.
    translated = True if distance > 0 else False
    return translated


def translate_cell(obj, axis):
    """Core function: translates cells along the given axis.

    This function is used during cell division so daughters have some space.

    :param bpy.data.objects['name'] obj: The Blender mesh to translate.
    :param tuple axis: The XYZ coordinates of the major axis of
        the Blender mesh to translate.
    :returns: None
    """

    # TODO: reasons to translate of such a distance??

    # Calculate the magnitude of the given axis
    magnitude = ((axis[0])**2 + (axis[1])**2 + (axis[2])**2)

    # Normalize the axis by dividing eagh term by the magnitude
    # Multiply each term by 0.1 so the cell is only translated a small distance
    # if this is not done, daughter cells are far away from each other
    new_coord = (0.1*axis[0]/magnitude, 0.1*axis[1]/magnitude, 0.1*axis[2]/magnitude)
    # translate the cell
    bpy.ops.transform.translate(value=new_coord)    


def separate_cell(obj):
    """Core function: splits one cell into two for cell division.

    The mother Blender mesh is split in two given a division plane.
    The division plane is specified by a point on this plane, 
    which is the center of mass of the mother, and the direction of this point, 
    given by the long axis of the mother mesh.
    By naming convention, the daughter cells inherit their mother name added with .001.

    :param bpy.data.objects['name'] obj: the Blender mesh to divide in two.
    :returns:
        - mother_name (:py:class:`str`) - The name of the mother cell.
        - daughter_name (:py:class:`str`) - The name of the daughter cell.
        - COM (:py:class:`tuple`) - The XYZ coordinates of the center of mass 
            of the mother cell.
        - major_axis (:py:class:`tuple`) - The XYZ coordinates of the long axis 
            of the mother cell.
    """

    # TODO allow division along a supplied axis.
    # If none supplied then divide along major axis

    # Get the cell as it is evaluated in Blender
    bpy.ops.object.mode_set(mode='OBJECT')
    dg = bpy.context.evaluated_depsgraph_get()
    obj = obj.evaluated_get(dg)
    # The mother cell name is the name of the cell currently being divided
    mother_name = obj.name
    # By Blender convention, duplicates of existing objects are given
    # the same name as the original object, but with ".001" added to the end
    # We use the naming convention to tentatively set the name of one daughter cell
    daughter_name = f"{mother_name}.001"
    # Get the cell's center of mass
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
    COM = obj.location
    # Enter edit mode
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    # Get the cell's major axis
    major_axis, tmp1, tmp2 = get_long_axis(obj)
    # The division plane is specified by a point on this plane (COM) 
    # and by the direction of the points on plane (major axis i.e. normal axis).  
    bpy.ops.mesh.bisect(plane_co=COM, plane_no=major_axis, use_fill=False, flip=False)
    # Enter object mode
    bpy.ops.object.mode_set(mode='OBJECT')
    # Now we separate the split mesh into two separate objects.
    # First, obtain the vertex coordinates of the bisected mother cell.
    obj = bpy.data.objects[mother_name]
    new_vertices = obj.data.vertices
    new_vert_coords = [(obj.matrix_world @ v.co) for v in new_vertices]
    new_vert_coords = np.asarray(new_vert_coords)

    # We will choose half of the vertices to be separated into a new object.
    # We create a list of the vertex indices we will use 
    # to create the daughter cell object
    separation_indices = []
    # We loop through each vertex and calculate the signed distance between the vertex
    # and the division plane (which is specified by the center of mass and major axis).
    # If the distance is greater than -0.05 (the vertex is on a specific
    # half of the cell),
    # we add the vertex's index i to our list of separation indices
    for i in range(len(new_vert_coords)):
        distance = mathutils.geometry.distance_point_to_plane(new_vert_coords[i],
                                                              COM, major_axis)
        if distance > -0.05:
            separation_indices.append(i)
    # Enter edit mode and deselect all vertices
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type="VERT")
    bpy.ops.mesh.select_all(action='DESELECT')
    # Enter object mode and loop through all the vertices.
    # If the vertex's index ix contained in the list of separation indices,
    # we select it for separation
    bpy.ops.object.mode_set(mode='OBJECT')
    for index in separation_indices:
        obj.data.vertices[index].select = True
    # Enter edit mode and separate the selected vertices as a new daughter cell.
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.separate(type='SELECTED')
    # Enter object mode and select only the "original" mother cell
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.data.objects[mother_name].select_set(False)
    bpy.data.objects[daughter_name].select_set(False)
    bpy.data.objects[mother_name].select_set(True)
    return mother_name, daughter_name, COM, major_axis


def get_line_angles(cell, p1, p2):

    # Calculate the direction vectors of the global axes
    x_dir = (Vector((1, 0, 0)) @ cell.matrix_world).normalized()
    y_dir = (Vector((0, 1, 0)) @ cell.matrix_world).normalized()
    z_dir = (Vector((0, 0, 1)) @ cell.matrix_world).normalized()

    # Get the vector representing the line
    line_vector = Vector((p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])).normalized()

    # Calculate the angles between the line and each axis
    # 90 is added to get the orthogonal
    x_angle = math.radians(math.degrees(math.acos(line_vector.dot(x_dir))) + 90)
    y_angle = math.radians(math.degrees(math.acos(line_vector.dot(y_dir))) + 90)
    z_angle = math.radians(math.degrees(math.acos(line_vector.dot(z_dir))))

    return x_angle, y_angle, z_angle


def divide_boolean(obj): 

    # Get COM
    p1 = get_centerofmass(obj)
    # Get the long axis in global coordinate from the origin
    p2 = get_long_axis_global(obj)    

    # Get mother modifiers and their values
    mother_modifiers = obj.modifiers

    # Get the values of the modifiers
    modifier_values = []
    for modifier in mother_modifiers:
        modifier_data = {}
        for attr in dir(modifier):
            if not attr.startswith("__") and not callable(getattr(modifier, attr)):
                modifier_data[attr] = getattr(modifier, attr)
        modifier_values.append(modifier_data)

    # Print the values of the modifiers
    for modifier_value in modifier_values:
        print(modifier_value)

    # Create a dictionary to store the modifier data
    '''modifier_data = {}
    # Loop through each modifier and store its properties in the dictionary
    for mod in modifiers:
        mod_data = {}
        for prop in dir(mod):
            if not prop.startswith("__") and not callable(getattr(mod, prop)):
                mod_data[prop] = getattr(mod, prop)
        modifier_data[mod.name] = mod_data'''

    # Center the long axis on the COM, and normalize vector
    # line = Vector((p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])).normalized()
    # Get the angles representing the orientation of the division plane
    x, y, z = get_line_angles(obj, p1, p2)
    # Get name of dividing cell
    mother_name = obj.name
    # Create division plane
    plane = bpy.ops.mesh.primitive_plane_add(enter_editmode=False, 
                                             size=100, 
                                             align='WORLD', 
                                             location=p1, 
                                             scale=(1, 1, 1), 
                                             rotation=(x, y, z))
    # Set plane as active object
    plane = bpy.context.active_object
    # Rename plane with specific cell name
    plane_name = f"{mother_name}_division plane"
    plane.name = plane_name

    # Add solidify modifier to the plane, add tickness to the plane
    bpy.ops.object.modifier_add(type='SOLIDIFY')
    solid_mod = plane.modifiers[-1]
    solid_mod.offset = 0
    solid_mod.thickness = 0.01
    bpy.ops.object.modifier_apply(modifier=solid_mod.name)

    Force.force_list = [
        force for force in Force.force_list if force.name != obj['adhesion force']
    ]
    strength = bpy.data.objects[obj['adhesion force']].field.strength
    falloff = bpy.data.objects[obj['adhesion force']].field.falloff_power
    force_collection = bpy.data.objects[obj['adhesion force']].users_collection[0].name
    if obj['adhesion force'] in bpy.data.objects:
        bpy.data.objects.remove(bpy.data.objects[obj['adhesion force']], do_unlink=True)

    # Add boolean modifier to the original object
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type='BOOLEAN')
    bool_mod = obj.modifiers[-1]
    bool_mod.operand_type = 'OBJECT'
    bool_mod.object = bpy.data.objects[plane_name]
    bool_mod.operation = 'DIFFERENCE'
    bool_mod.solver = 'EXACT'
    bpy.ops.object.modifier_apply(modifier=bool_mod.name)

    # Hide the plane
    bpy.data.objects[plane_name].hide_set(True)

    # Deselect all vertices in edit mode
    bpy.ops.object.mode_set(mode='EDIT')  
    bpy.ops.mesh.select_mode(type="VERT")
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

    # Select a vertex in object mode
    mesh = obj.data
    vertex = mesh.vertices[0]
    vertex.select = True

    # select all vertices linked to the selected vertex
    bpy.ops.object.mode_set(mode='EDIT') 
    bpy.ops.mesh.select_linked(delimit={'NORMAL'})
    # Separate the outer and inner parts of the mesh
    bpy.ops.mesh.separate(type='SELECTED')
    bpy.ops.object.mode_set(mode='OBJECT')  

    d1 = bpy.context.selected_objects[0]  # selects cell_003, to be renamed
    d1.name = f"{mother_name}.001"  # new cell

    d2 = bpy.context.scene.objects[mother_name]
    d2.name = f"{mother_name}.002"  # mother cell

    '''bpy.context.view_layer.objects.active = d1
    bpy.ops.mesh.primitive_round_cube_add(
        change=True, 
        radius=1, 
        size=(1, 1, 1), 
        arc_div=4, 
        lin_div=0, 
        div_type='CORNERS', 
        odd_axis_align=False, 
        no_limit=False
    )


    # edit to object update
    bpy.ops.object.mode_set(mode = 'EDIT')  
    bpy.ops.mesh.select_mode(type="VERT")
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.mesh.reveal()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.context.view_layer.update()
    
    bpy.context.view_layer.objects.active = d2
    bpy.ops.mesh.primitive_round_cube_add(
        change=True,
        radius=1,
        size=(1, 1, 1),
        arc_div=4,
        lin_div=0,
        div_type='CORNERS',
        odd_axis_align=False,
        no_limit=False
    )
    bpy.context.object.scale = (1,1,1)

    # edit to object update
    bpy.ops.object.mode_set(mode = 'EDIT')  
    bpy.ops.mesh.select_mode(type="VERT")
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.mesh.reveal()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.context.view_layer.update()'''

    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.objects.active = d1
    d1.select_set(True)
    bpy.ops.object.convert(target='MESH')
    bpy.context.view_layer.update()

    '''# Loop through each modifier in the modifier data dictionary 
    # and add it to the first empty object
    for mod_name, mod_data in modifier_data.items():
        mod = d1.modifiers.new(name=mod_name, type=mod_data['type'])
        for prop, value in mod_data.items():
            if hasattr(mod, prop):
                setattr(mod, prop, value)'''

    bpy.ops.object.select_all(action='DESELECT')
    d2.select_set(True)
    bpy.context.view_layer.objects.active = d2
    bpy.ops.object.convert(target='MESH')
    bpy.context.view_layer.update()

    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.update()
    # bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)

    return d1, d2, strength, force_collection, falloff, modifier_values


def apply_physics(obj, modifiers): 

    print('Starting modifiers declaration')

    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)

    # Add subsurface modifier to make smoother
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = 1
    # subdiv_mod = obj.modifiers[-1]
    # bpy.ops.object.modifier_apply(modifier=subdiv_mod.name)

    # Add cloth settings 
    bpy.ops.object.modifier_add(type='CLOTH')
    bpy.context.object.modifiers['Cloth'].settings.quality = 3
    bpy.context.object.modifiers['Cloth'].settings.air_damping = 10
    bpy.context.object.modifiers['Cloth'].settings.bending_model = 'ANGULAR'
    bpy.context.object.modifiers["Cloth"].settings.mass = 3
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1

    # Cloth > Stiffness 
    bpy.context.object.modifiers['Cloth'].settings.tension_stiffness = 1
    bpy.context.object.modifiers['Cloth'].settings.compression_stiffness = 1
    bpy.context.object.modifiers['Cloth'].settings.shear_stiffness = 1
    bpy.context.object.modifiers['Cloth'].settings.bending_stiffness = 1
    # Cloth > Damping
    bpy.context.object.modifiers['Cloth'].settings.tension_damping = 50
    bpy.context.object.modifiers['Cloth'].settings.compression_damping = 50
    bpy.context.object.modifiers['Cloth'].settings.shear_damping = 50
    bpy.context.object.modifiers['Cloth'].settings.bending_damping = 0.5
    # Cloth > Internal Springs
    bpy.context.object.modifiers['Cloth'].settings.use_internal_springs = True
    bpy.context.object.modifiers['Cloth'].settings.internal_spring_max_length = 1
    bpy.context.object.modifiers['Cloth'].settings.internal_spring_max_diversion = \
        0.785398
    bpy.context.object.modifiers['Cloth'].settings.internal_spring_normal_check = True
    bpy.context.object.modifiers['Cloth'].settings.internal_tension_stiffness = 10
    bpy.context.object.modifiers['Cloth'].settings.internal_compression_stiffness = 10
    bpy.context.object.modifiers['Cloth'].settings.internal_tension_stiffness_max = \
        10000
    bpy.context.object.modifiers['Cloth'].settings.internal_compression_stiffness_max \
        = 10000
    # Cloth > Pressure
    bpy.context.object.modifiers["Cloth"].settings.use_pressure = True
    bpy.context.object.modifiers['Cloth'].settings.uniform_pressure_force = 5
    bpy.context.object.modifiers['Cloth'].settings.use_pressure_volume = True
    bpy.context.object.modifiers['Cloth'].settings.target_volume = 0
    bpy.context.object.modifiers['Cloth'].settings.pressure_factor = 1
    bpy.context.object.modifiers['Cloth'].settings.fluid_density = 1.05
    # Cloth > Collisions
    bpy.context.object.modifiers['Cloth'].collision_settings.collision_quality = 3
    bpy.context.object.modifiers['Cloth'].collision_settings.use_collision = True
    bpy.context.object.modifiers['Cloth'].collision_settings.use_self_collision = True
    bpy.context.object.modifiers['Cloth'].collision_settings.self_friction = 50
    bpy.context.object.modifiers['Cloth'].collision_settings.friction = 80
    bpy.context.object.modifiers['Cloth'].collision_settings.self_distance_min = 0.001
    bpy.context.object.modifiers['Cloth'].collision_settings.distance_min = 0.001
    bpy.context.object.modifiers['Cloth'].collision_settings.self_impulse_clamp = 0

    # Collision
    bpy.ops.object.modifier_add(type='COLLISION')
    bpy.context.object.modifiers['Collision'].settings.use_culling = True
    bpy.context.object.modifiers['Collision'].settings.damping = 1
    bpy.context.object.modifiers['Collision'].settings.thickness_outer = 0.02
    bpy.context.object.modifiers['Collision'].settings.thickness_inner = 0.2
    bpy.context.object.modifiers['Collision'].settings.cloth_friction = 80

    # Remesh
    remesh_mod = obj.modifiers.new(name=f"Remesh_{obj.name}", type='REMESH')
    obj.modifiers[f"Remesh_{obj.name}"].name = f"Remesh_{obj.name}"
    remesh_mod.mode = 'SMOOTH'
    remesh_mod.octree_depth = 4
    remesh_mod.scale = 0.8
    remesh_mod.use_remove_disconnected = True
    remesh_mod.use_smooth_shade = False
    remesh_mod.show_in_editmode = True
    bpy.ops.object.modifier_move_to_index(modifier=f"Remesh1_{obj.name}", index=3)

    # Volume conservation
    '''bpy.ops.object.constraint_add(type='MAINTAIN_VOLUME')
    bpy.context.object.constraints["Maintain Volume"].volume = 1
    bpy.context.object.constraints["Maintain Volume"].owner_space = 'WORLD'
    bpy.context.object.constraints["Maintain Volume"].free_axis = 'SAMEVOL_Z'''

    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.update()


def apply_daugther_force(obj, strength, collection, falloff): 

    # strength = bpy.data.objects[d2['force']].field.strength
    # force_collection = bpy.data.objects[d2['force']].users_collection[0].name

    if obj['force'] in bpy.data.objects:
        bpy.data.objects.remove(bpy.data.objects[obj['force']], do_unlink=True)

    # Declare new forces 
    print(obj.name)
    make_force(f"{obj.name}_force", obj.name, strength, falloff, collection)
    # update the force name for its cell
    bpy.data.objects[obj.name]['force'] = f"{obj.name}_force"

    return 


def to_sphere(obj, rate): 
    print(f'{obj.name} under to_sphere')

    obj = bpy.context.active_object

    # hemi_key = obj.shape_key_add(name=f"{obj.name}_hemisphere", from_mix=False)
    bpy.context.object.active_shape_key_index = 0

    sphere_key = obj.shape_key_add(name=f"{obj.name}_sphere", from_mix=False)
    bpy.context.object.active_shape_key_index = 1
    bpy.context.view_layer.update()

    # Update the object and the shape key
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS', center='MEDIAN')
    '''# Calculate the average vertex position
    avg_pos = obj.location

    # Switch to Edit mode and select all vertices
    bpy.ops.object.mode_set(mode='EDIT')
    mesh = obj.data
    
    for vertex in mesh.vertices:
        vertex.select = True

    # Move each vertex to the surface of a sphere with center at the average position
    for vertex in mesh.vertices:
        vec = vertex.co - avg_pos
        vec.normalize()
        vertex.co = avg_pos + vec

    # Deselect all vertices
    for vertex in mesh.vertices:
        vertex.select = False

    # Switch back to Object mode and update the view layer
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.context.view_layer.update()'''

    '''# Get the active object and its mesh
    obj = bpy.context.active_object
    mesh = obj.data
    bpy.ops.object.mode_set(mode='EDIT')

    # Get the selected vertices in edit mode
    bm = bmesh.from_edit_mesh(mesh)
    verts = [v for v in bm.verts if v.select]

    # Calculate the centroid of the selected vertices
    centroid = get_centerofmass(obj)

    # Move each vertex towards the centroid until it lies on
    # the surface of a sphere with the same radius as the original hemisphere
    for v in verts:
        v.co = Vector(centroid) + (v.co - Vector(centroid)).normalized() * 0.1

    # Update the mesh and exit edit mode
    bmesh.update_edit_mesh(mesh)

    bpy.ops.mesh.reveal()
    bpy.ops.object.mode_set(mode = 'OBJECT')'''

    '''bpy.ops.object.mode_set(mode = 'EDIT')
    bpy.ops.mesh.select_mode(type="VERT")
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.reveal()

    bpy.ops.transform.tosphere(value=rate,
                               mirror=True,
                               use_proportional_edit=False,
                               proportional_edit_falloff='SMOOTH',
                               proportional_size=1,
                               use_proportional_connected=False,
                               use_proportional_projected=False)

    bpy.ops.mesh.reveal()
    bpy.ops.object.mode_set(mode = 'OBJECT')'''
    bpy.context.view_layer.update()

    sphere_key.slider_max = 1
    sphere_key.value = 1

    # bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
    bpy.ops.object.shape_key_remove(all=True, apply_mix=True)

    '''obj = bpy.context.active_object
    remesh_mod = obj.modifiers[-1]
    bpy.ops.object.modifier_apply(modifier=remesh_mod.name)'''

    obj = bpy.context.active_object

    # simplify topology 
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type="VERT")
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.decimate(ratio=0.3)
    bpy.ops.mesh.reveal()
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.context.view_layer.update()

    obj.select_set(False)


def remesh(obj): 

    if obj.name in bpy.data.objects:
        bpy.context.view_layer.objects.active = bpy.data.objects[obj.name]
        obj.select_set(True)
        print(obj.name)
        remesh_mod = obj.modifiers[-1]
        bpy.ops.object.modifier_apply(modifier=remesh_mod.name)

        bpy.ops.object.modifier_add(type='REMESH')
        bpy.context.object.modifiers["Remesh"].mode = 'SMOOTH'
        bpy.context.object.modifiers["Remesh"].octree_depth = 3
        bpy.context.object.modifiers["Remesh"].scale = 0.99
        bpy.context.object.modifiers["Remesh"].use_remove_disconnected = True
        bpy.context.object.modifiers["Remesh"].use_smooth_shade = True
        
    return 


def divide_new(obj, indices, line):

    com = get_centerofmass(obj)
    # Get the cell as it is evaluated in Blender
    obj = bpy.context.active_object
    bpy.context.view_layer.objects.active = obj

    # The mother cell name is the name of the cell currently being divided
    mother_name = obj.name
    # daughter_name = f"{mother_name}.001"

    print(f"Bisection normal: {line}")
                                                           
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type="VERT")
    # bpy.context.space_data.overlay.show_wireframes = True
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.bisect(plane_co=com,
                        plane_no=(line[0] + 0.005, line[1] + 0.005, line[2] + 0.005),
                        use_fill=True,
                        flip=False,
                        clear_outer=True
                        )
    
    # Select the active vertex    
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

    # Get the mesh data and select a random vertex
    mesh = obj.data
    vertex = mesh.vertices[0]
    vertex.select = True

    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_linked(delimit={'NORMAL'})
    bpy.ops.mesh.reveal()
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.mode_set(mode='EDIT')
    
    # Separate the outer and inner parts of the mesh
    bpy.ops.mesh.separate(type='SELECTED')
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.mesh.reveal()

    # bpy.ops.mesh.reveal()
    bpy.ops.object.mode_set(mode='OBJECT')  
    bpy.context.view_layer.update()

    # Get the selected objects (i.e., the newly separated meshes)
    selected_objects = bpy.context.selected_objects

    # Assign names to the separated meshes
    if len(selected_objects) == 2:
        print(selected_objects)
        selected_objects[0].name = f"{mother_name}.001"
        selected_objects[1].name = f"{mother_name}.002"
    else:
        print("Error: Expected 2 selected objects, found", len(selected_objects))

    d1 = bpy.data.objects[f"{mother_name}.001"]
    d2 = bpy.data.objects[f"{mother_name}.002"]

    com1 = get_centerofmass(d1)
    com2 = get_centerofmass(d2)

    print(f"New centers of mass: 1/ {com1}, 2/ {com2}")

    # Create forces for daughter cells
    '''make_force(f"{d1.name}_force", f"{d1.name}", -1000, 1, collection_name)
    make_force(f"{d2.name}_force", f"{d2.name}", -1000, 1, collection_name)'''

            
def divide(obj): 
    """Core function: divides a mother Blender mesh into two daughter Blender meshes. 
    
    The division plane is orthogonal to the major axis of the parent mesh. 
    Additionally, this function may translate daughter cells, 
    translating one mesh, filling the two holes, and retriangulating the meshes

    :param bpy.data.objects['mother_name'] obj: The soon-to-be mother Blender mesh. 
    :returns: 
        - daughter1 (bpy.data.objects['daughter1_name']) - The Blender mesh of 
            one of the daughter cells. 
        - daughter2 (bpy.data.objects['daughter2_name']) - The Blender mesh of 
            the other daughter cell. 
    """
    '''# Select mother cell
    obj.select_set(True)
    # Split mother's Blender mesh in two
    m_name, d_name, COM, major_axis = seperate_cell(obj)
    # Decides which cell gets translated
    translated = select_translate_cell(bpy.data.objects[m_name], COM, major_axis)

    if translated == True:
        translate_cell(bpy.data.objects[m_name], major_axis)
        bpy.data.objects[m_name].select_set(False)
    else:
        bpy.data.objects[d_name].select_set(True)
        bpy.data.objects[m_name].select_set(False)
        translate_cell(bpy.data.objects[d_name], major_axis)
        bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
        bpy.data.objects[d_name].select_set(False)

    bpy.data.objects[m_name].select_set(True)
    bpy.context.view_layer.objects.active = bpy.data.objects[m_name]
    repair_hole(bpy.data.objects[m_name])
    bpy.context.object.name = m_name + "0"
    daughter1 = Cell(m_name + "0", loc=bpy.context.object.location)
    daughter1.data['mother'] = m_name
    daughter1.data['daughters'] = ['none', 'none']
    bpy.data.objects[m_name + "0"].select_set(False)
    bpy.data.objects[d_name].select_set(True)
    bpy.context.view_layer.objects.active = bpy.data.objects[d_name]
    repair_hole(bpy.data.objects[d_name])
    bpy.context.object.name = m_name + "1"
    bpy.data.objects[m_name + "1"].select_set(False)
    daughter2 = Cell(bpy.context.object.name, loc=bpy.context.object.location)
    daughter2.data['mother'] = m_name
    daughter2.data['daughters'] = ['none', 'none']
    return daughter1, daughter2'''


def turn_off_physics():
    """Core function: turns physics off for the currently selected Blender object.

    The cloth physics for cells are turned off before division to avoid 
    irregular mesh behavior after. 

    :param: None
    :returns: None
    """
    bpy.ops.object.modifier_remove(modifier="Cloth")


def turn_on_physics():
    """Core function: turns physics on for the currently selected Blender object.
    
    Physics are turned back on after cell division occurs to avoid irregular mesh 
    behavior post-division.

    :param: None
    :returns: None
    """
    # Set all parameter values
    bpy.ops.object.modifier_add(type='CLOTH')
    bpy.context.object.modifiers["Cloth"].settings.bending_model = 'LINEAR'
    bpy.context.object.modifiers["Cloth"].settings.quality = 3
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1
    bpy.context.object.modifiers["Cloth"].settings.mass = 0.3
    bpy.context.object.modifiers["Cloth"].settings.air_damping = 10
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


def make_collection(name, type): 
    # allowed_type = ['cell', 'force', 'motion', 'global']
    # if type in allowed_type: 
    collection = bpy.data.collections.new(name)
    collection['type'] = type
    collection['object'] = 'collection'
    # link the collection to the scene for visualization 
    bpy.context.scene.collection.children.link(collection)

    '''if type == 'force': 
        # check if global collection is already created 
        if any(c.name == "Global force collection" for c in bpy.data.collections):
            print(f'Global force collection already created')
            global_force_collection = bpy.data.collections["Global force collection"]
            global_force_collection['type'] = 'global'
        else: 
            print(f'Global force collection created')
            global_force_collection = bpy.data.collections.new("Global force \
                collection")
            global_force_collection['type'] = 'global'
            bpy.context.scene.collection.children.link(global_force_collection)
            bpy.context.scene.collection.children.unlink(collection)
            global_force_collection.children.link(collection)
            return collection, global_force_collection
    else: '''

    return collection
    # else:
    #    print(f"Only the following types are allowed: {allowed_type}")

    '''if type == 'force' or type == 'motion':
        # check if global collection is already created
        if any(c.name == "Global force collection" for c in bpy.data.collections):
            print(f'Global force collection already created')
            global_force_collection = bpy.data.collections["Global force collection"]
            global_force_collection['type'] = 'global'
        else: 
            print(f'Global force collection created')
            global_force_collection = bpy.data.collections.new("Global force \
                collection")
            global_force_collection['type'] = 'global'
            bpy.context.scene.collection.children.link(global_force_collection)
            bpy.context.scene.collection.children.unlink(collection)
            global_force_collection.children.link(collection)
            return collection, global_force_collection
    else:
        return collection

else:
    print(f"{type} is not a valid collection type, should be in {allowed_type}")
    return'''


def make_cell(
    name,
    loc,
    type,
    remeshing=True,
    scale=(1, 1, 1),
    stiffness=1,
    material=("bubble", 0.007, 0.021, 0.3)
):
    # Function implementation
    # add mesh type as an cell argument? 
    """Core function: creates a Blender mesh corresponding to a Goo :class:`Cell`
    object.

    :param :class:`Cell` cell: The Goo :class:`Cell` object.
    :returns: None
    """

    collection = make_collection(f'{name}_collection', type=type)

    # Initialize Cell python object
    cell = Cell(name, loc, material, scale)

    '''bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=2,
                                          radius=cell.data['radius'],
                                          calc_uvs=True,
                                          enter_editmode=False,
                                          align='WORLD',
                                          location=loc,
                                          scale=(1, 1, 1))'''
    
    # Create mesh
    bpy.ops.mesh.primitive_round_cube_add(change=False,
                                          radius=cell.data['radius'],
                                          size=cell.data['size'],
                                          arc_div=cell.data['arcdiv'],
                                          lin_div=0,
                                          div_type='CORNERS',
                                          odd_axis_align=False,
                                          no_limit=False,
                                          location=cell.data['location'])

    # Give the Blender object the cell's name
    obj = bpy.context.object

    bpy.context.object.scale[0] = scale[0]
    bpy.context.object.scale[1] = scale[1]
    bpy.context.object.scale[2] = scale[2]

    bpy.context.object.name = name
    bpy.context.object['object'] = 'cell'
    bpy.context.object['past position'] = obj.location
    bpy.context.object['current position'] = obj.location
    bpy.context.object['displacement'] = 0 

    # bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS', center='MEDIAN')
    # bpy.context.view_layers.active = bpy.data.objects[cell.data['name']]

    # Smooth the mesh
    bpy.ops.object.select = True
    bpy.ops.object.shade_smooth()

    # Add subsurface modifier to make smoother
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = cell.data['subdiv']
    # subdiv_mod = obj.modifiers[-1]
    # bpy.ops.object.modifier_apply(modifier=subdiv_mod.name)

    # Add cloth settings 
    bpy.ops.object.modifier_add(type='CLOTH')
    bpy.context.object.modifiers['Cloth'].settings.quality = 3
    bpy.context.object.modifiers['Cloth'].settings.air_damping = 0
    bpy.context.object.modifiers['Cloth'].settings.bending_model = 'ANGULAR'
    bpy.context.object.modifiers["Cloth"].settings.mass = cell.data['vertex_mass']
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1

    # Cloth > Stiffness 
    bpy.context.object.modifiers['Cloth'].settings.tension_stiffness = stiffness
    bpy.context.object.modifiers['Cloth'].settings.compression_stiffness = stiffness
    bpy.context.object.modifiers['Cloth'].settings.shear_stiffness = stiffness
    bpy.context.object.modifiers['Cloth'].settings.bending_stiffness = stiffness
    # Cloth > Damping
    bpy.context.object.modifiers['Cloth'].settings.tension_damping = 50
    bpy.context.object.modifiers['Cloth'].settings.compression_damping = 50
    bpy.context.object.modifiers['Cloth'].settings.shear_damping = 50
    bpy.context.object.modifiers['Cloth'].settings.bending_damping = 0.5
    # Cloth > Internal Springs
    bpy.context.object.modifiers['Cloth'].settings.use_internal_springs = True
    bpy.context.object.modifiers['Cloth'].settings.internal_spring_max_length = 1
    bpy.context.object.modifiers['Cloth'].settings.internal_spring_max_diversion \
        = 0.785398
    bpy.context.object.modifiers['Cloth'].settings.internal_spring_normal_check = True
    bpy.context.object.modifiers['Cloth'].settings.internal_tension_stiffness = 1
    bpy.context.object.modifiers['Cloth'].settings.internal_compression_stiffness \
        = 10
    bpy.context.object.modifiers['Cloth'].settings.internal_tension_stiffness_max \
        = 10000
    bpy.context.object.modifiers['Cloth'].settings.internal_compression_stiffness_max \
        = 10000
    # Cloth > Pressure
    bpy.context.object.modifiers["Cloth"].settings.use_pressure = True
    bpy.context.object.modifiers['Cloth'].settings.uniform_pressure_force = 5
    bpy.context.object.modifiers['Cloth'].settings.use_pressure_volume = True
    bpy.context.object.modifiers['Cloth'].settings.target_volume = 1
    bpy.context.object.modifiers['Cloth'].settings.pressure_factor = 1
    bpy.context.object.modifiers['Cloth'].settings.fluid_density = 1.05
    # Cloth > Collisions
    bpy.context.object.modifiers['Cloth'].collision_settings.collision_quality = 3
    bpy.context.object.modifiers['Cloth'].collision_settings.use_collision = True
    bpy.context.object.modifiers['Cloth'].collision_settings.use_self_collision = True
    bpy.context.object.modifiers['Cloth'].collision_settings.self_friction = 50
    bpy.context.object.modifiers['Cloth'].collision_settings.friction = 80
    bpy.context.object.modifiers['Cloth'].collision_settings.self_distance_min = 0.001
    bpy.context.object.modifiers['Cloth'].collision_settings.distance_min = 0.001
    bpy.context.object.modifiers['Cloth'].collision_settings.self_impulse_clamp = 0

    # Collision
    bpy.ops.object.modifier_add(type='COLLISION')
    bpy.context.object.modifiers['Collision'].settings.use_culling = True
    bpy.context.object.modifiers['Collision'].settings.damping = 1
    bpy.context.object.modifiers['Collision'].settings.thickness_outer = 0.02
    bpy.context.object.modifiers['Collision'].settings.thickness_inner = 0.2
    bpy.context.object.modifiers['Collision'].settings.cloth_friction = 80

    '''bpy.ops.object.constraint_add(type='MAINTAIN_VOLUME')
    bpy.context.object.constraints["Maintain Volume"].volume = 1
    bpy.context.object.constraints["Maintain Volume"].owner_space = 'WORLD'
    bpy.context.object.constraints["Maintain Volume"].free_axis = 'SAMEVOL_Z'''

    if remeshing: 
        '''bpy.ops.object.modifier_add(type='REMESH')
        bpy.context.object.modifiers["Remesh"].mode = 'SMOOTH'
        bpy.context.object.modifiers["Remesh"].octree_depth = 3
        bpy.context.object.modifiers["Remesh"].scale = 0.99
        bpy.context.object.modifiers["Remesh"].use_remove_disconnected = True
        bpy.context.object.modifiers["Remesh"].use_smooth_shade = False'''
        # bpy.ops.object.modifier_move_to_index(modifier="Remesh", index=1)

        remesh_mod = obj.modifiers.new(name=f"Remesh_{obj.name}", type='REMESH')
        obj.modifiers[f"Remesh_{obj.name}"].name = f"Remesh_{obj.name}"
        remesh_mod.mode = 'SMOOTH'
        remesh_mod.octree_depth = 4
        remesh_mod.scale = 0.85
        remesh_mod.use_remove_disconnected = True
        remesh_mod.use_smooth_shade = True
        remesh_mod.show_in_editmode = True
        bpy.ops.object.modifier_move_to_index(modifier=f"Remesh1_{obj.name}", index=3)

        '''bpy.ops.object.modifier_add(type='REMESH')
        bpy.context.object.modifiers["Remesh"].mode = 'VOXEL'
        bpy.context.object.modifiers["Remesh"].voxel_size = 0.3
        bpy.context.object.modifiers["Remesh"].use_smooth_shade = True'''

    # add material, default is purple bubble
    bpy.context.view_layer.objects.active = bpy.data.objects[name]
    mat = add_material(material[0], 
                       float(material[1]), 
                       float(material[2]), 
                       float(material[3])
                       )
    bpy.context.active_object.data.materials.append(mat)

    # remove duplicate objects outside of the collection
    bpy.ops.collection.objects_remove_all()
    # Add the active cell to our specific collection 
    bpy.data.collections[collection.name].objects.link(bpy.data.objects[name])

    # handler = handler_class()


# Defines the Cell class
class Cell():
    '''Core class: creates Goo cell object. 

    :param str name: The name of the cell.
    :param tuple loc: The coordinates of the cell.   
    :returns None:
    '''
    
    def __init__(self, name, loc, material, scale=(1, 1, 1), mesh_type=""):
        '''
        Constructor method
        '''

        # The initialization function sets a cell data dictionary
        # for geometric parameters, physics parameters,
        # division information, and cell lineage
        self.data = {
            'ID': 0,
            'name': name,
            'radius': 1,
            'enter_editmode': False,
            'align': 'WORLD',
            'location': loc,
            'material': material,
            'size': (1, 1, 1),
            'scale': scale,
            'arcdiv': 8,  # changes vertices in object mode
            'subdiv': 1,  # changes vertices in edit mode (of the cloth)
            'vertex_mass': 1,
            'density': 1.05,
            'phi': 0,
            'theta': 0,
            'mother': 'none',
            'daughters': ['none', 'none'],
            'mesh_type': mesh_type
            }
        # The volume and mass are calculated from values in the data dictionary
        self.data['init_volume'] = ((4/3)*np.pi*(self.data['radius'])**3)
        self.data['mass'] = self.data['density']*self.data['init_volume']
        print(self.data['mass'])


def add_material(mat_name, r, g, b):
    """Core function: creates a soap bubble-like Blender material for use in rendering 
    cells.

    The material has a name that allows it to be shared across multiple cells. 

    :param str mat_name: The name of the material. 
    :param float r: The value of the red in RGB [0 to 1]. 
    :param float g: The value of the green value in RGB [0 to 1]. 
    :param float b: The value of the blue value in RGB [0 to 1]. 
    :returns: None
    """

    if bpy.data.materials.get(mat_name): 
        mat = bpy.data.materials.get(mat_name)
        mat.diffuse_color = (0.1, 0, 0, 0.8)  # viewport color
        mat.use_nodes = True
        mat.blend_method = 'BLEND'

        # get the material nodes
        nodes = mat.node_tree.nodes

        # clear all nodes to start clean
        for node in nodes:
            nodes.remove(node)

        # create principled node for main color
        node_main = nodes.new(type='ShaderNodeBsdfPrincipled')
        node_main.location = -200, 100
        node_main.inputs['Base Color'].default_value = (r, g, b, 0.8)
        node_main.inputs['Metallic'].default_value = 0.136
        node_main.inputs['Specular'].default_value = 0.500
        node_main.inputs['Specular Tint'].default_value = 0.555
        node_main.inputs['Roughness'].default_value = 0.318
        node_main.inputs['Anisotropic'].default_value = 0.041
        node_main.inputs['Anisotropic Rotation'].default_value = 0.048
        node_main.inputs['Sheen'].default_value = 0.052
        node_main.inputs['Sheen Tint'].default_value = 0.030
        node_main.inputs['Clearcoat'].default_value = 0.114
        node_main.inputs['Clearcoat Roughness'].default_value = 0.123
        node_main.inputs['IOR'].default_value = 1.450
        node_main.inputs['Transmission'].default_value = 0.882
        node_main.inputs['Transmission Roughness'].default_value = 0.0
        node_main.inputs['Alpha'].default_value = 0.414

        # create noise texture source
        node_noise = nodes.new(type="ShaderNodeTexNoise")
        node_noise.inputs['Scale'].default_value = 0.600
        node_noise.inputs['Detail'].default_value = 15.0
        node_noise.inputs['Roughness'].default_value = 0.500
        node_noise.inputs['Distortion'].default_value = 3.0

        # create HSV
        node_HSV = nodes.new(type="ShaderNodeHueSaturation")
        node_HSV.inputs['Hue'].default_value = 0.800
        node_HSV.inputs['Saturation'].default_value = 2.00
        node_HSV.inputs['Value'].default_value = 2.00
        node_HSV.inputs['Fac'].default_value = 1.00

        # create second principled node for random color variation
        node_random = nodes.new(type='ShaderNodeBsdfPrincipled')
        node_random.location = -200, -100
        node_random.inputs['Base Color'].default_value = (r, g, b, 1)
        node_random.inputs['Metallic'].default_value = 0.0
        node_random.inputs['Specular'].default_value = 0.500
        node_random.inputs['Specular Tint'].default_value = 0.0
        node_random.inputs['Roughness'].default_value = 0.482
        node_random.inputs['Anisotropic'].default_value = 0.0
        node_random.inputs['Anisotropic Rotation'].default_value = 0.0
        node_random.inputs['Sheen'].default_value = 0.0
        node_random.inputs['Sheen Tint'].default_value = 0.0
        node_random.inputs['Clearcoat'].default_value = 0.0
        node_random.inputs['Clearcoat Roughness'].default_value = 0.0
        node_random.inputs['IOR'].default_value = 1.450
        node_random.inputs['Transmission'].default_value = 1.0
        node_random.inputs['Transmission Roughness'].default_value = 0.0
        node_random.inputs['Alpha'].default_value = 0.555

        # create mix shader node
        node_mix = nodes.new(type='ShaderNodeMixShader')
        node_mix.location = 0, 0
        node_mix.inputs['Fac'].default_value = 0.079

        # create output node
        node_output = nodes.new(type='ShaderNodeOutputMaterial')
        node_output.location = 200, 0

        # link nodes
        links = mat.node_tree.links
        links.new(node_noise.outputs[1], node_HSV.inputs[4])  # link_noise_HSV
        links.new(node_HSV.outputs[0], node_random.inputs[0])  # link_HSV_random
        links.new(node_main.outputs[0], node_mix.inputs[1])  # link_main_mix
        links.new(node_random.outputs[0], node_mix.inputs[2])  # link_random_mix
        links.new(node_mix.outputs[0], node_output.inputs[0])  # link_mix_out

    else: 
        mat = bpy.data.materials.new(name=mat_name)
        mat.diffuse_color = (0.1, 0, 0, 0.8)  # viewport color
        mat.use_nodes = True
        mat.blend_method = 'BLEND'

        # get the material nodes
        nodes = mat.node_tree.nodes

        # clear all nodes to start clean
        for node in nodes:
            nodes.remove(node)

        # create principled node for main color
        node_main = nodes.new(type='ShaderNodeBsdfPrincipled')
        node_main.location = -200, 100
        node_main.inputs['Base Color'].default_value = (r, g, b, 0.8)
        node_main.inputs['Metallic'].default_value = 0.136
        node_main.inputs['Specular'].default_value = 0.500
        node_main.inputs['Specular Tint'].default_value = 0.555
        node_main.inputs['Roughness'].default_value = 0.318
        node_main.inputs['Anisotropic'].default_value = 0.041
        node_main.inputs['Anisotropic Rotation'].default_value = 0.048
        node_main.inputs['Sheen'].default_value = 0.052
        node_main.inputs['Sheen Tint'].default_value = 0.030
        node_main.inputs['Clearcoat'].default_value = 0.114
        node_main.inputs['Clearcoat Roughness'].default_value = 0.123
        node_main.inputs['IOR'].default_value = 1.450
        node_main.inputs['Transmission'].default_value = 0.882
        node_main.inputs['Transmission Roughness'].default_value = 0.0
        node_main.inputs['Alpha'].default_value = 0.414

        # create noise texture source
        node_noise = nodes.new(type="ShaderNodeTexNoise")
        node_noise.inputs['Scale'].default_value = 0.600
        node_noise.inputs['Detail'].default_value = 15.0
        node_noise.inputs['Roughness'].default_value = 0.500
        node_noise.inputs['Distortion'].default_value = 3.0

        # create HSV
        node_HSV = nodes.new(type="ShaderNodeHueSaturation")
        node_HSV.inputs['Hue'].default_value = 0.800
        node_HSV.inputs['Saturation'].default_value = 2.00
        node_HSV.inputs['Value'].default_value = 2.00
        node_HSV.inputs['Fac'].default_value = 1.00

        # create second principled node for random color variation
        node_random = nodes.new(type='ShaderNodeBsdfPrincipled')
        node_random.location = -200, -100
        node_random.inputs['Base Color'].default_value = (r, g, b, 1)
        node_random.inputs['Metallic'].default_value = 0.0
        node_random.inputs['Specular'].default_value = 0.500
        node_random.inputs['Specular Tint'].default_value = 0.0
        node_random.inputs['Roughness'].default_value = 0.482
        node_random.inputs['Anisotropic'].default_value = 0.0
        node_random.inputs['Anisotropic Rotation'].default_value = 0.0
        node_random.inputs['Sheen'].default_value = 0.0
        node_random.inputs['Sheen Tint'].default_value = 0.0
        node_random.inputs['Clearcoat'].default_value = 0.0
        node_random.inputs['Clearcoat Roughness'].default_value = 0.0
        node_random.inputs['IOR'].default_value = 1.450
        node_random.inputs['Transmission'].default_value = 1.0
        node_random.inputs['Transmission Roughness'].default_value = 0.0
        node_random.inputs['Alpha'].default_value = 0.555

        # create mix shader node
        node_mix = nodes.new(type='ShaderNodeMixShader')
        node_mix.location = 0, 0
        node_mix.inputs['Fac'].default_value = 0.079

        # create output node
        node_output = nodes.new(type='ShaderNodeOutputMaterial')
        node_output.location = 200, 0

        # link nodes
        links = mat.node_tree.links
        links.new(node_noise.outputs[1], node_HSV.inputs[4])  # link_noise_HSV
        links.new(node_HSV.outputs[0], node_random.inputs[0])  # link_HSV_random
        links.new(node_main.outputs[0], node_mix.inputs[1])  # link_main_mix
        links.new(node_random.outputs[0], node_mix.inputs[2])  # link_random_mix
        links.new(node_mix.outputs[0], node_output.inputs[0])  # link_mix_out

    return mat


class Force():
    """Core class: creates Goo force objects. 

    The class instantiates :class:`Force` objects that represent 
    adhesion forces between cells in Blender. 
    
    :param str force_name: The name of the force.
    :param str cell_name: The name of the cell. 
    :param float strength: The strength of the force.
    :param float falloff_power: The power of the falloff of the force. 
    :returns: None

    .. note:: ``falloff_power`` is a positive (:py:class:`float`). 
        By default, the type of the falloff is set to `SPHERE` 
        and its shape is set to `SURFACE`. 
    """

    force_list = []

    def __init__(self, 
                 force_name, 
                 cell_name, 
                 strength, 
                 falloff_power, 
                 collection_name, 
                 type
                 ):
        self.name = force_name
        self.strength = strength
        self.associated_cell = cell_name
        self.falloff_power = falloff_power
        self.collection = collection_name
        self.falloff_type = 'SPHERE'
        self.shape = 'SURFACE'
        self.type = type  # boolean
        if not self.type:
            Force.force_list.append(self)  # list of all the created forces

    def get_strength(self): 
        return self.strength
    
    def get_falloff(self): 
        return self.falloff_power

    def get_blender_force(self):
        obj = bpy.data.objects[self.name]
        return obj


def add_boundaries(loc, size, shape='BOX', type='REFLECTIVE', name='box'):

    if not isinstance(loc, tuple) or len(loc) != 3:
        raise ValueError(
            "Invalid 'loc' argument. It should be a tuple containing X, Y, and Z \
                coordinates."
            )

    if not all(isinstance(coord, (int, float)) for coord in loc):
        raise ValueError(
            "Invalid 'loc' coordinates. Coordinates should be integers or floats."
            )

    if not isinstance(size, tuple) or len(size) != 3:
        raise ValueError(
            "Invalid 'size' argument. It should be a tuple containing dimensions in X, \
                Y, and Z."
            )

    if not all(isinstance(dim, (int, float)) for dim in size):
        raise ValueError(
            "Invalid 'size' dimensions. Dimensions should be integers or floats."
            )

    if shape not in ['BOX', 'SPHERE']:
        raise ValueError(
            "Invalid 'shape' argument. Supported shapes are 'BOX' and 'SPHERE'."
            )

    if shape == 'BOX':
        bpy.ops.mesh.primitive_cube_add(enter_editmode=False,
                                        align='WORLD',
                                        location=loc
                                        )
        bpy.context.object.name = name
        bpy.context.object.scale[0] = -size[0]
        bpy.context.object.scale[1] = -size[1]
        bpy.context.object.scale[2] = -size[2]
        bpy.ops.object.modifier_add(type='COLLISION')
        bpy.ops.object.modifier_add(type='WIREFRAME')

    elif shape == 'SPHERE':
        bpy.ops.mesh.primitive_uv_sphere_add(radius=size[0],
                                             enter_editmode=False,
                                             align='WORLD',
                                             location=loc)
        bpy.context.object.name = name
        bpy.context.object.scale[0] = -size[0]
        bpy.context.object.scale[1] = -size[0]
        bpy.context.object.scale[2] = -size[0]
        bpy.ops.object.modifier_add(type='COLLISION')
        bpy.ops.object.modifier_add(type='WIREFRAME')


def add_motion(effector_name,
               strength,
               persistence=0,
               randomness=1,
               distribution='uniform',
               size=0.1):
    if (isinstance(effector_name, str)
            and isinstance(strength, (int, float))
            and isinstance(persistence, (int, float))
            and isinstance(randomness, (int, float))):

        cell_type = f"{bpy.data.objects[effector_name].users_collection[0]['type']}"
        force = make_force(force_name=f'motion_{effector_name}',
                           cell_name=f'{effector_name}',
                           type=cell_type,
                           strength=strength,
                           falloff=0,
                           motion=True,
                           min_dist=0,
                           max_dist=4)
        force = bpy.data.objects[force.name]
        print(force)

        force['persistence'] = persistence
        force['randomness'] = randomness
        force['distribution'] = distribution
        force['distribution size'] = size

        '''cell = bpy.data.objects[force.get('cell')]
        com = get_centerofmass(cell)
        dir_vec = (Vector(com) - bpy.data.objects[force.name].location)
        dir_vec.normalize()
        force['persistent direction'] = dir_vec

        # Create curve and delete one of the initial points
                
        # Define the list of coordinates
        coords = [(com[0], com[1], com[2])]

        # Create a new curve data block
        curve_cell_data = bpy.data.curves.new(name='cell motion', type='CURVE')
        curve_cell_data.dimensions = '3D'
        # Create a new spline and add it to the curve data block
        spline = curve_cell_data.splines.new(type='POLY')
        spline.points.add(len(coords)-1)
        # Set the coordinates for each control point in the spline
        for i, coord in enumerate(coords):
            x, y, z = coord
            spline.points[i].co = (x, y, z, 1)


        curve_force_data = bpy.data.curves.new(name='force motion', type='CURVE')
        curve_force_data.dimensions = '3D'
        # Create a new spline and add it to the curve data block
        spline_force = curve_force_data.splines.new(type='POLY')
        spline_force.points.add(len(coords)-1)
        # Set the coordinates for each control point in the spline
        for i, coord in enumerate(coords):
            x, y, z = coord
            spline_force.points[i].co = (x, y, z, 1)

        # Create a new object and link it to the scene
        curve_cell_obj = bpy.data.objects.new(f'{cell.name}_tracks', curve_cell_data)
        curve_force_obj = bpy.data.objects.new(f'{force.name}_force_tracks', 
                                                curve_force_data)

        #bpy.context.scene.collection.objects.link(curve_obj)
        cell.users_collection[0].objects.link(curve_cell_obj)
        cell.users_collection[0].objects.link(curve_force_obj)

        # Set the object to be in edit mode and select all points
        bpy.context.view_layer.objects.active = curve_cell_obj
        bpy.context.object.data.bevel_depth = 0.008
        # Assign the material to the curve object

        cell_curve_mat = add_material('cell_curve', 0, 1, 0)
        if curve_cell_obj.data.materials:
            curve_cell_obj.data.materials[0] = cell_curve_mat
        else:
            curve_cell_obj.data.materials.append(cell_curve_mat)

        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.curve.select_all(action='SELECT')
        # Update the curve display and exit edit mode
        bpy.ops.curve.reveal()
        bpy.ops.object.mode_set(mode='OBJECT')


        bpy.context.view_layer.objects.active = curve_force_obj
        bpy.context.object.data.bevel_depth = 0.008

        force_curve_mat = add_material('force_curve', 1, 0, 0)
        if curve_force_obj.data.materials:
            curve_force_obj.data.materials[0] = force_curve_mat
        else:
            curve_force_obj.data.materials.append(force_curve_mat)

        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.curve.select_all(action='SELECT')
        # Update the curve display and exit edit mode
        bpy.ops.curve.reveal()
        bpy.ops.object.mode_set(mode='OBJECT')'''

    else: 
        print('---- Types not supported')
        
        '''# collective motion on a whole cell type
        if type == 'collective':
            make_force(f'motion_{effector_name}',
                        f'{effector_name}',
                        f"{effector_name}",
                        strength, 0,
                        motion=True,
                        collective=True,
                        min_dist=0, max_dist=10)

        # individual motion for a single cell
        elif type == 'random' or type == 'persistent':
            make_force(f'motion_{effector_name}',
                        f'{effector_name}',
                        f"{bpy.data.objects[effector_name].users_collection[0]['type']}",
                        strength, 0,
                        motion=True,
                        collective =False,
                        min_dist=0, max_dist=10)'''
        
    '''if type == 'collection':
        if (bpy.data.collections[collection_name]
            and bpy.data.collections[collection_name]['type'] == 'cell'):
            for cell in bpy.data.collections.get(collection_name).all_objects:
                make_force(f'motion_{cell.name}',
                           f'{cell.name}',
                           strength=strength,
                           falloff=0,
                           motion=True,
                           min_dist=0,
                           max_dist=10)
    elif type == 'individual':
        for cell in bpy.data.collections.get(collection_name).all_objects:
            collection = make_collection(f'forces_{cell.name}', type = 'force')
            collection.children.link(cell['adhesion force'])
            collection.children.link(cell['motion force'])'''


def add_hetero_adhesion(cell_name, other_type_name, strength):

    hetero_collections = [
        coll for coll in bpy.data.collections if coll.get('type') in other_type_name
    ]

    force = make_force(force_name=f'heterotypic_{cell_name}_{other_type_name}',
                       cell_name=cell_name,
                       type=type,
                       strength=strength,
                       falloff=0,
                       motion=False)
    
    for coll in hetero_collections:
        coll.objects.link(bpy.data.objects[force.name])


def add_homo_adhesion(cell_name, strength):

    homo_collections = [
        coll for coll in bpy.data.collections
        if coll.get('type') 
        in bpy.data.objects[cell_name].users_collection[0].get("type")
    ]

    cell_type = f"{bpy.data.objects[cell_name].users_collection[0]['type']}"
    force = make_force(force_name=f'homotypic_{cell_name}',
                       cell_name=cell_name,
                       type=cell_type,
                       strength=strength,
                       falloff=0,
                       motion=False)
    
    for coll in homo_collections:
        coll.objects.link(bpy.data.objects[force.name])


def make_force(force_name,
               cell_name,
               type,
               strength,
               falloff,
               motion=False,
               min_dist=0.5,
               max_dist=1.5):
    """Core function: creates a Blender force from a Goo :class:`Force` object.

    :param :class:`Force` force: The Goo force object.
    :returns: None
    """
    collection = bpy.data.objects[cell_name].users_collection[0]
    force = Force(force_name, cell_name, strength, falloff, collection.name, motion)

    # Add a force object
    cell = force.associated_cell

    if not motion:
        bpy.ops.object.effector_add(
            type='FORCE',
            enter_editmode=False,
            align='WORLD',
            location=get_centerofmass(bpy.data.objects[cell]),
            scale=(1, 1, 1)
        )
        bpy.context.object['motion'] = False
        bpy.context.object.name = force_name
        bpy.data.objects.get(cell_name)["adhesion force"] = force_name

    elif motion:

        rand_coord = tuple(np.random.uniform(low=-0.05, high=0.05, size=(3,)))
        force_rand_coord = zip(get_centerofmass(bpy.data.objects[cell]), rand_coord)
        print(f"Random displacement new coord: "
              f"{tuple(map(sum, force_rand_coord))}")
        cell_com = get_centerofmass(bpy.data.objects[cell])
        bpy.ops.object.effector_add(type='FORCE',
                                    enter_editmode=False,
                                    align='WORLD',
                                    location=tuple(map(sum, zip(cell_com, rand_coord))),
                                    scale=(1, 1, 1))
        bpy.context.object['motion'] = True
        bpy.context.object.name = force_name
        bpy.data.objects.get(cell_name)["motion force"] = force_name
        collection.objects.link(bpy.data.objects[force_name])

    else:
        print('---- Unsupported forces ----')

    '''elif motion == True and collective == True: 
        bpy.ops.object.effector_add(type='FORCE',
                                enter_editmode=False,
                                align='WORLD',
                                location= tuple(np.random.uniform(low=-5, 
                                                                  high=5, 
                                                                  size=(3,))),
                                scale=(1, 1, 1))
        bpy.context.object['motion'] = True
        bpy.context.object.name = force_name
        for coll in same_collections: 
            for cell in coll.all_objects: 
                if cell.get('object') == 'cell': 
                    bpy.data.objects.get(cell.name)["motion force"] = force_name
            coll.objects.link(bpy.data.objects[force_name])'''

    # Add force parameters
    bpy.context.object.field.strength = force.strength
    bpy.context.object.field.use_max_distance = True
    bpy.context.object.field.use_min_distance = True
    bpy.context.object.field.distance_max = max_dist
    bpy.context.object.field.distance_min = min_dist
    bpy.context.object.name = force.name
    bpy.context.object.field.falloff_power = force.falloff_power

    bpy.context.object['cell'] = cell_name
 
    scene_collection = bpy.context.scene.collection
    scene_collection.objects.unlink(bpy.data.objects[force.name])

    return force


def setup_world(seed=None):
    """Auxiliary function: sets up the default values used for simulations in Goo 
    including units and rendering background. 

    :param seed: Seed value for random number generation (optional)
    :type seed: int or None
    :returns: None
    """
    # Set the random seed for reproducibility
    if seed is None:
        seed = int(datetime.now().timestamp())

    np.random.seed(seed)
    bpy.context.scene['seed'] = seed

    # Turn off gravity so cells don't fall in the simulation
    bpy.context.scene.use_gravity = False
    # Set units to the metric system
    bpy.context.scene.unit_settings.system = 'METRIC'
    bpy.context.scene.unit_settings.scale_length = 1
    bpy.context.scene.unit_settings.system_rotation = 'DEGREES'
    bpy.context.scene.unit_settings.length_unit = 'MICROMETERS'
    bpy.context.scene.unit_settings.mass_unit = 'MILLIGRAMS'
    bpy.context.scene.unit_settings.time_unit = 'SECONDS'
    bpy.context.scene.unit_settings.temperature_unit = 'CELSIUS'

    '''# Delete all existing collections 
    for collection in bpy.data.collections:  # loop through the existing collection
        # loop through objects in collection
        for objs in collection.objects:
            # delete existing objects in collection 
            bpy.data.objects.remove(objs)
        # Delete collection
        bpy.data.collections.remove(collection)'''
    
    # Delete all existing objects in the scene except cameras and lights
    for obj in bpy.context.scene.objects:
        if obj.type not in ['CAMERA', 'LIGHT']:
            bpy.data.objects.remove(obj)

    # Delete all existing collections 
    for collection in bpy.data.collections:
        # Delete collection
        bpy.data.collections.remove(collection)

    # Add an HDRI image for illumination
    add_world_HDRI()

    # Change the Viewport Shading to Rendered
    for area in bpy.data.screens[3].areas: 
        if area.type == 'VIEW_3D':
            for space in area.spaces: 
                if space.type == 'VIEW_3D':
                    space.shading.type = 'WIREFRAME'
                    space.overlay.show_overlays = False


def add_world_HDRI():
    """Auxilliary function: sets up Blender World properties for use in rendering.

    It adds an HDRI image for illumination. 

    :returns: None
    """
    C = bpy.context
    scn = C.scene

    # Get the environment node tree of the current scene
    node_tree = scn.world.node_tree
    tree_nodes = node_tree.nodes
    # Clear all nodes
    tree_nodes.clear()
    # Add Background node
    node_background = tree_nodes.new(type='ShaderNodeBackground')
    # Add Environment Texture node
    node_environment = tree_nodes.new('ShaderNodeTexEnvironment')
    # Load and assign the image to the node property
    scripts_paths = bpy.utils.script_paths()

    # Relative path- this file must be in same directory as blend file
    try:
        node_environment.image = bpy.data.images.load(
            scripts_paths[-1]+"/modules/goo/missile_launch_facility_01_4k.hdr")
    # If the user does not have this file in the right place, throw exception
    except Exception:
        print(sys.exc_info())
        print("WARNING FROM GOO: To enable proper rendering you must have")
        print("/modules/goo/missile_launch_facility_01_4k.hdr")
        print("in the right location")

    node_environment.location = -300, 0

    # Add Output node
    node_output = tree_nodes.new(type='ShaderNodeOutputWorld')
    node_output.location = 200, 0
    # Link all nodes
    links = node_tree.links
    links.new(node_environment.outputs["Color"], node_background.inputs["Color"])
    links.new(node_background.outputs["Background"], node_output.inputs["Surface"])
    # set film to transparent to hide background
    bpy.context.scene.render.film_transparent = True
    # change render preview mode
    # only updates windows in current tab, e.g. Sxripting but not Layout
    for area in bpy.context.screen.areas:
        if area.type == 'VIEW_3D':
            print("Update #D Viewport to 'RENDERED'")
            space = area.spaces.active
            if space.type == 'VIEW_3D':
                space.shading.type = 'WIREFRAME'


def render(file_path, start, end):
    """Auxilliary function: renders a simulation to create a set of still images that 
    can be made into a movie

    :param str file_path: The path of the folder used to store output images. 
    :param bpy.context.scene scene: The Blender current scene. 
    :param int start: The Blender starting frame. 
    :param int end: The Blender ending frame. 
    :returns: None
    """
    scene = bpy.context.scene
    # Set the image file format as PNG
    scene.render.image_settings.file_format = 'PNG'
    # Set the file path where the images will be saved
    old_fp = scene.render.filepath
    scene.render.filepath = file_path
    # Set the starting and ending frames for the simulation
    scene = bpy.context.scene
    scene.frame_start = start
    scene.frame_end = end
    # Set the handlers for the simulation
    handlers = bpy.app.handlers.frame_change_post.copy()
    bpy.app.handlers.frame_change_post.clear()
    # Loop through each frame
    for frame in range(scene.frame_start, scene.frame_end):
        # Set the frame
        bpy.context.scene.frame_set(frame)
        # Run each handler
        for func in handlers:
            func(scene)
        # Save the image
        file_name = "frame" + str(scene.frame_current)
        scene.render.filepath += file_name
        bpy.ops.render.render(write_still=True)
        scene.render.filepath = scene.render.filepath.removesuffix(file_name)
    scene.render.filepath = old_fp
    # Add each handler to the scene
    for func in handlers:
        bpy.app.handlers.frame_change_post.append(func)
    return


# collections are currently manually created in Blender - not used
def make_force_collections(master_collection, cell_types):
    """Auxilliary function: makes collections for forces to be stored in.

    :param bpy.context.view_layer.active_layer_collection master_collection: The 
    collection in which the force collections will be contained. 
    :param cell_types: list of active cell types
    :returns: None
    """
    bpy.context.view_layer.active_layer_collection = master_collection
    for type in cell_types:
        collection = bpy.context.blend_data.collections.new(name=type+"_forces")
        bpy.context.collection.children.link(collection)


def calculate_com(cell): 
    """Function to calculate the center of mass of a cell Blender object. 

    :param bpy.data.objects[cell.name] cell: cell
    :returns: The coordinate of the center of mass of the cell. 
    :rtype: Tuple(x,y,z)
    """
    bpy.context.view_layer.objects.active = cell
    dg = bpy.context.evaluated_depsgraph_get()
    cell_eval = cell.evaluated_get(dg)
    vertices = cell_eval.data.vertices
    vert_coords = np.asarray([(cell_eval.matrix_world @ v.co) for v in vertices])

    x = vert_coords[:, 0]
    y = vert_coords[:, 1]
    z = vert_coords[:, 2]
    COM = (np.mean(x), np.mean(y), np.mean(z))
    return COM


# not used
def get_contact_area_raycast():

    # Get the mesh objects
    obj1 = bpy.context.scene.objects["cell_A1"]
    obj2 = bpy.context.scene.objects["cell_A2"]

    # Get the mesh data for each object
    mesh_data1 = obj1.to_mesh()
    mesh_data2 = obj2.to_mesh()

    # Create BMesh objects for each mesh
    bm1 = bmesh.new()
    bm1.from_mesh(mesh_data1)

    bm2 = bmesh.new()
    bm2.from_mesh(mesh_data2)

    # Initialize lists to hold the contact vertices for each mesh
    contact_vertices1 = []
    contact_vertices2 = []

    # Loop over the vertices of the first mesh
    for v1 in bm1.verts:
        # Cast a ray from the vertex and check for intersection with the second mesh
        hit, loc, norm, face_index = obj2.ray_cast(
            obj1.matrix_world @ v1.co, 
            -obj1.matrix_world @ v1.normal
        )

        # If there is a hit and the distance is within a certain threshold, add the 
        # vertex to the contact list
        if hit and (loc - obj1.matrix_world @ v1.co).length < 0.1:
            contact_vertices1.append(v1.co)

    # Loop over the vertices of the second mesh
    for v2 in bm2.verts:
        # Cast a ray from the vertex and check for intersection with the first mesh
        hit, loc, norm, face_index = obj1.ray_cast(
            obj2.matrix_world @ v2.co, 
            -obj2.matrix_world @ v2.normal
        )

        # If there is a hit and the distance is within a certain threshold, add the 
        # vertex to the contact list
        if hit and (loc - obj2.matrix_world @ v2.co).length < 0.1:
            contact_vertices2.append(v2.co)

    # Compute the contact area and ratio
    contact_area1 = len(set(contact_vertices1))
    contact_area2 = len(set(contact_vertices1))
    total_area1 = len(bm1.faces)
    total_area2 = len(bm2.faces)
    contact_ratio1 = contact_area1 / total_area1
    contact_ratio2 = contact_area2 / total_area2

    # Print the results
    print("Contact area 1: ", contact_area1)
    print("Contact area 2: ", contact_area2)
    print("Object 1 total area: ", total_area1)
    print("Object 2 total area: ", total_area2)
    print("Object 1 contact ratio: ", contact_ratio1)
    print("Object 2 contact ratio: ", contact_ratio2)

    # Free the BMesh objects and mesh data
    bm1.free()
    bm2.free()
    bpy.data.meshes.remove(mesh_data1)
    bpy.data.meshes.remove(mesh_data2) 
    
    return 


'''def get_contact_area():
    """
    Calculate the contact area between two meshes using ray casting.

    Parameters:
    - mesh1: The first mesh as a Blender object.
    - mesh2: The second mesh as a Blender object.
    - threshold: The maximum distance between two vertices to consider them in contact.

    Returns:
    - The total contact area between the two meshes.
    """

    threshold = 0.03

    # Get cell objects
    mesh1 = bpy.data.objects['cell_A1']
    mesh2 = bpy.data.objects['cell_A2']

    # Evaluate the meshes to account for deformations.
    mesh1_eval = mesh1.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh2_eval = mesh2.evaluated_get(bpy.context.evaluated_depsgraph_get())

    # Get the vertices of the meshes as global coordinates.
    verts1 = [mesh1_eval.matrix_world @ v.co for v in mesh1_eval.data.vertices]
    verts2 = [mesh2_eval.matrix_world @ v.co for v in mesh2_eval.data.vertices]

    # Calculate the total surface area of each mesh.
    area1 = sum(p.area for p in mesh1_eval.data.polygons)
    area2 = sum(p.area for p in mesh2_eval.data.polygons)
    #print(f"total area: {area1}, {area2}")

    # Get the bounding boxes of the meshes and make them slightly bigger.
    #bbox1 = mesh1.bound_box
    #bbox2 = mesh2.bound_box
    #margin = 2  # Change this value to control the size of the margin.
    #bbox1 = [Vector(v) + Vector((margin, margin, margin)) for v in bbox1]
    #bbox2 = [Vector(v) + Vector((margin, margin, margin)) for v in bbox2]

    # Find the vertices in contact for each mesh.
    contact_verts1 = []
    contact_verts2 = []

    # Compute the contact area for each mesh.
    contact_area1 = 0.0
    contact_area2 = 0.0

    for v1 in verts1:
        for v2 in verts2:
            # Check if the distance between the vertices is below the threshold.
            if (v1 - v2).length < threshold:
                #print(f"coord v1, v2; {v1}, {v2} and length: {(v1 - v2).length}")
                contact_verts1.append(v1)
                contact_verts2.append(v2)

    for p in mesh1_eval.data.polygons:
        poly_verts = [
            mesh1_eval.matrix_world @ mesh1_eval.data.vertices[i].co 
            for i in p.vertices
        ]

        # Check if any of the polygon vertices are in contact.
        if any((v - v1).length < threshold for v in poly_verts 
                                            for v1 in contact_verts2):
            contact_area1 += p.area

    for p in mesh2_eval.data.polygons:
        poly_verts = [
            mesh2_eval.matrix_world @ mesh2_eval.data.vertices[i].co for i in p.vertices
            ]

        # Check if any of the polygon vertices are in contact.
        if any((v - v2).length < threshold for v in poly_verts 
                                            for v2 in contact_verts1):
            contact_area2 += p.area

    # Compute the contact area ratio for each mesh.
    ratio1 = contact_area1 / area1 if area1 > 0 else 0.0
    ratio2 = contact_area2 / area2 if area2 > 0 else 0.0

    return ratio1, ratio2
'''


def get_contact_area():
    """
    Calculate the contact ratio between cells in the scene.

    Parameters:
    - scene: The Blender scene object.

    Returns:
    - A dictionary containing cells in contact as keys and their contact ratio as 
        values.
    """

    threshold = 0.05
    contact_ratio_dict = {}
    contact_areas_dict = {}

    # Get all cell objects in the scene
    cell_objects = [
        obj for obj in bpy.context.scene.objects 
        if obj.get('object') is not None and obj.get('object') == 'cell'
    ]

    # Loop over cell pairs to compute contact ratio
    for i in range(len(cell_objects) - 1):
        mesh1 = cell_objects[i]
        for j in range(i + 1, len(cell_objects)):
            mesh2 = cell_objects[j]

            com_mesh1 = get_centerofmass(mesh1)
            com_mesh2 = get_centerofmass(mesh2)
            if (Vector(com_mesh1) - Vector(com_mesh2)).length < 2: 

                # Evaluate the meshes to account for deformations.
                mesh1_eval = mesh1.evaluated_get(bpy.context.evaluated_depsgraph_get())
                mesh2_eval = mesh2.evaluated_get(bpy.context.evaluated_depsgraph_get())

                # Get the vertices of the meshes as global coordinates.
                verts1 = [
                    mesh1_eval.matrix_world @ v.co 
                    for v in mesh1_eval.data.vertices
                ]
                verts2 = [
                    mesh2_eval.matrix_world @ v.co 
                    for v in mesh2_eval.data.vertices
                ]

                # Calculate the total surface area of each mesh.
                area1 = sum(p.area for p in mesh1_eval.data.polygons)
                area2 = sum(p.area for p in mesh2_eval.data.polygons)

                # Find the vertices in contact for each mesh.
                contact_verts1 = []
                contact_verts2 = []

                # Compute the contact area for each mesh.
                contact_area1 = 0.00
                contact_area2 = 0.00

                for v1 in verts1:
                    for v2 in verts2:
                        # Check if the distance between the vertices is below the 
                        # threshold.
                        if (v1 - v2).length < threshold:
                            contact_verts1.append(v1)
                            contact_verts2.append(v2)

                for p in mesh1_eval.data.polygons:
                    poly_verts = [
                        mesh1_eval.matrix_world @ mesh1_eval.data.vertices[i].co 
                        for i in p.vertices
                    ]

                    # Check if any of the polygon vertices are in contact.
                    if any((v - v1).length < threshold for v in poly_verts 
                           for v1 in contact_verts2):
                        contact_area1 += p.area

                for p in mesh2_eval.data.polygons:
                    poly_verts = [
                        mesh2_eval.matrix_world @ mesh2_eval.data.vertices[i].co 
                        for i in p.vertices
                    ]

                    # Check if any of the polygon vertices are in contact.
                    if any((v - v2).length < threshold for v in poly_verts 
                           for v2 in contact_verts1):
                        contact_area2 += p.area

                # Compute the contact area ratio for each mesh.
                ratio1 = contact_area1 / area1 if area1 > 0 else 0.0
                ratio2 = contact_area2 / area2 if area2 > 0 else 0.0

                # avg ratio because bidirectional contact areas
                bidir_ratio = (ratio1 + ratio2) / 2
                bidir_areas = (contact_area1 + contact_area2) / 2

                # Add contact ratio to the dictionary
                contact_ratio_dict[f"{mesh1.name}-{mesh2.name}"] = bidir_ratio
                contact_areas_dict[f"{mesh1.name}-{mesh2.name}"] = bidir_areas

            else: 
                contact_ratio_dict[f"{mesh1.name}-{mesh2.name}"] = 0
                contact_areas_dict[f"{mesh1.name}-{mesh2.name}"] = 0

    return contact_ratio_dict, contact_areas_dict


def separate_mesh(obj): 
    
    # Define the two tuple coordinates that define the line
    p1 = get_centerofmass(obj)
    p2 = get_long_axis_global(obj)
    
    # Get the vector representing the line
    line_vector = Vector((p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])).normalized()
    
    # Set the active object and get the evaluated mesh
    bpy.context.view_layer.objects.active = obj
    obj_eval = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh_eval = obj_eval.to_mesh()
    
    # List of intersection vertices and their indices
    intersection_verts = []
    intersection_indices = []
    
    # Get the edit-mode mesh and create a BMesh
    mesh = obj.data
    bm = bmesh.new()
    bm.from_mesh(mesh)

    # Ensure that the index table of the BMesh is up to date
    bm.verts.ensure_lookup_table()
    
    # Loop through the vertices in the evaluated mesh
    for vert_eval in mesh_eval.vertices:
        if vert_eval.index < len(mesh.vertices):
            # Get the vertex in the edit-mode mesh using the index
            vert = bm.verts[vert_eval.index]
            
            # Transform the vertex coordinate to global coordinates
            vert_coord = obj.matrix_world @ vert.co
            
            # Calculate the distance of the vertex from the plane
            distance = (vert_coord - Vector(p1)).dot(line_vector)
            
            # Check if the vertex is within the desired distance from the plane
            if abs(distance) < 0.05:
                intersection_verts.append(vert_coord)
                intersection_indices.append(vert.index)

    if intersection_verts: 
        # Create a new BMesh for each daughter mesh
        bm1 = bmesh.new()
        bm2 = bmesh.new()
        
        # Loop through the faces in the original mesh and add them to the appropriate 
        # daughter mesh
        for face in mesh.polygons:
            # Check if any of the vertices in the face are intersection vertices
            has_intersection = False
            for vert_index in face.vertices:
                if vert_index in intersection_indices:
                    has_intersection = True
                    break
            
            # Add the face to the appropriate daughter mesh
            if has_intersection:
                bm1.faces.new(
                    [bm1.verts[vert] for vert in [v.index for v in face.verts]]
                    )
            else:
                bm2.faces.new(
                    [bm2.verts[vert] for vert in [v.index for v in face.verts]]
                    )
        
        # Finish the BMeshes
        bm1.normal_update()
        bm1.to_mesh(obj.data)
        bm1.free()
        bm2.normal_update()
        bm2.to_mesh(obj.data)
        bm2.free()
    
    else:
        print("No intersection vertices found.")

    return


# not used in division   
def get_division_vertices(obj): 

    print(f"Dividing: {obj.name}")
    
    # Define the two tuple coordinates that define the line
    p1 = get_centerofmass(obj)
    p2 = get_long_axis_global(obj)
    print(f'COM: {p1}, long axis: {p2}')
    
    # Get the vector representing the line
    line_vector = Vector((p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])).normalized()
    
    # Set the active object and get the evaluated mesh
    bpy.context.view_layer.objects.active = obj
    obj_eval = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh_eval = obj_eval.to_mesh()
    
    # List of intersection vertices and their indices
    intersection_verts = []
    intersection_indices = []
    
    # Get the edit-mode mesh and create a BMesh
    mesh = obj.data
    bm = bmesh.new()
    bm.from_mesh(mesh)

    # Ensure that the index table of the BMesh is up to date
    bm.verts.ensure_lookup_table()
    
    # Loop through the vertices in the evaluated mesh
    for vert_eval in mesh_eval.vertices:
        if vert_eval.index < len(mesh.vertices):
            # Get the vertex in the edit-mode mesh using the index
            vert = bm.verts[vert_eval.index]
            
            # Transform the vertex coordinate to global coordinates
            vert_coord = obj.matrix_world @ vert.co
            
            # Calculate the distance of the vertex from the plane
            distance = (vert_coord - Vector(p1)).dot(line_vector)
            
            # Check if the vertex is within the desired distance from the plane
            if abs(distance) < 0.05:
                intersection_verts.append(vert_coord)
                intersection_indices.append(vert.index)

    print(f"Intersect coord: {len(intersection_verts)} vertices, {intersection_verts}")

    # Free the BMesh and the evaluated mesh
    bm.free()
    obj_eval.to_mesh_clear()
    
    return intersection_verts, intersection_indices, p1, line_vector


def add_vertices(obj):
    
    # Define the two tuple coordinates that define the line
    p1 = get_centerofmass(obj)
    p2 = get_long_axis_global(obj)
    print(f'COM: {p1}, long axis: {p2}')
    
    # Get the vector representing the line
    line_vector = Vector((p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])).normalized()
    
    # Set the active object and get the evaluated mesh
    bpy.context.view_layer.objects.active = obj
    obj_eval = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh_eval = obj_eval.to_mesh()
    
    # List of intersection vertices and their indices
    intersection_verts = []
    intersection_indices = []
    
    # Get the edit-mode mesh and create a BMesh
    mesh = obj.data
    bm = bmesh.new()
    bm.from_mesh(mesh)

    # Ensure that the index table of the BMesh is up to date
    bm.verts.ensure_lookup_table()
    
    # Loop through the vertices in the evaluated mesh
    for vert_eval in mesh_eval.vertices:
        if vert_eval.index < len(mesh.vertices):
            # Get the vertex in the edit-mode mesh using the index
            vert = bm.verts[vert_eval.index]
            
            # Transform the vertex coordinate to global coordinates
            vert_coord = obj.matrix_world @ vert.co
            
            # Calculate the distance of the vertex from the plane
            distance = (vert_coord - Vector(p1)).dot(line_vector)
            
            # Check if the vertex is within the desired distance from the plane
            if abs(distance) < 0.01:
                intersection_verts.append(vert_coord)
                intersection_indices.append(vert.index)

    print(f"Intersect coord: {len(intersection_verts)} vertices, {intersection_verts}")

    # Free the BMesh and the evaluated mesh
    bm.free()
    obj_eval.to_mesh_clear()

    # Switch to edit mode and select the intersection vertices
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type="VERT")

    # Retrieve the linked edges for the vertex object
    linked_edges = [e for e in mesh.edges if intersection_indices[0] in e.vertices]
    print(f"Linked edges index: {[edge.index for edge in linked_edges]}")
    edges_target = []

    # Get the vector representing the direction of the edge
    for edge in linked_edges: 
        edge_vector = (
            mesh.vertices[edge.vertices[1]].co - mesh.vertices[edge.vertices[0]].co
        ).normalized()
        # Get the vector representing the direction orthogonal to the division plane
        ortho_vector = line_vector.cross(edge_vector)
        # Get the vector representing the direction orthogonal to the division plane
        if abs(edge_vector.dot(ortho_vector)) < 0.01:
            edges_target.append(edge)
        else:
            print('No edge to further subdivide')
    print(f"Edges index being refined: {edges_target}")

    for edge in edges_target: 
        bpy.ops.mesh.loopcut_slide(
            MESH_OT_loopcut={"number_cuts": 3,
                             "smoothness": 0,
                             "falloff": 'INVERSE_SQUARE',
                             "object_index": 0,
                             "edge_index": edge.index,
                             "mesh_select_mode_init": (True, False, False)
                             },
            TRANSFORM_OT_edge_slide={
                "value": 0,
                "single_side": False,
                "use_even": False,
                "flipped": False,
                "use_clamp": True,
                "mirror": True,
                "snap": False,
                "snap_elements": {'INCREMENT'},
                "use_snap_project": False,
                "snap_target": 'CLOSEST',
                "use_snap_self": True,
                "use_snap_edit": True,
                "use_snap_nonedit": True,
                "use_snap_selectable_only": False,
                "use_snap_to_same_target": False,
                "snap_face_nearest_steps": 1,
                "snap_point": (0, 0, 0),
                "snap_align": False,
                "snap_normal": (0, 0, 0),
                "correct_uv": True,
                "release_confirm": True,
                "use_accurate": False})

    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

    # Set the active object and get the evaluated mesh
    bpy.context.view_layer.objects.active = obj
    obj_eval = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh_eval = obj_eval.to_mesh()
    
    # List of intersection vertices and their indices
    intersection_verts = []
    intersection_indices = []
    
    # Get the edit-mode mesh and create a BMesh
    mesh = obj.data
    bm = bmesh.new()
    bm.from_mesh(mesh)

    # Ensure that the index table of the BMesh is up to date
    bm.verts.ensure_lookup_table()
    
    # Loop through the vertices in the evaluated mesh
    for vert_eval in mesh_eval.vertices:
        if vert_eval.index < len(mesh.vertices):
            # Get the vertex in the edit-mode mesh using the index
            vert = bm.verts[vert_eval.index]
            
            # Transform the vertex coordinate to global coordinates
            vert_coord = obj.matrix_world @ vert.co
            
            # Calculate the distance of the vertex from the plane
            distance = (vert_coord - Vector(p1)).dot(line_vector)
            
            # Check if the vertex is within the desired distance from the plane
            if abs(distance) < 0.01:
                intersection_verts.append(vert_coord)
                intersection_indices.append(vert.index)
    
    # Free the BMesh and the evaluated mesh
    bm.free()
    obj_eval.to_mesh_clear()

    # Sort the indices to ensure that the vertices are selected in order
    intersection_indices.sort()

    return intersection_indices, line_vector


def constrict(obj, indices): 

    # Define the two tuple coordinates that define the line
    p1 = get_centerofmass(obj)
    p2 = get_long_axis_global(obj)
    print(f'COM: {p1}, long axis: {p2}')
    
    # Get the vector representing the line
    # line_vector = Vector((p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])).normalized()
    
    # Set the active object and get the evaluated mesh
    bpy.context.view_layer.objects.active = obj
    # obj_eval = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
    # mesh_eval = obj_eval.to_mesh()
    
    # List of intersection vertices and their indices
    # intersection_verts = []
    # intersection_indices = []
    
    # Get the edit-mode mesh and create a BMesh
    mesh = obj.data
    # bm = bmesh.new()
    # bm.from_mesh(mesh)

    # Ensure that the index table of the BMesh is up to date
    # bm.verts.ensure_lookup_table()
    
    '''# Loop through the vertices in the evaluated mesh
    for vert_eval in mesh_eval.vertices:
        if vert_eval.index < len(mesh.vertices):
            # Get the vertex in the edit-mode mesh using the index
            vert = bm.verts[vert_eval.index]
            
            # Transform the vertex coordinate to global coordinates
            vert_coord = obj.matrix_world @ vert.co
            
            # Calculate the distance of the vertex from the plane
            distance = (vert_coord - Vector(p1)).dot(line_vector)
            
            # Check if the vertex is within the desired distance from the plane
            if abs(distance) < 0.03:
                intersection_verts.append(vert_coord)
                intersection_indices.append(vert.index)

    print(f"Intersect coord: {len(intersection_verts)} vertices, {intersection_verts}")

    
    # Free the BMesh and the evaluated mesh
    bm.free()
    #bpy.data.meshes.remove(mesh_eval)
    obj_eval.to_mesh_clear()'''
    
    # Set the active object and get the evaluated mesh
    # bpy.context.view_layer.objects.active = obj
    # mesh = bpy.context.active_object.data

    # Switch to edit mode and select the intersection vertices
    # bpy.ops.object.mode_set(mode='EDIT')
    # bpy.ops.mesh.select_mode(type="VERT")

    '''first_vert = mesh.vertices[intersection_indices[0]]
    # Retrieve the linked edges for the vertex object
    linked_edges = [e for e in mesh.edges if intersection_indices[0] in e.vertices]
    print([edge.index for edge in linked_edges])
    edges_target = []

    # Get the vector representing the direction of the edge
    for edge in linked_edges:
        edge_vector = (
            mesh.vertices[edge.vertices[1]].co - mesh.vertices[edge.vertices[0]].co
            ).normalized()
        # Get the vector representing the direction orthogonal to the division plane
        ortho_vector = line_vector.cross(edge_vector)
        # Get the vector representing the direction orthogonal to the division plane
        if abs(edge_vector.dot(ortho_vector)) < 0.1:
            edges_target.append(edge)
        else:
            print('No edge to further subdivide')

    for edge in edges_target:
        bpy.ops.mesh.loopcut_slide(
            MESH_OT_loopcut={"number_cuts": 1,
                             "smoothness": 0,
                             "falloff": 'INVERSE_SQUARE',
                             "object_index": 0,
                             "edge_index": edge.index,
                             "mesh_select_mode_init": (True, False, False)
                             },
            TRANSFORM_OT_edge_slide={
                "value": 0,
                "single_side": False,
                "use_even": False,
                "flipped": False,
                "use_clamp": True,
                "mirror": True,
                "snap": False,
                "snap_elements": {'INCREMENT'},
                "use_snap_project": False,
                "snap_target": 'CLOSEST',
                "use_snap_self": True,
                "use_snap_edit": True,
                "use_snap_nonedit": True,
                "use_snap_selectable_only": False,
                "use_snap_to_same_target": False,
                "snap_face_nearest_steps": 1,
                "snap_point": (0, 0, 0),
                "snap_align": False,
                "snap_normal": (0, 0, 0),
                "correct_uv": True,
                "release_confirm": True,
                "use_accurate": False})'''

    # bpy.ops.mesh.select_all(action='DESELECT')
    # bpy.ops.object.mode_set(mode = 'OBJECT')

    '''# Set the active object and get the evaluated mesh
    bpy.context.view_layer.objects.active = obj
    obj_eval = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh_eval = obj_eval.to_mesh()
    
    # List of intersection vertices and their indices
    intersection_verts = []
    intersection_indices = []
    
    # Get the edit-mode mesh and create a BMesh
    mesh = obj.data
    bm = bmesh.new()
    bm.from_mesh(mesh)

    # Ensure that the index table of the BMesh is up to date
    bm.verts.ensure_lookup_table()
    
    # Loop through the vertices in the evaluated mesh
    for vert_eval in mesh_eval.vertices:
        if vert_eval.index < len(mesh.vertices):
            # Get the vertex in the edit-mode mesh using the index
            vert = bm.verts[vert_eval.index]
            
            # Transform the vertex coordinate to global coordinates
            vert_coord = obj.matrix_world @ vert.co
            
            # Calculate the distance of the vertex from the plane
            distance = (vert_coord - Vector(p1)).dot(line_vector)
            
            # Check if the vertex is within the desired distance from the plane
            if abs(distance) < 0.05:
                intersection_verts.append(vert_coord)
                intersection_indices.append(vert.index)
    
    # Free the BMesh and the evaluated mesh
    bm.free()
    #bpy.data.meshes.remove(mesh_eval)
    obj_eval.to_mesh_clear()'''

    # Sort the indices to ensure that the vertices are selected in order
    # intersection_indices.sort()
    for index in indices:
        mesh.vertices[index].select = True

    bpy.ops.object.mode_set(mode='EDIT') 

    '''bpy.ops.transform.shrink_fatten(value=-rate,
                                    use_even_offset=False,
                                    mirror=True,
                                    use_proportional_edit=False,
                                    proportional_edit_falloff='SMOOTH',
                                    proportional_size=1,
                                    use_proportional_connected=False,
                                    use_proportional_projected=False,
                                    release_confirm=True)'''
    
    bpy.ops.transform.resize(value=(0.1, 0.1, 0.1),
                             orient_type='GLOBAL',
                             orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                             orient_matrix_type='GLOBAL', mirror=True,
                             use_proportional_edit=False,
                             proportional_edit_falloff='SHARP',
                             proportional_size=1, use_proportional_connected=False,
                             use_proportional_projected=False)

    # Update the view
    bpy.ops.mesh.reveal()
    bpy.ops.object.mode_set(mode='OBJECT')
    # bpy.ops.object.mode_set(mode='EDIT')
    # bpy.ops.object.mode_set(mode='OBJECT')
    bpy.context.view_layer.update()


def get_division_angles(cell, alpha): 

    # Get the active object
    cell = bpy.context.active_object

    # Define the two tuple coordinates that define the line
    p1 = get_centerofmass(cell)
    p2 = get_long_axis_global(cell)

    '''# Calculate the direction vectors of the global axes
    x_dir = (Vector((1, 0, 0)) @ cell.matrix_world).normalized()
    y_dir = (Vector((0, 1, 0)) @ cell.matrix_world).normalized()
    z_dir = (Vector((0, 0, 1)) @ cell.matrix_world).normalized()'''

    # Get the vector representing the line
    line_vector = Vector((p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])).normalized()

    '''# Calculate the angles between the line and each axis
    x_angle = math.acos(line_vector.dot(x_dir))
    y_angle = math.acos(line_vector.dot(y_dir))
    z_angle = math.acos(line_vector.dot(z_dir))'''

    '''# Print the results
    print("COM:", p1)
    print("Long axis local:", p2)
    print("X angle:", math.degrees(x_angle))
    print("Y angle:", math.degrees(y_angle))
    print("Z angle:", math.degrees(z_angle))'''

    # Calculate a vector perpendicular to the line direction
    # perp_dir = Vector((-line_vector.y, line_vector.x, 0)).normalized()

    # Calculate a second perpendicular vector
    if abs(line_vector.z) < 0.999:
        up_dir = line_vector.cross(Vector((0, 0, 1)))
    else:
        up_dir = line_vector.cross(Vector((0, 1, 0)))
    up_dir.normalize()

    # Calculate a third perpendicular vector
    side_dir = line_vector.cross(up_dir)
    side_dir.normalize()
    # Calculate the quaternion that rotates the plane to align with the line direction
    rotation_matrix = Matrix((side_dir, up_dir, line_vector)).transposed()
    quat = rotation_matrix.to_quaternion()
    # Create the plane primitive
    '''plane = bpy.ops.mesh.primitive_plane_add(size=2,
                                             enter_editmode=False,
                                             location=p1)'''
    # Get a reference to the created plane object
    plane_obj = bpy.context.active_object
    # Set the plane's rotation
    plane_obj.rotation_mode = 'QUATERNION'
    plane_obj.rotation_quaternion = quat

    # Set the active object and get the evaluated mesh
    bpy.context.view_layer.objects.active = cell
    cell_eval = cell.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh_eval = cell_eval.to_mesh()

    # List of intersection vertices
    intersection_verts = []

    for vert in mesh_eval.vertices:
        if vert.index < len(cell.data.vertices):
            # Transform the vertex coordinate to global coordinates
            vert_coord = cell.matrix_world @ vert.co
            # Calculate the distance of the vertex from the plane
            distance = (vert_coord - Vector(p1)).dot(line_vector)
            # Check if the vertex is within the desired distance from the plane
            if abs(distance) < 0.05:
                intersection_verts.append(vert_coord)
                # Move the vertex towards the mesh center of mass
                vert_coord = (1 - alpha) * vert_coord + alpha * Vector(p1)
                # Transform the vertex coordinate back to local coordinates
                local_vert_coord = cell.matrix_world.inverted() @ vert_coord
                # Update the vertex coordinate in the mesh data
                cell.data.vertices[vert.index].co = local_vert_coord
    
    '''# Set the active object and get the evaluated mesh
    bpy.context.view_layer.objects.active = cell
    cell_eval = cell.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh_eval = cell_eval.to_mesh()

    # List of intersection vertices
    intersection_verts = []

    alpha = 0.99
    frames = 10  # Number of frames over which to smoothen the movement

    for vert in mesh_eval.vertices:
        if vert.index < len(cell.data.vertices):
            # Transform the vertex coordinate to global coordinates
            vert_coord = cell.matrix_world @ vert.co
            # Calculate the distance of the vertex from the plane
            distance = (vert_coord - Vector(p1)).dot(line_vector)
            # Check if the vertex is within the desired distance from the plane
            if abs(distance) < 0.05:
                intersection_verts.append(vert_coord)
                # Calculate the step size for each frame
                step_size = (Vector(p1) - vert_coord) / frames
                # Move the vertex towards the center of mass over the specified frames
                for i in range(frames):
                    vert_coord += step_size
                    # Transform the vertex coordinate back to local coordinates
                    local_vert_coord = cell.matrix_world.inverted() @ vert_coord
                    # Update the vertex coordinate in the mesh data
                    cell.data.vertices[vert.index].co = local_vert_coord'''

    print(f"Intersect coord: {len(intersection_verts)} vertices, {intersection_verts}")

    '''cell_mesh = cell.data

    # Create two lists to store the vertices for the two new meshes
    mesh1_verts = []
    mesh2_verts = []

    # Iterate over the vertices in the original mesh and 
    # add them to the appropriate list
    for vert in cell_mesh.vertices:
        if vert.co in intersection_verts:
            # If the vertex is an intersection vertex, add it to both lists
            mesh1_verts.append(vert.co)
            mesh2_verts.append(vert.co)
        elif vert.co[2] > p1[2]:
            # If the vertex is above the line, add it to mesh1_verts
            mesh1_verts.append(vert.co)
        else:
            # Otherwise, add it to mesh2_verts
            mesh2_verts.append(vert.co)

    # Create two new mesh objects
    mesh1 = bpy.data.meshes.new("Mesh1")
    mesh2 = bpy.data.meshes.new("Mesh2")

    # Create two new object objects, one for each mesh
    obj1 = bpy.data.objects.new("Object1", mesh1)
    obj2 = bpy.data.objects.new("Object2", mesh2)

    # Link the meshes to their respective objects and add the objects to the scene
    obj1.data = mesh1
    obj2.data = mesh2
    bpy.context.scene.collection.objects.link(obj1)
    bpy.context.scene.collection.objects.link(obj2)

    d1_com = get_centerofmass(obj1)
    d2_com = get_centerofmass(obj2)

    bpy.data.objects.remove(obj1, do_unlink=True)
    bpy.data.objects.remove(obj2, do_unlink=True)

    # Calculate the center of mass for each new mesh
    #d1_com = sum(mesh1_verts, Vector()) / len(mesh1_verts)
    #d2_com = sum(mesh2_verts, Vector()) / len(mesh2_verts)

    make_cell(name=f"{cell.name}_d1", 
              loc=d1_com, 
              collection=cell.users_collection[0].name)
    make_force(force_name=f"force_{cell.name}_d1", 
               cell_name=f"{cell.name}_d1", 
               strength=-1000, 
               falloff=1, 
               collection="A_Forces")
    make_cell(name=f"{cell.name}_d2", 
              loc=d2_com, 
              collection=cell.users_collection[0].name)
    make_force(force_name=f"force_{cell.name}_d2", 
               cell_name=f"{cell.name}_d2", 
               strength=-1000, 
               falloff=1, 
               collection="A_Forces")

    for force in Force.force_list: 
        if force.associated_cell == cell.name: 
            bpy.data.objects.remove(bpy.data.objects[force.name], do_unlink=True)

    # Delete the original object
    bpy.data.objects.remove(cell, do_unlink=True)'''

    # Create new vertices from intersection coordinates
    mesh = bpy.data.meshes.new('Intersection')
    obj = bpy.data.objects.new('Intersection', mesh)
    bpy.context.scene.collection.objects.link(obj)

    # Set the mesh data
    mesh.from_pydata(intersection_verts, [], [])
    mesh.update()

    # Create an empty object at p1
    # empty1 = bpy.ops.object.empty_add(type='PLAIN_AXES', location=p1)
    # Create an empty object at p2
    # empty2 = bpy.ops.object.empty_add(type='PLAIN_AXES', location=p2)
    # Create a new curve object and add a new spline to it
    curve_data = bpy.data.curves.new(name='MyCurve', type='CURVE')
    curve_data.dimensions = '3D'
    curve_data.resolution_u = 2
    curve_spline = curve_data.splines.new(type='POLY')
    curve_spline.points.add(2)
    # Set the points of the spline to the line points
    curve_spline.points[0].co = Vector(p1).to_4d()
    curve_spline.points[1].co = Vector(p2).to_4d()
    # Set the line thickness
    curve_data.bevel_depth = 0.03
    # Create a new object with the curve data
    curve_object = bpy.data.objects.new(name='MyCurveObject', object_data=curve_data)
    # Add the object to the scene
    bpy.context.scene.collection.objects.link(curve_object)  


class handler_class:
    # TODO document this class
    """
    A class for creating different types of handlers that trigger
    actions on ``Goo`` cells when certain criteria are met
    """

    # The initialization function specifies available cell types and associated
    # parameters like division rate, growth rate, and adhesion forces
    def __init__(self):

        self.seed = None
        self.cell_types = ['sphere', 'type1', 'type2']
        self.division_rates = {}
        self.growth_rates = {}
        self.adhesion_forces = {} 

        # Set parameter values for each cell type
        for type in self.cell_types:
            self.division_rates[type] = 0
            self.growth_rates[type] = 0
            self.adhesion_forces[type] = {}
            for i in self.cell_types:
                self.adhesion_forces[type][i] = 0

        # Set active (dividing) cell types
        self.active_cell_types = []

        # For adhesion handler
        self.forces = []

        # For motion handler
        self.random_motion_speed = 0

        # Data handler: for cell deformability 
        self.all_vertices = [[]]
        self.COMs = {}
        # Data handler: for total distance between cells - simulation stability
        self.frames = []
        self.distances_tot = []
        self.data_file_path = ''
        self.time = None
        self.times = defaultdict()
        self.absolute_time = [0]
        self.frame_interval = [None, None]
        self.strength = None
        self.falloff = None
        self.master_dict = None
        self.data_dict = {
            'Frames': [],
            'Distances': [],
            'Times': [],
            'Deformability': {},
            'Contact area': [],
            'Axis direction': [],
            'Axis length': [],
            'Tension': [],
            'Adhesion': [],
            'Volume': {}
        }
        # Data handler: for computational cost
        self.cell_number = 0
        # Data handler: for deformability
        self.displacement = defaultdict(list)
        self.previous_vertices = defaultdict(list)
        self.current_vertices = defaultdict(list)
        self.previous_com = defaultdict(list)
        self.sphere_centered = defaultdict(list)
        self.deformability = defaultdict(list)

        # Data handler: for phase diagrams
        self.stable_measure = []
        self.tension = []
        self.adhesion = []

        # For long axis
        self.vec_axis = defaultdict(list)
        self.len_axis = defaultdict(list)

        # Data handler: for contact area 
        self.contact_ratios = defaultdict(list)
        self.contact_areas = defaultdict(list)

        # For cell division
        self.division_rate = 0.0
        self.division_indices = defaultdict(list)
        self.cell_under_div = None
        self.daugthers = []

        # For cell growth 
        self.volumes = defaultdict(list)

        # For random motion
        self.seed = int()
        self.prev_frame = int()
        self.sorting_scores = defaultdict(dict)
        self.msd = defaultdict(list)
        self.speed = defaultdict(list)
        self.motion_path = defaultdict(list)

        # Data handler: flag in launch simulation
        self.data_flag = None

        return

    def launch_simulation(
            self,
            filepath,
            division_rate=100,
            start=1,
            end=250,
            motion_strength=-500,
            adhesion=True,
            growth=False,
            target=1,
            division=False,
            data=False,
            motility=False,
            remeshing=False,
            visualize=False,
            boundary=False):

        self.frame_interval = [start, end]
        bpy.context.scene.frame_set(start)
        bpy.context.scene.frame_start = start
        bpy.context.scene.frame_end = end

        self.data_file_path = filepath
        bpy.context.scene.render.filepath = filepath
        self.division_rate = division_rate
        self.motion_strength = motion_strength
        self.data_flag = data  # used to decide if data are computed

        # Set the end frame for all cloth simulation caches in the scene
        # To keep simulations running after 250 frames
        if bpy.context.scene.frame_current == start: 
            for collection in bpy.data.collections:
                # Loop through the objects existed in the collection 
                for obj in collection.objects:     
                    if obj.get('object') == 'cell':
                        obj.modifiers["Cloth"].point_cache.frame_start = start   
                        obj.modifiers["Cloth"].point_cache.frame_end = end   

        bpy.app.handlers.frame_change_pre.append(self.timing_init_handler)
        bpy.app.handlers.frame_change_post.clear()

        if adhesion:
            bpy.app.handlers.frame_change_post.append(self.adhesion_handler)
        if data:
            bpy.app.handlers.frame_change_post.append(self.data_export)
            bpy.app.handlers.frame_change_post.append(self.contact_area_handler)
        if growth:
            bpy.app.handlers.frame_change_post.append(
                lambda scene,
                depsgraph: self.growth_handler(scene,
                                               depsgraph,
                                               target)
                                               )
        if division:
            bpy.app.handlers.frame_change_post.append(self.division_handler)
        if motility:
            bpy.app.handlers.frame_change_post.append(self.motion_handler)
        if boundary:
            bpy.app.handlers.frame_change_post.append(self.boundary_handler)
        if remeshing:
            bpy.app.handlers.frame_change_post.append(self.remeshing_handler)
            # bpy.app.handlers.frame_change_post.append(self.select_dividing_cell)
        if visualize:
            bpy.app.handlers.frame_change_post.append(self.visualize_stretching)

        bpy.app.handlers.frame_change_post.append(self.timing_elapsed_handler)
        bpy.app.handlers.frame_change_post.append(self.stop_animation)

        # bpy.ops.screen.animation_play()

        return 

    # Member function to set division rate for a cell type
    def set_division_rate(self, cell_type, rate):
        """
        Sets division rate in the handler class that div_handler() can reference later
        :param cell_type: Name of cell type to apply this division rate to.
        Must be one of the active cell types. (String)
        :param rate: number of frames between each division (int)
        :return: None
        """

        self.division_rates[cell_type] = rate
        return

    # Member function to set growth rate for a cell type
    def set_growth_rate(self, cell_type, rate):
        """
        Sets rate rate in the handler class that growth_handler() can reference later
        :param cell_type: Name of cell type to apply this division rate to.
        Must be one of the active cell types. (String)
        :param rate: amount to change cell.modifiers["Cloth"].
        settings.shrink_min each frame.
        Should be between 0 and 1. (float)
        :return: None
        """

        if (rate > 0 and rate < 1):
            self.growth_rates[cell_type] = rate

        return

    # Member function to set adhesion forces between cell types 
    # currently not used
    def set_adhesion(self, type1, type2, force):
        """
        Sets a value adhesion_forces in the handler class that
        apply_forces() can reference later
        :param type1: Name of cell type that the force is attatched to.
        Must be one of the active cell types. (String)
        :param type1: Name of cell type that the force affects.
        Must be one of the active cell types. (String)
        :param force: the strenfgth of the force (int)
        :return: None
        """

        self.adhesion_forces[type1][type2] = force

        return

    # currently not used
    def apply_forces(self):
        """
        Add force fields to force collections and make them affect corresponding cell
        types
        :return: None
        """

        master_collection = bpy.context.view_layer.active_layer_collection
        for cell_type in self.active_cell_types:
            num_cells = len(bpy.data.collections[cell_type].objects)
            for i in range(num_cells):
                cell_name = bpy.data.collections[cell_type].objects[i].name
                cell = bpy.data.objects[cell_name]
                cell.modifiers["Cloth"].settings.effector_weights.collection = \
                    bpy.data.collections[cell_type+"_forces"]
                for affected_type in self.cell_types:
                    if self.adhesion_forces[cell_type][affected_type] != 0:
                        vl = bpy.context.view_layer
                        affected = vl.layer_collection.children[affected_type+"_forces"]
                        bpy.context.view_layer.active_layer_collection = affected
                        f = Force(cell_name+"_to_"+affected_type,
                                  cell_name,
                                  self.adhesion_forces[cell_type][affected_type])
                        make_force(f)
                        vl.active_layer_collection = master_collection
                        self.forces.append(f)

    '''def select_dividing_cell(self, scene, depsgraph): 

        for collection in bpy.data.collections: 
            # Exclude the objects in the force collections
            if collection['type'] == 'cell':
    
                # Loop through the objects existed in the collection 
                print(f"Collection under division: {collection.name_full}")
                cells = bpy.data.collections.get(collection.name_full).all_objects
                #cells.append(None)
                print(f"Cells under division: {[obj.name for obj in cells]}")

                # randomly choose a cell to divide
                self.cell_under_div = random.choice(cells)'''

    def remeshing_handler(self, scene, depsgraph): 
        remesh_frames = range(self.frame_interval[0], self.frame_interval[1], 1)[1:]
        # div_frames = range(self.frame_interval[0], self.frame_interval[1], 1)[1:]
        if scene.frame_current in remesh_frames:
            for collection in bpy.data.collections:
                # Exclude the objects in the force collections
                if collection['type'] == 'cell':
                    cells = bpy.data.collections.get(collection.name_full).all_objects
                    for cell in cells:
                        cell.select_set(True)
                    for cell in bpy.context.selected_objects:
                        cell.modifiers.new(
                            name=f"Remesh_tmp_{cell.name}",
                            type='REMESH'
                            )
                        remesh_mod = cell.modifiers.get(f"Remesh_tmp_{cell.name}")
                        remesh_mod.mode = 'SMOOTH'
                        remesh_mod.octree_depth = 3
                        remesh_mod.scale = 0.8
                        remesh_mod.use_remove_disconnected = True
                        remesh_mod.use_smooth_shade = False
                        remesh_mod.show_in_editmode = True
                        remesh_mod.show_in_editmode = True
                        bpy.ops.object.modifier_move_to_index(
                            modifier=f"Remesh_tmp_{cell.name}",
                            index=-1
                            )
                        bpy.context.view_layer.objects.active = \
                            bpy.data.objects[cell.name]      
                        bpy.ops.object.modifier_apply(
                            modifier=f"Remesh_tmp_{cell.name}"
                            )
                        # deselect all vertices in edit mode
                        bpy.ops.object.mode_set(mode='EDIT')
                        bpy.ops.mesh.select_mode(type="VERT")
                        bpy.ops.mesh.reveal()
                        bpy.ops.object.mode_set(mode='OBJECT')

    def division_handler(self, scene, despgraph):

        # rates = np.arange(0, 1, 0.05, dtype=np.float)
        div_frames = range(
            self.frame_interval[0],
            self.frame_interval[1],
            self.division_rate
            )[1:]
        '''sphere_frames = [
            num for start in div_frames 
            for num in range(start+5, start+len(rates))
            ]'''
        div_frames_physics = range(
            self.frame_interval[0] + 0, 
            self.frame_interval[1], 
            self.division_rate
            )[1:]
        '''div_frames_sphere = range(
            self.frame_interval[0] + 5, 
            self.frame_interval[1], 
            self.division_rate
            )[1:]'''
            
        if scene.frame_current in div_frames: 
            self.daugthers.clear()
            self.strength = None
            self.force_collection = None
            self.falloff = None

            for collection in bpy.data.collections: 
                # Exclude the objects in the force collections
                # if collection['type'] == 'cell':
                # Loop through the objects existed in the collection 
                print(f"Collection under division: {collection.name_full}")
                cells = [
                    obj for obj 
                    in bpy.data.collections.get(collection.name_full).all_objects 
                    if obj.get('object') == 'cell'
                    ]
                '''for cell in cells: 
                    print(cell.get('object'))
                    if cell.get('object') and cell['object'] == 'cell': 
                        #cells.append(None)'''
                print(f"Cells under division: {[obj.name for obj in cells]}")

                # randomly choose a cell to divide
                self.cell_under_div = random.choice(cells)

            if self.cell_under_div is not None: 
                # for obj in cells:   
                cell_name = self.cell_under_div.name
                print(f"-- Starting division of {cell_name} "
                      f"at frame {scene.frame_current}")
                d1, d2, tmp_strength, tmp_collection, tmp_falloff, mother_modifiers = \
                    divide_boolean(self.cell_under_div) 
                print(f"{tmp_strength}; {tmp_collection}; {tmp_falloff}")

                self.strength = tmp_strength
                self.force_collection = tmp_collection
                self.falloff = tmp_falloff
                self.modifiers = mother_modifiers
                print([mod for mod in mother_modifiers])

                self.daugthers.append(d1)     
                self.daugthers.append(d2)      

                # bpy.context.scene.frame_set(1)

                print(f"-- Finishing division of {cell_name} "
                      f"at frame {scene.frame_current}")

            else: 
                print('No cells under division')

        if scene.frame_current in div_frames_physics: 

            apply_physics(self.daugthers[0], self.modifiers)
            apply_physics(self.daugthers[1], self.modifiers)

            self.daugthers[0].modifiers["Cloth"].point_cache.frame_start = \
                self.frame_interval[0]   
            self.daugthers[0].modifiers["Cloth"].point_cache.frame_end = \
                self.frame_interval[1]  
            self.daugthers[1].modifiers["Cloth"].point_cache.frame_start = \
                self.frame_interval[0]   
            self.daugthers[1].modifiers["Cloth"].point_cache.frame_end = \
                self.frame_interval[1]  

            apply_daugther_force(
                self.daugthers[0], 
                self.strength, 
                self.force_collection, 
                self.falloff)
            apply_daugther_force(
                self.daugthers[1], 
                self.strength, 
                self.force_collection, 
                self.falloff)

        # if scene.frame_current in div_frames_sphere: 
            
        #    bpy.context.view_layer.objects.active = self.daugthers[0]
        #    to_sphere(self.daugthers[0], 1)

        #    bpy.context.view_layer.objects.active = self.daugthers[1]
        #    to_sphere(self.daugthers[1], 1)

            # remesh(self.daugthers[0])
            # remesh(self.daugthers[1])

        # CODE FOR TO SPHERE IMPLEMENTATION
        '''div_frames_sphere = range(
            self.frame_interval[0] + 2, 
            self.frame_interval[1], 
            self.division_rate
            [1:]
        frames = range(div_frames[0], div_frames[0] + len(rates), 1)

        if scene.frame_current in sphere_frames: 

            bpy.context.view_layer.objects.active = self.daugthers[0]
            to_sphere(
                self.daugthers[0], 
                rates[sphere_frames.index(scene.frame_current)]
                )
            #bpy.context.view_layer.objects.active = None

            bpy.context.view_layer.objects.active = self.daugthers[1]
            to_sphere(
                self.daugthers[1], 
                rates[sphere_frames.index(scene.frame_current)]
                )'''

    def get_division_vertices(self, scene, depsgraph): 

        cell = bpy.data.objects['cell_A1']

        if scene.frame_current == 100:
            print(f"Long axis: np {get_long_axis_np(cell)}; "
                  f"global {get_long_axis_global(cell)}")
            get_division_angles(cell, 0.9)

        '''if scene.frame_current == 102:
            print(f"Long axis: np {get_long_axis_np(cell)}; "
                  f"global {get_long_axis_global(cell)}")
            get_division_angles(cell, 0.6)

        if scene.frame_current == 104:
            get_division_angles(cell, 0.8)'''

        '''obj = bpy.data.objects['cell_A1']

        # Define the long axis as a vector
        # Change this vector to match your desired long axis
        long_axis = mathutils.Vector((0, 0, 2)) 

        # Define the center of mass as a point in world coordinates
        # Change this vector to match your desired center of mass
        center_of_mass = mathutils.Vector(get_centerofmass(obj)) 

        # Convert the center of mass to mesh-relative coordinates
        mat = obj.matrix_world.inverted()
        center_of_mass_rel = mat @ center_of_mass

        # Create a new plane to represent the division plane
        # Assume up direction is +Z
        normal = long_axis.cross(mathutils.Vector((0, 0, 2))) 
        normal.normalize()
        point = center_of_mass_rel
        division_plane = (point, normal)

        # Add new vertices at the intersection of the division plane and the mesh
        bm = bmesh.new()
        bm.from_mesh(obj.data)

        for edge in bm.edges:
            v1, v2 = edge.verts
            p1 = v1.co.copy()
            p2 = v2.co.copy()

            if (p1 - point).dot(normal) * (p2 - point).dot(normal) <= 0:
                # If the edge intersects the division plane
                intersection_point = mathutils.geometry.intersect_line_plane(p1, 
                                                                             p2, 
                                                                             point, 
                                                                             normal)
                new_vert = bm.verts.new(tuple(intersection_point))

        bm.to_mesh(obj.data)
        bm.free()'''

    def adjust_shrink_min(volume, target_volume):
        volume_deviation = (volume - target_volume) / target_volume
        shrink_adjustment = 0.01 * math.tanh(5 * volume_deviation)
        return shrink_adjustment                                                                             

    # Member function to handle cell growth
    def growth_handler(self, scene, depsgraph, target):
        for collection in bpy.data.collections: 
            # Exclude the objects in the force collections
            # if collection['type'] == 'cell':
            # Loop through the objects existed in the collection 
            cells = bpy.data.collections.get(collection.name_full).all_objects
            for cell in cells: 
                if (
                    cell.get('adhesion force') is not None or
                    cell.get('motion force') is not None or
                    (cell.get('adhesion force') is not None and 
                     cell.get('motion force') is not None)
                ):

                    volume = calculate_volume(cell)
                    target_volume = ((4/3)*np.pi*(1)**3)
                    cell['target volume'] = target_volume
                    cell['volume'] = volume
                    self.volumes[f"{cell.name}"].append(volume)

                    # TODO allow for multiple growth rate, with different functions
                    volume_deviation = (volume - target_volume) / target_volume
                    print(f"Current volume deviation: {volume_deviation}")
                    # shrink_adjustment = 0.001 * math.tanh(50 * volume_deviation)
                    shrink_adjustment = 0.001 * (math.exp(volume_deviation) - 1)
                    # shrink_adjustment = 0.01 * volume_deviation
                    cell.modifiers["Cloth"].settings.shrink_min += shrink_adjustment

                    # volume control
                    '''if volume < (target_volume - target_volume * 0.2): 
                        print(f"Current volume: {volume} is lower"
                              f"than target {target_volume} - 20%")
                        cell.modifiers["Cloth"].settings.shrink_min -= 0.01
                    elif volume < (target_volume - target_volume * 0.05): 
                        print(f"Current volume: {volume} is lower"
                              f"than target {target_volume} - 5%")
                        cell.modifiers["Cloth"].settings.shrink_min -= 0.001
                    elif volume < (target_volume - target_volume * 0.01): 
                        print(f"Current volume: {volume} is lower"
                              f"than target {target_volume} - 1%")
                        cell.modifiers["Cloth"].settings.shrink_min -= 0.0005
                    # shrink
                    elif volume > (target_volume + target_volume * 0.2): 
                        print(f"Current volume: {volume} is higher"
                              f"than target {target_volume} - 20%")
                        cell.modifiers["Cloth"].settings.shrink_min += 0.01
                    elif volume > (target_volume + target_volume * 0.05): 
                        print(f"Current volume: {volume} is higher"
                              f"than target {target_volume} - 5%")
                        cell.modifiers["Cloth"].settings.shrink_min += 0.001
                    elif volume < (target_volume + target_volume * 0.01): 
                        print(f"Current volume: {volume} is higher"
                              f"than target {target_volume} - 1%")
                        cell.modifiers["Cloth"].settings.shrink_min += 0.0005
                    else: 
                        print(f"Current volume: {volume} is within"
                              f"target {target_volume}")'''

    def boundary_handler(self, scene, depsgraph):
        # Get the 'box' object
        box_object = bpy.data.objects.get('box')
        print(box_object)

        if box_object:
            # Get the size and center from the object name
            box_dimensions = box_object.dimensions
            box_center = box_object.location

            # Calculate the boundaries of the box
            x_min = box_center.x - box_dimensions.x / 2
            x_max = box_center.x + box_dimensions.x / 2
            y_min = box_center.y - box_dimensions.y / 2
            y_max = box_center.y + box_dimensions.y / 2
            z_min = box_center.z - box_dimensions.z / 2
            z_max = box_center.z + box_dimensions.z / 2

            for collection in bpy.data.collections: 
                forces = bpy.data.collections.get(collection.name_full).all_objects
                for force in forces: 
                    if (force.get('motion') is not None and force.get('motion')): 

                        constrained_force_location = force.location
                        print(f'Before constrained: {force.location}')

                        # Reflect x-coordinate if it's outside the box boundaries
                        if constrained_force_location.x < x_min:
                            constrained_force_location.x = \
                                2 * x_min - constrained_force_location.x
                        elif constrained_force_location.x > x_max:
                            constrained_force_location.x = \
                                2 * x_max - constrained_force_location.x

                        # Reflect y-coordinate if it's outside the box boundaries
                        if constrained_force_location.y < y_min:
                            constrained_force_location.y = \
                                2 * y_min - constrained_force_location.y
                        elif constrained_force_location.y > y_max:
                            constrained_force_location.y = \
                                2 * y_max - constrained_force_location.y

                        # Reflect z-coordinate if it's outside the box boundaries
                        if constrained_force_location.z < z_min:
                            constrained_force_location.z = \
                                2 * z_min - constrained_force_location.z
                        elif constrained_force_location.z > z_max:
                            constrained_force_location.z = \
                                2 * z_max - constrained_force_location.z

                        force.location = constrained_force_location
                        print(f'After constrained: {force.location}')

                        # Check if the force location is outside the box
                        if not (x_min <= constrained_force_location.x <= x_max and
                                y_min <= constrained_force_location.y <= y_max and
                                z_min <= constrained_force_location.z <= z_max):
                            raise ValueError(
                                "Force location is outside the box boundaries."
                                )

    def motion_handler(self, scene, depsgraph): 
        
        cells = [
            obj for obj in bpy.data.objects 
            if "object" in obj.keys() and obj["object"] == "cell"
            ]
        print(cells)
        msd = dict()        
        
        for collection in bpy.data.collections: 
            forces = bpy.data.collections.get(collection.name_full).all_objects
            print(forces)
            for force in forces: 
                if (force.get('motion') is not None and force.get('motion')): 
                    print('Entering force loop')
                    # self.force_path[f'{force.name}_force_tracks'].append(tuple(force.location))
                    cell = bpy.data.objects[force.get('cell')]
                    com = get_centerofmass(cell)
                    cell['current position'] = com
                    disp = (Vector(cell.get('current position')) - 
                            Vector(cell.get('past position'))).length
                    cell['speed'] = disp
                    self.motion_path[f'{cell.name}'].append(tuple(com)[:3])
                    self.motion_path[f'{force.name}'].append(tuple(force.location)[:3])
                    self.speed[f'{cell.name}'].append(disp)
                    # Mean Squared Displacement for single particle
                    sq_displacement_cell = (
                        Vector(self.motion_path.get(cell.name)[-1]) - 
                        Vector(self.motion_path.get(cell.name)[0])
                        ).length_squared
                    sq_displacement_force = (
                        Vector(self.motion_path.get(force.name)[-1]) - 
                        Vector(self.motion_path.get(force.name)[0])
                        ).length_squared
                    self.msd[f'{cell.name}'].append(sq_displacement_cell)
                    self.msd[f'{force.name}'].append(sq_displacement_force)
                    cell['MSD'] = msd
                    force['MSD'] = msd

                    if force.get('distribution') == 'uniform': 
                        print('Distribution is uniform')
                        rand_coord = Vector(np.random.uniform(
                            low=-force['distribution size'],
                            high=force['distribution size'], 
                            size=(3,)
                        ))
                        new_loc = Vector(com) + rand_coord
                        print(new_loc)
                    elif force.get('distribution') == 'gaussian': 
                        rand_coord = Vector(np.random.normal(
                            loc=0, 
                            scale=force['distribution size'], 
                            size=(3,)
                        ))
                        new_loc = Vector(com) + rand_coord
                    else: 
                        print(f'{force.get("distribution")} '
                              f'is not a supported distribution')
                        continue
                    
                    force.location = new_loc

                    '''# constraint force field within the box

                    # Define the box's dimensions and center
                    # Define the box's length (half of the actual length)
                    box_length = 4.5  
                    box_center = Vector((-1.5, 0, 0))  # Define the center of the box

                    # Calculate the boundaries of the box
                    x_min = box_center.x - box_length
                    x_max = box_center.x + box_length
                    y_min = box_center.y - box_length
                    y_max = box_center.y + box_length
                    z_min = box_center.z - box_length
                    z_max = box_center.z + box_length

                    constrained_force_location = force.location
                    # Constrain the x-coordinate within the box
                    #force.location.x = max(x_min, min(x_max, force.location.x))
                    # Constrain the y-coordinate within the box
                    #force.location.y = max(y_min, min(y_max, force.location.y))
                    # Constrain the z-coordinate within the box
                    #force.location.z = max(z_min, min(z_max, force.location.z))

                    # Reflect x-coordinate if it's outside the box boundaries
                    if constrained_force_location.x < x_min:
                        constrained_force_location.x = \
                            2 * x_min - constrained_force_location.x
                    elif constrained_force_location.x > x_max:
                        constrained_force_location.x = \
                            2 * x_max - constrained_force_location.x

                    # Reflect y-coordinate if it's outside the box boundaries
                    if constrained_force_location.y < y_min:
                        constrained_force_location.y = \
                            2 * y_min - constrained_force_location.y
                    elif constrained_force_location.y > y_max:
                        constrained_force_location.y = \
                            2 * y_max - constrained_force_location.y

                    # Reflect z-coordinate if it's outside the box boundaries
                    if constrained_force_location.z < z_min:
                        constrained_force_location.z = \
                            2 * z_min - constrained_force_location.z
                    elif constrained_force_location.z > z_max:
                        constrained_force_location.z = \
                            2 * z_max - constrained_force_location.z

                    force.location = constrained_force_location'''

                    '''# cell tracks 
                    # Define the name of the curve object and the new point coordinates
                    # Find the curve object by name
                    curve_obj = bpy.data.objects[f'{cell.name}_tracks']

                    # Enter edit mode for the curve object
                    bpy.context.view_layer.objects.active = curve_obj
                    bpy.ops.object.mode_set(mode='EDIT')

                    # Add a new control point to the curve
                    spline = curve_obj.data.splines[0]
                    spline.points.add(1)
                    new_point_index = len(spline.points) - 1
                    spline.points[new_point_index].co = (*com, 1)

                    # Update the curve display and exit edit mode
                    bpy.ops.curve.reveal()
                    bpy.ops.object.mode_set(mode='OBJECT')

                    # force tracks
                    curve_force_obj = bpy.data.objects[f'{force.name}_force_tracks']

                    # Enter edit mode for the curve object
                    bpy.context.view_layer.objects.active = curve_force_obj
                    bpy.ops.object.mode_set(mode='EDIT')

                    # Add a new control point to the curve
                    spline = curve_force_obj.data.splines[0]
                    spline.points.add(1)
                    new_point_index = len(spline.points) - 1
                    spline.points[new_point_index].co = (*force.location, 1)

                    # Update the curve display and exit edit mode
                    bpy.ops.curve.reveal()
                    bpy.ops.object.mode_set(mode='OBJECT')'''

                    # cells will only adhere with other cells 
                    # that are in the same collection
                    obj_cell = bpy.data.objects[force.get('cell')]
                    cloth_modifier = obj_cell.modifiers["Cloth"]
                    cloth_modifier.settings.effector_weights.collection = collection
                    # update position of cell
                    cell['past position'] = com

                else: 
                    print('Motion handler is not working')
        
        '''if scene.frame_current == self.frame_interval[1] - 1:
        # Iterate over all spline objects in the scene
            tracks = [
                obj 
                for collection in bpy.data.collections 
                for obj in collection.objects 
                if obj.type == 'CURVE'
                ]
            
            for obj in tracks:
                    # Initialize variables
                    spline = obj.data.splines[0]
                    points = spline.points

                    # Iterate over control points and calculate displacements
                    for i in range(1, len(points)):
                        intial_point = points[0].co
                        prev_point = points[i-1].co
                        curr_point = points[i].co
                        # Mean Squared Displacement for single particle
                        displacement = (curr_point - intial_point).length_squared
                        # Single particle speed
                        speed = (curr_point - prev_point).length

                        msd = displacement / 1
                        obj['MSD'] = msd

                        self.msd[f'{obj.name}'].append(msd)
                        self.speed[f'{obj.name}'].append(speed)
                        self.motion_path[f'{obj.name}'].append(tuple(curr_point)[:3])'''
        if self.data_flag: 
            # Initialize a dictionary to store the sorting scores for each cell type
            cell_types = np.unique([coll['type'] for coll in bpy.data.collections])
            '''cells = [
                obj for obj in bpy.data.objects 
                if "object" in obj.keys() and obj["object"] == "cell"
                ]'''

            neighbors = {
                cell_type: {
                    cell.name: 0 
                    for cell in cells 
                    if cell.users_collection[0].get('type') == cell_type
                } 
                for cell_type in cell_types
            }

            neighbors_same_type = {
                cell_type: {
                    cell.name: 0 
                    for cell in cells 
                    if cell.users_collection[0].get('type') == cell_type
                } 
                for cell_type in cell_types
            }

            sorting_scores = {
                cell_type: {
                    cell.name: 0 
                    for cell in cells 
                    if cell.users_collection[0].get('type') == cell_type
                } 
                for cell_type in cell_types
            }

            sorting_scores_same_type = {
                cell_type: 0 
                for cell_type in cell_types
            }

            # loop for each cell
            for cell in cells:
                collection = cell.users_collection[0]
                for other_cell in cells:  
                    if cell is not other_cell: 
                        # get neighbors 
                        distance = (
                            bpy.data.objects[cell.get('adhesion force')].location - 
                            bpy.data.objects[other_cell.get('adhesion force')].location
                        ).length
                        if distance < 1.95:
                            # get number of neighbors for a specific cell
                            neighbors[collection.get('type')][cell.name] += 1

                            # get neighbors of same cell type
                            other_type = other_cell.users_collection[0].get('type')
                            if (
                                collection.get('type') is other_type
                            ):
                                neighbor_type = neighbors_same_type[
                                    collection.get('type')
                                    ]
                                neighbor_type[cell.name] += 1
                # sorting is null if cell has no neighbors                 
                if neighbors[collection.get('type')][cell.name] == 0: 
                    sorting_scores[collection.get('type')][cell.name] = 0 
                else: 
                    sorting_scores[collection.get('type')][cell.name] = (
                        neighbors_same_type[collection.get('type')][cell.name] / 
                        neighbors[collection.get('type')][cell.name]
                        )
                cell['sorting score'] = (
                    sorting_scores[collection.get('type')][cell.name]
                )
                print(neighbors)
                print(neighbors_same_type)
            
            # avg sorting score over cells among the same type
            for cell_type, cell_dict in sorting_scores.items():
                values = cell_dict.values()
                average = sum(values) / len(values)
                sorting_scores_same_type[cell_type] = average

            # Print the sorting scores for each cell type
            for cell_type, score in sorting_scores_same_type.items():
                for coll in [
                    coll 
                    for coll in bpy.data.collections 
                    if coll['type'] == cell_type
                ]:
                    coll['sorting score'] = score
                print(f"Cell Type: {cell}, Sorting Score: {score}")
            self.sorting_scores.update({scene.frame_current: sorting_scores_same_type})

    def set_scale(self, scale, cell_type):
        """
        Change the size of all cells of a certain type
        :param scale: value to set cell.modifiers["Cloth"].settings.shrink_min to.
        :param cell_type: Name of cell type to scale. (String)
        :return: None
        """

        num_cells = len(bpy.data.collections[cell_type].objects)
        for i in range(num_cells):
            cell_name = bpy.data.collections[cell_type].objects[i].name
            cell = bpy.data.objects[cell_name]
            cell.modifiers["Cloth"].settings.shrink_min = scale

    def contact_area_handler(self, scene, depsgraph): 

        if self.data_flag: 
            '''if (
                scene.frame_current 
                in range(self.frame_interval[0], self.frame_interval[1], 1)
            ):'''
            contact_ratio_dict, contact_areas_dict = get_contact_area()
            # Merge the dictionaries
            for key, value in contact_ratio_dict.items():
                if key not in self.contact_ratios:
                    # If the key is not present in the master dictionary, 
                    # create a new list with the current value
                    self.contact_ratios[key] = [value]
                else:
                    # If the key is already present in the master dictionary, 
                    # append the value to the existing list
                    self.contact_ratios[key].append(value)

            # Merge the dictionaries
            for key, value in contact_areas_dict.items():
                if key not in self.contact_areas:
                    # If the key is not present in the master dictionary, 
                    # create a new list with the current value
                    self.contact_areas[key] = [value]
                else:
                    # If the key is already present in the master dictionary, 
                    # append the value to the existing list
                    self.contact_areas[key].append(value)
            
                print(self.contact_areas)

    def adhesion_handler(self, scene, depsgraph):

        force_list = Force.force_list
        # print([force.name for force in force_list])

        for force in force_list:
            if not bpy.data.objects[force.name]['motion']: 
                assoc_cell = force.associated_cell
                bpy.context.view_layer.objects.active = bpy.data.objects[assoc_cell]
                # retrieve center of mass
                COM = get_centerofmass(bpy.data.objects[assoc_cell])
                # update the force location to its corresponding cell's center of mass
                bpy.data.objects[force.name].location = COM
                # cell type: 
                # cells will only adhere with cells from the same collection
                cloth_modifier = bpy.data.objects[assoc_cell].modifiers["Cloth"]
                cloth_settings = cloth_modifier.settings
                cloth_collection = bpy.data.objects[assoc_cell].users_collection[0]
                cloth_settings.effector_weights.collection = cloth_collection

                # bpy.data.objects[assoc_cell].modifiers["Cloth"].settings.effector_weights.
                # collection = bpy.data.objects[assoc_cell].users_collection[0]

    def set_random_motion_speed(self, motion_speed: float):

        self.random_motion_speed = motion_speed
        
    def data_export(self, scene, depsgraph): 

        self.seed = bpy.context.scene.get("seed")
        # var_to_export = ["sorting_scores", "times", "msd", "contact_ratios"]
        # Write the list at the end of the simulation
        if scene.frame_current == self.frame_interval[1]:
            with open(f"{self.data_file_path}_times.json", 'w') as write_file:
                write_file.write(json.dumps(self.times))
            with open(f"{self.data_file_path}_contact_ratios.json", 'w') as write_file:
                write_file.write(json.dumps(self.contact_ratios))
        # else: # data over simulation time
            with open(f"{self.data_file_path}_msd.json", 'w') as write_file:
                write_file.write(json.dumps(self.msd))
            with open(f"{self.data_file_path}_sorting_scores.json", 'w') as write_file:
                write_file.write(json.dumps(self.sorting_scores))
            with open(f"{self.data_file_path}_speed.json", 'w') as write_file:
                write_file.write(json.dumps(self.speed))
            with open(f"{self.data_file_path}_motion_path.json", 'w') as write_file:
                write_file.write(json.dumps(self.motion_path))
            with open(f"{self.data_file_path}_contact_areas.json", 'w') as write_file:
                write_file.write(json.dumps(self.contact_areas))
            with open(f"{self.data_file_path}_volumes.json", 'w') as write_file:
                write_file.write(json.dumps(self.volumes))

            if self.seed is not None: 
                with open(f"{self.data_file_path}_seed.json", 'w') as write_file:
                    write_file.write(json.dumps(self.seed))

            '''subprocess.run([
                "python",
                "C:\\Users\\anr9744\\Projects\\Goo\\scripts\\modules\\goo\\visualization.py",
                f"{self.data_file_path}"
            ])'''

    def visualize_stretching(self, scene, depsgraph):
        cells = [
            obj 
            for obj in bpy.data.objects 
            if "object" in obj.keys() and obj["object"] == "cell"
            ]

        for cell in cells:
            obj = bpy.data.objects.get(cell.name)

            # Make sure the object is a mesh
            if obj is not None and obj.type == 'MESH':
                # Duplicate the object and apply all modifiers
                mesh_eval = obj.evaluated_get(depsgraph)

                # Get the vertex coordinates of the original mesh
                vertices_original = [v.co for v in obj.data.vertices]

                # Calculate stretching ratio for each vertex
                stretching_ratios = []
                for vert_original, vert_evaluated in zip(
                    vertices_original, mesh_eval.data.vertices
                ):
                    stretching_ratio = (
                        (vert_evaluated.co - vert_original).length / 
                        vert_original.length
                    )
                    stretching_ratios.append(stretching_ratio)

                # Normalize the stretching ratios
                max_ratio = max(stretching_ratios)
                normalized_ratios = [
                    ratio / max_ratio if max_ratio != 0 else 0 
                    for ratio in stretching_ratios
                ]
                print(normalized_ratios)

                # Create a new vertex color layer if it doesn't exist
                color_layer = obj.data.vertex_colors.get("StretchingRatio")
                if color_layer is None:
                    color_layer = obj.data.vertex_colors.new(name="StretchingRatio")

                # Set vertex colors based on stretching ratio (gray to red gradient)
                for poly in obj.data.polygons:
                    for loop_index in poly.loop_indices:
                        vertex_index = obj.data.loops[loop_index].vertex_index
                        stretching_ratio = normalized_ratios[vertex_index]
                        color = (stretching_ratio, 1, 0, 1.0)  # Red-scale coloring
                        color_layer.data[loop_index].color = color

                # Create a new material and assign it to the object
                material = bpy.data.materials.new(name="StretchingMaterial")
                obj.data.materials.append(material)

                # Create a material slot for the object and link the material
                if obj.data.materials:
                    obj.data.materials[0] = material
                else:
                    obj.data.materials.append(material)

                # Set up vertex colors for the material
                material.use_nodes = True
                material.node_tree.nodes.clear()
                vertex_color_node = material.node_tree.nodes.new(
                    type='ShaderNodeVertexColor'
                )
                shader_node_output = material.node_tree.nodes.new(
                    type='ShaderNodeOutputMaterial'
                )
                material.node_tree.links.new(
                    vertex_color_node.outputs["Color"], 
                    shader_node_output.inputs["Surface"]
                )
                
                # Update the mesh to reflect the color changes
                obj.data.update()
            else:
                print("Object is not a mesh or not found.")

    # not used
    def data_handler(self, scene, depsgraph): 

        # initialization
        location_list = []
        com_list = []
        distances = set()
        total_dist = 0
        current_frame = bpy.data.scenes[0].frame_current

        print(f'Frame number: {current_frame}')
        print(self.data_file_path)

        ratio1, ratio2 = get_contact_area()
        print(f"Contact area: {ratio1}, {ratio2}")

        # get area with KD trees
        # areas, totals, ratios = mesh_contact_area_KD(0.03)
        # print(f"Contact area: {ratios[0]}, {ratios[1]}")

        # contact area 
        # contact_area1, surface_area1, ratio1 = compute_contact_area()
        # print(f"Contact ratios A1, A2: {contact_area1}, {surface_area1}, {ratio1}")
        # ratio = get_contact_area()
        '''self.contact_area.append(ratio1)'''

        # loop over each collection, then over each object
        for collection in bpy.data.collections:
            if 'Cells' in collection.name_full:
                # coll_name = collection.name_full

                # calculate the center of mass of cells in the scene
                for idx, cell in enumerate(collection.objects):
                    location_list.append((cell.name, cell.location))
                    bpy.context.view_layer.objects.active = cell
                    dg = bpy.context.evaluated_depsgraph_get()
                    cell_eval = cell.evaluated_get(dg)
                    vertices = cell_eval.data.vertices
                    vert_coords = np.asarray([
                        (cell_eval.matrix_world @ v.co) for v in vertices
                    ])
                    x = vert_coords[:, 0]
                    y = vert_coords[:, 1]
                    z = vert_coords[:, 2]
                    COM = (np.mean(x), np.mean(y), np.mean(z))
                    
                    # initialize the key for the cell name if the class dictionary 
                    com_list.append(COM)

                    # measure cell deformability between two frames
                    # idea: store coord of cell's COM, store coord of each vertex 
                    # of each cell's mesh at each frame
                    # at the end of the simulation, substract the COM's coord from each 
                    # vertex coord to isolate displacement caused by deformability
                    # at each frame, sum the 3D distance (displacement) between vertices 
                    # of current frame to previous frame 
                    # sum the overall displacement over the number of frame 

                for i in range(len(com_list)):
                    for j in range(len(com_list)):
                        # Calculate the Euclidean distance between the two coordinates
                        distance = math.sqrt(
                            sum(
                                [(a - b) ** 2 for a, b in zip(com_list[i], com_list[j])]
                            )
                        )
                        # Add the distance to the list of distances
                        distances.add(distance)

                total_dist = sum(distances)
                self.distances_tot.append(total_dist)
                self.frames.append(current_frame)

        # write the list at the end of the simulation
        if scene.frame_current == self.frame_interval[1]:

            self.data_dict['Frames'].append(self.frames)
            self.data_dict['Distances'].append(self.distances_tot)
            self.data_dict['Times'].append(self.times)

            # if file already exists, then merge new data to existing file
            # allows to save results over multiple runs in a single dict
            if os.path.isfile(f"{self.data_file_path}.json"): 
                with open(f"{self.data_file_path}.json", 'r') as f:
                    # data = f.read()
                    self.master_dict = json.load(f)

                    # append results of new run to the master dict
                    self.master_dict['Frames'].append(self.frames)
                    self.master_dict['Distances'].append(self.distances_tot)
                    self.master_dict['Times'].append(self.times)

                    # write master dict to file
                    with open(f"{self.data_file_path}.json", 'w') as write_file:
                        write_file.write(json.dumps(self.master_dict))

            # if file does not exist, create it with first dict
            else: 
                with open(f"{self.data_file_path}.json", 'a') as convert_file:
                    convert_file.write(json.dumps(self.data_dict))

            subprocess.run([
                "python",
                "C:\\Users\\anr9744\\Projects\\Goo\\scripts\\ \
                    modules\\goo\\visualization.py",
                f"{self.data_file_path}"
            ])        

    def timing_init_handler(self, scene, depsgraph): 

        if scene.frame_current == 2: 
            self.time = datetime.now()
            print(f'Render started for Frame 1 at: {self.time}')

    def timing_elapsed_handler(self, scene, depsgraph): 

        frame_written = scene.frame_current
        if frame_written not in [1, 2]: 
            elpased_time = datetime.now() - self.time
            elapsed_time_secs = elpased_time.seconds + elpased_time.microseconds/1000000
            self.times.update({frame_written: elapsed_time_secs})
            print('________________________________________________________')
            print(f"Render Started at:{self.time}")  
            print(f"Current frame: {frame_written}")
            print(f"Elapsed in seconds/microseconds:{elapsed_time_secs:.3f};"
                  f" {elapsed_time_secs:.1f}")

    def stop_animation(self, scene, depsgraph):

        # checks if the simulation has ended
        if scene.frame_current == self.frame_interval[1]:
            print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _')
            print("The simulation has ended.")
            # True enables the last frame not to be repeated
            bpy.ops.screen.animation_cancel(restore_frame=True) 
            # closes Blender then
            # bpy.ops.wm.quit_blender()
