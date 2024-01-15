from collections import defaultdict
import bpy
from mathutils import Vector, Matrix
import numpy as np
import sys
import math
from datetime import datetime
import json
import bmesh


"""
Authors: Antoine Ruzette, Sean Megason.
Developed and maintained by the Megason Lab, 2024. 

Sphynx docstring used for documentation. 
"""


def calculate_volume(obj):
    """Calculates the volume of the Blender mesh. 

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
    # print(f"Volume of {obj.name}: {abs(volume)}")
    # Free the bmesh
    bm.free()

    return volume


def get_centerofmass(obj): 
    """Calculates the center of mass of a mesh. 

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
    """
    Calculates the long axis of a 3D object in global coordinates.

    This function computes the long axis of a 3D object by determining the first
    eigenvector of its vertices' positions in global coordinates.

    :param bpy.data.objects['name'] obj: The Blender object for which to calculate
                                         the long axis.
    :returns: A tuple representing the coordinates of the long axis in global
              space as (x, y, z), and the length of the long axis.
    :rtype: tuple
    """

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

    # Convert the major axis to global coordinates
    major_axis_global = np.array(evaluated_object.matrix_world) @ np.append(major_axis, 
                                                                            1)

    # Find extreme points along the major axis
    min_point = np.min(vertex_coords.dot(major_axis))
    max_point = np.max(vertex_coords.dot(major_axis))
    length_long_axis = max_point - min_point

    return tuple(major_axis_global)[:3], length_long_axis


def get_long_axis(obj):
    """Calculates the long axis of a mesh. 

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
    """Retrieves the 2 angles associated with the long axis of a cell. 

    :param tuple axis: Coordinates of the long axis to the division plane 
                        as retrived by ``get_major_axis(obj)``.
    :returns: 
        - phi (:py:class:`tuple`) - Phi is the angle between the division axis \
            and the z-axis in spherical coordinates. 
        - theta (:py:class:`tuple`) - Theta is the angle projected on the xy plane \
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


def get_line_angles(cell, p1, p2):
    """
    Calculate the angles between a line and the global axes.

    This function calculates the angles between a line defined by two points
    (p1 and p2) and each of the global axes (x, y, z) of the given Blender cell.

    :param bpy.data.objects['name'] cell: The Blender object representing the cell.
    :param tuple p1: The coordinates of the first point of the line.
    :param tuple p2: The coordinates of the second point of the line.
    :returns: Angles between the line and each global axis in radians.
    :rtype: tuple
    """
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

 
def apply_modifiers(obj, which='all'): 
    """
    Apply modifiers to a Blender object.

    This function applies modifiers to a Blender object, either all modifiers or
    a specific one, based on the 'which' parameter.

    :param bpy.data.objects['name'] obj: The Blender object to which modifiers
                                         should be applied.
    :param str which: Specifies which modifier(s) to apply. Default is 'all'.
                     If 'all', all modifiers are applied; otherwise, specify the
                     name of the specific modifier to apply.
    """
    if which == 'all': 
        modifiers_to_apply = obj.modifiers

        for modifier in modifiers_to_apply:
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.modifier_apply(modifier=modifier.name)
    else:
        # Check if the specified modifier exists in the object's modifiers
        if which in obj.modifiers:
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.modifier_apply(modifier=which)
        else:
            print(f"Modifier '{which}' not found in object '{obj.name}'")

    return


def get_forces_in_collection(obj):
    """
    Get forces associated to the same cell.

    This function retrieves force objects within the same collection as the provided
    Blender object that have a tag 'object' with the value 'force'.

    :param bpy.data.objects['name'] obj: The Blender object whose collection
                                         is used for searching forces.
    :returns: A list of objects within the same collection tagged as 'force'.
    :rtype: list
    """
    forces = []

    for object in obj.users_collection[0].all_objects:
        if (
            object.get('object') == 'force'
        ): 
            forces.append(object)

    return forces


def get_cells_same_type(obj): 
    """
    Get cells of the same type.

    This function retrieves cells of the same type as the provided Blender
    object, and looks across different collections. Checks the `type` custom properties 
    of a Blender object. 

    :param bpy.data.objects['name'] obj: The Blender object representing a cell.
    :returns: A list of cell Blender object of the same type.
    :rtype: list
    """
    cells_same_type = []

    for coll in bpy.data.collections: 
        for cell in coll.all_objects: 
            if (
                cell.get('object') == 'cell' and 
                cell.users_collection[0].get('type') 
                == obj.users_collection[0].get('type') and 
                cell != obj
            ): 
                # Split the names into root and subsequent divisions
                # Naming by convention here
                cell_root = cell.name.split('.')[0:-1]
                obj_root = obj.name.split('.')[0:-1]

                # Compare roots to exclude cells with similar 
                # roots but different subsequent divisions
                if cell_root != obj_root:
                    cells_same_type.append(cell)

    return cells_same_type


def get_forces_same_type(obj): 
    """
    Get cells of the same type.

    This function retrieves forces of the same type as the provided Blender
    object, and looks across different collections. Checks the `type` custom properties 
    of a Blender object, but outputs the associated force. 

    :param bpy.data.objects['name'] obj: The Blender object representing a cell.
    :returns: A list of force Blender object of the same type.
    :rtype: list
    """
    forces_same_type = []

    for coll in bpy.data.collections: 
        for force in coll.all_objects: 
            if (
                force.get('object') == 'force' and 
                force.users_collection[0].get('type') 
                == obj.users_collection[0].get('type') and 
                force != obj
            ): 
                # Split the names into root and subsequent divisions
                # Naming by convention here
                cell_root = force.name.split('.')[0:-1]
                obj_root = obj.name.split('.')[0:-1]

                # Compare roots to exclude cells with similar 
                # roots but different subsequent divisions
                if cell_root != obj_root:
                    forces_same_type.append(force)

    return forces_same_type


def get_missing_adhesion_forces(obj1, obj2):
    """
    Get adhesion forces missing between two objects.

    This function identifies adhesion forces missing between two Blender objects.
    Adhesion forces are objects tagged with 'object' equal to 'force' and 'motion'
    equal to 0. It enforces daughter cells to inherit all adhesion forces 
    from their cell type.  

    :param bpy.data.objects['name'] obj1: The first Blender object.
    :param bpy.data.objects['name'] obj2: The second Blender object.
    :returns: A list of adhesion forces missing between the two objects.
    :rtype: list
    """
    forces_obj1 = {
        obj for obj in obj1.users_collection[0].all_objects
        if (
            obj.get('object') == 'force' and 
            obj.get('motion') == 0
        )
    }

    forces_obj2 = {
        obj for obj in obj2.users_collection[0].all_objects
        if (
            obj.get('object') == 'force' and 
            obj.get('motion') == 0
        )
    }

    missing_adhesion_forces = list(forces_obj1.symmetric_difference(forces_obj2))

    return missing_adhesion_forces


'''def gather_modifier_info(obj):
    modifiers_info = []
    for modifier in obj.modifiers:
        modifier_data = {'type': modifier.type, 'name': modifier.name}
        # Store modifier parameters
        for prop in modifier.bl_rna.properties:
            if prop.identifier != 'rna_type' and not prop.is_readonly:
                # Handling nested properties for collision_settings
                if prop.identifier == 'collision_settings': 
                    collision_settings = {}
                    for collision_prop in prop.type.bl_rna.properties:
                        if (collision_prop.identifier != 'rna_type' 
                                and not collision_prop.is_readonly):
                            collision_settings[collision_prop.identifier] = \
                                getattr(getattr(modifier, prop.identifier),
                                        collision_prop.identifier)
                    modifier_data[prop.identifier] = collision_settings
                # Handling cloth modifier settings
                elif prop.identifier == 'settings' and modifier.type == 'CLOTH':
                    cloth_settings = {}
                    for cloth_prop in prop.type.bl_rna.properties:
                        if (cloth_prop.identifier != 'rna_type' 
                                and not cloth_prop.is_readonly):
                            cloth_settings[cloth_prop.identifier] = \
                                getattr(getattr(modifier, prop.identifier), 
                                        cloth_prop.identifier)
                    modifier_data[prop.identifier] = cloth_settings
                else:
                    modifier_data[prop.identifier] = getattr(modifier, prop.identifier)
        modifiers_info.append(modifier_data)
    return modifiers_info'''


def store_subsurf_settings(obj):
    """
    Store Subdivision Surface (Subsurf) modifier settings.

    This function stores the settings of Subsurf modifiers applied to a Blender
    object, including subdivision type, levels, and quality.

    :param bpy.data.objects['name'] obj: The Blender cell object containing 
                                         Subsurf modifiers.
    :returns: A dictionary storing Subsurf modifier settings for each modifier.
    :rtype: dict
    """
    subsurf_settings = {}
    for modifier in obj.modifiers:
        if modifier.type == 'SUBSURF':
            subsurf_settings[modifier.name] = {}
            if 'subdivision_type' in dir(modifier):
                subsurf_settings[modifier.name]['subdivision_type'] = \
                    modifier.subdivision_type
            if 'levels' in dir(modifier):
                subsurf_settings[modifier.name]['levels'] = modifier.levels
            if 'quality' in dir(modifier):
                subsurf_settings[modifier.name]['quality'] = modifier.quality
    return subsurf_settings


def declare_subsurf_settings(obj, subsurf_settings):
    """
    Declare Subdivision Surface (Subsurf) modifiers with specified settings.

    This function adds Subsurf modifiers to a Blender object based on the provided
    settings dictionary.

    :param bpy.data.objects['name'] obj: The Blender cell object to which Subsurf
                                         modifiers will be added.
    :param dict subsurf_settings: A dictionary storing Subsurf modifier settings
                                  for each modifier.
    :return: None
    """
    for modifier_name, settings in subsurf_settings.items():
        new_modifier = obj.modifiers.new(name=modifier_name, type='SUBSURF')
        if 'subdivision_type' in settings:
            new_modifier.subdivision_type = settings['subdivision_type']
        if 'levels' in settings:
            new_modifier.levels = settings['levels']
        if 'quality' in settings:
            new_modifier.quality = settings['quality']

    return


def store_cloth_settings(obj):
    """
    Store settings of Cloth modifiers applied to a Blender object.

    This function stores the settings of Cloth modifiers applied to a Blender cell
    object, including various parameters related to cloth simulation.

    :param bpy.data.objects['name'] obj: The Blender cell object 
                                         containing Cloth modifiers.
    :returns: A dictionary storing Cloth modifier settings for each modifier.
    :rtype: dict
    """
    cloth_settings = {}
    for modifier in obj.modifiers:
        if modifier.type == 'CLOTH':
            cloth_settings[modifier.name] = {}
            settings = modifier.settings
            cloth_cell_settings = cloth_settings[modifier.name]
            
            # Using cloth_cell_settings for clarity
            cloth_cell_settings['quality'] = settings.quality
            cloth_cell_settings['air_damping'] = settings.air_damping
            cloth_cell_settings['bending_model'] = settings.bending_model
            cloth_cell_settings['mass'] = settings.mass
            cloth_cell_settings['time_scale'] = settings.time_scale
            cloth_cell_settings['tension_stiffness'] = settings.tension_stiffness
            cloth_cell_settings['compression_stiffness'] = \
                settings.compression_stiffness
            cloth_cell_settings['shear_stiffness'] = settings.shear_stiffness
            cloth_cell_settings['bending_stiffness'] = settings.bending_stiffness
            cloth_cell_settings['tension_damping'] = settings.tension_damping
            cloth_cell_settings['compression_damping'] = settings.compression_damping
            cloth_cell_settings['shear_damping'] = settings.shear_damping
            cloth_cell_settings['bending_damping'] = settings.bending_damping
            cloth_cell_settings['use_internal_springs'] = False
            cloth_cell_settings['internal_spring_max_length'] = \
                settings.internal_spring_max_length
            cloth_cell_settings['internal_spring_max_diversion'] = \
                settings.internal_spring_max_diversion
            cloth_cell_settings['internal_spring_normal_check'] = \
                settings.internal_spring_normal_check
            cloth_cell_settings['internal_tension_stiffness'] = \
                settings.internal_tension_stiffness
            cloth_cell_settings['internal_compression_stiffness'] = \
                settings.internal_compression_stiffness
            cloth_cell_settings['internal_tension_stiffness_max'] = \
                settings.internal_tension_stiffness_max
            cloth_cell_settings['internal_compression_stiffness_max'] = \
                settings.internal_compression_stiffness_max
            cloth_cell_settings['use_pressure'] = settings.use_pressure
            cloth_cell_settings['uniform_pressure_force'] = \
                settings.uniform_pressure_force
            cloth_cell_settings['use_pressure_volume'] = settings.use_pressure_volume
            cloth_cell_settings['target_volume'] = settings.target_volume
            cloth_cell_settings['pressure_factor'] = settings.pressure_factor
            cloth_cell_settings['fluid_density'] = settings.fluid_density

            cloth_cell_settings['collision_settings'] = {}
            collision_settings = modifier.collision_settings
            cloth_cell_settings['collision_settings']['collision_quality'] = \
                collision_settings.collision_quality
            cloth_cell_settings['collision_settings']['use_collision'] = \
                collision_settings.use_collision
            cloth_cell_settings['collision_settings']['use_self_collision'] = \
                collision_settings.use_self_collision
            cloth_cell_settings['collision_settings']['self_friction'] = \
                collision_settings.self_friction
            cloth_cell_settings['collision_settings']['friction'] = \
                collision_settings.friction
            cloth_cell_settings['collision_settings']['self_distance_min'] = \
                collision_settings.self_distance_min
            cloth_cell_settings['collision_settings']['distance_min'] = \
                collision_settings.distance_min
            cloth_cell_settings['collision_settings']['self_impulse_clamp'] = \
                collision_settings.self_impulse_clamp
    
    return cloth_settings


def declare_cloth_settings(obj, cloth_settings):
    """
    Declare Cloth modifiers with specified settings on a cell object.

    This function adds Cloth modifiers to a Blender cell object based on the provided
    settings dictionary.

    :param bpy.data.objects['name'] obj: The Blender cell object to which Cloth
                                         modifiers will be added.
    :param dict cloth_settings: A dictionary storing Cloth modifier settings
                               for each modifier.
    :return: None
    """
    for modifier_name, settings in cloth_settings.items():
        new_modifier = obj.modifiers.new(name=modifier_name, type='CLOTH')
        new_settings = new_modifier.settings
        for setting_name, value in settings.items():
            if setting_name != 'collision_settings':
                setattr(new_settings, setting_name, value)
            else:
                collision_settings = settings['collision_settings']
                new_collision_settings = new_modifier.collision_settings
                for coll_sett_name, coll_sett_value in collision_settings.items():
                    setattr(new_collision_settings, 
                            coll_sett_name, 
                            coll_sett_value)
                    
    return 


def store_collision_settings(obj):
    """
    Store settings of Collision modifiers applied to a cell object.

    This function stores the settings of Collision modifiers applied to a Blender cell
    object, including various collision-related parameters.

    :param bpy.data.objects['name'] obj: The Blender cell object 
                                         containing Collision modifiers.
    :returns: A dictionary storing Collision modifier settings for each modifier.
    :rtype: dict
    """
    collision_settings = {}
    for modifier in obj.modifiers:
        if modifier.type == 'COLLISION':
            collision_settings[modifier.name] = {}
            settings = modifier.settings
            collision_cell_settings = collision_settings[modifier.name]
            collision_cell_settings['use_culling'] = settings.use_culling
            collision_cell_settings['damping'] = settings.damping
            collision_cell_settings['thickness_outer'] = settings.thickness_outer
            collision_cell_settings['thickness_inner'] = settings.thickness_inner
            collision_cell_settings['cloth_friction'] = settings.cloth_friction
            # ... (other Collision settings)
    return collision_settings


def declare_collision_settings(obj, collision_settings):
    """
    Declare Collision modifiers with specified settings on a cell object.

    This function adds Collision modifiers to a Blender cell object 
    based on the provided settings dictionary.

    :param bpy.data.objects['name'] obj: The Blender cell object to which Collision
                                         modifiers will be added.
    :param dict collision_settings: A dictionary storing Collision modifier settings
                                    for each modifier.
    :return: None
    """
    for modifier_name, settings in collision_settings.items():
        new_modifier = obj.modifiers.new(name=modifier_name, type='COLLISION')
        new_settings = new_modifier.settings
        for setting_name, value in settings.items():
            setattr(new_settings, setting_name, value)
        
    return


def store_remesh_settings(obj):
    """
    Store settings of Remesh modifiers applied to a cell object.

    This function stores the settings of Remesh modifiers applied to a Blender cell
    object, including various remesh-related parameters.

    :param bpy.data.objects['name'] obj: The Blender cell object 
                                         containing Remesh modifiers.
    :returns: A dictionary storing Remesh modifier settings for each modifier.
    :rtype: dict
    """
    remesh_settings = {}
    for modifier in obj.modifiers:
        if modifier.type == 'REMESH':
            remesh_settings[modifier.name] = {}
            remesh_cell_settings = remesh_settings[modifier.name]
            remesh_cell_settings['mode'] = modifier.mode
            remesh_cell_settings['voxel_size'] = modifier.voxel_size
            remesh_cell_settings['adaptivity'] = modifier.adaptivity
            remesh_cell_settings['use_remove_disconnected'] = \
                modifier.use_remove_disconnected
            remesh_cell_settings['use_smooth_shade'] = modifier.use_smooth_shade
            remesh_cell_settings['show_in_editmode'] = modifier.show_in_editmode
    return remesh_settings


def declare_remesh_settings(obj, remesh_settings):
    for modifier_name, settings in remesh_settings.items():
        new_modifier = obj.modifiers.new(name=f"Remesh_{obj.name}", type='REMESH')
        for setting_name, value in settings.items():
            setattr(new_modifier, setting_name, value)


def divide_boolean(obj): 
    """
    Divide a cell object along its division plane using boolean operations.

    This function performs boolean division on a Blender cell object, creating two
    daughter cells. It also handles the declaration of daughter forces, 
    inheritance of properties, and separation of the mesh.

    :param bpy.data.objects['name'] obj: The Blender object to be divided.

    :return: 
             - bpy.data.objects['name']: First daughter cell.
             - bpy.data.objects['name']: Second daughter cell.
             - float: Adhesion strength.
             - float: Motion strength.
             - float: Falloff.
             - bool: Flag indicating if there was a mother adhesion force.
             - bool: Flag indicating if there was a mother motion force.
             - list: List of dictionaries containing values of the modifiers applied \
                    to the original object before division.   

    """
    # Initialize variables
    adhesion_strength = 0  # Default value
    motion_strength = 0  # Default value
    falloff = 0  # Default value
    mother_adhesion = False  # Default value
    mother_motion = False  # Default value

    # Get COM
    p1 = get_centerofmass(obj)
    # Get the long axis in global coordinate from the origin
    p2, _ = get_long_axis_global(obj)    

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

    # Get the angles representing the orientation of the division plane
    x, y, z = get_line_angles(obj, p1, p2)
    # Get name of dividing cell
    mother_name = obj.name
    # Create division plane
    plane = bpy.ops.mesh.primitive_plane_add(enter_editmode=False, 
                                             size=5, 
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
    solid_mod.thickness = 0.03
    bpy.ops.object.modifier_apply(modifier=solid_mod.name)

    Force.force_list = [
        force for force in Force.force_list if force.name != obj['adhesion force']
    ]
    
    # Remove mother adhesion force
    if 'adhesion force' in obj: 
        mother_adhesion = True
        adhesion_force = bpy.data.objects[obj['adhesion force']]
        adhesion_strength = adhesion_force.field.strength
        falloff = adhesion_force.field.falloff_power
        if obj['adhesion force'] in bpy.data.objects:
            bpy.data.objects.remove(
                bpy.data.objects[obj['adhesion force']], 
                do_unlink=True
            )

    # Remove mother motion force
    if 'motion force' in obj: 
        mother_motion = True
        motion_force = bpy.data.objects[obj['motion force']]
        motion_strength = motion_force.field.strength

        if obj['motion force'] in bpy.data.objects:
            bpy.data.objects.remove(
                bpy.data.objects[obj['motion force']], 
                do_unlink=True
            )

    # Add boolean modifier to the original object
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.modifier_add(type='BOOLEAN')
    bool_mod = obj.modifiers[-1]
    bool_mod.operand_type = 'OBJECT'
    bool_mod.object = bpy.data.objects[plane_name]
    bool_mod.operation = 'DIFFERENCE'
    bool_mod.solver = 'EXACT'

    # bpy.ops.object.modifier_apply(modifier="Cloth")
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

    d1 = bpy.context.selected_objects[0]  # to be renamed
    d1.name = f"{mother_name}.01"  # new cell

    d2 = bpy.context.scene.objects[mother_name]
    d2.name = f"{mother_name}.02"  # mother cell

    # Declare collections to contain daughter cells
    mother_collection = obj.users_collection[0]
    mother_type = mother_collection.get('type')
    d1_collection = make_collection(f'{d1.name}_collection', type=mother_type)
    d2_collection = make_collection(f'{d2.name}_collection', type=mother_type)

    bpy.ops.object.select_all(action='DESELECT')

    bpy.context.view_layer.objects.active = d1
    # remove duplicate objects outside of the collection
    bpy.ops.collection.objects_remove_all()
    # Add the active cell to our specific collection 
    bpy.data.collections[d1_collection.name].objects.link(d1)
 
    # Pass on the mother force
    d1.select_set(True)
    bpy.ops.object.convert(target='MESH')
    bpy.context.view_layer.update()
    # Adding a Remesh modifier to the converted mesh
    remesh_modifier = d1.modifiers.new(name='Remesh', type='REMESH')
    remesh_modifier.mode = 'VOXEL'
    remesh_modifier.voxel_size = 0.2  # microns
    remesh_modifier.adaptivity = 0
    remesh_modifier.use_remove_disconnected = True
    remesh_modifier.use_smooth_shade = True
    bpy.context.view_layer.update()
    bpy.ops.object.modifier_apply(modifier="Remesh")
    bpy.context.view_layer.update()

    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.objects.active = d2
    # remove duplicate objects outside of the collection
    bpy.ops.collection.objects_remove_all()
    # Add the active cell to our specific collection 
    bpy.data.collections[d2_collection.name].objects.link(d2)

    d2.select_set(True)
    bpy.ops.object.convert(target='MESH')
    bpy.context.view_layer.update()
    remesh_modifier = d2.modifiers.new(name='Remesh', type='REMESH')
    remesh_modifier.mode = 'VOXEL'
    remesh_modifier.voxel_size = 0.2  # microns
    remesh_modifier.adaptivity = 0
    remesh_modifier.use_remove_disconnected = True
    remesh_modifier.use_smooth_shade = True
    bpy.context.view_layer.update()
    bpy.ops.object.modifier_apply(modifier="Remesh")
    bpy.context.view_layer.update()

    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.update()

    # Retrieve the collection by name

    if mother_collection:
        # Remove the collection
        bpy.data.collections.remove(mother_collection)
    else:
        raise ValueError(f"Collection '{mother_collection}' not found")

    # Delete the plane if it exists
    if plane_name in bpy.data.objects:
        # Get the object reference
        obj = bpy.data.objects[plane_name]
        # Remove the object
        bpy.data.objects.remove(obj, do_unlink=True)
    else:
        print(f"Division plane: '{plane_name}' not found.")

    return d1, d2, adhesion_strength, motion_strength, \
        falloff, mother_adhesion, mother_motion, modifier_values


def disable_internal_springs(objects):
    """
    Disable internal springs for Cloth modifiers in a list of Blender cell objects.

    This function iterates through the provided list of Blender objects and, for each
    daughter cell object, checks if a Cloth modifier exists 
    and disables internal springs.

    :param list objects: A list of Blender cell objects.

    :return: None
    """
    # Disable internal springs for each object in the list
    for obj in objects:
        if obj is None:
            continue
        
        # Disable internal springs for daughter cells
        for daughter_cell in objects:
            if daughter_cell:
                modifiers = daughter_cell.modifiers
                
                # Check if a Cloth modifier exists and disable internal springs
                cloth_modifier = modifiers.get("Cloth")
                if cloth_modifier:
                    cloth_modifier.settings.use_internal_springs = False

    return


def apply_physics(obj, stiffness): 
    """
    Apply cloth physics and related modifiers to a Blender cell object.

    This function sets up cloth physics and associated modifiers for the 
    provided Blender object. It includes Subdivision Surface (Subsurf), 
    Cloth, Collision, and Remesh modifiers.

    :param bpy.data.objects['name'] obj: The Blender object to apply physics to.
    :param float stiffness: Stiffness value for Cloth modifier.

    :return: None
    """
    print('Starting modifiers declaration')
    
    # Select object and make it active
    bpy.ops.object.select_all(action='DESELECT')    
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)

    # Add subsurface modifier to make smoother
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].subdivision_type = 'CATMULL_CLARK'
    bpy.context.object.modifiers["Subdivision"].levels = 1
    bpy.context.object.modifiers["Subdivision"].quality = 3

    # Add cloth settings 
    bpy.ops.object.modifier_add(type='CLOTH')
    bpy.context.object.modifiers['Cloth'].settings.quality = 4
    bpy.context.object.modifiers['Cloth'].settings.air_damping = 10
    bpy.context.object.modifiers['Cloth'].settings.bending_model = 'ANGULAR'
    bpy.context.object.modifiers["Cloth"].settings.mass = 1
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1

    # Cloth > Stiffness 
    bpy.context.object.modifiers['Cloth'].settings.tension_stiffness = stiffness
    bpy.context.object.modifiers['Cloth'].settings.compression_stiffness = stiffness
    bpy.context.object.modifiers['Cloth'].settings.shear_stiffness = stiffness
    bpy.context.object.modifiers['Cloth'].settings.bending_stiffness = 1
    # Cloth > Damping
    bpy.context.object.modifiers['Cloth'].settings.tension_damping = 50
    bpy.context.object.modifiers['Cloth'].settings.compression_damping = 50
    bpy.context.object.modifiers['Cloth'].settings.shear_damping = 50
    bpy.context.object.modifiers['Cloth'].settings.bending_damping = 0.5
    # Cloth > Internal Springs
    bpy.context.object.modifiers['Cloth'].settings.use_internal_springs = False
    bpy.context.object.modifiers['Cloth'].settings.internal_spring_max_length = 1
    bpy.context.object.modifiers['Cloth'].settings.internal_spring_max_diversion = \
        0.785398
    bpy.context.object.modifiers['Cloth'].settings.internal_spring_normal_check = False
    bpy.context.object.modifiers['Cloth'].settings.internal_tension_stiffness = 10
    bpy.context.object.modifiers['Cloth'].settings.internal_compression_stiffness = 10
    bpy.context.object.modifiers['Cloth'].settings.internal_tension_stiffness_max = \
        10000
    bpy.context.object.modifiers['Cloth'].settings.internal_compression_stiffness_max \
        = 10000
    # Cloth > Pressure
    bpy.context.object.modifiers["Cloth"].settings.use_pressure = True
    bpy.context.object.modifiers['Cloth'].settings.uniform_pressure_force = 5
    bpy.context.object['previous pressure'] = \
        bpy.context.object.modifiers['Cloth'].settings.uniform_pressure_force
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
    bpy.context.object.modifiers['Cloth'].collision_settings.self_distance_min = 0.002
    bpy.context.object.modifiers['Cloth'].collision_settings.distance_min = 0.002
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
    remesh_mod.mode = 'VOXEL'
    remesh_mod.voxel_size = 0.2  # microns
    remesh_mod.adaptivity = 0
    # remesh_mod.mode = 'SMOOTH'
    # remesh_mod.octree_depth = 4
    # remesh_mod.scale = 0.75
    remesh_mod.use_remove_disconnected = True
    remesh_mod.use_smooth_shade = False
    remesh_mod.show_in_editmode = True
    # bpy.ops.object.modifier_move_to_index(modifier=f"Remesh1_{obj.name}", index=3)

    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.update()

    # Enable the cache of the daughter cloth simulations to match 
    # the simulation start and end frames
    obj.modifiers["Cloth"].point_cache.frame_start = bpy.context.scene.frame_start
    obj.modifiers["Cloth"].point_cache.frame_end = bpy.context.scene.frame_end


def cast_to_sphere(obj, factor=0.5):
    """
    Apply a Cast modifier to round up a cell object.

    This function applies a Cast modifier to the provided Blender cell object, 
    creating a rounded shape.

    :param bpy.data.objects['name'] obj: The Blender object to be rounded.
    :param float factor: The factor to adjust the level of rounding. Default is 0.5.

    :return: None
    """
    # Select the object to round up
    bpy.context.view_layer.objects.active = obj

    if obj is None:
        print(f"Object named '{obj}' not found.")
        return
    
    # Apply a "Cast" modifier to the duplicated object
    cast_modifier = obj.modifiers.new(name="Cast", type='CAST')
    cast_modifier.factor = factor  # Adjust the factor for the level of rounding
    
    # Apply the modifier to generate the rounded shape
    # bpy.ops.object.modifier_apply(modifier="Cast")


def to_sphere(obj): 
    """
    Forces a Blender cell object into a sphere shape.

    This function updates the object's origin, simplifies its topology, and returns
    the modified object.

    :param bpy.data.objects['name'] obj: The Blender object to be transformed \
        into a sphere.

    :return: None
    """
    obj = bpy.context.active_object

    # Update the object and the shape key
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS', center='MEDIAN')
    bpy.context.view_layer.update()

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

    return


def make_collection(name, type): 
    """
    Create and link a new Blender collection to the scene.

    This function creates a new collection with the given cell name and cell type,
    and links it to the current scene. Collections act as placeholder for cell types.

    :param str name: The name of the new collection.
    :param str type: The type of the new collection.

    :return: The newly created collection.
    :rtype: bpy.types.Collection
    """
    # Create a new collection
    collection = bpy.data.collections.new(name)
    # Cell type
    collection['type'] = type
    # Blender object type
    collection['object'] = 'collection'
    # link the collection to the scene for visualization 
    bpy.context.scene.collection.children.link(collection)

    return collection


def make_cell(
    name,
    loc,
    type,
    radius=1,
    remeshing=True,
    scale=(1, 1, 1),
    stiffness=1,
    material=("bubble", 0.007, 0.021, 0.3), 
    arcdiv=5,
    subdiv=3, 
    mesh='icosphere'
):

    """
    Creates a Blender cell object.

    This function creates a new Blender cell object with various settings,
    such as mesh type, scale, material, and modifiers.

    :param str name: The name of the new cell object.
    :param tuple loc: The location of the cell object.
    :param str type: The type of the new cell object.
    :param float radius: The radius of the cell object.
    :param bool remeshing: Whether to apply remeshing modifiers.
    :param tuple scale: The scale of the cell object.
    :param float stiffness: The stiffness of the cloth modifier.
    :param tuple material: The material properties for the cell.
    :param int arcdiv: The arc divisions for the primitive_round_cube_add function.
    :param int subdiv: The subdivision levels for the Subsurf modifier.
    :param str mesh: The type of mesh to create ('roundcube' or 'icosphere').

    :returns: None
    """

    collection = make_collection(f'{name}_collection', type=type)

    if mesh == 'roundcube': 
        # Create mesh
        bpy.ops.mesh.primitive_round_cube_add(change=False,
                                              radius=radius,
                                              size=scale,
                                              arc_div=arcdiv,
                                              lin_div=0,
                                              div_type='CORNERS',
                                              odd_axis_align=False,
                                              no_limit=False,
                                              location=loc
                                              )
        
    elif mesh == 'icosphere': 
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=2,
                                              radius=radius,
                                              calc_uvs=True,
                                              align='WORLD',
                                              location=loc,
                                              scale=scale
                                              )
    else:
        raise ValueError("Invalid mesh type. Goo currently supports \
                         'roundcube' and 'icosphere'")

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

    # Smooth the mesh
    bpy.ops.object.select = True
    bpy.ops.object.shade_smooth()

    # Add subsurface modifier to make smoother
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = subdiv

    # Add cloth settings 
    bpy.ops.object.modifier_add(type='CLOTH')
    bpy.context.object.modifiers['Cloth'].settings.quality = 4
    bpy.context.object.modifiers['Cloth'].settings.air_damping = 0
    bpy.context.object.modifiers['Cloth'].settings.bending_model = 'ANGULAR'
    bpy.context.object.modifiers["Cloth"].settings.mass = 1
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
    bpy.context.object['previous pressure'] = \
        bpy.context.object.modifiers['Cloth'].settings.uniform_pressure_force
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

    if remeshing: 

        remesh_mod = obj.modifiers.new(name=f"Remesh_{obj.name}", type='REMESH')
        obj.modifiers[f"Remesh_{obj.name}"].name = f"Remesh_{obj.name}"
        remesh_mod.mode = 'VOXEL'
        remesh_mod.voxel_size = 0.2  # microns
        remesh_mod.adaptivity = 0
        remesh_mod.use_remove_disconnected = True
        remesh_mod.use_smooth_shade = True
        remesh_mod.show_in_editmode = True
        bpy.ops.object.modifier_move_to_index(modifier=f"Remesh1_{obj.name}", index=3)

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

    return


def add_material(mat_name, r, g, b):
    """Creates a soap bubble-like Blender material for use in rendering 
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
    """Creates Force (Python) objects. 

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
        """
        Get the strength attribute of Force.

        :return: The strength of a Blender force.
        :rtype: float
        """
        return self.strength
    
    def get_falloff(self): 
        """
        Get the falloff attribute of Force.

        :return: The falloff of a Blender force.
        :rtype: float
        """
        return self.falloff_power

    def get_blender_force(self):
        """
        Get the Blender object of Force.

        :return: The Blender force obj.
        :rtype: float
        """
        obj = bpy.data.objects[self.name]
        return obj


def add_boundaries(loc, size, shape='BOX', type='REFLECTIVE', name='box'):
    """
    Adds reflective boundaries to the scene.

    The boundaries can be of the shape of sphere or boxes. 

    :param tuple loc: A tuple containing X, Y, and Z coordinates of the center of \
        the boundary.
    :param tuple size: A tuple containing dimensions in X, Y, and Z \
        for the boundary.
    :param str shape: The shape of the boundary. Supported values are \
        'BOX' and 'SPHERE'.
    :param str type: The type of boundary. (Unused parameter in the provided code)
    :param str name: The name to be assigned to the boundary object.

    :return: None
    """
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
               size=0.5):
    """
    Adds motion force to a cell.

    This function creates a motion force effector based on the specified parameters and 
    links it to the corresponding cell.

    :param str effector_name: The name of the cell object to which the motion \
        force will be applied.
    :param float strength: The strength of the motion force.
    :param float persistence: The persistence of the motion force (default is 0).
    :param float randomness: The randomness of the motion force.
    :param str distribution: The distribution type of the motion force \
        (default is 'uniform').
    :param float size: The size of the distribution for the motion force.

    :return: None
    """
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

    else: 
        print('---- Types not supported')
    
    return


def add_hetero_adhesion(cell_name, other_type_name, strength):
    """
    Adds heterotypic adhesion force between a cell and cells of another type.

    This function creates a heterotypic adhesion force between a cell and cells 
    of another specified type. The adhesion force is applied to the cell based 
    on the specified strength.

    :param str cell_name: The name of the cell to which the heterotypic adhesion force \
          will be applied.
    :param str other_type_name: The name of the other cell type to which the adhesion \
        force will be applied.
    :param float strength: The strength of the heterotypic adhesion force.

    Returns: None
    """
    hetero_collections = [
        coll for coll in bpy.data.collections if coll.get('type') in other_type_name
    ]

    force = make_force(force_name=f'heterotypic_{cell_name}_{other_type_name}',
                       cell_name=cell_name,
                       type=type,
                       strength=strength,
                       falloff=0.5,
                       motion=False)
    
    for coll in hetero_collections:
        coll.objects.link(bpy.data.objects[force.name])

    return


def add_homo_adhesion(cell_name, strength):
    """
    Adds homotypic adhesion force within a cell type.

    This function creates a homotypic adhesion force within a cell type. 
    The adhesion force is applied to the specified cell based on the specified 
    strength, and it affects other cells of the same type.

    :param str cell_name: The name of the cell to which the homotypic adhesion force \
        will be applied.
    :param float strength: The strength of the homotypic adhesion force.

    :return: None
    """
    homo_collections = [
        coll for coll in bpy.data.collections
        if coll.get('type') 
        in bpy.data.objects[cell_name].users_collection[0].get("type")
    ]
    print(f"Homotypic collections: {homo_collections}")

    cell_type = f"{bpy.data.objects[cell_name].users_collection[0]['type']}"
    force = make_force(force_name=f'homotypic_{cell_name}',
                       cell_name=cell_name,
                       type=cell_type,
                       strength=strength,
                       falloff=0.5,
                       motion=False)
    
    for coll in homo_collections:
        coll.objects.link(bpy.data.objects[force.name])

    return


def make_force(force_name,
               cell_name,
               type,
               strength,
               falloff=0.5,
               motion=False,
               min_dist=0,
               max_dist=1.2):
    """
    Creates a Blender force from a Goo :class:`Force` object.

    This function creates a Blender force object based on the provided parameters. 
    It adds the force to the associated Blender collection and sets various force 
    parameters such as strength, falloff, and motion.

    :param str force_name: The name of the force object.
    :param str cell_name: The name of the associated cell object.
    :param str type: The type of the cell.
    :param float strength: The strength of the force.
    :param float falloff: The falloff power of the force (default is 0.5).
    :param bool motion: A flag indicating whether the force is a motion force \
        (default is False).
    :param float min_dist: The minimum distance for the force (default is 0).
    :param float max_dist: The maximum distance for the force (default is 1.2).

    :return: The created Goo :class:`Force` object.
    :rtype: Force
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
        bpy.context.object.name = force_name
        bpy.data.objects.get(cell_name)["adhesion force"] = force_name

    elif motion:

        rand_coord = tuple(np.random.uniform(low=-0.05, high=0.05, size=(3,)))
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

    # Add force parameters
    bpy.context.object.field.strength = force.strength
    bpy.context.object.field.use_max_distance = True
    bpy.context.object.field.use_min_distance = True
    bpy.context.object.field.distance_max = max_dist
    bpy.context.object.field.distance_min = min_dist
    bpy.context.object.name = force.name
    bpy.context.object.field.falloff_power = force.falloff_power
    bpy.context.object['cell'] = cell_name
    bpy.context.object['object'] = 'force'

    scene_collection = bpy.context.scene.collection
    scene_collection.objects.unlink(bpy.data.objects[force.name])

    return force


def setup_world(seed=None):
    """Sets up the default values used for simulations in Goo 
    including units and rendering background. 

    :param seed: Seed value for random number generation (optional)
    :type seed: int or None
    :returns: None
    """

    # Add required add-ons
    addon_name = "add_mesh_extra_objects"
    # Check if the addon is not already enabled
    if addon_name not in bpy.context.preferences.addons:
        bpy.ops.preferences.addon_enable(module=addon_name)
        print(f"Addon '{addon_name}' has been enabled.")
    else:
        print(f"Addon '{addon_name}' is already enabled.")

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
    """Sets up Blender World properties for use in rendering.

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
    """Renders a simulation to create a set of still images that 
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


def calculate_com(cell): 
    """Calculates the center of mass of a Blender cell object. 

    :param bpy.data.objects[cell.name] cell: cell
    :return: The coordinate of the center of mass of the cell. 
    :rtype: tuple(x,y,z)
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


def get_contact_area():
    """
    Calculate the contact ratio between cells in the scene.

    This function computes the contact ratio between cell objects in the Blender scene. 
    Loops over all cells.
    It considers the contact area between each pair of cells and calculates the ratio 
    of the contact area to the total surface area of each cell.

    :param float threshold: The distance threshold to consider vertices in contact \
        (default is 0.05).

    :return: A dictionary containing cell pairs as keys and their corresponding \
        contact ratios and contact areas as values. The contact ratio represents \
            the ratio of the contact area to the total surface area of each cell. \
                The contact area represents the total area of contact between the cells.
    :rtype: dict, dict
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


def constrict(obj, indices): 
    """
    Constricts selected vertices along a line defined \
        by the center of mass and the long axis of the object.

    This function constrains the selected vertices of the given object along a line
    defined by the center of mass (COM) and the long axis of the object. 
    The selected vertices are resized to create a constriction effect.

    :param obj: The Blender cell object to constrict.
    :type obj: bpy.types.Object
    :param indices: A list of vertex indices to be constricted.
    :type indices: list[int]

    :return: None
    """
    # Define the two tuple coordinates that define the line
    p1 = get_centerofmass(obj)
    p2 = get_long_axis_global(obj)
    print(f'COM: {p1}, long axis: {p2}')

    # Set the active object and get the evaluated mesh
    bpy.context.view_layer.objects.active = obj
    
    # Get the edit-mode mesh and create a BMesh
    mesh = obj.data

    # Sort the indices to ensure that the vertices are selected in order
    for index in indices:
        mesh.vertices[index].select = True

    bpy.ops.object.mode_set(mode='EDIT') 
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
    bpy.context.view_layer.update()


def get_division_angles(cell, alpha): 

    # Get the active object
    cell = bpy.context.active_object

    # Define the two tuple coordinates that define the line
    p1 = get_centerofmass(cell)
    p2 = get_long_axis_global(cell)

    # Get the vector representing the line
    line_vector = Vector((p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])).normalized()

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
    
    print(f"Intersect coord: {len(intersection_verts)} vertices, {intersection_verts}")

    # Create new vertices from intersection coordinates
    mesh = bpy.data.meshes.new('Intersection')
    obj = bpy.data.objects.new('Intersection', mesh)
    bpy.context.scene.collection.objects.link(obj)

    # Set the mesh data
    mesh.from_pydata(intersection_verts, [], [])
    mesh.update()

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
    """
    A class for creating different types of handlers that trigger
    actions on ``Goo`` cells when certain criteria are met.

    The initialization function specifies available cell types and associated
    parameters like division rate, growth rate, and adhesion forces. It also
    initializes various data containers for simulation metrics.
    
    """

    # The initialization function specifies available cell types and associated
    # parameters like division rate, growth rate, and adhesion forces
    def __init__(self):
        """
        Initializes the handler class with default parameters and data containers.

        The initialization function specifies available cell types and associated
        parameters like division rate, growth rate, and adhesion forces. It also
        initializes various data containers for simulation metrics.

        :return: None
        """
        # For adhesion handler
        self.forces = []

        # Data handler: for total distance between cells - simulation stability
        self.data_file_path = ''
        self.time = None
        self.times = defaultdict()
        self.frame_cells = defaultdict()
        self.frame_interval = [None, None]
        self.mother_adhesion_strength = None
        self.falloff = None
        # self.master_dict = None

        # Data handler: for phase diagrams
        self.stable_measure = []
        self.tension = []
        self.adhesion = []

        # For long axis
        # self.vec_axis = defaultdict(list)
        # self.len_axis = defaultdict(list)

        # Data handler: for contact area 
        self.contact_ratios = defaultdict(list)
        self.contact_areas = defaultdict(list)

        # For cell division
        self.division_rate = 0.0
        self.cell_under_div = None
        self.daugthers = []

        # For cell growth 
        self.volumes = defaultdict(list)
        self.pressures = defaultdict(list)

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
            volume_scale=1,
            division_rate=50,
            start=1,
            end=250,
            motion_strength=-500,
            division_trigger='volume',
            adhesion=True,
            growth=False,
            division=False,
            data=False,
            motility=False,
            # remeshing=False,
            # visualize=False,
            boundary=False):
        """
        Launches the simulation with specified parameters.

        This method initializes and launches the simulation based on the provided 
        parameters. It sets up event handlers for various simulation aspects, 
        including adhesion, data export, growth, division, motility, and boundary 
        conditions.

        :param str filepath: The file path to save simulation data.
        :param float volume_scale: Scaling factor for cell volumes.
        :param int division_rate: The rate at which cells undergo division (in frames).
        :param int start: The starting frame of the simulation.
        :param int end: The ending frame of the simulation.
        :param int motion_strength: Strength of motion applied to cells.
        :param str division_trigger: The trigger for cell division ('volume' or \
              'random').
        :param bool adhesion: Enable/disable adhesion in the simulation.
        :param bool growth: Enable/disable cell growth in the simulation.
        :param bool division: Enable/disable cell division in the simulation.
        :param bool data: Enable/disable data export during the simulation.
        :param bool motility: Enable/disable cell motility in the simulation.
        :param bool boundary: Enable/disable boundary conditions in the simulation.

        :return: None
        """
        self.frame_interval = [start, end]
        bpy.context.scene.frame_set(start)
        bpy.context.scene.frame_start = start
        bpy.context.scene.frame_end = end

        self.data_file_path = filepath
        bpy.context.scene.render.filepath = filepath
        self.division_rate = division_rate
        self.motion_strength = motion_strength
        self.data_flag = data  # used to decide if data are computed
        self.division_trigger = division_trigger
        self.volume_scale = volume_scale
        self.unit_volume = ((4/3)*np.pi*(1)**3)

        # Set the end frame for all cloth simulation caches in the scene
        # To keep simulations running after the default 250 frames
        if bpy.context.scene.frame_current == start: 
            for collection in bpy.data.collections:
                # Loop through the objects existed in the collection 
                for obj in collection.objects:     
                    if obj.get('object') == 'cell':
                        obj.modifiers["Cloth"].point_cache.frame_start = start   
                        obj.modifiers["Cloth"].point_cache.frame_end = end   

        bpy.app.handlers.frame_change_post.clear()
        bpy.app.handlers.frame_change_post.append(self.timing_init_handler)

        if adhesion:
            bpy.app.handlers.frame_change_post.append(self.adhesion_handler)
        if data:
            bpy.app.handlers.frame_change_post.append(self.data_export_handler)
            bpy.app.handlers.frame_change_post.append(self.contact_area_handler)
        if growth:
            # bpy.app.handlers.frame_change_post.append(self.growth_handler)
            bpy.app.handlers.frame_change_post.append(self.growth_pressure_handler)
        if (division and division_trigger == 'random'):
            bpy.app.handlers.frame_change_post.append(self.division_handler_random)
        if (division and division_trigger == 'volume'):
            bpy.app.handlers.frame_change_post.append(self.division_handler_volume)            
        if motility:
            bpy.app.handlers.frame_change_post.append(self.motion_handler)
        if boundary:
            bpy.app.handlers.frame_change_post.append(self.boundary_handler)
        # if visualize:
        #    bpy.app.handlers.frame_change_post.append(self.visualize_stretching)

        bpy.app.handlers.frame_change_post.append(self.timing_elapsed_handler)
        bpy.app.handlers.frame_change_post.append(self.stop_animation)

        # bpy.ops.screen.animation_play()

        return 

    def remeshing_handler(self, scene, depsgraph): 
        """
        Handlers responsible for remeshing cells at every time step.

        It iterates over the cell objects in the scene, declares a remeshing modifier, 
        and adjusts the mesh properties.

        :param scene: The Blender scene object.
        :param depsgraph: The dependency graph object.

        :return: None
        """
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
                        remesh_mod.mode = 'VOXEL'
                        remesh_mod.voxel_size = 0.2  # microns
                        remesh_mod.adaptivity = 0 
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

    def division_handler_random(self, scene, despgraph):

        div_frames = range(
            self.frame_interval[0],
            self.frame_interval[1],
            self.division_rate
            )[1:]
        div_frames_physics = range(
            self.frame_interval[0] + 1, 
            self.frame_interval[1], 
            self.division_rate
            )[1:]
        div_frames_forces = range(
            self.frame_interval[0] + 5, 
            self.frame_interval[1], 
            self.division_rate
            )[1:]
                    
        if scene.frame_current in div_frames: 
            self.daugthers.clear()
            mother_adhesion_strength = None
            # self.force_collection = None
            self.falloff = None
            self.mother_adhesion = None
            self.mother_stiffness = None

            # Get cells that are candidates for division
            target_volume = self.unit_volume * self.volume_scale
            threshold_volume = self.volume_scale * self.unit_volume * 0.05

            candidate_cells = [
                obj for obj in bpy.context.scene.objects
                if (
                    obj.get('object') == 'cell' and 
                    obj.get('volume') is not None and 
                    abs(obj.get('volume') - target_volume) <= threshold_volume
                )
            ]

            print(f"Candidate cells: {candidate_cells}")

            if candidate_cells:
                # Select the largest cell among 5% error from target volume
                volumes = [cell['volume'] for cell in candidate_cells]
                largest_volume = np.max(volumes)
                largest_cell = candidate_cells[volumes.index(largest_volume)]
                # Replace with cell to divide    
                self.cell_under_div = largest_cell
            else:
                self.cell_under_div = None 
                
            if self.cell_under_div is not None: 
                # Get name of cell under division
                # cell_name = self.cell_under_div.name
                # Get the object by name
                dividing_cell = bpy.data.objects.get(self.cell_under_div.name)

                '''# If cell object exists
                if obj:
                    # Get cloth modifier and check if it exists
                    cloth_modifier = obj.modifiers.get('Cloth')
                    if cloth_modifier:
                        mother_stiffness = cloth_modifier.settings.tension_stiffness
                    else:
                        raise ValueError(f"No Cloth modifier found on {cell_name}")
                else:
                    raise ValueError(f"Object '{cell_name}' not found")'''
                
                apply_modifiers(obj=dividing_cell, which='all')

                d1, d2, mother_adhesion_strength, mother_motion_strength, \
                    tmp_falloff, mother_adhesion, mother_modifiers = \
                    divide_boolean(self.cell_under_div)

                # Pass on force info to daugther cells  
                # self.mother_adhesion_strength = mother_adhesion_strength
                # self.mother_motion_strength = mother_motion_strength
                # self.force_collection = tmp_collection
                # self.falloff = tmp_falloff
                self.mother_adhesion = mother_adhesion
                # self.mother_stiffness = mother_stiffness

                self.daugthers.append(d1)     
                self.daugthers.append(d2)

                print(f"-- Finishing division of {dividing_cell.name} "
                      f"at frame {scene.frame_current}")

            else: 
                print('No cells under division')

        if self.cell_under_div is not None: 

            if scene.frame_current in div_frames_physics: 

                # Toggle back on physics after mesh separation
                apply_physics(obj=self.daugthers[0], 
                              stiffness=self.mother_stiffness)
                apply_physics(obj=self.daugthers[1], 
                              stiffness=self.mother_stiffness)

                for cell in self.daugthers:
                    cells_same_type = get_cells_same_type(cell)
                    for cell_same_type in cells_same_type:
                        missing_forces = get_missing_adhesion_forces(cell_same_type, 
                                                                     cell)
                        for force in missing_forces:
                            cell.users_collection[0].objects.link(force)

            if scene.frame_current in div_frames_forces: 

                # Add adhesion forces to daugther cells
                if self.mother_adhesion: 
                    add_homo_adhesion(cell_name=self.daugthers[0].name, 
                                      strength=self.mother_adhesion_strength)
                    add_homo_adhesion(cell_name=self.daugthers[1].name, 
                                      strength=self.mother_adhesion_strength)

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
        
        '''for collection in bpy.data.collections:

            target_volume = self.unit_volume * self.volume_scale
            threshold_volume = self.volume_scale * self.unit_volume * 0.05

            candidate_cell = [
                obj for obj
                in collection.all_objects
                if (
                    obj.get('object') == 'cell' and 
                    'volume' in obj and (abs(obj.get('volume') - target_volume) 
                                         <= threshold_volume)
                )
            ]
        
            # Loop through only the cell objects in collections
            cells = [
                obj for obj 
                in bpy.data.collections.get(collection.name_full).all_objects 
                if obj.get('object') == 'cell'
            ]
            print(f"All cells: {cells}")

            # Filter cells with volume around 5% of the target volume
            target_volume = self.unit_volume * self.volume_scale
            threshold_volume = self.volume_scale * self.unit_volume * 0.05
            candidate_cell = [
                cell for cell in cells 
                if 'volume' in cell and (abs(cell['volume'] - target_volume) 
                                         <= threshold_volume)
            ]'''
        
    '''# poor results because creates overlapping meshes    
    def cast_modifier(self, scene, despgraph):
        # Iterate through all objects in the scene
        for obj in bpy.context.scene.objects:
            # Check if the object is a cell object (customize the condition as needed)
            if obj.get('object') == 'cell':  # Example condition for cell objects
                modifier_exists = False
                
                # Check if a Cast modifier already exists, if yes, update its factor
                for modifier in obj.modifiers:
                    if modifier.type == 'CAST':
                        modifier_exists = True
                        if modifier.factor < -0.5:
                            continue
                        else:
                            # Decrease factor by 0.1 (adjust as needed)
                            modifier.factor -= 0.001
                            break
                
                # If no Cast modifier exists, create a new one
                if not modifier_exists:
                    cast_modifier = obj.modifiers.new(name="Cast", type='CAST')
                    cast_modifier.factor = 0.0  # Start factor from 0'''

    def division_handler_volume(self, scene, despgraph):
        """
        Handler responsible for cell division based on volume criteria.

        This method is responsible for cell division based on volume criteria 
        during the simulation. It loop over all cells, identifies candidate 
        cells for division, selects the largest one, and then divides 
        into two daughter cells.

        :param scene: The Blender scene object.
        :param despgraph: The dependency graph object.

        :return: None
        """
    
        self.daugthers.clear()
        self.mother_adhesion_strength = None
        self.falloff = None
        self.mother_adhesion = None
        self.mother_motion = None
        self.mother_stiffness = None
        candidate_cells = []

        # Get cells that are candidates for division
        target_volume = self.unit_volume * self.volume_scale
        # threshold_volume = self.volume_scale * self.unit_volume * 0.05

        candidate_cells = [
            obj for obj in bpy.context.scene.objects
            if (
                obj.get('object') == 'cell' and 
                obj.get('volume') is not None and 
                obj.get('volume') >= target_volume
            )
        ]

        print(f"Candidate cells: {candidate_cells}")

        if candidate_cells:
            # print(f"Candidate cells: {candidate_cells}")
            # Select the largest cell among 5% error from target volume
            volumes = [cell['volume'] for cell in candidate_cells]
            largest_volume = np.max(volumes)
            largest_cell = candidate_cells[volumes.index(largest_volume)]
            print(f"Largest cell: {largest_cell}")
            # Replace with cell to divide    
            self.cell_under_div = largest_cell
            
        else:
            self.cell_under_div = None 

        if self.cell_under_div is not None: 

            print(f"Cell under division: {self.cell_under_div.name}")

            # Get name of cell under division
            # cell_name = self.cell_under_div.name
            # Get the object by name
            dividing_cell = bpy.data.objects.get(self.cell_under_div.name)

            # If cell object exists
            if dividing_cell:
                # Get cloth modifier and check if it exists
                cloth_modifier = dividing_cell.modifiers.get('Cloth')
                if cloth_modifier:
                    mother_stiffness = cloth_modifier.settings.tension_stiffness
                else:
                    raise ValueError(f"No Cloth modifier found on {dividing_cell.name}")
            else:
                raise ValueError(f"Object '{dividing_cell.name}' not found")
            
            mother_subsurf = store_subsurf_settings(dividing_cell)
            mother_cloth = store_cloth_settings(dividing_cell)
            mother_collision = store_collision_settings(dividing_cell)
            mother_remesher = store_remesh_settings(dividing_cell)
            print(f" -- Mother subsurf: {mother_subsurf}")
            print(f" -- Mother Cloth: {mother_cloth}")
            print(f" -- Mother collision: {mother_collision}")
            print(f" -- Mother remesher: {mother_remesher}")

            # Apply the whole modifer stack in sequential order    
            apply_modifiers(obj=dividing_cell, which='all')

            # Separation of mother in two distinct daughter cells
            d1, d2, mother_adhesion_strength, mother_motion_strength, \
                tmp_falloff, mother_adhesion, mother_motion, _ = \
                divide_boolean(dividing_cell)
            
            # Pass on force parameters to daugther cells  
            self.mother_adhesion_strength = mother_adhesion_strength
            self.mother_motion_strength = mother_motion_strength
            self.falloff = tmp_falloff
            self.mother_adhesion = mother_adhesion
            self.mother_motion = mother_motion
            self.mother_stiffness = mother_stiffness

            self.daugthers.append(d1)     
            self.daugthers.append(d2)

            # declare modifiers mother -> daughter
            '''declare_subsurf_settings(self.daugthers[1], mother_subsurf)
            declare_subsurf_settings(self.daugthers[0], mother_subsurf)

            declare_cloth_settings(self.daugthers[1], mother_cloth)
            declare_cloth_settings(self.daugthers[0], mother_cloth)            

            declare_collision_settings(self.daugthers[1], mother_collision)
            declare_collision_settings(self.daugthers[0], mother_collision)            

            declare_remesh_settings(self.daugthers[1], mother_remesher)
            declare_remesh_settings(self.daugthers[0], mother_remesher)'''         

            # Toggle back on physics after mesh separation
            apply_physics(obj=self.daugthers[0], 
                          stiffness=self.mother_stiffness)
            apply_physics(obj=self.daugthers[1], 
                          stiffness=self.mother_stiffness)

            # Add adhesion forces to daugther cells
            if self.mother_adhesion: 
                print('Mother adhesion detected, passing it on')
                add_homo_adhesion(cell_name=self.daugthers[0].name, 
                                  strength=self.mother_adhesion_strength)
                add_homo_adhesion(cell_name=self.daugthers[1].name, 
                                  strength=self.mother_adhesion_strength)

                for cell in self.daugthers:
                    cells_same_type = get_cells_same_type(cell)
                    for cell_same_type in cells_same_type:
                        missing_forces = get_missing_adhesion_forces(cell_same_type, 
                                                                     cell)
                        print(f"Missing adhesion forces: {missing_forces}")
                        for force in missing_forces:
                            cell.users_collection[0].objects.link(force)

            if self.mother_motion: 
                print('Mother motion detected, passing it on')
                add_motion(effector_name=self.daugthers[0].name, 
                           strength=self.mother_motion_strength)
                add_motion(effector_name=self.daugthers[1].name, 
                           strength=self.mother_motion_strength)
                
            print(f"-- Finishing division of {dividing_cell.name} "
                  f"at frame {scene.frame_current}")
    
    def growth_pressure_handler(self, scene, depsgraph):
        """
        Handler responsible for cell growth based on cell's pressure.

        This method calculates the volume and target volume of each cell, 
        adjusts the uniform pressure force based on the deviation from the 
        target volume, and updates the pressure settings accordingly.

        :param scene: The Blender scene object.
        :param depsgraph: The dependency graph object.

        :return: None
        """
        for collection in bpy.data.collections: 
            cells = bpy.data.collections.get(collection.name_full).all_objects
            for cell in cells: 
                if cell.modifiers.get('Cloth'):
                    volume = calculate_volume(cell)
                    target_volume = self.volume_scale * self.unit_volume
                    cell['target volume'] = target_volume
                    cell['volume'] = volume
                    self.volumes[f"{cell.name}"].append(volume)
                    self.pressures[f"{cell.name}"].append(
                        cell.modifiers["Cloth"].settings.uniform_pressure_force
                    )
                    volume_deviation = (target_volume - volume) / target_volume
                    print(
                        f"Volume deviation: {volume_deviation}; "
                        f"Target volume: {target_volume}; "
                        f"Volume: {volume}"
                    )
                    # Define the pressure-volume relationship 
                    previous_pressure = cell.get('previous pressure')
                    # Modify pressure based on volume deviation
                    # delta_pressure = 0.001 * (math.exp(volume_deviation) - 1)
                    # print(f"Delta pressure: {delta_pressure}")
                    new_pressure = previous_pressure + (volume_deviation
                                                        * previous_pressure*0.03)
                    print(f"New pressure: {new_pressure}")
                    # Update the cloth pressure settings
                    cell.modifiers["Cloth"].settings.uniform_pressure_force = \
                        new_pressure
                    cell['previous pressure'] = new_pressure                                                    

    # Member function to handle cell growth
    def growth_handler(self, scene, depsgraph):
        """
        Handler responsible for cell growth based on mesh scaling.

        This method calculates the volume and target volume of each cell, 
        adjusts the shrink_min parameter of the Cloth modifier based on 
        the deviation from the target volume, resulting in mesh scaling 
        and mimicking cell growth.

        :param scene: The Blender scene object.
        :param depsgraph: The dependency graph object.

        :return: None
        """
        for collection in bpy.data.collections: 
            # Loop through the cell objects in collections, excluding forces
            cells = bpy.data.collections.get(collection.name_full).all_objects
            for cell in cells: 
                if (
                    # cell.get('adhesion force') is not None or
                    # cell.get('motion force') is not None or
                    cell.modifiers.get('Cloth')
                ):

                    volume = calculate_volume(cell)
                    # target_volume = ((4/3)*np.pi*(1)**3)
                    target_volume = self.volume_scale * self.unit_volume
                    cell['target volume'] = target_volume
                    cell['volume'] = volume
                    self.volumes[f"{cell.name}"].append(volume)

                    # TODO allow for multiple growth rate, with different functions
                    volume_deviation = ((volume - target_volume) / target_volume)
                    # (f"Current volume deviation: {volume_deviation}")
                    # shrink_adjustment = 0.001 * math.tanh(50 * volume_deviation)
                    shrink_adjustment = 0.01 * (math.exp(volume_deviation) - 1)
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
        """
        Handler responsible for enforcing reflective boundaries on motion forces.

        This method ensures that motion forces are constrained within the boundaries 
        of a specified box or sphere. If a force location exceeds the box boundaries, 
        it is reflected back into the box.

        :param scene: The Blender scene object.
        :param depsgraph: The dependency graph object.

        :return: None
        """
        # Get the 'box' object, per convention
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
        # print(cells)
        msd = dict()        
        
        for collection in bpy.data.collections: 
            forces = bpy.data.collections.get(collection.name_full).all_objects
            # print(forces)
            for force in forces: 
                if (force.get('motion') is not None and force.get('motion')): 
                    # print('Entering force loop')
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
                        # print('Distribution is uniform')
                        rand_coord = Vector(np.random.uniform(
                            low=-force['distribution size'],
                            high=force['distribution size'], 
                            size=(3,)
                        ))
                        new_loc = Vector(com) + rand_coord
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
                # print(f"Cell Type: {cell}, Sorting Score: {score}")
            self.sorting_scores.update({scene.frame_current: sorting_scores_same_type})

    def contact_area_handler(self, scene, depsgraph): 

        if self.data_flag: 
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
            
                # print(self.contact_areas)

    def adhesion_handler(self, scene, depsgraph):
        """
        Handler responsible for updating adhesion forces between cells.

        This method updates the location and parameters of adhesion forces 
        between cells. Adhesion forces are adapted to the size of the 
        corresponding cells, and the force location is set to the center 
        of mass of each associated cell.

        :param scene: The Blender scene object.
        :param depsgraph: The dependency graph object.

        :return: None
        """
        force_list = Force.force_list

        for force in force_list:
            if not bpy.data.objects[force.name]['motion']: 
                # Retrieve objects of interest
                assoc_cell = force.associated_cell
                cell_obj = bpy.data.objects[assoc_cell]
                force_obj = bpy.data.objects.get(force.name)
                bpy.context.view_layer.objects.active = cell_obj

                # Update the force location to its corresponding cell's center of mass
                COM = get_centerofmass(cell_obj)
                force_obj.location = COM
                # Adapt the radius of the force field to match with cell size
                _, len_long_axis = get_long_axis_global(cell_obj)                
                force_obj.field.distance_max = (len_long_axis / 2) + 0.2

                # cell type: 
                # cells will only adhere with cells from the same collection
                cloth_modifier = bpy.data.objects[assoc_cell].modifiers["Cloth"]
                cloth_settings = cloth_modifier.settings
                cloth_collection = bpy.data.objects[assoc_cell].users_collection[0]
                cloth_settings.effector_weights.collection = cloth_collection

    def data_export_handler(self, scene, depsgraph): 
        """
        Handler responsible for exporting simulation data.

        This method writes various simulation data to JSON files at the end of 
        the simulation.

        :param scene: The Blender scene object.
        :param depsgraph: The dependency graph object.

        :return: None
        """
        self.seed = bpy.context.scene.get("seed")
        # Write the list at the end of the simulation
        if scene.frame_current == self.frame_interval[1]:
            with open(f"{self.data_file_path}_times.json", 'w') as write_file:
                write_file.write(json.dumps(self.times))
            with open(f"{self.data_file_path}_frame_cells.json", 'w') as write_file:
                write_file.write(json.dumps(self.frame_cells))
            with open(f"{self.data_file_path}_contact_ratios.json", 'w') as write_file:
                write_file.write(json.dumps(self.contact_ratios))
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
            with open(f"{self.data_file_path}_pressures.json", 'w') as write_file:
                write_file.write(json.dumps(self.pressures))
            if self.seed is not None: 
                with open(f"{self.data_file_path}_seed.json", 'w') as write_file:
                    write_file.write(json.dumps(self.seed))

            '''subprocess.run([
                "python",
                "C:\\Users\\anr9744\\Projects\\Goo\\scripts\\modules\\goo\\visualization.py",
                f"{self.data_file_path}"
            ])'''

    # not used
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
                # print(normalized_ratios)

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

    def timing_init_handler(self, scene, depsgraph): 
        """
        Initializes timing at the start of rendering.

        This method is called when the rendering starts for Frame 1.

        :param scene: The Blender scene object.
        :param depsgraph: The dependency graph object.

        :return: None
        """
        if scene.frame_current == 2: 
            self.time = datetime.now()
            print(f'Render started for Frame 1 at: {self.time}')

    def timing_elapsed_handler(self, scene, depsgraph): 
        """
        Handler responsible for timing each frame during rendering.

        This method is called for each frame after the rendering is completed. 
        It is used to calculate the computational cost of each time step during 
        rendering.

        :param scene: The Blender scene object.
        :param depsgraph: The dependency graph object.

        :return: None
        """
        frame_written = scene.frame_current
        if frame_written not in [1, 2]: 
            elpased_time = datetime.now() - self.time
            elapsed_time_secs = elpased_time.seconds + elpased_time.microseconds/1000000
            self.times.update({frame_written: elapsed_time_secs})
            print('________________________________________________________')
            # print(f"Render Started at:{self.time}")  
            print(f"Current frame: {frame_written}")
            print(f"Elapsed in seconds/microseconds:{elapsed_time_secs:.3f};"
                  f" {elapsed_time_secs:.1f}")
            
        cells = [
            obj 
            for obj in bpy.data.objects 
            if "object" in obj.keys() and obj["object"] == "cell"
            ]
        
        for cell in cells: 
            if cell.name not in self.frame_cells:
                if frame_written in [1, 2]:
                    self.frame_cells[cell.name] = [1]  # Initialize list if not present
                else:
                    self.frame_cells[cell.name] = []  # Initialize list if not present
            self.frame_cells[cell.name].append(frame_written)

            if frame_written == self.frame_interval[1] - 1:
                break 

    def stop_animation(self, scene, depsgraph):
        """
        Handler responsible to stop the simulation.

        This method is called for each frame after the rendering is completed, 
        and is executed only at the last time step. 
        May be toggled to automatically quit Blender when the simulation is completed, 
        which is needed when scanning parameter spaces. 

        :param scene: The Blender scene object.
        :param depsgraph: The dependency graph object.

        :return: None
        """
        # checks if the simulation has ended
        if scene.frame_current == self.frame_interval[1]:
            print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _')
            print("The simulation has ended.")
            # True enables the last frame not to be repeated
            bpy.ops.screen.animation_cancel(restore_frame=True) 
            # closes Blender then
            # bpy.ops.wm.quit_blender()
