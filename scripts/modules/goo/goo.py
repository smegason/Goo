# goo.py - This is the Goo library. It contains all helper functions for goo
# goo is licensed under BSDv2

import bpy
import bmesh
import mathutils
import numpy as np
import sys

"""
Refactored by Antoine Ruzette, November 2022.

Note on comments: 
Core functions are required for ``Goo`` to correctly run. 
Auxilliary function are not required for ``Goo`` to run but are used for supporting tasks such as retrieving metrics - typically they can be removed. 

Sphynx docstring was enhanced by adding types, and written in a more sphynx-ic way. 
"""

# TODO: explain differences between edit and object mode
# TODO: explain Blender's modifiers as a data type

def calculate_volume(obj):
    """Auxilliary function: calculates the volume of the Blender mesh. 

    In order to retrieve the mesh as it is currrently evaluated - including 
    the effect of all modifiers - in the simulation, its corresponding evaluated 
    ID is obtained from the dependency graph for the current context. 
    .. seealso:: [Blender API Documentation > ``evaluated_get(depsgraph)``](https://docs.blender.org/api/current/bpy.types.ID.html?highlight=evaluated_get#bpy.types.ID.evaluated_get)

    :param bpy.data.objects['name'] obj: The Blender mesh.
    :returns: The volume of the mesh. 
    :rtype: float

    .. note:: The function may return a negative volume - see ``calc_volume(signed=True)``. 

    """

    # We need to get the cell as it is evaluated in the simulation.
    # To do this, we fetch its dependency graph and obtain the
    # evaluated cell (denoted as obj_eval here)
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    # We use this graph to obtain a new mesh from which
    # we can calculate the volume using Blender's bmesh library
    mesh_from_eval = obj_eval.to_mesh()
    bm = bmesh.new()
    bm.from_mesh(mesh_from_eval)
    volume = bm.calc_volume(signed=True)
    print(type(volume))
    return volume


def get_long_axis(obj):
    """Core function: calculates the long axis of a mesh. 

    This function calculates the first eigenvector of the vertices in the mesh, 
    which corresponds to the long axis.

    :param bpy.data.objects['name'] obj: The Blender mesh.
    :returns: The coordinates of the long axis of the mesh as a tuple(x, y, z) which gives direction from the origin (0, 0, 0). 
    :rtype: tuple
    """

    # We need to get the cell as it is evaluated in the simulation.
    # To do this, we fetch its dependency graph and obtain the
    # evaluated cell (denoted as obj_eval here)
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    # We obtain the (x, y, z) coordinates of the vertices in the
    # evaluated cell
    vertices = obj_eval.data.vertices
    vert_coords = [(obj_eval.matrix_world @ v.co) for v in vertices]
    vert_coords = np.asarray(vert_coords)

    # We separate the x, y, and z coordinates into their own arrays.
    # We also subtract the mean of each dimension from the corresponding
    # array values. This is part of the PCA algorithm.
    x = vert_coords[:, 0]
    x = x - np.mean(x)
    y = vert_coords[:, 1]
    y = y - np.mean(y)
    z = vert_coords[:, 2]
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
    major_axis = (major_x, major_y, major_z)
    return major_axis


def get_division_angles(axis):
    """Auxilliary function: retrieves the 2 angles associated with the long axis of a cell. 

    :param tuple axis: Coordinates of the long axis to the division plane - as retrived by ``get_major_axis(obj)``.
    :returns: 
        - phi (:py:class:`tuple`) - Phi is the angle between the division axis and the z-axis in spherical coordinates. 
        - theta (:py:class:`tuple`) - Theta is the angle projected on the xy plane in cartesian 3D coordinates. 
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

# TODO: to be removed 
'''
def calculate_contact_area(obj):  # obj = bpy.context.object
    """
    Auxilliary function: calculates the contact area of a cell with all surrounding objects
    by looking for flat areas of the mesh. Not always a good assumption
    :param obj: the blender object (mesh) -> obj = bpy.data.objects['object name']
    :returns: contact_area
    """

    # We need to get the cell as it is evaluated in the simulation.
    # To do this, we fetch its dependency graph and obtain the
    # evaluated cell (denoted as obj_eval here)
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    # We use this graph to obtain a new mesh from which
    # we can calculate the angles between adjacent faces of the cell
    mesh_from_eval = obj_eval.to_mesh()
    bm = bmesh.new()
    bm.from_mesh(mesh_from_eval)

    # The contact area variable is initialized to zero.
    contact_area = 0
    # All faces on the mesh are initially deselected
    for face in bm.faces:
        face.select = False

    # We loop over all the edges in the mesh
    for edge in bm.edges:
        # We find the angle between the two adjacent faces
        # of the edge. If the angle is less than 0.01, the
        # adjacent faces are either flat or concave, so this
        # face is presumed to be in contact with another object
        # and is selected
        angle = edge.calc_face_angle_signed()
        if angle < 0.01:
            for face in edge.link_faces:
                face.select = True

    # We loop through the selected faces and add the area of
    # each face to the contact area variable
    for face in bm.faces:
        if face.select:
            contact_area += face.calc_area()

    return contact_area
'''

def repair_hole(obj):

    """Core function: repair the holes created by the bissection of the Blender mesh. 
    
    This function adds a new face to the cell after division (bissection - in two equal parts) 
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
    """Auxilliary function: Retrieves the number of vertices, edges and faces of a Blender object

    :param bpy.data.objects['name'] obj: The Blender mesh. 
    :returns: The number of vertices, faces and edges of ``obj``. 
    :rtype: list of float
    """
    vertices = obj.data.vertices
    edges = obj.data.edges
    faces = obj.data.polygons

    print(f'Vertices: {len(vertices)}; Edges: {len(edges)}; Faces: {len(faces)}')
    
    return [vertices, edges, faces]


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
    :param tuple major_axis: The major axis of the mother mesh specified by XYZ coordinates from the origin. 
    :returns: A boolean flag indicating which cell has to be translated in ``translate_cell(obj, axis)``. 
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
    :param tuple axis: The XYZ coordinates of the major axis of the Blender mesh to translate.
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


def seperate_cell(obj):
    """Core function: splits one cell into two for cell division. 
    
    The mother Blender mesh is split in two given a division plane. 
    The division plane is specified by a point on this plane, that is the center of mass of the mother,
    and the direction of this point, given by the long axis of the mother mesh. 
    By naming convention, the daughter cells inherit their mother name added with .001. 

    :param bpy.data.objects['name'] obj: the Blender mesh to divide in two. 
    :returns: 
        - mother_name (:py:class:`str`) - The name of the mother cell. 
        - daughter_name (:py:class:`str`) - The name of the daughter cell. 
        - COM (:py:class:`tuple`) - The XYZ coordinates of the center of mass of the mother cell. 
        - major_axis (:py:class:`tuple`) - The XYZ coordinates of the long axis of the mother cell. 
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
    major_axis = get_long_axis(obj)
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
    # We create a list of the vertex indices we will use to create the daughter cell object
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

def divide(obj): 
    """Core function: divides a mother Blender mesh into two daughter Blender meshes. 
    
    The division plane is orthogonal to the major axis of the parent mesh. 
    Additionally, this function may translate daughter cells, 
    translating one mesh, filling the two holes, and retriangulating the meshes

    :param bpy.data.objects['mother_name'] obj: The soon-to-be mother Blender mesh. 
    :returns: 
        - daughter1 (bpy.data.objects['daughter1_name']) - The Blender mesh of one of the daughter cells. 
        - daughter2 (bpy.data.objects['daughter2_name']) - The Blender mesh of the other daughter cell. 
    """
    # Select mother cell
    obj.select_set(True)
    # Split mother's Blender mesh in two
    m_name, d_name, COM, major_axis = seperate_cell(obj)
    # Decides which ce
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
    return daughter1, daughter2

def turn_off_physics():
    """Core function: turns physics off for the currently selected Blender object.

    The cloth physics for cells are turned off before division to avoid irregular mesh behavior after. 

    :param: None
    :returns: None
    """
    bpy.ops.object.modifier_remove(modifier="Cloth")


def turn_on_physics():
    """Core function: turns physics on for the currently selected Blender object.
    
    Physics are turned back on after cell division occurs to avoid irregular mesh behavior post-division.

    :param: None
    :returns: None
    """
    # Set all parameter values
    bpy.ops.object.modifier_add(type='CLOTH')
    bpy.context.object.modifiers["Cloth"].settings.bending_model = 'LINEAR'
    bpy.context.object.modifiers["Cloth"].settings.quality = 5
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


# TODO: to be removed as redundant with the div_handler of the handler class
'''
def mitosis_handler(scene):
    """
    Mitosis handler is a function triggered every frame to check if all cells
    need to divide
    Currently just based on volume
    :param scene: Blender syntax
    :param tree: ???
    :return: None
    """
    # TODO check docstring
    # TODO should merge mitosis_handler and div_handler into one function.
    # When to divide should be based on a cell property that allows for volume
    # or time based trigger and allows for some noise in timing

    # Find the number of cells by checking the length of the collection
    # where they are kept in the Blender simulation
    num_cells = len(bpy.data.collections["Cells"].objects)
    cell_tree = None
    # Loop through each cell
    for i in range(num_cells):
        # Make that cell the active object in the simulation
        bpy.context.view_layer.objects.active = bpy.data.collections["Cells"].objects[i]
        cell = bpy.context.view_layer.objects.active
        # Calculate the cell volume
        volume = calculate_volume(cell)
        # If the volume is above a certain threshold, divide
        if volume > .3:
            cell_tree = divide(cell, cell_tree)
    cell_tree.show()
'''

# TODO: to be removed as redundant 
'''
# Remove this function. Take try/except code to make_cell function
def make_mesh(cell):
    """
    Makes a Blender mesh object for the given ``Goo`` cell
    :param cell: the ``Goo`` cell to make a Blender mesh from
    :return: None
    """
    # Attempt to add a RoundCube mesh
    try:
        bpy.ops.mesh.primitive_round_cube_add(change=False,
                                              radius=cell.data['radius'],
                                              size=cell.data['size'],
                                              arc_div=cell.data['arcdiv'],
                                              lin_div=0, div_type='CORNERS',
                                              odd_axis_align=False, no_limit=False,
                                              location=cell.data['location'])
    # If the user does not have this object enabled, throw exception
    except Exception:
        print("Exception=" + sys.exc_info())
        print("To enable RoundCube creation for Cells you must go \
               to Edit->Preferences->AddOns->Add Mesh:ExtraObjects\
               and check the box to enable it")

    # Get the cell name
    bpy.context.object.name = cell.data['name']

    # Name the Blender mesh the cell's name
    bpy.context.view_layer.objects.active = bpy.data.objects[cell.data['name']]

    # Add a subdivision surface for a smoother shape
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = cell.data['subdiv']
'''

def make_cell(cell, force = None):
    # add mesh type as an cell argument? 
    """Core function: creates a Blender mesh corresponding to a Goo :class:`Cell` object. 

    :param :class:`Cell` cell: The Goo :class:`Cell` object. 
    :returns: None
    """

    # Add a round_cube mesh
    if cell.data['mesh_type'] == "round_cube" or "":
        print("Making mesh as a Round Cube")
        try:
            bpy.ops.mesh.primitive_round_cube_add(change=False,
                                                  radius=cell.data['radius'],
                                                  size=cell.data['size'],
                                                  arc_div=cell.data['arcdiv'],
                                                  lin_div=0,
                                                  div_type='CORNERS',
                                                  odd_axis_align=False,
                                                  no_limit=False,
                                                  location=cell.data['location'])
        except Exception:
            print(sys.exc_info())
            print(f"To enable RoundCube creation for Cells you must go to \n \
                    Edit->Preferences->AddOns->Add Mesh:ExtraObjects and check the box to enable it")
    
    else: 
        print("Making mesh as Ico Sphere")
        try:
            bpy.ops.mesh.primitive_ico_sphere_add(enter_editmode=False,
                                                  align='WORLD',
                                                  location=cell.data['location'],
                                                  scale=cell.data['scale'],
                                                  radius=cell.data['radius'])

        except Exception:
            print(sys.exc_info())
            print(f"{cell.data['mesh_type']} is not recognised as a valid mesh type.")        

    # Give the Blender object the cell's name
    obj = bpy.context.object
    bpy.context.object.name = cell.data['name']
    bpy.context.view_layer.objects.active = bpy.data.objects[cell.data['name']]

    # Smooth the mesh
    bpy.ops.object.select = True
    bpy.ops.object.shade_smooth()

    # Add subsurface modifier to make smoother
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = cell.data['subdiv']

    # Add cloth settings for physics
    bpy.ops.object.modifier_add(type='CLOTH')
    bpy.context.object.modifiers["Cloth"].settings.bending_model = 'LINEAR'
    bpy.context.object.modifiers["Cloth"].settings.quality = 5
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1
    bpy.context.object.modifiers["Cloth"].settings.mass = cell.data['vertex_mass']
    obj.modifiers["Cloth"].settings.air_damping = cell.data['air_damping']
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
    bpy.context.object.collision.use_culling = True
    bpy.context.object.collision.damping = 0.579821
    bpy.context.object.collision.thickness_outer = 0.02
    bpy.context.object.collision.thickness_inner = 0.2
    bpy.context.object.collision.cloth_friction = 5

    '''    
    if force == None: 
        # FIXME: add forces to daughter cells!!
        bpy.ops.object.forcefield_toggle()
        bpy.context.object.field.type = 'FORCE'
        bpy.context.object.field.strength = -20
        bpy.context.object.field.shape = 'POINT'
    elif type(force) == Force(): 
        cell = force.associated_cell # Link forces to cells in Python/Goo
        bpy.ops.object.forcefield_toggle()
        bpy.context.object.field.type = force.type
        bpy.context.object.field.name = force.name
        bpy.context.object.field.falloff_power = force.falloff_power
        bpy.context.object.field.strength = force.strength
        bpy.context.object.field.falloff_type = 'SPHERE'
        bpy.context.object.field.shape = 'SURFACE'

        # add scale (1,1,1)

    '''

    # add material to cell based on name of material
    material = bpy.data.materials.get(cell.data['material'])
    if (material):
        bpy.context.active_object.data.materials.append(material)
    else:
        print(f"No material: {cell.data['material']}")

# TODO: to be removed - not used
'''
def delete_cell(cell):
    """
    Auxilliary function: deletes a Goo cell's Blender object

    :param cell: the Goo cell to delete
    :return: None
    """
    obj = cell.get_blender_object()
    bpy.data.objects.remove(obj, do_unlink=True)
'''

# Defines the Cell class
class Cell():
    '''Core class: creates Goo cell object. 

    :param str name: The name of the cell.
    :param tuple loc: The coordinates of the cell.   
    :returns None:
    '''
    def __init__(self, name, loc, mesh_type='', material=''):
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
            'scale': (1, 1, 1),
            'arcdiv': 8,
            'subdiv': 2,
            'vertex_mass': 0.3,
            'density': 1.0,
            'edges_pull': 0.0,
            'edges_bend': 0.2,
            'air_damping': 10,
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
            'mesh_type': mesh_type
        }
        # The volume and mass are calculated from values in the data dictionary
        self.data['init_volume'] = ((4/3)*np.pi*(self.data['radius'])**3)
        self.data['mass'] = self.data['density']*self.data['init_volume']

    # Member function for obtaining the Blender object corresponding to the cell
    def get_blender_object(self):
        obj = bpy.data.objects[self.data["name"]]
        return obj


# TODO: to be removed as redundant with the standalone function
'''
    # Member function to divide
    def divide(self):
        obj = self.get_blender_object()
        obj.select_set(True)
        m_name, d_name, COM, major_axis = seperate_cell(obj)
        print(major_axis)
        val = select_translate_cell(bpy.data.objects[m_name], COM, major_axis)
        if val == 0:
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
        return daughter1, daughter2
'''

def add_material_cell(mat_name, r, g, b):
    """Core function: creates a soap bubble-like Blender material for use in rendering cells.

    The material has a name that allows it to be shared across multiple cells. 

    :param str mat_name: The name of the material. 
    :param float r: The value of the red in RGB [0 to 1]. 
    :param float g: The value of the green value in RGB [0 to 1]. 
    :param float b: The value of the blue value in RGB [0 to 1]. 
    :returns: None
    """

    # check whether the material already exists
    if bpy.data.materials.get(mat_name):
        mat = bpy.data.materials[mat_name]
    else:
        # create the material
        mat = bpy.data.materials.new(mat_name)

    mat.diffuse_color = (1.0, 1.0, 1.0, 1.0)  # viewport color
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
    node_main.inputs['Base Color'].default_value = (r, g, b, 1)
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
    node_random.inputs['Base Color'].default_value = (r, g, b, 1)  # RGB or HSV?? todo
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
    def __init__(self, force_name, cell_name, strength, falloff_power):
        self.name = force_name
        self.strength = strength
        self.associated_cell = cell_name
        self.falloff_power = falloff_power
        self.falloff_type = 'SPHERE'
        self.shape = 'SURFACE'


def make_force(force):
    """Core function: creates a Blender force from a Goo :class:`Force` object. 

    :param :class:`Force` force: The Goo force object. 
    :returns: None
    """
    # Add a force object
    cell = force.associated_cell
    bpy.ops.object.effector_add(type='FORCE',
                                enter_editmode=False,
                                align='WORLD',
                                location=bpy.data.objects[cell].location,
                                scale=(1, 1, 1))

    # Add force parameters
    bpy.context.object.field.strength = force.strength
    bpy.context.object.name = force.name
    bpy.context.object.field.falloff_power = force.falloff_power
    # add links for falloff_type and shape

# not used - under development
'''
def initialize_cell_sheet():
    bpy.ops.mesh.primitive_grid_add(size=2,
                                    enter_editmode=False,
                                    align='WORLD',
                                    location=(0, 0, 0),
                                    scale=(8.28558, 8.28558, 8.28558))
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.transform.shear(value=-0.5,
                            orient_axis='Z',
                            orient_axis_ortho='Y',
                            orient_type='GLOBAL',
                            orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                            orient_matrix_type='GLOBAL',
                            mirror=True,
                            use_proportional_edit=False,
                            proportional_edit_falloff='SMOOTH',
                            proportional_size=1,
                            use_proportional_connected=False,
                            use_proportional_projected=False,
                            release_confirm=True)
    bpy.ops.object.mode_set(mode='OBJECT')

    bpy.ops.transform.resize(value=(8.28558, 8.28558, 8.28558),
                             orient_type='GLOBAL',
                             orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                             orient_matrix_type='GLOBAL',
                             mirror=False,
                             use_proportional_edit=False,
                             proportional_edit_falloff='SMOOTH',
                             proportional_size=1,
                             use_proportional_connected=False,
                             use_proportional_projected=False)

    c = Cell("cell", loc=(0, 0, 0))
    make_cell(c)
    bpy.data.objects["Grid"].select_set(False)
    bpy.data.objects["cell"].select_set(False)
    bpy.data.objects["cell"].select_set(True)
    bpy.data.objects["Grid"].select_set(True)
    bpy.ops.object.parent_set(type='OBJECT')
    bpy.context.object.instance_type = 'VERTS'
    bpy.ops.object.duplicates_make_real()
    bpy.ops.outliner.item_activate(deselect_all=True)
    bpy.data.objects["Grid"].select_set(True)
    bpy.ops.object.delete(use_global=False)


def initialize_cell_shell():
    bpy.ops.mesh.primitive_ico_sphere_add(radius=1,
                                          enter_editmode=False,
                                          align='WORLD',
                                          location=(0, 0, 0),
                                          scale=(4.16825, 4.16825, 4.16825))
    c = Cell("cell", loc=(0, 0, 0))
    make_cell(c)
    bpy.data.objects["Icosphere"].select_set(False)
    bpy.data.objects["cell"].select_set(False)
    bpy.data.objects["cell"].select_set(True)
    bpy.data.objects["Icosphere"].select_set(True)
    bpy.ops.object.parent_set(type='OBJECT')
    bpy.context.object.instance_type = 'VERTS'
    bpy.ops.object.duplicates_make_real()
    bpy.ops.outliner.item_activate(deselect_all=True)
    bpy.data.objects["Icosphere"].select_set(True)
    bpy.ops.object.delete(use_global=False)
    bpy.data.objects["cell"].select_set(True)
    bpy.ops.object.delete(use_global=False)


def initialize_solid_tissue():  # Work in Progress
    c1 = Cell("cell1", loc=(0, 0, 0))
    make_cell(c1)
    c2 = Cell("cell2", loc=(1, 1, -1))
    make_cell(c2)
    c3 = Cell("cell3", loc=(-1, 1, 1))
    make_cell(c3)
    c4 = Cell("cell4", loc=(-1, -1, -1))
    make_cell(c4)
    c5 = Cell("cell5", loc=(1, -1, 1))
    make_cell(c5)
    bpy.data.objects["cell1"].select_set(True)
    bpy.data.objects["cell2"].select_set(True)
    bpy.data.objects["cell3"].select_set(True)
    bpy.data.objects["cell4"].select_set(True)
    bpy.data.objects["cell5"].select_set(True)
    bpy.ops.object.duplicate_move(
       OBJECT_OT_duplicate={"linked": False,
                            "mode": 'TRANSLATION'},
       TRANSFORM_OT_translate={"value": (-0, -2.33693, -0),
                               "orient_axis_ortho": 'X',
                               "orient_type": 'GLOBAL',
                               "orient_matrix": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                               "orient_matrix_type": 'GLOBAL',
                               "constraint_axis": (False, True, False),
                               "mirror": False,
                               "use_proportional_edit": False,
                               "proportional_edit_falloff": 'SMOOTH',
                               "proportional_size": 1,
                               "use_proportional_connected": False,
                               "use_proportional_projected": False,
                               "snap": False,
                               "snap_target": 'CLOSEST',
                               "snap_point": (0, 0, 0),
                               "snap_align": False,
                               "snap_normal": (0, 0, 0),
                               "gpencil_strokes": False,
                               "cursor_transform": False,
                               "texture_space": False,
                               "remove_on_cancel": False,
                               "view2d_edge_pan": False,
                               "release_confirm": False,
                               "use_accurate": False,
                               "use_automerge_and_split": False})
'''
def setup_world():
    """Auxilliary function: sets up the default values used for simulations in Goo 
    including units and rendering background. 

    :returns: None
    """
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

    # Addn an HDRI image for illumination
    add_world_HDRI()


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
            print("update view 3d to rendered")
            space = area.spaces.active
            if space.type == 'VIEW_3D':
                space.shading.type = 'RENDERED'


def render(file_path, scene, start, end):
    """Auxilliary function: renders a simulation to create a set of still images that can be made into a movie

    :param str file_path: The path of the folder used to store output images. 
    :param bpy.context.scene scene: The Blender current scene. 
    :param int start: The Blender starting frame. 
    :param int end: The Blender ending frame. 
    :returns: None
    """
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

# collections are currently manually created in Blender - this function is not being used
def make_force_collections(master_collection, cell_types):
    """Auxilliary function: makes collections for forces to be stored in.

    :param bpy.context.view_layer.active_layer_collection master_collection: The collection in which the force
    collections will be contained. 
    :param cell_types: list of active cell types
    :returns: None
    """
    bpy.context.view_layer.active_layer_collection = master_collection
    for type in cell_types:
        collection = bpy.context.blend_data.collections.new(name=type+"_forces")
        bpy.context.collection.children.link(collection)


class handler_class:
    # TODO document this class
    """
    A class for creating different types of handlers that trigger
    actions on ``Goo`` cells when certain criteria are met
    """

    # The initialization function specifies available cell types and associated
    # parameters like division rate, growth rate, and adhesion forces
    def __init__(self):

        self.cell_types = ['sphere', 'type1', 'type2']
        self.division_rates = {}
        self.growth_rates = {}
        self.adhesion_forces = {}  # dictionary of dictionaries
        # ex: {"sphere": {"sphere": 100, "type1": 100, "type2": 100},
        #      "type1": {"sphere": 200, "type1": 200, "type2": 200},
        #      "type2": {"sphere": 300, "type1": 300, "type2": 300}}
        # in the example above, spheres exert 100 force on other
        # spheres, 100 force on type 1, and 100 on type  2
        # could change values to have different amounts of force on different cell types
        # Set parameter values for each cell type
        for type in self.cell_types:
            self.division_rates[type] = 0
            self.growth_rates[type] = 0
            self.adhesion_forces[type] = {}
            for i in self.cell_types:
                self.adhesion_forces[type][i] = 0
        # Set active (dividing) cell types
        # add active types to know what collections to divide
        self.active_cell_types = []
        self.forces = []
        return

    # Member function to set division rate for a cell type
    def set_division_rate(self, cell_type, rate):
        # assume 60 frames per second
        # rate is in divisions per second
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
        if(rate > 0 and rate < 1):
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

    # Member function to set the division handler for cells
    def div_handler(self, scene, depsgraph):
        # Loop through active cell types
        for cell_type in self.active_cell_types:
            # Get the number of cells of that type
            num_cells = len(bpy.data.collections[cell_type].objects)
            # Get the currrernt frame number
            current_frame = bpy.data.scenes[0].frame_current
            if self.division_rates[cell_type] == 0:
                continue
            # Divide based on the division rate
            if current_frame % self.division_rates[cell_type] == 0:
                print("DIVIDING CELLS")

                # Loop over all the cells of a type
                for i in range(num_cells):
                    # Get the cell name
                    cell_name = bpy.data.collections[cell_type].objects[i].name

                    # Get the corresponding Blender object
                    cell = bpy.data.objects[cell_name]

                    # get scale before division
                    # scale = cell.modifiers["Cloth"].settings.shrink_min

                    # Select the Blender object and make it the active object
                    bpy.data.objects[cell_name].select_set(True)
                    bpy.context.view_layer.objects.active = bpy.data.objects[cell_name]
                    print(cell_name)

                    # Turn off the cloth physics for this cell. This mitigates
                    # irregular mesh behavior after division
                    turn_off_physics()

                    # Divide the cell
                    d1, d2 = divide(cell)

                    # Select the first daughter cell only
                    vl = bpy.context.view_layer
                    bpy.data.objects[d1.data["name"]].select_set(True)
                    bpy.data.objects[d2.data["name"]].select_set(False)
                    vl.objects.active = bpy.data.objects[d1.data["name"]]

                    # Turn on the physics for this daughter cell
                    turn_on_physics()

                    # Select the second daughter cell only
                    bpy.data.objects[d1.data["name"]].select_set(False)
                    bpy.data.objects[d2.data["name"]].select_set(True)
                    vl.objects.active = bpy.data.objects[d2.data["name"]]

                    # Turn on the physics for this daughter cell
                    turn_on_physics()
                    bpy.data.objects[d2.data["name"]].select_set(False)

    # Member function to handle cell growth
    def growth_handler(self, scene, depsgraph):
        for cell_type in self.active_cell_types:
            num_cells = len(bpy.data.collections[cell_type].objects)
            for i in range(num_cells):
                cell_name = bpy.data.collections[cell_type].objects[i].name
                cell = bpy.data.objects[cell_name]
                cell.modifiers["Cloth"].settings.shrink_min -= 0.01

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

    def adhesion_handler(self, scene, depsgraph):
        for force in self.forces:
            assoc_cell = force.associated_cell
            bpy.context.view_layer.objects.active = bpy.data.objects[assoc_cell]
            dg = bpy.context.evaluated_depsgraph_get()
            cell_eval = bpy.data.objects[assoc_cell].evaluated_get(dg)
            vertices = cell_eval.data.vertices
            vert_coords = [(cell_eval.matrix_world @ v.co) for v in vertices]
            vert_coords = np.asarray(vert_coords)

            x = vert_coords[:, 0]
            y = vert_coords[:, 1]
            z = vert_coords[:, 2]
            COM = (np.mean(x), np.mean(y), np.mean(z))
            bpy.data.objects[force.name].location = COM