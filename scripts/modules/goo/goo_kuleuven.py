import bpy
import bmesh
import mathutils
import numpy as np
import sys


# Before all, please note that there are two kinds of object in the follwing codes,
# one is Python object created by Goo package, another is the Blender mesh object.

# Note: Open the cmd to run blender, so that you can see the print information and error 
# from cmd, it is more clear.
def setup_world():
    """
    Sets up the default values used for simulations in Goo
    including units and rendering background

    :return: None
    """
    # print("-------------setup world--------------")
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

    # Clear frames
    bpy.app.handlers.frame_change_post.clear()

    # Delete all existing collections 
    for collection in bpy.data.collections:  # loop through the existing collection
        # loop through objects in collection
        for objs in collection.objects:
            # delete existing objects in collection 
            bpy.data.objects.remove(objs)
        # Delete collection
        bpy.data.collections.remove(collection)


def add_world_HDRI():
    """
    Sets up Blender World properties for use in rendering most importantly
    adds an HDRI image for illumination

    :return: None
    """
    print("add world HDRI")
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


class Cell():
    # It is a class to make cell object
    def __init__(self, name_string, cell_type, radius, shape, loc, material=""):
        # The initialization function sets all parameters of a cell
        self.data = {
            'name': name_string,
            'cell_type': cell_type,
            'radius': radius,
            'location': loc,
            'material': material,
            'shape': shape,
            'mother': 'none',
            'daughters': ['none', 'none'],
        }
        # The volume and mass are calculated from values in the data dictionary
        self.data['init_volume'] = ((4/3)*np.pi*(self.data['radius'])**3)


def calculate_contact_area(obj):  # obj = bpy.context.object
    """
    Calculates the contact area of a cell with all surrounding objects
    by looking for flat areas of the mesh. Not always a good assumption

    :param obj: the Blender mesh ??

    :return: contact_area
    """
    # TODO check docstring - I don't understand this function. Doesn't it
    # need two cells to calculate the contact area??? Contact area with what?

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


def render(file_path, scene, start, end):
    """
    Renders a simulation to create a set of still images that can be made into a movie

    :param file_path: path for storing outputted images
    :param scene: the Blender scene
    :param start: the Blender start key frame
    :param end: the Blender end key frame

    :return: None
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


def make_cell(cell):
    """
    Makes a Blender mesh object for the given Goo cell Python object

    :param cell: the Goo cell Python object used to make a Blender mesh

    :return: None
    """
    # Use cell_type to name a new collection in Blender,
    # and link the new Blender mesh under this collection
    make_new_collection = True
    for collection in bpy.data.collections:
        if collection.name == cell.data['cell_type']:
            bpy.context.view_layer.active_layer_collection = \
                bpy.context.view_layer.layer_collection.children[cell.data['cell_type']]
            make_new_collection = False

    if make_new_collection is True:
        collection = bpy.data.collections.new(cell.data['cell_type'])
        bpy.context.scene.collection.children.link(collection)
        bpy.context.view_layer.active_layer_collection = \
            bpy.context.view_layer.layer_collection.children[cell.data['cell_type']]

    # Create a ico_sphere mesh object if the shape of cell is "sphere"
    if cell.data['shape'] == "sphere":
        try:
            bpy.ops.mesh.primitive_ico_sphere_add(enter_editmode=False,
                                                  align='WORLD',
                                                  location=cell.data['location'],
                                                  scale=(1, 1, 1),
                                                  radius=cell.data['radius'])

        except Exception:
            print(sys.exc_info())
            print("To enable RoundCube creation for Cells you must go to ")
            print("Edit->Preferences->AddOns->Add Mesh:ExtraObjects and ")
            print("check the box to enable it")
            return

    # Give the Blender object the cell's name
    obj = bpy.context.object
    bpy.context.object.name = cell.data['name']
    bpy.context.view_layer.objects.active = bpy.data.objects[cell.data['name']]

    # Smooth the mesh
    bpy.ops.object.select = True
    bpy.ops.object.shade_smooth()

    # Add subsurface modifier to make smoother
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = 2

    # Add cloth settings for physics
    bpy.ops.object.modifier_add(type='CLOTH')
    obj.modifiers['Cloth'].settings.quality = 6
    obj.modifiers['Cloth'].settings.air_damping = 5
    obj.modifiers['Cloth'].settings.bending_model = 'ANGULAR'
    obj.modifiers['Cloth'].settings.tension_stiffness = 1
    obj.modifiers['Cloth'].settings.compression_stiffness = 1
    obj.modifiers['Cloth'].settings.shear_stiffness = 1
    obj.modifiers['Cloth'].settings.bending_stiffness = 1
    obj.modifiers['Cloth'].settings.tension_damping = 25
    obj.modifiers['Cloth'].settings.compression_damping = 25
    obj.modifiers['Cloth'].settings.shear_damping = 25
    obj.modifiers['Cloth'].settings.bending_damping = 0.5
    obj.modifiers['Cloth'].settings.use_internal_springs = True
    obj.modifiers['Cloth'].settings.internal_spring_max_length = 0
    obj.modifiers['Cloth'].settings.internal_spring_max_diversion = 0.785398
    obj.modifiers['Cloth'].settings.internal_spring_normal_check = True
    obj.modifiers['Cloth'].settings.internal_tension_stiffness = 0
    obj.modifiers['Cloth'].settings.internal_compression_stiffness = 1
    obj.modifiers['Cloth'].settings.uniform_pressure_force = 5
    obj.modifiers['Cloth'].settings.use_pressure_volume = True
    obj.modifiers['Cloth'].settings.target_volume = 1
    obj.modifiers['Cloth'].settings.pressure_factor = 1
    obj.modifiers['Cloth'].settings.fluid_density = 0
    obj.modifiers['Cloth'].collision_settings.collision_quality = 6
    obj.modifiers['Cloth'].collision_settings.use_collision = True
    obj.modifiers['Cloth'].collision_settings.use_self_collision = True
    obj.modifiers['Cloth'].collision_settings.self_friction = 5
    obj.modifiers['Cloth'].collision_settings.self_distance_min = 0.015
    obj.modifiers['Cloth'].collision_settings.self_impulse_clamp = 0

    # add Collision modifier for physics
    bpy.ops.object.modifier_add(type='COLLISION')
    # bpy.context.object.collision.use_culling = False
    obj.modifiers['Collision'].settings.damping = 0.579821
    obj.modifiers['Collision'].settings.thickness_outer = 0.02
    obj.modifiers['Collision'].settings.thickness_inner = 0.2
    obj.modifiers['Collision'].settings.cloth_friction = 5
    obj.modifiers['Collision'].settings.use_culling = True

    # add material to cell based on name of material
    material = bpy.data.materials.get(cell.data['material'])
    if (material):
        bpy.context.active_object.data.materials.append(material)
    else:
        print("No material ", cell.data['material'])

    # Add Cast modifier to control the shape of Blender object
    # Factor 1 will force the Blender object into a sphere
    bpy.ops.object.modifier_add(type='CAST')
    bpy.context.object.modifiers["Cast"].factor = 1
    bpy.context.object.modifiers["Cast"].show_in_editmode = True
    bpy.context.object.modifiers["Cast"].show_on_cage = True


def delete_cell(cell):
    """
    Deletes a Goo cell's Blender object

    :param cell: the Goo cell to delete

    :return: None
    """
    obj = cell.get_blender_object()
    bpy.data.objects.remove(obj, do_unlink=True)


def get_major_axis(obj):
    """
    Calculates the major (long) axis of a mesh by calculating the first eigenvector
    of the vertices in the mesh. This can be used to choose the axis of division

    :param obj: the blender object (mesh) -> obj = bpy.data.objects['object name']

    :return: major axis as (x, y, z) which gives direction from origin (0, 0, 0)
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
    # index 1 is used here for simulating the blast-yolk system for
    # the first 3 division
    # index -1 will be used for the purely divison simulation
    major_x, major_y, major_z = eigenvecs[:, sort_indices[1]]
    major_axis = (major_x, major_y, major_z)
    # print("this is majir axis", major_axis)
    return major_axis


def seperate_cell(obj):
    """
    Splits one cell into two for cell division. Currently divides
    along the major axis but should allow a supplied axis

    :param obj: the Blender mesh to divide in two

    :return: mother_name, daughter_name, COM, major_axis
    """
    # TODO allow division along a supplied axis.
    # If none supplied then divide along major axis

    # We need to get the cell as it is evaluated in the simulation.
    # To do this, we fetch its dependency graph and obtain the
    # evaluated cell (denoted as obj_eval here)
    bpy.ops.object.mode_set(mode='OBJECT')
    dg = bpy.context.evaluated_depsgraph_get()
    obj = obj.evaluated_get(dg)
    # The mother cell name is the name of the cell currently being divided
    mother_name = obj.name
    # By Blender convention, duplicates of existing objects are given
    # the same name as the original object, but with ".001" added to the end
    # We use the naming convention to tentatively set the name of one
    # daughter cell
    daughter_name = mother_name + ".001"
    # Get the cell's center of mass
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
    COM = obj.location
    # Go into edit mode
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    # Get the cell's major axis
    major_axis = get_major_axis(obj)
    # Get the division angles, if desired
    # phi, theta = get_division_angles(major_axis)
    # print(mother_name + " PHI: " + str(phi) + ", THETA: " + str(theta))
    # Bisect the cell along the plane specified by the center of mass (point)
    # and major axis (normal vectorr)
    bpy.ops.mesh.bisect(plane_co=COM, plane_no=major_axis, use_fill=False, flip=False)
    # Go into object mode
    bpy.ops.object.mode_set(mode='OBJECT')
    # Now we separate the split mesh into two separate objects.
    # First, obtain the vertex coordinates of the bisected mother cell.
    obj = bpy.data.objects[mother_name]
    new_vertices = obj.data.vertices
    new_vert_coords = [(obj.matrix_world @ v.co) for v in new_vertices]
    new_vert_coords = np.asarray(new_vert_coords)

    # We will choose half of the vertices to be separated into a new object.
    # We create a list of the vertex indices we will use to create the
    # daughter cell object
    separation_indices = []
    # We loop through each vertex and calculate the signed distance between the vertex
    # and the division plane (which is specified by the center of mass and major axis).
    # If the distance is greater than -0.05 (the vertex is on a specific
    # half of the cell),
    # we add the vertex's index i to our list of separation indices
    for i in range(len(new_vert_coords)):
        distance = mathutils.geometry.distance_point_to_plane(new_vert_coords[i],
                                                              COM, major_axis)
        if distance > -0.0001:
            separation_indices.append(i)
    # We go back into edit mode and deselect all vertices
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type="VERT")
    bpy.ops.mesh.select_all(action='DESELECT')
    # We return to object mode and loop through all the vertices.
    # If the vertex's index ix contained in the list of separation indices,
    # we select it for separation
    bpy.ops.object.mode_set(mode='OBJECT')
    for index in separation_indices:
        obj.data.vertices[index].select = True
    # We go back into edit mode and separate the selected vertices
    # as a new daughter cell.
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.separate(type='SELECTED')
    # We return to object mode and select only the "original" mother cell
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.data.objects[mother_name].select_set(False)
    bpy.data.objects[daughter_name].select_set(False)
    bpy.data.objects[mother_name].select_set(True)
    return mother_name, daughter_name, COM, major_axis


def calculate_volume(obj):
    """
    Calculates volume of a cell

    :param obj: the blender object (mesh) -> obj = bpy.data.objects['object name']

    :return: volume
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
    return volume


def vertices_control(obj):
    """
    Control the vertices number under EDIT mode for each daughter cell
    In this case, vertices of each cell should always be 100
    
    :param obj: The blender mesh whose vertices need to be controled

    :return: none
    
    """
    # Get vertices number of the obj
    vertices = obj.data.vertices
    vert_coords = [(obj.matrix_world @ v.co) for v in vertices]
    vert_array = np.asarray(vert_coords)
    vertices_nb = len(vert_array)

    # Get into EDIT mode, and SELCET all mesh
    # Print the vertices number before correcting
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    # print("Initial Vertices:", vertices_nb)

    # If the vertices number is not equal to 100, increase the vertices first,
    # then decrease the vertices number to the value 100
    if vertices_nb != 100:
        # Increase the vertices number, but increased vertices number cannot
        # be expected
        bpy.ops.mesh.subdivide(number_cuts=1, smoothness=0, fractal=0.1)

        # Go to the OBJECT mode, and calculate the current vertices after increasing
        bpy.ops.object.mode_set(mode='OBJECT')
        vertices = obj.data.vertices
        vert_coords = [(obj.matrix_world @ v.co) for v in vertices]
        vert_array = np.asarray(vert_coords)
        vertices_nb = len(vert_array)
        # print("First corrected number:", vertices_nb)

        # Go back to the EDIT mode, decrease the vertices number, now decreased vertices
        # can be controlled
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.decimate(ratio=100/vertices_nb)

        # DESELECT all mesh and go to the OBJECT mode, calculate the final vertices number
        # and print it out
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        vertices = obj.data.vertices
        vert_coords = [(obj.matrix_world @ v.co) for v in vertices]
        vert_array = np.asarray(vert_coords)
        vertices_nb = len(vert_array)
        # print("final corrected number:", vertices_nb)


def volume_control(obj, new_volume):
    """
    Control the volume of daughter cells after division.

    :param obj: the Blender mesh to change the volume
    :new_volume: the volume that this object will have

    :return: none
    """
    # Calculate the current volume of obj
    volume = calculate_volume(obj)

    # Go to the EDIT mode and SELECT all mesh
    # change the volume of obj from the current volume to the new volume
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.print3d_scale_to_volume(
        volume_init=volume,
        volume=new_volume)

    # DESELCT all mesh and go back to the OBJECT mode
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')


def divide(obj, cell):
    """
    Divides a cell in two along its major axis by splitting the parent mesh in two,
    translating one mesh, filling the two holes, and retriangulating the meshes

    :param obj: the Blender mesh to divide
    :param cell: the Python object containing the information of mother cell

    :return: daughter1, daughter2
    """
    # Select the obj you want to divide
    obj.select_set(True)
    # Print the theoretical mother volume based on the Python object data
    # print("theoretical mother volume: ", cell.data['init_volume'])

    # Seperate the cell, and get the mother name, daughter name, COM, and major_axis
    m_name, d_name, COM, major_axis = seperate_cell(obj)

    val = select_translate_cell(bpy.data.objects[m_name], COM, major_axis)

    # Obtain radius, shape, cell type of the mother cell
    m_radius = cell.data['radius']
    m_shape = cell.data['shape']
    m_type = cell.data['cell_type']

    # Daughter 1 is the Blender object still use the mother name
    # Select daughter 1 and repair the hole after division
    d1 = bpy.data.objects[m_name]
    d1.select_set(True)
    bpy.context.view_layer.objects.active = d1
    repair_hole(d1)

    # Rename daughter 1
    bpy.context.object.name = m_name + "0"

    # Make a cell Python object of daughter 1 to store information
    daughter1 = Cell(
        bpy.context.object.name,
        cell_type=m_type,
        radius=m_radius*pow(0.5, 1/3),
        shape=m_shape,
        loc=bpy.context.object.location)
    daughter1.data['mother'] = m_name
    daughter1.data['daughters'] = ['none', 'none']

    # Control the vertices of daughter 1
    vertices_control(d1)
    # Deselect the daughter 1
    bpy.data.objects[m_name + "0"].select_set(False)

    # Daughter 2 is the Blender object still use the d_name (mother name + 001)
    # Select daughter 2 and repair the hole after division
    d2 = bpy.data.objects[d_name]
    d2.select_set(True)
    bpy.context.view_layer.objects.active = d2
    repair_hole(d2)

    # Rename the daughter 2
    bpy.context.object.name = m_name + "1"

    # Make a cell Python object of daughter 2 to store information
    daughter2 = Cell(
        bpy.context.object.name,
        cell_type=m_type,
        radius=m_radius*pow(0.5, 1/3),
        shape=m_shape,
        loc=bpy.context.object.location)
    daughter2.data['mother'] = m_name
    daughter2.data['daughters'] = ['none', 'none']

    # Control the vertices of daughter 1 and deselect the daughter 2
    vertices_control(d2)
    bpy.data.objects[m_name + "1"].select_set(False)

    # Translate cell based on the value of 'val'

    d1_name = daughter1.data['name']
    d2_name = daughter2.data['name']
    if val == 0:
        bpy.data.objects[d2_name].select_set(False)
        bpy.data.objects[d1_name].select_set(True)
        translate_cell(bpy.data.objects[d1_name], major_axis)
        bpy.data.objects[d1_name].select_set(False)
    else:
        bpy.data.objects[d1_name].select_set(False)
        bpy.data.objects[d2_name].select_set(True)
        translate_cell(bpy.data.objects[d2_name], major_axis)
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
        bpy.data.objects[d2_name].select_set(False)
    bpy.data.objects[d2_name].select_set(False)
    bpy.data.objects[d1_name].select_set(True)

    return daughter1, daughter2


def repair_hole(obj):
    """
    Adds a new face to the cell after division (splitting) of the mesh and
    regularizes the mesh so that the mesh is evenly covered with faces

    :param obj: the Blender mesh

    :return: None
    """

    # bpy.context.view_layer.objects.active = obj
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


def translate_cell(obj, axis):
    """
    Moves cell along the given axis.
    Used during cell division so daughters have some space

    :param obj: the Blender mesh to move -> obj = bpy.data.objects['object name']
    :param axis: the axis along which to move it

    :return: None
    """

    # Calculate the magnitude of the given axis
    magnitude = ((axis[0])**2 + (axis[1])**2 + (axis[2])**2)

    # Calculate the radius of the current cell and translated based 
    # on its radius
    v = calculate_volume(obj)
    r = np.cbrt(3*v/4*np.pi)
    t = (r/2.5)*0.1
    # Normalize the axis by dividing eagh term by the magnitude.
    # Multiply each term by 0.1 so the cell is only translated a small distance
    new_vec = (t*axis[0]/magnitude, t*axis[1]/magnitude, t*axis[2]/magnitude)
    
    # translate the cell
    for i in range(3):
        if np.abs(bpy.context.object.location[i]) < 0.001:
            bpy.context.object.location[i] = 0
        bpy.context.object.location[i] = new_vec[i] + bpy.context.object.location[i]
        if np.abs(bpy.context.object.location[i]) < 0.001:
            bpy.context.object.location[i] = 0


def select_translate_cell(obj, COM, division_axis):
    """
    This function is called after dividing the mother cell.
    We compare the center of mass of the original mother cell to
    that of the corresponding daughter cell.
    We calculate the signed distance between the daughter cell's
    center of mass and the original plane of division.
    If the distance is positive, the function returns 0
    and this daughter cell will be translated later.
    If the distance is negative, the function returns 1
    and the other daughter cell will be translated later

    :param obj: the Blender mesh to divide in two
    :param COM: center of mass
    :param division_axis: division axis of cell (cleavage plane is
    perpendicular to division axis)

    :return: 0 if the daughter cell is a positive distance away from
    the plane of division and 1 otherwise
    """

    # Get the mother cell's center of mass
    com = np.copy(COM)
    # print("OLD COM: " + str(COM))
    # Set the daughter cell's location to its center of mass,
    # and set this value as the new center of mass
    # (denoted new_COM)
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
    new_COM = obj.location
    # print("NEW COM: " + str(new_COM))
    # print(com)

    # Calculate the signed distance between the new center of mass
    # and the plane of division
    distance = mathutils.geometry.distance_point_to_plane(new_COM, com, division_axis)

    # If the signed distance is greater than 0, return 0
    # print(distance)
    if distance > 0:
        return 0
    # Else, if the signed distance is less than 0, return 1
    else:
        return 1


def turn_off_physics():
    """
    Turns physics off for the currently selected Blender object.
    The cloth physics for cells are turned off before division to
    avoid irregular mesh behavior after

    :return: None
    """
    bpy.ops.object.modifier_remove(modifier="Cloth")


def turn_on_physics():
    """
    Turns physics on for the currently selected Blender object.
    Physics are turned back on after cell division occurs to avoid
    irregular mesh behavior post-division.

    :return: None
    """
    # TODO should these values be cell properties?

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


def add_material_cell(mat_name, r, g, b):
    """
    Creates a soap bubble like Blender material for use in rendering cells.
    The material has a name that allows it to be shared for multiple cells

    :param mat_name: a text name for the material
    :param r: red value [0 to 1]
    :param g: green value [0 to 1]
    :param b: blue value [0 to 1]

    :return: None
    """
    print("add cell material " + mat_name)

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
    """
    A class for representing forces between cells in Goo. Used for cell adhesion

    :param force_name: Name of force
    :param cell_name: Name of cell
    :param strength: Strength of force
    """
    def __init__(self, force_name, cell_name, strength, falloff):
        self.name = force_name
        self.strength = strength
        self.associated_cell = cell_name
        self.falloff_power = falloff

    # Member function that gets the Blender object corresponding to a Goo Force object
    def get_blender_force(self):
        obj = bpy.data.objects[self.name]
        return obj


def make_force(force):
    """
    Makes a Blender force from the Goo force

    :param force: a Goo force

    :return: None
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


def make_force_collections(master_collection, cell_types):
    """
    Make collections for forces to be stored in.

    :param master_collection: The collection that the force
    collections will be contained
    :param cell_types: list of active cell types

    :return: None
    """
    bpy.context.view_layer.active_layer_collection = master_collection
    for type in cell_types:
        collection = bpy.context.blend_data.collections.new(name=type+"_forces")
        bpy.context.collection.children.link(collection)


class handler_class:
    """
    A class for creating different types of handlers that trigger
    actions on Goo cells when certain criteria are met
    """

    # The initialization function specifies available cell types and associated
    # parameters like division rate, growth rate, and adhesion forces
    def __init__(self):
        # all the cell types
        self.cell_types = {'blast', 'yolk'}
        # all types are actived to growth, divide, generate force...
        self.active_cell_types = set()
        # python objects only for cells that aim for division 
        self.active_cells = []
        # dictionary {'cell type', value}, ex: {'blast', 60}
        self.division_rates = {}
        # dictionary {'cell type', value}, ex: {'blast', 10}
        self.growth_rates = {}
        # one of the parameters to define the force
        self.falloff = 0
        
        # dictionary of dictionaries
        # ex: {"blast": {"blast": 100, "yolk": 100},
        #      "yolk": {"blast": 200, "yolk": 200,}}
        # in the example above, blast exert 100 force on other
        # blast, 100 force on yolk
        # could change values to have different amounts of force on different cell types


        # initial parameter values to zero for each cell type
        self.adhesion_forces = {}
        for type in self.cell_types:
            self.division_rates[type] = 0
            self.growth_rates[type] = 0
            self.adhesion_forces[type] = {}
            for i in self.cell_types:
                self.adhesion_forces[type][i] = 0
        
        # collection of forces 
        self.forces = []

        return

    # Member function to add a new cell types
    def add_cell_type(self, cell_type):
        """
        Add a new cell type in the handler class
        :param cell_type: Name of cell type. 
        Must be one of the active cell types. (String)
        :return: None
        """
        self.cell_types.add(cell_type)

        # Reinitial the division, growth rate for this specif cell type
        # reset of all the adhesion_force dictionary
        self.division_rates[cell_type] = 0
        self.growth_rates[cell_type] = 0
        self.adhesion_forces[cell_type] = {}
        for type in self.cell_types:
            for i in self.cell_types:
                self.adhesion_forces[type][i] = 0

        return
    # Member function to remove a existing cell types
    def remove_cell_type(self, cell_type):
        # TODO 
        # 1. throw an error if the cell type is not in set of cell_types
        # 2. idealy, the defiend value of other cell types will not be influence 

        """
        Remove a cell type in the handler class already, if not throw
        an error. If yes, also remove this cell type in the set of 
        active_cell_types also.
        :param cell_type: Name of cell type. (String)
        :return: None
        """
        self.remove_active_cell_type(cell_type)
        self.cell_types.remove(cell_type)

        # Reinitial the division, growth rate for rest cell types
        # reset of all the adhesion_force dictionary
        self.division_rates[cell_type] = 0
        self.growth_rates[cell_type] = 0
        self.adhesion_forces[cell_type] = {}
        for type in self.cell_types:
            for i in self.cell_types:
                self.adhesion_forces[type][i] = 0

        return

    # Member function to active a specific cell type that already exist
    def add_active_cell_type(self, active_cell_type):
        # TODO throws an error if the "active_cell_type" is not in the
        # set of cell_types
        """
        Active a cell type, make it into the set of active_cell_types

        :param active_cell_type: One of cell types that exists in the
        set of cell_types (String) 

        :return: None
        """
        self.active_cell_types.add(active_cell_type)
        return
    
    # Member function to active a specific cell type that already exist
    # in the set of active_cell_types
    def remove_active_cell_type(self, active_cell_type):
        # TODO throws an error if the "active_cell_type" is not in the
        # set of active_cell_types
        """
        Deactive a cell type, delete it in the set of active_cell_types
        
        :param active_cell_type: One of cell types that exists in the
        set of active_cell_types (String) 
        :return: None
        """

        self.active_cell_types.remove(active_cell_type)
        return

    # Member function to add a python object(cell) into the active_cells,
    # which is a list that stores the python object of division mesh 
    def add_cell(self, cell):
        """
        Add the pyhton object of a cell into the handler.active_cells,
        which allows us to get the corresponding python object of 
        every mesh for better volume control and growth control
        Also cell type of this object will be automatically added into
        set of cell_types and set of active_cell_types

        :param cell: python object of a cell, which belongs to
        Cell class(Cell)

        :return: None
        """
        # if the cell type of this cell is not in the set of cell_type
        # it's cell type will be automatically added into this set
        if cell.data['cell_type'] not in self.cell_types:
            self.add_cell_type(cell.data['cell_type'])
        # same for the set of active_cell_type
        if cell.data['cell_type'] not in self.active_cell_types:
            self.add_active_cell_type(cell.data['cell_type'])
        # add this cell into list of active_cells
        self.active_cells.append(cell)
        return

    # Member function to remove a python object(cell) into the handler
    def remove_cell(self, cell):
        """
        Remove the pyhton object of a cell into the handler.active_cells,
        :param cell: python object of a cell, which belongs to
        Cell class(Cell)

        :return: None
        """
        self.active_cells.remove(cell)
        return

     # Member function to set growth rate for a cell type
    def set_growth_rate(self, cell_type, rate):
        """
        Sets rate rate in the handler class that growth_handler() can reference later

        :param cell_type: Name of cell type to apply this division rate to.
        Must be one of the active cell types. (String)
        :param rate: For every given frame number, tge cell will growth once 

        :return: None
        """
        self.growth_rates[cell_type] = rate
        return
    
    # Member function to handle cell growth
    def growth_handler(self, scene, depsgraph):
        # TODO this function can apply for all cell types
        """
        Increase the volume of a cell type that will experience division

        :return: None
        """
        current_frame = scene.frame_current

        for cell in self.active_cells:
            cell_name = cell.data['name']
            # get the defined growth rate
            growth_rate = self.growth_rates[cell.data['cell_type']]
            # get the defined division rate
            division_rate = self.division_rates[cell.data['cell_type']]
            if current_frame % division_rate != 0:
                if current_frame % division_rate % growth_rate == 0:
                    cell = bpy.data.objects[cell_name]
                    cell.select_set(True)
                    # print("----------")
                    # print("current_frame: ", current_frame)
                    # print("cell: ", cell_name)
                    # print("old volume: ", calculate_volume(cell))

                    # get the pyhton object of a cell
                    for c in self.active_cells:
                        if c.data['name'] == cell_name :
                            # find the corresponding mesh in blender
                            bpy.context.view_layer.objects.active = bpy.data.objects[cell_name]
                            bpy.ops.object.select = True
                            turn_off_physics()

                            # change the volume of the mesh based on the calculation by python object
                            volume_control(cell, c.data['init_volume'] * 1.01)
                            c.data['init_volume'] = c.data['init_volume'] * 1.01
                            c.data['radius'] = \
                                pow(c.data['init_volume'] * 3 / 4 / np.pi, 1/3)
                            # print("new volume: ", calculate_volume(cell))
                            turn_on_physics()
                            bpy.context.active_object.select_set(False)
                    cell.select_set(False)

    # Member function to define the division rate for a specif cell type
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

    # Member function to define the fallof parameter for the force objects
    # which will be generated during the division
    def set_falloff(self,value):
        """
        Sets division rate in the handler class that div_handler() can reference later
        :param value: value for the falloff strength for the force (int)

        :return: None
        """
        self.falloff = value
        return

    # Member function to set the division handler for cells and generate
    # corresponding forces between different cell types if they have been defined
    def div_handler(self, scene, depsgraph):
        current_frame = scene.frame_current
        # collect the new generated cell objects during the division
        new_cells = []

        # Loop through active cell types
        for cell_type in self.active_cell_types:
            # Get the number of cells of that type
            num_cells = len(bpy.data.collections[cell_type].objects)
            # Divide based on the division rate
            division_rate = self.division_rates[cell_type]
            if division_rate == 0:
                continue
            if current_frame % division_rate == 0:
                # print("-------DIVIDING CELLS------")
                for i in range(num_cells):
                    # Get the cell name
                    cell_name = bpy.data.collections[cell_type].objects[i].name

                    # Get the corresponding Blender object
                    Blender_cell = bpy.data.objects[cell_name]

                    # Select the Blender object and make it the active object
                    bpy.data.objects[cell_name].select_set(True)
                    bpy.context.view_layer.objects.active = bpy.data.objects[cell_name]
                    # print(cell_name)

                    # Turn off the cloth physics for this cell. This mitigates
                    # irregular mesh behavior after division
                    turn_off_physics()

                    for cell in self.active_cells:
                        if cell_name == cell.data['name']:
                            global d1, d2
                            d1, d2 = divide(Blender_cell, cell)
                            self.remove_cell(cell)
                            new_cells.append(d1)
                            new_cells.append(d2)

                    blender_d1 = bpy.data.objects[d1.data["name"]]
                    blender_d2 = bpy.data.objects[d2.data["name"]]

                    # Select the first daughter cell only
                    vl = bpy.context.view_layer
                    blender_d1.select_set(True)
                    blender_d2.select_set(False)
                    vl.objects.active = blender_d1
                    
                    turn_on_physics()

                    # Select the second daughter cell only
                    blender_d1.select_set(False)
                    blender_d2.select_set(True)
                    vl.objects.active = blender_d2

                    # Turn on the physics for this daughter cell
                    turn_on_physics()

                    blender_d2.select_set(False)
                    # print("active cells number: ", len(self.active_cells))

                # print("----------Force---------")
                self.apply_forces(self.falloff)

        # if current_frame % division_rate == 0:
            if len(bpy.data.collections[cell_type].objects) > num_cells:
                for cell in new_cells:
                    self.add_cell(cell)

            # print("Now cells in active_cells: ")
        # for cell in self.active_cells:
                # print(cell.data['name'], ", volume: ", cell.data['init_volume'])    
        return

    # Member function to set adhesion forces between cell types
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

    # Member function to generate forces between different cell types
    # if they have been redefined
    def apply_forces(self, falloff):
        """
        Add force fields to force collections and make them affect corresponding cell
        types

        :param falloff: Set the falloff value for all forces are going to be generated
        :return: None
        """
        master_collection = bpy.context.view_layer.active_layer_collection
        # Find the number of all cells in set of active cell type
        for cell_type in self.active_cell_types:
            num_cells = len(bpy.data.collections[cell_type].objects)
            # loop all cells
            for i in range(num_cells):
                # set the effector collection parameter in the "CLOTH" modifier
                # corresponding to its cell type
                cell_name = bpy.data.collections[cell_type].objects[i].name
                cell = bpy.data.objects[cell_name]
                cell.modifiers["Cloth"].settings.effector_weights.collection = \
                    bpy.data.collections[cell_type+"_forces"]

                # generated forces between different cell types based on the adhesion_force    
                for affected_type in self.active_cell_types:
                    if self.adhesion_forces[cell_type][affected_type] != 0:
                        vl = bpy.context.view_layer
                        affected = vl.layer_collection.children[affected_type+"_forces"]
                        bpy.context.view_layer.active_layer_collection = affected
                        f = Force(cell_name+'_to_'+affected_type,
                                  cell_name,
                                  self.adhesion_forces[cell_type][affected_type], falloff)
                        make_force(f)
                        vl.active_layer_collection = master_collection
                        # add all generated force objects into the forces list
                        self.forces.append(f)
                    
                        
    # Member function to remove forces that belong to mother cell(which will 
    # dispear after new daughter cells are produced) and also make the each 
    # force moves as the corresponding cell moves
    def adhesion_handler(self, scene, depsgraph):
        """
        A handler that focuses on force number and force loction.
        It will make the number of force the same as the number of cells,
        and moves the force as the corresponding cell moves

        """
        current_cells_names = []
        current_forces = []

        # get all force objects store in the list of forces
        for f in self.forces:
            current_forces.append(f)
        # get all the cell (actice cell type) names and store in the list
        for cell_type in self.active_cell_types:
            num_cells = len(bpy.data.collections[cell_type].objects)
            for i in range(num_cells):
                cell_name = bpy.data.collections[cell_type].objects[i].name
                current_cells_names.append(cell_name)
                # print(cell_name)

        # print("How many forces: ", len(self.forces))

        # if the corresponding cell of a force does not exist in the colletion
        # this force will be removed, others stay the same 
        for force in self.forces:
            if force.associated_cell not in current_cells_names:
                # print(force.name, ' is removed')
                bpy.data.objects.remove(
                                 bpy.data.objects[force.get_blender_force().name],
                                 do_unlink=True)
                current_forces.remove(force)
            else:
                # print(force.name, ' is not removed')

                # calculate the location of each cell
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
                # set the location of the force exactly the same as the cell that
                # force attaches to
                bpy.data.objects[force.name].location = COM
        
        # after deletion, put exsiting forces into the list of forces
        self.forces = current_forces