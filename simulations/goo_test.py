# goo.py - This is the Goo library. It contains all helper functions for goo
# goo is licensed under BSDv2

import sys
import bpy
import numpy as np

sys.path.append('C:\\Users\\anr9744\\Projects\\Goo\\scripts\\modules\\goo')
import simulations.cell as cell
import simulations.force as force


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


'''def mitosis_handler(scene):
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
    cell_tree.show()'''


'''# Remove this function. Take try/except code to make_cell function
def make_mesh(cell):
    """
    Makes a Blender mesh object for the given Goo cell

    :param cell: the Goo cell to make a Blender mesh from

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

'''def make_cell(cell):
    """
    Makes a Blender mesh object for the given Goo cell Python object

    :param cell: the Goo cell Python object to make a Blender mesh from

    :return: None
    """
    # TODO make_mesh and make_cell are redundant. Need to merge.

    # Making cell
    print('Making cell')
    print(cell.data['flavor'] == "")
    # if cell.data['flavor'] == "round_cube" or None:
    # Add a round_cube mesh
    if cell.data['flavor'] == "round_cube" or "":
        print("Making Round Cube")
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
            print("To enable RoundCube creation for Cells you must go to ")
            print("Edit->Preferences->AddOns->Add Mesh:ExtraObjects and ")
            print("check the box to enable it")
            return           
            
    # Ico sphere
    # elif cell.data['flavor'] == "ico_sphere":
    else:
        print("Making Ico Sphere")
        try:
            bpy.ops.mesh.primitive_ico_sphere_add(enter_editmode=False,
                                                  align='WORLD',
                                                  location=cell.data['location'],
                                                  scale=cell.data['scale'],
                                                  radius=cell.data['radius'])

        except Exception:
            print(sys.exc_info())
            print("Make sure you spell the correct name")
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
    material = bpy.data.materials.get(cell.data['material'])
    if (material):
        bpy.context.active_object.data.materials.append(material)
    else:
        print("No material ", cell.data['material'])
'''

def delete_cell(cell):
    """
    Deletes a Goo cell's Blender object

    :param cell: the Goo cell to delete

    :return: None
    """
    obj = cell.get_blender_object()
    bpy.data.objects.remove(obj, do_unlink=True)


'''def add_material_cell(mat_name, r, g, b):
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
    links.new(node_mix.outputs[0], node_output.inputs[0])  # link_mix_out'''

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
    '''
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
    '''
    c = cell.Cell("cell", loc=(0, 0, 0))
    cell.make_cell(c)
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
    c = cell.Cell("cell", loc=(0, 0, 0))
    cell.make_cell(c)
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
    c1 = cell.Cell("cell1", loc=(0, 0, 0))
    cell.make_cell(c1)
    c2 = cell.Cell("cell2", loc=(1, 1, -1))
    cell.make_cell(c2)
    c3 = cell.Cell("cell3", loc=(-1, 1, 1))
    cell.make_cell(c3)
    c4 = cell.Cell("cell4", loc=(-1, -1, -1))
    cell.make_cell(c4)
    c5 = cell.Cell("cell5", loc=(1, -1, 1))
    cell.make_cell(c5)
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


def setup_world():
    """
    Sets up the default values used for simulations in Goo
    including units and rendering background

    :return: None
    """
    print("setup world")
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

# TODO: document and see if it works 
# idea is to clear the object hierarchy before running a simulation    
def delete_hierarchy(parent_obj_name):
        bpy.ops.object.select_all(action='DESELECT')
        obj = bpy.data.objects[parent_obj_name]
        obj.animation_data_clear()
        names = set()
        # Go over all the objects in the hierarchy like @zeffi suggested:
        def get_child_names(obj):
            for child in obj.children:
                names.add(child.name)
                if child.children:
                    get_child_names(child)

        get_child_names(obj)
        names.add(parent_obj_name)
        print(names)
        objects = bpy.data.objects
        
        # Remove the animation from the all the child objects
        if names:
            for child_name in names:
                bpy.data.objects[child_name].animation_data_clear()
                objects[child_name].select_set(state=True)
                bpy.data.objects.remove(objects[child_name])
            print ("Successfully deleted object")
        else:
            print ("Could not delete object")


class handler_class:
    # TODO document this class
    """
    A class for creating different types of handlers that trigger
    actions on Goo cells when certain criteria are met
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
                        f = force.Force(cell_name+"_to_"+affected_type,
                                  cell_name,
                                  self.adhesion_forces[cell_type][affected_type])
                        force.make_force(f)
                        vl.active_layer_collection = master_collection
                        self.forces.append(f)

    # Member function to set the division handler for cells - not used 
    '''def div_handler(self, scene, depsgraph):
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
                    bpy.data.objects[d2.data["name"]].select_set(False)'''

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
