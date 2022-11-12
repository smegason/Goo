#script to simulate yolk-cell interface
import sys
from goo import goo
import bpy
import numpy as np

goo.setup_world()
goo.add_material_cell("CellGreen", 0.007, 0.300, 0.005)
goo.add_material_cell("CellBlue", 0.007, 0.021, 0.300)

cell = goo.Cell(name_string="Yolk", loc=(0, 0, 0), material="CellGreen")
goo.make_cell(cell)

##is the 'Cell' class from original goo => 'cell' to 'blast' 
##IDEA: create method to adjust the size of the of each cell, instead of creating a new class blast (How??)
    ## add extra scale and radius parameter in the arguments? (similar to loc, material, flavour)
class Blast():
    def __init__(self, name_string, loc, material="", flavor=""):
        # The initialization function sets a cell data dictionary
        # for geometric parameters, physics parameters,
        # division information, and cell lineage
        self.data = {
            'ID': 0,
            'name': name_string,
##adapted radius: smaller radius
            'radius': 1.5,
            'enter_editmode': False,
            'align': 'WORLD',
            'location': loc,
            'material': material,
            'size': (1, 1, 1),
##adapted scale to flatten the sphere in the Z-direction
            'scale': (1, 1, 0.7),
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
            'flavor': flavor
        }
        # The volume and mass are calculated from values in the data dictionary
        self.data['init_volume'] = ((4/3)*np.pi*(self.data['radius'])**3)
        self.data['mass'] = self.data['density']*self.data['init_volume']

def make_blast(Blast):
    if blast.data['flavor'] == "round_cube" or "":
            print("Making Round Cube")
            try:
                bpy.ops.mesh.primitive_round_cube_add(change=False,
                                                      radius=blast.data['radius'],
                                                      size=blast.data['size'],
                                                      arc_div=blast.data['arcdiv'],
                                                      lin_div=0,
                                                      div_type='CORNERS',
                                                      odd_axis_align=False,
                                                      no_limit=False,
                                                      location=blast.data['location'])
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
                                                    location=blast.data['location'],
                                                    scale=blast.data['scale'],
                                                    radius=blast.data['radius'])

        except Exception:
            print(sys.exc_info())
            print("Make sure you spell the correct name")
            return

    # Give the Blender object the cell's name
    obj = bpy.context.object
    bpy.context.object.name = blast.data['name']
    bpy.context.view_layer.objects.active = bpy.data.objects[blast.data['name']]

    # Smooth the mesh
    bpy.ops.object.select = True
    bpy.ops.object.shade_smooth()

    # Add subsurface modifier to make smoother
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = blast.data['subdiv']

    # Add cloth settings for physics
    bpy.ops.object.modifier_add(type='CLOTH')
    bpy.context.object.modifiers["Cloth"].settings.bending_model = 'LINEAR'
    bpy.context.object.modifiers["Cloth"].settings.quality = 5
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1
    bpy.context.object.modifiers["Cloth"].settings.mass = blast.data['vertex_mass']
#    obj.modifiers["Cloth"].settings.air_damping = blast.data['air_damping']
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
    material = bpy.data.materials.get(blast.data['material'])
    if (material):
        bpy.context.active_object.data.materials.append(material)
    else:
        print("No material ", blast.data['material'])
        
        
##Location sets the blast on top of the yolk        
blast = Blast(name_string="Blast", loc=(0, 0, 0.6), material="CellBlue", flavor = "")
make_blast(blast)

bpy.data.objects['Blast'].select_set(True)
bpy.ops.object.modifier_add(type='BOOLEAN')
bpy.context.object.modifiers["Boolean"].operation = 'DIFFERENCE'
bpy.context.object.modifiers["Boolean"].object = bpy.data.objects['Yolk']
bpy.context.object.modifiers["Boolean"].solver = 'FAST'
bpy.ops.object.modifier_apply(modifier="Boolean")