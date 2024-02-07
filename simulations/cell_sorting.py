from goo import goo
from importlib import reload
import bpy
reload(goo)
goo.setup_world(seed=1)

# Cell A
goo.make_cell("cell_A1", loc=(0, 0, 0), type="cellsA")
goo.make_cell("cell_A2", loc=(0, 3, -2), type="cellsA")
goo.make_cell("cell_A3", loc=(-1, 5, 0), type="cellsA")
goo.make_cell("cell_A4", loc=(-1, -5, 0), type="cellsA")
goo.make_cell("cell_A5", loc=(-1, 5, -3), type="cellsA")
goo.make_cell("cell_A6", loc=(-5, 0, 0), type="cellsA")
goo.make_cell("cell_A7", loc=(2, 4, 0), type="cellsA")
goo.make_cell("cell_A8", loc=(4, 0, 2), type="cellsA")
goo.make_cell("cell_A9", loc=(-1, -5, 0), type="cellsA")
goo.make_cell("cell_A10", loc=(1, -2, 3), type="cellsA")

# Force A
homoA = 2000
goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_homo_adhesion('cell_A2', -homoA)
goo.add_homo_adhesion('cell_A3', -homoA)
goo.add_homo_adhesion('cell_A4', -homoA)
goo.add_homo_adhesion('cell_A5', -homoA)
goo.add_homo_adhesion('cell_A6', -homoA)
goo.add_homo_adhesion('cell_A7', -homoA)
goo.add_homo_adhesion('cell_A8', -homoA)
goo.add_homo_adhesion('cell_A9', -homoA)
goo.add_homo_adhesion('cell_A10', -homoA)
# goo.add_motion('cell_A1', -motion)

bpy.ops.object.effector_add(type='TURBULENCE',
                            enter_editmode=False, 
                            align='WORLD', 
                            location=(0, 0, 0), 
                            scale=(1, 1, 1))
name = 'random motion'
bpy.context.object.name = name
bpy.context.object.field.shape = 'POINT'
bpy.context.object.field.strength = -20000.0
bpy.context.object.field.size = 6
bpy.context.object.field.noise = 0
bpy.context.object.field.seed = 92
bpy.context.object.field.apply_to_location = True
bpy.context.object.field.apply_to_rotation = True
bpy.context.object.field.use_min_distance = True
bpy.context.object.field.distance_min = 0
bpy.context.object.field.use_max_distance = True
bpy.context.object.field.distance_max = 20
bpy.context.object.field.flow = 0
bpy.context.object.field.wind_factor = 0
bpy.context.scene.collection.objects.link(bpy.data.objects[name])

goo.add_sphere_boundaries(loc=(0, 0, 0), radius=10)

# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=10000,  # default, 250
                           filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/division/20240102_division_tests/test01", 
                           adhesion=True,  # default, True
                           data=True,  # default, False
                           growth=True, 
                           motility=False, 
                           division=False, 
                           target_volume=50,
                           growth_rate=1, 
                           growth_type='linear'
                           )
