#script to simulate first four divisions up to 16 cells 

# Import Libraries and setup world
from goo import goo_kuleuven as modifiedGoo
import importlib
import bpy

importlib.reload(modifiedGoo)

modifiedGoo.setup_world()

#setup material

modifiedGoo.add_material_cell("CellGreen", 0.007, 0.300, 0.005)
modifiedGoo.add_material_cell("purple", 0.1, 0, 0.8)

# create a python object named 'blast'
blast=modifiedGoo.Cell(name_string="blast",
                       cell_type="blast",
                       loc=(0,0,1.2),
                       radius=0.63,
                       shape='sphere',
                       material="CellGreen")
                       
# make it as a mesh in blender
modifiedGoo.make_cell(blast)

# create a python object named 'yolk'
yolk = modifiedGoo.Cell(name_string="yolk",
                       cell_type="yolk",
                       loc=(0, 0, 0),
                       radius=1,
                       shape='sphere',
                       material="purple")
# make it as a mesh 
modifiedGoo.make_cell(yolk)

# create a collection to all store the forces that only blast will be affected
f_collection = bpy.data.collections.new("blast_forces")
# link the collection to the scene
bpy.context.scene.collection.children.link(f_collection)
bpy.context.view_layer.active_layer_collection = \
            bpy.context.view_layer.layer_collection.children["blast_forces"]
            
# create a force that interact with blast to blast     
blast_force=modifiedGoo.Force('blast_to_blast',
                                  'blast',
                                  -100, 1)
# make it as a blender object
modifiedGoo.make_force(blast_force)

# create a force that interact with yolk to blast   
yolk_blast_force=modifiedGoo.Force('yolk_to_blast',
                                  'yolk',
                                  -800, 1)
# make it as a blender object
modifiedGoo.make_force(yolk_blast_force)

# create a collection to all store the forces that only yolk will be affected
f_collection = bpy.data.collections.new("yolk_forces")
# link the collection to the scene
bpy.context.scene.collection.children.link(f_collection)
bpy.context.view_layer.active_layer_collection = \
            bpy.context.view_layer.layer_collection.children["yolk_forces"]

# create a force that interact with blast to yolk  
yolk_force=modifiedGoo.Force('blast_to_yolk',
                                  'blast',
                                  -100, 1)
# make it as a blender object                               
modifiedGoo.make_force(yolk_force)

# create a handler class
handlers = modifiedGoo.handler_class()
# add the 'blast' python object into handler
handlers.add_cell(blast)
# make the 'yolk' cell type activated
handlers.add_active_cell_type('yolk')
# add the blast_force into the list of forces
handlers.forces.append(blast_force)
# add the yolk_force into the list of forces
handlers.forces.append(yolk_force)

# set falloff, division rate, growth rate, and adhesion rate 
handlers.set_falloff(1)
handlers.set_division_rate('blast', 30)
handlers.set_division_rate('yolk', 0)
handlers.set_growth_rate('blast', 30)
handlers.set_adhesion('blast', 'blast', -100)
handlers.set_adhesion('blast', 'yolk', -100)

# Add handlers for division, growth and adhesion
bpy.app.handlers.frame_change_post.append(handlers.div_handler)
bpy.app.handlers.frame_change_post.append(handlers.growth_handler)
bpy.app.handlers.frame_change_post.append(handlers.adhesion_handler)

# make the 'yolk' object only will be influenced by the forces under the "yolk_ forces" collection
bpy.context.view_layer.objects.active = bpy.data.objects['yolk']
bpy.ops.object.select = True
bpy.context.object.modifiers["Cloth"].settings.effector_weights.collection = bpy.data.collections["yolk_forces"]
bpy.context.active_object.select_set(False)

# close and reopen the "CLOTH" modifier for the 'blast' to make blast and yolk start to interacti
# which is weired but we do have any other solution to make them interact, which they should be 
bpy.context.view_layer.objects.active = bpy.data.objects['blast']
bpy.ops.object.select = True
modifiedGoo.turn_off_physics()
modifiedGoo.turn_on_physics()
bpy.context.active_object.select_set(False)

bpy.context.view_layer.objects.active = bpy.data.objects['blast']
bpy.ops.object.select = True
bpy.context.object.modifiers["Cloth"].settings.effector_weights.collection = bpy.data.collections["blast_forces"]
bpy.context.object.modifiers["Cloth"].settings.bending_model = 'ANGULAR'
bpy.context.object.modifiers["Cloth"].settings.quality = 5
bpy.context.object.modifiers["Cloth"].settings.time_scale = 1
bpy.context.object.modifiers["Cloth"].settings.mass = 0.3
bpy.context.object.modifiers["Cloth"].settings.air_damping = 5

bpy.context.object.modifiers["Cloth"].settings.tension_stiffness = 1
bpy.context.object.modifiers["Cloth"].settings.compression_stiffness = 1
bpy.context.object.modifiers["Cloth"].settings.shear_stiffness = 1
bpy.context.object.modifiers["Cloth"].settings.bending_stiffness = 1

bpy.context.object.modifiers["Cloth"].settings.tension_damping = 25
bpy.context.object.modifiers["Cloth"].settings.compression_damping = 25
bpy.context.object.modifiers["Cloth"].settings.shear_damping = 25
bpy.context.object.modifiers["Cloth"].settings.bending_damping = 0.5

# Cloth > Pressure
bpy.context.object.modifiers['Cloth'].settings.uniform_pressure_force = 5
bpy.context.object.modifiers['Cloth'].settings.use_pressure_volume = True
bpy.context.object.modifiers['Cloth'].settings.target_volume = 1
bpy.context.object.modifiers['Cloth'].settings.pressure_factor = 1
bpy.context.object.modifiers['Cloth'].settings.fluid_density = 0

# Cloth > Collisions
bpy.context.object.modifiers['Cloth'].collision_settings.collision_quality = 6
bpy.context.object.modifiers['Cloth'].collision_settings.use_collision = True
bpy.context.object.modifiers['Cloth'].collision_settings.use_self_collision = True
bpy.context.object.modifiers['Cloth'].collision_settings.self_friction = 5
bpy.context.object.modifiers['Cloth'].collision_settings.self_distance_min = 0.015
bpy.context.object.modifiers['Cloth'].collision_settings.self_impulse_clamp = 0

# Collision
bpy.context.object.modifiers['Collision'].settings.damping = 0.579821
bpy.context.object.modifiers['Collision'].settings.thickness_outer = 0.02
bpy.context.object.modifiers['Collision'].settings.thickness_inner = 0.2
bpy.context.object.modifiers['Collision'].settings.cloth_friction = 5
bpy.context.object.modifiers['Collision'].settings.use_culling = True

# Cloth > Internal Springs
bpy.context.object.modifiers['Cloth'].settings.use_internal_springs = True
bpy.context.object.modifiers['Cloth'].settings.internal_spring_max_length = 0
bpy.context.object.modifiers['Cloth'].settings.internal_spring_max_diversion = 0.785398
bpy.context.object.modifiers['Cloth'].settings.internal_spring_normal_check = True
bpy.context.object.modifiers['Cloth'].settings.internal_tension_stiffness = 1
bpy.context.object.modifiers['Cloth'].settings.internal_compression_stiffness = 1
bpy.context.object.modifiers['Cloth'].settings.internal_tension_stiffness_max = 5
bpy.context.object.modifiers['Cloth'].settings.internal_compression_stiffness_max = 5

bpy.context.active_object.select_set(False)
