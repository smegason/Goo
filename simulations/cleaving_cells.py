# benchmark_two_cells.py - simulates cell doublets

from goo import goo
import bpy

goo.setup_world()

#====================== Colors =========================
goo.add_material("green", 0, 0.1, 0)
goo.add_material("red", 0.1, 0, 0)


#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection("A_Cells")
# Define cell A1

goo.make_cell("cell_A1", loc = (0,0,0), material = "green", collection = "A_Cells")


#================== Simulation setup ==================
#goo.simulation_stiffness(1, 1, 1, 1)

handlers = goo.handler_class() 
handlers.frame_interval = [1,100]

handlers.active_cell_types += ['A_Cells']
handlers.set_division_rate('A_Cells', 20)
handlers.set_growth_rate('A_Cells', 20)

# Set the current frame to the start frame
bpy.context.scene.frame_set(handlers.frame_interval[0])
# Set the animation to play from start to end frame
bpy.context.scene.frame_start = handlers.frame_interval[0]
bpy.context.scene.frame_end = handlers.frame_interval[1]

handlers.data_file_path = "C:\\Users\\anr9744\\Projects\\Goo\\data\\doublet_test_axis"
#bpy.app.handlers.frame_change_pre.append(handlers.timing_init_handler)

bpy.app.handlers.frame_change_post.clear()

bpy.app.handlers.frame_change_post.append(handlers.div_handler)
bpy.app.handlers.frame_change_post.append(handlers.growth_handler)


#bpy.app.handlers.frame_change_post.append(handlers.adhesion_handler)
#bpy.app.handlers.frame_change_pre.append(handlers.data_handler)

#bpy.app.handlers.frame_change_post.append(handlers.timing_elapsed_handler)
bpy.app.handlers.frame_change_post.append(handlers.stop_animation)




bpy.ops.screen.animation_play()
