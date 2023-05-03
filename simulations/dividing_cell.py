# Import Libraries and setup world
from goo import goo
import importlib
import bpy

importlib.reload(goo)

goo.setup_world()

#====================== Colors =========================
goo.add_material("green", 0, 0.1, 0)
goo.add_material("red", 0.1, 0, 0)


#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection("A_Cells")
# Define cell A1

goo.make_cell("cell_A1", loc = (0,0,0), material = "green", collection = "A_Cells")

#================== Force A Collection ==================
# Force parameters; strength and decay power
force_strength = 1
force_falloff = 1

#float(sys.argv[-1])
# Create a collection for force A
goo.make_collection("A_Forces")
# Define force A1
fA1 = goo.make_force("force_A1", "cell_A1", force_strength, force_falloff, "A_Forces")

#================== Simulation setup ==================
goo.simulation_stiffness(1, 1, 1, 1)
            
# Add handlers for division, growth and adhesion
handlers = goo.handler_class()
handlers.active_cell_types += ['A_Cells']
handlers.set_division_rate('A_Cells', 50)
handlers.set_growth_rate('A_Cells', 50)
handlers.forces = [fA1]
bpy.app.handlers.frame_change_post.append(handlers.div_handler)
bpy.app.handlers.frame_change_post.append(handlers.adhesion_division_handler)
bpy.app.handlers.frame_change_post.append(handlers.growth_handler)
