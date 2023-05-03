# benchmark_two_cells.py - simulates cell doublets

from goo import goo
import bpy
import argparse
import sys

parser = argparse.ArgumentParser()
#parser.add_argument('--falloff', type=float, required=True)
parser.add_argument('--adhesion', type=float, required=True)
parser.add_argument('--tension', type=float, required=True)
args, unknown = parser.parse_known_args(sys.argv[sys.argv.index("--") + 1:])

goo.setup_world()


#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection("A_Cells")
# Define cell A1
goo.make_cell("cell_A1", loc = (3,0,0), stiffness = float(args.tension), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A2", loc = (5,0,0), stiffness = float(args.tension), collection = "A_Cells")


#================== Force A Collection ==================
# Force parameters; strength and decay power
force_strength = args.adhesion
force_falloff = 1

# Create a collection for force A
goo.make_collection("A_Forces")
# Define force A1
goo.make_force("force_A1", "cell_A1", force_strength, force_falloff, "A_Forces")
# Define force A2
goo.make_force("force_A2", "cell_A2", force_strength, force_falloff, "A_Forces")

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, 
                           end = 20, 
                           filepath = "C:\\Users\\anr9744\\Projects\\Goo\\data\\TEST11_doublets", 
                           adhesion = True, 
                           data = False)
