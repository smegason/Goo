import bpy
import sys
import os

# Set the path to the directory containing the goo module
directory_path = "/home/runner/work/Goo/Goo/scripts/modules"

# Add the directory to sys.path to enable module import
if directory_path not in sys.path:
    sys.path.append(directory_path)

# Import the goo module and access its functions
try:
    from goo import goo  # Replace 'goo_function' with the actual function name you want to use
except ImportError:
    print("Error: Unable to import 'goo' module")

# Test using the imported function
try:
    goo.setup_world()  # Replace 'goo_function' with the actual function you want to use
except Exception as e:
    print("Error:", e)
