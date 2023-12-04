import sys

# Set the path to the directory containing the goo module
directory_path = "/home/runner/work/Goo/Goo/scripts/modules"

# Add the directory to sys.path to enable module import
if directory_path not in sys.path:
    sys.path.append(directory_path)

# Import the goo module and access its functions
try:
    from goo import goo 
except ImportError:
    print("Error: Unable to import 'goo' module")

# Test using the imported function
try:
    goo.setup_world(seed=1)
    '''goo.make_cell(
        name='cell_A1',
        loc=(0, 0, 0), 
        type='cellsA'        
    )

    print(bpy.data.objects['cellA1'].users_collection[0])'''
    
except Exception as e:
    print("Error:", e)
