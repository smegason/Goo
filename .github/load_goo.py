import bpy

# Check if bpy can be imported
try:
    import bpy
except ImportError:
    print("Error: Unable to import 'bpy'")

# Check if sys can be imported
try:
    import sys
except ImportError:
    print("Error: Unable to import 'sys'")

# Check if goo can be imported and set the script path
try:    
    # Access the preferences
    prefs = bpy.context.preferences
    # Set the scripts path
    scripts_path = "/home/runner/work/Goo/Goo/scripts/"
    prefs.filepaths.script_directory = scripts_path

    from goo import goo
    
except ImportError:
    print("Error: Unable to import 'goo'")
