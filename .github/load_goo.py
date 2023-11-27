import bpy
print('bpy imported')

# Set the path to your scripts folder
scripts_path = "/home/runner/work/Goo/Goo/scripts/"
print(scripts_path)

# Access the preferences
prefs = bpy.context.preferences

# Set the scripts path
prefs.filepaths.script_directory = scripts_path

from goo import goo
