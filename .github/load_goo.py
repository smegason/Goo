import bpy


# Set the path to your scripts folder
scripts_path = "/path/to/your/scripts"

# Access the preferences
prefs = bpy.context.preferences

# Set the scripts path
prefs.filepaths.script_directory = scripts_path
