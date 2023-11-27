
import os 

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

# Set the path to the directory
directory_path = "/home/runner/work/Goo/Goo/scripts/modules/goo"

# Check if the directory exists
if os.path.exists(directory_path):
    # Perform an ls in the directory
    try:
        files = os.listdir(directory_path)
        print("Files in the directory:")
        for file in files:
            print(file)
    except Exception as e:
        print(f"Error: {e}")
else:
    print("Error: Directory does not exist")
