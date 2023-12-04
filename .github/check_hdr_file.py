import os

file_path = '/home/runner/blender/3.0/scripts/modules/\
    goo/missile_launch_facility_01_4k.hdr'

if os.path.exists(file_path):
    print("File exists")
else:
    print("File does not exist")