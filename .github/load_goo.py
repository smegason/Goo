try:
    import bpy
except ImportError:
    print("Error: Unable to import 'bpy'")

try:
    import sys
except ImportError:
    print("Error: Unable to import 'sys'")

try:
    sys.path.append("/home/runner/work/Goo/Goo/scripts/modules/goo/goo.py")
    from goo import goo
except ImportError:
    print("Error: Unable to import 'goo'")
