#test_goo.py
#this file is read by pytest to test all functions in Goo  xxx

# importing module
#import sys
#sys.path.append('scripts/modules')

# importing required module
#import goo
import bpy

import os
os.chdir("scripts/modules/goo")
print (os.listdir())   # For some weird reason goo.py is not in this directory but other .py files are e.g. goo2.py

import importlib.util
spec = importlib.util.spec_from_file_location("goo", "goo2.py")
goo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(goo)
#print(foo.var)

def test_sample_func():
    print("test sample function------------")

    version = bpy.app.version_string
    print (version)
 
    print ("AAAAAA")
    #print (dir(goo))
    print ("BBBBBBB")
    ret = goo.sample_func()
    print (ret)
    print("Test sample func..........")

    assert True  


def test_calculate_volume():
    #obj = ???? need to get Blender object to test with somehow
    #volume = calculate_volume(obj)
    #assert volume = ??
    print ("test calculate volume")

    assert True

    
