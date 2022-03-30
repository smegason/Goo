#test_goo.py
#this file is read by pytest to test all functions in Goo  xxx

# importing module
#import sys
#sys.path.append('scripts/modules')

#import importlib.util
#spec = importlib.util.spec_from_file_location("goo", "../scripts/modules/goo/goo.py")
#foo = importlib.util.module_from_spec(spec)
#spec.loader.exec_module(foo)
#print(foo.var)

# importing required module
#import goo
import bpy

import os
cwd = os.getcwd()

def test_sample_func():
    print("test sample function------------")

    version = bpy.app.version_string
    print (version)
 
    print ("AAAAAA")
    #print (dir(goo))
    print ("BBBBBBB")
   # ret = goo.sample_func()
   # print (ret)
    print("Test sample func..........")

    print (cwd)
    print("FFFFFF")
    assert True  


def test_calculate_volume():
    #obj = ???? need to get Blender object to test with somehow
    #volume = calculate_volume(obj)
    #assert volume = ??
    print ("test calculate volume")

    assert True

    
