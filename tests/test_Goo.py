#test_goo.py
#this file is read by pytest to test all functions in Goo  xxx

# importing module
import sys
  
# appending a path
sys.path.append('scripts/modules')
  
# importing required module
import goo
import bpy

def test_sample_func():
    print("test sample function------------")

    version = bpy.app.version_string
    print (version)
 
    print ("AAAAAA")
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

    
