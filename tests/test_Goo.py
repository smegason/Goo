#test_goo.py
#this file is read by pytest to test all functions in Goo cccccccc

#import Goo

class classA():
    def __init__(self):
        self.var1 = 0
        self.var2 = "Hello"

A = classA()


def test_sample_func():
    print("test samp[le function")
    print ("A")
    print(dir(classA))

    #print ("Goo class:")
    #print (dir(Goo))
    
    #Goo.ghdew()

    #Goo.sample_func()
    #print (ret)
    print("Test sample func")
    assert 1 == 1   

def test_calculate_volume():
    #obj = ???? need to get Blender object to test with somehow
    #volume = calculate_volume(obj)
    #assert volume = ??
    print ("test calculate volume")

    assert True

    
