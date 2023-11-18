def test_uppercase():
    assert "loud noises".upper() == "LOUD NOISES"


def test_reversed():
    assert list(reversed([1, 2, 3, 4])) == [4, 3, 2, 1]


def test_some_primes():
    assert 37 in {
        num
        for num in range(2, 50)
        if not any(num % div == 0 for div in range(2, num))
    }


def test_make_cell_with_attributes():
    from goo import goo
    import bpy

    # Define attributes for the cell
    name = "MyCell"
    location = (1.0, 2.0, 3.0)
    obj_type = "CustomType"

    # Check if 'name' is a string and 'location' is a tuple
    assert isinstance(name, str), (f"'name' should be a string"
                                   f" but got {type(name)}")
    assert isinstance(location, tuple), (f"'location' should be a tuple"
                                         f" but got {type(location)}")
    assert isinstance(obj_type, str), (f"'obj_type' should be a string"
                                       f" but got {type(obj_type)}")

    # Create the cell using the make_cell function with specified attributes
    goo.make_cell(name=name, loc=location, type=obj_type)

    # Get the object by its name from Blender's data objects
    obj = bpy.data.objects.get(name)

    # Assert that the object with the specified name exists
    assert obj is not None, f"Object with name '{name}' was not found"

    # Assert that the object matches the specified attributes
    assert obj.name == name, (f"Expected object name '{name}'"
                              f" but got '{obj.name}'")
    assert tuple(obj.location) == location, (f"Expected location '{location}'"
                                             f" but got '{obj.location}'")

    # Check if the 'type' property exists in the collection and matches as expected
    collection = obj.users_collection[0] if obj.users_collection else None

    # First assertion for property existence
    assert collection is not None, f"No collection associated with object '{name}'"

    # Second assertion for property value
    assert collection.get('type') == obj_type, \
        (f"Expected object type '{obj_type}'"
         f" but got '{collection.get('type')}'"
         f" or 'type' property not found in the collection")
