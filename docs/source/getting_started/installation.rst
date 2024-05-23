.. _installation:

Installation
============

Goo runs within Blender, so you will first need to download Blender 3.6 LTS from `here <https://www.blender.org/download/lts/3-6/>`__.
The package was developed under Blender 3.6 LTS, so we recommend using this version. We aim to maintain Goo for future Blender LTS versions but might have a lag form their releases. 

Once you have Blender installed:

1. Download the Goo library from `GitHub <https://github.com/smegason/Goo>`__. 

2. In Blender, Edit > Preferences > File Paths > Scripts: add <your_root_path>/Goo/scripts

Dependencies
------------

Goo has the following dependencies.

- python 3.7 or newer
- numpy_
- matplotlib_
- bpy_
- bmesh_
- mathutils_
- scipy_

.. _numpy: http://www.numpy.org/
.. _bpy: https://docs.blender.org/api/current/info_advanced_blender_as_bpy.html
.. _bmesh: https://docs.blender.org/api/current/bmesh.html
.. _pandas: http://pandas.pydata.org/
.. _matplotlib: https://matplotlib.org/
.. _json: https://docs.python.org/3/library/json.html
.. _mathutils: https://pypi.org/project/mathutils/
.. _scipy: https://pypi.org/project/scipy/


Install dependencies in Blender
------------------------------------

Blender comes with its own Python interpreter, which is isolated from the system's Python environment. 
Goo requires a few additional Python packages that need to be installed directly into Blender's Python environment. 
Replace the Blender executable path in the sript and run ``python install_dependencies_blender.py`` to install the dependencies. 


.. literalinclude:: ../../install_dependencies_blender.py
   :language: python