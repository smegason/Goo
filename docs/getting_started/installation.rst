.. _installation:

Installation
============

Goo runs within Blender, so you will first need to download Blender 3.6 LTS from `here <https://www.blender.org/download/lts/3-6/>`__.
The package was developed under Blender 3.6 LTS, so we recommend using this version. We aim to maintain Goo for future Blender LTS versions but might have a lag form their releases. 

Once you have Blender installed:

1. Download the Goo library from `GitHub <https://github.com/smegason/Goo>`__. 

2. In Blender, Edit > Preferences > File Paths > Scripts: add <your_root_path>/Goo/scripts


Simulation scripts
------------------

Examples of simulation script can be found in the /simulations folder, located `here <https://github.com/smegason/Goo/tree/main/simulations>`__. 
Once you get a good grasp of the library, you will be able to write your own Goo scripts and specify lots of initial conditions for your simulations of cells. 

Goo scripts typically get ran from Blender's scripting tab, though they can be ran from Visual Studio Code directly using the `developer's extension <https://marketplace.visualstudio.com/items?itemName=JacquesLucke.blender-development>`__ developed by Jacques Lucke. 

Dependencies
------------

Goo has the following dependencies.

- python 3.7 or newer
- numpy_
- matplotlib_
- bpy_
- bmesh_
- mathutils_

.. _numpy: http://www.numpy.org/
.. _bpy: https://docs.blender.org/api/current/info_advanced_blender_as_bpy.html
.. _bmesh: https://docs.blender.org/api/current/bmesh.html
.. _pandas: http://pandas.pydata.org/
.. _matplotlib: https://matplotlib.org/
.. _json: https://docs.python.org/3/library/json.html
.. _mathutils: https://pypi.org/project/mathutils/
