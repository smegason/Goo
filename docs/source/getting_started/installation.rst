.. _installation:

Installation
============

Goo runs within Blender, so you will first need to download Blender 4.1 from `here <https://www.blender.org/download/lts/4-1/>`__.
We aim to maintain Goo for future Blender LTS versions but might have a lag from their releases. Currently, we support Blender 3.3 and above.  

Once you have Blender installed:

1. Download the Goo library from `GitHub <https://github.com/smegason/Goo>`__. Download the zip file of the latest release, and unzip it in an empty folder. \n Alternatively, it can be downloaded from the command line as follows:

   .. code-block:: bash

      mkdir Goo
      cd Goo
      wget <Goo-latest-release>.tar.gz
      tar -xvf <Goo-latest-release>.tar.gz

2. In Blender, go to `Edit > Preferences > File Paths > Scripts` and add `<your_root_path>/Goo/scripts`.

Dependencies
------------

- python 3.7 or newer
- numpy_
- scipy_
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
.. _scipy: https://pypi.org/project/scipy/

Install Dependencies in Blender
-------------------------------

Blender comes with its own Python interpreter, which is isolated from the system's Python environment. 
Goo requires a few additional Python packages that need to be installed directly into Blender's Python environment. 

To install Goo's dependencies: 

1. Find the paths of the Blender executable and its Python interpreter.

   For macOS, it is usually in the Applications folder, e.g., `/Applications/Blender.app/Contents/MacOS/Blender` and `/Applications/Blender.app/Contents/Resources/3.6/python/bin/python3.10`. For Windows, it is usually in the Program Files folder.

2. Create a new environment using Blender's Python interpreter:

   .. code-block:: bash

      /Applications/Blender.app/Contents/Resources/4.1/python/bin/python3.11 -m venv blender_venv

3. Activate the environment:

   .. code-block:: bash

      source blender_venv/bin/activate

4. Install the dependencies:

   .. code-block:: bash

      /Applications/Blender.app/Contents/Resources/3.6/python/bin/python3.10 -m pip install setuptools numpy scipy sphinx sphinx_copybutton furo typing_extensions

5. Check that the dependencies are installed:

   .. code-block:: bash

      /Applications/Blender.app/Contents/Resources/4.1/python/bin/python3.11 -m pip list

6. Launch Blender:

   .. code-block:: bash

      /Applications/Blender.app/Contents/MacOS/Blender
