.. _installation:

Installation
============

1. Install Blender

Goo runs within Blender, so you will first need to download Blender 4.0 from `here <https://www.blender.org/download/lts/4-1/>`__.
Goo currently support Blender 3.3 to 4.0. We aim to maintain Goo for future Blender LTS versions but might have a slight lag from their releases. 

.. note::

   Goo is not currently not compatible with Blender 4.1 because of a dependency clash with the `RoadRunner Simulation Engine <https://libroadrunner.readthedocs.io/en/latest/index.html>`__. Use Blender 4.0 while we are working towards implementing biocircuitery in Goo. 


2. Download the Goo library from `GitHub <https://github.com/smegason/Goo>`__. Download the zip file of the latest release, and unzip it in an empty folder. \n Alternatively, it can be downloaded from the command line as follows:

   .. code-block:: bash

      mkdir Goo
      cd Goo
      wget <Goo-latest-release>.tar.gz
      tar -xvf <Goo-latest-release>.tar.gz

2. In Blender, go to `Edit > Preferences > File Paths > Scripts` and add `<your_root_path>/Goo/scripts`.

Dependencies
------------

- python 3.7 or newer
- bpy_ (bundled with Blender)
- numpy_
- scipy_
- antimony_
- libroadrunner_
- h5py_


.. _bpy: https://docs.blender.org/api/current/info_advanced_blender_as_bpy.html
.. _numpy: http://www.numpy.org/
.. _scipy: https://scipy.org/
.. _antimony: https://tellurium.readthedocs.io/en/latest/antimony.html
.. _libroadrunner: https://www.libroadrunner.org/
.. _h5py: https://www.h5py.org/

Install dependencies in Blender
------------------------------------

Blender comes with its own Python interpreter, which is isolated from the system's Python environment. 
Goo requires some additional packages that must installed and then be exposed to Blender.

MacOS/Linux
^^^^^^^^^^^

For MacOS and Linux, Goo comes packaged with a Makefile to streamline dependency installation. To use it:

1. Set the Makefile variables `BLENDER_PATH` and `BPY_PATH` to the paths of the Blender executable and its Python interpreter, respectively.

   For macOS, it is usually in the Applications folder, e.g.: 

   .. code-block:: bash

      BLENDER_PATH = /Applications/Blender.app/Contents/MacOS/Blender
      BPY_PATH = /Applications/Blender.app/Contents/Resources/4.0/python/bin/python3.10

2. Set up the Blender environment:

   .. code-block:: bash

      make setup
   
3. Allow Blender to access the dependency folder:

   In Blender, go to `Edit > Preferences > File Paths > Scripts` and add `<your_root_path>/Goo/hook/scripts`.

Windows
^^^^^^^

For Windows, the setup must be done manually.

1. Find the paths of the Blender executable and its Python interpreter.

   These are usually found in Program Files, e.g., `C:\\Program Files\\Blender Foundation\\Blender 4.0\\Blender.exe` and `C:\\Program Files\\Blender Foundation\\Blender 4.0\\4.0\\python\\bin\\python.exe`.

2. Create a new virtual environment using Blender's Python interpreter:

   .. code-block:: bash

      C:\\Program Files\\Blender Foundation\\Blender 4.0\\4.0\\python\\bin\\python.exe -m venv .blender_env


3. Install dependencies into the virtual environment

   .. code-block:: bash

      .blender_env\\bin\\python.exe -m pip install -r requirements.txt

4. Check that the dependencies are installed:

   .. code-block:: bash

      .blender_env\\bin\\python.exe -m pip list

5. Create a "hook" folder that enables the installed packages to be exposed to Blender.

   .. code-block:: bash

      mkdir hook\\scripts\\modules
      xcopy .blender_venv\\lib\\python3.10\\site-packages\\* hook\\scripts\\modules /E /H /I

6. Allow Blender to access the dependency folder:

   In Blender, go to `Edit > Preferences > File Paths > Scripts` and add `<your_root_path>/Goo/hook/scripts`.