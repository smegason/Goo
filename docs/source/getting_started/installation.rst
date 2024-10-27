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
- numpy_
- scipy_
- bpy_ (bundled with Blender)

.. _numpy: http://www.numpy.org/
.. _bpy: https://docs.blender.org/api/current/info_advanced_blender_as_bpy.html
.. _scipy: https://pypi.org/project/scipy/

Install dependencies in Blender
------------------------------------

MacOS/Linux
------------

Blender comes with its own Python interpreter, which is isolated from the system's Python environment. 
Goo requires minimal dependencies that need to be installed directly into Blender's Python environment. 

.. note::

   Blender's Python interpreter changes from Python 3.10 to 3.11 between 3.x and 4.x versions. Make sure to adapt your path accordingly. 


To install Goo's dependencies from the terminal: 

1. Find the paths of the Blender executable and its Python interpreter.

   For macOS, it is usually in the Applications folder, e.g., `/Applications/Blender.app/Contents/MacOS/Blender` and `/Applications/Blender.app/Contents/Resources/4.1/python/bin/python3.11`.

2. Create a new environment using Blender's Python interpreter:

   .. code-block:: bash

      /Applications/Blender.app/Contents/Resources/4.0/python/bin/python3.11 -m venv blender_venv

3. Activate the environment:

   .. code-block:: bash

      source blender_venv/bin/activate

4. Install the dependencies:

   .. code-block:: bash

      /Applications/Blender.app/Contents/Resources/4.0/python/bin/python3.11 -m pip install numpy scipy typing_extensions

5. Check that the dependencies are installed:

   .. code-block:: bash

      /Applications/Blender.app/Contents/Resources/4.0/python/bin/python3.11 -m pip list

6. Launch Blender by running this command in your terminal:

   .. code-block:: bash

      /Applications/Blender.app/Contents/MacOS/Blender



Windows
------------

Blender comes with its own Python interpreter, which is isolated from the system's Python environment. 
Goo requires a few additional Python packages that need to be installed directly into Blender's Python environment. 

To install Goo's dependencies from a terminal: 

1. Find the paths of the Blender executable and its Python interpreter.

   For Windows, it is usually in the Program Files, e.g., `C:\\Program Files\\Blender Foundation\\Blender 4.0\\Blender.exe` and `C:\\Program Files\\Blender Foundation\\Blender 4.0\\4.0\\python\\bin\\python.exe`.

2. Create a new virtual environment using Blender's Python interpreter:

   .. code-block:: bash

      C:\\Program Files\\Blender Foundation\\Blender 4.0\\4.0\\python\\bin\\python.exe -m venv .blender_env


3. Activate the environment:

   .. code-block:: bash

      .blender_env\\Scripts\\activate

4. Install the dependencies:

   .. code-block:: bash

      C:\\Program Files\\Blender Foundation\\Blender 4.0\\4.0\\python\\bin\\python.exe -m pip install numpy scipy typing_extensions

5. Check that the dependencies are installed:

   .. code-block:: bash

      C:\\Program Files\\Blender Foundation\\Blender 4.0\\4.0\\python\\bin\\python.exe -m pip list

6. Launch Blender from within the activated virtual environment:

   .. code-block:: bash

       C:\\Program Files\\Blender Foundation\\Blender 4.0\\Blender.exe


Blender supports virtual environment and the installed packages will be available to use for scripting in Blender. 
