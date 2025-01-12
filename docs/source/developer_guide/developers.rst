Developer's guide
====================

We recommend using Visual Studio Code (VSCode) as code editor. 

Contributing to the codebase
-------------------------------

1. Fork the `github repository <https://github.com/smegason/Goo>`__. 

2. Clone your fork repository within VSCode. 

3. Install dependencies: follow the steps of the `installation guide <https://smegason.github.io/Goo/getting_started/installation.html>`__. 

4. Install the `Blender developer's extension <https://marketplace.visualstudio.com/items?itemName=JacquesLucke.blender-development>`__ developed by Jacques Lucke inside VSCode. 

Installing new packages
-----------------------

1. Activate the virtual environment created during installation and install desired packages:
    .. code-block:: bash

        .blender_env\\Scripts\\activate
        pip install <package>

2. Update the `hook/modules` folder:
    macOS:

    .. code-block:: bash

        make update_modules

    Windows:

    .. code-block:: bash

        rmdir /S /Q hook
        mkdir hook\\scripts\\modules
        xcopy .blender_venv\\lib\\python3.10\\site-packages\\* hook\\scripts\\modules /E /H /I
 

Testing your code
--------------------

We developed a suite of tests to ensure the correct functioning of the codebase. They do not cover all possible cases, but they are a good starting point. Also, they do not test the realism and accuracy of the physics simulations. 

1. We use pytest from within Blender's Python interpreter to run tests. Make sure to install the `pytest` library in your Blender Python environment.

2. A suite of test cases can be found inside the `/tests` folder.

3. To run the suite of tests, open a terminal in VSCode and run the following command from the codebase root directory. Make sure to replace the path to the Blender executable with the correct path on your system.
    macOS: 

    .. code-block:: bash

        make test

    or

    .. code-block:: bash

        /Applications/Blender-4.0.app/Contents/MacOS/Blender --background  --python-expr "import pytest; pytest.main(['-v', './tests'])"

    Windows: 

    .. code-block:: bash

        "C:\Program Files\Blender Foundation\Blender 3.0\blender.exe" --background  --python-expr "import pytest; pytest.main(['-v', './tests'])"

4. All tests should pass before submitting a pull request.

5. If you are adding new features, please add tests to cover the new functionality. See the existing test cases for examples.


Update the documentation
-----------------------------

1. Install dependencies about outlined above

2. In your local Goo codebase

.. code-block:: bash

    cd docs

3. Modify source files: make modifications to the source .rst files or to the codebase. Code documentation follows `sphynx <https://www.sphinx-doc.org/en/master/>`__ syntax. 

4. Running the build
Due to some dependencies from Blender (i.e. `bpy`), the build has to be ran from Blender using its Python interpreter. The docs directory contains a modified `Makefile` so that it manually uses Blender executable. It should be ready-to-use for MacOS users. Other platform users should modify the 'SPHINXBUILD' variable so that it points towards a Blender executable file path in `Makefile`. 

.. literalinclude:: ../../Makefile
   :language: make

Then, running the build will update the documentation in the `/build` directory. It will also automatically update the library documentation using Sphynx's `automodule`.

.. code-block:: bash

    make clean # clean up the build directory
    make html # populate the build directory 

5. Move the content of the `docs/build/html`` directory to the root of the `gh-pages` branch using a `tmp/` folder on your local machine. Then push the changes to the remote repository with: 
    
.. code-block:: bash

    git add .
    git commit -m "Update docs"
    git push -f