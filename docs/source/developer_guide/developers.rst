Developer's guide
====================

Install libraries
-----------------------------

1. Fork the github repository

2. Set-up environment

.. code-block:: bash

    conda create --name goo-dev --file requirements.txt

.. code-block:: bash

    conda activate goo-dev

3. Install the Blender `developer's extension <https://marketplace.visualstudio.com/items?itemName=JacquesLucke.blender-development>`__ developed by Jacques Lucke. 

4. Blender test suites will be implemented in the future. In the meantime, create a PR and the Goo team wil review it. Feel free to reach out to sean_megason AT hms.harvard.edu to get involved. 

Update the documentation
-----------------------------

1. Install dependencies about outlined above

2. In your local Goo codebase

.. code-block:: bash

    cd docs

3. Modify source files: make modifications to the source .rst files or to the codebase. Code documentation follows `sphynx <https://www.sphinx-doc.org/en/master/>`__ syntax. 

4. Running the build
Due to some dependencies from Blender (eg `bpy`), the build has to be ran from Blender using its Python interpreter. The docs directory contains a modified `Makefile` so that it manually uses Blender executable. It should be ready-to-use for MacOS users. Other platform users should modify the 'SPHINXBUILD' variable so that it points towards a Blender executable file path in `Makefile`. 

.. literalinclude:: ../../Makefile
   :language: make

Then, running the build will update the documentation in the `/build` directory.

.. code-block:: bash

    make clean # clean up the build directory
    make html # populate the build directory 

5. Push the `docs/build/html`` directory to the `gh-pages` branch
```bash
cd docs/build/html
git add .
git commit -m "Update docs"
git remote add origin
git push -f origin gh-pages
```