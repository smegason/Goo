.. _installation:

Running Goo on O2
====================

O2 is Harvard Medical School's computing cluster. If you don't have an O2 account, follow RC's `instructions <https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1918304257/How+to+request+or+retain+an+O2+account>`__.

1. Log in to O2

.. code-block:: bash

    ssh <HMS-ID>@o2.hms.harvard.edu

2. Install Blender on O2 
Ranit Karmakar (HMS IAC) has created a bash script to directy install Blender executable on O2. It will only allow to run Blender headless, which is enough if rendering is not requried. Note that `EVEE` is not an available rendering engine on Linux clusters thus it is recommended to use `Cycles`. 

.. literalinclude:: ../examples/setup_blender.sh
   :language: bash

3. Install libmvec (vector math library) from Glibc 2.22w