.. _example_scripts:

Example scripts
========================

Simulation scripts
------------------

Examples of simulation script can be found in the /simulations folder, located `here <https://github.com/smegason/Goo/tree/main/simulations>`__. 
Once you get a good grasp of the library, you will be able to write your own Goo scripts and specify lots of initial conditions for your simulations of cells. 

Goo scripts typically get ran from Blender's scripting tab, though they can be ran from Visual Studio Code directly using the `developer's extension <https://marketplace.visualstudio.com/items?itemName=JacquesLucke.blender-development>`__ developed by Jacques Lucke. 

Adherent cell doublets
------------------------

.. literalinclude:: ../examples/2_randomly_moving_doublets.py
   :language: python


Random cell mixing inside a spherical volume
------------------------------------------------

.. literalinclude:: ../examples/3_volume_dividing_cells.py
   :language: python


Dividing cells based on volume
-----------------------------------

.. literalinclude:: ../examples/5_random_cell_mixing.py
   :language: python