
User's Guide
============

Installation
-------------
1. Install the latest edition of Blender from https://www.blender.org/download/
2. Clone the Goo repository or download and unzip
3. In Blender, go to `Edit>Preferences`, then go to the 'Add-ons' tab and enable the checkbox next to `Add Mesh: Extra Objects`

   .. image:: ../img/blender_edit_preferences.png
   .. image:: ../img/blender_add_mesh.png

4. Then, in `Edit>Preferences`, go to the `File Paths` tab and add the `<location>/Goo/scripts/` folder to `Scripts`. `<location>` should be replaced with the root leading up to where you cloned the Goo repository. You may need to close and re-open Blender afther this change.

   .. image:: ../img/blender_add_path.png

Usage
-----
In a new `General` Project within Blender, delete the default cube by left-clicking on the cube, type `X`, and then `return`

In the Scripting tab of Blender, use the desired cell functions

Example script (Create cell)::

   from goo import goo
   goo.setup_world()
   cell = goo.Cell(name_string = "Cell_", loc = (0, 0, 0))
   goo.make_cell(cell)

Click the play button in the scripting tab of Blender and you should see a Goo cell appear
  
Note, Blender uses its own built in Python interpreter which may be different than another Python instance you've installed.

Examples
--------
Look in the Simulations folder in the Github repository for Goo to find example scripts for running simulations. These can be loaded from Blender in the scripts tab and then run to create the simulation.