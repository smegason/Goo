.. _example_scripts:

Example scripts
========================

Writing your first simulation script 
------------------------------------

Examples of simulation script can be found in the `simulations/` folder, located `here <https://github.com/smegason/Goo/tree/main/simulations>`__. 
Once you get a good grasp of the library, you will be able to write your own Goo scripts and specify lots of initial conditions for your simulations of cells. 

Goo scripts typically get ran from Blender's scripting tab, though they can be ran from Visual Studio Code directly using the `developer's extension <https://marketplace.visualstudio.com/items?itemName=JacquesLucke.blender-development>`__ developed by Jacques Lucke. 

Simulation scripts should start with the `reset_scene()` to start with a fresh scene - it removes all objects contained in the scene and resets simulation parameters to default. 

.. code-block:: python

   import goo

   goo.reset_scene()

Then, you should declare cell types, along with some of their physical parameters such as surface stiffness, homotypic and heterotypic adhesion strength. 
Stiffer cells will lead to less deformable cells. Higher adhesion strength will lead to larger cell deformation from their initial spherical shape. 
The ratio of stiffness to adhesion strength controls the resulting patterns cells make up (Garner, Tsai and Megason, 2022).
Parameters are dimensionless in Goo.

.. code-block:: python

   cellsA = goo.CellType("A")
   cellsB = goo.CellType("B")
   cellsC = goo.CellType("C")

   cellsA.homo_adhesion_strength = 50
   cellsA.stiffness = 10 # ratio = 5
   cellsB.homo_adhesion_strength = 250
   cellsB.stiffness = 5 # ratio = 50
   cellsC.homo_adhesion_strength = 500
   cellsC.stiffness = 1 # ratio = 500


Cell types can then be populated with actual cell instances using `create_cell()`. 
The initial location, size and shape of each cell need to be specified - or set randomly. A material (i.e. a color) can also be given, which comes handy to visualize cell types or track individual cells with a color. 

.. code-block:: python

   cellsA.create_cell("A1", (+1.75, -5, 0), color=(0.5, 0, 0), size=1.6)
   cellsA.create_cell("A2", (-1.75, -5, 0), color=(0.5, 0, 0), size=1.6)

   cellsB.create_cell("B1", (+1.75, 0, 0), color=(0, 0.5, 0), size=1.6)
   cellsB.create_cell("B2", (-1.75, 0, 0), color=(0, 0.5, 0), size=1.6)

   cellsC.create_cell("C1", (+1.75, 5, 0), color=(0, 0, 0.5), size=1.6)
   cellsC.create_cell("C2", (-1.75, 5, 0), color=(0, 0, 0.5), size=1.6)


Until now, cells remain static meshes - no cell behaviours have been modelled yet. To do so, you need to call the simulator, 
which handles simulating the cell physics and solving sets of ODEs over time for gene regulation circuitry. 
The simulator takes lists of cell types and reaction-diffusion systems to be included in the simulation. If not included, objects will remain static.
This is also where the total simulation time is specified, as well as the time step used in simulations.

.. note::

   Goo uses two simulation engines: one for cell mechanics on meshes and one for discrete molecular reactions on KD-trees and gene regulatory circuitry for each cell. 
   Molecular processes happen at faster time scale than cell mechanics; therefore `physics_dt` typically needs be set at least 10 times larger than `molecular_dt`.

The `setup_world()` function always needs be defined. It sets some general parameters (e.g. turning gravity off), units and length scales. 

.. code-block:: python

   sim = goo.Simulator([cellsA, cellsB, cellsC], time=300, physics_dt=1)
   sim.setup_world(seed=2024)

Lastly, cell behaviors can be modelled in a modular fashion using handlers. By a cell behavior, we mean a cellular function such as cell division, adhesion, motility, signal sensing and transduction.
The role of a handler is to execute a function over time. They are executed sequentially at every time step.  
Goo extends Blender towards agent-based simulations of cell mechanics, molecular reactions and gene regulatory networks.
Handlers take some parameters as arguments to control e.g. the rate of division based on time or cell volume. 

.. code-block:: python

   sim.add_handlers(
      [
         goo.GrowthPIDHandler(target_volume=50), # in cubic micrometers
         goo.AdhesionLocationHandler(), # no parameters
         goo.RandomMotionHandler(goo.ForceDist.GAUSSIAN, max_strength=750)
      ]
   )

.. note::
   
   The full list of handlers - cell behaviors the library currently supports - can be found in the codebase documentation. 


Put together, this script models three cell doublets of varying stiffness and homotypic adhesion strength, all moving at the same strength for 300 minutes. Heterotypic adhesion is ignored here.


.. literalinclude:: ../examples/2_randomly_moving_doublets.py
   :language: python

Running this script in Blender should produce the following simulation:

.. video:: ../examples/adherent_cell_doublets0001-0300.mp4
   :width: 810

.. .. raw:: html

..     <div style="max-width: 100%; height: auto;">
..         <video style="width: 100%; height: auto;" controls>
..             <source src="/Users/antoine/Harvard/MegasonLab/github/Goo/docs/source/examples/adherent_cell_doublets0001-0300.mp4" type="video/mp4">
..             Your browser does not support the video tag.
..         </video>
..     </div>

   
Cell division based a threshold volume
------------------------------------------------

.. literalinclude:: ../examples/3_volume_dividing_cells.py
   :language: python

Running this script in Blender should produce the following simulation:

.. video:: ../examples/division_cells0001-0150.mp4
   :width: 810


.. .. raw:: html

..    <div style="max-width: 100%; height: auto;">
..       <video style="width: 100%; height: auto;" controls>
..          <source src="/Users/antoine/Harvard/MegasonLab/github/Goo/docs/source/examples/division_cells0001-0150.mp4" type="video/mp4">
..          Your browser does not support the video tag.
..       </video>
..    </div>


Random cell motion in a confined sphere
------------------------------------------

.. literalinclude:: ../examples/5_random_cell_mixing.py
   :language: python