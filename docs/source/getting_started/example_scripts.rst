.. _example_scripts:

Example scripts
=================

Writing your first simulation script 
--------------------------------------

Examples of simulation script can be found in the `/simulations/` folder, located `here <https://github.com/smegason/Goo/tree/main/simulations>`__. 
Goo extends Blender towards agent-based simulations of cell mechanics, molecular reactions and gene regulatory networks.
Once you get a good grasp of the library, you will be able to write your own Goo scripts and specify lots of initial conditions for your simulations of cells. 

Goo scripts typically get ran from Blender's scripting tab, though they can be ran from Visual Studio Code directly using the `developer's extension <https://marketplace.visualstudio.com/items?itemName=JacquesLucke.blender-development>`__ developed by Jacques Lucke. 

Simulation scripts should start with the `reset_scene()` to start with a fresh scene - it removes all objects contained in the scene and resets simulation parameters to default. 

.. code-block:: python

   import goo

   goo.reset_modules()
   goo.reset_scene()

Then, you should declare cell types, along with some of their physical parameters such as surface stiffness, homotypic and heterotypic adhesion strength. 
Stiffer cells will lead to less deformable cells. Higher adhesion strength will lead to larger cell deformation from their initial spherical shape. 
The ratio of stiffness to adhesion strength controls the resulting patterns cells make up (Garner, Tsai and Megason, 2022).
Parameters are dimensionless in Goo.

.. code-block:: python

   cellsA = goo.CellType("A", target_volume=70, pattern="simple")
   cellsB = goo.CellType("B", target_volume=70, pattern="simple")
   cellsC = goo.CellType("C", target_volume=70, pattern="simple")

   cellsA.homo_adhesion_strength = 1
   cellsA.stiffness = 1  # ratio = 1
   cellsB.homo_adhesion_strength = 500
   cellsB.stiffness = 1  # ratio = 500
   cellsC.homo_adhesion_strength = 2000
   cellsC.stiffness = 1  # ratio = 2000


Cell types can then be populated with actual cell instances using `create_cell()`. 
The initial location, size and shape of each cell need to be specified, or set randomly. A material (i.e. a color) can also be given, which comes handy to visualize cell types or track individual cells with a color. 

.. code-block:: python

   cellsA.create_cell("A1", (-20, +1.75, 0), color=(1, 1, 0), size=1.6)
   cellsA.create_cell("A2", (-20, -1.75, 0), color=(1, 1, 0), size=1.6)

   cellsB.create_cell("B1", (0, +1.75, 0), color=(0, 1, 1), size=1.6)
   cellsB.create_cell("B2", (0, -1.75, 0), color=(0, 1, 1), size=1.6)

   cellsC.create_cell("C1", (20, +1.75, 0), color=(1, 0, 1), size=1.6)
   cellsC.create_cell("C2", (20, -1.75, 0), color=(1, 0, 1), size=1.6)


Until now, cells remain static meshes–no cell behaviours have been modelled yet. To do so, you need to call the simulator, 
which handles simulating the cell physics and solving sets of ODEs over time for gene regulation circuitry. 
The simulator takes lists of cell types and reaction-diffusion systems to be included in the simulation. If not included, objects will remain static.
This is also where the total simulation time is specified, as well as the time step used in simulations.

.. note::

   Goo uses two simulation engines: one for cell mechanics on meshes and one for discrete molecular reactions on KD-trees and gene regulatory circuitry for each cell. 
   Molecular processes happen at faster time scale than cell mechanics; therefore `physics_dt` typically needs be set at least 10 times larger than `molecular_dt`.

The `setup_world()` function always needs be defined. It sets some general parameters (e.g. turning gravity off), units and length scales. 

.. code-block:: python

   sim = goo.Simulator([cellsA, cellsB, cellsC], time=180, physics_dt=1)
   sim.setup_world(seed=2024)

Lastly, cell behaviors–cell division, adhesion, motility, signal sensing and transduction–can be modelled in a modular fashion using handlers.
Handlers execute functions sequentially at every time step. They can take some parameters as arguments to control e.g. the rate of division based on cell volume. 
In this example, we model cell growth, homotypic adhesion, volume-based division (target volume of 50 :math:`\mu m^3` with a std.dev. of 1) and gaussian random motion. 

.. code-block:: python

   sim.add_handlers(
      [
         goo.GrowthPIDHandler(),                                           # in um3
         goo.RecenterHandler(),
         goo.SizeDivisionHandler(goo.BisectDivisionLogic, mu=60, sigma=1), # in um3
         goo.RandomMotionHandler(goo.ForceDist.GAUSSIAN, strength=500)
      ]
   )

.. note::
   
   The full list of handlers–cell behaviors the library currently supports–can be found in the codebase documentation. 


Put together, this script models three cell doublets of varying stiffness and homotypic adhesion strength, all moving at the same strength for ~150 minutes. Heterotypic adhesion is ignored here.


.. literalinclude:: ../examples/1_growing_doublets.py
   :language: python

Running this script in Blender should produce the following simulation:

.. video:: ../examples/1_growing_doublets.mp4
   :width: 810
   :loop:

   
Cell division based on volume
------------------------------

.. literalinclude:: ../examples/2_volume_dividing_cells.py
   :language: python

Running this script in Blender should produce the following simulation:

.. video:: ../examples/division_cells0001-0150.mp4
   :width: 810
   :loop:
   

Random cell mixing
-------------------

.. literalinclude:: ../examples/3_random_cell_mixing.py
   :language: python

Running this script in Blender should produce the following simulation:

.. video:: ../examples/3_random_cell_mixing0001-0200.mp4
   :width: 810
   :loop:

Blinking cells
----------------

In this example, we model cells that blink on and off, representing their gene concentration that is controlled by a gene regulatory network akin to a repressilator (Elowitz et al., 2000). 
To do so, we need to define a set of genes and their interactions. Here, gene X repress gene Z, gene Z repress gene Y, and gene Y repress gene X, and all genes degrades at a constant rate.
Other than that, cells still homotypically adhere, move, grow and divide. 

.. code-block:: python

   x = Gene("x")
   y = Gene("y")
   z = Gene("z")

   network1 = GeneRegulatoryNetwork()
   network1.load_circuits(
      DegFirstOrder(x, 0.1),
      DegFirstOrder(y, 0.1),
      DegFirstOrder(z, 0.1),
      ProdRepression(y, x, kcat=0.4, n=3),
      ProdRepression(z, y, kcat=0.4, n=3),
      ProdRepression(x, z, kcat=0.4, n=3),
   )
   cell.grn = network1
   cell.metabolites = {x: 2, y: 0.1, z: 0.1}

When put all together, this is the script that models blinking cells:

.. literalinclude:: ../examples/4_blinking_cells.py
   :language: python


Running this script in Blender should produce the following simulation: 
