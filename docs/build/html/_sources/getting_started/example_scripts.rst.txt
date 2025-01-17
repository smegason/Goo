.. _example_scripts:

Example scripts
=================

1. Writing your first script: a tutorial
--------------------------------------------

Examples of simulation script can be found in the `/simulations/` folder, located `here <https://github.com/smegason/Goo/tree/main/simulations>`__. 
Goo extends Blender towards agent-based simulations of cell mechanics, molecular reactions and gene regulatory networks.
With Goo, you can create custom scripts to simulate various cellular phenomena by specifying initial conditions and parameters.

Running scripts
^^^^^^^^^^^^^^^^

Goo scripts typically get ran from Blender's scripting tab, though they can be ran from Visual Studio Code directly using the `developer's extension <https://marketplace.visualstudio.com/items?itemName=JacquesLucke.blender-development>`__ developed by Jacques Lucke. 

Initialization
^^^^^^^^^^^^^^^
All simulation scripts should begin with `goo.reset_modules()` and `reset_scene()` to ensure a clean starting environment. 
The latter removes all objects from the Blender scene and resets simulation parameters to their default values.

.. code-block:: python

   import goo

   goo.reset_modules()
   goo.reset_scene()

Next, define cell types along with their physical properties, such as surface stiffness and adhesion strength. 
Stiffer cells are less deformable, while higher adhesion strength increases deformation from their initial spherical shape. 
The ratio of stiffness to adhesion strength influences resulting cell patterns (Garner, Tsai, and Megason, 2022). 
Most parameters in Goo are dimensionless.

Defining cell types
^^^^^^^^^^^^^^^^^^^^

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


Creating cells
^^^^^^^^^^^^^^^

Populate cell types with individual cells using `create_cell()`. Specify their initial location, size, shape, and optional material color for visualization.

.. code-block:: python

   cellsA.create_cell("A1", (-20, +1.75, 0), color=(1, 1, 0), size=1.6)
   cellsA.create_cell("A2", (-20, -1.75, 0), color=(1, 1, 0), size=1.6)

   cellsB.create_cell("B1", (0, +1.75, 0), color=(0, 1, 1), size=1.6)
   cellsB.create_cell("B2", (0, -1.75, 0), color=(0, 1, 1), size=1.6)

   cellsC.create_cell("C1", (20, +1.75, 0), color=(1, 0, 1), size=1.6)
   cellsC.create_cell("C2", (20, -1.75, 0), color=(1, 0, 1), size=1.6)


Setting up the simulator
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To introduce cell behaviors like adhesion, motility and division, use the simulator. It handles the simmulation of cell physics and solving sets of ODEs over time for genetic circuitry. 
Define total simulation time, time step (`physics_dt` for mechanics, `molecular_dt` for reactions), and which cell types to include in the simulation.
If not included, objects will remain static.

.. note::

   Goo uses two simulation engines: one for cell mechanics on meshes and one for discrete molecular reactions on KD-trees and gene regulatory circuitry for each cell. 
   Molecular processes happen at faster time scale than cell mechanics; therefore `physics_dt` typically needs be set at least 10 times larger than `molecular_dt`.

The `setup_world()` function always needs be defined, and it's best practice to always set a random seed for reproducibility. It sets some general parameters (e.g. turning gravity off), units and length scales. 

.. code-block:: python

   sim = goo.Simulator([cellsA, cellsB, cellsC], time=180, physics_dt=1)
   sim.setup_world(seed=2024)


Appending handlers to the simulator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Handlers modularly define cell behavior. They execute functions sequentially at every time step. They can take some parameters as arguments to control e.g. the rate of division based on cell volume. 
Add handlers to the simulator to control these aspects. For example: these lines model cell growth, homotypic adhesion, volume-based division (target volume of 50 :math:`\mu m^3` with a std.dev. of 1) and gaussian random motion. 

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

2. Modeling growing cell doublets with differential adhesion
---------------------------------------------------------------

This example simulates three cell doublets of varying stiffness and adhesion strength over ~150 minutes, ignoring heterotypic adhesion.

.. literalinclude:: ../examples/1_growing_doublets.py
   :language: python

Running this script in Blender produces the following simulation:

.. video:: ../examples/1_growing_doublets.mp4
   :width: 740
   :loop:


3. Cell division based on volume growth
------------------------------------------

Simulates cells dividing based on volume. 

.. admonition:: Goo script
   :class: dropdown

   .. literalinclude:: ../examples/2_volume_dividing_cells.py
      :language: python


.. .. literalinclude:: ../examples/2_volume_dividing_cells.py
..    :language: python

Running this script in Blender produces the following simulation:

.. video:: ../examples/2_volume_dividing_cells0001-0250.mp4
   :width: 740
   :loop:
   

4. Random cell mixing 
------------------------

Simulates random motion and mixing of cells. 


.. admonition:: Goo script
   :class: dropdown

   .. literalinclude:: ../examples/3_random_cell_mixing.py
      :language: python

Running this script in Blender produces the following simulation:

.. video:: ../examples/3_random_cell_mixing0001-0200.mp4
   :width: 740
   :loop:


5. Blinking cells and the repressilator 
-----------------------------------------

This example introduces the simulation of gene regulation networks.

We model cells that blink on and off, representing their gene concentration that is controlled by a gene regulatory network akin to a repressilator (Elowitz et al., 2000). 
To do so, we need to define a set of genes and their interactions. Here, gene X repress gene Z, gene Z repress gene Y, and gene Y repress gene X, and all genes degrades following first order equation. That is, the degradation of a gene concentration proceeds at a rate that is linearly proportional to the gene concentration itself. 
We declare two repressilator networks, one with a faster production rate than the other i.e. higher `kcat`, leading to a shorter period of oscillation. Cells are then populated with these networks and initial gene concentrations.

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

   network2 = GeneRegulatoryNetwork()
   network2.load_circuits(
      DegFirstOrder(x, 0.1),
      DegFirstOrder(y, 0.1),
      DegFirstOrder(z, 0.1),
      ProdRepression(y, x, kcat=0.4, n=3),
      ProdRepression(z, y, kcat=0.4, n=3),
      ProdRepression(x, z, kcat=0.4, n=3),
   )

   cell1.grn = network1
   cell1.metabolites = {x: 2, y: 0.1, z: 0.1}
   cell2.grn = network2
   cell2.metabolites = {x: 2, y: 0.1, z: 0.1}


When put all together, this is the script that models blinking cells:

.. admonition:: Goo script
   :class: dropdown

   .. literalinclude:: ../examples/4_blinking_cells.py
      :language: python


Running this script in Blender produces the following simulation: 

.. video:: ../examples/4_blinking_cells0001-0200.mp4
   :width: 740
   :loop:



6. Controlling cell motility with gene concentration
------------------------------------------------------

In this example, the strength of the motion of a cell is controlled by the concentration of gene x, 
which itself is modeled as part of a repressilator regulatory network. The concentration of gene x oscillates; 
therefore, the speed of the cell oscillates between exploratory (fast motility) and stationary (low motility) periods.  

The script to model this behavior: 

.. admonition:: Goo script
   :class: dropdown

   .. literalinclude:: ../examples/5_motion_controlled_by_gene_conc.py
      :language: python


Running this script in Blender produces the following simulation: 

.. video:: ../examples/5_motion_controlled_by_gene_conc0001-1000.mp4
   :width: 740
   :loop: