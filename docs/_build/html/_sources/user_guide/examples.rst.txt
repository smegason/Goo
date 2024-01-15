.. _examples:

Examples
========

We provide a couple script examples to start simulating cells, tissues and early embryos and their properties including cell stiffness, 
pressure, adhesion strength (homotypic and heterotypic), motion speed (random walk), growth rate, target volume triggering division and type. 

Cleaving cells
---------------
Simulating growing and dividing cells takes little over 20 lines of code with Goo. This examples shows a single cell dividing and growing until reaching about 100 cells after 250 time steps. Tested on MacBook Pro M2 Max and NVIDIA RTX4090. 

.. code-block:: python


    from goo import goo
    from importlib import reload
    reload(goo)
    goo.setup_world()

    # Cell A
    goo.make_cell("cell_A1", loc=(0, 0, 0), type="cellsA")

    # Simulation setup
    handlers = goo.handler_class()
    handlers.launch_simulation(start=1,  # default, 1
                            end=250,  # default, 250
                            filepath="<your_file_path_including_format_extension>", 
                            adhesion=False,  # default, True
                            growth=True, 
                            motility=False, 
                            division=True, 
                            volume_scale=2
                            )


Homotypically adhering doublets
---------------------------------
This example highlights cell-cell adhesion in Goo. Doublets are declared, then their homotypic forces are declared at different strengths. Homotypic interactions are declared by the `type` argument in `make_cell()`. 

.. code-block:: python

    from goo import goo
    from importlib import reload
    reload(goo)
    goo.setup_world()

    # Cells A 
    goo.make_cell("cell_A1", loc=(-1, 0, 0), type="cellsA")
    goo.make_cell("cell_A2", loc=(1, 0, 0), type="cellsA")
    goo.make_cell("cell_B3", loc=(-1, 2.2, 0), material=('green', 0, .1, 0), type="cellsB")
    goo.make_cell("cell_B4", loc=(1, 2.2, 0), material=('green', 0, .1, 0), type="cellsB")
    goo.make_cell("cell_C5", loc=(-1, 4.4, 0), material=('red', .1, 0, 0), type="cellsC")
    goo.make_cell("cell_C6", loc=(1, 4.4, 0), material=('red', .1, 0, 0), type="cellsC")

    # Forces A 
    homoA = 750
    homoB = 1500
    homoC = 3000

    goo.add_homo_adhesion('cell_A1', -homoA)
    goo.add_homo_adhesion('cell_A2', -homoA)
    goo.add_homo_adhesion('cell_B3', -homoB)
    goo.add_homo_adhesion('cell_B4', -homoB)
    goo.add_homo_adhesion('cell_C5', -homoC)
    goo.add_homo_adhesion('cell_C6', -homoC)

    # Simulation setup
    handlers = goo.handler_class()
    handlers.launch_simulation(start=1,
                               end=500,
                               filepath="<your_file_path_including_format_extension>", 
                               adhesion=True,
                               data=False,
                               growth=True
                               )