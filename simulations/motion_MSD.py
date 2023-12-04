from goo import goo
import sys
from importlib import reload


def main(seed):
    reload(goo)
    goo.setup_world(seed=seed)

    # Cells A 
    # Define cell A1
    goo.make_cell("cell_A1", loc=(0, 0, 0), type="cellsA")

    # Forces A
    # Define force A1
    goo.add_motion('cell_A1', -200, distribution='uniform', size=5)

    # Simulation setup 
    handlers = goo.handler_class()
    handlers.launch_simulation(start=1,  # default, 1x
                               end=5000,  # default, 250
                               filepath=f"C:\\tmp\\differential_motion\\single_cell_motion_seed_1to100_uni10\\data_{seed}", 
                               adhesion=True,  # default, True
                               data=True,  # default, False
                               growth=True, 
                               division=False, 
                               motility=True
                               )


if __name__ == "__main__":
    # Check if the "--seed" argument is provided
    if "--seed" in sys.argv:
        seed_idx = sys.argv.index("--seed") + 1
        seed = int(sys.argv[seed_idx])

        # Call the main function with the provided seed
        main(seed)
    else:
        print("Error: Seed argument is missing. Please provide the --seed argument.")
