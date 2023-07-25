import subprocess
import numpy as np

#multiple = False

'''if multiple: 

    for adhesion in [*range(-900, -1100, -5)]: 
        for tension in [*np.arange(0, 10, 0.2)]: 
            subprocess.run(
                ["C:\\Program Files\\Blender Foundation\\Blender 3.3\\blender.exe",
                "--python",
                "C:\\Users\\anr9744\\Projects\\Goo\\simulations\\three_cells_benchmark copy.py",
                "--",
                "--adhesion", f"{adhesion}", 
                "--tension", f"{tension}"])

else: '''

for seed in range(1, 101, 1): 
    subprocess.run(
        ["C:\\Program Files\\Blender Foundation\\Blender 3.3\\blender.exe",
        "--python",
        "C:\\Users\\anr9744\\Projects\\Goo\\simulations\\motion_MSD.py",
        "--",
        "--seed", f"{seed}"])


# additional "--" is required to separate the arguments used by blender and python
# "--falloff", f"{falloff}"])

#for i in range(1): 
#    subprocess.run(["C:\\Program Files\\Blender Foundation\\Blender 3.3\\blender.exe", "--python", "C:\\Users\\anr9744\\Projects\\Goo\\simulations\\benchmark_two_cells.py"])

#subprocess.run(["C:\\Program Files\\Blender Foundation\\Blender 3.3\\blender.exe", "--python", "C:\\Users\\anr9744\\Projects\\Goo\\simulations\\benchmark_two_cells.py"])
