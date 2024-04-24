import subprocess
import numpy as np
from enum import Enum
import random


# Define Blender path enumeration
class BlenderPath(Enum):
    GPU_COMP_ROOM = "C:\\Program Files\\Blender Foundation\\Blender 3.3\\blender.exe"
    PERSONAL_MACOS = "/Applications/Blender.app/Contents/MacOS/Blender"
    WINDOWS_PATH = "C:\\Program Files\\Blender Foundation\\Blender 3.3\\blender.exe"


# Set the desired Blender path
blender_path = BlenderPath.PERSONAL_MACOS.value  # Change this to choose the desired path

# Rest of your code
for adhesion in [0, 1000]: # 2000, 4000, 6000, 8000, 10000]:
    for motion in [5000]: #, 1000, 2000, 3000, 4000, 5000]: 

        seed = int(np.random.uniform(0, 1000))
            
        subprocess.run(
            [blender_path,
             "--python", "simulations/scan_sorting.py", 
             "--", 
             "--seed", f"{seed}", 
             "--adhesion", f"{adhesion}", 
             "--motion", f"{motion}"])
