import shutil
import sys


def which_blender_by_os():
    """Get the expected Blender executable location by operative system."""
    return (
        shutil.which("Blender")
        if sys.platform == "darwin"
        else ("blender.exe" if "win" in sys.platform else shutil.which("blender"))
    )