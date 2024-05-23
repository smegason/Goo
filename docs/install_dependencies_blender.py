import os
import subprocess

# Define the Blender Python path as a global variable
BLENDER_PYTHON_PATH = "/Applications/Blender.app/Contents/Resources/3.3/python/bin/python3.10"  # here for MacOS


def get_blender_python_path():

    if not os.path.exists(BLENDER_PYTHON_PATH):
        raise FileNotFoundError("Blender's Python interpreter not found at: " + BLENDER_PYTHON_PATH)
    
    return BLENDER_PYTHON_PATH


def install_package(python_executable, package):
    subprocess.check_call([python_executable, "-m", "ensurepip"])
    subprocess.check_call([python_executable, "-m", "pip", "install", "--upgrade", "pip"])
    subprocess.check_call([python_executable, "-m", "pip", "install", package])


def main():
    try:
        blender_python = get_blender_python_path()
        packages = ["scipy", "numpy", "pandas", "matplotlib", "sphinx", "sphinx_rtd_theme", "sphinx_copybutton", "furo", "typing_extensions"]

        for package in packages:
            print(f"Installing {package} using {blender_python}")
            install_package(blender_python, package)

        print("Dependencies have been successfully installed in Blender's Python environment.")

    except Exception as e:
        print("An error occurred: ", str(e))


if __name__ == "__main__":
    main()
