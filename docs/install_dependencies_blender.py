import os
import subprocess
import site
import sys

# Define the Blender Python path as a global variable
BL_PYTHON_PATH = "/Applications/Blender.app/Contents/Resources/3.6/python/bin/python3.10"


def get_blender_python_path():
    if not os.path.exists(BL_PYTHON_PATH):
        raise FileNotFoundError("Python interpreter not found at: " + BL_PYTHON_PATH)
    return BL_PYTHON_PATH


def list_installed_packages(python_executable):
    try:
        result = subprocess.run([python_executable, "-m", "pip", "list"], 
                                capture_output=True, text=True)
        print("Currently installed packages:")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error listing installed packages:", e)
        raise


def install_package(python_executable, package):
    subprocess.check_call([python_executable, "-m", "ensurepip"])
    subprocess.check_call([python_executable, "-m", "pip", "install", "--upgrade", "pip"])
    subprocess.check_call([python_executable, "-m", "pip", "install", package])


def main():
    try:
        blender_python = get_blender_python_path()
        list_installed_packages(blender_python)

        packages = ["numpy", 
                    "scipy",
                    "sphinx", 
                    "sphinx_copybutton", 
                    "furo", 
                    "typing_extensions"]

        for package in packages:
            print(f"Installing {package} using {blender_python}")
            install_package(blender_python, package)

        print("Dependencies have been successfully installed in"
              "Blender's Python environment.")

    except Exception as e:
        print("An error occurred: ", str(e))


if __name__ == "__main__":
    main()
