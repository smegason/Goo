#!/bin/bash

# Creates a Goo folder if it does not exist yet
BL_PATH="$HOME/Goo"

if [ ! -d "$BL_PATH" ]; then
    # if the folder does not exist, create it
    mkdir -p "$BL_PATH"
    echo "Folder created."
else
    echo "Folder already exists."
fi

# Define the name of the software
SOFTWARE_NAME="Goo"
SOFTWARE_URL="https://builder.blender.org/download/daily/blender-3.6.12-stable+v36.626a6b1c6799-linux.x86_64-release.tar.xz"
FILE_NAME="blender-3.6.12-stable+v36.626a6b1c6799-linux.x86_64-release.tar.xz"

# Get Blender 3.6 LTS
echo "Downloading Blender"
# Download the software
wget -P "$BL_PATH" "$SOFTWARE_URL"

# Check if the download was successful
if [ $? -ne 0 ]; then
    echo "Download failed. Exiting."
    exit 1
fi

cd "$BL_PATH"

# Untar the tarball 
echo "Extracting Blender"
tar -xf "$FILE_NAME"

# Verify extraction
EXTRACTED_DIR="$BL_PATH/blender-3.6.12-stable+v36.626a6b1c6799-linux.x86_64-release"
if [ -d "$EXTRACTED_DIR" ]; then
    echo "Extraction completed successfully."
else
    echo "Extraction failed. Exiting."
    exit 1
fi

echo "Removing installer"
rm "$FILE_NAME"

# Make Blender executable
BLENDER_EXEC="$EXTRACTED_DIR/blender"
chmod u+x "$BLENDER_EXEC"

# Create the bin directory in the home directory if it does not exist
BIN_DIR="$HOME/bin"
if [ ! -d "$BIN_DIR" ]; then
    mkdir -p "$BIN_DIR"
    echo "Created $BIN_DIR directory."
else
    echo "$BIN_DIR directory already exists."
fi

# Create a symbolic link to the Blender executable in ~/bin
ln -sf "$BLENDER_EXEC" "$BIN_DIR/blender"

# Ensure ~/bin is in the PATH
if [[ ":$PATH:" != *":$HOME/bin:"* ]]; then
    echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
    export PATH="$HOME/bin:$PATH"
    echo "Added $HOME/bin to PATH."
else
    echo "$HOME/bin is already in the PATH."
fi

# Check if Blender runs
echo "Checking if Blender runs"
blender --version
if [ $? -ne 0 ]; then
    echo "Blender failed to run. Please check the installation."
    exit 1
else
    echo "Blender is installed and running successfully."
fi

echo "Blender is set up and ready to use. You can run it with the command: blender"
