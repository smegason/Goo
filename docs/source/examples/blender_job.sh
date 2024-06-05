#!/bin/bash
BLENDER_PATH=/path/to/blender  # Update this path
BLEND_FILE=/path/to/your/project.blend  # Update this path
OUTPUT_DIR=/path/to/output/directory  # Update this path

$BLENDER_PATH/blender -b $BLEND_FILE -o $OUTPUT_DIR/frame_##### -a
