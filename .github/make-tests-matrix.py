jobs = [
    {
        "version": "3.0.0",
        "version_x_y": "3.0",
        "sha": "released",
        "download_url": "https://download.blender.org/release/Blender3.0/blender-3.0.0-linux-x64.tar.xz",
    }
    # {'version': '', 'version_x_y': '', 'download_url': ''},
]

if __name__ == "__main__":
    matrix = {"include": jobs}
    print(f"matrix={matrix}")
