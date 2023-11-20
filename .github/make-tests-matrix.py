# import re
# from urllib import request


jobs = [
    {
        "version": "3.3.0",
        "version_x_y": "3.3",
        "sha": "released",
        "download_url": "https://download.blender.org/release/Blender3.3/blender-3.3.0-linux-x64.tar.xz",
    },
    {
        "version": "3.3.12",
        "version_x_y": "3.3",
        "sha": "released",
        "download_url": "https://download.blender.org/release/Blender3.3/blender-3.3.12-linux-x64.tar.xz",
    },
    # {'version': '', 'version_x_y': '', 'download_url': ''},
]


'''def get_daily_builds(jobs: list):
    resp = request.urlopen("https://builder.blender.org/download/daily/")
    page = resp.read().decode("utf-8")
    regex_pattern = (
        r"(https://builder.blender.org/download/daily/"
        r"blender-(((?:3|4)\.\d)\.\d-\w+)\+\S{1,6}\.(\S{12})-"
        r"linux\.x86_64-release\.tar\.xz)"
    )   
    releases = re.findall(
        regex_pattern,
        page,
    )
    for release in releases:
        new_job = {
            "version": release[1],
            "version_x_y": release[2],
            "download_url": release[0],
            "sha": release[3],
        }
        if new_job["version"] in [job["version"] for job in jobs]:
            continue

        if (
            new_job["version"].removesuffix("-stable")
            not in [job["version"] for job in jobs]
        ):
            jobs.append(new_job)'''


if __name__ == "__main__":
    # get_daily_builds(jobs)
    matrix = {"include": jobs}
    print(f"matrix={matrix}")
