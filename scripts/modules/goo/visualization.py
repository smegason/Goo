import matplotlib.pyplot as plt
import matplotlib
import json
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.pyplot import xticks
import sys
import os 


with open(f"{sys.argv[1]}.json", 'r') as f:
    print(f"TEST: {sys.argv[1]}.json")
    #data = f.read()
    master_dict = json.load(f)
#len(master_dict['Times'][0])

if len(master_dict["Distances"]) > 1: 

    fig, ax = plt.subplots()
    enum = list(range(len(master_dict["Distances"])))
    times = master_dict['Times']
    distances = master_dict['Distances']

    hsv2rgb = lambda hue: mcolors.hsv_to_rgb([hue,0.9,0.7])
    hues = np.linspace(0, 0.7, len(enum))
    colors = [hsv2rgb(hue) for hue in hues]

    for idx, n in enumerate(enum):
        ax.plot(times[idx], distances[idx], color=colors[idx])

    # reverse colormap to allow for negative plotting
    cmap = mcolors.LinearSegmentedColormap.from_list('my_cmap', colors)
    norm = mcolors.Normalize(min(enum), max(enum))
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    clbr = fig.colorbar(sm)
    clbr.set_label('Irerations')

    plt.ylabel('Total distance between cells')
    plt.xlabel('Times')
    plt.legend(title = r'$falloff\_power = 0;   strength = -1000$')
    ax.grid(False)
    locs, labels = xticks()

    plt.savefig(f"{sys.argv[1]}_TEST.png", dpi=500)

else: 

    fig, ax = plt.subplots()
    enum = list(range(len(master_dict["Distances"])))
    print(enum)
    times = master_dict['Times']
    distances = master_dict['Distances']

    hsv2rgb = lambda hue: mcolors.hsv_to_rgb([hue,0.9,0.7])
    hues = np.linspace(0, 0.7, len(enum))
    colors = [hsv2rgb(hue) for hue in hues]

    ax.plot(times[0], distances[0], color=colors[0])

    plt.ylabel('Total distance between cells')
    plt.xlabel('Times')
    plt.legend(title = r'$falloff\_power = 0;   strength = -1000$')
    ax.grid(False)
    locs, labels = xticks()

    plt.savefig(f"{sys.argv[1]}.png", dpi=500)
