import matplotlib.pyplot as plt
import matplotlib
import json
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.pyplot import xticks
import sys
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

with open(f"{sys.argv[1]}.json", 'r') as f:
    print(f"Data locator: {sys.argv[1]}.json")
    master_dict = json.load(f)

# plot when .json file contains data for multiple runs (then plots a legend)
if len(master_dict["Distances"]) > 1: 

#------------------------ DISTANCE BETWEEN CELLS, MULTIPLE -------------------------------#
    
    if master_dict.get("Distances") is not None and master_dict.get("Times") is not None:
        fig, ax = plt.subplots()
        enum = list(range(len(master_dict["Distances"])))
        times_tmp = master_dict['Times']
        times = [[time / 1000000 for time in times_tmp[i]] for i in enum]
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
        clbr.set_label('irerations')

        plt.ylabel(r'distance between cells [$\mu m$]')
        plt.xlabel(r'time [s]')
        #plt.legend(title = r'$falloff\_power = 0;   strength = -1000$')
        ax.grid(False)
        locs, labels = xticks()

        plt.savefig(f"{sys.argv[1]}_multiple.png", dpi=500)


#--------------------------------- PHASE DIAGRAMS ----------------------------------------#
    
    if master_dict.get("Distances") is not None and master_dict.get("Tension") is not None and master_dict.get("Adhesion") is not None: 
        tension = master_dict['Tension']
        adhesion = master_dict['Adhesion']
        steady_distance = [dist[-1] for dist in master_dict['Distances']]

        data_array = np.array(steady_distance).reshape(len(np.unique(adhesion)), len(np.unique(tension)))

        # Create the heatmap using imshow()
        fig, ax = plt.subplots()
        im = ax.imshow(data_array)

        # Set the tick labels for the x and y axes
        ax.set_yticks(np.arange(len(np.unique(adhesion))))
        ax.set_xticks(np.arange(len(np.unique(tension))))
        ax.set_yticklabels(np.unique(np.abs(adhesion)))
        ax.set_xticklabels(np.unique(tension))

        # Rotate the tick labels and set their alignment
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        # Add a colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = ax.figure.colorbar(im, cax=cax)
        cbar.set_label(r'distance between cells [$\mu m$]')


        # Add labels and title to the plot
        plt.ylabel(r'adhesion strength [$-$]')
        plt.xlabel(r'tension stiffness [$-$]')
        fig.tight_layout()

        # save image
        plt.savefig(f"{sys.argv[1]}_heatmap.png", dpi=500)

# plot data for single run
else: 

#----------------------------- DISTANCE BETWEEN CELLS ------------------------------------#
    data_needed = ["Distances", "Times", "Frames"]
    if all(key in master_dict for key in data_needed) and all(master_dict[key] for key in data_needed):
    
        fig, ax1 = plt.subplots(1, figsize = (9, 3.7))

        enum = list(range(len(master_dict["Distances"])))
        times_tmp = master_dict['Times']
        times = [[time / 1000000 for time in times_tmp[i]] for i in enum]
        distances = master_dict['Distances']
        frames = master_dict['Frames']

        # mapping from time to frame
        df = pd.DataFrame({'x': times[0], 'y': frames[0]})
        # fit a polynomial curve on forward data
        poly3_forward = np.poly1d(np.polyfit(df.x, df.y, deg = 5))
        # fit a polynomial curve on reverse data
        poly3_reverse = np.poly1d(np.polyfit(df.y, df.x, deg = 5))

        # declare colors for multiple plots
        hsv2rgb = lambda hue: mcolors.hsv_to_rgb([hue,0.9,0.7])
        hues = np.linspace(0, 0.7, len(enum))
        colors = [hsv2rgb(hue) for hue in hues]

        # plot data
        ax1.plot(times[0], distances[0], color=colors[0])

        # secondary axis: frame number
        ax2 = ax1.secondary_xaxis('top', functions=(poly3_forward, poly3_reverse))
        ax2.set_xlabel(r'frame')
        ax2.tick_params(axis='both', which='both', length=0)

        # layout
        plt.ylabel(r'distance between cells [$\mu m$]')
        plt.xlabel(r'time [s]')
        plt.setp(ax1.get_xticklabels(), visible=True)
        ax1.grid(False)
        locs, labels = xticks()

        # save image
        plt.savefig(f"{sys.argv[1]}_single.png", dpi=500)

#---------------------------------- DISPLACEMENT -----------------------------------------#
    data_needed = ["Displacement"]
    if all(key in master_dict for key in data_needed) and all(master_dict[key] for key in data_needed):

        enum = ('x','y','z')

        hsv2rgb = lambda hue: mcolors.hsv_to_rgb([hue,0.9,0.7])
        hues = np.linspace(0, 0.7, len(enum))
        colors = [hsv2rgb(hue) for hue in hues]


        for cell_name, data in master_dict.get('Displacement').items():
            
            fig, ax1 = plt.subplots(1, figsize = (9, 4))

            x = [d[0] for d in data]
            y = [d[1] for d in data]
            z = [d[2] for d in data]
            
            # plot
            ax1.plot(times[0], x, color = colors[0], linestyle='dashed', marker = 'o', ms = 3, label = f'{cell_name}, x')    
            ax1.plot(times[0], y, color = colors[1], linestyle = 'solid', marker = 'v', ms = 3, label = f'{cell_name}, y')    
            ax1.plot(times[0], z, color = colors[2], linestyle = 'dotted', marker = 's', ms = 2, label = f'{cell_name}, z')   

            # secondary axis: frame number
            ax2 = ax1.secondary_xaxis('top', functions=(poly3_forward, poly3_reverse))
            ax2.set_xlabel(r'frame')
            ax2.tick_params(axis='both', which='both', length=0) 

            # layout
            plt.xlabel(r'time [s]')
            plt.ylabel(r'displacement [$\mu m$]')
            plt.setp(ax1.get_xticklabels(), visible=True)
            ax1.grid(False)
            locs, labels = xticks()
            plt.legend()

            # save image
            plt.savefig(f"{sys.argv[1]}_displacement_{cell_name}.png", dpi=500)

#--------------------------------- DEFORMABILITY ----------------------------------------#
    data_needed = ["Deformability"]
    if all(key in master_dict for key in data_needed) and all(master_dict[key] for key in data_needed):       

        for cell_name, data in master_dict.get('Deformability').items():
            
            fig, ax1 = plt.subplots(1, figsize = (9, 4))

            x = [d[0] for d in data]
            y = [d[1] for d in data]
            z = [d[2] for d in data]
            
            # plot
            ax1.plot(times[0], x, linestyle='dashed', marker = 'o', ms = 3, label = f'{cell_name}, x')    
            ax1.plot(times[0], y, linestyle = 'solid', marker = 'v', ms = 3, label = f'{cell_name}, y')    
            ax1.plot(times[0], z, linestyle = 'dotted', marker = 's', ms = 2, label = f'{cell_name}, z')    

            # secondary axis: frame number
            ax2 = ax1.secondary_xaxis('top', functions=(poly3_forward, poly3_reverse))
            ax2.set_xlabel(r'frame')
            ax2.tick_params(axis='both', which='both', length=0) 

            # layout
            plt.xlabel(r'time [s]')
            plt.ylabel(r'shape deformation [$\mu m$]')
            plt.setp(ax1.get_xticklabels(), visible=True)
            ax1.grid(False)
            locs, labels = xticks()
            plt.legend()

            # save image
            plt.savefig(f"{sys.argv[1]}_deformability_{cell_name}.png", dpi=500)


#--------------------------------- CONTACT AREA ----------------------------------------#

    data_needed = ["Contact area", "Times", "Frames"]
    if all(key in master_dict for key in data_needed) and all(master_dict[key] for key in data_needed):

        fig, ax1 = plt.subplots(1, figsize = (9, 3.7))

        enum = list(range(len(master_dict["Contact area"])))
        times_tmp = master_dict['Times']
        times = [[time / 1000000 for time in times_tmp[i]] for i in enum]
        areas = master_dict['Contact area']
        frames = master_dict['Frames']

        # declare colors for multiple plots
        hsv2rgb = lambda hue: mcolors.hsv_to_rgb([hue,0.9,0.7])
        hues = np.linspace(0, 0.7, len(enum))
        colors = [hsv2rgb(hue) for hue in hues]

        # plot data
        ax1.plot(times[0], areas[0], color=colors[0])

        # secondary axis: frame number
        ax2 = ax1.secondary_xaxis('top', functions=(poly3_forward, poly3_reverse))
        ax2.set_xlabel(r'frame')
        ax2.tick_params(axis='both', which='both', length=0)

        # layout
        plt.ylabel(r'contact area ratio [$A_0/A_n$]')
        plt.xlabel(r'time [s]')
        plt.setp(ax1.get_xticklabels(), visible=True)
        ax1.grid(False)
        locs, labels = xticks()

        # save image
        plt.savefig(f"{sys.argv[1]}_contact.png", dpi=500)

#--------------------------------- VOLUME ----------------------------------------#
    data_needed = ["Volume", "Times", "Frames"]
    if all(key in master_dict and len(master_dict[key]) > 0 for key in data_needed):

        enum = list(range(len(master_dict["Times"])))
        times_tmp = master_dict['Times']
        times = [[time / 1000000 for time in times_tmp[i]] for i in enum]
        volumes = master_dict['Volume']
        frames = master_dict['Frames']

        # mapping from time to frame
        df = pd.DataFrame({'x': times[0], 'y': frames[0]})
        # fit a polynomial curve on forward data
        poly3_forward = np.poly1d(np.polyfit(df.x, df.y, deg = 5))
        # fit a polynomial curve on reverse data
        poly3_reverse = np.poly1d(np.polyfit(df.y, df.x, deg = 5))


        # declare colors for multiple plots
        hsv2rgb = lambda hue: mcolors.hsv_to_rgb([hue,0.9,0.7])
        hues = np.linspace(0, 0.7, len(enum))
        colors = [hsv2rgb(hue) for hue in hues]

        for cell_name, data in master_dict.get('Volume').items():
        
            fig, ax1 = plt.subplots(1, figsize = (9, 4))

            # plot
            ax1.plot(times[0][:-1], volumes.get(f'{cell_name}'), color = colors[0], label = f'{cell_name}')    

            # Add horizontal lines with the first and last values of volumes
            first_vol = volumes.get(f'{cell_name}')[0]
            last_vol = volumes.get(f'{cell_name}')[-1]
            ax1.axhline(y=first_vol, linestyle='--', color='black', linewidth=1)
            ax1.axhline(y=last_vol, linestyle='--', color='black', linewidth=1)
            
            # Add ticks to the y-axis at the locations of the horizontal lines
            yticks = ax1.get_yticks().tolist()
            yticks.append(first_vol)
            yticks.append(last_vol)
            yticks = sorted(list(set(yticks)))
            ax1.set_yticks(yticks)

            # secondary axis: frame number
            ax2 = ax1.secondary_xaxis('top', functions=(poly3_forward, poly3_reverse))
            ax2.set_xlabel(r'frame')
            ax2.tick_params(axis='both', which='both', length=0) 

            # layout
            plt.xlabel(r'time [s]')
            plt.ylabel(r'volume [$\mu m ^ 3$]')
            plt.setp(ax1.get_xticklabels(), visible=True)
            ax1.grid(False)
            locs, labels = xticks()
            plt.legend()

            # save image
            plt.savefig(f"{sys.argv[1]}_volume_{cell_name}.png", dpi=500)


'''
#--------------------------------- LONG AXIS ----------------------------------------#

    for cell_name, data in master_dict.get('Axis length').items():
        
        fig, ax1 = plt.subplots(1, figsize = (9, 4))

        # plot
        ax1.plot(times[0], data, color = colors[0], label = f'{cell_name}')    

        # secondary axis: frame number
        ax2 = ax1.secondary_xaxis('top', functions=(poly3_forward, poly3_reverse))
        ax2.set_xlabel(r'frame')
        ax2.tick_params(axis='both', which='both', length=0) 

        # layout
        plt.xlabel(r'time [s]')
        plt.ylabel(r'long axis length [$\mu m$]')
        plt.setp(ax1.get_xticklabels(), visible=True)
        ax1.grid(False)
        locs, labels = xticks()
        plt.legend()

        # save image
        plt.savefig(f"{sys.argv[1]}_axis_length_{cell_name}.png", dpi=500)


'''
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
