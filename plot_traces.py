#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

def main():

    timestamp = "1545167676"

    connections = [True,False]
    targets = ['PY', 'HTC', 'RTC', 'FS', 'IN', 'RE', 'TC']
    left_options = [True,False]

    heatmaps(targets, connections, left_options, timestamp)
    voltage_traces(targets, connections, left_options, timestamp)
    voltage_traces(targets, connections, left_options, timestamp, average = True)

def get_output_stem(timestamp):

    return "output/figures/" + timestamp

def heatmaps(targets, connections, left_options, timestamp):

    for target in targets:
        for connected in connections:
            for left in left_options:

                fig = plt.figure(figsize = (11, 8.5))

                stem = get_filename_stem(target, timestamp, connected, left)
                volt_filename = stem + "volt.npy"

                data = np.load(volt_filename)

                #normalize plot height
                num_neurons = data.shape[0]
                aspect = 20000 / num_neurons

                #determine maximum magnitude for centering colormap around 0
                m = 0
                for x in data.flat:
                    if abs(x) > m:
                        m = abs(x)

                plt.imshow(data, aspect = aspect, cmap = 'seismic', vmin = -m, vmax = m)
                plt.colorbar(shrink = 0.35)
                plt.xlabel('Time')
                plt.ylabel('Neuron index')

                output_dir = get_output_stem(timestamp) + "/heatmaps/"
                if(connected):
                    output_dir = output_dir + "connected/"
                else:
                    output_dir = output_dir + "unconnected/"

                if(left):
                    output_dir = output_dir + "left/"
                else:
                    output_dir = output_dir + "right/"

                os.makedirs(output_dir, exist_ok = True)

                output_filename = output_dir + target

                if(connected):
                    output_filename = output_filename + "_connected"
                else:
                    output_filename = output_filename + "_unconnected"
                if(left):
                    output_filename = output_filename + "_left"
                else:
                    output_filename = output_filename + "_right"


                output_filename = output_filename + "_heatmap.pdf"
                print(output_filename)
                fig.savefig(output_filename)
                plt.close()






def get_filename_stem(target, timestamp, connected, is_left):

    if(connected):
        stem = "output/" + timestamp + "/connected/"
    else:
        stem = "output/" + timestamp + "/unconnected/"

    if(is_left):
        name = target + "L_"
    else:
        name = target + "R_"

    return stem + name



def voltage_traces(targets, connections, left_options, timestamp, average = False):

    plt.rcParams.update({'font.size':6})

    for target in targets:

        nrows = len(connections)
        ncols = len(left_options)

        f, axarr = plt.subplots(2,2)

        subplot_counter = 1
        i = 0
        for connected in connections:
            j = 0

            for left in left_options:

                ax = axarr[i,j]
                stem = get_filename_stem(target, timestamp, connected, left)
                time_filename = stem + "time.npy"
                volt_filename = stem + "volt.npy"

                time = np.load(time_filename)
                volt = np.load(volt_filename)

                subplot_code = (100 * nrows) + (10 * ncols) + subplot_counter

                if average:

                    num_neurons = volt.shape[0]
                    num_points = volt.shape[1]
                    values = np.empty(shape = (volt.shape[1]))
                    for k in range(num_points):

                        s = 0
                        for l in range(num_neurons):
                            s = s + volt[l,k]
                        values[k] = s / num_neurons


                    ax.plot(time, values, linewidth = 0.5, color = "#3333cc")

                else:
                    for x in range(volt.shape[0]):
                        ax.plot(time, volt[x], linewidth = 0.5, color = "#3333cc",
                                alpha = 0.5)

                xlabel = "Time (seconds)"
                ylabel = "Voltage (volts)"

                ax.set(xlabel = xlabel, ylabel = ylabel)

                if(j == 0):
                    if(connected):
                        label = "Connected"
                    else:
                        label = "Not Connected"
                    ax.annotate(label, xy=(-0.65,0.5), xycoords=("axes fraction",
                        "axes fraction"), weight = "bold")

                if(i == 0):
                    if(left):
                        label = "Left"
                    else:
                        label = "Right"

                    ax.annotate(label, xy=(0.45,1.1), xycoords=("axes fraction",
                        "axes fraction"), weight = "bold")
                subplot_counter = subplot_counter + 1
                j = j + 1
            i = i + 1

        f.tight_layout(rect=[0.15,0,1,0.9])

        title = target

        if(average):
            title = title + " (Average across %d neurons)" % volt.shape[0]

        f.suptitle(title, x = 0.6, fontsize = 15)

        output_dir = get_output_stem(timestamp) + "/voltage_traces/"
        if(average):
            output_dir = output_dir + "averaged/"
        os.makedirs(output_dir, exist_ok = True)

        output_filename = output_dir + target

        if(average):
            output_filename = output_filename + "_average"

        output_filename = output_filename + "_voltage_trace.pdf"

        print(output_filename)
        f.savefig(output_filename, pad_inches = 10)

if __name__ == "__main__":
    main()




