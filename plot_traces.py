#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def main():

    plt.rcParams.update({'font.size':6})
    connections = [True,False]
    targets = ['PY', 'HTC', 'RTC', 'FS', 'IN', 'RE', 'TC']
    #targets = ['PY']
    left_options = [True,False]

    for target in targets:

        plot(target, connections, left_options)

def plot(target, connections, left_options):

    nrows = len(connections)
    ncols = len(left_options)

    f, axarr = plt.subplots(2,2)

    subplot_counter = 1
    i = 0
    for connected in connections:
        j = 0

        for left in left_options:

            ax = axarr[i,j]

            if(connected):
                stem = "output/connected/"
            else:
                stem = "output/unconnected/"

            if(left):
                name = target + "L_"
            else:
                name = target + "R_"

            time_file = name + "time.npy"
            volt_file = name + "volt.npy"

            time = np.load(stem + time_file)
            volt = np.load(stem + volt_file)

            subplot_code = (100 * nrows) + (10 * ncols) + subplot_counter

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
    f.suptitle(target, x = 0.6, fontsize = 15)

    output_dir = "output/figures/run_2/"
    filename = output_dir + target + ".png"
    print(filename)
    f.savefig(filename, pad_inches = 10)

if __name__ == "__main__":
    main()




