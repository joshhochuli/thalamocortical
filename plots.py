#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

def main():

    timestamp = "test"

    connections = [True,False]
    targets = ['PY', 'HTC', 'RTC', 'FS', 'IN', 'RE', 'TC']
    left_options = [True,False]

    #heatmaps(targets, connections, left_options, timestamp)
    #voltage_traces(targets, connections, left_options, timestamp)
    #voltage_traces(targets, connections, left_options, timestamp, average = True)
    #percentage(targets, connections, left_options, timestamp)
    percentage_overlayed(targets, connections, left_options, timestamp)

#name directory structure and filename redundantly in case they get separated
def generate_output_filename(timestamp, target, left, connected, plot_type):

    output_dir = "output/figures/" + timestamp + '/' + plot_type + '/'
    filename = target

    if(connected):
        output_dir = output_dir + "connected/"
        filename = filename + "_connected"
    else:
        output_dir = output_dir + "unconnected/"
        filename = filename + "_unconnected"

    if(left):
        output_dir = output_dir + "left/"
        filename = filename + "_left_"
    else:
        output_dir = output_dir + "right/"
        filename = filename + "_right_"

    filename = filename + plot_type + ".pdf"

    os.makedirs(output_dir, exist_ok = True)

    return output_dir + filename

def percentage(targets, connections, left_options, timestamp):

    plt.rcParams.update({'font.size':6})
    fs = 12

    for target in targets:

        nrows = len(connections)
        ncols = len(left_options)

        f, axarr = plt.subplots(2,2)

        subplot_counter = 1
        i = 0
        for connected in connections:
            j = 0

            for left in left_options:

                ind_fig = plt.figure(figsize = (7,5))
                ax = axarr[i,j]
                stem = get_filename_stem(target, timestamp, connected, left)
                volt_filename = stem + "volt.npy"

                volt = np.load(volt_filename)

                subplot_code = (100 * nrows) + (10 * ncols) + subplot_counter

                '''
                ax.plot(time, volt[x], linewidth = 0.5, color = "#3333cc",
                        alpha = 0.5)

                plt.plot(time, volt[x], linewidth = 0.5)
                '''


                plot_data = []

                n_neurons = volt.shape[0]
                total_time = volt.shape[1]

                for time in range(total_time):
                    s = 0
                    for neuron in range(n_neurons):
                        if volt[neuron, time] > -0.051:
                            s = s + 1
                    plot_data.append(float(s) / n_neurons)


                x = range(len(plot_data))
                plt.plot(x,plot_data)
                ax.plot(x,plot_data)
                plt.ylim([0,1])
                ax.set_ylim(0,1)
                plt.fill_between(x,0, plot_data)
                ax.fill_between(x,0, plot_data)

                ind_filename = generate_output_filename(timestamp, target,
                    connected, left, 'percent_activity')


                xlabel = "Time (seconds)"
                ylabel = "% Active Neurons"

                plt.xlabel(xlabel, fontsize = fs)
                plt.ylabel(ylabel, fontsize = fs)

                #individual plot title
                title = target
                if(connected):
                    title = title + " (Connected)"
                else:
                    title = title + " (Not Connected)"
                if(left):
                    title = title + " (Left)"
                else:
                    title = title + " (Right)"

                plt.title(title, fontsize = fs)
                print(ind_filename)
                ind_fig.savefig(ind_filename)
                plt.close(ind_fig)

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

        f.suptitle(title, x = 0.6, fontsize = 15)


        output_dir = "output/figures/" + timestamp
        filename = target
        output_dir = output_dir + "/percent_activity/"
        filename = filename + "_percent_activity.pdf"

        output_dir = output_dir + "grouped/"
        os.makedirs(output_dir, exist_ok = True)

        output_filename = output_dir + filename

        print(output_filename)

        f.savefig(output_filename, pad_inches = 10)
        plt.close()


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

                output_filename = generate_output_filename(timestamp, target,
                        connected, left, 'heatmap')

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
    fs = 12

    for target in targets:

        #find axis limits
        #don't think it's possible to calculate on the fly
        target_min = 100000
        target_max = -100000

        for connected in connections:
            for left in left_options:
                stem = get_filename_stem(target, timestamp, connected, left)
                volt_filename = stem + "volt.npy"
                volt = np.load(volt_filename)
                mi = volt.min()
                ma = volt.max()
                if mi < target_min:
                    target_min = mi
                if ma > target_max:
                    target_max = ma


        nrows = len(connections)
        ncols = len(left_options)

        f, axarr = plt.subplots(2,2)

        subplot_counter = 1
        i = 0
        for connected in connections:
            j = 0

            for left in left_options:

                ind_fig = plt.figure(figsize = (7,5))
                ax = axarr[i,j]
                stem = get_filename_stem(target, timestamp, connected, left)
                time_filename = stem + "time.npy"
                volt_filename = stem + "volt.npy"

                time = np.load(time_filename)
                volt = np.load(volt_filename)

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
                    plt.plot(time, values, linewidth = 1)

                    ind_filename = generate_output_filename(timestamp, target,
                        connected, left, 'average_voltage_trace')

                else:
                    for x in range(volt.shape[0]):
                        ax.plot(time, volt[x], linewidth = 0.5, color = "#3333cc",
                                alpha = 0.5)

                        plt.plot(time, volt[x], linewidth = 0.5)

                    ind_filename = generate_output_filename(timestamp, target,
                        connected, left, 'voltage_trace')


                xlabel = "Time (seconds)"
                ylabel = "Voltage (volts)"

                plt.xlabel(xlabel, fontsize = fs)
                plt.ylabel(ylabel, fontsize = fs)

                #individual plot title
                title = target
                if(connected):
                    title = title + " (Connected)"
                else:
                    title = title + " (Not Connected)"
                if(left):
                    title = title + " (Left)"
                else:
                    title = title + " (Right)"
                if(average):
                    title = title + " Average across %d neurons" % volt.shape[0]

                plt.title(title, fontsize = fs)
                print(ind_filename)
                ind_fig.savefig(ind_filename)
                plt.close(ind_fig)

                ax.set(xlabel = xlabel, ylabel = ylabel)

                buff = (target_max - target_min) * 0.1
                ax.set_ylim(target_min - buff, target_max + buff)

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


        output_dir = "output/figures/" + timestamp
        filename = target
        if(average):
            output_dir = output_dir + "/average_voltage_trace/"
            filename = filename + "_average_voltage_trace.pdf"
        else:
            output_dir = output_dir + "/voltage_trace/"
            filename = filename + "_voltage_trace.pdf"

        output_dir = output_dir + "grouped/"
        os.makedirs(output_dir, exist_ok = True)

        output_filename = output_dir + filename

        print(output_filename)

        f.savefig(output_filename, pad_inches = 10)
        plt.close()

def percentage_overlayed(targets, connections, left_options, timestamp):

    plt.rcParams.update({'font.size':6})
    fs = 12

    for target in targets:

        f, axarr = plt.subplots(2,1)

        subplot_counter = 1
        i = 0
        for connected in connections:

            for left in left_options:

                if(left):
                    label = "Left"
                else:
                    label = "Right"

                ind_fig = plt.figure(figsize = (7,5))
                ax = axarr[i]
                stem = get_filename_stem(target, timestamp, connected, left)
                volt_filename = stem + "volt.npy"

                volt = np.load(volt_filename)

                '''
                ax.plot(time, volt[x], linewidth = 0.5, color = "#3333cc",
                        alpha = 0.5)

                plt.plot(time, volt[x], linewidth = 0.5)
                '''


                plot_data = []

                n_neurons = volt.shape[0]
                total_time = volt.shape[1]

                for time in range(total_time):
                    s = 0
                    for neuron in range(n_neurons):
                        if volt[neuron, time] > -0.051:
                            s = s + 1
                    plot_data.append(float(s) / n_neurons)


                x = range(len(plot_data))
                ax.plot(x,plot_data, label = label)
                ax.set_ylim(-0.1,1.1)
                ax.fill_between(x,0, plot_data, alpha = 0.5)



                xlabel = "Time (milliseconds)"
                ylabel = "% Active Neurons"

                ax.set(xlabel = xlabel, ylabel = ylabel)

                if(connected):
                    label = "Connected"
                else:
                    label = "Not Connected"
                ax.annotate(label, xy=(-0.25,0.5), xycoords=("axes fraction",
                    "axes fraction"), weight = "bold")

                ax.legend(loc = 1)

            i = i + 1

        f.tight_layout(rect=[0.15,0,1,0.9])

        title = target

        f.suptitle(title, x = 0.6, fontsize = 15)

        output_dir = "output/figures/" + timestamp
        filename = target
        output_dir = output_dir + "/percent_activity/"
        filename = filename + "_percent_activity.pdf"

        output_dir = output_dir + "grouped/"
        os.makedirs(output_dir, exist_ok = True)

        output_filename = output_dir + filename

        print(output_filename)

        f.savefig(output_filename, pad_inches = 10)
        plt.close(f)



if __name__ == "__main__":
    main()




