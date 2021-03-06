#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

def main():

    parent_dir = "2s_3"

    comparison = ["unconnected", "physiological", "fully_connected"]
    targets = ['PY', 'HTC', 'RTC', 'FS', 'IN', 'RE', 'TC']
    left_options = [True,False]

    voltage_traces_overlayed(targets, comparison, left_options, parent_dir,
            average = True)
    heatmaps(targets, comparison, left_options, parent_dir)
    voltage_traces(targets, comparison, left_options, parent_dir)
    voltage_traces(targets, comparison, left_options, parent_dir, average = True)
    percentage(targets, comparison, left_options, parent_dir)
    percentage_overlayed(targets, comparison, left_options, parent_dir)
    voltage_traces_overlayed(targets, comparison, left_options, parent_dir)

#name directory structure and filename redundantly in case they get separated
#also makes directories required to write files
def generate_output_filename(parent_dir, target, left, choice, plot_type):

    output_dir = "output/figures/" + parent_dir + '/' + plot_type + '/'
    filename = target

    output_dir = output_dir + choice +"/"
    filename = filename + "_" + choice

    if(left):
        output_dir = output_dir + "left/"
        filename = filename + "_left_"
    else:
        output_dir = output_dir + "right/"
        filename = filename + "_right_"

    filename = filename + plot_type + ".pdf"

    os.makedirs(output_dir, exist_ok = True)

    return output_dir + filename

def dirname_to_title(dirname):

    s = dirname.replace("_", " ")
    s = s.title()
    return s

def percentage(targets, comparison, left_options, parent_dir):

    plt.rcParams.update({'font.size':6})
    fs = 12

    for target in targets:
        f, axarr = plt.subplots(len(comparison),2)

        subplot_counter = 1
        i = 0
        for choice in comparison:
            j = 0

            for left in left_options:

                ind_fig = plt.figure(figsize = (7,5))
                ax = axarr[i,j]
                stem = get_filename_stem(target, parent_dir, choice, left)
                volt_filename = stem + "volt.npy"
                time_filename = stem + "time.npy"

                voltage_data = np.load(volt_filename)
                time_data = np.load(time_filename)

                plot_data = []

                n_neurons = voltage_data.shape[0]
                total_time = voltage_data.shape[1]

                for time in range(total_time):
                    s = 0
                    for neuron in range(n_neurons):

                        #voltage threshold is completely arbitrary at the moment
                        if voltage_data[neuron, time] > -0.051:
                            s = s + 1

                    plot_data.append(float(s) / n_neurons)

                x = range(len(plot_data))


                #similar operations for individual plot and grouped plots
                plt.plot(time_data,plot_data)
                ax.plot(time_data,plot_data)

                plt.ylim([0,1])
                ax.set_ylim(0,1)

                plt.fill_between(time_data,0, plot_data)
                ax.fill_between(time_data,0, plot_data)

                ind_filename = generate_output_filename(parent_dir, target,
                    left, choice, 'percent_activity')

                xlabel = "Time (seconds)"
                ylabel = "% Active Neurons"

                plt.xlabel(xlabel, fontsize = fs)
                plt.ylabel(ylabel, fontsize = fs)

                #individual plot title
                title = target
                title = title + " (%s)" % dirname_to_title(choice) 
                if(left):
                    title = title + " (Left)"
                else:
                    title = title + " (Right)"

                plt.title(title, fontsize = fs)
                print(ind_filename)
                ind_fig.savefig(ind_filename)
                plt.close(ind_fig)

                ax.set(xlabel = xlabel, ylabel = ylabel)


                #outer y-axis labels
                if(j == 0):
                    label = dirname_to_title(choice) 

                    ax.annotate(label, xy=(-0.65,0.5), xycoords=("axes fraction",
                        "axes fraction"), weight = "bold")

                #outer x-axis labels
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

        #fragile, keep upright
        f.tight_layout(rect=[0.15,0,1,0.9])

        title = target

        f.suptitle(title, x = 0.6, fontsize = 15)

        output_dir = "output/figures/" + parent_dir
        filename = target
        output_dir = output_dir + "/percent_activity/"
        filename = filename + "_grouped_percent_activity.pdf"

        output_dir = output_dir + "grouped/"
        os.makedirs(output_dir, exist_ok = True)

        output_filename = output_dir + filename

        print(output_filename)

        f.savefig(output_filename, pad_inches = 10)
        plt.close()


def heatmaps(targets, comparison, left_options, parent_dir):

    for target in targets:
        for choice in comparison:
            for left in left_options:

                fig = plt.figure(figsize = (11, 8.5))

                stem = get_filename_stem(target, parent_dir, choice, left)
                volt_filename = stem + "volt.npy"
                time_filename = stem + "time.npy"

                voltage_data = np.load(volt_filename)
                time_data = np.load(time_filename)

                #normalize plot height
                num_neurons = voltage_data.shape[0]
                magic_size = 20000
                aspect = magic_size / num_neurons

                #determine maximum magnitude for centering colormap around 0
                m = 0
                for x in voltage_data.flat:
                    if abs(x) > m:
                        m = abs(x)

                plt.imshow(voltage_data, aspect = aspect, cmap = 'seismic', vmin = -m, vmax = m)
                plt.colorbar(shrink = 0.35)
                plt.xlabel('Time (seconds)')
                plt.ylabel('Neuron index')

                title = "%s (%s)" % (target, dirname_to_title(choice))
                plt.title(title)

                #shifty indexing to get last time value to be labeled
                pos = np.arange(0,voltage_data.shape[1], 9999)
                lab = []
                for val in pos:
                    lab.append("%.1f" % time_data[val])


                plt.xticks(pos, lab)

                output_filename = generate_output_filename(parent_dir, target,
                        left, choice, 'heatmap')

                print(output_filename)

                fig.savefig(output_filename)
                plt.close()






def get_filename_stem(target, parent_dir, choice, is_left):

    stem = "output/" + parent_dir + "/%s/" % choice

    if(is_left):
        name = target + "L_"
    else:
        name = target + "R_"

    return stem + name

def voltage_traces(targets, comparison, left_options, parent_dir, average = False):

    plt.rcParams.update({'font.size':6})
    fs = 12

    for target in targets:

        #find axis limits
        #don't think it's possible to calculate on the fly
        target_min = 100000
        target_max = -100000

        for choice in comparison:
            for left in left_options:
                stem = get_filename_stem(target, parent_dir, choice, left)
                volt_filename = stem + "volt.npy"
                volt = np.load(volt_filename)
                mi = volt.min()
                ma = volt.max()
                if mi < target_min:
                    target_min = mi
                if ma > target_max:
                    target_max = ma

        f, axarr = plt.subplots(len(comparison),2)

        subplot_counter = 1
        i = 0
        for choice in comparison:
            j = 0

            for left in left_options:

                ind_fig = plt.figure(figsize = (7,5))
                ax = axarr[i,j]
                stem = get_filename_stem(target, parent_dir, choice, left)
                time_filename = stem + "time.npy"
                volt_filename = stem + "volt.npy"

                time = np.load(time_filename)
                volt = np.load(volt_filename)


                #ax is in group, plt is individual
                #you can't say 'fig.plt()', but I don't think this is proper

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

                    ind_filename = generate_output_filename(parent_dir, target,
                        left, choice, 'average_voltage_trace')

                else:
                    for x in range(volt.shape[0]):
                        ax.plot(time, volt[x], linewidth = 0.5, color = "#3333cc",
                                alpha = 0.5)

                        plt.plot(time, volt[x], linewidth = 0.5, alpha = 0.5)

                    ind_filename = generate_output_filename(parent_dir, target,
                        left, choice, 'voltage_trace')


                xlabel = "Time (seconds)"
                ylabel = "Voltage (volts)"

                plt.xlabel(xlabel, fontsize = fs)
                plt.ylabel(ylabel, fontsize = fs)

                #individual plot title
                title = target

                title = title + " (%s)" % dirname_to_title(choice)

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
                    label = dirname_to_title(choice)

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

        #fragile, handle with care
        f.tight_layout(rect=[0.15,0,1,0.9])

        title = target

        if(average):
            title = title + " (Average across %d neurons)" % volt.shape[0]

        f.suptitle(title, x = 0.6, fontsize = 15)


        output_dir = "output/figures/" + parent_dir
        filename = target
        if(average):
            output_dir = output_dir + "/average_voltage_trace/"
            filename = filename + "_grouped_average_voltage_trace.pdf"
        else:
            output_dir = output_dir + "/voltage_trace/"
            filename = filename + "_grouped_voltage_trace.pdf"

        output_dir = output_dir + "grouped/"
        os.makedirs(output_dir, exist_ok = True)

        output_filename = output_dir + filename

        print(output_filename)

        f.savefig(output_filename, pad_inches = 10)
        plt.close()

def percentage_overlayed(targets, comparison, left_options, parent_dir):

    plt.rcParams.update({'font.size':6})
    fs = 12

    for target in targets:

        f, axarr = plt.subplots(len(comparison),1)

        subplot_counter = 1
        i = 0
        for choice in comparison:

            for left in left_options:

                if(left):
                    label = "Left"
                else:
                    label = "Right"

                ax = axarr[i]
                stem = get_filename_stem(target, parent_dir, choice, left)
                volt_filename = stem + "volt.npy"
                time_filename = stem + "time.npy"

                time_data = np.load(time_filename)
                voltage_data = np.load(volt_filename)

                plot_data = []

                n_neurons = voltage_data.shape[0]
                total_samples = voltage_data.shape[1]

                for time in range(total_samples):
                    s = 0
                    for neuron in range(n_neurons):

                        #threshold is currently arbitrary
                        if voltage_data[neuron, time] > -0.051:
                            s = s + 1
                    plot_data.append(float(s) / n_neurons)


                x = range(len(plot_data))
                ax.plot(time_data,plot_data, label = label, linewidth = 0.5)
                ax.set_ylim(-0.1,1.1)
                ax.fill_between(time_data,0, plot_data, alpha = 0.5)



                xlabel = "Time (seconds)"
                ylabel = "% Active Neurons"

                ax.set(xlabel = xlabel, ylabel = ylabel)

                label = dirname_to_title(choice)

                ax.annotate(label, xy=(-0.25,0.5), xycoords=("axes fraction",
                    "axes fraction"), weight = "bold")

                ax.legend(loc = 1)

            i = i + 1

        #fragile, be gentle
        f.tight_layout(rect=[0.15,0,1,0.9])

        title = target

        f.suptitle(title, x = 0.6, fontsize = 15)

        output_dir = "output/figures/" + parent_dir
        filename = target
        output_dir = output_dir + "/percent_activity/"
        filename = filename + "_overlayed_percent_activity.pdf"

        output_dir = output_dir + "overlayed/"
        os.makedirs(output_dir, exist_ok = True)

        output_filename = output_dir + filename

        print(output_filename)

        f.savefig(output_filename, pad_inches = 10)
        plt.close(f)


def voltage_traces_overlayed(targets, comparison, left_options, parent_dir, average = False):

    plt.rcParams.update({'font.size':6})
    fs = 12

    for target in targets:

        f, axarr = plt.subplots(len(comparison),1)

        subplot_counter = 1
        i = 0

        for choice in comparison:

            #find axis limits
            #don't think it's possible to calculate on the fly
            target_min = 100000
            target_max = -100000

            for left in left_options:
                stem = get_filename_stem(target, parent_dir, choice, left)
                volt_filename = stem + "volt.npy"
                volt = np.load(volt_filename)
                volt = volt[:,5000:]

                n_neurons = volt.shape[0]
                time_steps = volt.shape[1]

                averaged = np.empty(shape = (time_steps))
                for m in range(time_steps):
                    s = 0
                    for n in range(n_neurons):
                        s = s + volt[n,m]
                    averaged[m] = s / n_neurons

                mi = averaged.min()
                ma = averaged.max()

                if mi < target_min:
                    target_min = mi
                if ma > target_max:
                    target_max = ma

            if(average):
                for_cross_corr = []

            for left in left_options:

                if(left):
                    col = 'red'
                    label = "Left"
                else:
                    label = "Right"
                    col = 'blue'

                ax = axarr[i]
                stem = get_filename_stem(target, parent_dir, choice, left)
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


                    ax.plot(time, values, linewidth = 0.2, color = col, label =
                            label)

                    for_cross_corr.append(values)

                else:
                    for x in range(volt.shape[0]):

                        #only label first line drawn to avoid redundant labels
                        #surely a better option exists
                        if x == 0:
                            ax.plot(time, volt[x], linewidth = 0.5,
                                    alpha = 0.5, color = col, label = label)
                        else:
                            ax.plot(time, volt[x], linewidth = 0.5,
                                    alpha = 0.5, color = col)

                xlabel = "Time (seconds)"
                ylabel = "Voltage (volts)"

                ax.set(xlabel = xlabel, ylabel = ylabel)

                buff = (target_max - target_min) * 0.1
                ax.set_ylim(target_min - buff, target_max + buff)

                ax.legend(loc = 1)

                label = dirname_to_title(choice)

                ax.annotate(label, xy=(-0.3,0.5), xycoords=("axes fraction",
                    "axes fraction"), weight = "bold")

                subplot_counter = subplot_counter + 1

            i = i + 1

            if(average):
                a = for_cross_corr[0]
                b = for_cross_corr[1]

                za = (a - a.mean()) / a.std()
                zb = (b - b.mean()) / b.std()

                cc = np.correlate(za,zb) / a.shape[0]

                ax.annotate("CC: %.2f" % cc, xy = (-0.3,0.4), xycoords = 
                    ("axes fraction", "axes fraction"))

        #still fragile
        f.tight_layout(rect=[0.15,0,1,0.9])

        title = target

        if(average):
            title = title + " (Average across %d neurons)" % volt.shape[0]

        f.suptitle(title, x = 0.6, fontsize = 15)

        output_dir = "output/figures/" + parent_dir
        filename = target
        if(average):
            output_dir = output_dir + "/average_voltage_trace/"
            filename = filename + "_overlayed_average_voltage_trace.pdf"
        else:
            output_dir = output_dir + "/voltage_trace/"
            filename = filename + "_overlayed_voltage_trace.pdf"

        output_dir = output_dir + "overlayed/"
        os.makedirs(output_dir, exist_ok = True)

        output_filename = output_dir + filename

        print(output_filename)

        f.savefig(output_filename, pad_inches = 10)
        plt.close()

if __name__ == "__main__":
    main()




