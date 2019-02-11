#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import signal
import os
import glob
import math
import scipy.io

import random

def main():

    parent_dir = "diff_stim_both_connected"

    frequencies = ["0hz", "10hz", "20hz"]
    comparison = []

    for freq_a in frequencies:
        for freq_b in frequencies:
            comp = freq_a + "_" + freq_b
            comparison.append(comp)

    targets = ['PY', 'FS']
    left_options = [True,False]

    file_dic = get_file_dic(targets, parent_dir, comparison,
            left_options)

    trace_spd_coherence(file_dic, targets, comparison, left_options, parent_dir)
    rasters(file_dic, targets, comparison, left_options, parent_dir)
    #voltage_traces(targets, comparison, left_options, parent_dir, average = True)
    '''
    voltage_traces_overlayed(targets, comparison, left_options, parent_dir)
    percentage_overlayed(targets, comparison, left_options, parent_dir)
    percentage(targets, comparison, left_options, parent_dir)
    voltage_traces_overlayed(targets, comparison, left_options, parent_dir)
    heatmaps(targets, comparison, left_options, parent_dir) ''' 
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
                title = title + " (%s)" % dirname_totitle(choice) 
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

    stems = glob.glob("output/" + parent_dir + "/" + choice + "/*/")

    if(is_left):
        name = target + "L_"
    else:
        name = target + "R_"

    stems = [stem + name for stem in stems]

    return stems

def get_file_key(target, parent_dir, choice, is_left):

    s = target + parent_dir + choice + str(is_left)
    return s

def plot_spd(ax, volt, left):

    #use copy for trimming so originals aren't altered
    volt = np.copy(volt)

    #trim out first half for SPD
    length = volt.size
    volt = volt[int(length * 0.5):]

    if(left):
        color = "red"
        label = "Left"
    else:
        color = "blue"
        label = "Right"

    volt = volt - volt.mean()
    ps = np.abs(np.fft.fft(volt))**2

    time_step = .02 * .001

    freqs = np.fft.fftfreq(volt.size, time_step)
    idx = np.argsort(freqs)

    ax.plot(freqs[idx], ps[idx], color = color, label = label)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude")
    ax.set_title("SPD (only from second half)")
    ax.legend()
    ax.set_xlim((0,50))
    #ax.set_ylim((0,2.5 * (10**10)))

def plot_raster(ax, volt, time, left):

    cutoff = volt.mean() + (1.5 * volt.std())

    raster_data = volt > cutoff
    raster_data = raster_data.astype(int)

    if(left):
        cmap = "Reds"
        title = "Left"
    else:
        cmap = "Blues"
        title = "Right"
    ax.imshow(raster_data, aspect = 1000, cmap = cmap)

    #shifty indexing to get last time value to be labeled
    step = volt.shape[1] // 5
    pos = np.arange(0,volt.shape[1], step - 1)
    lab = []
    for val in pos:
        lab.append("%.1f" % time[val])

    ax.set_title(title)
    ax.set_xticks(pos)
    ax.set_xticklabels(lab)
    ax.set_xlabel("Time (seconds)")
    ax.set_ylabel("Neuron Index")

def plot_trace(ax, volt, time, left):

    #use copy for trimming so originals aren't altered
    volt = np.copy(volt)
    time = np.copy(time)

    #trim out first 20%
    length = time.shape[0]
    time = time[int(length * 0.2):]
    volt = volt[int(length * 0.2):]

    if(left):
        color = "red"
        label = "Left"

    else:
        color = "blue"
        label = "Right"

    ax.plot(time, volt, color = color, label = label, linewidth = 0.5)
    xlabel = "Time (seconds)"
    ylabel = "Voltage (millivolts)"
    ax.legend(loc = 1)
    ax.set(xlabel = xlabel, ylabel = ylabel)
    ax.set_title("Voltage Trace")
    ax.set_ylim([-60,-20])

def plot_coherence(ax, x1, x2):
    
    '''
    x1 = np.linspace(0,100000, 10000000)
    x1 = np.sin(2 * np.pi * x1) + np.sin(5 * 2 * np.pi * x1)
    x2 = np.linspace(0,100000, 10000000)
    x2 = np.sin(2 * np.pi * x2) + np.sin(5 * 2 * np.pi * x2)

    noise = np.random.normal(0,0.1, x1.size)
    x1 = np.add(x1, noise) 
    noise = np.random.normal(0,0.1, x1.size)
    x2 = np.add(x2, noise) 


    (x,y) = signal.coherence(x1, x2, fs = 100)



    plt.plot(x1[0:100])
    plt.plot(x2[0:100])
    plt.figure()
    plt.semilogy(x,y**(1/2))
    plt.show()
    exit()

    ''' 
    x1 = x1.mean(axis = 0)
    x2 = x2.mean(axis = 0)
    plt.plot(x1)
    plt.plot(x2)

    plt.figure()


    '''
    x1 = signal.decimate(x1, 100)
    x2 = signal.decimate(x2, 100)
    '''
    print(x1.size)
    print(x1.size / 500)
    x1 = signal.resample(x1, int(x1.size / 500))
    x2 = signal.resample(x2, int(x2.size / 500))


    np.savetxt('/home/josh/tmp/signal_1_downsampled.dat', x1)
    np.savetxt('/home/josh/tmp/signal_2_downsampled.dat', x2)

    '''
    plt.plot(x1)
    plt.plot(x2)
    plt.show()
    '''


    #copy for trimming
    x1 = np.copy(x1)
    x2 = np.copy(x2)

    #trim out first 50%
    x1 = x1[x1.size // 2:]
    x2 = x2[x2.size // 2:]
    fs = 1 / (0.02 * 0.001)

    fs = fs / 500

    (x,y) = signal.coherence(x1, x2, fs = fs, detrend =
            "linear", nperseg = fs)

    #100hz is normal
    #take one value for every 500

    ax.semilogy(x,y, linewidth = 1)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power")
    ax.set_title("Coherence (only from second half)")
    #ax.set_xlim((0,100))

#read in traces and compute average of each neuron across all runs
def get_file_dic(targets, parent_dir, comparison, left_options):

    file_dic = {}

    num_runs = len(get_filename_stem(targets[0], parent_dir, comparison[0],
        left_options[0]))

    #read in all files and average
    for target in targets:
        for choice in comparison:
            for left in left_options:
                stems = get_filename_stem(target, parent_dir, choice, left)

                key = get_file_key(target, parent_dir, choice, left)

                #determine size of array from first file
                volt_filename = stems[0] + "volt.npy"
                volt = np.load(volt_filename)

                volt = volt * 1000

                length = volt.shape[1]
                n_neurons = volt.shape[0]


                per_neuron_array = np.zeros(shape = (length,))

                for i, stem in enumerate(stems):
                    volt_filename = stem + "volt.npy"
                    volt = np.load(volt_filename)

                    volt = volt * 1000
                    per_neuron_array = np.add(per_neuron_array, volt)

                per_neuron_array = per_neuron_array / num_runs

                file_dic[key] = per_neuron_array

                time_filename = stems[0] + "time.npy"
                time = np.load(time_filename)

    print("Done reading.")
    return file_dic

def trace_spd_coherence(file_dic, targets, comparison, left_options, parent_dir):

    num_runs = len(get_filename_stem(targets[0], parent_dir, comparison[0],
        left_options[0]))

    plt.rcParams.update({'font.size':10})

    for target in targets:

        fig = plt.figure(figsize = (20, 20))


        outer = fig.add_gridspec(len(comparison), 1, wspace = 0.1, hspace = 1)

        for i, choice in enumerate(comparison):

            inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec =
                    outer[i], wspace= 0.1, hspace = 0.1)

            ax = plt.Subplot(fig, inner[0])
            ax2 = plt.Subplot(fig, inner[1])

            fig.add_subplot(ax)
            fig.add_subplot(ax2)

            for left in left_options:
                stem = get_filename_stem(target, parent_dir, choice, left)
                time_filename = stem[0] + "time.npy"

                time = np.load(time_filename)
                key = get_file_key(target, parent_dir, choice, left)
                volt = file_dic[key]
                avg_volt = volt.mean(axis = 0)

                plot_trace(ax, avg_volt, time, left)
                plot_spd(ax2, avg_volt, left)

                choice_label = choice.split("_")
                new = []
                for item in choice_label:
                    if item == "0hz":
                        new.append("None")
                    else:
                        new.append(item)

                choice_label = "/".join(new)
                ax.annotate(choice_label, xy=(-0.2,0.5), xycoords=("axes fraction",
                    "axes fraction"), weight = "bold")

        title = target

        for filetype in [".pdf", ".png"]:
            filename = target
            output_dir = "output/figures/" + parent_dir

            output_dir = output_dir + "/trace_spd_coherence/"
            filename = filename + "_trace_spd_coherence" + filetype

            os.makedirs(output_dir, exist_ok = True)

            output_filename = output_dir + filename

            print(output_filename)

            titledict = {'fontweight':'bold',
                         'fontsize':20,
                         'horizontalalignment':'center'
                        }

            subtextdict = {'fontsize':10,
                           'horizontalalignment':'center'
                          }

            fig.text(s = "Cortex with symmetrical projections", x = 0.5, y =
                    0.93, fontdict = titledict)

            fig.text(s = "LRP's are same probability as intra-hemisphere connections", x = 0.5, y = 0.92, fontdict = subtextdict)
            fig.text(s = "Stimulation amplitude chosen arbitrarily to get visible effects", x = 0.5, y = 0.91, fontdict = subtextdict)
            fig.text(s = "\"10hz/20hz\" is 10hz left stimulation, 20hz right stimulation", x = 0.5, y = 0.90, fontdict = subtextdict)

            fig.savefig(output_filename, bbox_inches = "tight")

        plt.close()


def rasters(file_dic, targets, comparison, left_options, parent_dir):

    plt.rcParams.update({'font.size':10})

    for target in targets:

        fig = plt.figure(figsize = (15,10))

        outer = fig.add_gridspec(len(comparison),1, wspace = 0.2, hspace = 0.5)

        for i, choice in enumerate(comparison):

            inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec =
                    outer[i], wspace= 0.1, hspace = 0.05)

            for j, left in enumerate(left_options):

                stem = get_filename_stem(target, parent_dir, choice, left)[0]
                time_filename = stem + "time.npy"

                time = np.load(time_filename)
                volt = file_dic[get_file_key(target, parent_dir, choice, left)]
                avg_volt = volt.mean(axis = 0)

                ax = plt.Subplot(fig, inner[j])
                if(left):
                    ax.annotate(choice, xy=(-0.6,0.5), xycoords=("axes fraction",
                     "axes fraction"), weight = "bold")

                plot_raster(ax, volt, time, left)

                fig.add_subplot(ax)

        for filetype in [".pdf", ".png"]:
            filename = target

            output_dir = "output/figures/" + parent_dir

            output_dir = output_dir + "/raster/"
            filename = filename + "_raster" + filetype

            os.makedirs(output_dir, exist_ok = True)

            output_filename = output_dir + filename

            print(output_filename)

            fig.savefig(output_filename, bbox_inches = "tight")

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

                if(average):
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

                '''
                za = (a - a.mean()) / a.std()
                zb = (b - b.mean()) / b.std()

                cc = np.correlate(za,zb) / a.shape[0]

                ax.annotate("CC: %.2f" % cc, xy = (-0.3,0.3), xycoords = 
                    ("axes fraction", "axes fraction"))
                '''
                f, c = signal.coherence(a,b)
                plt.figure()
                plt.plot(f,c)
                plt.show()
                ax.annotate("C: %.2f" % c, xy = (-0.3,0.3), xycoords = 
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

