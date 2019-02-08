#!/usr/bin/env python
from brian2 import *
import time
import os
import argparse

def main():

    timestamp = str(int(time.time()))

    num_runs = 1

    prefs.codegen.target = 'numpy'
    BrianLogger.log_level_debug()

    duration = 1*second

    frequencies = [0,5,10,15,20]

    for frequency in frequencies:

        for i in range(num_runs):

            print("RUN: %d" % i)

            run_name = str(frequency) + "hz"

            cortex_unconnected(run_name, i, duration, frequency)


def cortex_unconnected(run_name, i, duration, frequency):

    name = "cortex_cortex_left_sin"

    stem = "output/" + name + "/" + run_name + "/"
    print(stem)

    #Time constants
    t_step = 0.02*ms

    t_d = t_step*numpy.arange(duration/t_step)/second
    Sgn = TimedArray(sin(frequency*2*pi*t_d)*50*pA, dt = t_step)

    #Number of Each Cell Type(sqrt of # Thalamus Cells)
    nPY = 80
    nFS = 20

    #Main Izhikevich Equation, PY FS
    shared_PY_FS = '''
    Cm                                                                                          :farad
    k                                                                                           :siemens/volt
    vt                                                                                          :volt
    vr                                                                                          :volt
    du/dt = a*(b*((v - vr)/mV)**un*(v > vb) - u)                                                :amp
    a                                                                                           :Hz
    b                                                                                           :amp
    vb                                                                                          :volt
    un                                                                                          :1
    vth                                                                                         :volt
    c                                                                                           :volt
    d                                                                                           :amp
    sg                                                                                          :1
    iStm                                                                                        :amp
    '''

    right_PY_FS = '''
    dv/dt = (k*(v - vr)*(v - vt) - u - iSyn + iStm)/Cm                                          :volt
    iPYR                                                                                        :amp
    iFSR                                                                                        :amp
    iSyn = iPYR + iFSR                                                                          :amp
    LFP = abs(iPYR + iFSR)                                                                      :amp
    ''' + shared_PY_FS

    left_PY_FS = '''
    dv/dt = (k*(v - vr)*(v - vt) - u - iSyn + iStm + Sgn(t))/Cm                                 :volt
    iPYL                                                                                        :amp
    iFSL                                                                                        :amp
    iSyn = iPYL + iFSL                                                                          :amp
    LFP = abs(iPYL + iFSL)                                                                      :amp
    ''' + shared_PY_FS


    FSRg = NeuronGroup(nFS, right_PY_FS, method = 'rk4', dt = t_step, threshold = 'v>= vth', reset = '''
    v = c
    u += d
    ''', name = "FSRg")

    PYRg = NeuronGroup(nPY, right_PY_FS, method = 'rk4', dt = t_step, threshold = 'v>= vth', reset = '''
    v = c
    u += d
    ''', name = "PYRg")

    FSLg = NeuronGroup(nFS, left_PY_FS, method = 'rk4', dt = t_step, threshold = 'v>= vth', reset = '''
    v = c
    u += d
    ''', name = "FSLg")

    PYLg = NeuronGroup(nPY, left_PY_FS, method = 'rk4', dt = t_step, threshold = 'v>= vth', reset = '''
    v = c
    u += d
    ''', name = "PYLg")

    #PYR Group Declaration
    PYRg.Cm   = PYLg.Cm   = 100*pF
    PYRg.k    = PYLg.k    = 0.7*pA/(mvolt*mvolt)
    PYRg.vr   = PYLg.vr   = -60*mV
    PYRg.vt   = PYLg.vt   = -40*mV
    PYRg.a    = PYLg.a    =  0.03*kHz
    PYRg.b    = PYLg.b    = -2*pA
    PYRg.vb   = PYLg.vb   = -200*volt
    PYRg.un   = PYLg.un   = 1
    PYRg.vth  = PYLg.vth  = 35*mV
    PYRg.c    = PYLg.c    = -50*mV
    PYRg.d    = PYLg.d    = 100*pA
    PYRg.v    = PYLg.v    = -60*mV
    PYRg.iStm = PYLg.iStm = 79*pA
    PYRg.sg   = PYLg.sg   = 1

    PYRg.u = PYRg.b*(PYRg.v-PYRg.vr)/mV
    PYLg.u = PYLg.b*(PYLg.v-PYLg.vr)/mV

    #FS Group Declaration
    FSRg.Cm   = FSLg.Cm   = 20*pF
    FSRg.k    = FSLg.k    = 1*pA/(mvolt*mvolt)
    FSRg.vr   = FSLg.vr   = -55*mV
    FSRg.vt   = FSLg.vt   = -40*mV
    FSRg.a    = FSLg.a    = 0.2*kHz
    FSRg.b    = FSLg.b    = 0.025*pA
    FSRg.un   = FSLg.un   = 3
    FSRg.vth  = FSLg.vth  = 25*mV
    FSRg.c    = FSLg.c    = -45*mV
    FSRg.d    = FSLg.d    = 0*pA
    FSRg.v    = FSLg.v    = -55*mV
    FSRg.u    = FSLg.u    = 0*pA
    FSRg.iStm = FSLg.iStm = 60*pA

    FSRg.vb = FSRg.vr
    FSLg.vb = FSLg.vr

    #Chemical Synapse, GABA_A, AMPA
    CSa = '''
    D_i                                                                                         :1
    t_spike                                                                                     :second
    g_syn                                                                                       :siemens
    e_syn                                                                                       :volt
    alpha                                                                                       :Hz
    beta                                                                                        :Hz
    T = 0.5*(t-t_spike < 2.3*ms)*(t-t_spike > 2*ms)*(t_spike > 0*ms)                            :1
    D = 1 - (1 - D_i*(1 - 0.07))*exp(-(t - t_spike)/700/ms)                                     :1
    ds/dt = T*alpha*(1-s) - beta*s                                                              :1 (clock-driven)
    iCS = D*g_syn*s*(v_post - e_syn)                                                            :amp
    '''

    #Chemical Synapse, NMDA
    CSb = '''
    D_i                                                                                         :1
    t_spike                                                                                     :second
    g_syn                                                                                       :siemens
    alpha                                                                                       :Hz
    beta                                                                                        :Hz
    T = 0.5*(t-t_spike < 2.3*ms)*(t-t_spike > 2*ms)*(t_spike > 0*ms)                            :1
    ds/dt = T*alpha*(1-s) - beta*s                                                              :1 (clock-driven)
    D = 1 - (1 - D_i*(1 - 0.07))*exp(-(t - t_spike)/700/ms)                                     :1
    B = 1/(1+exp(-(v_post + 25*mV)/12.5/mV))                                                    :1
    iCS = D*B*g_syn*s*v_post                                                                    :amp
    '''

    preCS = '''
    t_spike = t
    D_i = D
    '''

    iPYRa_eqs = '''
    iPYRa_post = iCS                                                                             :amp (summed)
    '''

    iPYLa_eqs = '''
    iPYLa_post = iCS                                                                             :amp (summed)
    '''

    iPYRb_eqs = '''
    iPYRb_post = iCS                                                                             :amp (summed)
    '''

    iPYLb_eqs = '''
    iPYLb_post = iCS                                                                             :amp (summed)
    '''

    #Modular Synapse Equations, Izhikevich
    #Chemical Synapse, GABA_A, AMDA
    PYFS_cs = '''
    tau                                                                                         :second
    e_syn                                                                                       :volt
    g_spike                                                                                     :siemens
    dg/dt = -g/tau                                                                              :siemens (clock-driven)
    iCS = g*(v_post - e_syn)                                                                    :amp
    '''

    prePYFS = '''
    g += g_spike
    '''

    iPYR_eqs = '''
    iPYR_post = iCS                                                                             :amp (summed)
    '''

    iPYL_eqs = '''
    iPYL_post = iCS                                                                             :amp (summed)
    '''

    iFSR_eqs = '''
    iFSR_post = iCS                                                                             :amp (summed)
    '''

    iFSL_eqs = '''
    iFSL_post = iCS                                                                             :amp (summed)
    '''

    PYR_PYR_cs = Synapses(PYRg,PYRg,(PYFS_cs + iPYR_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "PYR_PYR_cs")
    PYR_FSR_cs = Synapses(FSRg,PYRg,(PYFS_cs + iFSR_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "PYR_FSR_cs")

    PYL_PYL_cs = Synapses(PYLg,PYLg,(PYFS_cs + iPYL_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "PYL_PYL_cs")
    PYL_FSL_cs = Synapses(FSLg,PYLg,(PYFS_cs + iFSL_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "PYL_FSL_cs")

    FSR_PYR_cs = Synapses(PYRg,FSRg,(PYFS_cs + iPYR_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "FSR_PYR_cs")
    FSR_FSR_cs = Synapses(FSRg,FSRg,(PYFS_cs + iFSR_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "FSR_FSR_cs")

    FSL_PYL_cs = Synapses(PYLg,FSLg,(PYFS_cs + iPYL_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "FSL_PYL_cs")
    FSL_FSL_cs = Synapses(FSLg,FSLg,(PYFS_cs + iFSL_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "FSL_FSL_cs")

    #ipsilateral connection probabilities
    PY_PY_p = 0.5
    PY_FS_p = 0.8
    FS_FS_p = 0.8

    PYR_PYR_cs.connect(p = '0.5')
    PYR_FSR_cs.connect('abs(i*(nPY-1.0)/(nFS-1.0) - j) < nPY*0.2', p = '0.8')

    PYL_PYL_cs.connect(p = '0.5')
    PYL_FSL_cs.connect('abs(i*(nPY-1.0)/(nFS-1.0) - j) < nPY*0.2', p = '0.8')

    FSR_PYR_cs.connect(i = PYR_FSR_cs.j[:], j = PYR_FSR_cs.i[:])
    FSR_FSR_cs.connect('abs(i-j) < nFS/2.0*(i!=j)', p = '0.8')

    FSL_PYL_cs.connect(i = PYL_FSL_cs.j[:], j = PYL_FSL_cs.i[:])
    FSL_FSL_cs.connect('abs(i-j) < nFS/2.0*(i!=j)', p = '0.8')

    PYR_PYR_cs.tau = 2*ms
    PYR_PYR_cs.e_syn = 0*mV
    PYR_PYR_cs.g_spike = 0.3*nS
    PYR_FSR_cs.tau = 10*ms
    PYR_FSR_cs.e_syn = -70*mV
    PYR_FSR_cs.g_spike = 0.3*nS

    PYL_PYL_cs.tau = 2*ms
    PYL_PYL_cs.e_syn = 0*mV
    PYL_PYL_cs.g_spike = 0.3*nS
    PYL_FSL_cs.tau = 10*ms
    PYL_FSL_cs.e_syn = -70*mV
    PYL_FSL_cs.g_spike = 0.3*nS

    FSR_PYR_cs.tau = 2*ms
    FSR_PYR_cs.e_syn = 0*mV
    FSR_PYR_cs.g_spike = 0.4*nS
    FSR_FSR_cs.tau = 10*ms
    FSR_FSR_cs.e_syn = -70*mV
    FSR_FSR_cs.g_spike = 0.03*nS

    FSL_PYL_cs.tau = 2*ms
    FSL_PYL_cs.e_syn = 0*mV
    FSL_PYL_cs.g_spike = 0.4*nS
    FSL_FSL_cs.tau = 10*ms
    FSL_FSL_cs.e_syn = -70*mV
    FSL_FSL_cs.g_spike = 0.03*nS

    PYR_volt = StateMonitor(PYRg, 'v', record = True)
    PYL_volt = StateMonitor(PYLg, 'v', record = True)
    FSL_volt = StateMonitor(FSLg,'v', record = True)
    FSR_volt = StateMonitor(FSRg,'v', record = True)

    print("Setup complete.")
    run(duration)


    os.makedirs(stem, exist_ok = True)
    stem = stem + str(i) + "/"
    print(stem)

    os.makedirs(stem, exist_ok = True)

    np.save(stem + "PYR_time.npy", PYR_volt.t)
    np.save(stem + "PYR_volt.npy", PYR_volt.v)
    np.save(stem + "PYL_time.npy", PYL_volt.t)
    np.save(stem + "PYL_volt.npy", PYL_volt.v)

    np.save(stem + "FSR_time.npy", FSR_volt.t)
    np.save(stem + "FSR_volt.npy", FSR_volt.v)
    np.save(stem + "FSL_time.npy", FSL_volt.t)
    np.save(stem + "FSL_volt.npy", FSL_volt.v)

    settings_file = stem + "settings.txt"
    s = open(settings_file, "w+")
    s.write("duration: %s\n" % duration)

    print("Output written to directory: %s" % run_name)
if __name__ == "__main__":
    main()
