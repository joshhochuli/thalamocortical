#!/usr/bin/env python
from brian2 import *
import time
import os
def main():
    
    connected = True

    #seed(555)

    prefs.codegen.target = 'numpy'
    BrianLogger.log_level_debug()

    left = []
    right = []
    other = []

    #Global Variables Across Thalamus Cells
    #Capacitance
    C = 1*uF/cm**2

    #Leak Variables
    gL = 0.01*mS/cm**2          #average g_L value +- 0.0025, uniformly distributed
    e_KL = -90*mV

    #Ca Variables
    F = 96489*coulomb/mole
    w = 0.5*um
    z = 2
    ca_Rest = 5e-5*mM

    #iCa Variables
    R = 8.32441*joule/(mole*kelvin)
    T = 309.15*kelvin
    eca0 = R*T/(2*F)
    ca_0 = 2*mM

    #iCat Variables
    g_Cat = 2.1*mS/cm**2

    #iNa Variables
    g_Na = 90*mS/cm**2
    e_Na = 50*mV

    #iK Variables
    g_K = 10*mS/cm**2
    e_K = -90*mV                        #-80*mV in INaK, -90*mV in paper

    #iH Variables       HTC
    e_H = -43*mV

    #iAHP Variables
    e_AHP = -90*mV                      #e_K

    #iCan Variables
    e_Can = 10*mV

    #Number of Each Cell Type(sqrt of # Thalamus Cells)
    nHTC = 7
    nRTC = 12
    nIN = 8
    nRE = 10
    nPY = 80
    nFS = 20

    #Time constants
    duration = 0.1*second
    t_step = 0.02*ms

    #tACs Signal Array
    t_d = t_step*numpy.arange(duration/t_step)/second

    #Modular Equation Expressions
    #Main HH Equation, HTC, RTC
    main_TC = '''
    dv/dt=(-g_L*(v-e_L)-g_KL*(v-e_KL)-iCa-iCal-iCat-iNa-iK-iH-iAHP-iCan-iSyn+iStm)/C            :volt
    e_L                                                                                         :volt
    g_KL                                                                                        :siemens/meter**2
    g_L                                                                                         :siemens/meter**2
    '''

    #Main HH Equation, IN
    main_IN = '''
    dv/dt=(-g_L*(v-e_L)-g_KL*(v-e_KL)-iCa-iNa-iK-iH-iAHP-iCan-iSyn+iStm)/C                      :volt
    e_L                                                                                         :volt
    g_KL                                                                                        :siemens/meter**2
    g_L                                                                                         :siemens/meter**2
    '''

    #Main HH Equation, RE
    main_RE = '''
    dv/dt=(-g_L*(v-e_L)-g_KL*(v-e_KL)-iCa-iNa-iK-iAHP-iCan-iSyn+iStm)/C                         :volt
    e_L                                                                                         :volt
    g_KL                                                                                        :siemens/meter**2
    g_L                                                                                         :siemens/meter**2
    '''

    #Ca_TC, HTC, RTC
    Ca_TC ='''
    tau_Ca = 10*ms                                                                              :second
    drive = -(iCa+iCat+iCal)/(z*F*w)*(-(iCa+iCat+iCal)/(z*F*w)>=0*mM/second)                    :mM/second
    dca/dt = drive + (ca_Rest - ca)/tau_Ca                                                      :mM
    '''

    #Ca_RE_IN, IN, RE
    Ca_RE_IN ='''
    tau_Ca                                                                                      :second
    drive = -iCa/(z*F*w)*(-iCa/(z*F*w)>=0*mM/second)                                            :mM/second
    dca/dt = drive + (ca_Rest - ca)/tau_Ca                                                      :mM
    '''

    #iCa, current 0, T0_TC
    i_Ca = '''
    g_Ca                                                                                        :siemens/meter**2
    e_Ca = eca0*log(ca_0/ca)                                                                    :volt
    m0_inf = 1/(1+exp(-(v+34*mV)/(6.2*mV)))                                                     :1
    tau_m0 = 0.612*ms+1*ms/(exp(-(v+107*mV)/(16.7*mV))+exp((v-8.2*mV)/(18.2*mV)))               :second
    dm0/dt = 4.6*(m0_inf-m0)/tau_m0                                                             :1
    h0_inf = 1/(1+exp((v+58*mV)/(4*mV)))                                                        :1
    tau_h0a = (1*ms)*(exp((v+442*mV)/(66.6*mV)))*(v < -55*mV)                                   :second
    tau_h0b = 28*ms+(1*ms)*(exp(-(v-3*mV)/(10.5*mV)))*(v >= -55*mV)                             :second
    tau_h0 = tau_h0a + tau_h0b                                                                  :second
    dh0/dt = 3.7*(h0_inf-h0)/tau_h0                                                             :1
    iCa = g_Ca*m0**2*h0*(v-e_Ca)                                                                :amp/meter**2
    '''

    #iCa_RE, current 0
    i_Ca_RE = '''
    g_Ca = 1.3*mS/cm**2                                                                         :siemens/meter**2
    e_Ca = eca0*log(ca_0/ca)                                                                    :volt
    m0_inf = 1/(1+exp(-(v+55*mV)/(7.4*mV)))                                                     :1
    tau_m0 = 3*ms + 1*ms/(exp(-(v+105*mV)/(15*mV))+exp((v+30*mV)/(10*mV)))                      :second
    dm0/dt = 6.9*(m0_inf-m0)/tau_m0                                                             :1
    h0_inf = 1/(1+exp((v+83*mV)/(5*mV)))                                                        :1
    tau_h0 = 85*ms + 1*ms/(exp(-(v+410*mV)/(50*mV))+exp((v+51*mV)/(4*mV)))                      :second
    dh0/dt = 3.7*(h0_inf-h0)/tau_h0                                                             :1
    iCa = g_Ca*m0**2*h0*(v-e_Ca)                                                                :amp/meter**2
    '''

    #iCat, current 1, IT0_TC
    i_Cat = '''
    m1_inf = 1/(1+exp(-(v+62*mV)/(6.2*mV)))                                                     :1
    tau_m1 = 0.612*ms+1*ms/(exp(-(v+135*mV)/(16.7*mV))+exp((v+19.8*mV)/(18.2*mV)))              :second
    dm1/dt = 4.6*(m1_inf-m1)/tau_m1                                                             :1
    h1_inf = 1/(1+exp((v+86*mV)/(4*mV)))                                                        :1
    tau_h1a = (1*ms)*(exp((v+470*mV)/(66.6*mV)))*(v < -83*mV)                                   :second
    tau_h1b = 28*ms+(1*ms)*(exp(-(v+25*mV)/(10.5*mV)))*(v >= -83*mV)                            :second
    tau_h1 = tau_h1a + tau_h1b                                                                  :second
    dh1/dt = 3.7*(h1_inf-h1)/tau_h1                                                             :1
    iCat = g_Cat*m1**2*h1*(v-e_Ca)                                                              :amp/meter**2
    '''

    #iCal, current 2, ICaL_TC
    i_Cal = '''
    g_Cal                                                                                       :siemens/meter**2
    m2_inf = 1/(1+exp(-(v+10*mV)/(4*mV)))                                                       :1
    tau_m2 = 0.4*ms + 0.7*ms/(exp(-(v+5*mV)/(15*mV))+exp((v+5*mV)/(15*mV)))                     :second
    dm2/dt = 4.6*(m2_inf-m2)/tau_m2                                                             :1
    h2_inf = 1/(1+exp((v+25*mV)/(2*mV)))                                                        :1
    tau_h2 = 300*ms+(100*ms)/(exp(-(v+40*mV)/(9.5*mV))+exp((v+40*mV)/(9.5*mV)))                 :second
    dh2/dt = 3.7*(h2_inf-h2)/tau_h2                                                             :1
    iCal = g_Cal*m2**2*h2*(v-e_Ca)                                                              :amp/meter**2
    '''

    #iNa, current 3, INaK
    i_Na = '''
    v_SH                                                                                        :volt
    a_m3 = (0.32*ms**-1*mV**-1)*(v-v_SH-13*mV)/(1-exp(-(v-v_SH-13*mV)/(4*mV)))                  :Hz
    b_m3 = (-0.28*mV**-1*ms**-1)*(v-v_SH-40*mV)/(1-exp((v-v_SH-40*mV)/(5*mV)))                  :Hz
    dm3/dt = (a_m3*(1-m3)-b_m3*m3)                                                              :1
    a_h3 = (0.128*ms**-1)*exp(-(v-v_SH-17*mV)/(18*mV))                                          :Hz
    b_h3 = (4*ms**-1)/(1+exp(-(v-v_SH-40*mV)/(5*mV)))                                           :Hz
    dh3/dt = (a_h3*(1-h3)-b_h3*h3)                                                              :1
    iNa = g_Na*m3**3*h3*(v-e_Na)                                                                :amp/meter**2
    '''

    #iK, current 4, INaK
    i_K = '''
    tm                                                                                          :1
    a_m4 = tm*(0.032*ms**-1*mV**-1)*(v-v_SH-15*mV)/(1-exp(-(v-v_SH-15*mV)/(5*mV)))              :Hz
    b_m4 = tm*(0.5*ms**-1)*exp(-(v-v_SH-10*mV)/(40*mV))                                         :Hz
    dm4/dt = (a_m4*(1-m4)-b_m4*m4)                                                              :1
    iK = g_K*m4**4*(v-e_K)                                                                      :amp/meter**2
    '''

    #iH, current 5, Ih_TC
    i_H = '''
    g_H                                                                                         :siemens/meter**2
    m5_inf = 1/(1+exp((v+75*mV)/(5.5*mV)))                                                      :1
    tau_m5 = 1*ms/(exp(-0.086/mV*v-14.59)+exp(0.0701/mV*v-1.87))                                :second
    dm5/dt = (m5_inf-m5)/(tau_m5)                                                               :1
    iH = g_H*m5*(v-e_H)                                                                         :amp/meter**2
    '''

    #iAHP, current 6, Iahp2
    i_AHP = '''
    g_AHP                                                                                       :siemens/meter**2
    m6_inf = 48*ca**2/(48*ca**2+0.09*mM**2)                                                     :1
    tau_m6 = (1*ms*mM**2)/(48*ca**2+0.09*mM**2)                                                 :second
    dm6/dt = (m6_inf - m6)/(tau_m6)                                                             :1
    iAHP = g_AHP*m6*(v-e_AHP)                                                                   :amp/meter**2
    '''

    #iCan, current 7, Ican_TC
    i_Can = '''
    g_Can                                                                                       :siemens/meter**2
    mCa = ca/(0.2*mM+ca)                                                                        :1
    m7_inf = 1/(1+exp(-(v+43*mV)/(5.2*mV)))                                                     :1
    tau_m7 = 1.6*ms+(2.7*ms)/(exp(-(v+55*mV)/(15*mV))+exp((v+55*mV)/(15*mV)))                   :second
    dm7/dt = (m7_inf - m7)/(tau_m7)                                                             :1
    iCan = g_Can*mCa*m7*(v-e_Can)                                                               :amp/meter**2
    '''

    #iStim, General, Time-Dependent Stimulus
    i_Stim = '''
    n_Stim                                                                                      :meter**2
    mag                                                                                         :amp
    iStm = (t > 1*second)*(t < 2*second)*mag/n_Stim                                             :amp/meter**2
    '''

    #iSyn, Synapse Currents, RE
    i_Syn = '''
    igj                                                                                         :amp
    igjHR                                                                                       :amp
    iHTCa                                                                                       :amp
    iHTCb                                                                                       :amp
    iRTCLa                                                                                      :amp
    iRTCRa                                                                                      :amp
    iRTCLb                                                                                      :amp
    iRTCRb                                                                                      :amp
    iIN                                                                                         :amp
    iRE                                                                                         :amp
    iPYLa                                                                                       :amp
    iPYRa                                                                                       :amp
    iPYLb                                                                                       :amp
    iPYRb                                                                                       :amp
    iRND                                                                                        :amp
    iGLUT = iHTCa + iHTCb + iRTCLa + iRTCLb + iRTCRa + iRTCRb + iPYLa + iPYRb + iRND            :amp
    iGABA = iIN + iRE                                                                           :amp
    iSyn = (igj + igjHR + iGLUT + iGABA)/n_Stim                                                 :amp/meter**2
    '''

    #Location equations
    loc = '''
    x                                                                                           :1
    y                                                                                           :1
    gjb                                                                                         :1
    '''

    #Main Izhikevich Equation, PY FS
    main_PY_FS = '''
    dv/dt = (k*(v - vr)*(v - vt) - u - iSyn + iStm + sg*Sgn(t))/Cm                              :volt
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
    iRTCLa                                                                                      :amp
    iRTCRa                                                                                      :amp
    iPYL                                                                                        :amp
    iPYR                                                                                        :amp
    iFS                                                                                         :amp
    iFSL                                                                                        :amp
    iFSR                                                                                        :amp
    iRND                                                                                        :amp
    iStm                                                                                        :amp
    iSyn = iRTCLa + iRTCRa + iPYL + iPYR + iFS + iRND                                           :amp
    LFP = abs(iRTCLa + iRTCRa + iPYL + iPYR) + abs(iFS)                                         :amp
    '''

    #Cell Group Declarations
    #TC Equations
    TC_eqs = (main_TC + Ca_TC + i_Ca + i_Cat + i_Cal + i_Na + i_K + i_H + i_AHP + i_Can + i_Syn + i_Stim + loc)

    #IN Equations
    IN_eqs = (main_IN + Ca_RE_IN + i_Ca + i_Na + i_K + i_H + i_AHP + i_Can + i_Syn + i_Stim)

    #RE Equations
    RE_eqs = (main_RE + Ca_RE_IN + i_Ca_RE + i_Na + i_K + i_AHP + i_Can + i_Syn + i_Stim + loc)

    #Declare TCL Group
    TCLg = NeuronGroup(nHTC**2+nRTC**2,TC_eqs,method = 'rk4',dt = t_step,threshold =
            'v > 0*mV',refractory = 3*ms, name = "TCLg")
    TCLg.e_L = -70*mV
    TCLg.g_KL = 0*mS/cm**2
    TCLg.g_L = gL + (rand(nHTC**2+nRTC**2) - 0.5)*0.005*mS/cm**2
    TCLg.tm  = 0.25
    TCLg.g_H = 0.01*mS/cm**2
    TCLg.n_Stim = 2.9e-4*cm**2
    left.append(TCLg)

    #Declare TCR Group
    TCRg = NeuronGroup(nHTC**2+nRTC**2,TC_eqs,method = 'rk4',dt = t_step,threshold =
            'v > 0*mV',refractory = 3*ms, name = "TCRg")
    TCRg.e_L = -70*mV
    TCRg.g_KL = 0*mS/cm**2
    TCRg.g_L = gL + (rand(nHTC**2+nRTC**2) - 0.5)*0.005*mS/cm**2
    TCRg.tm  = 0.25
    TCRg.g_H = 0.01*mS/cm**2
    TCRg.n_Stim = 2.9e-4*cm**2
    right.append(TCRg)

    #Declare HTCL Group
    HTCLg = Subgroup(source = TCLg, start = 0, stop = nHTC**2, name = "HTCLg")
    HTCLg.g_Ca = 3*mS/cm**2
    HTCLg.g_Cal = 0.5*mS/cm**2
    HTCLg.v_SH = -30*mV
    HTCLg.g_AHP = 0.3*mS/cm**2
    HTCLg.g_Can = 0.5*mS/cm**2
    HTCLg.x = 'i%nHTC'
    HTCLg.y = '(i - x)/nHTC'
    left.append(HTCLg)

    #Declare HTCR Group
    HTCRg = Subgroup(source = TCRg, start = 0, stop = nHTC**2, name = "HTCRg")
    HTCRg.g_Ca = 3*mS/cm**2
    HTCRg.g_Cal = 0.5*mS/cm**2
    HTCRg.v_SH = -30*mV
    HTCRg.g_AHP = 0.3*mS/cm**2
    HTCRg.g_Can = 0.5*mS/cm**2
    HTCRg.x = 'i%nHTC'
    HTCRg.y = '(i - x)/nHTC'
    right.append(HTCRg)

    #Declare RTCL Group
    RTCLg = Subgroup(source = TCLg, start = nHTC**2,  stop = len(TCLg), name = "RTCLg")
    RTCLg.g_Ca = 0.6*mS/cm**2
    RTCLg.g_Cal = 0.3*mS/cm**2
    RTCLg.v_SH = -40*mV
    RTCLg.g_AHP = 0.1*mS/cm**2
    RTCLg.g_Can = 0.6*mS/cm**2
    RTCLg.x = 'i%nRTC'
    RTCLg.y = '(i - x)/nRTC'
    RTCLg.gjb = rand(nRTC*nRTC) < 0.2
    left.append(RTCLg)

    #Declare RTCR Group
    RTCRg = Subgroup(source = TCRg, start = nHTC**2,  stop = len(TCRg), name = "RTCRg")
    RTCRg.g_Ca = 0.6*mS/cm**2
    RTCRg.g_Cal = 0.3*mS/cm**2
    RTCRg.v_SH = -40*mV
    RTCRg.g_AHP = 0.1*mS/cm**2
    RTCRg.g_Can = 0.6*mS/cm**2
    RTCRg.x = 'i%nRTC'
    RTCRg.y = '(i - x)/nRTC'
    RTCRg.gjb = rand(nRTC*nRTC) < 0.2
    right.append(RTCRg)

    #Declare INL Group
    INLg = NeuronGroup(nIN*nIN,IN_eqs,method = 'rk4',dt = t_step,threshold = 'v > 0*mV',
            refractory = 3*ms, name = "INLg")
    INLg.e_L = -60*mV
    INLg.g_KL = 0.02*mS/cm**2
    INLg.g_L = gL + (rand(nIN*nIN) - 0.5)*0.005*mS/cm**2
    INLg.tau_Ca = 10*ms
    INLg.g_Ca = 2.5*mS/cm**2
    INLg.v_SH = -30*mV
    INLg.tm  = 0.25
    INLg.g_H = 0.05*mS/cm**2
    INLg.g_AHP = 0.2*mS/cm**2
    INLg.g_Can = 0.1*mS/cm**2
    INLg.n_Stim = 1.7e-4*cm**2
    left.append(INLg)

    #Declare INR Group
    INRg = NeuronGroup(nIN*nIN,IN_eqs,method = 'rk4',dt = t_step,threshold = 'v > 0*mV',
            refractory = 3*ms, name = "INRg")
    INRg.e_L = -60*mV
    INRg.g_KL = 0.02*mS/cm**2
    INRg.g_L = gL + (rand(nIN*nIN) - 0.5)*0.005*mS/cm**2
    INRg.tau_Ca = 10*ms
    INRg.g_Ca = 2.5*mS/cm**2
    INRg.v_SH = -30*mV
    INRg.tm  = 0.25
    INRg.g_H = 0.05*mS/cm**2
    INRg.g_AHP = 0.2*mS/cm**2
    INRg.g_Can = 0.1*mS/cm**2
    INRg.n_Stim = 1.7e-4*cm**2
    right.append(INRg)

    #Declare REL Group
    RELg = NeuronGroup(nRE*nRE,RE_eqs,method = 'rk4',dt = t_step,threshold = 'v > 0*mV',
            refractory = 3*ms, name = "RELg")
    RELg.e_L = -60*mV
    RELg.g_KL = 0.01*mS/cm**2
    RELg.g_L = gL + (rand(nRE*nRE) - 0.5)*0.005*mS/cm**2
    RELg.tau_Ca = 100*ms
    RELg.v_SH = -40*mV
    RELg.tm  = 1
    RELg.g_AHP = 0.2*mS/cm**2
    RELg.g_Can = 0.2*mS/cm**2
    RELg.n_Stim = 1.45e-4*cm**2
    RELg.x = 'i%nRE'
    RELg.y = '(i - x)/nRE'
    RELg.gjb = rand(nRE*nRE) < 0.2
    left.append(RELg)

    #Declare RER Group
    RERg = NeuronGroup(nRE*nRE,RE_eqs,method = 'rk4',dt = t_step,threshold = 'v > 0*mV',
            refractory = 3*ms, name = "RERg")
    RERg.e_L = -60*mV
    RERg.g_KL = 0.01*mS/cm**2
    RERg.g_L = gL + (rand(nRE*nRE) - 0.5)*0.005*mS/cm**2
    RERg.tau_Ca = 100*ms
    RERg.v_SH = -40*mV
    RERg.tm  = 1
    RERg.g_AHP = 0.2*mS/cm**2
    RERg.g_Can = 0.2*mS/cm**2
    RERg.n_Stim = 1.45e-4*cm**2
    RERg.x = 'i%nRE'
    RERg.y = '(i - x)/nRE'
    RERg.gjb = rand(nRE*nRE) < 0.2
    right.append(RERg)

    #Declare PY FS Groups
    FSLg = NeuronGroup(nFS, main_PY_FS, method = 'rk4', dt = t_step, threshold = 'v>= vth', reset = '''
    v = c
    u += d
    ''', name = "FSLg")
    left.append(FSLg)

    PYLg = NeuronGroup(nPY, main_PY_FS, method = 'rk4', dt = t_step, threshold = 'v>= vth', reset = '''
    v = c
    u += d
    ''', name = "PYLg")
    left.append(PYLg)


    FSRg = NeuronGroup(nFS, main_PY_FS, method = 'rk4', dt = t_step, threshold = 'v>= vth', reset = '''
    v = c
    u += d
    ''', name = "FSRg")
    right.append(FSRg)

    PYRg = NeuronGroup(nPY, main_PY_FS, method = 'rk4', dt = t_step, threshold = 'v>= vth', reset = '''
    v = c
    u += d
    ''', name = "PYRg")
    right.append(PYRg)



    #PYL Group Declaration
    PYLg.Cm = (100 + 0.1*randn(nPY))*pF
    PYLg.k = 0.7*pA/(mvolt*mvolt)
    PYLg.vr = -60*mV + .1*randn(nPY)*mV
    PYLg.vt = -40*mV + .1*randn(nPY)*mV
    PYLg.a = (0.03 + 0.001*randn(nPY))*kHz
    PYLg.b = (-2 + 0.01*randn(nPY))*pA
    PYLg.vb = -200*volt
    PYLg.un = 1
    PYLg.vth = 35*mV + .1*randn(nPY)*mV
    PYLg.c = -50*mV + 0.1*randn(nPY)*mV
    PYLg.d = (100 + 0.1*randn(nPY))*pA
    PYLg.v -60*mV + .1*randn(nPY)*mV
    PYLg.u = PYLg.b*(PYLg.v-PYLg.vr)/mV
    PYLg.iStm = 79*pA
    PYLg.sg = 1
    left.append(PYLg)

    #PYR Group Declaration
    PYRg.Cm = (100 + 0.1*randn(nPY))*pF
    PYRg.k = 0.7*pA/(mvolt*mvolt)
    PYRg.vr = -60*mV + .1*randn(nPY)*mV
    PYRg.vt = -40*mV + .1*randn(nPY)*mV
    PYRg.a = (0.03 + 0.001*randn(nPY))*kHz
    PYRg.b = (-2 + 0.01*randn(nPY))*pA
    PYRg.vb = -200*volt
    PYRg.un = 1
    PYRg.vth = 35*mV + .1*randn(nPY)*mV
    PYRg.c = -50*mV + 0.1*randn(nPY)*mV
    PYRg.d = (100 + 0.1*randn(nPY))*pA
    PYRg.v = -60*mV + .1*randn(nPY)*mV
    PYRg.u = PYRg.b*(PYRg.v-PYRg.vr)/mV
    PYRg.iStm = 79*pA
    PYRg.sg = 1
    right.append(PYRg)


    #FS Group Declaration
    FSLg.Cm = (20 + 0.1*randn(nFS))*pF
    FSLg.k = (1 + 0.01*randn(nFS))*pA/(mvolt*mvolt)
    FSLg.vr = -55*mV + 0.1*randn(nFS)*mV
    FSLg.vt = -40*mV + 0.1*randn(nFS)*mV
    FSLg.a = (0.2 + 0.001*randn(nFS))*kHz
    FSLg.b = (0.025 + 0.0001*randn(nFS))*pA
    FSLg.vb = FSLg.vr
    FSLg.un = 3
    FSLg.vth = 25*mV + 0.1*randn(nFS)*mV
    FSLg.c = -45*mV + 0.1*randn(nFS)*mV
    FSLg.d = (0 + 0.001*randn(nFS))*pA
    FSLg.v = -55*mV + 0.1*randn(nFS)*mV
    FSLg.u = (0 + 0.1*randn(nFS))*pA
    FSLg.iStm = 60*pA
    left.append(FSLg)

    #FS Group Declaration
    FSRg.Cm = (20 + 0.1*randn(nFS))*pF
    FSRg.k = (1 + 0.01*randn(nFS))*pA/(mvolt*mvolt)
    FSRg.vr = -55*mV + 0.1*randn(nFS)*mV
    FSRg.vt = -40*mV + 0.1*randn(nFS)*mV
    FSRg.a = (0.2 + 0.001*randn(nFS))*kHz
    FSRg.b = (0.025 + 0.0001*randn(nFS))*pA
    FSRg.vb = FSRg.vr
    FSRg.un = 3
    FSRg.vth = 25*mV + 0.1*randn(nFS)*mV
    FSRg.c = -45*mV + 0.1*randn(nFS)*mV
    FSRg.d = (0 + 0.001*randn(nFS))*pA
    FSRg.v = -55*mV + 0.1*randn(nFS)*mV
    FSRg.u = (0 + 0.1*randn(nFS))*pA
    FSRg.iStm = 60*pA
    right.append(FSRg)

    #Constant Spiking Neuron for Noise
    RNDg = NeuronGroup(1, model = ' v = 0 :1', method = 'rk4', dt = t_step,
            threshold = 'v < 1', name = "RNDg")
    other.append(RNDg)


    #Declare Initial Conditions Thalamus Cells
    #Global Variables
    ca0 = 5e-5*mM
    v0 = -60*mV

    #iCa, HTC, RTC, IN
    m00 = 1/(1+exp(-(v0+34*mV)/(6.2*mV)))
    h00 = 1/(1+exp((v0+58*mV)/(4*mV)))

    #iCa, RE
    m00_RE = 1/(1+exp(-(v0+55*mV)/(7.4*mV)))
    h00_RE = 1/(1+exp((v0+83*mV)/(5*mV)))

    #iCat HTC, RTC
    m10 = 1/(1+exp(-(v0+62*mV)/(6.2*mV)))
    h10 = 1/(1+exp((v0+86*mV)/(4*mV)))

    #iCal HTC, RTC
    m20 = 1/(1+exp(-(v0+10*mV)/(4*mV)))
    h20 = 1/(1+exp((v0+25*mV)/(2*mV)))

    #iNa, HTC, IN
    m30a_a = (0.32*ms**-1*mV**-1)*(v0+30*mV-13*mV)/(1-exp(-(v0+30*mV-13*mV)/(4*mV)))
    m30b_a = (-0.28*mV**-1*ms**-1)*(v0+30*mV-40*mV)/(1-exp((v0+30*mV-40*mV)/(5*mV)))
    m30_a = m30a_a/(m30a_a + m30b_a)
    h30a_a = (0.128*ms**-1)*exp(-(v0+30*mV-17*mV)/(18*mV))
    h30b_a = (4*ms**-1)/(1+exp(-(v0+30*mV-40*mV)/(5*mV)))
    h30_a = h30a_a/(h30a_a + h30b_a)

    #iNa, RTC, RE
    m30a_b = (0.32*ms**-1*mV**-1)*(v0+40*mV-13*mV)/(1-exp(-(v0+40*mV-13*mV)/(4*mV)))
    m30b_b = (-0.28*mV**-1*ms**-1)*(v0+40*mV-40*mV)/(1-exp((v0+40*mV-40*mV)/(5*mV)))
    m30_b = m30a_b/(m30a_b + m30b_b)
    h30a_b = (0.128*ms**-1)*exp(-(v0+40*mV-17*mV)/(18*mV))
    h30b_b = (4*ms**-1)/(1+exp(-(v0+40*mV-40*mV)/(5*mV)))
    h30_b = h30a_b/(h30a_b + h30b_b)

    #iK, HTC, IN
    m40a_a = (0.032*ms**-1*mV**-1)*(v0+30*mV-15*mV)/(1-exp(-(v0+30*mV-15*mV)/(5*mV)))
    m40b_a = (0.5*ms**-1)*exp(-(v0+30*mV-10*mV)/(40*mV))
    m40_a = m40a_a/(m40a_a + m40b_a)

    #iK, RTC, RE
    m40a_b = (0.032*ms**-1*mV**-1)*(v0+40*mV-15*mV)/(1-exp(-(v0+40*mV-15*mV)/(5*mV)))
    m40b_b = (0.5*ms**-1)*exp(-(v0+40*mV-10*mV)/(40*mV))
    m40_b = m40a_b/(m40a_b + m40b_b)

    #iH
    m50 = 1/(1+exp((v0+75*mV)/(5.5*mV)))

    #iAHP
    m60 = 48*ca0**2/(48*ca0**2+0.09*mM**2)

    #iCan
    m70 = 1/(1+exp(-(v0+43*mV)/(5.2*mV)))

    #HTC Initial Conditions
    int_HTCg = {'v':[v0],'ca':[ca0],'m0':[m00],'h0':[h00],'m1':[m10],'h1':[h10],'m2':[m20],'h2':[h20],'m3':[m30_a],'h3':[h30_a],'m4':[m40_a],'m5':[m50],'m6':[m60],'m7':[m70]}
    HTCLg.set_states(int_HTCg)
    HTCRg.set_states(int_HTCg)

    #RTC Initial Conditions
    int_RTCg = {'v':[v0],'ca':[ca0],'m0':[m00],'h0':[h00],'m1':[m10],'h1':[h10],'m2':[m20],'h2':[h20],'m3':[m30_b],'h3':[h30_b],'m4':[m40_b],'m5':[m50],'m6':[m60],'m7':[m70]}
    RTCLg.set_states(int_RTCg)
    RTCRg.set_states(int_RTCg)

    #IN Initial Conditions
    int_INg = {'v':[v0],'ca':[ca0],'m0':[m00],'h0':[h00],'m3':[m30_a],'h3':[h30_a],'m4':[m40_a],'m5':[m50],'m6':[m60],'m7':[m70]}
    INLg.set_states(int_INg)
    INRg.set_states(int_INg)


    #RE Initial Conditions
    int_REg = {'v':[v0],'ca':[ca0],'m0':[m00_RE],'h0':[h00_RE],'m3':[m30_b],'h3':[h30_b],'m4':[m40_b],'m6':[m60],'m7':[m70]}
    RELg.set_states(int_REg)
    RERg.set_states(int_REg)

    #Modular Synapse Equations, Biophysical
    #Gap Junction
    GJ = '''
    rgap                                                                                        :ohm
    iGJ = (v_post - v_pre)/rgap                                                                 :amp
    '''

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

    #Random Spiking Synapse, HTC, RTC, IN, RE
    iRND = '''
    tau                                                                                         :second
    t_spike                                                                                     :second
    g_spike                                                                                     :siemens
    d                                                                                           :1
    g                                                                                           :siemens
    iRND_post = g*v_post                                                                        :amp (summed)
    '''
    preRND = '''
    d = (rand() < 0.003)*(t-t_spike > 3*ms)
    t_spike = t*d + t_spike*(1-d)
    g = (g + g_spike*d)*exp(-(t-t_spike)/tau)
    '''

    #->Thalamus Equation Expressions
    igj_eqs = '''
    igj_post = iGJ                                                                              :amp (summed)
    '''
    igjHR_eqs = '''
    igjHR_post = iGJ                                                                            :amp (summed)
    '''
    iHTCa_eqs = '''
    iHTCa_post = iCS                                                                            :amp (summed)
    '''
    iHTCb_eqs = '''
    iHTCb_post = iCS                                                                            :amp (summed)
    '''

    iRTCLa_eqs = '''
    iRTCLa_post = iCS                                                                            :amp (summed)
    '''
    iRTCLb_eqs = '''
    iRTCLb_post = iCS                                                                            :amp (summed)
    '''

    iRTCRa_eqs = '''
    iRTCRa_post = iCS                                                                            :amp (summed)
    '''
    iRTCRb_eqs = '''
    iRTCRb_post = iCS                                                                            :amp (summed)
    '''

    iIN_eqs = '''
    iIN_post = iCS                                                                              :amp (summed)
    '''
    iRE_eqs = '''
    iRE_post = iCS                                                                              :amp (summed)
    '''

    iPYLa_eqs = '''
    iPYLa_post = iCS                                                                             :amp (summed)
    '''

    iPYRa_eqs = '''
    iPYRa_post = iCS                                                                             :amp (summed)
    '''


    iPYLb_eqs = '''
    iPYLb_post = iCS                                                                             :amp (summed)
    '''

    iPYRb_eqs = '''
    iPYRb_post = iCS                                                                             :amp (summed)
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

    #Random Spiking Synapse, PY, FS
    RND_PYFS = '''
    std                                                                                         :amp
    t_spike                                                                                     :second
    randi                                                                                       :1
    iRND_post = std*randi                                                                       :amp(summed)
    '''
    preRND_PYFS = '''
    randi = randi*(t-t_spike<0.5*ms)+randn()*(1-(t-t_spike<0.5*ms))
    t_spike = t_spike*(t-t_spike<0.5*ms)+t*(1-(t-t_spike<0.5*ms))
    '''

    #->PYFS Equation Expressions
    iPYL_eqs = '''
    iPYL_post = iCS                                                                             :amp (summed)
    '''

    iPYR_eqs = '''
    iPYR_post = iCS                                                                             :amp (summed)
    '''

    iFS_eqs = '''
    iFS_post = iCS                                                                              :amp (summed)
    '''
    iFSL_eqs = '''
    iFSL_post = iCS                                                                             :amp (summed)
    '''
    iFSR_eqs = '''
    iFSR_post = iCS                                                                             :amp (summed)
    '''
    #Synapse Group Declarations

    #Intra-hemisphere


    #->HTCL Synapses
    HTCL_HTCL_gj = Synapses(HTCLg,HTCLg,(GJ + igj_eqs), name = "HTCL_HTCL_gj")
    HTCL_RTCL_gj = Synapses(RTCLg,HTCLg,(GJ + igjHR_eqs), name = "HTCL_RTCL_gj")
    HTCL_REL_cs = Synapses(RELg,HTCLg,(CSa + iRE_eqs),on_pre = preCS,
            method = 'rk4',dt = t_step, name = "HTCL_REL_cs")
    HTCL_RND_cs = Synapses(RNDg, HTCLg, iRND, on_pre = preRND,
            method = 'rk4', dt = t_step, name = "HTCL_RND_cs")
    #->HTCR Synapses
    HTCR_HTCR_gj = Synapses(HTCRg,HTCRg,(GJ + igj_eqs), name = "HTCR_HTCR_gj")
    HTCR_RTCR_gj = Synapses(RTCRg,HTCRg,(GJ + igjHR_eqs), name = "HTCR_RTCT_gj")
    HTCR_RER_cs = Synapses(RERg,HTCRg,(CSa + iRE_eqs),on_pre = preCS,
            method = 'rk4',dt = t_step, name = "HTCR_RER_cs")
    HTCR_RND_cs = Synapses(RNDg, HTCRg, iRND, on_pre = preRND,
            method = 'rk4', dt = t_step, name = "HTCR_RND_cs")

    #->RTCL Synapses
    RTCL_HTCL_gj = Synapses(HTCLg,RTCLg,(GJ + igjHR_eqs), name = "RTCL_HTCL_gj")
    RTCL_INL_cs = Synapses(INLg,RTCLg,(CSa + iIN_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "RTCL_INL_cs")
    RTCL_REL_cs = Synapses(RELg,RTCLg,(CSa + iRE_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "RTCL_REL_cs")
    RTCL_PYL_csa = Synapses(PYLg, RTCLg, (CSa + iPYLa_eqs), on_pre = preCS,
            method = 'rk4', dt = t_step, name = "RTCL_PYL_csa")
    RTCL_PYL_csb = Synapses(PYLg, RTCLg, (CSb + iPYLb_eqs), on_pre = preCS,
            method = 'rk4', dt = t_step, name = "RTCL_PYL_csb")
    RTCL_RND_cs = Synapses(RNDg, RTCLg, iRND, on_pre = preRND,
            method = 'rk4', dt = t_step, name = "RTCL_RND_cs")
    #->RTCR Synapses
    RTCR_HTCR_gj = Synapses(HTCRg,RTCRg,(GJ + igjHR_eqs))
    RTCR_INR_cs = Synapses(INRg,RTCRg,(CSa + iIN_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "RCTR_INR_cs")
    RTCR_RER_cs = Synapses(RERg,RTCRg,(CSa + iRE_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "RTCR_RER_cs")
    RTCR_PYR_csa = Synapses(PYRg, RTCRg, (CSa + iPYRa_eqs), on_pre = preCS,
            method = 'rk4', dt = t_step, name = "RTCR_PYR_csa")
    RTCR_PYR_csb = Synapses(PYRg, RTCRg, (CSb + iPYRb_eqs), on_pre = preCS,
            method = 'rk4', dt = t_step, name = "RTCR_PYR_csb")
    RTCR_RND_cs = Synapses(RNDg, RTCRg, iRND, on_pre = preRND,
            method = 'rk4', dt = t_step, name = "RTCR_RND_cs")

    #->INL Synapses
    INL_HTCL_csa = Synapses(HTCLg,INLg,(CSa + iHTCa_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "INL_HTCL_csa")
    INL_HTCL_csb = Synapses(HTCLg,INLg,(CSb + iHTCb_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "INL_HTCL_csb")
    INL_REL_cs = Synapses(RELg,INLg,(CSa + iRE_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "INL_REL_cs")
    INL_RND_cs = Synapses(RNDg, INLg, iRND, on_pre = preRND,
            method = 'rk4', dt = t_step, name = "INL_RND_cs")
    #->INR Synapses
    INR_HTCR_csa = Synapses(HTCRg,INRg,(CSa + iHTCa_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "INR_HTCR_csa")
    INR_HTCR_csb = Synapses(HTCRg,INRg,(CSb + iHTCb_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "INR_HTCR_csb")
    INR_RER_cs = Synapses(RERg,INRg,(CSa + iRE_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "INR_RER_cs")
    INR_RND_cs = Synapses(RNDg, INRg, iRND, on_pre = preRND,
            method = 'rk4', dt = t_step, name = "INR_RND_cs")

    #->REL Synapses
    REL_REL_gj = Synapses(RELg,RELg,(GJ + igj_eqs))
    REL_HTCL_csa = Synapses(HTCLg,RELg,(CSa + iHTCa_eqs),on_pre = preCS,
            method = 'rk4', dt = t_step, name = "REL_HTCL_csa")
    REL_HTCL_csb = Synapses(HTCLg,RELg,(CSb + iHTCb_eqs),on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "REL_HTCL_csb")
    REL_RTCL_csa = Synapses(RTCLg,RELg,(CSa + iRTCLa_eqs),on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "REL_RTCL_csa")
    REL_RTCL_csb = Synapses(RTCLg,RELg,(CSb + iRTCLb_eqs),on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "REL_RTCL_csb")
    REL_REL_cs = Synapses(RELg,RELg,(CSa + iRE_eqs),on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "REL_REL_cs")
    REL_PYL_csa = Synapses(PYLg, RELg, (CSa + iPYLa_eqs), on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "REL_PYL_csa")
    REL_PYL_csb = Synapses(PYLg, RELg, (CSb + iPYLb_eqs), on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "REL_PYL_csb")
    REL_RND_cs = Synapses(RNDg, RELg, iRND, on_pre = preRND, 
            method = 'rk4', dt = t_step, name = "REL_RND_cs")
    #->RER Synapses
    RER_RER_gj = Synapses(RERg,RERg,(GJ + igj_eqs))
    RER_HTCR_csa = Synapses(HTCRg,RERg,(CSa + iHTCa_eqs),on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "RER_HTCR_csa")
    RER_HTCR_csb = Synapses(HTCRg,RERg,(CSb + iHTCb_eqs),on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "RER_HTCR_csb")
    RER_RTCR_csa = Synapses(RTCRg,RERg,(CSa + iRTCRa_eqs),on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "RER_RTCR_csa")
    RER_RTCR_csb = Synapses(RTCRg,RERg,(CSb + iRTCRb_eqs),on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "RER_RTCR_csb")
    RER_RER_cs = Synapses(RERg,RERg,(CSa + iRE_eqs),on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "RER_RER_cs")
    RER_PYR_csa = Synapses(PYRg, RERg, (CSa + iPYRa_eqs), on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "RER_PYR_csa")
    RER_PYR_csb = Synapses(PYRg, RERg, (CSb + iPYRb_eqs), on_pre = preCS, 
            method = 'rk4', dt = t_step, name = "RER_PYR_csb")
    RER_RND_cs = Synapses(RNDg, RERg, iRND, on_pre = preRND, 
            method = 'rk4', dt = t_step, name = "RER_RND_cs")

    #->PYL Synapses
    PYL_RTCL_cs = Synapses(RTCLg,PYLg,(PYFS_cs + iRTCLa_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "PYL_RTCL_cs")
    PYL_PYL_cs = Synapses(PYLg,PYLg,(PYFS_cs + iPYL_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "PYL_PYL_cs")
    PYL_FSL_cs = Synapses(FSLg,PYLg,(PYFS_cs + iFS_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "PYL_FSL_cs")
    PYL_RND = Synapses(RNDg,PYLg,RND_PYFS,on_pre = preRND_PYFS, 
            method = 'rk4', dt = t_step, name = "PYL_RND")
    #->PY Synapses
    PYR_RTCR_cs = Synapses(RTCRg,PYRg,(PYFS_cs + iRTCRa_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "PYR_RTCR_cs")
    PYR_PYR_cs = Synapses(PYRg,PYRg,(PYFS_cs + iPYR_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "PYR_PYR_cs")
    PYR_FSR_cs = Synapses(FSRg,PYRg,(PYFS_cs + iFS_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "PYR_FSR_cs")
    PYR_RND = Synapses(RNDg,PYRg,RND_PYFS,on_pre = preRND_PYFS, 
            method = 'rk4', dt = t_step, name = "PYR_RND")

    #->FSL Synapses
    FSL_RTCL_cs = Synapses(RTCLg,FSLg,(PYFS_cs + iRTCLa_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "FSL_RTCL_cs")
    FSL_PYL_cs = Synapses(PYLg,FSLg,(PYFS_cs + iPYL_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "FSL_PYL_cs")
    FSL_FSL_cs = Synapses(FSLg,FSLg,(PYFS_cs + iFS_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "FSL_FSL_cs")
    FSL_RND = Synapses(RNDg,FSLg,RND_PYFS,on_pre = preRND_PYFS, 
            method = 'rk4', dt = t_step, name = "FSL_RND")

    #->FSR Synapses
    FSR_RTCR_cs = Synapses(RTCRg,FSRg,(PYFS_cs + iRTCRa_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "FSR_RTCT_cs")
    FSR_PYR_cs = Synapses(PYRg,FSRg,(PYFS_cs + iPYR_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "FSR_PYR_cs")
    FSR_FSR_cs = Synapses(FSRg,FSRg,(PYFS_cs + iFS_eqs),on_pre = prePYFS, 
            method = 'rk4', dt = t_step, name = "FSR_FSR_cs")
    FSR_RND = Synapses(RNDg,FSRg,RND_PYFS,on_pre = preRND_PYFS, 
            method = 'rk4', dt = t_step, name = "FSR_RND")

    #Declare Connections

    #->HTCL
    HTCL_HTCL_gj.connect('(x_pre - x_post)**2+(y_pre-y_post)**2 < 4.01*(i < j)' , p = '0.3')
    HTCL_HTCL_gj.connect(i = HTCL_HTCL_gj.j[:], j = HTCL_HTCL_gj.i[:])
    HTCL_RTCL_gj.connect('(x_pre*(nHTC-1.0)/(nRTC-1.0)-x_post)**2+(y_pre*(nHTC-1.0)/(nRTC-1.0)-y_post)**2 < 4.01*gjb_pre', p = '0.3')
    HTCL_REL_cs.connect(p = '0.2')
    HTCL_RND_cs.connect()
    #->HTCR
    HTCR_HTCR_gj.connect('(x_pre - x_post)**2+(y_pre-y_post)**2 < 4.01*(i < j)' , p = '0.3')
    HTCR_HTCR_gj.connect(i = HTCR_HTCR_gj.j[:], j = HTCR_HTCR_gj.i[:])
    HTCR_RTCR_gj.connect('(x_pre*(nHTC-1.0)/(nRTC-1.0)-x_post)**2+(y_pre*(nHTC-1.0)/(nRTC-1.0)-y_post)**2 < 4.01*gjb_pre', p = '0.3')
    HTCR_RER_cs.connect(p = '0.2')
    HTCR_RND_cs.connect()

    #->RTCL
    RTCL_HTCL_gj.connect(i = HTCL_RTCL_gj.j[:], j = HTCL_RTCL_gj.i[:])
    RTCL_INL_cs.connect(p = '0.3')
    RTCL_REL_cs.connect(p = '0.2')
    RTCL_PYL_csa.connect(p = '0.23')
    RTCL_PYL_csb.connect(i = RTCL_PYL_csa.i[:], j = RTCL_PYL_csa.j[:])
    RTCL_RND_cs.connect()
    #->RTCR
    RTCR_HTCR_gj.connect(i = HTCR_RTCR_gj.j[:], j = HTCR_RTCR_gj.i[:])
    RTCR_INR_cs.connect(p = '0.3')
    RTCR_RER_cs.connect(p = '0.2')
    RTCR_PYR_csa.connect(p = '0.23')
    RTCR_PYR_csb.connect(i = RTCR_PYR_csa.i[:], j = RTCR_PYR_csa.j[:])
    RTCR_RND_cs.connect()

    #->INL
    INL_HTCL_csa.connect(p = '0.3')
    INL_HTCL_csb.connect(i = INL_HTCL_csa.i[:], j = INL_HTCL_csa.j[:])
    INL_REL_cs.connect(p = '0.05')
    INL_RND_cs.connect()
    #->INR
    INR_HTCR_csa.connect(p = '0.3')
    INR_HTCR_csb.connect(i = INR_HTCR_csa.i[:], j = INR_HTCR_csa.j[:])
    INR_RER_cs.connect(p = '0.05')
    INR_RND_cs.connect()

    #->REL
    REL_REL_gj.connect('(x_pre - x_post)**2+(y_pre-y_post )**2 < 4.01*(gjb_pre*i < gjb_post*j)', p = '0.3')
    REL_REL_gj.connect(i = REL_REL_gj.j[:], j = REL_REL_gj.i[:])
    REL_HTCL_csa.connect(p = '0.2')
    REL_HTCL_csb.connect(i = REL_HTCL_csa.i[:], j = REL_HTCL_csa.j[:])
    REL_RTCL_csa.connect(p = '0.2')
    REL_RTCL_csb.connect(i = REL_RTCL_csa.i[:], j = REL_RTCL_csa.j[:])
    REL_REL_cs.connect(p = '0.2')
    REL_PYL_csa.connect(p = '0.3')
    REL_PYL_csb.connect(i = REL_PYL_csa.i[:], j = REL_PYL_csa.j[:])
    REL_RND_cs.connect()
    #->RER
    RER_RER_gj.connect('(x_pre - x_post)**2+(y_pre-y_post )**2 < 4.01*(gjb_pre*i < gjb_post*j)', p = '0.3')
    RER_RER_gj.connect(i = RER_RER_gj.j[:], j = RER_RER_gj.i[:])
    RER_HTCR_csa.connect(p = '0.2')
    RER_HTCR_csb.connect(i = RER_HTCR_csa.i[:], j = RER_HTCR_csa.j[:])
    RER_RTCR_csa.connect(p = '0.2')
    RER_RTCR_csb.connect(i = RER_RTCR_csa.i[:], j = RER_RTCR_csa.j[:])
    RER_RER_cs.connect(p = '0.2')
    RER_PYR_csa.connect(p = '0.3')
    RER_PYR_csb.connect(i = RER_PYR_csa.i[:], j = RER_PYR_csa.j[:])
    RER_RND_cs.connect()

    #->PYL
    PYL_RTCL_cs.connect(p = '0.04')
    PYL_PYL_cs.connect(p = '0.5')
    PYL_FSL_cs.connect('abs(i*(nPY-1.0)/(nFS-1.0) - j) < nPY*0.2', p = '0.8')
    PYL_RND.connect()
    #->PYR
    PYR_RTCR_cs.connect(p = '0.04')
    PYR_PYR_cs.connect(p = '0.5')
    PYR_FSR_cs.connect('abs(i*(nPY-1.0)/(nFS-1.0) - j) < nPY*0.2', p = '0.8')
    PYR_RND.connect()

    #->FS
    FSL_RTCL_cs.connect(p = '0.02')
    FSL_PYL_cs.connect(i = PYL_FSL_cs.j[:], j = PYL_FSL_cs.i[:])
    FSL_FSL_cs.connect('abs(i-j) < nFS/2.0*(i!=j)', p = '0.8')
    FSL_RND.connect()
    #->FS
    FSR_RTCR_cs.connect(p = '0.02')
    FSR_PYR_cs.connect(i = PYR_FSR_cs.j[:], j = PYR_FSR_cs.i[:])
    FSR_FSR_cs.connect('abs(i-j) < nFS/2.0*(i!=j)', p = '0.8')
    FSR_RND.connect()

    #Declare Synapse Variables

    #->HTCL
    HTCL_HTCL_gj.rgap = 100*Mohm
    HTCL_RTCL_gj.rgap = 300*Mohm
    HTCL_REL_cs.D_i = 1.07
    HTCL_REL_cs.g_syn = 3*nS
    HTCL_REL_cs.e_syn = -80*mV
    HTCL_REL_cs.alpha = 10.5/ms
    HTCL_REL_cs.beta = 0.166/ms
    HTCL_RND_cs.t_spike = -3*ms
    HTCL_RND_cs.g_spike = 1.5*nS
    HTCL_RND_cs.tau = 5*ms
    #->HTCR
    HTCR_HTCR_gj.rgap = 100*Mohm
    HTCR_RTCR_gj.rgap = 300*Mohm
    HTCR_RER_cs.D_i = 1.07
    HTCR_RER_cs.g_syn = 3*nS
    HTCR_RER_cs.e_syn = -80*mV
    HTCR_RER_cs.alpha = 10.5/ms
    HTCR_RER_cs.beta = 0.166/ms
    HTCR_RND_cs.t_spike = -3*ms
    HTCR_RND_cs.g_spike = 1.5*nS
    HTCR_RND_cs.tau = 5*ms

    #->RTCL
    RTCL_HTCL_gj.rgap = 300*Mohm
    RTCL_INL_cs.D_i = 1.07
    RTCL_INL_cs.g_syn = 3*nS
    RTCL_INL_cs.e_syn = -80*mV
    RTCL_INL_cs.alpha = 10.5/ms
    RTCL_INL_cs.beta = 0.166/ms
    RTCL_REL_cs.D_i = 1.07
    RTCL_REL_cs.g_syn = 3*nS
    RTCL_REL_cs.e_syn = -80*mV
    RTCL_REL_cs.alpha = 10.5/ms
    RTCL_REL_cs.beta = 0.166/ms
    RTCL_PYL_csa.D_i = 1.07
    RTCL_PYL_csa.g_syn = 4*nS
    RTCL_PYL_csa.e_syn = 0*mV
    RTCL_PYL_csa.alpha = 0.94/ms
    RTCL_PYL_csa.beta = 0.18/ms
    RTCL_PYL_csb.D_i = 1.07
    RTCL_PYL_csb.g_syn = 2*nS
    RTCL_PYL_csb.alpha = 1/ms
    RTCL_PYL_csb.beta = 0.0067/ms
    RTCL_RND_cs.t_spike = -3*ms
    RTCL_RND_cs.g_spike = 1.5*nS
    RTCL_RND_cs.tau = 5*ms
    #->RTCR
    RTCR_HTCR_gj.rgap = 300*Mohm
    RTCR_INR_cs.D_i = 1.07
    RTCR_INR_cs.g_syn = 3*nS
    RTCR_INR_cs.e_syn = -80*mV
    RTCR_INR_cs.alpha = 10.5/ms
    RTCR_INR_cs.beta = 0.166/ms
    RTCR_RER_cs.D_i = 1.07
    RTCR_RER_cs.g_syn = 3*nS
    RTCR_RER_cs.e_syn = -80*mV
    RTCR_RER_cs.alpha = 10.5/ms
    RTCR_RER_cs.beta = 0.166/ms
    RTCR_PYR_csa.D_i = 1.07
    RTCR_PYR_csa.g_syn = 4*nS
    RTCR_PYR_csa.e_syn = 0*mV
    RTCR_PYR_csa.alpha = 0.94/ms
    RTCR_PYR_csa.beta = 0.18/ms
    RTCR_PYR_csb.D_i = 1.07
    RTCR_PYR_csb.g_syn = 2*nS
    RTCR_PYR_csb.alpha = 1/ms
    RTCR_PYR_csb.beta = 0.0067/ms
    RTCR_RND_cs.t_spike = -3*ms
    RTCR_RND_cs.g_spike = 1.5*nS
    RTCR_RND_cs.tau = 5*ms

    #->INL
    INL_HTCL_csa.D_i = 1.07
    INL_HTCL_csa.g_syn = 6*nS
    INL_HTCL_csa.e_syn = 0*mV
    INL_HTCL_csa.alpha = 0.94/ms
    INL_HTCL_csa.beta = 0.18/ms
    INL_HTCL_csb.D_i = 1.07
    INL_HTCL_csb.g_syn = 3*nS
    INL_HTCL_csb.alpha = 1/ms
    INL_HTCL_csb.beta = 0.0067/ms
    INL_REL_cs.D_i = 1.07
    INL_REL_cs.g_syn = 1*nS
    INL_REL_cs.e_syn = -80*mV
    INL_REL_cs.alpha = 10.5/ms
    INL_REL_cs.beta = 0.166/ms
    INL_RND_cs.t_spike = -3*ms
    INL_RND_cs.g_spike = 1.5*nS
    INL_RND_cs.tau = 5*ms
    #->INR
    INR_HTCR_csa.D_i = 1.07
    INR_HTCR_csa.g_syn = 6*nS
    INR_HTCR_csa.e_syn = 0*mV
    INR_HTCR_csa.alpha = 0.94/ms
    INR_HTCR_csa.beta = 0.18/ms
    INR_HTCR_csb.D_i = 1.07
    INR_HTCR_csb.g_syn = 3*nS
    INR_HTCR_csb.alpha = 1/ms
    INR_HTCR_csb.beta = 0.0067/ms
    INR_RER_cs.D_i = 1.07
    INR_RER_cs.g_syn = 1*nS
    INR_RER_cs.e_syn = -80*mV
    INR_RER_cs.alpha = 10.5/ms 
    INR_RER_cs.beta = 0.166/ms
    INR_RND_cs.t_spike = -3*ms
    INR_RND_cs.g_spike = 1.5*nS
    INR_RND_cs.tau = 5*ms

    #->REL
    REL_REL_gj.rgap = 300*Mohm
    REL_HTCL_csa.D_i = 1.07
    REL_HTCL_csa.g_syn = 4*nS
    REL_HTCL_csa.e_syn = 0*mV
    REL_HTCL_csa.alpha = 0.94/ms
    REL_HTCL_csa.beta = 0.18/ms
    REL_HTCL_csb.D_i = 1.07 
    REL_HTCL_csb.g_syn = 2*nS
    REL_HTCL_csb.alpha = 1/ms
    REL_HTCL_csb.beta = 0.0067/ms
    REL_RTCL_csa.D_i = 1.07
    REL_RTCL_csa.g_syn = 4*nS
    REL_RTCL_csa.e_syn = 0*mV
    REL_RTCL_csa.alpha = 0.94/ms
    REL_RTCL_csa.beta = 0.18/ms
    REL_RTCL_csb.D_i = 1.07
    REL_RTCL_csb.g_syn = 2*nS
    REL_RTCL_csb.alpha = 1/ms
    REL_RTCL_csb.beta = 0.0067/ms
    REL_REL_cs.D_i = 1.07
    REL_REL_cs.g_syn = 1*nS
    REL_REL_cs.e_syn = -70*mV
    REL_REL_cs.alpha = 10.5/ms
    REL_REL_cs.beta = 0.166/ms
    REL_PYL_csa.D_i = 1.07
    REL_PYL_csa.g_syn = 4*nS
    REL_PYL_csa.e_syn = 0*mV
    REL_PYL_csa.alpha = 0.94/ms
    REL_PYL_csa.beta = 0.18/ms
    REL_PYL_csb.D_i = 1.07
    REL_PYL_csb.g_syn = 2*nS
    REL_PYL_csb.alpha = 1/ms
    REL_PYL_csb.beta = 0.0067/ms
    REL_RND_cs.t_spike = -3*ms
    REL_RND_cs.g_spike = 1.5*nS
    REL_RND_cs.tau = 5*ms
    #->RER
    RER_RER_gj.rgap = 300*Mohm
    RER_HTCR_csa.D_i = 1.07
    RER_HTCR_csa.g_syn = 4*nS
    RER_HTCR_csa.e_syn = 0*mV
    RER_HTCR_csa.alpha = 0.94/ms
    RER_HTCR_csa.beta = 0.18/ms
    RER_HTCR_csb.D_i = 1.07
    RER_HTCR_csb.g_syn = 2*nS
    RER_HTCR_csb.alpha = 1/ms
    RER_HTCR_csb.beta = 0.0067/ms
    RER_RTCR_csa.D_i = 1.07
    RER_RTCR_csa.g_syn = 4*nS
    RER_RTCR_csa.e_syn = 0*mV
    RER_RTCR_csa.alpha = 0.94/ms
    RER_RTCR_csa.beta = 0.18/ms
    RER_RTCR_csb.D_i = 1.07
    RER_RTCR_csb.g_syn = 2*nS
    RER_RTCR_csb.alpha = 1/ms
    RER_RTCR_csb.beta = 0.0067/ms
    RER_RER_cs.D_i = 1.07
    RER_RER_cs.g_syn = 1*nS
    RER_RER_cs.e_syn = -70*mV
    RER_RER_cs.alpha = 10.5/ms
    RER_RER_cs.beta = 0.166/ms
    RER_PYR_csa.D_i = 1.07
    RER_PYR_csa.g_syn = 4*nS
    RER_PYR_csa.e_syn = 0*mV
    RER_PYR_csa.alpha = 0.94/ms
    RER_PYR_csa.beta = 0.18/ms
    RER_PYR_csb.D_i = 1.07
    RER_PYR_csb.g_syn = 2*nS
    RER_PYR_csb.alpha = 1/ms
    RER_PYR_csb.beta = 0.0067/ms
    RER_RND_cs.t_spike = -3*ms
    RER_RND_cs.g_spike = 1.5*nS
    RER_RND_cs.tau = 5*ms

    #->PYL
    PYL_RTCL_cs.tau = 2*ms
    PYL_RTCL_cs.e_syn = 0*mV
    PYL_RTCL_cs.g_spike = 0.3*nS
    PYL_PYL_cs.tau = 2*ms
    PYL_PYL_cs.e_syn = 0*mV
    PYL_PYL_cs.g_spike = 0.3*nS
    PYL_FSL_cs.tau = 10*ms
    PYL_FSL_cs.e_syn = -70*mV
    PYL_FSL_cs.g_spike = 0.3*nS
    PYL_RND.std = 0.1*pA
    PYL_RND.randi = randn()
    #->PYR
    PYR_RTCR_cs.tau = 2*ms
    PYR_RTCR_cs.e_syn = 0*mV
    PYR_RTCR_cs.g_spike = 0.3*nS
    PYR_PYR_cs.tau = 2*ms
    PYR_PYR_cs.e_syn = 0*mV
    PYR_PYR_cs.g_spike = 0.3*nS
    PYR_FSR_cs.tau = 10*ms
    PYR_FSR_cs.e_syn = -70*mV
    PYR_FSR_cs.g_spike = 0.3*nS
    PYR_RND.std = 0.1*pA
    PYR_RND.randi = randn()

    #->FSL
    FSL_RTCL_cs.tau = 2*ms
    FSL_RTCL_cs.e_syn = 0*mV
    FSL_RTCL_cs.g_spike = 0.4*nS
    FSL_PYL_cs.tau = 2*ms
    FSL_PYL_cs.e_syn = 0*mV
    FSL_PYL_cs.g_spike = 0.4*nS
    FSL_FSL_cs.tau = 10*ms
    FSL_FSL_cs.e_syn = -70*mV
    FSL_FSL_cs.g_spike = 0.03*nS
    FSL_RND.std = 0.1*pA
    FSL_RND.randi = randn()
    #->FSR
    FSR_RTCR_cs.tau = 2*ms
    FSR_RTCR_cs.e_syn = 0*mV
    FSR_RTCR_cs.g_spike = 0.4*nS
    FSR_PYR_cs.tau = 2*ms
    FSR_PYR_cs.e_syn = 0*mV
    FSR_PYR_cs.g_spike = 0.4*nS
    FSR_FSR_cs.tau = 10*ms
    FSR_FSR_cs.e_syn = -70*mV
    FSR_FSR_cs.g_spike = 0.03*nS
    FSR_RND.std = 0.1*pA
    FSR_RND.randi = randn()

    if(connected):

        PYL_PYR_cs = Synapses(PYRg,PYLg,(PYFS_cs + iPYR_eqs),on_pre = prePYFS,
                method = 'rk4', dt = t_step, name = "PYL_PYR")
        PYL_PYR_cs.connect(p = '0.5')
        PYL_PYR_cs.tau = 2*ms
        PYL_PYR_cs.e_syn = 0*mV
        PYL_PYR_cs.g_spike = 0.3*nS

        PYR_PYL_cs = Synapses(PYLg,PYRg,(PYFS_cs + iPYL_eqs),on_pre = prePYFS,
                method = 'rk4', dt = t_step, name = "PYR_PYL")
        PYR_PYL_cs.connect(p = '0.5')
        PYR_PYL_cs.tau = 2*ms
        PYR_PYL_cs.e_syn = 0*mV
        PYR_PYL_cs.g_spike = 0.3*nS

        PYL_FSR_cs = Synapses(FSRg,PYLg,(PYFS_cs + iFSR_eqs),on_pre = prePYFS,
                method = 'rk4', dt = t_step, name = "PYL_FSR")
        PYL_FSR_cs.connect('abs(i*(nPY-1.0)/(nFS-1.0) - j) < nPY*0.2', p = '0.8')
        PYL_FSR_cs.tau = 10*ms
        PYL_FSR_cs.e_syn = -70*mV
        PYL_FSR_cs.g_spike = 0.3*nS

        PYR_FSL_cs = Synapses(FSLg,PYRg,(PYFS_cs + iFSL_eqs),on_pre = prePYFS,
                method = 'rk4', dt = t_step, name = "PYR_FSL")
        PYR_FSL_cs.connect('abs(i*(nPY-1.0)/(nFS-1.0) - j) < nPY*0.2', p = '0.8')
        PYR_FSL_cs.tau = 10*ms
        PYR_FSL_cs.e_syn = -70*mV
        PYR_FSL_cs.g_spike = 0.3*nS

        FSR_PYL_cs = Synapses(PYLg,FSRg,(PYFS_cs + iPYL_eqs),on_pre = prePYFS,
                method = 'rk4', dt = t_step, name = "FSR_PYL")
        FSR_PYL_cs.connect(i = PYL_FSR_cs.j[:], j = PYL_FSR_cs.i[:])
        FSR_PYL_cs.tau = 2*ms
        FSR_PYL_cs.e_syn = 0*mV
        FSR_PYL_cs.g_spike = 0.4*nS

        FSL_PYR_cs = Synapses(PYRg,FSLg,(PYFS_cs + iPYR_eqs),on_pre = prePYFS,
                method = 'rk4', dt = t_step, name = "FSL_PYR")
        FSL_PYR_cs.connect(i = PYR_FSL_cs.j[:], j = PYR_FSL_cs.i[:])
        FSL_PYR_cs.tau = 2*ms
        FSL_PYR_cs.e_syn = 0*mV
        FSL_PYR_cs.g_spike = 0.4*nS


        PYL_RTCR_cs = Synapses(RTCRg,PYLg,(PYFS_cs + iRTCRa_eqs),on_pre = prePYFS,
                method = 'rk4', dt = t_step, name = "PYL_RTCR_cs")
        PYL_RTCR_cs.connect(p = '0.04')
        PYL_RTCR_cs.tau = 2*ms
        PYL_RTCR_cs.e_syn = 0*mV
        PYL_RTCR_cs.g_spike = 0.3*nS


        PYR_RTCL_cs = Synapses(RTCLg,PYRg,(PYFS_cs + iRTCLa_eqs),on_pre = prePYFS, 
                method = 'rk4', dt = t_step, name = "PYR_RTCL_cs")
        PYR_RTCL_cs.connect(p = '0.04')
        PYR_RTCL_cs.tau = 2*ms
        PYR_RTCL_cs.e_syn = 0*mV
        PYR_RTCL_cs.g_spike = 0.3*nS



        FSR_RTCL_cs = Synapses(RTCLg,FSRg,(PYFS_cs + iRTCLa_eqs),on_pre = prePYFS, 
                method = 'rk4', dt = t_step, name = "FSR_RTCL_cs")
        FSR_RTCL_cs.connect(p = '0.02')
        FSR_RTCL_cs.tau = 2*ms
        FSR_RTCL_cs.e_syn = 0*mV
        FSR_RTCL_cs.g_spike = 0.4*nS


        FSL_RTCR_cs = Synapses(RTCRg,FSLg,(PYFS_cs + iRTCRa_eqs),on_pre = prePYFS, 
                method = 'rk4', dt = t_step, name = "FSL_RTCR_cs")
        FSL_RTCR_cs.connect(p = '0.02')
        FSL_RTCR_cs.tau = 2*ms
        FSL_RTCR_cs.e_syn = 0*mV
        FSL_RTCR_cs.g_spike = 0.4*nS


        RTCR_PYL_csa = Synapses(PYLg, RTCRg, (CSa + iPYLa_eqs), on_pre = preCS, 
                method = 'rk4', dt = t_step, name = "RTCR_PYL_csa")
        RTCR_PYL_csa.connect(p = '0.23')
        RTCR_PYL_csa.D_i = 1.07
        RTCR_PYL_csa.g_syn = 4*nS
        RTCR_PYL_csa.e_syn = 0*mV
        RTCR_PYL_csa.alpha = 0.94/ms
        RTCR_PYL_csa.beta = 0.18/ms

        RTCL_PYR_csa = Synapses(PYRg, RTCLg, (CSa + iPYRa_eqs), on_pre = preCS, 
                method = 'rk4', dt = t_step, name = "RTCL_PYR_csa")
        RTCL_PYR_csa.connect(p = '0.23')
        RTCL_PYR_csa.D_i = 1.07
        RTCL_PYR_csa.g_syn = 4*nS
        RTCL_PYR_csa.e_syn = 0*mV
        RTCL_PYR_csa.alpha = 0.94/ms
        RTCL_PYR_csa.beta = 0.18/ms

        RTCR_PYL_csb = Synapses(PYLg, RTCRg, (CSb + iPYLb_eqs), on_pre = preCS, 
                method = 'rk4', dt = t_step, name = "RTCR_PYL_csb")
        RTCR_PYL_csb.connect(i = RTCR_PYL_csa.i[:], j = RTCR_PYL_csa.j[:])
        RTCR_PYL_csb.D_i = 1.07
        RTCR_PYL_csb.g_syn = 2*nS
        RTCR_PYL_csb.alpha = 1/ms
        RTCR_PYL_csb.beta = 0.0067/ms

        RTCL_PYR_csb = Synapses(PYRg, RTCLg, (CSb + iPYRb_eqs), on_pre = preCS, 
                method = 'rk4', dt = t_step, name = "RTCL_PYR_csb")
        RTCL_PYR_csb.connect(i = RTCL_PYR_csa.i[:], j = RTCL_PYR_csa.j[:])
        RTCL_PYR_csb.D_i = 1.07
        RTCL_PYR_csb.g_syn = 2*nS
        RTCL_PYR_csb.alpha = 1/ms
        RTCL_PYR_csb.beta = 0.0067/ms

        RER_PYL_csa = Synapses(PYLg, RERg, (CSa + iPYLa_eqs), on_pre = preCS, 
                method = 'rk4', dt = t_step, name = "RER_PYL_csa")
        RER_PYL_csa.connect(p = '0.3')
        RER_PYL_csa.D_i = 1.07
        RER_PYL_csa.g_syn = 4*nS
        RER_PYL_csa.e_syn = 0*mV
        RER_PYL_csa.alpha = 0.94/ms
        RER_PYL_csa.beta = 0.18/ms

        REL_PYR_csa = Synapses(PYRg, RELg, (CSa + iPYRa_eqs), on_pre = preCS, 
                method = 'rk4', dt = t_step, name = "REL_PYR_csa")
        REL_PYR_csa.connect(p = '0.3')
        REL_PYR_csa.D_i = 1.07
        REL_PYR_csa.g_syn = 4*nS
        REL_PYR_csa.e_syn = 0*mV
        REL_PYR_csa.alpha = 0.94/ms
        REL_PYR_csa.beta = 0.18/ms


        RER_PYL_csb = Synapses(PYLg, RERg, (CSb + iPYLb_eqs), on_pre = preCS, 
                method = 'rk4', dt = t_step, name = "RER_PYL_csb")
        RER_PYL_csb.connect(i = RER_PYL_csa.i[:], j = RER_PYL_csa.j[:])
        RER_PYL_csb.D_i = 1.07
        RER_PYL_csb.g_syn = 2*nS
        RER_PYL_csb.alpha = 1/ms
        RER_PYL_csb.beta = 0.0067/ms

        REL_PYR_csb = Synapses(PYRg, RELg, (CSb + iPYRb_eqs), on_pre = preCS, 
                method = 'rk4', dt = t_step, name = "REL_PYR_csb")
        REL_PYR_csb.connect(i = REL_PYR_csa.i[:], j = REL_PYR_csa.j[:])
        REL_PYR_csb.D_i = 1.07
        REL_PYR_csb.g_syn = 2*nS
        REL_PYR_csb.alpha = 1/ms
        REL_PYR_csb.beta = 0.0067/ms

    #State Monitors
    TCR_volt = StateMonitor(TCRg,'v', record = True)
    TCL_volt = StateMonitor(TCLg,'v', record = True)

    PYL_volt = StateMonitor(PYLg, 'v', record = True)
    PYR_volt = StateMonitor(PYRg, 'v', record = True)

    HTCR_volt = StateMonitor(HTCRg,'v', record = True)
    HTCL_volt = StateMonitor(HTCLg,'v', record = True)

    RTCR_volt = StateMonitor(RTCRg,'v', record = True)
    RTCL_volt = StateMonitor(RTCLg,'v', record = True)

    FSR_volt = StateMonitor(FSRg,'v', record = True)
    FSL_volt = StateMonitor(FSLg,'v', record = True)

    RER_volt = StateMonitor(RERg,'v', record = True)
    REL_volt = StateMonitor(RELg,'v', record = True)

    INR_volt = StateMonitor(INRg,'v', record = True)
    INL_volt = StateMonitor(INLg,'v', record = True)

    TCR_volt = StateMonitor(TCRg,'v', record = True)
    TCL_volt = StateMonitor(TCLg,'v', record = True)


    #Spike Monitors
    TC_spike = SpikeMonitor(TCLg, record = True)
    IN_spike = SpikeMonitor(INLg, record = True)
    RE_spike = SpikeMonitor(RELg, record = True)
    #PYFS_spike = SpikeMonitor(PYFSLg, record = True)

    Sgn = TimedArray(0*sin(20*pi*t_d)*25*pA, dt = t_step)

    print("Setup complete.")
    run(duration)

    t = int(time.time())

    base = "output/" + str(t) + "/"
    os.mkdir(base)

    if(connected):
        stem = base + "connected/"
    else:
        stem = base + "unconnected/"

    os.mkdir(stem)


    np.save(stem + "PYL_time.npy", PYL_volt.t)
    np.save(stem + "PYL_volt.npy", PYL_volt.v)
    np.save(stem + "PYR_time.npy", PYR_volt.t)
    np.save(stem + "PYR_volt.npy", PYR_volt.v)

    np.save(stem + "INL_time.npy", INL_volt.t)
    np.save(stem + "INL_volt.npy", INL_volt.v)
    np.save(stem + "INR_time.npy", INR_volt.t)
    np.save(stem + "INR_volt.npy", INR_volt.v)

    np.save(stem + "REL_time.npy", REL_volt.t)
    np.save(stem + "REL_volt.npy", REL_volt.v)
    np.save(stem + "RER_time.npy", RER_volt.t)
    np.save(stem + "RER_volt.npy", RER_volt.v)

    np.save(stem + "RTCL_time.npy", RTCL_volt.t)
    np.save(stem + "RTCL_volt.npy", RTCL_volt.v)
    np.save(stem + "RTCR_time.npy", RTCR_volt.t)
    np.save(stem + "RTCR_volt.npy", RTCR_volt.v)

    np.save(stem + "HTCL_time.npy", HTCL_volt.t)
    np.save(stem + "HTCL_volt.npy", HTCL_volt.v)
    np.save(stem + "HTCR_time.npy", HTCR_volt.t)
    np.save(stem + "HTCR_volt.npy", HTCR_volt.v)

    np.save(stem + "TCL_time.npy", TCL_volt.t)
    np.save(stem + "TCL_volt.npy", TCL_volt.v)
    np.save(stem + "TCR_time.npy", TCR_volt.t)
    np.save(stem + "TCR_volt.npy", TCR_volt.v)

    np.save(stem + "FSL_time.npy", FSL_volt.t)
    np.save(stem + "FSL_volt.npy", FSL_volt.v)
    np.save(stem + "FSR_time.npy", FSR_volt.t)
    np.save(stem + "FSR_volt.npy", FSR_volt.v)



if __name__ == "__main__":
    main()
