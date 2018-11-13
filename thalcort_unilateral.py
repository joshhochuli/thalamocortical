from brian2 import *
seed(555)

#Global Variables Across Thalamus Cells
#Capacitance
C = 1*uF/cm**2

#Leak Variables
gL = 0.01*mS/cm**2		#average g_L value +- 0.0025, uniformly distributed	
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
e_K = -90*mV 			#-80*mV in INaK, -90*mV in paper

#iH Variables	HTC
e_H = -43*mV

#iAHP Variables		
e_AHP = -90*mV			#e_K

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
duration = 0.05*second
t_step = 0.02*ms

#tACs Signal Array
t_d = t_step*numpy.arange(duration/t_step)/second

#Modular Equation Expressions
#Main HH Equation, HTC, RTC
main_TC = '''
dv/dt=(-g_L*(v-e_L)-g_KL*(v-e_KL)-iCa-iCal-iCat-iNa-iK-iH-iAHP-iCan-iSyn+iStm)/C:volt
e_L																				:volt
g_KL																			:siemens/meter**2
g_L																				:siemens/meter**2
'''

#Main HH Equation, IN
main_IN = '''
dv/dt=(-g_L*(v-e_L)-g_KL*(v-e_KL)-iCa-iNa-iK-iH-iAHP-iCan-iSyn+iStm)/C			:volt
e_L																				:volt
g_KL																			:siemens/meter**2
g_L																				:siemens/meter**2
'''

#Main HH Equation, RE
main_RE = '''
dv/dt=(-g_L*(v-e_L)-g_KL*(v-e_KL)-iCa-iNa-iK-iAHP-iCan-iSyn+iStm)/C				:volt
e_L																				:volt
g_KL																			:siemens/meter**2
g_L																				:siemens/meter**2
'''

#Ca_TC, HTC, RTC
Ca_TC ='''
tau_Ca = 10*ms																	:second					
drive = -(iCa+iCat+iCal)/(z*F*w)*(-(iCa+iCat+iCal)/(z*F*w)>=0*mM/second)		:mM/second	
dca/dt = drive + (ca_Rest - ca)/tau_Ca											:mM
'''

#Ca_RE_IN, IN, RE
Ca_RE_IN ='''
tau_Ca																			:second					
drive = -iCa/(z*F*w)*(-iCa/(z*F*w)>=0*mM/second)								:mM/second	
dca/dt = drive + (ca_Rest - ca)/tau_Ca											:mM
'''

#iCa, current 0, T0_TC
i_Ca = '''
g_Ca																			:siemens/meter**2				
e_Ca = eca0*log(ca_0/ca) 														:volt
m0_inf = 1/(1+exp(-(v+34*mV)/(6.2*mV))) 										:1		
tau_m0 = 0.612*ms+1*ms/(exp(-(v+107*mV)/(16.7*mV))+exp((v-8.2*mV)/(18.2*mV)))	:second	
dm0/dt = 4.6*(m0_inf-m0)/tau_m0													:1
h0_inf = 1/(1+exp((v+58*mV)/(4*mV)))											:1
tau_h0a = (1*ms)*(exp((v+442*mV)/(66.6*mV)))*(v < -55*mV)						:second
tau_h0b = 28*ms+(1*ms)*(exp(-(v-3*mV)/(10.5*mV)))*(v >= -55*mV)					:second
tau_h0 = tau_h0a + tau_h0b														:second		
dh0/dt = 3.7*(h0_inf-h0)/tau_h0	 												:1
iCa = g_Ca*m0**2*h0*(v-e_Ca) 													:amp/meter**2
'''

#iCa_RE, current 0
i_Ca_RE = '''
g_Ca = 1.3*mS/cm**2																:siemens/meter**2 
e_Ca = eca0*log(ca_0/ca) 														:volt
m0_inf = 1/(1+exp(-(v+55*mV)/(7.4*mV))) 										:1		
tau_m0 = 3*ms + 1*ms/(exp(-(v+105*mV)/(15*mV))+exp((v+30*mV)/(10*mV)))			:second	
dm0/dt = 6.9*(m0_inf-m0)/tau_m0													:1
h0_inf = 1/(1+exp((v+83*mV)/(5*mV)))											:1	
tau_h0 = 85*ms + 1*ms/(exp(-(v+410*mV)/(50*mV))+exp((v+51*mV)/(4*mV)))			:second		
dh0/dt = 3.7*(h0_inf-h0)/tau_h0	 												:1
iCa = g_Ca*m0**2*h0*(v-e_Ca) 													:amp/meter**2
'''

#iCat, current 1, IT0_TC
i_Cat = '''
m1_inf = 1/(1+exp(-(v+62*mV)/(6.2*mV))) 										:1		
tau_m1 = 0.612*ms+1*ms/(exp(-(v+135*mV)/(16.7*mV))+exp((v+19.8*mV)/(18.2*mV)))	:second	
dm1/dt = 4.6*(m1_inf-m1)/tau_m1													:1
h1_inf = 1/(1+exp((v+86*mV)/(4*mV)))											:1
tau_h1a = (1*ms)*(exp((v+470*mV)/(66.6*mV)))*(v < -83*mV)						:second
tau_h1b = 28*ms+(1*ms)*(exp(-(v+25*mV)/(10.5*mV)))*(v >= -83*mV)				:second
tau_h1 = tau_h1a + tau_h1b														:second		
dh1/dt = 3.7*(h1_inf-h1)/tau_h1	 												:1
iCat = g_Cat*m1**2*h1*(v-e_Ca) 													:amp/meter**2
'''

#iCal, current 2, ICaL_TC
i_Cal = '''
g_Cal 																			:siemens/meter**2 
m2_inf = 1/(1+exp(-(v+10*mV)/(4*mV))) 											:1		
tau_m2 = 0.4*ms + 0.7*ms/(exp(-(v+5*mV)/(15*mV))+exp((v+5*mV)/(15*mV)))			:second	
dm2/dt = 4.6*(m2_inf-m2)/tau_m2													:1
h2_inf = 1/(1+exp((v+25*mV)/(2*mV)))											:1
tau_h2 = 300*ms+(100*ms)/(exp(-(v+40*mV)/(9.5*mV))+exp((v+40*mV)/(9.5*mV)))		:second	
dh2/dt = 3.7*(h2_inf-h2)/tau_h2	 												:1
iCal = g_Cal*m2**2*h2*(v-e_Ca) 													:amp/meter**2
''' 

#iNa, current 3, INaK
i_Na = '''
v_SH																			:volt
a_m3 = (0.32*ms**-1*mV**-1)*(v-v_SH-13*mV)/(1-exp(-(v-v_SH-13*mV)/(4*mV)))		:Hz 	
b_m3 = (-0.28*mV**-1*ms**-1)*(v-v_SH-40*mV)/(1-exp((v-v_SH-40*mV)/(5*mV)))		:Hz	
dm3/dt = (a_m3*(1-m3)-b_m3*m3)													:1	
a_h3 = (0.128*ms**-1)*exp(-(v-v_SH-17*mV)/(18*mV))								:Hz		
b_h3 = (4*ms**-1)/(1+exp(-(v-v_SH-40*mV)/(5*mV)))								:Hz	
dh3/dt = (a_h3*(1-h3)-b_h3*h3)													:1		
iNa = g_Na*m3**3*h3*(v-e_Na) 													:amp/meter**2
'''

#iK, current 4, INaK
i_K = '''
tm																				:1					
a_m4 = tm*(0.032*ms**-1*mV**-1)*(v-v_SH-15*mV)/(1-exp(-(v-v_SH-15*mV)/(5*mV)))	:Hz	
b_m4 = tm*(0.5*ms**-1)*exp(-(v-v_SH-10*mV)/(40*mV))								:Hz		
dm4/dt = (a_m4*(1-m4)-b_m4*m4)													:1		
iK = g_K*m4**4*(v-e_K) 															:amp/meter**2
'''

#iH, current 5, Ih_TC
i_H = '''
g_H																				:siemens/meter**2 
m5_inf = 1/(1+exp((v+75*mV)/(5.5*mV))) 											:1
tau_m5 = 1*ms/(exp(-0.086/mV*v-14.59)+exp(0.0701/mV*v-1.87))					:second
dm5/dt = (m5_inf-m5)/(tau_m5) 													:1
iH = g_H*m5*(v-e_H) 															:amp/meter**2
'''

#iAHP, current 6, Iahp2
i_AHP = '''
g_AHP																			:siemens/meter**2
m6_inf = 48*ca**2/(48*ca**2+0.09*mM**2)											:1
tau_m6 = (1*ms*mM**2)/(48*ca**2+0.09*mM**2)										:second
dm6/dt = (m6_inf - m6)/(tau_m6) 												:1
iAHP = g_AHP*m6*(v-e_AHP) 														:amp/meter**2
'''

#iCan, current 7, Ican_TC
i_Can = '''
g_Can																			:siemens/meter**2
mCa = ca/(0.2*mM+ca) 															:1
m7_inf = 1/(1+exp(-(v+43*mV)/(5.2*mV)))											:1
tau_m7 = 1.6*ms+(2.7*ms)/(exp(-(v+55*mV)/(15*mV))+exp((v+55*mV)/(15*mV)))		:second
dm7/dt = (m7_inf - m7)/(tau_m7) 												:1
iCan = g_Can*mCa*m7*(v-e_Can) 													:amp/meter**2
'''

#iStim, General, Time-Dependent Stimulus
i_Stim = '''
n_Stim 																			:meter**2
mag																				:amp
iStm = (t > 1*second)*(t < 2*second)*mag/n_Stim									:amp/meter**2
'''

#iSyn, Synapse Currents, RE
i_Syn = '''
igj																				:amp
igjHR																			:amp
iHTCa																			:amp
iHTCb																			:amp
iRTCa																			:amp
iRTCb																			:amp
iIN																				:amp
iRE																				:amp
iPYa																			:amp
iPYb																			:amp
iRND																			:amp
iGLUT = iHTCa + iHTCb + iRTCa + iRTCb + iPYa + iPYb + iRND						:amp 
iGABA = iIN + iRE																:amp
iSyn = (igj + igjHR + iGLUT + iGABA)/n_Stim										:amp/meter**2
'''

#Location equations
loc = '''
x 																				:1
y 																				:1
gjb																				:1
'''

#Main Izhikevich Equation, PY FS
main_PY_FS = '''
dv/dt = (k*(v - vr)*(v - vt) - u - iSyn + iStm + sg*Sgn(t))/Cm					:volt
Cm																				:farad
k																				:siemens/volt
vt																				:volt
vr																				:volt
du/dt = a*(b*((v - vr)/mV)**un*(v > vb) - u)									:amp
a 																				:Hz
b 																				:amp
vb																				:volt
un																				:1
vth																				:volt
c 																				:volt
d 																				:amp
sg																				:1
iRTCa																			:amp
iPY																				:amp
iFS																				:amp
iRND																			:amp
iStm																			:amp
iSyn = iRTCa + iPY + iFS + iRND													:amp
LFP = abs(iRTCa + iPY) + abs(iFS)												:amp
'''

#Cell Group Declarations
#TC Equations
TC_eqs = (main_TC + Ca_TC + i_Ca + i_Cat + i_Cal + i_Na + i_K + i_H + i_AHP + i_Can + i_Syn + i_Stim + loc)

#IN Equations
IN_eqs = (main_IN + Ca_RE_IN + i_Ca + i_Na + i_K + i_H + i_AHP + i_Can + i_Syn + i_Stim)

#RE Equations
RE_eqs = (main_RE + Ca_RE_IN + i_Ca_RE + i_Na + i_K + i_AHP + i_Can + i_Syn + i_Stim + loc)

#Declare TC Group
TCg = NeuronGroup(nHTC**2+nRTC**2,TC_eqs,method = 'rk4',dt = t_step,threshold = 'v > 0*mV',refractory = 3*ms)
TCg.e_L = -70*mV
TCg.g_KL = 0*mS/cm**2
TCg.g_L = gL + (rand(nHTC**2+nRTC**2) - 0.5)*0.005*mS/cm**2
TCg.tm  = 0.25
TCg.g_H = 0.01*mS/cm**2
TCg.n_Stim = 2.9e-4*cm**2

#Declare HTC Group
HTCg = TCg[:nHTC**2]
HTCg.g_Ca = 3*mS/cm**2
HTCg.g_Cal = 0.5*mS/cm**2
HTCg.v_SH = -30*mV
HTCg.g_AHP = 0.3*mS/cm**2 
HTCg.g_Can = 0.5*mS/cm**2
HTCg.x = 'i%nHTC'
HTCg.y = '(i - x)/nHTC'

#Declare RTC Group
RTCg = TCg[nHTC**2:]
RTCg.g_Ca = 0.6*mS/cm**2
RTCg.g_Cal = 0.3*mS/cm**2
RTCg.v_SH = -40*mV
RTCg.g_AHP = 0.1*mS/cm**2 
RTCg.g_Can = 0.6*mS/cm**2
RTCg.x = 'i%nRTC'
RTCg.y = '(i - x)/nRTC'
RTCg.gjb = rand(nRTC*nRTC) < 0.2

#Declare IN Group
INg = NeuronGroup(nIN*nIN,IN_eqs,method = 'rk4',dt = t_step,threshold = 'v > 0*mV',refractory = 3*ms)
INg.e_L = -60*mV
INg.g_KL = 0.02*mS/cm**2
INg.g_L = gL + (rand(nIN*nIN) - 0.5)*0.005*mS/cm**2
INg.tau_Ca = 10*ms
INg.g_Ca = 2.5*mS/cm**2
INg.v_SH = -30*mV
INg.tm  = 0.25
INg.g_H = 0.05*mS/cm**2 
INg.g_AHP = 0.2*mS/cm**2 
INg.g_Can = 0.1*mS/cm**2
INg.n_Stim = 1.7e-4*cm**2

#Declare RE Group
REg = NeuronGroup(nRE*nRE,RE_eqs,method = 'rk4',dt = t_step,threshold = 'v > 0*mV',refractory = 3*ms)
REg.e_L = -60*mV
REg.g_KL = 0.01*mS/cm**2
REg.g_L = gL + (rand(nRE*nRE) - 0.5)*0.005*mS/cm**2
REg.tau_Ca = 100*ms
REg.v_SH = -40*mV
REg.tm  = 1
REg.g_AHP = 0.2*mS/cm**2 
REg.g_Can = 0.2*mS/cm**2
REg.n_Stim = 1.45e-4*cm**2
REg.x = 'i%nRE'
REg.y = '(i - x)/nRE'
REg.gjb = rand(nRE*nRE) < 0.2

#Declare PY FS Group
PYFSg = NeuronGroup(nPY+nFS, main_PY_FS, method = 'rk4', dt = t_step, threshold = 'v>= vth', reset = '''
v = c
u += d
''')

#PY Group Declaration
PYg = PYFSg[:nPY]
PYg.Cm = (100 + 0.1*randn(nPY))*pF
PYg.k = 0.7*pA/(mvolt*mvolt)
PYg.vr = -60*mV + .1*randn(nPY)*mV
PYg.vt = -40*mV + .1*randn(nPY)*mV
PYg.a = (0.03 + 0.001*randn(nPY))*kHz
PYg.b = (-2 + 0.01*randn(nPY))*pA
PYg.vb = -200*volt
PYg.un = 1
PYg.vth = 35*mV + .1*randn(nPY)*mV
PYg.c = -50*mV + 0.1*randn(nPY)*mV
PYg.d = (100 + 0.1*randn(nPY))*pA
PYg.v = -60*mV + .1*randn(nPY)*mV
PYg.u = PYg.b*(PYg.v-PYg.vr)/mV
PYg.iStm = 79*pA
PYg.sg = 1

#FS Group Declaration
FSg = PYFSg[nPY:]
FSg.Cm = (20 + 0.1*randn(nFS))*pF
FSg.k = (1 + 0.01*randn(nFS))*pA/(mvolt*mvolt)
FSg.vr = -55*mV + 0.1*randn(nFS)*mV
FSg.vt = -40*mV + 0.1*randn(nFS)*mV
FSg.a = (0.2 + 0.001*randn(nFS))*kHz
FSg.b = (0.025 + 0.0001*randn(nFS))*pA
FSg.vb = FSg.vr
FSg.un = 3
FSg.vth = 25*mV + 0.1*randn(nFS)*mV
FSg.c = -45*mV + 0.1*randn(nFS)*mV
FSg.d = (0 + 0.001*randn(nFS))*pA
FSg.v = -55*mV + 0.1*randn(nFS)*mV
FSg.u = (0 + 0.1*randn(nFS))*pA
FSg.iStm = 60*pA

#Constant Spiking Neuron for Noise
RNDg = NeuronGroup(1, model = ' v = 0 :1', method = 'rk4', dt = t_step, threshold = 'v < 1')


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
HTCg.set_states(int_HTCg)

#RTC Initial Conditions
int_RTCg = {'v':[v0],'ca':[ca0],'m0':[m00],'h0':[h00],'m1':[m10],'h1':[h10],'m2':[m20],'h2':[h20],'m3':[m30_b],'h3':[h30_b],'m4':[m40_b],'m5':[m50],'m6':[m60],'m7':[m70]}
RTCg.set_states(int_RTCg)

#IN Initial Conditions
int_INg = {'v':[v0],'ca':[ca0],'m0':[m00],'h0':[h00],'m3':[m30_a],'h3':[h30_a],'m4':[m40_a],'m5':[m50],'m6':[m60],'m7':[m70]}
INg.set_states(int_INg)

#RE Initial Conditions
int_REg = {'v':[v0],'ca':[ca0],'m0':[m00_RE],'h0':[h00_RE],'m3':[m30_b],'h3':[h30_b],'m4':[m40_b],'m6':[m60],'m7':[m70]}
REg.set_states(int_REg)

#Modular Synapse Equations, Biophysical
#Gap Junction
GJ = '''
rgap 										 									:ohm
iGJ = (v_post - v_pre)/rgap														:amp
'''

#Chemical Synapse, GABA_A, AMPA
CSa = '''
D_i 																			:1
t_spike																			:second
g_syn 																			:siemens
e_syn																			:volt
alpha 																			:Hz
beta 																			:Hz
T = 0.5*(t-t_spike < 2.3*ms)*(t-t_spike > 2*ms)*(t_spike > 0*ms)				:1
D = 1 - (1 - D_i*(1 - 0.07))*exp(-(t - t_spike)/700/ms)							:1
ds/dt = T*alpha*(1-s) - beta*s 													:1 (clock-driven)
iCS = D*g_syn*s*(v_post - e_syn)												:amp
'''

#Chemical Synapse, NMDA
CSb = '''
D_i 																			:1
t_spike																			:second
g_syn 																			:siemens
alpha 																			:Hz
beta 																			:Hz
T = 0.5*(t-t_spike < 2.3*ms)*(t-t_spike > 2*ms)*(t_spike > 0*ms)				:1
ds/dt = T*alpha*(1-s) - beta*s 													:1 (clock-driven)
D = 1 - (1 - D_i*(1 - 0.07))*exp(-(t - t_spike)/700/ms)							:1
B = 1/(1+exp(-(v_post + 25*mV)/12.5/mV))										:1
iCS = D*B*g_syn*s*v_post														:amp
'''

preCS = '''
t_spike = t
D_i = D
'''

#Random Spiking Synapse, HTC, RTC, IN, RE
iRND = '''
tau 																			:second
t_spike																			:second
g_spike																			:siemens
d 																				:1
g																				:siemens
iRND_post = g*v_post															:amp (summed)
'''
preRND = '''
d = (rand() < 0.003)*(t-t_spike > 3*ms)
t_spike = t*d + t_spike*(1-d)
g = (g + g_spike*d)*exp(-(t-t_spike)/tau)
'''

#->Thalamus Equation Expressions
igj_eqs = '''
igj_post = iGJ																	:amp (summed)
'''
igjHR_eqs = '''
igjHR_post = iGJ																:amp (summed)
'''
iHTCa_eqs = '''
iHTCa_post = iCS																:amp (summed)
'''
iHTCb_eqs = '''
iHTCb_post = iCS																:amp (summed)
'''
iRTCa_eqs = '''
iRTCa_post = iCS																:amp (summed)
'''
iRTCb_eqs = '''
iRTCb_post = iCS																:amp (summed)
'''
iIN_eqs = '''
iIN_post = iCS																	:amp (summed)
'''
iRE_eqs = '''
iRE_post = iCS																	:amp (summed)
'''
iPYa_eqs = '''
iPYa_post = iCS																	:amp (summed)
'''
iPYb_eqs = '''
iPYb_post = iCS																	:amp (summed)
'''

#Modular Synapse Equations, Izhikevich
#Chemical Synapse, GABA_A, AMDA
PYFS_cs = '''
tau 																			:second 
e_syn																			:volt
g_spike																			:siemens
dg/dt = -g/tau																	:siemens (clock-driven)
iCS = g*(v_post - e_syn)														:amp
'''
prePYFS = '''
g += g_spike
'''

#Random Spiking Synapse, PY, FS
RND_PYFS = '''
std 																			:amp
t_spike																			:second
randi																			:1
iRND_post = std*randi															:amp(summed)
'''
preRND_PYFS = '''
randi = randi*(t-t_spike<0.5*ms)+randn()*(1-(t-t_spike<0.5*ms))
t_spike = t_spike*(t-t_spike<0.5*ms)+t*(1-(t-t_spike<0.5*ms))
'''

#->PYFS Equation Expressions
iPY_eqs = '''
iPY_post = iCS																	:amp (summed)
'''
iFS_eqs = '''
iFS_post = iCS																	:amp (summed)
'''

#Synapse Group Declarations

#->HTC Synapses
HTC_HTC_gj = Synapses(HTCg,HTCg,(GJ + igj_eqs))
HTC_RTC_gj = Synapses(RTCg,HTCg,(GJ + igjHR_eqs))
HTC_RE_cs = Synapses(REg,HTCg,(CSa + iRE_eqs),on_pre = preCS,method = 'rk4',dt = t_step)
HTC_RND_cs = Synapses(RNDg, HTCg, iRND, on_pre = preRND, method = 'rk4', dt = t_step)

#->RTC Synapses
RTC_HTC_gj = Synapses(HTCg,RTCg,(GJ + igjHR_eqs))
RTC_IN_cs = Synapses(INg,RTCg,(CSa + iIN_eqs),on_pre = preCS, method = 'rk4', dt = t_step)
RTC_RE_cs = Synapses(REg,RTCg,(CSa + iRE_eqs),on_pre = preCS, method = 'rk4', dt = t_step)
RTC_PY_csa = Synapses(PYg, RTCg, (CSa + iPYa_eqs), on_pre = preCS, method = 'rk4', dt = t_step)
RTC_PY_csb = Synapses(PYg, RTCg, (CSb + iPYb_eqs), on_pre = preCS, method = 'rk4', dt = t_step)
RTC_RND_cs = Synapses(RNDg, RTCg, iRND, on_pre = preRND, method = 'rk4', dt = t_step)

#->IN Synapses
IN_HTC_csa = Synapses(HTCg,INg,(CSa + iHTCa_eqs),on_pre = preCS, method = 'rk4', dt = t_step)
IN_HTC_csb = Synapses(HTCg,INg,(CSb + iHTCb_eqs),on_pre = preCS, method = 'rk4', dt = t_step)
IN_RE_cs = Synapses(REg,INg,(CSa + iRE_eqs),on_pre = preCS, method = 'rk4', dt = t_step)
IN_RND_cs = Synapses(RNDg, INg, iRND, on_pre = preRND, method = 'rk4', dt = t_step)

#->RE Synapses
RE_RE_gj = Synapses(REg,REg,(GJ + igj_eqs))
RE_HTC_csa = Synapses(HTCg,REg,(CSa + iHTCa_eqs),on_pre = preCS, method = 'rk4', dt = t_step)
RE_HTC_csb = Synapses(HTCg,REg,(CSb + iHTCb_eqs),on_pre = preCS, method = 'rk4', dt = t_step)
RE_RTC_csa = Synapses(RTCg,REg,(CSa + iRTCa_eqs),on_pre = preCS, method = 'rk4', dt = t_step)
RE_RTC_csb = Synapses(RTCg,REg,(CSb + iRTCb_eqs),on_pre = preCS, method = 'rk4', dt = t_step)
RE_RE_cs = Synapses(REg,REg,(CSa + iRE_eqs),on_pre = preCS, method = 'rk4', dt = t_step)
RE_PY_csa = Synapses(PYg, REg, (CSa + iPYa_eqs), on_pre = preCS, method = 'rk4', dt = t_step)
RE_PY_csb = Synapses(PYg, REg, (CSb + iPYb_eqs), on_pre = preCS, method = 'rk4', dt = t_step)
RE_RND_cs = Synapses(RNDg, REg, iRND, on_pre = preRND, method = 'rk4', dt = t_step)

#->PY Synapses
PY_RTC_cs = Synapses(RTCg,PYg,(PYFS_cs + iRTCa_eqs),on_pre = prePYFS, method = 'rk4', dt = t_step)
PY_PY_cs = Synapses(PYg,PYg,(PYFS_cs + iPY_eqs),on_pre = prePYFS, method = 'rk4', dt = t_step)
PY_FS_cs = Synapses(FSg,PYg,(PYFS_cs + iFS_eqs),on_pre = prePYFS, method = 'rk4', dt = t_step)
PY_RND = Synapses(RNDg,PYg,RND_PYFS,on_pre = preRND_PYFS, method = 'rk4', dt = t_step)

#->FS Synapses
FS_RTC_cs = Synapses(RTCg,FSg,(PYFS_cs + iRTCa_eqs),on_pre = prePYFS, method = 'rk4', dt = t_step)
FS_PY_cs = Synapses(PYg,FSg,(PYFS_cs + iPY_eqs),on_pre = prePYFS, method = 'rk4', dt = t_step)
FS_FS_cs = Synapses(FSg,FSg,(PYFS_cs + iFS_eqs),on_pre = prePYFS, method = 'rk4', dt = t_step) 
FS_RND = Synapses(RNDg,FSg,RND_PYFS,on_pre = preRND_PYFS, method = 'rk4', dt = t_step)
FS_RND = Synapses(RNDg,FSg,RND_PYFS,on_pre = preRND_PYFS, method = 'rk4', dt = t_step)


#Declare Connections

#->HTC
HTC_HTC_gj.connect('(x_pre - x_post)**2+(y_pre-y_post)**2 < 4.01*(i < j)' , p = '0.3')
HTC_HTC_gj.connect(i = HTC_HTC_gj.j[:], j = HTC_HTC_gj.i[:]) 
HTC_RTC_gj.connect('(x_pre*(nHTC-1.0)/(nRTC-1.0)-x_post)**2+(y_pre*(nHTC-1.0)/(nRTC-1.0)-y_post)**2 < 4.01*gjb_pre', p = '0.3')
HTC_RE_cs.connect(p = '0.2')
HTC_RND_cs.connect()

#->RTC
RTC_HTC_gj.connect(i = HTC_RTC_gj.j[:], j = HTC_RTC_gj.i[:])
RTC_IN_cs.connect(p = '0.3')
RTC_RE_cs.connect(p = '0.2')
RTC_PY_csa.connect(p = '0.23')
RTC_PY_csb.connect(i = RTC_PY_csa.i[:], j = RTC_PY_csa.j[:])
RTC_RND_cs.connect()

#->IN
IN_HTC_csa.connect(p = '0.3')
IN_HTC_csb.connect(i = IN_HTC_csa.i[:], j = IN_HTC_csa.j[:])
IN_RE_cs.connect(p = '0.05')
IN_RND_cs.connect()

#->RE
RE_RE_gj.connect('(x_pre - x_post)**2+(y_pre-y_post )**2 < 4.01*(gjb_pre*i < gjb_post*j)', p = '0.3')
RE_RE_gj.connect(i = RE_RE_gj.j[:], j = RE_RE_gj.i[:]) 
RE_HTC_csa.connect(p = '0.2')
RE_HTC_csb.connect(i = RE_HTC_csa.i[:], j = RE_HTC_csa.j[:])
RE_RTC_csa.connect(p = '0.2')
RE_RTC_csb.connect(i = RE_RTC_csa.i[:], j = RE_RTC_csa.j[:])
RE_RE_cs.connect(p = '0.2')
RE_PY_csa.connect(p = '0.3')
RE_PY_csb.connect(i = RE_PY_csa.i[:], j = RE_PY_csa.j[:])
RE_RND_cs.connect()

#->PY
PY_RTC_cs.connect(p = '0.04')
PY_PY_cs.connect(p = '0.5')
PY_FS_cs.connect('abs(i*(nPY-1.0)/(nFS-1.0) - j) < nPY*0.2', p = '0.8')
PY_RND.connect()

#->FS
FS_RTC_cs.connect(p = '0.02')
FS_PY_cs.connect(i = PY_FS_cs.j[:], j = PY_FS_cs.i[:])
FS_FS_cs.connect('abs(i-j) < nFS/2.0*(i!=j)', p = '0.8')
FS_RND.connect()


#Declare Synapse Variables

#->HTC
HTC_HTC_gj.rgap = 100*Mohm
HTC_RTC_gj.rgap = 300*Mohm
HTC_RE_cs.D_i = 1.07
HTC_RE_cs.g_syn = 3*nS
HTC_RE_cs.e_syn = -80*mV
HTC_RE_cs.alpha = 10.5/ms
HTC_RE_cs.beta = 0.166/ms
HTC_RND_cs.t_spike = -3*ms
HTC_RND_cs.g_spike = 1.5*nS
HTC_RND_cs.tau = 5*ms
#->RTC

RTC_HTC_gj.rgap = 300*Mohm
RTC_IN_cs.D_i = 1.07
RTC_IN_cs.g_syn = 3*nS
RTC_IN_cs.e_syn = -80*mV
RTC_IN_cs.alpha = 10.5/ms
RTC_IN_cs.beta = 0.166/ms
RTC_RE_cs.D_i = 1.07
RTC_RE_cs.g_syn = 3*nS
RTC_RE_cs.e_syn = -80*mV
RTC_RE_cs.alpha = 10.5/ms
RTC_RE_cs.beta = 0.166/ms
RTC_PY_csa.D_i = 1.07
RTC_PY_csa.g_syn = 4*nS
RTC_PY_csa.e_syn = 0*mV
RTC_PY_csa.alpha = 0.94/ms
RTC_PY_csa.beta = 0.18/ms
RTC_PY_csb.D_i = 1.07
RTC_PY_csb.g_syn = 2*nS
RTC_PY_csb.alpha = 1/ms
RTC_PY_csb.beta = 0.0067/ms
RTC_RND_cs.t_spike = -3*ms
RTC_RND_cs.g_spike = 1.5*nS
RTC_RND_cs.tau = 5*ms

#->IN
IN_HTC_csa.D_i = 1.07
IN_HTC_csa.g_syn = 6*nS
IN_HTC_csa.e_syn = 0*mV
IN_HTC_csa.alpha = 0.94/ms
IN_HTC_csa.beta = 0.18/ms
IN_HTC_csb.D_i = 1.07
IN_HTC_csb.g_syn = 3*nS
IN_HTC_csb.alpha = 1/ms
IN_HTC_csb.beta = 0.0067/ms
IN_RE_cs.D_i = 1.07
IN_RE_cs.g_syn = 1*nS
IN_RE_cs.e_syn = -80*mV
IN_RE_cs.alpha = 10.5/ms 
IN_RE_cs.beta = 0.166/ms
IN_RND_cs.t_spike = -3*ms
IN_RND_cs.g_spike = 1.5*nS
IN_RND_cs.tau = 5*ms

#->RE
RE_RE_gj.rgap = 300*Mohm
RE_HTC_csa.D_i = 1.07
RE_HTC_csa.g_syn = 4*nS
RE_HTC_csa.e_syn = 0*mV
RE_HTC_csa.alpha = 0.94/ms
RE_HTC_csa.beta = 0.18/ms
RE_HTC_csb.D_i = 1.07 
RE_HTC_csb.g_syn = 2*nS
RE_HTC_csb.alpha = 1/ms
RE_HTC_csb.beta = 0.0067/ms
RE_RTC_csa.D_i = 1.07
RE_RTC_csa.g_syn = 4*nS
RE_RTC_csa.e_syn = 0*mV
RE_RTC_csa.alpha = 0.94/ms
RE_RTC_csa.beta = 0.18/ms
RE_RTC_csb.D_i = 1.07
RE_RTC_csb.g_syn = 2*nS
RE_RTC_csb.alpha = 1/ms
RE_RTC_csb.beta = 0.0067/ms
RE_RE_cs.D_i = 1.07 
RE_RE_cs.g_syn = 1*nS
RE_RE_cs.e_syn = -70*mV
RE_RE_cs.alpha = 10.5/ms
RE_RE_cs.beta = 0.166/ms
RE_PY_csa.D_i = 1.07
RE_PY_csa.g_syn = 4*nS
RE_PY_csa.e_syn = 0*mV
RE_PY_csa.alpha = 0.94/ms
RE_PY_csa.beta = 0.18/ms
RE_PY_csb.D_i = 1.07
RE_PY_csb.g_syn = 2*nS
RE_PY_csb.alpha = 1/ms
RE_PY_csb.beta = 0.0067/ms
RE_RND_cs.t_spike = -3*ms
RE_RND_cs.g_spike = 1.5*nS
RE_RND_cs.tau = 5*ms

#->PY
PY_RTC_cs.tau = 2*ms
PY_RTC_cs.e_syn = 0*mV
PY_RTC_cs.g_spike = 0.3*nS
PY_PY_cs.tau = 2*ms
PY_PY_cs.e_syn = 0*mV
PY_PY_cs.g_spike = 0.3*nS
PY_FS_cs.tau = 10*ms
PY_FS_cs.e_syn = -70*mV
PY_FS_cs.g_spike = 0.3*nS
PY_RND.std = 0.1*pA
PY_RND.randi = randn()

#->FS
FS_RTC_cs.tau = 2*ms
FS_RTC_cs.e_syn = 0*mV
FS_RTC_cs.g_spike = 0.4*nS
FS_PY_cs.tau = 2*ms
FS_PY_cs.e_syn = 0*mV
FS_PY_cs.g_spike = 0.4*nS
FS_FS_cs.tau = 10*ms
FS_FS_cs.e_syn = -70*mV
FS_FS_cs.g_spike = 0.03*nS
FS_RND.std = 0.1*pA
FS_RND.randi = randn()

#State Monitors
TC_volt = StateMonitor(TCg,'v', record = True)
PY_lfp = StateMonitor(PYg, 'LFP', record = True)

#Spike Monitors
TC_spike = SpikeMonitor(TCg, record = True)
IN_spike = SpikeMonitor(INg, record = True)
RE_spike = SpikeMonitor(REg, record = True)
PYFS_spike = SpikeMonitor(PYFSg, record = True)

Sgn = TimedArray(0*sin(20*pi*t_d)*25*pA, dt = t_step)

run(duration)
