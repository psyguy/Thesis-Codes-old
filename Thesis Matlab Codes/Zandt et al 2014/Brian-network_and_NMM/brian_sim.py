import time
print('Preparing network')
tic = time.time()
import makepops2 as mp

from brian import *

## parameters
simdur = 500*ms
Ncellex = 100
Ncellin = 99
P=0.5           # connection probability (all synaptic populations)
Lambda_rate = 600/second
gamma = 1.0/(array([45,45,34,45,34]))
g0    = array([0.0054, 0.5*0.0057/Ncellex/P, 0.5*1.32*6*0.0057/Ncellex/P, 15*5*0.5/1.32*0.015/Ncellin/P, 0.5*0.015/Ncellin/P])

std_cellparam = [0.005, 0.005];    # standard deviation of the additional synaptic conductance for both populations
offsetg = 0.0;            	       # add additional constant synaptic input to the neurons
gsynex = randn(1,Ncellex)*std_cellparam[0]
gsynex = gsynex - mean(gsynex) + offsetg;
gsynin = randn(1,Ncellin)*std_cellparam[1];
gsynin = gsynin - mean(gsynin) + offsetg;

## generate populations and connections
# generate poisson input population
parsPoisson = {'rate1': Lambda_rate, 'rate2': Lambda_rate/1.5, 'StepTime': 0.5*second}
PoissonCells = mp.makepop('PoissonStep',Ncellex,parsPoisson)
# generate Excitatory population
# use dictionary to pass parameter values
parsHHe = {'Iinput': 0*uA/cm**2,'Cm': 10.0 *uF/cm**2, 'g_na': 100.0*msiemens/cm**2,'g_naL': 0.0175*msiemens/cm**2,
'g_k': 40.0*msiemens/cm**2,'g_kL':0.05*msiemens/cm**2,'g_clL':0.05*msiemens/cm**2,'phi': 3.0/msecond,
'E_k':-95*mV,'E_na': 53*mV,'E_cl': -82*mV,'gsyntest': gsynex*msiemens/cm**2,'Esyntest': 50* mV,'Vm': -65*mV,'n':0.07,'h':0.97,
'Esyn_ext':50 * mV, 'Esyn_e': 50 * mV, 'Esyn_i': -82 * mV}
HHe = mp.makepop('HH',Ncellex,parsHHe)
@network_operation
def HHe_spike():    # define that Vm needs to go below a certain value to generate the next spike, to prevent generation of spikes from a cell in depolarization block
    HHe.spiked = 0
    HHe.spiked[(HHe.was_above_threshold==0) & (HHe.Vm > 0*mV)] = 1
    HHe.was_above_threshold[HHe.Vm < -10*mV] = 0
    HHe.was_above_threshold[HHe.Vm > 0*mV] = 1
#generate Inhibitory population
parsHHi = parsHHe
parsHHi['gsyntest'] = gsynin*msiemens/cm**2
HHi = mp.makepop('HH',Ncellin,parsHHe)
@network_operation
def HHi_spike():    # define that Vm needs to go below a certain value to generate the next spike, to prevent generation of spikes from a cell in depolarization block
    HHi.spiked = 0
    HHi.spiked[(HHi.was_above_threshold==0) & (HHi.Vm > 0*mV)] = 1
    HHi.was_above_threshold[HHi.Vm < -10*mV] = 0
    HHi.was_above_threshold[HHi.Vm > 0*mV] = 1

synmodel = '''  dg/dt=g1 : msiemens/cm**2
                dg1/dt = -2*gamma*g1 - gamma*gamma*g : msiemens/cm**2/second
                gamma : 1/second
                g0 : msiemens/cm**2'''

# synmodel = '''  dg/dt=g1 : msiemens/cm**2
#                 dg1/dt = -2*gamma*g1 - gamma*gamma*g : msiemens/cm**2/second
#                 Isyn = g*(Esyn-Vm_post) : uA/cm**2
#                 Esyn : mV
#                 gamma : 1/second
#                 g0 : msiemens/cm**2'''
                
#dg/dt=-g/(80*msecond) : msiemens/cm**2
#           Isyn = g*(Esyn-Vm_post) : uA/cm**2
#           Esyn : mV
#           w : msiemens/cm**2'''

# ext -> e synapses
Synexte = Synapses(PoissonCells,HHe,model=synmodel,pre='g1+=e*gamma*g0')
HHe.gsyn_ext = Synexte.g
Synexte[:,:]='i==j'
#Synexte.Esyn = 50 * mV
Synexte.g0[:] = g0[0]*msiemens/cm**2
Synexte.gamma[:] = gamma[0]/msecond

# e -> e synapses
Synee = Synapses(HHe,HHe,model=synmodel,pre='g1+=e*gamma*g0')
HHe.gsyn_e = Synee.g
ConMat2 = rand(Ncellex,Ncellex)<P
Synee[:,:]='ConMat2[j,i]';# Synee.Esyn = 50 * mV;
Synee.g0[:] = g0[1]*msiemens/cm**2
Synee.gamma[:] = gamma[1]/msecond

# e -> i synapses
Synei = Synapses(HHe,HHi,model=synmodel,pre='g1+=e*gamma*g0')
HHi.gsyn_e = Synei.g
w = 0.001 * msiemens/cm**2;
ConMat3 = rand(Ncellin,Ncellex)<P
Synei[:,:]='ConMat3[j,i]';# Synei.Esyn = 50 * mV;
Synei.g0[:] = g0[2]*msiemens/cm**2
Synei.gamma[:] = gamma[2]/msecond

# i -> e synapses
Synie = Synapses(HHi,HHe,model=synmodel,pre='g1+=e*gamma*g0')
HHe.gsyn_i = Synie.g
w = 0.001 * msiemens/cm**2;
ConMat4 = rand(Ncellex,Ncellin)<P
Synie[:,:]='ConMat4[j,i]';# Synie.Esyn = -82 * mV;
Synie.g0[:] = g0[3]*msiemens/cm**2
Synie.gamma[:] = gamma[3]/msecond

# i -> i synapses
Synii = Synapses(HHi,HHi,model=synmodel,pre='g1+=e*gamma*g0')
HHi.gsyn_i = Synii.g
w = 0.001 * msiemens/cm**2;
ConMat5 = rand(Ncellin,Ncellin)<P
Synii[:,:]='ConMat5[j,i]';# Synii.Esyn = -82 * mV;
Synii.g0[:] = g0[4]*msiemens/cm**2
Synii.gamma[:] = gamma[4]/msecond

# generate Brian network and monitors
sm_ext = SpikeMonitor(PoissonCells,record=True)
sm_e = SpikeMonitor(HHe,record=True)
sm_i = SpikeMonitor(HHi,record=True)
synstate_ext = StateMonitor(HHe,'gsyn_ext',record=True)
synstate_ee = StateMonitor(HHe,'gsyn_e',record=True)
synstate_ie = StateMonitor(HHe,'gsyn_i',record=True)
synstate_ei = StateMonitor(HHi,'gsyn_e',record=True)
synstate_ii = StateMonitor(HHi,'gsyn_i',record=True)

net = MagicNetwork(verbose=True)
net.add(HHe,HHi,PoissonCells,Synexte,Synee,Synei, Synie, Synii, sm_ext, sm_e, sm_i, synstate_ext,synstate_ee,synstate_ei,synstate_ie,synstate_ii)

## run simulation
toc = time.time(); print('preparation took ' + str(toc-tic) + ' seconds')
print('Simulating...')
tic = time.time()
net.run(simdur)
toc = time.time(); print('Simulation took ' + str(toc-tic) + ' seconds')
#raster_plot(sm_e)
#figure()
#plot(synstate_ie.times/msecond, transpose(synstate_ie[0:10]) /(msiemens/cm**2))
#show()

# save spike times of neurons, and conductions of synapses (or Isyn of neurons?)

## save the network data
import scipy.io as sio
Spoisson = concatenate(sm_ext.spiketimes.values())/msecond # save spike timings of external input in ms
SynState = transpose(concatenate((synstate_ext[:],synstate_ee[:],synstate_ei[:],synstate_ie[:],synstate_ii[:]))/(msiemens/cm**2))
Spikes_e = array(sm_e.spikes)[:,::-1] # switch columns
Spikes_i = array(sm_i.spikes)[:,::-1] # switch columns
Spikes_i[:,1] = Spikes_i[:,1] + Ncellex # create unique cell numbers for i cells (right above numbers of e cells)
Spikes = concatenate((Spikes_e,Spikes_i))
sio.savemat('Network_brian.mat',{'Ncellex':double(Ncellex),'Ncellin':double(Ncellin),'ConMat2':double(ConMat2), 'ConMat3':double(ConMat3),'ConMat4':double(ConMat4),'ConMat5':double(ConMat5), 'gsynex': gsynex, 'gsynin': gsynin, 'offsetg':offsetg,
'g0':g0, 'gamma':gamma*1e3, 't_brian':synstate_ee.times/msecond,'Spoisson':Spoisson, 'SynState':SynState, 'Spikes': Spikes})