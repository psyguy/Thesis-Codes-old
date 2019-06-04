from brian import *
def makepop(celltypeid,Ncell,pardict):

    if celltypeid == 'HH': # HH neuron
        pop = NeuronGroup(Ncell, HHeqs(), freeze=True, compile=False, threshold='spiked==1')        

        for key in pardict.keys():
            setattr(pop,key,pardict[key])

    elif celltypeid == 'PoissonStep':
        rate1 = pardict['rate1']; rate2 = pardict['rate2']; StepTime = pardict['StepTime']
        pop = PoissonGroup(Ncell,rates=lambda t:rate1+(rate2-rate1)*(t>StepTime))
        #for key in pardict.keys():
        #    setattr(pop,key,pardict[key])
    else:
        pop = 'Population undefined'
        Warning('Population undefined')
    return(pop)
    
    
    
def HHeqs():
    eqs = Equations('''
    dVm/dt = (1/Cm)*(-Ina -Ik -Icl +Iinput + Iconstsyn + Isyn_ext + Isyn_e + Isyn_i)        : mV
    dn/dt = phi*(alpha_n*(1-n)-beta_n*n)                                    : 1
    dh/dt = phi*(alpha_h*(1-h)-beta_h*h)                                    : 1
    
    alpha_n =0.01 * (Vm/mV+34.0)/( 1.0 - exp(-0.1 * (Vm/mV+34.0)) ) : 1
    beta_n = 0.125 * exp(-(Vm/mV+44.0)/80.0)                     : 1
    alpha_m = 0.1 * (Vm/mV+30.0)/( 1.0 - exp(-0.1 * (Vm/mV+30.0)) ) : 1 
    beta_m = 4.0 * exp(-(Vm/mV+55.0)/18.0)                       : 1
    alpha_h = 0.07 * exp(-(Vm/mV+44.0)/20.0)                     : 1
    beta_h = 1.0/( 1.0 + exp(-0.1 * (Vm/mV+14.0)) )              : 1
    m_inf = alpha_m/(alpha_m + beta_m)                        : 1
    
    Ina = g_na*(m_inf**3)*h*(Vm-E_na) + g_naL*(Vm-E_na)       : uA/cm**2
    Ik =  (g_k*n**4)*(Vm-E_k) + g_kL*(Vm-E_k)                 : uA/cm**2
    Icl = g_clL*(Vm-E_cl)                                     : uA/cm**2
    Iconstsyn = gsyntest*(Esyntest-Vm)                        : uA/cm**2
    Isyn_ext = gsyn_ext*(Esyn_ext-Vm) : uA/cm**2
    Isyn_e   = gsyn_e*(Esyn_e  -Vm)   : uA/cm**2
    Isyn_i   = gsyn_i*(Esyn_i  -Vm)   : uA/cm**2
    ''')
    eqs += 'spiked : 1'
    eqs += 'was_above_threshold : 1'
    eqs += 'Iinput : uA/cm**2'
    eqs += 'Cm : uF/cm**2'
    eqs += 'g_na : msiemens/cm**2'
    eqs += 'g_naL : msiemens/cm**2'
    eqs += 'g_k : msiemens/cm**2'
    eqs += 'g_kL : msiemens/cm**2'
    eqs += 'g_clL : msiemens/cm**2'
    eqs += 'phi : 1/msecond'
    eqs += 'E_k : mV'
    eqs += 'E_na : mV'
    eqs += 'E_cl : mV'
    eqs += 'gsyntest : msiemens/cm**2'
    eqs += 'Esyntest : mV'
    eqs += 'gsyn_ext : msiemens/cm**2'
    eqs += 'gsyn_e   : msiemens/cm**2'
    eqs += 'gsyn_i   : msiemens/cm**2'
    eqs += 'Esyn_ext : mV'
    eqs += 'Esyn_e   : mV'
    eqs += 'Esyn_i   : mV'
    
    return(eqs)