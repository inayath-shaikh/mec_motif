TITLE Stellate synapse
COMMENT
Acker, C. D., Kopell, N., & White, J. A. (2003). Synchronization of strongly coupled excitatory neurons: Relating network behavior to biophysics. Journal of Computational Neuroscience, 15(1), 71–90. https://doi.org/10.1023/A:1024474819512
ENDCOMMENT


UNITS {
    (mV)=(millivolt)
    (uS) = (microsiemens)
    (nA) = (nanoamp)
}

NEURON {
    POINT_PROCESS stell_syn
    NONSPECIFIC_CURRENT isyn 
    RANGE  gsynbar,bsyn,asyn,msyn
    POINTER vpre

      
}

PARAMETER {
   

    gsynbar = 0.14e-4 (uS)
    esyn = 0 (mV)
    asyn = 100e-1 (1/s)
    bsyn = 0.33e-1 (1/s)

}   

ASSIGNED {
    vpre (mV)
    v (mV)
    gsyn (uS)
    isyn (nA)

} 


STATE {
    msyn
}


BREAKPOINT { 
    SOLVE states METHOD cnexp

    gsyn = gsynbar*msyn
    isyn = gsyn*(v-esyn)
    
}
UNITSOFF
INITIAL {

    msyn =  asyn*F(vpre)/(asyn*F(vpre)+bsyn)
    

}

DERIVATIVE states {

    msyn' = asyn*F(vpre)*(1-msyn) - bsyn*msyn

    
}



FUNCTION F(vpre) {
    F = (1+tanh(vpre/4))/2
}
UNITSON



NET_RECEIVE (weight (uS)){
  msyn=msyn+weight
}