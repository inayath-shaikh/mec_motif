TITLE itheta
COMMENT
Theta oscillations gate the transmission of reliable sequences in the medial entorhinal cortex
Arun Neru, Collins Assisi
ENDCOMMENT


UNITS {
    (mV)=(millivolt)
    (S) = (siemens)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX i_theta
    NONSPECIFIC_CURRENT itheta 
    RANGE gnabar
    GLOBAL Amp,vthresh,omega
    

      
}

PARAMETER {
    Amp = 2e-4 (S/cm2)
    vthresh = -80 (mV)
    omega = 0.01 (1/ms)
    phi=0






}   

ASSIGNED {
        v (mV)
        itheta (mA/cm2)

} 




BREAKPOINT { 

	itheta = Amp*sin(2*3.14*omega*t+phi): *(v-vthresh) :phi not phase
    

}
