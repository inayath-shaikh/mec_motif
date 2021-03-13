TITLE Stellate cells mechanism
COMMENT
Acker, C. D., Kopell, N., & White, J. A. (2003). Synchronization of strongly coupled excitatory neurons: Relating network behavior to biophysics. Journal of Computational Neuroscience, 15(1), 71–90. https://doi.org/10.1023/A:1024474819512
ENDCOMMENT


UNITS {
    (mV)=(millivolt)
    (S) = (siemens)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX stellate_mech
    USEION na READ ena WRITE ina 
    USEION k READ ek WRITE ik 
    NONSPECIFIC_CURRENT il 
    NONSPECIFIC_CURRENT ih 
    RANGE gnabar,gkbar,gna,gk, gl,el,gnap,ghbar,gh,eh,ena
    GLOBAL amna,bmna,ahna,bhna,ank,bnk ,mhfinf,mhftau,mhsinf,mhstau,mnapinf,mnaptau
    

      
}

PARAMETER {
    gnap = 0.0005 (S/cm2)
    gnabar = 0.052 (S/cm2) 
    gkbar = 0.011 (S/cm2)

    ghbar = 0.0015 (S/cm2)
    ena = 55 (mV) :reset by neuron
    ek = -90 (mV) :reset by neuron
    el = -65 (mV)
    gl = 0.0005 (S/cm2)
    eh = -20 (mV)
    mnaptau = 0.15 (ms)







}   

ASSIGNED {
        v (mV)

	    gna (S/cm2)
	    gk (S/cm2)
        gh (S/cm2)

        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        ih (mA/cm2)

    amna (1/ms) bmna (1/ms) ahna (1/ms) bhna (1/ms) ank (1/ms) bnk (1/ms) 
    mhfinf mhsinf  mnapinf   
    mhstau (ms)  mhftau (ms)   
} 


STATE {
    mna hna mnap nk mhf mhs 
}


BREAKPOINT { 
    SOLVE states METHOD cnexp

    gna = (gnabar*mna*mna*mna*hna) + (gnap*mnap)
	ina = gna*(v - ena)
    gk = (gkbar*nk*nk*nk*nk)
	ik = gk*(v - ek)      
    il = gl*(v - el)
    gh = ghbar*(0.65*mhf+0.35*mhs)
    ih = gh*(v-eh)
    



}
INITIAL {
    rates(v)
    mna = amna/(amna+bmna)
    hna = ahna/(ahna+bhna)
    mnap = mnapinf

    nk = ank/(ank+bnk)

    mhf = mhfinf
    mhs = mhsinf
    

}

DERIVATIVE states {
    rates(v)

    mna' = amna*(1-mna) - bmna*mna
    hna' = ahna*(1-hna) - bhna*hna
    mnap' = (mnapinf-mnap)/mnaptau

    nk' = ank*(1-nk) - bnk*nk

    mhf' = (mhfinf-mhf)/(mhftau)
    mhs' = (mhsinf-mhs)/(mhstau)
    


}

PROCEDURE rates(v (mV)) {

    TABLE amna, bmna,ahna,bhna,ank,bnk,mhfinf,mhftau,mhsinf,mhstau,mnapinf FROM -100 TO 100 WITH 200 
    UNITSOFF
    amna = (.1)*vtrap(-(v+23),10)
    bmna = 4*exp(-(v+48)/18)


    ahna = 0.07*exp(-(v+37)/20)
    bhna = 1/(exp(-0.1*(v+7))+1)


    mnapinf = 1/(1+exp(-(v+38)/6.5))

 
    ank = 0.01*vtrap(-(v+27),10)
    bnk = 0.125 * exp(-(v+37)/80)


    mhfinf = 1/(1+exp((v+79.2)/9.78))
    mhftau = (0.51 / ((exp((v-1.7)/10)) + exp(-(v+340)/52))) + 1 

    mhsinf = 1/((1+exp((v+2.83)/15.9))^58)
    mhstau = (5.6/ ((exp((v-1.7)/14)) + exp(-(v+260)/43))) + 1
    


}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
