TITLE A-type (fast inactivating) Kv channel
COMMENT

Original Mod Files:
Original names 'bgka.mod' or 'borgka.mod'
Santhakumar V, Aradi I, Soltesz I (2005) J Neurophysiol 93:437-53 
https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=51781&file=/dentategyrusnet2005/bgka.mod
Migliore M, Cook EP, Jaffe DB, Turner DA, Johnston D (1995) J Neurophysiol 73:1157-68 
https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=3263&file=/ca3_db/borgka.mod

Current version by A. Hanuschkin <AH, 2011> for:
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

ENDCOMMENT

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
        ek (mV)
	celsius =6.3	(degC)  : default Neuron value (implictly used)
	gkabar	=.01 	(mho/cm2)
        vhalfn	=-33.6  (mV)
        vhalfl	=-83   	(mV)
        a0n	=0.02   (/ms)
        a0l     =0.08   (/ms)
        zetan	=-3    	(1)
        zetal	=4    	(1)
        gmn	=0.6   	(1)
        gml	=1   	(1)
}


NEURON {
	SUFFIX ka
	USEION k READ ek WRITE ik
        RANGE gkabar,gka, ik
        GLOBAL ninf,linf,taul,taun
}

STATE {
	n
        l
}

INITIAL {
        rates(v)
        n=ninf
        l=linf
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        linf      
        taul
        taun
        gka
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gka = gkabar*n*l
	ik = gka*(v-ek)
}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpl(v(mV)) {
  alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betl(v(mV)) {
  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states { 
        rates(v)
        n' = (ninf - n)/taun
        l' = (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,q10
        q10=3^((celsius-30)/10)
        a = alpn(v)
        ninf = 1/(1 + a)
        taun = betn(v)/(q10*a0n*(1+a))
        a = alpl(v)
        linf = 1/(1+ a)
        taul = betl(v)/(q10*a0l*(1 + a))
}

