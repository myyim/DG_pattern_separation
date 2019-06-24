TITLE ichan2.mod combination of Nav and Kv channels for Hodgin-Huxely-type action potential mechanism
 
COMMENT

Original Mod File:
Original name 'ichan2.mod'
Santhakumar V, Aradi I, Soltesz I (2005) J Neurophysiol 93:437-53 
https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=51781&file=%2fdentategyrusnet2005%2fichan2.mod
Morgan RJ, Soltesz I (2008) Proc Natl Acad Sci U S A 105:6179-84
Morgan RJ, Santhakumar V, Soltesz I (2007) Prog Brain Res 163:639-58
Cutsuridis V, Cobb S, Graham BP (2009) Hippocampus 20(3):423-46 

Current version by A. Hanuschkin <AH, 2011> for:
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

Changes in current versus original version:
 - added a tonic (leak) GABAA conductance to be modified during epilepsy (see Young CC, Stegen M, Bernard R, Muller M, Bischofberger J, Veh RW, Haas CA, Wolfart J (2009) J Physiol 587:4213-4233)
 - checked, simplified (reduced to single k Ion, Ekf=Eks=Ek) & commented by A. Hanuschkin 2011,2012

Mod File history:
I_Na: (equivalent to Santhakumar et al 2005)
* modified parameters (shifted voltage dependence by 68mV) 
  Aradi I and Soltesz I (2002) J Physiol. 538(Pt 1):227-51.
* dynamics (NOTE in this paper (at least) sign error in alpha_m, beta_m and alpha_n! (1-exp()) -> should be (exp()-1) <ah>)
Aradi and Soltesz (2002)
* modified from 
  Yuen GL, Durand D, (1991) Neuroscience 41(2-3):411-23.
* Aradi and Soltesz (2002) parameters =  
  Aradi I, Holmes WR (1999) J Comput Neurosci 6:215-35
parameters = 
  Yuen and Durand (1991) parameters + V shift of 16mV

I_Kf: (equivalent to Santhakumar et al 2005)
* modified parameters Aradi and Soltesz (2002)/Aradi & Holmes (1999) (shifted voltage dependence by 65mV) 
* NOTE typo in formular beta_n in Aradi and Soltesz (2002): beta_n = 0.264/exp((v-22)/4) -> should be beta_n = 0.264/exp((v-22)/40) <ah>
* modified parameters from Yuen and Durand (1991)

I_Ks: (equivalent to Santhakumar et al 2005)
* modified parameters Aradi & Holmes (1999) (shifted voltage dependence by 65mV) 

I_leak: (equivalent to Santhakumar et al 2005)

I_GABAA: (tonic GABAA leak (see above), added in Yim et al (2015))
* replicated from I_leak
 
A. Hanuschkin(c) 2011,2012

ENDCOMMENT
 
UNITS {
        (mA) =(milliamp)
        (mV) =(millivolt)
        (uF) = (microfarad)
	(molar) = (1/liter)
	(nA) = (nanoamp)
	(mM) = (millimolar)
	(um) = (micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
NEURON { 
: "suffix marks the mechanism to be distributed and whose variables & parameters are identified in hoc by a particular suffix" The Neuron Book Chap 9.5
SUFFIX ichan2

: ION usage block
USEION na READ ena WRITE ina VALENCE 1			: Na current
USEION k READ ek WRITE ik  				: K current
NONSPECIFIC_CURRENT il, igabaa 				: leak current

: range variable definition block,
: i.e. variables that might change with space along a compartment / could be declared global in this case
RANGE gnatbar, gkfbar, gksbar				: gbar values for Na, K(slow/fast)
RANGE gl, el, ina, ik, il, ggabaa, igabaa, egabaa	: gbar and reversal poti for leak current	
}

: The INDEPENDENT statement was omitted; INDEPENDENT statements are irrelevant to NEURON. http://www.neuron.yale.edu/phpbb/viewtopic.php?f=16&t=2351 
: INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
: Variables whose values are normally specified by the user are parameters, and are declared in a PARAMETER block.
: Variables in the parameter section will have global scope
PARAMETER {						
        :v (mV) 
        :celsius = 6.3 (degC)
        :dt (ms) 

        gnatbar (mho/cm2)   				: Na (gbar and reversal poti)
        ena  	(mV)	
		
	gkfbar 	(mho/cm2)				: K  (gbar(slow/fast), reversal is ek)
	gksbar = 0 (mho/cm2)	                        : init to 0 (not included in BC, HIPP and MC) <ah>
        ek     	(mV)                      

	gl 	(mho/cm2)    				: leak (gbar and reversal poti)
 	el 	(mV)

	ggabaa 	(mho/cm2)    				: GABAA (gbar and reversal poti)
 	egabaa 	(mV)
}


: "If a model involves differential equations [..] their dependent variables or unknowns are to be listed in the STATE block" The Neuron Book Chap 9.5
STATE {
	m h nf ns
}

: The ASSIGNED block is used to declare two kinds of variables
: 1) those given values outside the mod file (variables potentially available to every mechanism (e.g. v, celsius,t..)
: 2) left hand side of assignment statements (unknowns in set of equations, dependent variables in differential euqtions ...)
ASSIGNED {		
: 1)
        v (mV) 
        celsius (degC)
        dt (ms) 
	
: 2) 
        gna (mho/cm2) 					: Na
        ina (mA/cm2)
	
        gkf (mho/cm2)					: K
        gks (mho/cm2)
	ik (mA/cm2)

	il (mA/cm2)					: leak 

	igabaa (mA/cm2)					: GABAA 

	minf hinf nfinf nsinf				: left hand side of differential equations
 	mtau (ms) htau (ms) nftau (ms) nstau (ms)	: and other assignment variables
	mexp hexp nfexp nsexp
	q10
} 

: This block is evaluated every time step. 
BREAKPOINT {
	SOLVE states					: here the state variables are updated 
        gna = gnatbar*m*m*m*h  			: calculated g at timepoint t
        gkf = gkfbar*nf*nf*nf*nf
        gks = gksbar*ns*ns*ns*ns

        ina = gna*(v - ena)				: calculated currents flowing
       	ik = gkf*(v-ek) + gks*(v-ek)
	il = gl*(v-el)
	igabaa = ggabaa*(v-egabaa)
}
 
: UNITSOFF
 
: Called from Neuron during initializing the model
INITIAL {
	trates(v)
	
	m = minf
	h = hinf
        nf = nfinf
        ns = nsinf
	
	VERBATIM
	return 0;
	ENDVERBATIM
}

: discreticed versions of the differential equations, hence a PROCEDURE and not DERIVATIVE block
PROCEDURE states() {	: Computes state variables m, h, and n 
        trates(v)	: at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        nf = nf + nfexp*(nfinf-nf)
        ns = ns + nsexp*(nsinf-ns)
        VERBATIM
        return 0;
        ENDVERBATIM
}

: moved this to assign block <ah> 
: LOCAL q10

:Computes rate and other constants at current v.
PROCEDURE rates(v) {  
        LOCAL  alpha, beta, sum
        q10 = 3^((celsius - 6.3)/10)
                :"m" sodium activation system - act and inact cross at -40	: shifted by 68mV compared to in Aradi 1999/2002
	alpha = -0.3*vtrap((v+60-17),-5)		: in Aradi 1999: alpha = -0.3*vtrap((v-25),-5); in Aradi 2002: alpha = 0.3*vtrap((v-25),-5) <ah>
	beta = 0.3*vtrap((v+60-45),5)			: in Aradi 1999: beta = 0.3*vtrap((v-53),5);  in Aradi 2002:  beta = -0.3*vtrap((v-53),5) <ah>
	sum = alpha+beta        
	mtau = 1/sum      minf = alpha/sum
                :"h" sodium inactivation system		: shifted by 68mV compared to in Aradi 1999/2002
	alpha = 0.23/exp((v+60+5)/20)			: in Aradi 1999/2002:  alpha = 0.23/exp((v-3)/20) <ah>
	beta = 3.33/(1+exp((v+60-47.5)/-10))		: in Aradi 1999/2002:  beta = 3.33/(1+exp((v-55.5)/-10)) <ah>
	sum = alpha+beta
	htau = 1/sum 
        hinf = alpha/sum 


             :"ns" sKDR activation system		: shifted by 65mV compared to Aradi 1999 <ah>
        alpha = -0.028*vtrap((v+65-35),-6)		: in Aradi 1999: alpha = -0.028*vtrap((v-35),-6) 
	beta = 0.1056/exp((v+65-10)/40)			: in Aradi 1999: beta = 0.1056/exp((v-10)/40)   
	sum = alpha+beta        			
	nstau = 1/sum      nsinf = alpha/sum		
            :"nf" fKDR activation system		: shifted by 65mV compared to Aradi 1999/2002 <ah>
        alpha = -0.07*vtrap((v+65-47),-6)		: in Aradi 1999: alpha = -0.07*vtrap((v-47),-6); in Aradi 2002: alpha = 0.07*vtrap((v-47),-6) <ah>
	beta = 0.264/exp((v+65-22)/40)			: in Aradi 1999/2002: beta = 0.264/exp((v-22)/40)  // probably typo in Aradi & Soltez 2002 there: beta = 0.264/exp((v-22)/4) <ah>
	sum = alpha+beta        
	nftau = 1/sum      nfinf = alpha/sum
}

: Computes rate and other constants at current v. 
PROCEDURE trates(v) {  
	LOCAL tinc
        : TABLE minf, mexp, hinf, hexp, nfinf, nfexp, nsinf, nsexp, mtau, htau, nftau, nstau   : <ah>
	: DEPEND dt, celsius FROM -100 TO 100 WITH 200					       : <ah>
                           
	rates(v)	: not consistently executed from here if usetable_hh == 1
			: so don't expect the tau values to be tracking along with
			: the inf values in hoc

        tinc = -dt * q10
        mexp = 1 - exp(tinc/mtau)
        hexp = 1 - exp(tinc/htau)
	nfexp = 1 - exp(tinc/nftau)
	nsexp = 1 - exp(tinc/nstau)
}

:Traps for 0 in denominator of rate eqns.
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
:UNITSON

