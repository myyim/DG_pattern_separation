TITLE SK channel (small conductance, calcium-activated potassium channel)
COMMENT

Original Mod File:
Original name 'gskch.mod', gsk granule
Santhakumar et al. (2005)
https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=51781&file=/dentategyrusnet2005/gskch.mod

Current version by A. Hanuschkin <AH, 2011> for:
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

Changes in current versus original version:
Correction: use of correct dynamics (see rate() lines: 95-101)

Further Mod File history:
- gsk granule
- modified from Aradi I, Holmes WR (1999) J Comput Neurosci 6:215-35

ENDCOMMENT

UNITS {
        (molar) = (1/liter)
        (mM)    = (millimolar)
	(mA)	= (milliamp)
	(mV)	= (millivolt)
}

NEURON {
	SUFFIX sk
	USEION k READ ek WRITE ik 
	USEION nca READ ncai VALENCE 2
	USEION lca READ lcai VALENCE 2
	USEION tca READ tcai VALENCE 2
	RANGE gsk, gskbar, qinf, qtau, ik
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	celsius=6.3 (degC)
	v		(mV)
	dt		(ms)
	gskbar  (mho/cm2)
	ek	(mV)
	cai (mM)
	ncai (mM)
	lcai (mM)
	tcai (mM)
}

STATE { q }

ASSIGNED {
	ik (mA/cm2) gsk (mho/cm2) qinf qtau (ms) qexp
}


BREAKPOINT {          :Computes i=g*q^2*(v-ek)
	SOLVE state
        gsk = gskbar * q*q
	ik = gsk * (v-ek)
}

UNITSOFF

INITIAL {
	cai = ncai + lcai + tcai	
	rate(cai)
	q=qinf
	VERBATIM
	ncai = _ion_ncai;
	lcai = _ion_lcai;
	tcai = _ion_tcai;
	ENDVERBATIM
}


PROCEDURE state() {  :Computes state variable q at current v and dt.
	cai = ncai + lcai + tcai
	rate(cai)
	q = q + (qinf-q) * qexp
	VERBATIM
	return 0;
	ENDVERBATIM
}

LOCAL q10
PROCEDURE rate(cai) {  :Computes rate and other constants at current v.
	LOCAL alpha, beta, tinc
	q10 = 3^((celsius - 6.3)/10) : q10=1 for 6.3 celcius
		:"q" activation system

        : this is the correct dynamics <AH>
	alpha = 0.00246/exp((12*log10(cai)+28.48)/-4.5)
	beta = 0.006/exp((12*log10(cai)+60.4)/35)
	qtau = 1 / (alpha + beta)
	qinf = alpha * qtau
	tinc = -dt*q10
	qexp = 1 - exp(tinc/qtau)*q10
}

UNITSON
