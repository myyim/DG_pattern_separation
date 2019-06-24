TITLE BK channel (big conductance, calcium-activated potassium channel)
COMMENT

Original Mod File:
Original name 'CaBK.mod'
Santhakumar V, Aradi I, Soltesz I (2005) J Neurophysiol 93:437-53 
https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=51781&file=/dentategyrusnet2005/CaBK.mod
Original name 'cagk.mod'
Migliore M, Cook EP, Jaffe DB, Turner DA, Johnston D (1995) J Neurophysiol 73:1157-68 
https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=3263&file=/ca3_db/cagk.mod

Current version by A. Hanuschkin <AH, 2011> for:
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

Further Mod File history:
Moczydlowski and Latorre (1983) J. Gen. Physiol. 82

ENDCOMMENT


UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


NEURON {
	SUFFIX bk
	USEION nca READ ncai VALENCE 2
	USEION lca READ lcai VALENCE 2
	USEION tca READ tcai VALENCE 2
	USEION k READ ek WRITE ik
	RANGE gkbar,gkca, ik
	GLOBAL oinf, otau
}

UNITS {
	FARADAY = (faraday)  (kilocoulombs)
	R = 8.313424 (joule/degC)
}

PARAMETER {
	celsius		(degC)
	v		(mV)
	gkbar=.01	(mho/cm2)	: Maximum Permeability
	cai = 5.e-5	(mM)
	ek		(mV)

	d1 = .84
	d2 = 1.
	k1 = .48e-3	(mM)
	k2 = .13e-6	(mM)
	abar = .28	(/ms)
	bbar = .48	(/ms)
        st=1            (1)
	lcai		(mV)
	ncai		(mV)
	tcai		(mV)
}

ASSIGNED {
	ik		(mA/cm2)
	oinf
	otau		(ms)
        gkca          	(mho/cm2)
}

INITIAL {
	cai= ncai + lcai + tcai
        rate(v,cai)
        o=oinf
}

STATE {	o }		: fraction of open channels

BREAKPOINT {
	SOLVE state METHOD cnexp
	gkca = gkbar*o^st
	ik = gkca*(v - ek)
}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
: Please note that cai was not assiged here in the original Santhakumar (2005) version (which we used). It should be cai = ncai + lcai + tcai, as noted by
: Morgan RJ, Santhakumar V, Soltesz I (2007) Prog Brain Res 163:639-58
: See
: https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=124513&file=/dentate_gyrus/CaBK.mod
	rate(v, cai)
	o' = (oinf - o)/otau
}

FUNCTION alp(v (mV), c (mM)) (1/ms) { :callable from hoc
	alp = c*abar/(c + exp1(k1,d1,v))
}

FUNCTION bet(v (mV), c (mM)) (1/ms) { :callable from hoc
	bet = bbar/(1 + c/exp1(k2,d2,v))
}

FUNCTION exp1(k (mM), d, v (mV)) (mM) { :callable from hoc
	exp1 = k*exp(-2*d*FARADAY*v/R/(273.15 + celsius))
}

PROCEDURE rate(v (mV), c (mM)) { :callable from hoc
	LOCAL a
	a = alp(v,c)
	otau = 1/(a + bet(v, c))
	oinf = a*otau
}

