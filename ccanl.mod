COMMENT

Original Mod Files:
Santhakumar V, Aradi I, Soltesz I (2005) J Neurophysiol 93:437-53 
https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=51781&file=/dentategyrusnet2005/ccanl.mod
Morgan RJ, Soltesz I (2008) Proc Natl Acad Sci U S A 105:6179-84
Morgan RJ, Santhakumar V, Soltesz I (2007) Prog Brain Res 163:639-58
Dyhrfjeld-Johnsen J, Santhakumar V, Morgan RJ, Huerta R, Tsimring L, Soltesz I (2007) J Neurophysiol 97:1566-87 https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=124513&file=/dentate_gyrus/ccanl.mod

Current version by A. Hanuschkin <AH, 2011> for:
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

Mod File history:
calcium accumulation into a volume of area*depth next to the membrane with a decay (time constant tau) to resting level
given by caiinf	updating ECa for CaN.mod, CaT.mod and CaL.mod <ah 2011>

Warning by Ted Carnevale 2015:
The expression that this mechanism uses to calculate the contribution of ica to the rate of change of calcium concentration in the shell is 
-ica*(1e7)/(depth*FARADAY)
but it should really be
-ica*(1e7)/(depth*2*FARADAY)
because the valence of ca is 2.  The result of this omission is that the mechanism behaves as if the shell is only 1/2 as thick as the value specified by the depth parameter.

ENDCOMMENT

NEURON {
	SUFFIX ccanl
USEION nca READ ncai, inca, enca WRITE enca, ncai VALENCE 2
USEION lca READ lcai, ilca, elca WRITE elca, lcai VALENCE 2
USEION tca READ tcai, itca, etca WRITE etca, tcai VALENCE 2
RANGE caiinf, catau, cai, ncai, lcai,tcai, eca, elca, enca, etca
}

UNITS {
        (mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (milli/liter)
	(mA) = (milliamp)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}

INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

PARAMETER {
        celsius = 6.3 (degC)
	depth = 200 (nm)	: assume volume = area*depth
	catau = 9 (ms)
	caiinf = 50.e-6 (mM)	: takes precedence over cai0_ca_ion
			: Do not forget to initialize in hoc if different
			: from this default.
	cao = 2 (mM)
	ica (mA/cm2)
	inca (mA/cm2)
	ilca (mA/cm2)
	itca (mA/cm2)
	cai= 50.e-6 (mM)
}

ASSIGNED {
	enca (mV)
	elca (mV)
	etca (mV)
	eca (mV)
}

STATE {
	ncai (mM)
	lcai (mM)
	tcai (mM)
}

INITIAL {
	VERBATIM
	ncai = _ion_ncai;
	lcai = _ion_lcai;
	tcai = _ion_tcai; 
	ENDVERBATIM
	ncai=caiinf/3
	lcai=caiinf/3
	tcai=caiinf/3
	cai = caiinf	
	eca = ktf() * log(cao/caiinf)	
	enca = eca
	elca = eca
	etca = eca
}


BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
	cai = ncai+lcai+tcai 	
	eca = ktf() * log(cao/cai)	
	enca = eca : communicate new Ca reversal potential
	elca = eca : communicate new Ca reversal potential
	etca = eca : communicate new Ca reversal potential
}

DERIVATIVE integrate {
ncai' = -(inca)/depth/FARADAY * (1e7) + (caiinf/3 - ncai)/catau
lcai' = -(ilca)/depth/FARADAY * (1e7) + (caiinf/3 - lcai)/catau
tcai' = -(itca)/depth/FARADAY * (1e7) + (caiinf/3 - tcai)/catau
}

FUNCTION ktf() (mV) {
	ktf = (1000)*R*(celsius +273.15)/(2*FARADAY)
} 
