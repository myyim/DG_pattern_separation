### Analysis of DG network data ###
# This Python code creates a scatter plot of output vs input sim scores.
# Enter the idname

# ModelDB file along with publication:
# Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
# http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

# modified and augmented by
# Man Yi Yim / 2015
# Alexander Hanuschkin / 2011

import pylab
import numpy
import matplotlib as mpl
font_size = 12	# 20
mpl.rcParams['axes.titlesize'] = font_size+10
mpl.rcParams['xtick.labelsize'] = font_size+6
mpl.rcParams['ytick.labelsize'] = font_size+6
mpl.rcParams['axes.labelsize'] = font_size+8
mpl.rcParams['legend.fontsize'] = font_size
mpl.rcParams['font.size'] = font_size+10

sid = 0
sidm = 12
nid = 6
step = 2
#delta = numpy.array([5,10,20,30])/2.	# half the width of the kernel base
delta = numpy.array([20])/2.
dur = 200.	# duration of the simulation
bin = 1.
ninput = 100	# number of input sources
ncell = 500	# number of neurons (GCs)
#idname = "-pp10-gaba1-kir1-st0"

pylab.figure(figsize=(15,16))
for k in range(4):	#number to show
	pylab.subplot(4,1,k+1)
	stimin = numpy.loadtxt('StimIn-'+str(step*k)+'-'+str(nid)+idname+'.txt')
	for j in range(stimin.size/2):
		if stimin.size == 2:
			pylab.plot(stimin[0]*numpy.ones(2),numpy.array([stimin[1]-0.4,stimin[1]+0.4]),'k',linewidth=2)
		else:
			pylab.plot(stimin[j,0]*numpy.ones(2),numpy.array([stimin[j,1]-0.4,stimin[j,1]+0.4]),'k',linewidth=2)
	pylab.xticks([])
	pylab.yticks([])
	pylab.xlim(0.,dur)
	pylab.ylabel(str(step*k+1))
	pylab.axis([0,dur,-0.5,19-0.5])
	if k == 0:
		pylab.title('Pattern'+', '+idname)		
	if k == 12:
		pylab.xticks(range(0,201,100))
		pylab.xlabel('Time (ms)')
pylab.xticks(range(0,201,50))
pylab.savefig('PP_'+idname+'.eps')

###############
pylab.figure(figsize=(15,16))
for k in range(4):
	pylab.subplot(4,1,k+1)
	dgsp = numpy.loadtxt('DGsp-'+str(step*k)+'-'+str(nid)+idname+'.txt')
	for j in range(dgsp.size/2):
		if dgsp.size == 2:
			pylab.plot(dgsp[0]*numpy.ones(2),numpy.array([dgsp[1]-0.4,dgsp[1]+0.4]),color,linewidth=2)
		else:
			pylab.plot(dgsp[j,0]*numpy.ones(2),numpy.array([dgsp[j,1]-0.4,dgsp[j,1]+0.4]),'b',linewidth=2)		# color for GC
	pylab.xticks([])
	pylab.yticks([])
	pylab.xlim(0.,dur)
	pylab.ylabel(str(step*k+1))
	pylab.axis([0,dur,-0.5,500-0.5])
	if k == 0:
		pylab.title('GC'+', '+idname)		
	if k == 12:
		pylab.xticks(range(0,201,100))
		pylab.xlabel('Time (ms)')
pylab.xticks(range(0,201,50))
pylab.savefig('GCoutput_'+idname+'.eps')
###############
pylab.show()
