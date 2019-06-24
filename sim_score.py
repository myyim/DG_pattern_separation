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
delta = numpy.array([20])/2.
dur = 200.	# duration of the simulation
bin = 1.
ninput = 100	# number of input sources
ncell = 500	# number of neurons (GCs)
#idname ='-pp10-gaba1-kir1-st0'

def corr_score(file1,file2,delta,bin=1.,dur=100.,ncell=500):
	"""Similarity score by correlation coefficient. The spike trains are convolved with a triangular kernel."""
	d1 = numpy.loadtxt(file1)
	d2 = numpy.loadtxt(file2)
	x = numpy.zeros(int(ncell*dur/bin))
	y = numpy.zeros(int(ncell*dur/bin))
	for j in range(ncell):
		if d1.size == 2:
			s1 = numpy.array(d1[0]*(d1[1]==j))
		else:
			s1 = d1[d1[:,1]==j,0]
		if d2.size == 2:
			s2 = numpy.array(d2[0]*(d2[1]==j))
		else:
			s2 = d2[d2[:,1]==j,0]
		kern = numpy.append(numpy.arange(delta/bin),numpy.arange(delta/bin,-1,-1))
		ts1,dump = pylab.histogram(s1,numpy.arange(0.,dur+bin,bin))
		ts2,dump = pylab.histogram(s2,numpy.arange(0.,dur+bin,bin))
		x[j*dur/bin:(j+1)*dur/bin] = numpy.convolve(ts1,kern,'same')
		y[j*dur/bin:(j+1)*dur/bin] = numpy.convolve(ts2,kern,'same')
        x = x - pylab.mean(x)
        y = y - pylab.mean(y)
        cor = sum(x*y)/(len(x)*pylab.std(x)*pylab.std(y))
        return cor

pylab.figure()
fin1 = 'StimIn-'+str(sid)+'-'+str(nid)+idname+'.txt'
fout1 = 'GCsp-'+str(sid)+'-'+str(nid)+idname+'.txt'
inscore = numpy.empty((delta.size,sidm*(sidm+1)/2+1))
outscore = numpy.empty((delta.size,sidm*(sidm+1)/2+1))
for k in range(delta.size):
	inscore[k,0] = corr_score(fin1,fin1,delta[k],bin,dur,ninput)
	outscore[k,0] = corr_score(fout1,fout1,delta[k],bin,dur,ncell)
count = 1
for j in range(0,sidm):
	fin1 = 'StimIn-'+str(j)+'-'+str(nid)+idname+'.txt'
	fout1 = 'GCsp-'+str(j)+'-'+str(nid)+idname+'.txt'
	for m in range(j+1,sidm+1):
		fin2 = 'StimIn-'+str(m)+'-'+str(nid)+idname+'.txt'
		fout2 = 'GCsp-'+str(m)+'-'+str(nid)+idname+'.txt'
		for k in range(delta.size):
			inscore[k,count] = corr_score(fin1,fin2,delta[k],bin,dur,ninput)
			outscore[k,count] = corr_score(fout1,fout2,delta[k],bin,dur,ncell)
		count += 1

for k in range(delta.size):
	pylab.plot(inscore[k,:],outscore[k,:],'ko',label='width='+str(2*delta[k]))
pylab.plot([0,1],[0,1],'k--',linewidth=2)
pylab.axis([-0.05,1.05,-0.05,1.05])
pylab.xlabel('Input similarity score')
pylab.ylabel('Output similarity score')
pylab.savefig('sim_score_'+idname+'.eps')
f_write = open('sim_score'+idname+'.txt','w')
for j in range(inscore.size):
	f_write.write(str(inscore[0][j])+'\t'+str(outscore[0][j])+'\n')
f_write.close()
pylab.show()
