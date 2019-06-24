### Analysis of DG network data ###
# This Python code extracts and plots the inputs to a selected GC.
# Enter the idname and cellID

# ModelDB file along with publication:
# Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
# http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

# modified and augmented by
# Man Yi Yim / 2015
# Alexander Hanuschkin / 2011

import pylab
import numpy
import matplotlib as mpl

# Font size in the figure
font_size = 12	# 20
mpl.rcParams['axes.titlesize'] = font_size+10
mpl.rcParams['xtick.labelsize'] = font_size+6
mpl.rcParams['ytick.labelsize'] = font_size+6
mpl.rcParams['axes.labelsize'] = font_size+8
mpl.rcParams['legend.fontsize'] = font_size
mpl.rcParams['font.size'] = font_size+10

# Parameters
tstop = 200	# duration of the simulation to be analysed
is_all_Vt = 1	# 1 if the dgvtfile contains all GCs
cellID = 8
npp = 100	# number of PP inputs
sid = 0		# PP_box_start_
nid = 6		# number of active PPs
#idname = "-pp10-gaba1-kir1-st0"
print 'duration = ' + str(tstop)
print 'npp = ' + str(npp)
print 'idname = ' + idname

# Import data
stiminfile = 'StimIn-'+str(sid)+'-'+str(nid)+idname+'.txt'
dgspfile = 'DGsp-'+str(sid)+'-'+str(nid)+idname+'.txt'
dgvtfile = 'DGVt-'+str(sid)+'-'+str(nid)+idname+'.txt'
con_file = 'DGNC-'+str(sid)+'-'+str(nid)+idname+'.txt'
con_file_id = 'id_DGNC.txt'

f2 = open(stiminfile,'r')
if f2.readline() == '':
	stimin = numpy.empty(0)
else:
	stimin = numpy.loadtxt(stiminfile)
f3 = open(dgspfile,'r')
if f3.readline() == '':
	dgsp = numpy.empty(0)
else:
	dgsp = numpy.loadtxt(dgspfile)
dgvt = numpy.loadtxt(dgvtfile)
g1 = open('130408_GC_and_input.txt','w')
g2 = open('130408_GC_mp.txt','w')

def extract_conn():
	"""extract_conn() extracts connectivity and prints the global IDs for Santhakumar model by default.
	Modify the number in this code if you change the network size or ordering.
	ID 0-499: GC
	ID 500-505: BC
	ID 506-520: MC
	ID 521-526: HC
	ID 527-626: PP
	ID >=527: BG	"""
	f_read = open(con_file)
	f_write = open(con_file_id,'w')
	for line in f_read:
		if line[0] == 'G':
			j = 13
			while line[j] != ']':
				j += 1
			pre = int(line[12:j])
		elif line[0] == 'B':
			j = 12
			while line[j] != ']':
				j += 1
			pre = 500 + int(line[11:j])
		elif line[0] == 'M':
			j = 11
			while line[j] != ']':
				j += 1
			pre = 506 + int(line[10:j])
		elif line[0] == 'H':
			j = 10
			while line[j] != ']':
				j += 1
			pre = 521 + int(line[9:j])
		elif line[0] == 'N':
			j = 9 #12
			while line[j] != ']':
				j += 1
			if line[7] == 'B':	# correlated PP
				pre = 527 + int(line[11:j])
			elif line[7] == '1':	# uncorrelated PP
				pre = 527 + int(line[11:j])
			elif line[7] == '[':	# BG
				pre = 527 + npp + int(line[8:j])
		else:
			print 'Please check the file: ' + con_file
		j += 2
		if line[j] == 'G':
			j += 13
			k = j
			while line[j] != ']':
				j += 1
			post = int(line[k-1:j])
		elif line[j] == 'B':
			j += 12
			k = j
			while line[j] != ']':
				j += 1
			post = 500 + int(line[k-1:j])
		elif line[j] == 'M':
			j += 11
			k = j
			while line[j] != ']':
				j += 1
			post = 506 + int(line[k-1:j])
		elif line[j] == 'H':
			j += 10
			k = j
			while line[j] != ']':
				j += 1
			post = 521 + int(line[k-1:j])
		else:
			print 'Please check the file: ' + con_file
		f_write.write(str(pre)+'\t'+str(post)+'\n')
	f_write.close()

#extract_conn()
conn = numpy.loadtxt(con_file_id)
pylab.figure(figsize=(10,8))
pylab.subplot(211)
ind = conn[conn[:,1]==cellID,0]
ind.sort()
j = 0
ystring = ['G'+str(cellID)]
if dgsp.size == 2:
	if dgsp[1] == cellID:
		pylab.plot(dgsp[0]*numpy.ones(2),numpy.array([j-0.4,j+0.4]),'b',linewidth=2)
elif dgsp.size > 2:
    if any(dgsp[:,1]==cellID):
        for spk in dgsp[dgsp[:,1]==cellID,0]:
            pylab.plot(spk*numpy.ones(2),numpy.array([j-0.4,j+0.4]),'b',linewidth=2)
spktime = dgsp[dgsp[:,1]==cellID,0]
for kk in range(spktime.size):
	g1.write(str(spktime[kk])+'\t'+str(cellID))
	g1.write('\n')
for k in ind:
	spktime = numpy.empty(0)
	if k < 527:
		if dgsp.size == 2 and dgsp[1] == k:
			spktime = numpy.append(spktime,dgsp[0])
		elif dgsp.size > 2:
			spktime = numpy.append(spktime,dgsp[dgsp[:,1]==k,0])
		if spktime.size > 0:
			j += 1
			if k < 500:
				color = 'b'
				ystring.append('G'+str(int(k))+r'$\to$G'+str(cellID))
			elif k < 506:
				color = 'g'
				ystring.append('B'+str(int(k-500))+r'$\to$G'+str(cellID))
			elif k < 521:
				color = 'r'
				ystring.append('M'+str(int(k-506))+r'$\to$G'+str(cellID))
			elif k < 527:
				color = 'c'
				ystring.append('H'+str(int(k-521))+r'$\to$G'+str(cellID))
			for m in spktime:
				pylab.plot(m*numpy.ones(2),numpy.array([j-0.4,j+0.4]),color,linewidth=2)
	else:
		if stimin.size == 2 and stimin[1] == k-527:
			spktime = numpy.append(spktime,stimin[0])
		elif stimin.size > 2:
			spktime = numpy.append(spktime,stimin[stimin[:,1]==k-527,0])
		if spktime.size > 0:
			j += 1
			ystring.append('PP'+str(int(k-527))+r'$\to$G'+str(cellID))
			for m in spktime:
				pylab.plot(m*numpy.ones(2),numpy.array([j-0.4,j+0.4]),'k',linewidth=2)
	for kk in range(spktime.size):
		g1.write(str(spktime[kk])+'\t'+str(int(k)))
		g1.write('\n')
pylab.yticks(numpy.arange(0.,j+1,1),ystring)
pylab.axis([0,tstop,-0.5,j+0.5])

pylab.subplot(212)
pylab.plot(dgvt[:,0],dgvt[:,cellID+1],'b',linewidth=2)
for kk in range(dgvt[:,0].size):
	g2.write(str(dgvt[kk,0])+'\t'+str(dgvt[kk,cellID+1]))
	g2.write('\n')
g1.close()
g2.close()
pylab.ylabel('Voltage (mV)')
pylab.xlabel('Time (ms)')
pylab.xlim(0,tstop)
pylab.subplots_adjust(bottom=0.08,left=0.145)
pylab.savefig('DG_Fig1C.eps')
pylab.show()
