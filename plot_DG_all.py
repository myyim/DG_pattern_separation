### Analysis of DG network data ###
# This Python code plots DG neurons' activity.
# Enter the idname, sid and nid.

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
selected_Vt = numpy.array([135,157,185,224,239])	# used if is_all_Vt is 0; put 5 IDs of GCs to print
npp = 100	# number of PP inputs
sid = 2		# PP_box_start_
nid = 6		# number of active PPs
#idname = "-pp10-gaba1-kir1-st0"
print 'duration = ' + str(tstop)
print 'npp = ' + str(npp)
print 'idname = ' + idname

# Import data
#bginfile = 'BGIn-'+str(sid)+'-'+str(nid)+idname+'.txt'
stiminfile = 'StimIn-'+str(sid)+'-'+str(nid)+idname+'.txt'
dgspfile = 'DGsp-'+str(sid)+'-'+str(nid)+idname+'.txt'
dgvtfile = 'DGVt-'+str(sid)+'-'+str(nid)+idname+'.txt'
con_file = 'DGNC-'+str(sid)+'-'+str(nid)+idname+'.txt'
con_file_id = 'id_DGNC.txt'

# Load files
#f1 = open(bginfile,'r')
#if f1.readline() == '':
#	bgin = numpy.empty(0)
#else:
#	bgin = numpy.loadtxt(bginfile)
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

def conn_mat(is_weighted = 1):
	"""conn_mat(is_weighted = 1) prints out the connectivity matrix according to con_file_id
	if is_weighted = 0, the matrix element is either connected (white) or unconnected (black)
	if is_weighted = 1, the matrix element is the connection strength that has to be entered manually.
	"""
	conn = numpy.loadtxt(con_file_id)
	ncell = 527
	npp = 100
	# connection strength
	gc2gc = 0
	gc2bc = 4.7e-3
	gc2mc = 0.2e-3
	gc2hc = 0.5e-3
	bc2gc = 1.6e-3
	bc2bc = 7.6e-3
	bc2mc = 1.5e-3
	bc2hc = 0
	mc2gc = 0.3e-3
	mc2bc = 0.3e-3
	mc2mc = 0.5e-3
	mc2hc = 0.2e-3
	hc2gc = 0.5e-3
	hc2bc = 0.5e-3
	hc2mc = 1.5e-3
	hc2hc = 0
	pp2gc = 0.02
	pp2bc = 0
	if is_weighted:
		cmat = numpy.zeros((ncell+npp,ncell))
		for j in range(conn.size/2):
			if conn[j,0] < 500:
				if conn[j,1] < 500:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = gc2gc
				elif conn[j,1] < 506:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = gc2bc
				elif conn[j,1] < 521:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = gc2mc
				elif conn[j,1] < 527:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = gc2hc
			elif conn[j,0] < 506:
				if conn[j,1] < 500:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = bc2gc
				elif conn[j,1] < 506:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = bc2bc
				elif conn[j,1] < 521:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = bc2mc
				elif conn[j,1] < 527:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = bc2hc
			elif conn[j,0] < 521:
				if conn[j,1] < 500:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = mc2gc
				elif conn[j,1] < 506:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = mc2bc
				elif conn[j,1] < 521:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = mc2mc
				elif conn[j,1] < 527:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = mc2hc
			elif conn[j,0] < 527:
				if conn[j,1] < 500:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = hc2gc
				elif conn[j,1] < 506:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = hc2bc
				elif conn[j,1] < 521:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = hc2mc
				elif conn[j,1] < 527:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = hc2hc
			elif conn[j,0] < 527 + npp:
				if conn[j,1] < 500:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = pp2gc
				elif conn[j,1] < 506:
					cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = pp2bc
		# The first figure is the connectivity matrix of projections from HC, MC and BC to all
		pylab.figure()
		pylab.pcolor(cmat[npp:ncell+npp-500,:],cmap=pylab.cm.spectral)
		pylab.colorbar()
		pylab.title('Connectivity matrix')
		pylab.axis([0,ncell,0,27])
		pylab.xlabel('Post')
		pylab.ylabel('Pre')
		pylab.yticks([3,13.5,24],('HC','MC','BC'))
		pylab.xticks([250,503,513.5,523],('GC','BC','MC','HC'))
		pylab.savefig('120412_conn_mat_to_all.png')
		# The second figure is the connectivity matrix of projections from all to HC, MC and BC
		pylab.figure()
		pylab.pcolor(cmat[npp:,500:],cmap=pylab.cm.spectral)
		pylab.colorbar()
		pylab.title('Connectivity matrix')
		pylab.axis([0,27,0,ncell])
		pylab.xlabel('Post')
		pylab.ylabel('Pre')
		pylab.yticks([3,13.5,24,177],('HC','MC','BC','GC'))
		pylab.xticks([3,13.5,23],('BC','MC','HC'))
		pylab.savefig('120412_conn_mat_from_all.png')
	else:
		cmat = numpy.ones((ncell+npp,ncell))
		for j in range(conn.size/2):
			if conn[j,0] < ncell+npp:	# background input is not considered
				cmat[ncell+npp-conn[j,0]-1,conn[j,1]] = 0
		# This figure is the connectivity matrix of projections from PP and all neurons to all neurons
		pylab.figure(figsize=(16,12))
		pylab.pcolor(cmat,cmap=pylab.cm.gray)
		pylab.title('Connectivity matrix: white means connected')
		pylab.axis([0,ncell,0,ncell+npp])
		pylab.xlabel('Post')
		pylab.ylabel('Pre')
		pylab.yticks([50,103,113.5,124,377],('PP','HC','MC','BC','GC'))
		pylab.xticks([250,503,513.5,523],('GC','B','M','H'))
		pylab.savefig('120412_conn_mat.png')

# Load connectivity
extract_conn()	# extract the connectivity profile 'con_file_id' for further analysis
# conn_mat(1)	# print out the connectivity matrix figure
conn = numpy.loadtxt(con_file_id)
pylab.figure(figsize=(10,15))

# Plot PP activity
pylab.subplot(411)
pylab.plot([-1],[-1],linewidth=2)
pylab.plot([-1],[-1],linewidth=2)
pylab.plot([-1],[-1],linewidth=2)
pylab.plot([-1],[-1],linewidth=2)
pylab.plot([-1],[-1],'k',linewidth=2)
pylab.legend(['GC','BC','MC','HC','PP'])
for j in range(stimin.size/2):
	if stimin.size == 2:
		pylab.plot(stimin[0]*numpy.ones(2),numpy.array([stimin[1]-0.4,stimin[1]+0.4]),'k',linewidth=2)
	else:
		pylab.plot(stimin[j,0]*numpy.ones(2),numpy.array([stimin[j,1]-0.4,stimin[j,1]+0.4]),'k',linewidth=2)
#pylab.xticks(range(0,201,25))
pylab.axis([0,tstop,-0.5,20-0.5])
pylab.ylabel('PP')
pylab.title(idname)

# Plot DG activity
pylab.subplot(412)
if dgsp.size > 0:
	for j in range(dgsp.size/2):
		if dgsp.size == 2:
			if dgsp[1] < 500:
				color = 'b'
			elif dgsp[1] < 506:
				color = 'g'
			elif dgsp[1] < 521:
				color = 'r'
			elif dgsp[1] < 527:
				color = 'c'
			pylab.plot(dgsp[0]*numpy.ones(2),numpy.array([dgsp[1]-0.4,dgsp[1]+0.4]),color,linewidth=2)
		else:
			if dgsp[j,1] < 500:
				color = 'b'
			elif dgsp[j,1] < 506:
				color = 'g'
			elif dgsp[j,1] < 521:
				color = 'r'
			elif dgsp[j,1] < 527:
				color = 'c'
			pylab.plot(dgsp[j,0]*numpy.ones(2),numpy.array([dgsp[j,1]-0.4,dgsp[j,1]+0.4]),color,linewidth=2)
#pylab.xticks(range(0,201,25))
pylab.axis([0,tstop,-0.5,527-0.5])
pylab.ylabel('DG')

# Plot spiking GCs and their inputs
nspk = 0
if dgsp.size == 2 and dgsp[1] < 500:
	nspk = 1
elif dgsp.size > 2:
	nspk = numpy.unique(dgsp[dgsp[:,1]<500,1]).size
if is_all_Vt:
    spkgc = numpy.arange(5)
    if dgsp.size > 0:
		if dgsp.size == 2 and dgsp[1] < 500:
			spkgc = dgsp[1]
		elif dgsp.size > 2:
			spkgc = numpy.unique(dgsp[dgsp[:,1]<500,1])
		if spkgc.size > 5:
			spkgc = spkgc[:5]
		elif spkgc.size < 5:
			for j in range(5-spkgc.size):
				spkgc = numpy.append(spkgc,[j])	# append a GC
else:
	spkgc = selected_Vt
spkgc.sort()
print 'selected neurons = ' + str(spkgc)

# Rasterplot of a few GCs and their inputs
pylab.subplot(413)
for j in range(dgsp.size/2):
	if dgsp.size == 2:
		if any(spkgc == dgsp[1]):
			n = pylab.argwhere(spkgc==dgsp[1])[0]
			pylab.plot(dgsp[0]*numpy.ones(2),numpy.array([n*2-0.4,n*2+0.4])+1,'b',linewidth=2)
	else:
		if any(spkgc == dgsp[j,1]):
			n = pylab.argwhere(spkgc==dgsp[j,1])[0]
			pylab.plot(dgsp[j,0]*numpy.ones(2),numpy.array([n*2-0.4,n*2+0.4])+1,'b',linewidth=2)
for j in range(spkgc.size):
	ind = conn[conn[:,1]==spkgc[j],0]
	for k in ind:
		if k < 527:
			for l in range(dgsp.size/2):
				if dgsp.size == 2:
					if dgsp[1] == k:
						if k < 500:
							color = 'b'
						elif k < 506:
							color = 'g'
						elif k < 521:
							color = 'r'
						elif k < 527:
							color = 'c'
						pylab.plot(dgsp[0]*numpy.ones(2),numpy.array([j*2-0.4,j*2+0.4])+2,color,linewidth=2)
				else:
					if dgsp[l,1] == k:
						if k < 500:
							color = 'b'
						elif k < 506:
							color = 'g'
						elif k < 521:
							color = 'r'
						elif k < 527:
							color = 'c'
						pylab.plot(dgsp[l,0]*numpy.ones(2),numpy.array([j*2-0.4,j*2+0.4])+2,color,linewidth=2)
		else:
			for l in range(stimin.size/2):
				if stimin.size == 2:
					if stimin[1] == k-527:
						pylab.plot(stimin[0]*numpy.ones(2),numpy.array([j*2-0.4,j*2+0.4])+2,'k',linewidth=2)
				else:
					if stimin[l,1] == k-527:
						pylab.plot(stimin[l,0]*numpy.ones(2),numpy.array([j*2-0.4,j*2+0.4])+2,'k',linewidth=2)
#			for l in range(bgin.size/2):
#				if bgin.size == 2:
#					if bgin[1] == k-527-npp:
#						pylab.plot(bgin[0]*numpy.ones(2),numpy.array([j*2-0.4,j*2+0.4])+2,'0.5',linewidth=2)
#				else:
#					if bgin[l,1] == k-527-npp:
#						pylab.plot(bgin[l,0]*numpy.ones(2),numpy.array([j*2-0.4,j*2+0.4])+2,'0.5',linewidth=2)
pylab.plot([0,tstop],[2.5,2.5],'k')
pylab.plot([0,tstop],[4.5,4.5],'k')
pylab.plot([0,tstop],[6.5,6.5],'k')
pylab.plot([0,tstop],[8.5,8.5],'k')
pylab.yticks(numpy.arange(1.5,10,2),spkgc.astype(int))
#pylab.xticks(range(0,201,25))
pylab.axis([0,tstop,0.5,spkgc.size*2+0.5])
pylab.ylabel('GCs & input')	

# Plot membrane potentials
pylab.subplot(414)
for j in range(spkgc.size):
#	if spkgc[j] > 5:
	pylab.plot(dgvt[:,0],dgvt[:,spkgc[j]+1],'b',linewidth=2)
#	else:
#		pylab.plot(dgvt[:,0],dgvt[:,spkgc[j]+1],'m',linewidth=2)
#	if any(dgsp[:,1] == spkgc[j]):
#		pylab.plot(dgvt[:,0],dgvt[:,spkgc[j]+1],'b',linewidth=2)
#	else:
#		pylab.plot(dgvt[:,0],dgvt[:,spkgc[j]+1],'m',linewidth=2)
pylab.ylabel('Potential (mV)')
#pylab.xticks(range(0,201,25))
pylab.xlabel('Time (ms) - '+str(nspk)+' GCs spike')
#pylab.xlabel('Time (ms)')
pylab.xlim(0,tstop)
pylab.savefig('DG-'+str(sid)+'-'+str(nid)+idname+'.eps')
pylab.show()
