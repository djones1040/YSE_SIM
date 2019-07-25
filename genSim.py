#!/usr/bin/env python
from __future__ import print_function

# current sim list, diff 1-day fracs
# python genSim.py --gcadence 3 --rcadence 9 --icadence 9 --zcadence 9 --roffset 0 --ioffset 3 --zoffset 6 --sim -i snanainputs/yse_gr_gi_gz_10per.input --ztfoffset 0.875 --onedayfrac 0.1
# python genSim.py --gcadence 3 --rcadence 9 --icadence 9 --zcadence 9 --roffset 0 --ioffset 3 --zoffset 6 --sim -i snanainputs/yse_gr_gi_gz_20per.input --ztfoffset 0.875 --onedayfrac 0.2
# python genSim.py --gcadence 3 --rcadence 9 --icadence 9 --zcadence 9 --roffset 0 --ioffset 3 --zoffset 6 --sim -i snanainputs/yse_gr_gi_gz_30per.input --ztfoffset 0.875 --onedayfrac 0.3



# gr, gi, gz every 3 days
# ZTF the day after
# python genSim.py --gcadence 3 --rcadence 9 --icadence 9 --zcadence 9 --roffset 0 --ioffset 3 --zoffset 6 --dosim -i snanainputs/yse_gr_gi_gz.input --ztfoffset 0.875 --onedayfrac 1.0

# gr, gi, gz every 3 days
# python genSim.py --gcadence 3 --rcadence 9 --icadence 9 --zcadence 9 --roffset 0 --ioffset 3 --zoffset 6 --dosim -i snanainputs/yse_gr_gi_gz.input
# griz + ZTF
# python genSim.py --dosim -i snanainputs/yse_ztf.input
# python genSim.py --dosim -i snanainputs/yse_ztf_gri.input --zcadence 1500
# python genSim.py --gcadence 3 --rcadence 9 --icadence 9 --zcadence 9 --roffset 0 --ioffset 3 --zoffset 6 --dosim -i snanainputs/yse_ztf_gr_gi_gz.input
# python genSim.py --dosim -i snanainputs/yse_ztf_perfect.input --perfect

# rsync -avz *.py djones1741@midway.rcc.uchicago.edu:/project/rkessler/djones/YSE/
# rsync -avz djones1741@midway.rcc.uchicago.edu:/project/rkessler/djones/YSE/YSE/* YSE/
# rsync -avz djones1741@midway.rcc.uchicago.edu:/project/rkessler/SN/SNDATA_ROOT/SIM/yse_YOUNG/* $SNDATA_ROOT/SIM/yse_YOUNG/
# rsync -avz djones1741@midway.rcc.uchicago.edu:/project/rkessler/SN/SNDATA_ROOT/SIM/yse/* $SNDATA_ROOT/SIM/yse/
# rsync -avz djones1741@midway.rcc.uchicago.edu:/project/rkessler/SN/SNDATA_ROOT/SIM/yse* $SNDATA_ROOT/SIM/
# rsync -avz djones1741@midway.rcc.uchicago.edu:/project/rkessler/SN/SNDATA_ROOT/SIM/yse/* $SNDATA_ROOT/SIM/yse/

import numpy as np
import random
from astropy.time import Time

haleakaladict = {'01':0.323,'02':0.357,'03':0.355,'04':0.333,'05':0.355,'06':0.267,
				 '07':0.226,'08':0.226,'09':0.367,'10':0.323,'11':0.333,'12':0.258}
palomardict = {'01':0.452,'02':0.607,'03':0.452,'04':0.300,'05':0.161,'06':0.067,
			   '07':0.097,'08':0.161,'09':0.167,'10':0.226,'11':0.300,'12':0.387}

def mkGoodWeatherList(simlibfile = 'PS1MD.simlib',
					  obslistfile='weather/yse_goodweather.list'):

	fin = open(simlibfile,'r')
	mjdlist = []
	for line in fin:
		if line.startswith('S:'):
			mjdlist += [int(float(line.split()[1]))]
	fin.close()
	mjdlist = np.sort(np.unique(mjdlist))

	fout = open(obslistfile,'w')
	for m in mjdlist:
		if m > 55317:
			print(m+2920,file=fout)
	fout.close()
	
class mkSimlibs:
	def __init__(self):
		pass

	def add_options(self, parser=None, usage=None, config=None):
		import optparse
		if parser == None:
			parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

		parser.add_option(
			'-g','--gcadence', default=3,type="float",
			help='cadence in g (default=%default)')
		parser.add_option(
			'-r','--rcadence', default=3,type="float",
			help='cadence in r (default=%default)')
		parser.add_option(
			'-i','--icadence', default=3,type="float",
			help='cadence in i (default=%default)')
		parser.add_option(
			'-z','--zcadence', default=3,type="float",
			help='cadence in z (default=%default)')
		parser.add_option(
			'--goffset', default=0,type="float",
			help='cadence offset in g (default=%default)')
		parser.add_option(
			'--roffset', default=0,type="float",
			help='cadence offset in r (default=%default)')
		parser.add_option(
			'--ioffset', default=0,type="float",
			help='cadence offset in i (default=%default)')
		parser.add_option(
			'--zoffset', default=0,type="float",
			help='cadence offset in z (default=%default)')
		parser.add_option(
			'--ztfoffset', default=-0.125,type="float",
			help='cadence offset for ZTF (default=%default)')
		parser.add_option(
			'--onedayfrac', default=0.1,type="float",
			help='fraction of survey with a one-day cadence, w/ offsets also divided by 3 for this subset (default=%default)')
		parser.add_option(
			'--simlibfile', default='snanasimlibs/yse.simlib',type="string",
			help='name of simlib file (default=%default)')
		parser.add_option(
			'-i','--inputfile', default='snanainputs/yse.input',type="string",
			help='name of simlib file (default=%default)')
		parser.add_option('--batchtmpl',default='/home/djones1741/djones/SBATCH_sandyb.TEMPLATE',
						  type="string",help='cluster batch template')
		parser.add_option(
			'-s','--sim', default=False,action="store_true",
			help='run SNANA simulation, if set (default=%default)')
		parser.add_option(
			'--fit', default=False,action="store_true",
			help='run SNANA fitting, if set (default=%default)')
		parser.add_option(
			'--analyze', default=False,action="store_true",
			help='run basic analysis tools, if set (default=%default)')
		parser.add_option(
			'--perfect', default=False,action="store_true",
			help='generate perfect high cadence simulations, if set (default=%default)')

		return(parser)

	def mksimlib(self,gcadence,rcadence,icadence,zcadence,
				 goffset,roffset,ioffset,zoffset,
				 simlibfile,simperfect=False,ztfsim=True,
				 ztf_offset=None,onedayfrac=0):
		# 3 day cadence, gr gi gz
		glines,rlines,ilines,zlines = getlines()
		mjd_goodweather = np.loadtxt('weather/yse_goodweather.list',unpack=True)

		# ZTF is 3 hours earlier
		#ztf_offset = -0.125

		if not ztfsim:
			print('Warning : not simulating ZTF!!')
		
		def mjd_to_month(mjd):
			date = Time(mjd,format='mjd').isot
			return date.split('-')[1]

		if not simperfect:
			mjd = np.arange(58240,59517,1)
		else:
			mjd = np.arange(58240,59335,1)

		# nominal survey
		count = 0; nightcount = -1; usednightcount = 0; ps1count = 0
		simliblines = []
		for m in mjd[0:int(len(mjd)*(1-onedayfrac))]:
			nightcount += 1

			if ztf_offset < 0:
				randval = random.uniform(0, 1)
				if ztfsim and (randval > palomardict[mjd_to_month(m)] or simperfect):
					if not nightcount % 3 or simperfect:
						iLine = random.sample(range(len(glines)),1)[0]
						for lines in [rlines]:
							line = lines[iLine][0]
							linemjd = line.split()[1]
							lineid = line.split()[2]
							simliblines += [line.replace(linemjd,'%.2f'%(m+ztf_offset)).replace(' %i '%int(lineid),' %i '%count).replace(' r ',' X ')]
							simliblines += [line.replace(linemjd,'%.2f'%(m+ztf_offset)).replace(' %i '%int(lineid),' %i '%count).replace(' r ',' Y ')]
							count += 2

			# replace random weather with PS1 MD observations
			# includes not much data at bright time
			randval = random.uniform(0, 1)
			if randval > haleakaladict[mjd_to_month(m)] or simperfect:
				iLine = random.sample(range(len(glines)),1)[0]
				if (not (nightcount - goffset) % gcadence and nightcount - goffset >= 0) or simperfect:
					for lines in [glines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						simliblines += [line.replace(linemjd,'%.2f'%m).replace(' %i '%int(lineid),' %i '%count)]
						count += 1
						ps1count += 1
				if (not (nightcount - roffset) % rcadence and nightcount - roffset >= 0) or simperfect:
					for lines in [rlines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						simliblines += [line.replace(linemjd,'%.2f'%m).replace(' %i '%int(lineid),' %i '%count)]
						count += 1
						ps1count += 1
				if (not (nightcount - ioffset) % icadence and nightcount - ioffset >= 0) or simperfect:
					for lines in [ilines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						simliblines += [line.replace(linemjd,'%.2f'%m).replace(' %i '%int(lineid),' %i '%count)]
						count += 1
						ps1count += 1
				if (not (nightcount - zoffset) % zcadence and nightcount - zoffset >= 0) or simperfect:
					for lines in [zlines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						simliblines += [line.replace(linemjd,'%.2f'%m).replace(' %i '%int(lineid),' %i '%count)]
						count += 1
						ps1count += 1
						
				if not nightcount % gcadence or not nightcount % rcadence or not \
				   nightcount % icadence or not nightcount % zcadence:
					usednightcount += 1

			if ztf_offset > 0:
				randval = random.uniform(0, 1)
				if ztfsim and (randval > palomardict[mjd_to_month(m)] or simperfect):
					if not nightcount % 3 or simperfect:
						iLine = random.sample(range(len(glines)),1)[0]
						for lines in [rlines]:
							line = lines[iLine][0]
							linemjd = line.split()[1]
							lineid = line.split()[2]
							simliblines += [line.replace(linemjd,'%.2f'%(m+ztf_offset)).replace(' %i '%int(lineid),' %i '%count).replace(' r ',' X ')]
							simliblines += [line.replace(linemjd,'%.2f'%(m+ztf_offset)).replace(' %i '%int(lineid),' %i '%count).replace(' r ',' Y ')]
							count += 2

					
		# one-day survey
		for m in mjd[int(len(mjd)*(1-onedayfrac)):]:
			nightcount += 1

			if ztf_offset < 0:
				randval = random.uniform(0, 1)
				if ztfsim and (randval > palomardict[mjd_to_month(m)] or simperfect):
					if not nightcount % 3 or simperfect:
						iLine = random.sample(range(len(glines)),1)[0]
						for lines in [rlines]:
							line = lines[iLine][0]
							linemjd = line.split()[1]
							lineid = line.split()[2]
							simliblines += [line.replace(linemjd,'%.2f'%(m+ztf_offset)).replace(' %i '%int(lineid),' %i '%count).replace(' r ',' X ')]
							simliblines += [line.replace(linemjd,'%.2f'%(m+ztf_offset)).replace(' %i '%int(lineid),' %i '%count).replace(' r ',' Y ')]
							count += 2
							
			# replace random weather with PS1 MD observations
			# includes not much data at bright time
			randval = random.uniform(0, 1)
			if randval > haleakaladict[mjd_to_month(m)] or simperfect:
				iLine = random.sample(range(len(glines)),1)[0]
				if (not (nightcount - goffset/3)%1 and nightcount - goffset/3 >= 0) or simperfect:
					for lines in [glines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						simliblines += [line.replace(linemjd,'%.2f'%m).replace(' %i '%int(lineid),' %i '%count)]
						count += 1
						ps1count += 1
				if (not (nightcount - roffset/3)%1 and nightcount - roffset/3 >= 0) or simperfect:
					for lines in [rlines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						simliblines += [line.replace(linemjd,'%.2f'%m).replace(' %i '%int(lineid),' %i '%count)]
						count += 1
						ps1count += 1
				if (not (nightcount - ioffset/3)%1 and nightcount - ioffset/3 >= 0) or simperfect:
					for lines in [ilines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						simliblines += [line.replace(linemjd,'%.2f'%m).replace(' %i '%int(lineid),' %i '%count)]
						count += 1
						ps1count += 1
				if (not (nightcount - zoffset/3)%1 and nightcount - zoffset/3 >= 0) or simperfect:
					for lines in [zlines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						simliblines += [line.replace(linemjd,'%.2f'%m).replace(' %i '%int(lineid),' %i '%count)]
						count += 1
						ps1count += 1
						
				if not nightcount % gcadence or not nightcount % rcadence or not \
				   nightcount % icadence or not nightcount % zcadence:
					usednightcount += 1

			if ztf_offset > 0:
				randval = random.uniform(0, 1)
				if ztfsim and (randval > palomardict[mjd_to_month(m)] or simperfect):
					if not nightcount % 3 or simperfect:
						iLine = random.sample(range(len(glines)),1)[0]
						for lines in [rlines]:
							line = lines[iLine][0]
							linemjd = line.split()[1]
							lineid = line.split()[2]
							simliblines += [line.replace(linemjd,'%.2f'%(m+ztf_offset)).replace(' %i '%int(lineid),' %i '%count).replace(' r ',' X ')]
							simliblines += [line.replace(linemjd,'%.2f'%(m+ztf_offset)).replace(' %i '%int(lineid),' %i '%count).replace(' r ',' Y ')]
							count += 2

				
		fout = open(simlibfile,'w')
		print(simlibheader%count,file=fout)
		for l in simliblines:
			print(l,file=fout)
		print(simlibfooter,file=fout)
		fout.close()

		surveyarea = 2.*3600./20*7/(ps1count/float(usednightcount))*(np.pi/180.)**2.
		return surveyarea
	
		
	def mkinput(self,gcadence,rcadence,icadence,zcadence,inputfile,simlibfile,surveyarea,simperfect=False,batchtmpl=None):

		filtstr = ''
		for filt,cadence in zip('grizXY',[gcadence,rcadence,icadence,zcadence,3,3]):
			if cadence: filtstr += filt

		fout = open(inputfile.replace('.','_ia.'),'w')
		print(iainputheader%(
			inputfile.split('/')[-1].split('.')[0]+'_ia',simlibfile,filtstr,surveyarea),file=fout)
		if simperfect:
			print('GENPERFECT: 1',file=fout)
		fout.close()

		fout = open(inputfile.replace('.','_nonia.'),'w')
		print(noniainputheader%(
			inputfile.split('/')[-1].split('.')[0]+'_nonia',simlibfile,filtstr,surveyarea),file=fout)
		if simperfect:
			print('GENPERFECT: 1',file=fout)
		fout.close()

		#fout = open(inputfile.replace('.','_young.'),'w')
		#print(noniainputheader_young%(
		#	inputfile.split('/')[-1].split('.')[0]+'_young',simlibfile,filtstr,surveyarea),file=fout)
		#fout.close()

		if simperfect:
			zmax = 0.2
			cutwinstr = 'GENOPT: CUTWIN_SNRMAX 50.0 grizXY 2 -20. 80.'
			ratestr = 'GENOPT(NON1A): DNDZ_PEC1A POWERLAW 1E-5 2.15\nGENOPT(NON1A): DNDZ POWERLAW 1E-4   4.5'
		else:
			zmax = 0.5
			cutwinstr = ''
			ratestr = ''
		
		fout = open(inputfile.replace('.','_MASTER.'),'w')
		genversion = inputfile.split('/')[-1].split('.')[0]
		plasticc = plasticcmodelstr%(
			genversion,genversion,genversion,genversion,
			genversion,genversion,genversion,genversion,
			genversion,genversion,genversion,genversion,
			genversion)
		
		print(mastertmpl%(
			batchtmpl,
			inputfile.split('/')[-1].split('.')[0],inputfile.split('/')[-1].split('.')[0],
			cutwinstr,ratestr,plasticc,inputfile.replace('.','_ia.'),inputfile.replace('.','_nonia.'),zmax,
			inputfile.split('/')[-1].split('.')[0]),file=fout)
		fout.close()	

		return genversion
		
def getlines(simlibfile = 'snanasimlibs/found_yse.simlib'):


	glines,rlines,ilines,zlines = [],[],[],[]
	gmjd = rmjd = imjd = zmjd = np.array([])
	glib = rlib = ilib = zlib = np.array([])
	fin = open(simlibfile,'r')
	lines = fin.readlines()
	fin.close()
	for line in lines:
		if line.startswith('LIBID:'):
			libid = int(line.split()[1])
		elif line.startswith('S:'):
			line = line.replace('\n','')
			filt = line.split()[3]
			if filt == 'g':
				gmjd = np.append(gmjd,[float(line.split()[1])])
				glib = np.append(glib,[libid])
				glines.append([line])
			elif filt == 'r':
				rmjd = np.append(rmjd,[float(line.split()[1])])
				rlib = np.append(rlib,[libid])
				rlines.append([line])
			elif filt == 'i':
				imjd = np.append(imjd,[float(line.split()[1])])
				ilib = np.append(ilib,[libid])
				ilines.append([line])
			elif filt == 'z':
				zmjd = np.append(zmjd,[float(line.split()[1])])
				zlib = np.append(zlib,[libid])
				zlines.append([line])

	outlinesg,outlinesr,outlinesi,outlinesz = [],[],[],[]
	for gm,gl,gli in zip(gmjd,glib,glines):
		rMatch = np.where((np.abs(gm-rmjd) < 0.01) & (rlib == gl))[0]
		iMatch = np.where((np.abs(gm-imjd) < 0.01) & (ilib == gl))[0]
		zMatch = np.where((np.abs(gm-zmjd) < 0.01) & (zlib == gl))[0]
		if len(rMatch) == 1 and len(iMatch) == 1 and len(zMatch) == 1:
			outlinesg.append(gli)
			outlinesr.append(rlines[int(rMatch)])
			outlinesi.append(ilines[int(iMatch)])
			outlinesz.append(zlines[int(zMatch)])
			
	return(outlinesg,outlinesr,outlinesi,outlinesz)

iainputheader = """
GENVERSION: %s         # simname
GENSOURCE:  RANDOM   
GENMODEL:   SALT2.Guy10_UV2IR
GENPREEFIX: YSE_IA

SIMLIB_FILE: %s # simlib file

CIDOFF: 500
KCOR_FILE:  $PS1_ROOT/kcor/ZTF/kcor_PS1_ZTF_none.fits
APPLY_SEARCHEFF_OPT: 0

EXPOSURE_TIME_FILTER: g 1.0
EXPOSURE_TIME_FILTER: r 1.0
EXPOSURE_TIME_FILTER: i 1.0
EXPOSURE_TIME_FILTER: z 1.0

GENMAG_SMEAR_MODELNAME: G10
# selection criteria for generation
GENFILTERS:       %s

GENSIGMA_SEARCH_PEAKMJD:  1.0         # sigma-smearing for  SEARCH_PEAKMJD (days)

GENRANGE_PEAKMJD:  58240  59617
SOLID_ANGLE: %.3f # 0.148 # 1 field, 7 sq degreees *7
# baseline for 4 filters should be 630 degrees (0.192 steradians)

SEARCHEFF_PIPELINE_FILE:  SEARCHEFF_PIPELINE_YSE.DAT
SEARCHEFF_PIPELINE_LOGIC_FILE:  SEARCHEFF_PIPELINE_LOGIC_YSE.DAT

GENRANGE_REDSHIFT:  0.001    0.5
GENSIGMA_REDSHIFT:  0.000001
DNDZ: POWERLAW  2.6E-5  2.2
GENRANGE_TREST:   -100.0    80.0     # rest epoch relative to peak (days)

GENMEAN_RV:         3.1               # mean RV to generate

OPT_MWEBV: 1

RANSEED: 128473       # random number seed

# smear flags: 0=off, 1=on
SMEARFLAG_FLUX:    1  # photo-stat smearing of signal, sky, etc ...
SMEARFLAG_ZEROPT:  1  # smear zero-point with zptsig

# SEARCHEFF_SPEC_FILE:  SEARCHEFF_SPEC_MF_final.DAT #this is the big one for now.
APPLY_CUTWIN_OPT:     1
CUTWIN_NEPOCH:   5 -5.              # require 5 epochs (no S/N requirement)
CUTWIN_TRESTMIN: -20  10
CUTWIN_TRESTMAX:   9  40
CUTWIN_MWEBV:      0 .20

FORMAT_MASK:  2 # terse format
CUTWIN_SNRMAX:   5.0 grizXY 2 -20. 80.  # require 1 of griz with S/N > 5

GENMEAN_SALT2x1:     0.703
GENRANGE_SALT2x1:   -5.0  +4.0     # x1 (stretch) range
GENSIGMA_SALT2x1:    2.15  0.472      # bifurcated sigmas

GENMEAN_SALT2c:     -0.04
GENRANGE_SALT2c:   -0.4   0.4     # color range
GENSIGMA_SALT2c:    0.033   0.125     # bifurcated sigmas

# SALT2 alpha and beta

GENMEAN_SALT2ALPHA:   0.14
GENMEAN_SALT2BETA:   3.1

# cosmological params for lightcurve generation and redshift distribution
OMEGA_MATTER:  0.3
OMEGA_LAMBDA:  0.7
W0_LAMBDA:    -1.00
H0:           70.0   

SIMGEN_DUMP:  10  CID  Z  PEAKMJD S2c S2x1 SNRMAX MAGT0_r MAGT0_g MJD_TRIGGER NON1A_INDEX
"""

noniainputheader = """
GENVERSION: %s         # simname
GENSOURCE:  RANDOM   
GENMODEL:   NON1A
GENPREEFIX: YSE_IA

SIMLIB_FILE: %s # simlib file

CIDOFF: 500
KCOR_FILE:  $PS1_ROOT/kcor/ZTF/kcor_PS1_ZTF_none.fits
APPLY_SEARCHEFF_OPT: 0

EXPOSURE_TIME_FILTER: g 1.0
EXPOSURE_TIME_FILTER: r 1.0
EXPOSURE_TIME_FILTER: i 1.0
EXPOSURE_TIME_FILTER: z 1.0

SEARCHEFF_PIPELINE_FILE:  SEARCHEFF_PIPELINE_YSE.DAT
SEARCHEFF_PIPELINE_LOGIC_FILE:  SEARCHEFF_PIPELINE_LOGIC_YSE.DAT

GENMAG_SMEAR_MODELNAME: G10
# selection criteria for generation
GENFILTERS:       %s

GENSIGMA_SEARCH_PEAKMJD:  1.0         # sigma-smearing for  SEARCH_PEAKMJD (days)

GENRANGE_PEAKMJD:  58240  59617
SOLID_ANGLE: %.3f # 0.148 # 1 field, 7 sq degreees *7
# baseline for 4 filters should be 630 degrees (0.192 steradians)

GENRANGE_REDSHIFT:  0.001    0.5
GENSIGMA_REDSHIFT:  0.000001
DNDZ: POWERLAW 5E-5   4.5
DNDZ_PEC1A: POWERLAW  2.6E-5  2.2
GENRANGE_TREST:   -100.0    80.0     # rest epoch relative to peak (days)

GENMEAN_RV:         3.1               # mean RV to generate

OPT_MWEBV: 1

RANSEED: 128473       # random number seed

# smear flags: 0=off, 1=on
SMEARFLAG_FLUX:    1  # photo-stat smearing of signal, sky, etc ...
SMEARFLAG_ZEROPT:  1  # smear zero-point with zptsig

# SEARCHEFF_SPEC_FILE:  SEARCHEFF_SPEC_MF_final.DAT #this is the big one for now.
APPLY_CUTWIN_OPT:     1
CUTWIN_NEPOCH:   5 -5.              # require 5 epochs (no S/N requirement)
CUTWIN_TRESTMIN: -20  10
CUTWIN_TRESTMAX:   9  40
CUTWIN_MWEBV:      0 .20

FORMAT_MASK:  2 # terse format
CUTWIN_SNRMAX:   5.0 grizXY 2 -20. 80.  # require 1 of griz with S/N > 5

GENMEAN_SALT2x1:     0.703
GENRANGE_SALT2x1:   -5.0  +4.0     # x1 (stretch) range
GENSIGMA_SALT2x1:    2.15  0.472      # bifurcated sigmas
#GENSIGMA_SALT2x1:    0.1  0.1      # bifurcated sigmas

GENMEAN_SALT2c:     -0.04
GENRANGE_SALT2c:   -0.4   0.4     # color range
GENSIGMA_SALT2c:    0.033   0.125     # bifurcated sigmas

# SALT2 alpha and beta

GENMEAN_SALT2ALPHA:   0.14
GENMEAN_SALT2BETA:   3.1

# cosmological params for lightcurve generation and redshift distribution
OMEGA_MATTER:  0.3
OMEGA_LAMBDA:  0.7
W0_LAMBDA:    -1.00
H0:           70.0   

SIMGEN_DUMP:  10  CID  Z  PEAKMJD S2c S2x1 SNRMAX MAGT0_r MAGT0_g MJD_TRIGGER NON1A_INDEX

INPUT_FILE_INCLUDE: LFs/SIMGEN_INCLUDE_NON1A_J17-beforeAdjust.INPUT
"""

noniainputheader_young = """
GENVERSION: %s         # simname
GENSOURCE:  RANDOM   
GENMODEL:   NON1A
GENPREEFIX: YSE_IA

SIMLIB_FILE: %s # simlib file

CIDOFF: 500
KCOR_FILE:  $PS1_ROOT/kcor/ZTF/kcor_PS1_ZTF_none.fits
APPLY_SEARCHEFF_OPT: 0

SEARCHEFF_PIPELINE_FILE:  SEARCHEFF_PIPELINE_YSE.DAT
SEARCHEFF_PIPELINE_LOGIC_FILE:  SEARCHEFF_PIPELINE_LOGIC_YSE.DAT

EXPOSURE_TIME_FILTER: g 1.0
EXPOSURE_TIME_FILTER: r 1.0
EXPOSURE_TIME_FILTER: i 1.0
EXPOSURE_TIME_FILTER: z 1.0

GENMAG_SMEAR_MODELNAME: G10
# selection criteria for generation
GENFILTERS:       %s

GENSIGMA_SEARCH_PEAKMJD:  1.0         # sigma-smearing for  SEARCH_PEAKMJD (days)

GENRANGE_PEAKMJD:  58240  59617
SOLID_ANGLE: %.3f # 0.148 # 1 field, 7 sq degreees *7
# baseline for 4 filters should be 630 degrees (0.192 steradians)

GENRANGE_REDSHIFT:  0.001    0.5
GENSIGMA_REDSHIFT:  0.000001
DNDZ: POWERLAW2 5E-5   4.5   0.0   0.8    # rate = R0(1+z)^Beta for z<0.8
DNDZ: POWERLAW2 5.44E-4  0.0   0.8   9.1  # rate = constant for z>0.8
DNDZ_PEC1A: POWERLAW  2.6E-5  2.2
GENRANGE_TREST:   -100.0    80.0     # rest epoch relative to peak (days)

GENMEAN_RV:         3.1               # mean RV to generate

OPT_MWEBV: 1

RANSEED: 128473       # random number seed

# smear flags: 0=off, 1=on
SMEARFLAG_FLUX:    1  # photo-stat smearing of signal, sky, etc ...
SMEARFLAG_ZEROPT:  1  # smear zero-point with zptsig

# SEARCHEFF_SPEC_FILE:  SEARCHEFF_SPEC_MF_final.DAT #this is the big one for now.
APPLY_CUTWIN_OPT:     1
CUTWIN_NEPOCH:   5 -5.              # require 5 epochs (no S/N requirement)
CUTWIN_TRESTMIN: -20  10
CUTWIN_TRESTMAX:   9  40
CUTWIN_MWEBV:      0 .20

FORMAT_MASK:  2 # terse format
CUTWIN_SNRMAX:   5.0 grizXY 2 -20. 80.  # require 1 of griz with S/N > 5

GENMEAN_SALT2x1:     0.703
GENRANGE_SALT2x1:   -5.0  +4.0     # x1 (stretch) range
GENSIGMA_SALT2x1:    2.15  0.472      # bifurcated sigmas
#GENSIGMA_SALT2x1:    0.1  0.1      # bifurcated sigmas

GENMEAN_SALT2c:     -0.04
GENRANGE_SALT2c:   -0.4   0.4     # color range
GENSIGMA_SALT2c:    0.033   0.125     # bifurcated sigmas

# SALT2 alpha and beta

GENMEAN_SALT2ALPHA:   0.14
GENMEAN_SALT2BETA:   3.1

# cosmological params for lightcurve generation and redshift distribution
OMEGA_MATTER:  0.3
OMEGA_LAMBDA:  0.7
W0_LAMBDA:    -1.00
H0:           70.0   

SIMGEN_DUMP:  10  CID  Z  PEAKMJD S2c S2x1 SNRMAX MAGT0_r MAGT0_g MJD_TRIGGER NON1A_INDEX

PATH_NON1ASED: /project/rkessler/djones/YSE/LFs/NON1A
INPUT_FILE_INCLUDE: SIMGEN_INCLUDE_NON1A_YOUNGSN.INPUT
"""


simlibheader = """SURVEY:		  PS1MD
FILTERS:		grizXY
PSF_UNIT:		ARCSEC_FWHM
SKYSIG_UNIT:	ADU_PER_SQARCSEC

# Assume instrument parameters for GROUND
# Assume SKYMAG(  2700.) = 23.80 mag/asec^2
# Assume SKYMAG(  3714.) = 22.70 mag/asec^2
# Assume SKYMAG(  4790.) = 21.00 mag/asec^2
# Assume SKYMAG(  6220.) = 20.40 mag/asec^2
# Assume SKYMAG(  7544.) = 19.40 mag/asec^2
# Assume SKYMAG(  8679.) = 18.10 mag/asec^2
# Assume SKYMAG( 10095.) = 17.90 mag/asec^2

FLUXERR_COR:  grizXY  -0.90      4.7703  1.0526  1.0000  0.8650  4.7703  1.0526
FLUXERR_COR:  grizXY  -0.70      5.3170  1.0000  1.0000  1.3651  5.3170  1.0000
FLUXERR_COR:  grizXY  -0.50      2.3181  1.7215  1.3091  0.7351  2.3181  1.7215
FLUXERR_COR:  grizXY  -0.30      3.3645  1.8474  1.1026  0.7446  3.3645  1.8474
FLUXERR_COR:  grizXY  -0.10      3.1311  2.0894  1.0315  1.0230  3.1311  2.0894
FLUXERR_COR:  grizXY   0.10      3.3423  2.1468  1.6439  1.0871  3.3423  2.1468
FLUXERR_COR:  grizXY   0.30      3.0158  2.8774  1.2655  1.1967  3.0158  2.8774
FLUXERR_COR:  grizXY   0.50      3.1574  2.2492  2.1569  1.2198  3.1574  2.2492
FLUXERR_COR:  grizXY   0.70      2.3857  2.4199  1.5264  1.2647  2.3857  2.4199
FLUXERR_COR:  grizXY   0.90      2.2685  2.1145  1.4510  1.2045  2.2685  2.1145
FLUXERR_COR:  grizXY   1.10      2.0670  1.8696  1.3235  1.0674  2.0670  1.8696
FLUXERR_COR:  grizXY   1.30      1.8561  1.7418  1.1890  1.1428  1.8561  1.7418
FLUXERR_COR:  grizXY   1.50      1.4960  1.4025  1.1767  1.0346  1.4960  1.4025
FLUXERR_COR:  grizXY   1.70      1.4736  1.2986  1.0032  1.2882  1.4736  1.2986
FLUXERR_COR:  grizXY   1.90      1.1383  1.1056  1.0742  1.0000  1.1383  1.1056
FLUXERR_COR:  grizXY   2.10      1.6320  1.0000  1.0000  1.0000  1.6320  1.0000
FLUXERR_COR:  grizXY   2.30      1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
FLUXERR_COR:  grizXY   2.50      1.0000  1.0000  1.0000  1.0000  1.0000  1.0000
FLUXERR_COR:  grizXY   2.70      1.0000  1.0000  1.0000  1.0000  1.0000  1.0000

BEGIN LIBGEN

#--------------------------------------------
LIBID:		  1		# cadence from 2016W^@
RA:	   37.665487		DECL:	42.235969	  MWEBV:  0.059
NOBS:	%i		PIXSIZE:  0.500		REDSHIFT: 0.01925	  PEAKMJD: 57417.953
SUBSURVEY: PS1MD
"""

simlibfooter = """
#  MJD	   IDUM	 BAND  GAIN RDNOISE	 SKYSIG	   PSF1 PSF2 PSFRAT	   ZP	ZPERR
END_LIBID:		1


END_OF_SIMLIB:
"""

mastertmpl = """
BATCH_INFO:  sbatch  %s 20

# nominal generation

GENVERSION: %s
# GENOPT(NON1A): DNDZ_PEC1A POWERLAW 3.9E-6 2.15
GENOPT(NON1A): DNDZ_PEC1A POWERLAW 5.5E-6 2.15
GENOPT(NON1A): PATH_NON1ASED /project/rkessler/djones/YSE/LFs/NON1A
# new rates, allowing SN Iax to be 31%%

GENVERSION: %s_YOUNG
GENOPT(NON1A): PATH_NON1ASED /project/rkessler/djones/YSE/LFs/NON1A
GENOPT(NON1A): INPUT_FILE_INCLUDE LFs/SIMGEN_INCLUDE_NON1A_YOUNGSN.INPUT
%s
%s

%s # PLASTICC Models

ENDLIST_GENVERSION:

NGEN_UNIT:  0.05  SEASONS

# specify sim-input files for snlc_sim.exe
SIMGEN_INFILE_Ia: %s
SIMGEN_INFILE_NONIa: %s

# define required global items to ensure uniformity among all jobs
H0: 70
ZRANGE:      0.002  %.1f
GENPREFIX:   %s          # prefix of all data filenames
FORMAT_MASK: 48           # 2=TERSE    16=RanCID  32=FITS-FORMAT
RESET_CIDOFF: 2

RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349
RANSEED: 12349

"""

plasticcmodelstr = """
# 91BG from S. Gonzalez-Gaitan and Felipe Lagos
#  (more templates than J17, and stretch-color correlation)
GENVERSION: %s_PLASTICC_MODEL67_SNIa-91bg
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_SNIa-91bg.INPUT
GENOPT: GENTYPE 67
GENOPT: SEARCHEFF_SPEC_SCALE 1.0


# SNIax from Saurabh
GENVERSION: %s_PLASTICC_MODEL52_SNIax
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_SNIax.INPUT
GENOPT: GENTYPE 52
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# Superluminous SN:  SLSN-I
GENVERSION:  %s_PLASTICC_MODEL95_SLSN-I
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_SLSN-I-MOSFIT.INPUT
GENOPT: GENTYPE 95
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# pair instability SN: PISN
GENVERSION:  %s_PLASTICC_MODEL99_PISN
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_PISN-MOSFIT.INPUT
GENOPT: GENTYPE 99
GENOPT: SEARCHEFF_SPEC_FILE ZERO

# Intermediate Luminosity Optical Transients (ILOT)
GENVERSION:  %s_PLASTICC_MODEL99_ILOT
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_ILOT-MOSFIT.INPUT
GENOPT: GENTYPE 99
GENOPT: SEARCHEFF_SPEC_FILE ZERO

# Ca Rich Transients (CART)
GENVERSION:  %s_PLASTICC_MODEL99_CART
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_CART-MOSFIT.INPUT
GENOPT: GENTYPE 99
GENOPT: SEARCHEFF_SPEC_FILE ZERO

# TDE
#GENVERSION:  %s_PLASTICC_MODEL15_TDE
#GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_TDE-MOSFIT.INPUT
#GENOPT: GENTYPE 15
##GENOPT: SEARCHEFF_SPEC_FILE ZERO
#GENOPT: SEARCHEFF_SPEC_SCALE 1.00

# MOSFIT-IIn
GENVERSION:  %s_PLASTICC_MODEL42_SNIIn
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_SNIIn-MOSFIT.INPUT
GENOPT: GENTYPE 42  SIMLIB_NREPEAT 2
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# MOSFIT-Ibc
GENVERSION:  %s_PLASTICC_MODEL62_SNIbc-MOSFIT
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_SNIbc-MOSFIT.INPUT
GENOPT: GENTYPE 62   SIMLIB_NREPEAT 4
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# Core collapse Type II using pca (5->12 on May 9 2018)
# for end-of-challenge model release
GENVERSION: %s_PLASTICC_MODEL42_SNII-NMF
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_SNII-NMF.INPUT
GENOPT: GENTYPE 42  SIMLIB_NREPEAT 4
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# - - - - - - - - - - - - - - - - - - - - - - - - -
# legacy NON1ASED
# Core collapse Type II from K10 templates
GENVERSION: %s_PLASTICC_MODEL42_SNII-Templates
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_SNII-Templates.INPUT
GENOPT: GENTYPE 42   SIMLIB_NREPEAT 4
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# NON1ASED-Ibc
# Core collapse Type Ibc from K10 templates
GENVERSION: %s_PLASTICC_MODEL62_SNIbc-Templates
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_SNIbc-Templates.INPUT
GENOPT: GENTYPE 62   SIMLIB_NREPEAT 4
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# - - - - - - - -
# Type Ia SN
GENVERSION: %s_PLASTICC_MODEL90_SNIa-SALT2
GENOPT: INPUT_FILE_INCLUDE $PLASTICC_ROOT/SIMGEN/SIMGEN_INCLUDE_SNIa-SALT2.INPUT
GENOPT: GENTYPE 90
GENOPT: SEARCHEFF_SPEC_SCALE 0.5

"""

if __name__ == "__main__":

	import os
	import optparse

	mks = mkSimlibs()

	usagestring = 'getSim.py <options>'
	parser = mks.add_options(usage=usagestring)
	options,  args = parser.parse_args()

	if not options.sim or not options.fit:
		surveyarea = mks.mksimlib(options.gcadence,options.rcadence,options.icadence,
								  options.zcadence,options.goffset,options.roffset,
								  options.ioffset,options.zoffset,
								  options.simlibfile.replace('.simlib','_%s.simlib'%options.inputfile.split('/')[-1].split('.')[0]),
								  simperfect=options.perfect,
								  ztf_offset=options.ztfoffset,
								  onedayfrac=options.onedayfrac)
		genversion = mks.mkinput(options.gcadence,options.rcadence,options.icadence,options.zcadence,
								 options.inputfile,options.simlibfile.replace('.simlib','_%s.simlib'%options.inputfile.split('/')[-1].split('.')[0]),
								 surveyarea,simperfect=options.perfect,batchtmpl=options.batchtmpl)
	
	if options.sim:
		#os.system('sim_SNmix.pl %s'%options.inputfile.replace('.','_MASTER.'))
		os.system('rm -r SIMLOGS_%s'%genversion)
		simtext = os.popen('sim_SNmix.pl %s'%options.inputfile.replace('.','_MASTER.')).read()
		import pdb; pdb.set_trace()
		
		# check for job completion
		job_complete=False
		while not job_complete:
			time.sleep(30)
			
		from sim_serializer import serialize
		from sim_serializer.validutils.io import save_compressed
		serialize.main(genversion,verbose=True)
		serialize.main('%s_YOUNG'%genversion,verbose=True)		
		fulldatadict = {}
		for versionsuffix in ['PLASTICC_MODEL67_SNIa-91bg','PLASTICC_MODEL52_SNIax','PLASTICC_MODEL95_SLSN-I',
							  'PLASTICC_MODEL99_ILOT','PLASTICC_MODEL99_CART','PLASTICC_MODEL42_SNIIn',
							  'PLASTICC_MODEL62_SNIbc-Templates','PLASTICC_MODEL62_SNIbc-Templates',
							  'PLASTICC_MODEL42_SNII-NMF','PLASTICC_MODEL42_SNII-Templates',
							  'PLASTICC_MODEL62_SNIbc-Templates','PLASTICC_MODEL90_SNIa-SALT2']:
			datadict = serialize.main('%s_%s'%(genversion,plasticcsuffix),verbose=True)
			for k in datadict.keys():
				fulldatadict[k] = datadict[k]
		save_compressed(table, '%s_PLASTICC.pkl.gz')
		
	if options.fit:
		fittext = os.popen('split_and_fit.pl YSE.nml').read()
		
	if options.analyze:
		pass

