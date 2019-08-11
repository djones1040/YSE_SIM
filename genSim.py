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
import time
from sim_serializer import serialize
from sim_serializer.validutils.io import save_compressed

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
			'--exptime', default=15,type="float",
			help='exposure time in seconds (default=%default)')
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

	def skynoisefrommaglim(self,maglim,zpt,areascale=3):

		skysig = 1/np.sqrt(np.pi)/areascale*np.sqrt(0.25*((0.4*10**(0.4*(zpt-maglim)) + 5)**2. - 25))

		return skysig
	
	def mksimlib(self,gcadence,rcadence,icadence,zcadence,
				 goffset,roffset,ioffset,zoffset,
				 simlibfile,simperfect=False,ztfsim=True,
				 ztf_offset=None,onedayfrac=0,exptime=15):
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
			mjd = np.arange(58240,58440,1)
			#mjd = np.arange(58240,59000,1)

		# mag lims with moon phase from Sofie
		limmjdg,limmagg,limmjdr,limmagr,limmjdi,limmagi = getmjdmaglims()
			
		# nominal survey
		count = 0; nightcount = -1; usednightcount = 0; ps1count = 0
		simliblines = []
		for m in mjd:
			nightcount += 1

			maglimg = np.interp(m,limmjdg,limmagg)
			maglimr = np.interp(m,limmjdr,limmagr)
			maglimi = np.interp(m,limmjdi,limmagi)

			
			if ztf_offset < 0:
				randval = random.uniform(0, 1)
				if ztfsim and (randval > palomardict[mjd_to_month(m)] or simperfect):
					if not nightcount % 3 or simperfect:
						iLine = random.sample(range(len(glines)),1)[0]
						for lines in [rlines]:
							line = lines[iLine][0]
							linemjd = line.split()[1]
							lineid = line.split()[2]
							skysig = line.split()[6]
							linezpt = line.split()[10]
							newskysigg = self.skynoisefrommaglim(maglimg,float(linezpt),areascale=3)
							newskysigr = self.skynoisefrommaglim(maglimr,float(linezpt),areascale=3)
							lineparts = line.split()
							lineparts[1] = '%.2f'%(m+ztf_offset); lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysigg; lineparts[3] = 'X'
							simlibline = " ".join(lineparts)
							simliblines += [simlibline]
							lineparts[1] = '%.2f'%(m+ztf_offset); lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysigr; lineparts[3] = 'Y'
							simlibline = " ".join(lineparts)
							simliblines += [simlibline]

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
						lineid = int(line.split()[2])
						skysig = line.split()[6]
						linezpt = line.split()[10]
						newskysig = self.skynoisefrommaglim(maglimg,float(linezpt),areascale=3)
						lineparts = line.split()
						lineparts[1] = '%.2f'%m; lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysig
						simlibline = " ".join(lineparts)
						simliblines += [simlibline]
						
						count += 1
						ps1count += 1
				if (not (nightcount - roffset) % rcadence and nightcount - roffset >= 0) or simperfect:
					for lines in [rlines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						skysig = line.split()[6]
						linezpt = line.split()[10]
						newskysig = self.skynoisefrommaglim(maglimr,float(linezpt),areascale=2.1)
						lineparts = line.split()
						lineparts[1] = '%.2f'%m; lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysig
						simlibline = " ".join(lineparts)
						simliblines += [simlibline]
						
						count += 1
						ps1count += 1
				if (not (nightcount - ioffset) % icadence and nightcount - ioffset >= 0) or simperfect:
					for lines in [ilines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						skysig = line.split()[6]
						linezpt = line.split()[10]
						newskysig = self.skynoisefrommaglim(maglimi,float(linezpt),areascale=2.1)
						lineparts = line.split()
						lineparts[1] = '%.2f'%m; lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysig
						simlibline = " ".join(lineparts)
						simliblines += [simlibline]
						line = lines[iLine][0]
						
						count += 1
						ps1count += 1
				if (not (nightcount - zoffset) % zcadence and nightcount - zoffset >= 0) or simperfect:
					for lines in [zlines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						skysig = line.split()[6]
						linezpt = line.split()[10]
						lineparts = line.split()
						lineparts[1] = '%.2f'%m; lineparts[2] = '%i'%count
						simlibline = " ".join(lineparts)
						simliblines += [simlibline]

						count += 1
						ps1count += 1
						
				if not nightcount % gcadence - goffset or not nightcount % rcadence - roffset or not \
				   nightcount % icadence - ioffset or not nightcount % zcadence - zoffset:
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
							skysig = line.split()[6]
							linezpt = line.split()[10]
							newskysigg = self.skynoisefrommaglim(maglimg,float(linezpt),areascale=3)
							newskysigr = self.skynoisefrommaglim(maglimr,float(linezpt),areascale=3)
							lineparts = line.split()
							lineparts[1] = '%.2f'%(m+ztf_offset); lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysigg; lineparts[3] = 'X'
							simlibline = " ".join(lineparts)
							simliblines += [simlibline]
							lineparts[1] = '%.2f'%(m+ztf_offset); lineparts[2] = '%i'%(count+1); lineparts[6] = '%.2f'%newskysigr; lineparts[3] = 'Y'
							simlibline = " ".join(lineparts)
							simliblines += [simlibline]
							
							count += 2


		fout = open(simlibfile,'w')
		print(simlibheader,file=fout)
		fout.close()

		surveyarea = 1.6*3600./(exptime+12)*0.76*7/(ps1count/float(usednightcount))*(np.pi/180.)**2.*(1-onedayfrac)
		nfields = int(surveyarea/(7*(np.pi/180.)**2.))

		fout = open(simlibfile,'a')
		for i in range(nfields):
			print("""LIBID:		  %i		# cadence from 2016W^@
RA:	   37.665487		DECL:	42.235969	  MWEBV:  0.059
NOBS:	%i		PIXSIZE:  0.500		REDSHIFT: 0.01925	  PEAKMJD: 57417.953
SUBSURVEY: PS1MD   FIELD: NORMAL"""%(i+1,count),file=fout)
			for l in simliblines:
				print(l,file=fout)
			print("""#  MJD	   IDUM	 BAND  GAIN RDNOISE	 SKYSIG	   PSF1 PSF2 PSFRAT	   ZP	ZPERR
END_LIBID:		1

		""",file=fout)
		fout.close()
		
		# one-day survey
		count = 0; nightcount = -1; usednightcount = 0; ps1count = 0
		simliblines_oneday = []
		for m in mjd:
			nightcount += 1

			maglimg = np.interp(m,limmjdg,limmagg)
			maglimr = np.interp(m,limmjdr,limmagr)
			maglimi = np.interp(m,limmjdi,limmagi)
			
			if ztf_offset < 0:
				randval = random.uniform(0, 1)
				if ztfsim and (randval > palomardict[mjd_to_month(m)] or simperfect):
					if not nightcount % 3 or simperfect:
						iLine = random.sample(range(len(glines)),1)[0]
						for lines in [rlines]:
							line = lines[iLine][0]
							linemjd = line.split()[1]
							lineid = line.split()[2]
							skysig = line.split()[6]
							linezpt = line.split()[10]
							newskysigg = self.skynoisefrommaglim(maglimg,float(linezpt),areascale=3)
							newskysigr = self.skynoisefrommaglim(maglimr,float(linezpt),areascale=3)
							lineparts = line.split()
							lineparts[1] = '%.2f'%(m+ztf_offset); lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysigg; lineparts[3] = 'X'
							simlibline = " ".join(lineparts)
							simliblines_oneday += [simlibline]
							lineparts[1] = '%.2f'%(m+ztf_offset); lineparts[2] = '%i'%(count+1); lineparts[6] = '%.2f'%newskysigr; lineparts[3] = 'Y'
							simlibline = " ".join(lineparts)
							simliblines_oneday += [simlibline]

							count += 2
							
			# replace random weather with PS1 MD observations
			# includes not much data at bright time
			randval = random.uniform(0, 1)
			if randval > haleakaladict[mjd_to_month(m)] or simperfect:
				iLine = random.sample(range(len(glines)),1)[0]
				if (not (nightcount - goffset/3)%(gcadence/3) and nightcount - goffset/3 >= 0) or simperfect:
					for lines in [glines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = int(line.split()[2])
						skysig = line.split()[6]
						linezpt = line.split()[10]
						newskysig = self.skynoisefrommaglim(maglimg,float(linezpt),areascale=3)
						lineparts = line.split()
						lineparts[1] = '%.2f'%m; lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysig
						simlibline = " ".join(lineparts)
						simliblines_oneday += [simlibline]
						
						count += 1
						ps1count += 1
				if (not (nightcount - roffset/3)%(rcadence/3) and nightcount - roffset/3 >= 0) or simperfect:
					for lines in [rlines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						skysig = line.split()[6]
						linezpt = line.split()[10]
						newskysig = self.skynoisefrommaglim(maglimr,float(linezpt),areascale=2.1)
						lineparts = line.split()
						lineparts[1] = '%.2f'%m; lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysig
						simlibline = " ".join(lineparts)
						simliblines_oneday += [simlibline]
						
						count += 1
						ps1count += 1
				if (not (nightcount - ioffset/3)%(icadence/3) and nightcount - ioffset/3 >= 0) or simperfect:
					for lines in [ilines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						skysig = line.split()[6]
						linezpt = line.split()[10]
						newskysig = self.skynoisefrommaglim(maglimi,float(linezpt),areascale=2.1)
						lineparts = line.split()
						lineparts[1] = '%.2f'%m; lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysig
						simlibline = " ".join(lineparts)
						simliblines_oneday += [simlibline]
						
						count += 1
						ps1count += 1
				if (not (nightcount - zoffset/3)%(zcadence/3) and nightcount - zoffset/3 >= 0) or simperfect:
					for lines in [zlines]:
						line = lines[iLine][0]
						linemjd = line.split()[1]
						lineid = line.split()[2]
						skysig = line.split()[6]
						linezpt = line.split()[10]
						lineparts = line.split()
						lineparts[1] = '%.2f'%m; lineparts[2] = '%i'%count
						simlibline = " ".join(lineparts)
						simliblines_oneday += [simlibline]

						count += 1
						ps1count += 1

				if not nightcount % gcadence - goffset or not nightcount % rcadence - roffset or not \
				   nightcount % icadence - ioffset or not nightcount % zcadence - zoffset:
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
							skysig = line.split()[6]
							linezpt = line.split()[10]
							newskysigg = self.skynoisefrommaglim(maglimg,float(linezpt),areascale=3)
							newskysigr = self.skynoisefrommaglim(maglimr,float(linezpt),areascale=3)
							lineparts = line.split()
							lineparts[1] = '%.2f'%(m+ztf_offset); lineparts[2] = '%i'%count; lineparts[6] = '%.2f'%newskysigg; lineparts[3] = 'X'
							simlibline = " ".join(lineparts)
							simliblines_oneday += [simlibline]
							lineparts[1] = '%.2f'%(m+ztf_offset); lineparts[2] = '%i'%(count+1); lineparts[6] = '%.2f'%newskysigr; lineparts[3] = 'Y'
							simlibline = " ".join(lineparts)
							simliblines_oneday += [simlibline]
						
							count += 2
							

		surveyarea_oneday = 1.6*3600./(exptime+12)*0.76*7/(ps1count/float(usednightcount))*(np.pi/180.)**2.*onedayfrac
		nfields_oneday = surveyarea_oneday/(7*(np.pi/180.)**2.)

		fout = open(simlibfile,'a')
		for i in range(int(nfields_oneday)):
			print("""LIBID:		  %i		# cadence from 2016W^@
RA:	   37.665487		DECL:	42.235969	  MWEBV:  0.059
NOBS:	%i		PIXSIZE:  0.500		REDSHIFT: 0.01925	  PEAKMJD: 57417.953
SUBSURVEY: PS1MD   FIELD: ONEDAY"""%(i+nfields,count),file=fout)			
			for l in simliblines_oneday:
				print(l,file=fout)
			print("""#  MJD	   IDUM	 BAND  GAIN RDNOISE	 SKYSIG	   PSF1 PSF2 PSFRAT	   ZP	ZPERR
END_LIBID:		2
""",file=fout)
		print(simlibfooter,file=fout)
		fout.close()

		#surveyarea = 2.*3600./20*7/(ps1count/float(usednightcount))*(np.pi/180.)**2.

		# 12s overhead, 16% of a telescope, 0.76 detector area, 7 sq deg
		
		return surveyarea/(1-onedayfrac),surveyarea/(1-onedayfrac)
	
		
	def mkinput(self,gcadence,rcadence,icadence,zcadence,inputfile,simlibfile,
				surveyarea,surveyarea_oneday,simperfect=False,batchtmpl=None,exptime=15):

		expadj = exptime/15.
		
		filtstr = ''
		for filt,cadence in zip('grizXY',[gcadence,rcadence,icadence,zcadence,3,3]):
			if cadence: filtstr += filt

		fout = open(inputfile.replace('.','_ia.'),'w')
		print(iainputheader%(
			inputfile.split('/')[-1].split('.')[0]+'_ia',simlibfile,
			expadj,expadj,expadj,expadj,filtstr,surveyarea),file=fout)
		if simperfect:
			print('GENPERFECT: 1',file=fout)
		fout.close()

		fout = open(inputfile.replace('.','_nonia.'),'w')
		print(noniainputheader%(
			inputfile.split('/')[-1].split('.')[0]+'_nonia',simlibfile,
			expadj,expadj,expadj,expadj,
			filtstr,surveyarea),file=fout)
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
			ratestr = 'GENOPT(NON1A): DNDZ_PEC1A POWERLAW 1E-5 2.15\nGENOPT(NON1A): DNDZ POWERLAW 1E-4	 4.5'
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

def getmjdmaglims():

	fin = open('weather/predictions_g.txt','r')
	lines = fin.readlines()
	targetgmjd = np.array(lines[0].replace('\n','').split()).astype(float)
	targetglim = np.array(lines[1].replace('\n','').split()).astype(float)
	fin.close()

	fin = open('weather/predictions_r.txt','r')
	lines = fin.readlines()
	targetrmjd = np.array(lines[0].replace('\n','').split()).astype(float)
	targetrlim = np.array(lines[1].replace('\n','').split()).astype(float)
	fin.close()

	fin = open('weather/predictions_i.txt','r')
	lines = fin.readlines()
	targetimjd = np.array(lines[0].replace('\n','').split()).astype(float)
	targetilim = np.array(lines[1].replace('\n','').split()).astype(float)
	fin.close()

	return targetgmjd,targetglim,targetrmjd,targetrlim,targetimjd,targetilim

iainputheader = """
GENVERSION: %s		   # simname
GENSOURCE:	RANDOM	 
GENMODEL:	SALT2.Guy10_UV2IR
GENPREEFIX: YSE_IA

SIMLIB_FILE: %s # simlib file

CIDOFF: 500
KCOR_FILE:	$PS1_ROOT/kcor/ZTF/kcor_PS1_ZTF_none.fits
APPLY_SEARCHEFF_OPT: 0

EXPOSURE_TIME_FILTER: g %.1f
EXPOSURE_TIME_FILTER: r %.1f
EXPOSURE_TIME_FILTER: i %.1f
EXPOSURE_TIME_FILTER: z %.1f

GENMAG_SMEAR_MODELNAME: G10
# selection criteria for generation
GENFILTERS:		  %s

GENSIGMA_SEARCH_PEAKMJD:  1.0		  # sigma-smearing for	SEARCH_PEAKMJD (days)

GENRANGE_PEAKMJD:  58240  59617
GENRANGE_RA: -999 +999
GENRANGE_DECL: -999 +999
SOLID_ANGLE: %.3f # 0.148 # 1 field, 7 sq degreees *7
# baseline for 4 filters should be 630 degrees (0.192 steradians)

SEARCHEFF_PIPELINE_FILE:  SEARCHEFF_PIPELINE_YSE.DAT
SEARCHEFF_PIPELINE_LOGIC_FILE:	SEARCHEFF_PIPELINE_LOGIC_YSE.DAT

GENRANGE_REDSHIFT:	0.001	 0.5
GENSIGMA_REDSHIFT:	0.000001
GENRANGE_TREST:	  -100.0	80.0	 # rest epoch relative to peak (days)

GENMEAN_RV:			3.1				  # mean RV to generate

OPT_MWEBV: 1

RANSEED: 128473		  # random number seed

# smear flags: 0=off, 1=on
SMEARFLAG_FLUX:	   1  # photo-stat smearing of signal, sky, etc ...
SMEARFLAG_ZEROPT:  1  # smear zero-point with zptsig

# SEARCHEFF_SPEC_FILE:	SEARCHEFF_SPEC_MF_final.DAT #this is the big one for now.
APPLY_CUTWIN_OPT:	  1
CUTWIN_NEPOCH:	 5 -5.				# require 5 epochs (no S/N requirement)
CUTWIN_TRESTMIN: -20  10
CUTWIN_TRESTMAX:   9  40
CUTWIN_MWEBV:	   0 .20

FORMAT_MASK:  2 # terse format
CUTWIN_SNRMAX:	 5.0 grizXY 2 -20. 80.	# require 1 of griz with S/N > 5

GENMEAN_SALT2x1:	 0.703
GENRANGE_SALT2x1:	-5.0  +4.0	   # x1 (stretch) range
GENSIGMA_SALT2x1:	 2.15  0.472	  # bifurcated sigmas

GENMEAN_SALT2c:		-0.04
GENRANGE_SALT2c:   -0.4	  0.4	  # color range
GENSIGMA_SALT2c:	0.033	0.125	  # bifurcated sigmas

# SALT2 alpha and beta

GENMEAN_SALT2ALPHA:	  0.14
GENMEAN_SALT2BETA:	 3.1

# cosmological params for lightcurve generation and redshift distribution
OMEGA_MATTER:  0.3
OMEGA_LAMBDA:  0.7
W0_LAMBDA:	  -1.00
H0:			  70.0	 

SIMGEN_DUMP:  8	 CID  Z	 PEAKMJD SNRMAX MAGT0_r MAGT0_g MJD_TRIGGER NON1A_INDEX
"""

noniainputheader = """
GENVERSION: %s		   # simname
GENSOURCE:	RANDOM	 
GENMODEL:	NON1A
GENPREEFIX: YSE_IA

SIMLIB_FILE: %s # simlib file

CIDOFF: 500
KCOR_FILE:	$PS1_ROOT/kcor/ZTF/kcor_PS1_ZTF_none.fits
APPLY_SEARCHEFF_OPT: 0

EXPOSURE_TIME_FILTER: g %.1f
EXPOSURE_TIME_FILTER: r %.1f
EXPOSURE_TIME_FILTER: i %.1f
EXPOSURE_TIME_FILTER: z %.1f

SEARCHEFF_PIPELINE_FILE:  SEARCHEFF_PIPELINE_YSE.DAT
SEARCHEFF_PIPELINE_LOGIC_FILE:	SEARCHEFF_PIPELINE_LOGIC_YSE.DAT

GENMAG_SMEAR_MODELNAME: G10
# selection criteria for generation
GENFILTERS:		  %s

GENSIGMA_SEARCH_PEAKMJD:  1.0		  # sigma-smearing for	SEARCH_PEAKMJD (days)

GENRANGE_PEAKMJD:  58240  59617
GENRANGE_RA: -999 +999
GENRANGE_DECL: -999 +999
SOLID_ANGLE: %.3f # 0.148 # 1 field, 7 sq degreees *7
# baseline for 4 filters should be 630 degrees (0.192 steradians)

GENRANGE_REDSHIFT:	0.001	 0.5
GENSIGMA_REDSHIFT:	0.000001
GENRANGE_TREST:	  -100.0	80.0	 # rest epoch relative to peak (days)

GENMEAN_RV:			3.1				  # mean RV to generate

OPT_MWEBV: 1

RANSEED: 128473		  # random number seed

# smear flags: 0=off, 1=on
SMEARFLAG_FLUX:	   1  # photo-stat smearing of signal, sky, etc ...
SMEARFLAG_ZEROPT:  1  # smear zero-point with zptsig

# SEARCHEFF_SPEC_FILE:	SEARCHEFF_SPEC_MF_final.DAT #this is the big one for now.
APPLY_CUTWIN_OPT:	  1
CUTWIN_NEPOCH:	 5 -5.				# require 5 epochs (no S/N requirement)
CUTWIN_TRESTMIN: -20  10
CUTWIN_TRESTMAX:   9  40
CUTWIN_MWEBV:	   0 .20

FORMAT_MASK:  2 # terse format
CUTWIN_SNRMAX:	 5.0 grizXY 2 -20. 80.	# require 1 of griz with S/N > 5

GENMEAN_SALT2x1:	 0.703
GENRANGE_SALT2x1:	-5.0  +4.0	   # x1 (stretch) range
GENSIGMA_SALT2x1:	 2.15  0.472	  # bifurcated sigmas
#GENSIGMA_SALT2x1:	  0.1  0.1		# bifurcated sigmas

GENMEAN_SALT2c:		-0.04
GENRANGE_SALT2c:   -0.4	  0.4	  # color range
GENSIGMA_SALT2c:	0.033	0.125	  # bifurcated sigmas

# SALT2 alpha and beta

GENMEAN_SALT2ALPHA:	  0.14
GENMEAN_SALT2BETA:	 3.1

# cosmological params for lightcurve generation and redshift distribution
OMEGA_MATTER:  0.3
OMEGA_LAMBDA:  0.7
W0_LAMBDA:	  -1.00
H0:			  70.0	 

SIMGEN_DUMP:  8	 CID  Z	 PEAKMJD SNRMAX MAGT0_r MAGT0_g MJD_TRIGGER NON1A_INDEX

INPUT_FILE_INCLUDE: LFs/SIMGEN_INCLUDE_NON1A_J17-beforeAdjust.INPUT
"""

noniainputheader_young = """
GENVERSION: %s		   # simname
GENSOURCE:	RANDOM	 
GENMODEL:	NON1A
GENPREEFIX: YSE_IA

SIMLIB_FILE: %s # simlib file

CIDOFF: 500
KCOR_FILE:	$PS1_ROOT/kcor/ZTF/kcor_PS1_ZTF_none.fits
APPLY_SEARCHEFF_OPT: 0

SEARCHEFF_PIPELINE_FILE:  SEARCHEFF_PIPELINE_YSE.DAT
SEARCHEFF_PIPELINE_LOGIC_FILE:	SEARCHEFF_PIPELINE_LOGIC_YSE.DAT

EXPOSURE_TIME_FILTER: g %.1f
EXPOSURE_TIME_FILTER: r %.1f
EXPOSURE_TIME_FILTER: i %.1f
EXPOSURE_TIME_FILTER: z %.1f

GENMAG_SMEAR_MODELNAME: G10
# selection criteria for generation
GENFILTERS:		  %s

GENSIGMA_SEARCH_PEAKMJD:  1.0		  # sigma-smearing for	SEARCH_PEAKMJD (days)

GENRANGE_PEAKMJD:  58240  59617
GENRANGE_RA: -999 +999
GENRANGE_DECL: -999 +999
SOLID_ANGLE: %.3f # 0.148 # 1 field, 7 sq degreees *7
# baseline for 4 filters should be 630 degrees (0.192 steradians)

GENRANGE_REDSHIFT:	0.001	 0.5
GENSIGMA_REDSHIFT:	0.000001
GENRANGE_TREST:	  -100.0	80.0	 # rest epoch relative to peak (days)

GENMEAN_RV:			3.1				  # mean RV to generate

OPT_MWEBV: 1

RANSEED: 128473		  # random number seed

# smear flags: 0=off, 1=on
SMEARFLAG_FLUX:	   1  # photo-stat smearing of signal, sky, etc ...
SMEARFLAG_ZEROPT:  1  # smear zero-point with zptsig

# SEARCHEFF_SPEC_FILE:	SEARCHEFF_SPEC_MF_final.DAT #this is the big one for now.
APPLY_CUTWIN_OPT:	  1
CUTWIN_NEPOCH:	 5 -5.				# require 5 epochs (no S/N requirement)
CUTWIN_TRESTMIN: -20  10
CUTWIN_TRESTMAX:   9  40
CUTWIN_MWEBV:	   0 .20

FORMAT_MASK:  2 # terse format
CUTWIN_SNRMAX:	 5.0 grizXY 2 -20. 80.	# require 1 of griz with S/N > 5

GENMEAN_SALT2x1:	 0.703
GENRANGE_SALT2x1:	-5.0  +4.0	   # x1 (stretch) range
GENSIGMA_SALT2x1:	 2.15  0.472	  # bifurcated sigmas
#GENSIGMA_SALT2x1:	  0.1  0.1		# bifurcated sigmas

GENMEAN_SALT2c:		-0.04
GENRANGE_SALT2c:   -0.4	  0.4	  # color range
GENSIGMA_SALT2c:	0.033	0.125	  # bifurcated sigmas

# SALT2 alpha and beta

GENMEAN_SALT2ALPHA:	  0.14
GENMEAN_SALT2BETA:	 3.1

# cosmological params for lightcurve generation and redshift distribution
OMEGA_MATTER:  0.3
OMEGA_LAMBDA:  0.7
W0_LAMBDA:	  -1.00
H0:			  70.0	 

SIMGEN_DUMP:  8	 CID  Z	 PEAKMJD SNRMAX MAGT0_r MAGT0_g MJD_TRIGGER NON1A_INDEX

PATH_NON1ASED: /project/rkessler/SURVEYS/YSE/USERS/djones/YSE_SIM/LFs/NON1A
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

FLUXERR_COR: XYgriz -0.90 2.6708 1.9649 2.5373 1.8667 1.9486 1.0957
FLUXERR_COR: XYgriz -0.48 3.4204 2.2901 3.2494 2.1756 2.3005 1.2364
FLUXERR_COR: XYgriz -0.06 3.7884 2.5821 3.5990 2.4530 2.5564 1.2882
FLUXERR_COR: XYgriz  0.37 3.3087 3.0493 3.1433 2.8968 2.8876 1.3509
FLUXERR_COR: XYgriz  0.79 2.6774 2.3056 2.5435 2.1903 2.3163 1.3317
FLUXERR_COR: XYgriz  1.21 2.0638 1.8177 1.9606 1.7268 1.8920 1.3963
FLUXERR_COR: XYgriz  1.63 1.5051 1.5201 1.4298 1.4441 1.4120 1.3013
FLUXERR_COR: XYgriz  2.06 1.0000 1.0000 0.9500 0.9500 1.0500 1.0000
FLUXERR_COR: XYgriz  2.48 1.0000 1.0000 0.9500 0.9500 1.0500 1.0000

BEGIN LIBGEN

#--------------------------------------------
"""

simlibfooter = """

END_OF_SIMLIB:
"""

mastertmpl = """
BATCH_INFO:	 sbatch	 %s 20

# nominal generation

GENVERSION: %s
GENOPT(1A): DNDZ POWERLAW  2.6E-5  2.2
GENOPT(NON1A): DNDZ POWERLAW 5E-5	4.5
GENOPT(NON1A): DNDZ_PEC1A POWERLAW	2.6E-5	2.2
GENOPT(NON1A): PATH_NON1ASED /project/rkessler/SURVEYS/YSE/USERS/djones/YSE_SIM/LFs/NON1A
# new rates, allowing SN Iax to be 31%%

GENVERSION: %s_YOUNG
GENOPT(1A): DNDZ POWERLAW  2.6E-5  2.2
GENOPT(NON1A): DNDZ POWERLAW 5E-5	4.5
GENOPT(NON1A): DNDZ_PEC1A POWERLAW	2.6E-5	2.2
GENOPT(NON1A): PATH_NON1ASED /project/rkessler/SURVEYS/YSE/USERS/djones/YSE_SIM/LFs/NON1A
GENOPT(NON1A): INPUT_FILE_INCLUDE LFs/SIMGEN_INCLUDE_NON1A_YOUNGSN.INPUT
%s
%s

%s # PLASTICC Models

ENDLIST_GENVERSION:

NGEN_UNIT:	0.014  SEASONS
# 0.286 seasons (1 year)/20 jobs

# specify sim-input files for snlc_sim.exe
SIMGEN_INFILE_Ia: %s
SIMGEN_INFILE_NONIa: %s

# define required global items to ensure uniformity among all jobs
H0: 70
ZRANGE:		 0.002	%.1f
GENPREFIX:	 %s			 # prefix of all data filenames
FORMAT_MASK: 48			  # 2=TERSE	   16=RanCID  32=FITS-FORMAT
RESET_CIDOFF: 2
PATH_SNDATA_SIM:  $SCRATCH_SIMDIR

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
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_SNIa-91bg.INPUT
GENOPT: GENTYPE 67
GENOPT: SEARCHEFF_SPEC_SCALE 1.0


# SNIax from Saurabh
GENVERSION: %s_PLASTICC_MODEL52_SNIax
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_SNIax.INPUT
GENOPT: GENTYPE 52
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# Superluminous SN:	 SLSN-I
GENVERSION:	 %s_PLASTICC_MODEL95_SLSN-I
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_SLSN-I-MOSFIT.INPUT
GENOPT: GENTYPE 95
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# pair instability SN: PISN
#GENVERSION:  %s_PLASTICC_MODEL99_PISN
#GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_PISN-MOSFIT.INPUT
#GENOPT: GENTYPE 99
#GENOPT: SEARCHEFF_SPEC_FILE ZERO

# Intermediate Luminosity Optical Transients (ILOT)
GENVERSION:	 %s_PLASTICC_MODEL99_ILOT
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_ILOT-MOSFIT.INPUT
GENOPT: GENTYPE 99
GENOPT: SEARCHEFF_SPEC_FILE ZERO

# Ca Rich Transients (CART)
GENVERSION:	 %s_PLASTICC_MODEL99_CART
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_CART-MOSFIT.INPUT
GENOPT: GENTYPE 99
GENOPT: SEARCHEFF_SPEC_FILE ZERO

# TDE
#GENVERSION:  %s_PLASTICC_MODEL15_TDE
#GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_TDE-MOSFIT.INPUT
#GENOPT: GENTYPE 15
##GENOPT: SEARCHEFF_SPEC_FILE ZERO
#GENOPT: SEARCHEFF_SPEC_SCALE 1.00

# MOSFIT-IIn
GENVERSION:	 %s_PLASTICC_MODEL42_SNIIn
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_SNIIn-MOSFIT.INPUT
GENOPT: GENTYPE 42	SIMLIB_NREPEAT 2
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# MOSFIT-Ibc
GENVERSION:	 %s_PLASTICC_MODEL62_SNIbc-MOSFIT
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_SNIbc-MOSFIT.INPUT
GENOPT: GENTYPE 62	 SIMLIB_NREPEAT 4
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# Core collapse Type II using pca (5->12 on May 9 2018)
# for end-of-challenge model release
GENVERSION: %s_PLASTICC_MODEL42_SNII-NMF
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_SNII-NMF.INPUT
GENOPT: GENTYPE 42	SIMLIB_NREPEAT 4
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# - - - - - - - - - - - - - - - - - - - - - - - - -
# legacy NON1ASED
# Core collapse Type II from K10 templates
GENVERSION: %s_PLASTICC_MODEL42_SNII-Templates
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_SNII-Templates.INPUT
GENOPT: GENTYPE 42	 SIMLIB_NREPEAT 4
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# NON1ASED-Ibc
# Core collapse Type Ibc from K10 templates
GENVERSION: %s_PLASTICC_MODEL62_SNIbc-Templates
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_SNIbc-Templates.INPUT
GENOPT: GENTYPE 62	 SIMLIB_NREPEAT 4
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

# - - - - - - - -
# Type Ia SN
GENVERSION: %s_PLASTICC_MODEL90_SNIa-SALT2
GENOPT: INPUT_FILE_INCLUDE $LSST_USERS/djones/PLASTICC/SIMGEN/SIMGEN_INCLUDE_SNIa-SALT2.INPUT
GENOPT: GENTYPE 90
GENOPT: SEARCHEFF_SPEC_SCALE 1.0

"""

if __name__ == "__main__":

	import os
	import optparse

	mks = mkSimlibs()

	usagestring = 'getSim.py <options>'
	parser = mks.add_options(usage=usagestring)
	options,  args = parser.parse_args()

	if not options.sim or not options.fit:
		surveyarea,surveyarea_oneday = mks.mksimlib(options.gcadence,options.rcadence,options.icadence,
													options.zcadence,options.goffset,options.roffset,
													options.ioffset,options.zoffset,
													options.simlibfile.replace('.simlib','_%s.simlib'%options.inputfile.split('/')[-1].split('.')[0]),
													simperfect=options.perfect,
													ztf_offset=options.ztfoffset,
													onedayfrac=options.onedayfrac,exptime=options.exptime)
		genversion = mks.mkinput(options.gcadence,options.rcadence,options.icadence,options.zcadence,
								 options.inputfile,options.simlibfile.replace('.simlib','_%s.simlib'%options.inputfile.split('/')[-1].split('.')[0]),
								 surveyarea,surveyarea_oneday,simperfect=options.perfect,batchtmpl=options.batchtmpl,exptime=options.exptime)
	
	if options.sim:
		#os.system('rm -r SIMLOGS_%s'%genversion)
		#os.system('sim_SNmix.pl %s'%options.inputfile.replace('.','_MASTER.'))

		# check for job completion
		print('waiting for job to finish...')
		job_complete=False
		while not job_complete:
			time.sleep(15)
			simtext = os.popen('squeue --user=djones1741 --format="%.30j"').read()

			job_complete = True
			for line in simtext.split('\n'):
				if '%s_'%genversion in line: job_complete = False
		print('starting serialization step')
		serialize.main(genversion,verbose=True,filters='grizXY',dirpath='$SCRATCH_SIMDIR')
		os.system('cp $SCRATCH_SIMDIR/%s/%s.DUMP dump/'%(genversion,genversion))

		serialize.main('%s_YOUNG'%genversion,verbose=True,filters='grizXY',dirpath='$SCRATCH_SIMDIR')
		os.system('cp $SCRATCH_SIMDIR/%s_YOUNG/%s_YOUNG.DUMP dump/'%(genversion,genversion))

		fulldatadict = {}
		for versionsuffix in ['PLASTICC_MODEL67_SNIa-91bg','PLASTICC_MODEL52_SNIax','PLASTICC_MODEL95_SLSN-I',
							  'PLASTICC_MODEL99_ILOT','PLASTICC_MODEL99_CART','PLASTICC_MODEL42_SNIIn',
							  'PLASTICC_MODEL62_SNIbc-Templates','PLASTICC_MODEL62_SNIbc-Templates',
							  'PLASTICC_MODEL42_SNII-NMF','PLASTICC_MODEL42_SNII-Templates',
							  'PLASTICC_MODEL62_SNIbc-Templates','PLASTICC_MODEL90_SNIa-SALT2']:
			try:
				datadict = serialize.main('%s_%s'%(genversion,versionsuffix),verbose=True,save=False,
										  filters='grizXY',dirpath='$SCRATCH_SIMDIR')
				for k in datadict.keys():
					fulldatadict[k] = datadict[k]
			except Exception as e:
				print('error for %s: %s'%(versionsuffix,e))
			os.system('cp $SCRATCH_SIMDIR/%s_%s/%s_%s.DUMP dump/'%(genversion,versionsuffix,genversion,versionsuffix))
			os.system('cat $SCRATCH_SIMDIR/SIM/%s_%s/%s_%s.DUMP >> dump/%s_PLASTICC.dump'%(genversion,versionsuffix,genversion,versionsuffix,genversion))

		save_compressed(fulldatadict, '%s_PLASTICC.pkl.gz'%genversion)
		
	if options.fit:
		fittext = os.popen('split_and_fit.pl YSE.nml').read()
		
	if options.analyze:
		pass

