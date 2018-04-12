#!/usr/bin/env python
from __future__ import print_function
import gzip
import pickle
import pandas as pd

# python YSFoM.py -d $SNDATA_ROOT/SIM/yse_gr_gi_gz/yse_gr_gi_gz.DUMP -y $SNDATA_ROOT/SIM/yse_gr_gi_gz_YOUNG/yse_gr_gi_gz_YOUNG.DUMP -f YSE/yse_gr_gi_gz/FITOPT000.FITRES -d $SNDATA_ROOT/SIM/yse_gr_gi_gz/yse_gr_gi_gz.DUMP
# python YSFoM.py -d $SNDATA_ROOT/SIM/yse_ztf_gr_gi_gz/yse_ztf_gr_gi_gz.DUMP -y $SNDATA_ROOT/SIM/yse_ztf_gr_gi_gz_YOUNG/yse_ztf_gr_gi_gz_YOUNG.DUMP -f YSE/yse_ztf_gr_gi_gz/FITOPT000.FITRES -d $SNDATA_ROOT/SIM/yse_ztf_gr_gi_gz/yse_ztf_gr_gi_gz.DUMP
# python YSFoM.py -d $SNDATA_ROOT/SIM/yse_ztf/yse_ztf.DUMP -y $SNDATA_ROOT/SIM/yse_ztf_YOUNG/yse_ztf_YOUNG.DUMP -f YSE/yse_ztf/FITOPT000.FITRES -d $SNDATA_ROOT/SIM/yse_ztf/yse_ztf.DUMP
# python YSFoM.py -d $SNDATA_ROOT/SIM/yse_ztf_gri/yse_ztf_gri.DUMP -y $SNDATA_ROOT/SIM/yse_ztf_gri_YOUNG/yse_ztf_gri_YOUNG.DUMP -f YSE/yse_ztf_gri/FITOPT000.FITRES -d $SNDATA_ROOT/SIM/yse_ztf_gri/yse_ztf_gri.DUMP

# days before max that explosion happens for young templates
explosiondict = {601:-12, # 1993J, IIb
				 602:-23, # 2011dh, IIb, this light curve looks weird
				 603:-30, # 2008ax, IIb
				 604:-7, # 2013am, IIP
				 605:-8, # 2013ej, IIP
				 606:-6, # 2012aw, IIP
				 607:-8, # 2015bh, IIn
				 608:-27.0, # 2011ht, IIn - ok now, no more break because I fudged the max
				 609:-14, # 2009ip, IIn
				 610:-16, # iPTF13bvn, Ib
				 611:-12, # 2007Y, Ib
				 612:-9, # 2007gr, Ic
				 613:-6, # 2013ge, Ic
				 614:-16, # 2011fe, Ia
				 615:-13, # 2014J, Ia, possible issues here
				 616:-10} # 2005hk, Ia

explosiondict = {605: -8.0, 604: -7.0, 607: -8.0, 606: -6.1,
				 601: -15.1, 603: -38.0, 602: -40.1, 612: -9.1,
				 613: -6.1, 610: -16.1, 611: -12.1, 609: -14.1,
				 608: -42.1, 614: -17.0, 617: -14.0, 618: -19.1,
				 620: -94.1, 621: -64.0, 616: -10.0, 619: -15.1,
				 0:-15,
				 216: -27.1, 217: -18.0,
				 214: -14.1, 215: -14.0,
				 212: -19.0, 213: -14.0,
				 210: -13.0, 211: -16.1,
				 218: -26.0, 219: -29.0,
				 104: -20.1, 227: -14.0,
				 103: -23.1, 234: -28.0,
				 235: -15.0, 230: -23.1,
				 231: -21.1, 232: -14.0,
				 233: -21.1, 224: -16.0,
				 21: -11.1, 22: -23.0,
				 2: -7.1, 403: -29.1,
				 402: -30.0, 401: -30.0,
				 509: -29.0, 506: -15.1,
				 220: -16.1, 502: -14.1,
				 503: -17.0, 201: -22.1,
				 203: -29.0, 202: -21.1,
				 205: -14.0, 204: -14.0,
				 207: -25.1, 206: -28.1,
				 209: -28.1, 208: -14.1,
				 400: -17.0, 229: -13.1,
				 228: -25.1, 102: -21.0,
				 226: -12.1, 225: -28.0,
				 101: -21.1, 223: -14.1,
				 222: -17.1, 221: -15.1,
				 105: -24.1, 0:-15}

explosiondict_orig = {216: -27.1, 217: -18.0,
					  214: -14.1, 215: -14.0,
					  212: -19.0, 213: -14.0,
					  210: -13.0, 211: -16.1,
					  218: -26.0, 219: -29.0,
					  104: -20.1, 227: -14.0,
					  103: -23.1, 234: -28.0,
					  235: -15.0, 230: -23.1,
					  231: -21.1, 232: -14.0,
					  233: -21.1, 224: -16.0,
					  21: -11.1, 22: -23.0,
					  2: -7.1, 403: -29.1,
					  402: -30.0, 401: -30.0,
					  509: -29.0, 506: -15.1,
					  220: -16.1, 502: -14.1,
					  503: -17.0, 201: -22.1,
					  203: -29.0, 202: -21.1,
					  205: -14.0, 204: -14.0,
					  207: -25.1, 206: -28.1,
					  209: -28.1, 208: -14.1,
					  400: -17.0, 229: -13.1,
					  228: -25.1, 102: -21.0,
					  226: -12.1, 225: -28.0,
					  101: -21.1, 223: -14.1,
					  222: -17.1, 221: -15.1,
					  105: -24.1, 0:-15}

noniadict = {0:'Ia',201:'IIP',204:'IIP',208:'IIP',
			 210:'IIP',213:'IIP',214:'IIP',
			 215:'IIP',216:'IIP',219:'IIP',
			 220:'IIP',221:'IIP',222:'IIP',
			 223:'IIP',224:'IIP',225:'IIP',
			 226:'IIP',227:'IIP',228:'IIP',
			 229:'IIP',230:'IIP',231:'IIP',
			 232:'IIP',233:'IIP',235:'IIP',
			 206:'IIn',209:'IIn',2:'IIL',
			 103:'Ib',104:'Ib',105:'Ib',
			 202:'Ib',203:'Ib',212:'Ib',
			 234:'Ib',21:'Ic',22:'Ic',
			 101:'Ic',102:'Ic',205:'Ic',
			 207:'Ic',211:'Ic',217:'Ic',
			 218:'Ic',400:'IIb',401:'IIb',
			 402:'IIb',403:'IIb',502:'Ia-91bg',
			 503:'Ia-91bg',506:'Ia-91bg',509:'Ia-91bg',
			 601:'IIb',602:'IIb',603:'IIb',
			 604:'IIP',605:'IIP',606:'IIP',
			 607:'IIn',608:'IIn',609:'IIn',
			 610:'Ib',611:'Ib',612:'Ic',
			 613:'Ic',614:'Ia',615:'Ia',
			 616:'Iax',617:'Ia',618:'Ia',
			 619:'Ia',620:'SLSN',621:'TDe'}

namedict = {0:'SALT2',
			601:'1993J',
			602:'2011dh',
			603:'2008ax',
			604:'2013am',
			605:'2013ej',
			606:'2012aw',
			607:'2015bh',
			608:'2011ht',
			609:'2009ip',
			610:'iPTF13bvn',
			611:'2007Y',
			612:'2007gr',
			613:'2013ge',
			614:'2011fe',
			615:'2014J',
			616:'2005hk',
			617:'2012fr',
			618:'2017cbv',
			619:'2018gv',
			620:'2015bn',
			621:'PS1-10jh'}


import glob
import numpy as np
import os
import pylab as plt
import snana

def perfectlcplots():
	import glob

	lcfiles = glob.glob('/usr/local/SNDATA_ROOT/SIM/youngsn_perfect/youngsn_perfect_SN*DAT')

	for l in lcfiles:
		sn = snana.SuperNova(l)
		if sn.SIM_NON1a != '621   (NONIA index)': continue
		plt.clf()
		for f in sn.FILTERS:
			plt.plot(sn.MJD[sn.FLT == f],sn.FLUXCAL[sn.FLT == f],'o',label=f)
		plt.xlabel('MJD'); plt.legend()
		plt.ylabel('Flux')
		plt.title(sn.SIM_COMMENT.split('/')[-1].split('NEW')[0])
		plt.savefig('%s_SIM.png'%sn.SIM_COMMENT.split('/')[-1].split('NEW')[0])
		break

def getExplosionDate(simdir='/usr/local/SNDATA_ROOT/SIM/youngsn_perfect/'):
	import snana
	files = glob.glob('%s/*DAT'%simdir)
	expdict = {}
	for f in files:
		sn = snana.SuperNova(f)
		if sn.SIM_NON1a.split()[0] in expdict.keys(): continue
		plt.clf()
		for f in sn.FILTERS:
			plt.plot(sn.MJD[sn.FLT == f],sn.FLUXCAL[sn.FLT == f],'o')
		expdict[sn.SIM_NON1a.split()[0]] = np.round(np.min(sn.MJD[sn.FLUXCAL > 1])-float(sn.SIM_PEAKMJD.split()[0]),decimals=1)

	print(expdict)
			
def synphot(x,spc,pb,zp,plot=False,oplot=False,allowneg=False):
	import numpy as np
	x1=pb.split('/')
	nx=len(x1)
	xfin=x1[nx-1]
	fp='/'
	for i in range(nx-1): fp=fp+x1[i]+'/'


	import pysynphot.spectrum as S
	sp = S.Vega
	mag = zp - 2.5 * np.log10( synflux(x,spc,pb,plot=plot,oplot=oplot,
										 allowneg=allowneg))
	vegamag = zp - 2.5 * np.log10( synflux(x,sp(x),pb,plot=plot,oplot=oplot,
										 allowneg=allowneg))

	return(mag)

def synflux(x,spc,pb,plot=False,oplot=False,allowneg=False):
	import numpy as np

	nx = len(x)
	pbphot = 1
	pbx,pby = np.loadtxt(pb,unpack=True)

	npbx = len(pbx)
	if (len(pby) != npbx):
		print(' pbs.wavelength and pbs.response have different sizes')

	if nx == 1 or npbx == 1:
		print('warning! 1-element array passed, returning 0')
		return(spc[0]-spc[0])

	diffx = x[1:nx]-x[0:nx-1]
	diffp = pbx[1:npbx]-pbx[0:npbx-1]

	if (np.min(diffx) <= 0) or (np.min(diffp) <= 0):
		print('passed non-increasing wavelength array')

	if x[0] > pbx[0]:
		print("spectrum doesn''t go blue enough for passband!")

	if x[nx-1] < pbx[npbx-1]:
		print("spectrum doesn''t go red enough for passband!")

	g = np.where((x >= pbx[0]) & (x <= pbx[npbx-1]))  # overlap range

	pbspl = np.interp(x[g],pbx,pby)#,kind='cubic')

	if not allowneg:
		pbspl = pbspl
		col = np.where(pbspl < 0)[0]
		if len(col): pbspl[col] = 0

	if (pbphot): pbspl *= x[g]


	res = np.trapz(pbspl*spc[g],x[g])/np.trapz(pbspl,x[g])

	return(res)
		
def getBMax():
	sedfiles = glob.glob('LFs/NON1A/youngsn/*10jh.SED')
	filtfile = '/Users/David/Dropbox/research/manglefig/json/filters/keplercam_B.dat'
	for s in sedfiles:
		if 'NEW' in s: continue
		phase,wave,flux = np.loadtxt(s,unpack=True)
		maxphase,bmags = [],[]
		for p in np.unique(phase):
			bmag = synphot(wave[phase == p],flux[phase == p],filtfile,0)
			maxphase += [p]; bmags += [bmag]

		maxphase = np.array(maxphase)[bmags == min(bmags)][0]

		print('max phase is %.2f for file %s'%(maxphase,s))
		fout = open(s.replace('.SED','NEW.SED'),'w')
		
		for p,w,f,i in zip(phase,wave,flux,range(len(phase))):
			#maxphase = -15
			if ('2015bn' in s or 'PS1-10jh' in s) and p % 2: continue
			if '2015bn' in s: f *= 1e-5
			print('  %.4f  %.2f  %8.5e'%(p-maxphase,w,f),file=fout)
		fout.close()
			
class YSFoM:
	def __init__(self):
		pass

	def add_options(self, parser=None, usage=None, config=None):
		import optparse
		if parser == None:
			parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

		parser.add_option(
			'-s','--simname', default='yse',type="string",
			help='name of sim (default=%default)')
		parser.add_option(
			'-d','--dumpfile', default='$SNDATA_ROOT/SIM/yse/yse.DUMP',type="string",
			help='DUMP file (default=%default)')
		parser.add_option(
			'-y','--youngdumpfile', default='$SNDATA_ROOT/SIM/yse_YOUNG/yse_YOUNG.DUMP',type="string",
			help='DUMP file for young SN (default=%default)')
		parser.add_option(
			'-f','--fitresfile', default='YSE/yse/FITOPT000.FITRES',type="string",
			help='name of sim (default=%default)')
		parser.add_option(
			'--histplot', default='ysehist.png',type="string",
			help='name of histogram plot (default=%default)')
		parser.add_option(
			'--lcplot', default='yselc.png',type="string",
			help='name of light curve plot (default=%default)')
		parser.add_option(
			'--yseztfplot', default='yseztf.png',type="string",
			help='name of YSE/ZTF plot (default=%default)')
		parser.add_option(
			'--pkmagcut', default=None,type="float",
			help='maximum peak mag (default=%default)')
		parser.add_option(
			'--zcut', default=None,type="float",
			help='maximum z (default=%default)')
		parser.add_option(
			'--serialize', default=False,action="store_true",
			help='convert SNANA sims to pickle files, if set (default=%default)')

		return parser

	def earlydetect(self,youngdumpfile,debug=False,pkmagcut=None,zcut=None,pklfile=None):
		import snana
		import pylab as plt

		CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
			np.loadtxt(os.path.expandvars(youngdumpfile),unpack=True,
					   skiprows=7,usecols=[1,2,3,4,5,6,7,8,9,10])
		NON1A_INDEX = NON1A_INDEX.astype(int)

		if pkmagcut:
			ipk = np.where(MAGT0_r < pkmagcut)[0]
			CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
				CID[ipk],Z[ipk],PEAKMJD[ipk],S2c[ipk],S2x1[ipk],SNRMAX[ipk],\
				MAGT0_r[ipk],MAGT0_g[ipk],MJD_TRIGGER[ipk],NON1A_INDEX[ipk]
		if zcut:
			iz = np.where(Z < zcut)[0]
			CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
				CID[iz],Z[iz],PEAKMJD[iz],S2c[iz],S2x1[iz],SNRMAX[iz],\
				MAGT0_r[iz],MAGT0_g[iz],MJD_TRIGGER[iz],NON1A_INDEX[iz]

		X, y = read_data(pklfile)
		onedaycount,twodaycount = 0,0
		onedayidx,twodayidx = [],[]
		onedayz,twodayz = [],[]
		for cd,mt,pkm,ni,z in zip(CID,MJD_TRIGGER,PEAKMJD,NON1A_INDEX,Z):
			#mt = np.min(sn.MJD[sn.FLUXCAL/sn.FLUXCALERR > 5])
			mjdmin = 1e8
			sn = X.loc[cd]
			for filt in 'grizXY':
				mjd = sn['mjd_%s'%filt]
				mjd_goodsnr = mjd[sn['fluxcal_%s'%filt]/sn['fluxcalerr_%s'%filt] > 5]
				if len(mjd_goodsnr) and np.min(mjd_goodsnr) < mjdmin: mjdmin = np.min(mjd_goodsnr)
			
			if ni != 0 and mjdmin - pkm <= explosiondict[ni]+0.5:
				onedaycount += 1
				onedayidx += [noniadict[ni]]
				onedayz += [z]
			elif debug and ni == 621 and ni != 0 and z < 0.2: #and mt - pkm <= explosiondict[ni]+1.5: #ni != 0 and mt != 1000000.0 and noniadict[ni] == 'IIn':
				#sn = snana.SuperNova(simname='yse_gr_gi_gz',simdir=youngdumpfile.split('/')[-2],snid=cd)
				plt.clf()
				plt.ion()
				print('%i days before max'%(mt-pkm))
				for f in 'griz':
					plt.errorbar(sn['mjd_%s'%f],sn['fluxcal_%s'%f],
								 yerr=sn['fluxcalerr_%s'%f],fmt='o',label=f)
				plt.title('%i, %s, %s'%(ni,noniadict[ni],namedict[ni]))
				plt.axvline(pkm,label='peak',color='k')
				plt.axvline(mt,label='disc',color='r')
				plt.legend()
				plt.show()
				import pdb; pdb.set_trace()
			if ni != 0 and mjdmin - pkm <= explosiondict[ni]+1.5:
				twodaycount += 1
				twodayidx += [noniadict[ni]]
				twodayz += [z]
		print('%i SNe discovered w/i a day of explosion'%onedaycount)
		print('%i SNe discovered w/i two days of explosion'%twodaycount)
		nonia_idx_oneday,counts_oneday = np.unique(onedayidx,return_counts=True)
		nonia_idx_twoday,counts_twoday = np.unique(twodayidx,return_counts=True)
		
		iibct,iinct,iipct,ibct,icct,iact,iaxct,slsnct,tdect = 0,0,0,0,0,0,0,0,0
		for ix,ct in zip(nonia_idx_oneday,counts_oneday):
			if ix == 'IIb': iibct += ct
			elif ix == 'IIn': iinct += ct
			elif ix == 'IIP': iipct += ct
			elif ix == 'Ib': ibct += ct
			elif ix == 'Ic': icct += ct
			elif ix == 'Ia': iact += ct
			elif ix == 'Iax': iaxct += ct
			elif ix == 'SLSN': slsnct += ct
			elif ix == 'TDe': tdect += ct
			
		iibct_2day,iinct_2day,iipct_2day,ibct_2day,icct_2day,iact_2day,iaxct_2day,slsnct_2day,tdect_2day = 0,0,0,0,0,0,0,0,0
		for ix,ct in zip(nonia_idx_twoday,counts_twoday):
			if ix == 'IIb': iibct_2day += ct
			elif ix == 'IIn': iinct_2day += ct
			elif ix == 'IIP': iipct_2day += ct
			elif ix == 'Ib': ibct_2day += ct
			elif ix == 'Ic': icct_2day += ct
			elif ix == 'Ia': iact_2day += ct
			elif ix == 'Iax': iaxct_2day += ct
			elif ix == 'SLSN': slsnct_2day += ct
			elif ix == 'TDe': tdect_2day += ct

		ax = plt.subplot(151)
		ax.hist(twodayz,label='%i SNe w/ 2 days'%twodaycount)
		ax.hist(onedayz,label='%i SNe w/ 1 day'%onedaycount)
		ax.legend()
		ax.set_xlabel('$z$',fontsize=15)
		ax.set_xlabel('# of SNe',fontsize=15)
		
		print('Detected w/i 1 day of explosion: %i IIb, %i IIn, %i IIP, %i Ib, %i Ic, %i Ia, %i Iax, %i SLSN, %i TDe'%(
			iibct,iinct,iipct,ibct,icct,iact,iaxct,slsnct,tdect))
		print('Detected w/i 2 days of explosion: %i IIb, %i IIn, %i IIP, %i Ib, %i Ic, %i Ia, %i Iax, %i SLSN, %i TDe'%(
			iibct_2day,iinct_2day,iipct_2day,ibct_2day,icct_2day,iact_2day,iaxct_2day,slsnct_2day,tdect_2day))

		labels = 'IIb','IIn','IIP','Ib','Ic','Ia','Iax','SLSN-I','TDe'
		sizes = [iibct,iinct,iipct,ibct,icct,iact,iaxct,slsnct,tdect]
		colors = ['#1BE7FF','#6EEB83','#E4FF1A','#E8AA14','#FF5714','#731963','C0','C1','C2']
 
		# Plot
		ax = plt.subplot(152)
		ax.set_title('SN types detected w/i \n1 day of explosion\n (ZTF or PS1 detection)')
		def absolute_value(val):
			a  = np.round(val/100.*np.sum(sizes), 0)
			return '%i'%a
		ax.pie(sizes, labels=labels, colors=colors,
				autopct=absolute_value, shadow=False, startangle=140)
		
		
	def sn_plus_cosmo(self,dumpfile,fitresfile,pkmagcut=None,zcut=None):
		import getmu
		from txtobj import txtobj
		fr = txtobj(fitresfile,fitresheader=True)
		fr = getmu.getmu(fr,sigint=None)
		iMuerr = np.where(fr.muerr < 0.15)[0]
		for k in fr.__dict__.keys():
			fr.__dict__[k] = fr.__dict__[k][iMuerr]
		
		CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
			np.loadtxt(os.path.expandvars(dumpfile),unpack=True,
					   skiprows=7,usecols=[1,2,3,4,5,6,7,8,9,10])
		NON1A_INDEX = NON1A_INDEX.astype(int)

		if pkmagcut and not zcut:
			ipk = np.where(MAGT0_r < pkmagcut)[0]
			CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
				CID[ipk],Z[ipk],PEAKMJD[ipk],S2c[ipk],S2x1[ipk],SNRMAX[ipk],\
				MAGT0_r[ipk],MAGT0_g[ipk],MJD_TRIGGER[ipk],NON1A_INDEX[ipk]

		if zcut and not pkmagcut:
			iz = np.where(Z < zcut)[0]
			CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
				CID[iz],Z[iz],PEAKMJD[iz],S2c[iz],S2x1[iz],SNRMAX[iz],\
				MAGT0_r[iz],MAGT0_g[iz],MJD_TRIGGER[iz],NON1A_INDEX[iz]

		if zcut and pkmagcut:
			ipkz = np.where((MAGT0_r < pkmagcut) | (Z < zcut))[0]
			CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
				CID[ipkz],Z[ipkz],PEAKMJD[ipkz],S2c[ipkz],S2x1[ipkz],SNRMAX[ipkz],\
				MAGT0_r[ipkz],MAGT0_g[ipkz],MJD_TRIGGER[ipkz],NON1A_INDEX[ipkz]
			
		if zcut or pkmagcut:
			ifr,rpk = [],[]
			for i,j in zip(fr.CID,range(len(fr.CID))):
				if i in CID:
					ifr += [j]
					rpk += [MAGT0_r[i == CID][0]]
			for k in fr.__dict__.keys():
				fr.__dict__[k] = fr.__dict__[k][ifr]
		else:
			ifr,rpk = [],[]
			for i,j in zip(fr.CID,range(len(fr.CID))):
				if i in CID:
					ifr += [j]
					rpk += [MAGT0_r[i == CID][0]]

		ax = plt.subplot(153)
		zbins = np.arange(0,0.4,0.02)
		ax.hist(Z,label='All SNe',bins=zbins)
		ax.hist(fr.zHD[fr.SIM_TYPE_INDEX == 1],label='SNe Ia, $\sigma_{\mu} < 0.15$',bins=zbins)
		ax.legend()
		ax.set_xlabel('$z$',fontsize=15)
		ax.set_xlabel('# of SNe',fontsize=15)
		print('%i SNe pass cuts'%len(fr.CID[fr.SIM_TYPE_INDEX == 1]))
		print('median redshift = %.3f'%np.median(fr.zHD))

		ax.set_title('%i SNe Ia (after cuts), $z_{med}$ = %.3f, \nN(z < 0.1) = %i'%(len(fr.CID),np.median(fr.zHD),len(fr.CID[fr.zHD < 0.1])))
		
		ax = plt.subplot(154)
		magbins = np.arange(15,22,0.5)
		ax.hist(MAGT0_r,label='All SNe',bins=magbins)
		ax.hist(np.array(rpk)[np.where(fr.SIM_TYPE_INDEX == 1)[0]],label='SNe Ia, $\sigma_{\mu} < 0.15$',bins=magbins)
		ax.legend()
		ax.set_xlabel('$r_{pk}$',fontsize=15)
		ax.set_xlabel('# of SNe',fontsize=15)
		ax.set_title('%i SNe, %i SNe Ia (after cuts)'%(len(MAGT0_r),len(fr.CID)))
		ax.set_xlim([15,22])

		count,idx,z = 0,[],[]
		for cd,mt,pkm,ni,z in zip(CID,MJD_TRIGGER,PEAKMJD,NON1A_INDEX,Z):			
			count += 1
			idx += [noniadict[ni]]
			z += [z]
		nonia_idx,counts = np.unique(idx,return_counts=True)
		
		iibct,iinct,iipct,iilct,ibct,icct,iact,iabgct,iaxct,tdect,slsnct = 0,0,0,0,0,0,0,0,0,0,0
		for ix,ct in zip(nonia_idx,counts):
			if ix == 'IIb': iibct += ct
			elif ix == 'IIn': iinct += ct
			elif ix == 'IIP': iipct += ct
			elif ix == 'Ib': ibct += ct
			elif ix == 'Ic': icct += ct
			elif ix == 'Ia': iact += ct
			elif ix == 'IIL': iilct += ct
			elif ix == 'Ia-91bg': iabgct += ct
			elif ix == 'Iax': iaxct += ct
			elif ix == 'TDe': tdect += ct
			elif ix == 'SLSN': slsnct += ct
			else:
				import pdb; pdb.set_trace()
			
		print('Detected by YSE: %i IIb, %i IIn, %i IIP, %i IIL, %i Ib, %i Ic, %i Ia %i Ia-91bg, %i Iax, %i SLSN, %i TDe'%(
			iibct,iinct,iipct,iilct,ibct,icct,iact,iabgct,iaxct, slsnct,tdect))

		labels = 'IIb','IIn','IIP','IIL','Ib','Ic','Ia', '91bg', 'Iax', 'SLSN-I', 'TDe'
		sizes = [iibct,iinct,iipct,iilct,ibct,icct,iact,iabgct,iaxct,slsnct,tdect]
		colors = ['#1BE7FF','#6EEB83','#E4FF1A','#E8AA14','#FF5714','#731963','g','r','C0','C1','C3']
 
		# Plot
		ax = plt.subplot(155)
		ax.set_title('%i SNe from YSE'%np.sum(sizes))
		def absolute_value(val):
			a  = np.round(val/100.*np.sum(sizes), 0)
			return '%i'%a

		ax.pie(sizes, labels=labels, colors=colors,
				autopct=absolute_value, shadow=False, startangle=140)

		
	def youngLCplots(self,youngdumpfile,nplotsfortype=2,pklfile=None):
		CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
			np.loadtxt(os.path.expandvars(youngdumpfile),unpack=True,
					   skiprows=7,usecols=[1,2,3,4,5,6,7,8,9,10])
		NON1A_INDEX = NON1A_INDEX.astype(int)
		X, y = read_data(pklfile)
		
		twodaycount,plotcount = 0,0
		twodayidx = []
		for cd,mt,pkm,ni,z in zip(CID,MJD_TRIGGER,PEAKMJD,NON1A_INDEX,Z):
			#mt = np.min(sn.MJD[sn.FLUXCAL/sn.FLUXCALERR > 5])

			mjdmin = 1e8
			sn = X.loc[cd]
			for filt in 'griz':
				mjd = sn['mjd_%s'%filt]
				mjd_goodsnr = mjd[sn['fluxcal_%s'%filt]/sn['fluxcalerr_%s'%filt] > 5]
				if len(mjd_goodsnr) and np.min(mjd_goodsnr) < mjdmin: mjdmin = np.min(mjd_goodsnr)
			
			if ni != 0 and mjdmin - pkm <= explosiondict[ni]+1.5:
				twodayidx += [noniadict[ni]]
				if len(np.where(np.array(twodayidx) == noniadict[ni])[0]) > nplotsfortype:
					continue
				plotcount += 1
				ax = plt.subplot(4,4,plotcount)
				sn = X.loc[cd]
				for f in 'grizXY':
					ax.errorbar(sn['mjd_%s'%f],sn['fluxcal_%s'%f],
								yerr=sn['fluxcalerr_%s'%f],fmt='o',label=f)

				ax.text(0.5,0.9,'%i, %s, %s'%(ni,noniadict[ni],namedict[ni]),
						ha='center',transform=ax.transAxes,
						bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.5})
				ax.axvline(pkm,label='peak',color='k')
				ax.axvline(mt,label='disc',color='r')
				if plotcount == 1: ax.legend()
		return

	def sampleLCplots(self,dumpfile,nplotsfortype=2,pklfile=None):

		CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
			np.loadtxt(os.path.expandvars(dumpfile),unpack=True,
					   skiprows=7,usecols=[1,2,3,4,5,6,7,8,9,10])
		NON1A_INDEX = NON1A_INDEX.astype(int)

		X, y = read_data(pklfile)
		
		twodaycount,plotcount = 0,0
		idx = []
		for cd,mt,pkm,ni,z in zip(CID,MJD_TRIGGER,PEAKMJD,NON1A_INDEX,Z):
			if mt > 999999: continue
			
			idx += [noniadict[ni]]
			if len(np.where(np.array(idx) == noniadict[ni])[0]) > nplotsfortype:
				continue
			plotcount += 1
			ax = plt.subplot(5,5,plotcount)
			sn = X.loc[cd]
			for f in 'grizXY':
				ax.errorbar(sn['mjd_%s'%f],sn['fluxcal_%s'%f],
							yerr=sn['fluxcalerr_%s'%f],fmt='o',label=f)

			ax.text(0.5,0.8,'%i, %s, \nz=%.2f, g$_{pk}$=%.1f'%(
				ni,noniadict[ni],z,-2.5*np.log10(np.max(sn['fluxcal_g']))+27.5),
					ha='center',transform=ax.transAxes,
					bbox={'facecolor':'1.0','edgecolor':'1.0','alpha':0.5})
			ax.axvline(pkm,label='peak',color='k')
			ax.axvline(mt,label='disc',color='r')
			ax.set_ylim(bottom=-250)
			if plotcount == 1: ax.legend()

	def ZTF_YSE(self,youngdumpfile,debug=False,pkmagcut=None,zcut=None,pklfile=None):
		import snana
		import pylab as plt

		CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
			np.loadtxt(os.path.expandvars(youngdumpfile),unpack=True,
					   skiprows=7,usecols=[1,2,3,4,5,6,7,8,9,10])
		NON1A_INDEX = NON1A_INDEX.astype(int)

		if pkmagcut:
			ipk = np.where(MAGT0_r < pkmagcut)[0]
			CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
				CID[ipk],Z[ipk],PEAKMJD[ipk],S2c[ipk],S2x1[ipk],SNRMAX[ipk],\
				MAGT0_r[ipk],MAGT0_g[ipk],MJD_TRIGGER[ipk],NON1A_INDEX[ipk]
		if zcut:
			iz = np.where(Z < zcut)[0]
			CID,Z,PEAKMJD,S2c,S2x1,SNRMAX,MAGT0_r,MAGT0_g,MJD_TRIGGER,NON1A_INDEX = \
				CID[iz],Z[iz],PEAKMJD[iz],S2c[iz],S2x1[iz],SNRMAX[iz],\
				MAGT0_r[iz],MAGT0_g[iz],MJD_TRIGGER[iz],NON1A_INDEX[iz]

		X, y = read_data(pklfile)

		onedaycount,twodaycount = 0,0
		onedayidx,twodayidx = [],[]
		onedayz,twodayz = [],[]
		for cd,mt,pkm,ni,z in zip(CID,MJD_TRIGGER,PEAKMJD,NON1A_INDEX,Z):
			mjdmin = 1e8
			sn = X.loc[cd]
			for filt in 'grizXY':
				mjd = sn['mjd_%s'%filt]
				mjd_goodsnr = mjd[sn['fluxcal_%s'%filt]/sn['fluxcalerr_%s'%filt] > 5]
				if len(mjd_goodsnr) and np.min(mjd_goodsnr) < mjdmin: mjdmin = np.min(mjd_goodsnr)

			#import pdb; pdb .set_trace()			
			if ni != 0 and mjdmin - pkm <= explosiondict[ni]+0.5:

				#sn = snana.SuperNova(simname='YSE',simdir=youngdumpfile.split('/')[-2],snid=cd)

				idx = np.where(X.snid == cd)[0]
				sn = X.iloc[idx]
				#yseobs = np.min(np.abs(sn.MJD[(sn.FLT != 'X') & (sn.FLT != 'Y')]-(pkm+explosiondict[ni])))
				#ztfobs = np.min(np.abs(sn.MJD[(sn.FLT == 'X') | (sn.FLT == 'Y')]-(pkm+explosiondict[ni])))
				yseobs = np.min(np.abs(np.array(sn['mjd'])[0][(np.array(sn['FLT'])[0] != 'X') &
															  (np.array(sn['FLT'])[0] != 'Y')]-(pkm+explosiondict[ni])))
				ztfobs = np.min(np.abs(np.array(sn['mjd'])[0][(np.array(sn['FLT'])[0] == 'X') |
															  (np.array(sn['FLT'])[0] == 'Y')]-(pkm+explosiondict[ni])))

				if yseobs < 0.5 and ztfobs < 0.5:
					onedaycount += 1
					onedayidx += [noniadict[ni]]
					onedayz += [z]
				
			if ni != 0 and mt - pkm <= explosiondict[ni]+1.5:
				#sn = snana.SuperNova(simname='YSE',simdir=youngdumpfile.split('/')[-2],snid=cd)
				idx = np.where(X.snid == cd)[0]
				sn = X.iloc[idx]

				#yseobs = np.min(np.abs(sn.MJD[(sn.FLT != 'X') & (sn.FLT != 'Y')]-(pkm+explosiondict[ni])))
				#ztfobs = np.min(np.abs(sn.MJD[(sn.FLT == 'X') | (sn.FLT == 'Y')]-(pkm+explosiondict[ni])))
				yseobs = np.min(np.abs(np.array(sn['mjd'])[0][(np.array(sn['FLT'])[0] != 'X') &
															  (np.array(sn['FLT'])[0] != 'Y')]-(pkm+explosiondict[ni])))
				ztfobs = np.min(np.abs(np.array(sn['mjd'])[0][(np.array(sn['FLT'])[0] == 'X') |
															  (np.array(sn['FLT'])[0] == 'Y')]-(pkm+explosiondict[ni])))

				if yseobs < 1.5 and ztfobs < 1.5:
					twodaycount += 1
					twodayidx += [noniadict[ni]]
					twodayz += [z]


		nonia_idx_oneday,counts_oneday = np.unique(onedayidx,return_counts=True)
		nonia_idx_twoday,counts_twoday = np.unique(twodayidx,return_counts=True)


		iibct,iinct,iipct,ibct,icct,iact,iaxct,slsnct,tdect = 0,0,0,0,0,0,0,0,0
		for ix,ct in zip(nonia_idx_oneday,counts_oneday):
			if ix == 'IIb': iibct += ct
			elif ix == 'IIn': iinct += ct
			elif ix == 'IIP': iipct += ct
			elif ix == 'Ib': ibct += ct
			elif ix == 'Ic': icct += ct
			elif ix == 'Ia': iact += ct
			elif ix == 'Iax': iaxct += ct
			elif ix == 'SLSN': slsnct += ct
			elif ix == 'TDe': tdect += ct
			
		iibct_2day,iinct_2day,iipct_2day,ibct_2day,icct_2day,iact_2day,iaxct_2day,slsnct_2day,tdect_2day = 0,0,0,0,0,0,0,0,0
		for ix,ct in zip(nonia_idx_twoday,counts_twoday):
			if ix == 'IIb': iibct_2day += ct
			elif ix == 'IIn': iinct_2day += ct
			elif ix == 'IIP': iipct_2day += ct
			elif ix == 'Ib': ibct_2day += ct
			elif ix == 'Ic': icct_2day += ct
			elif ix == 'Ia': iact_2day += ct
			elif ix == 'Iax': iaxct_2day += ct
			elif ix == 'SLSN': slsnct_2day += ct
			elif ix == 'TDe': tdect_2day += ct

		ax = plt.subplot(121)
		ax.hist(twodayz,label='%i SNe w/ 2 days, ZTF+PS1'%twodaycount)
		ax.hist(onedayz,label='%i SNe w/ 1 day, ZTF+PS1'%onedaycount)
		ax.legend()
		ax.set_xlabel('$z$',fontsize=15)
		ax.set_xlabel('# of SNe',fontsize=15)

		labels = 'IIb','IIn','IIP','Ib','Ic','Ia','Iax','SLSN-I','TDe'
		sizes = [iibct,iinct,iipct,ibct,icct,iact,iaxct,slsnct,tdect]
		colors = ['#1BE7FF','#6EEB83','#E4FF1A','#E8AA14','#FF5714','#731963','C0','C1','C2']
 
		# Plot
		ax = plt.subplot(122)
		ax.set_title('SN types detected w/i \n1 day of explosion \nby BOTH ZTF and PS1')
		def absolute_value(val):
			a  = np.round(val/100.*np.sum(sizes), 0)
			return '%i'%a
		ax.pie(sizes, labels=labels, colors=colors,
				autopct=absolute_value, shadow=False, startangle=140)
			

def read_data(filename):
	"""Read data from pickled file to a pandas dataframe"""
	with gzip.open(filename, 'rb') as f:
		data = pickle.load(f)

	X = to_dataframe(data)
	y = pd.get_dummies(X.type == 0, prefix='SNIa', drop_first=True)
	X = X.drop(columns=['type'])

	return X, y

def to_dataframe(data):
	"""Converts from a python dictionary to a pandas dataframe"""
	for idx in data:
		sn = data[idx]
		for filt in 'grizXY':
			sn['mjd_%s' % filt] = np.array(sn[filt]['mjd'])
			sn['fluxcal_%s' % filt] = np.array(sn[filt]['fluxcal'])
			sn['fluxcalerr_%s' % filt] = np.array(sn[filt]['fluxcalerr'])
#			if sn['header']['SIM_NON1a'] != 0:
			sn['days_from_explosion_%s' % filt] = \
				sn['mjd_%s' % filt] - sn['header']['SIM_PEAKMJD'] - explosiondict[sn['header']['SIM_NON1a']]
#		else:
#			sn['days_from_explosion_%s' % filt] = np.array([-99]*len(sn['mjd_%s' % filt]))

			del sn[filt]
		sn['mjd'] = np.concatenate((sn['mjd_g'],sn['mjd_r'],sn['mjd_i'],sn['mjd_z'],sn['mjd_X'],sn['mjd_Y'],))
		sn['FLT'] = np.concatenate((['g']*len(sn['mjd_g']),['r']*len(sn['mjd_r']),
									['i']*len(sn['mjd_i']),['z']*len(sn['mjd_z']),
									['X']*len(sn['mjd_X']),['Y']*len(sn['mjd_Y']),))

		sn.update(sn['header'])
		del sn['header']

	return pd.DataFrame.from_dict(data, orient='index')

if __name__ == "__main__":

	import os
	import optparse
	ys = YSFoM()

	usagestring = 'YSFoM.py <options>'
	parser = ys.add_options(usage=usagestring)
	options,  args = parser.parse_args()

	if options.serialize:
		os.system('serialize_lsst_model.py %s'%os.path.dirname(os.path.expandvars(options.youngdumpfile)))
		os.system('serialize_lsst_model.py %s'%os.path.dirname(os.path.expandvars(options.dumpfile)))
		
	if 'hi':
		plt.rcParams['figure.figsize'] = (25,5)
		plt.clf()
		ys.earlydetect(options.youngdumpfile,pkmagcut=options.pkmagcut,zcut=options.zcut,
					   pklfile='%s.pkl.gz'%os.path.dirname(os.path.expandvars(options.youngdumpfile)).split('/')[-1])
		ys.sn_plus_cosmo(options.dumpfile,options.fitresfile,pkmagcut=options.pkmagcut,zcut=options.zcut)
		plt.savefig(options.histplot,dpi = 600)

		plt.close()
		plt.rcParams['figure.figsize'] = (10,10)
		ys.youngLCplots(options.youngdumpfile,
						pklfile='%s.pkl.gz'%os.path.dirname(os.path.expandvars(options.youngdumpfile)).split('/')[-1])
		plt.savefig('%s_2day.png'%options.lcplot.split('.')[0],dpi = 600)
		plt.close()
		plt.rcParams['figure.figsize'] = (10,10)
		ys.sampleLCplots(options.dumpfile,
						 pklfile='%s.pkl.gz'%os.path.dirname(os.path.expandvars(options.dumpfile)).split('/')[-1])
		plt.savefig(options.lcplot,dpi = 600)

	# summary figs and
	plt.close()
	plt.rcParams['figure.figsize'] = (10,5)
	ys.ZTF_YSE(options.youngdumpfile,pkmagcut=options.pkmagcut,zcut=options.zcut,
			   pklfile='%s.pkl.gz'%os.path.dirname(os.path.expandvars(options.youngdumpfile)).split('/')[-1])
	plt.savefig(options.yseztfplot,dpi = 600)
	
# metrics:
# SNe < 1 day old
# SNe < 2 day old
# total SNe (Ia and CC)
# cosmologically useful SNe Ia
# median redshift (all, and cosmo-Ia)
