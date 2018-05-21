#!/usr/bin/env python
# D. Jones - 3/28/18
# color-color plots to identify colors of young SNe
# python FindYoungSN.py --simname yse_ztf --empirical

from YSFoM import read_data
import numpy as np
import pylab as plt
from fill_between import fill_between_steps
LSST_FILTERS = 'grizXY'

explosiondict = {605: -8.0, 604: -7.0, 607: -8.0, 606: -6.1,
				 601: -15.1, 603: -38.0, 602: -40.1, 612: -9.1,
				 613: -6.1, 610: -16.1, 611: -12.1, 609: -14.1,
				 608: -42.1, 614: -17.0, 617: -14.0, 618: -19.1,
				 620: -94.1, 621: -64.0, 616: -10.0, 619: -15.1,
				 0:-15, 615:-13}

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


def plot_lightcurves(idx, X, ax):

    flux = []
    for f in LSST_FILTERS:
        if len(X.iloc[idx]['fluxcal_' + f]) > 0:
            flux.append(max(X.iloc[idx]['fluxcal_' + f]))

    for id_f, f in enumerate(LSST_FILTERS):


        #ax = axes[id_f // 3, id_f % 3]
        ax.errorbar(X.iloc[idx]['mjd_%s' % f] - X.iloc[idx]['pkmjd'],
                    X.iloc[idx]['fluxcal_%s' % f]/max(flux),
                    X.iloc[idx]['fluxcalerr_%s' % f]/max(flux),
                    fmt='o')
        ax.set_xlabel('MJD')
        ax.set_ylabel('Calibrated flux')
        ax.set_title('%s-band' % f)

class YoungSNPlots:
	def __init__(self):
		pass

	def add_options(self, parser=None, usage=None, config=None):
		import optparse
		if parser == None:
			parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

		parser.add_option(
			'-s','--simname', default='yse_ztf_perfect',type="string",
			help='name of sim (default=%default)')
		parser.add_option(
			'-t','--sntype', default='all',type="string",
			help='type of SN to plot (default=%default)')			
		parser.add_option(
			'--serialize', default=False,action="store_true",
			help='convert SNANA sims to pickle files, if set (default=%default)')
		parser.add_option(
			'--perfect', default=False,action="store_true",
			help='look at early light curves for the perfect simulations, if set (default=%default)')
		parser.add_option(
			'--empirical', default=False,action="store_true",
			help='test empirical light curve cuts, if set (default=%default)')

		return(parser)

	def risehist_perfect(self,pklfile,sntype='Ia'):
		X, y = read_data(pklfile)

		plt.clf()
		plt.ion()
		plt.rcParams['figure.figsize'] = (13,7)
		axg = plt.subplot(141); axr = plt.subplot(142); axi = plt.subplot(143); axz = plt.subplot(144)

		gdays1,gdays2,gdays3,gdays4 = np.array([]),np.array([]),np.array([]),np.array([])
		rdays1,rdays2,rdays3,rdays4 = np.array([]),np.array([]),np.array([]),np.array([])
		idays1,idays2,idays3,idays4 = np.array([]),np.array([]),np.array([]),np.array([])
		zdays1,zdays2,zdays3,zdays4 = np.array([]),np.array([]),np.array([]),np.array([])
		for i in range(len(X.snid)):
			if X.iloc[i]['z'] > 0.1: continue
			if (noniadict[X.iloc[i]['SIM_NON1a']] == sntype or sntype == 'all') and X.iloc[i]['SIM_NON1a'] != 0:
				print(X.iloc[i]['SIM_NON1a'],X.iloc[i]['z'])
				sn = X.iloc[i]
				iPlot0 = (sn['days_from_explosion_g'] >= -1) & \
						 (sn['days_from_explosion_g'] <= 0) & \
						 (sn['fluxcal_i']/sn['fluxcalerr_i'] > 10)
				iPlot1 = (sn['days_from_explosion_g'] >= 0) & \
						 (sn['days_from_explosion_g'] <= 1) & \
						 (sn['fluxcal_i']/sn['fluxcalerr_i'] > 10)
				iPlot2 = (sn['days_from_explosion_g'] >= 1) & \
				         (sn['days_from_explosion_g'] <= 2) & \
						 (sn['fluxcal_i']/sn['fluxcalerr_i'] > 10)
				iPlot3 = (sn['days_from_explosion_g'] >= 2) & \
				         (sn['days_from_explosion_g'] <= 3) & \
						 (sn['fluxcal_i']/sn['fluxcalerr_i'] > 10)
				iPlot4 = (sn['days_from_explosion_g'] > 3) & \
						 (sn['fluxcal_i']/sn['fluxcalerr_i'] > 10)


				gmag = -2.5*np.log10(sn['fluxcal_g'])+27.5
				rmag = -2.5*np.log10(sn['fluxcal_r'])+27.5
				imag = -2.5*np.log10(sn['fluxcal_i'])+27.5
				zmag = -2.5*np.log10(sn['fluxcal_z'])+27.5
				gdays1 = np.append(gdays1,gmag[iPlot0]-gmag[iPlot1])
				gdays2 = np.append(gdays2,gmag[iPlot1]-gmag[iPlot2])
				gdays3 = np.append(gdays3,gmag[iPlot2]-gmag[iPlot3])
				for day in sn['days_from_explosion_g'][iPlot4]:
					gdays4 = np.append(gdays4,gmag[sn['days_from_explosion_g'] == day-1] - gmag[sn['days_from_explosion_g'] == day])

				rdays1 = np.append(rdays1,rmag[iPlot0]-rmag[iPlot1])
				rdays2 = np.append(rdays2,rmag[iPlot1]-rmag[iPlot2])
				rdays3 = np.append(rdays3,rmag[iPlot2]-rmag[iPlot3])
				for day in sn['days_from_explosion_g'][iPlot4]:
					rdays4 = np.append(rdays4,rmag[sn['days_from_explosion_g'] == day-1] - rmag[sn['days_from_explosion_g'] == day])

				idays1 = np.append(idays1,imag[iPlot0]-imag[iPlot1])
				idays2 = np.append(idays2,imag[iPlot1]-imag[iPlot2])
				idays3 = np.append(idays3,imag[iPlot2]-imag[iPlot3])
				for day in sn['days_from_explosion_g'][iPlot4]:
					idays4 = np.append(idays4,imag[sn['days_from_explosion_g'] == day-1] - imag[sn['days_from_explosion_g'] == day])

				if len(zmag) == len(gmag):
					zdays1 = np.append(zdays1,zmag[iPlot0]-zmag[iPlot1])
					zdays2 = np.append(zdays2,zmag[iPlot1]-zmag[iPlot2])
					zdays3 = np.append(zdays3,zmag[iPlot2]-zmag[iPlot3])
					for day in sn['days_from_explosion_g'][iPlot4]:
						zdays4 = np.append(zdays4,zmag[sn['days_from_explosion_g'] == day-1] - zmag[sn['days_from_explosion_g'] == day])
				#elif len(zmag) == len(gmag)-1:
				#	zdays1 = np.append(zdays1,zmag[iPlot0[:-1]]-zmag[iPlot1[:-1]])
				#	zdays2 = np.append(zdays2,zmag[iPlot1[:-1]]-zmag[iPlot2[:-1]])
				#	zdays3 = np.append(zdays3,zmag[iPlot2[:-1]]-zmag[iPlot3[:-1]])
				#	for day in sn['days_from_explosion_g'][iPlot4[:-1]]:
				#		zdays4 = np.append(zdays4,zmag[sn['days_from_explosion_g'][:-1] == day-1] - zmag[sn['days_from_explosion_g'][:-1] == day])

		gdays1 = gdays1[gdays1 == gdays1]
		gdays2 = gdays2[gdays2 == gdays2]
		gdays3 = gdays3[gdays3 == gdays3]
		gdays4 = gdays4[gdays4 == gdays4]

		rdays1 = rdays1[rdays1 == rdays1]
		rdays2 = rdays2[rdays2 == rdays2]
		rdays3 = rdays3[rdays3 == rdays3]
		rdays4 = rdays4[rdays4 == rdays4]

		idays1 = idays1[idays1 == idays1]
		idays2 = idays2[idays2 == idays2]
		idays3 = idays3[idays3 == idays3]
		idays4 = idays4[idays4 == idays4]

		zdays1 = zdays1[zdays1 == zdays1]
		zdays2 = zdays2[zdays2 == zdays2]
		zdays3 = zdays3[zdays3 == zdays3]
		zdays4 = zdays4[zdays4 == zdays4]

		
		magbins = np.linspace(np.min([np.min(gdays1),np.min(gdays2),np.min(gdays3),np.min(gdays4)]),
							  np.max([np.max(gdays1),np.max(gdays2),np.max(gdays3),np.max(gdays4)]),20)
		axg.hist(gdays4[gdays4 == gdays4],bins=magbins,label='>3 days',color='C0')
		axg.hist(gdays3[gdays3 == gdays3],bins=magbins,label='<3 days',color='C1',histtype='step',lw=2,ls='--')
		axg.hist(gdays2[gdays2 == gdays2],bins=magbins,label='<2 days',color='C2',histtype='step',lw=2,ls='-.')
		axg.hist(gdays1[gdays1 == gdays1],bins=magbins,label='<1 days',color='C3',histtype='step',lw=2,ls='-')

		magbins = np.linspace(np.min([np.min(rdays1),np.min(rdays2),np.min(rdays3),np.min(rdays4)]),
							  np.max([np.max(rdays1),np.max(rdays2),np.max(rdays3),np.max(rdays4)]),20)
		axr.hist(rdays4[rdays4 == rdays4],bins=magbins,label='<4 days',color='C0')
		axr.hist(rdays3[rdays3 == rdays3],bins=magbins,label='<3 days',color='C1',histtype='step',lw=2,ls='--')
		axr.hist(rdays2[rdays2 == rdays2],bins=magbins,label='<2 days',color='C2',histtype='step',lw=2,ls='-.')
		axr.hist(rdays1[rdays1 == rdays1],bins=magbins,label='<1 days',color='C3',histtype='step',lw=2,ls='-')
		
		magbins = np.linspace(np.min([np.min(idays1),np.min(idays2),np.min(idays3),np.min(idays4)]),
							  np.max([np.max(idays1),np.max(idays2),np.max(idays3),np.max(idays4)]),20)
		axi.hist(idays4[idays4 == idays4],bins=magbins,label='<4 days',color='C0')
		axi.hist(idays3[idays3 == idays3],bins=magbins,label='<3 days',color='C1',histtype='step',lw=2,ls='--')
		axi.hist(idays2[idays2 == idays2],bins=magbins,label='<2 days',color='C2',histtype='step',lw=2,ls='-.')
		axi.hist(idays1[idays1 == idays1],bins=magbins,label='<1 days',color='C3',histtype='step',lw=2,ls='-')

		magbins = np.linspace(np.min([np.min(zdays1),np.min(zdays2),np.min(zdays3),np.min(zdays4)]),
							  np.max([np.max(zdays1),np.max(zdays2),np.max(zdays3),np.max(zdays4)]),20)
		axz.hist(zdays4[zdays4 == zdays4],bins=magbins,label='<4 days',color='C0')
		axz.hist(zdays3[zdays3 == zdays3],bins=magbins,label='<3 days',color='C1',histtype='step',lw=2,ls='--')
		axz.hist(zdays2[zdays2 == zdays2],bins=magbins,label='<2 days',color='C2',histtype='step',lw=2,ls='-.')
		axz.hist(zdays1[zdays1 == zdays1],bins=magbins,label='<1 days',color='C3',histtype='step',lw=2,ls='-')

		
		for ax,flt in zip([axg,axr,axi,axz],'griz'):
			ax.set_xlabel('%s rise (mag/day)'%flt)
			ax.set_ylabel('N$_{SNe}$')
			ax.set_yscale('log')
		axg.legend()

		
		plt.show()
		import pdb; pdb.set_trace()
		
	def colorcolor_perfect(self,pklfile,sntype='Ia'):
		X, y = read_data(pklfile)

		plt.clf()
		plt.ion()
		plt.rcParams['figure.figsize'] = (13,7)
		ax1 = plt.axes([0.1,0.1,0.3,0.6]); ax2 = plt.axes([0.55,0.1,0.3,0.6])
		ax1histy = plt.axes([0.4,0.1,0.07,0.6],sharey=ax1); ax1histx = plt.axes([0.1,0.7,0.3,0.12],sharex=ax1)
		ax2histy = plt.axes([0.85,0.1,0.07,0.6],sharey=ax2); ax2histx = plt.axes([0.55,0.7,0.3,0.12],sharex=ax2)
		ax1histx.xaxis.tick_top(); ax2histx.xaxis.tick_top()
		ax1histy.yaxis.tick_right(); ax2histy.yaxis.tick_right()
		
		cax = plt.axes([0.1,0.9,0.42,0.05])
		ax1.set_ylabel('$g - r$')
		ax1.set_xlabel('$r - i$')
		ax1.set_title('SN %s'%sntype.split('SN')[-1])
		ax2.set_ylabel('$r - i$')
		ax2.set_xlabel('$i - z$')
		ax2.set_title('SN %s'%sntype.split('SN')[-1])


		rmi,gmr,imz = np.array([]),np.array([]),np.array([])
		rmidays,gmrdays,imzdays = np.array([]),np.array([]),np.array([])
		count = 0
		for i in range(len(X.snid)):
			if X.iloc[i]['z'] > 0.1: continue
			if (noniadict[X.iloc[i]['SIM_NON1a']] == sntype or sntype == 'all') and X.iloc[i]['SIM_NON1a'] != 0:
				print(X.iloc[i]['SIM_NON1a'],X.iloc[i]['z'])
				iPlot = (X.iloc[i]['days_from_explosion_g'] >= 0) & \
				        (X.iloc[i]['days_from_explosion_g'] < 12) & \
						(X.iloc[i]['fluxcal_i']/X.iloc[i]['fluxcalerr_i'] > 10)
				mag_r = -2.5*np.log10(X.iloc[i]['fluxcal_r'][iPlot]) + 27.5
				mag_g = -2.5*np.log10(X.iloc[i]['fluxcal_g'][iPlot]) + 27.5
				mag_i = -2.5*np.log10(X.iloc[i]['fluxcal_i'][iPlot]) + 27.5
				if len(X.iloc[i]['fluxcal_z']) == len(X.iloc[i]['fluxcal_g'])-1:
					mag_z = -2.5*np.log10(X.iloc[i]['fluxcal_z'][iPlot[:-1]]) + 27.5

				ax1.plot(mag_r - mag_i,
						 mag_g - mag_r,'.')
				sc = ax1.scatter(mag_r - mag_i,mag_g - mag_r,
								 c=X.iloc[i]['days_from_explosion_g'][iPlot],s=10,zorder=30,cmap='Paired')
				rmi = np.append(rmi,mag_r - mag_i)
				rmidays = np.append(rmidays,X.iloc[i]['days_from_explosion_g'][iPlot])
				gmr = np.append(gmr,mag_g - mag_r)
				gmrdays = np.append(gmrdays,X.iloc[i]['days_from_explosion_g'][iPlot])
				
				if len(X.iloc[i]['fluxcal_z']) == len(X.iloc[i]['fluxcal_g'])-1:				
					ax2.plot(mag_i - mag_z,
							 mag_r - mag_i,'.')
					sc = ax2.scatter(mag_i - mag_z,mag_r - mag_i,
									 c=X.iloc[i]['days_from_explosion_g'][iPlot],s=10,zorder=30,cmap='Paired')

					imz = np.append(imz,mag_i - mag_z)
					imzdays = np.append(imzdays,X.iloc[i]['days_from_explosion_z'][iPlot[:-1]])
					
				if count == 0:
					cb = plt.colorbar(sc,cax=cax,orientation='horizontal')
					cax.set_xlabel('Days from Explosion for SN %s'%sntype,labelpad=-50)
				count += 1
				plt.show()

		gmrdays = gmrdays[gmr == gmr]
		gmr = gmr[(gmr == gmr)]
		imzdays = imzdays[imz == imz]
		imz = imz[(imz == imz)]

		rmihistbins = np.linspace(np.min(rmi),np.max(rmi),20)
		ax1histx.hist(rmi,bins=rmihistbins,alpha=0.5,color=plt.cm.Paired(12/12.),zorder=8)
		ax1histx.hist(rmi[rmidays <= 2],bins=rmihistbins,alpha=0.5,color=plt.cm.Paired(2/12.),zorder=8)
		ax1histx.hist(rmi[rmidays <= 1],bins=rmihistbins,alpha=0.5,color=plt.cm.Paired(1/12.),zorder=10)
		
		gmrhistbins = np.linspace(np.min(gmr[gmr == gmr]),np.max(gmr[gmr == gmr]),20)
		ax1histy.hist(gmr,bins=gmrhistbins,alpha=0.5,color=plt.cm.Paired(12/12.),zorder=8,orientation='horizontal')
		ax1histy.hist(gmr[gmrdays <= 2],bins=gmrhistbins,alpha=0.5,color=plt.cm.Paired(2/12.),zorder=8,orientation='horizontal')
		ax1histy.hist(gmr[gmrdays <= 1],bins=gmrhistbins,alpha=0.5,color=plt.cm.Paired(1/12.),zorder=10,orientation='horizontal')

		rmihistbins = np.linspace(np.min(rmi[rmi == rmi]),np.max(rmi[rmi == rmi]),20)
		ax2histy.hist(rmi,bins=rmihistbins,alpha=0.5,color=plt.cm.Paired(12/12.),zorder=8,orientation='horizontal')
		ax2histy.hist(rmi[rmidays <= 2],bins=rmihistbins,alpha=0.5,color=plt.cm.Paired(2/12.),zorder=8,orientation='horizontal')
		ax2histy.hist(rmi[rmidays <= 1],bins=rmihistbins,alpha=0.5,color=plt.cm.Paired(1/12.),zorder=10,orientation='horizontal')
		
		imzhistbins = np.linspace(np.min(imz[imz == imz]),np.max(imz[imz == imz]),20)
		ax2histx.hist(imz,bins=imzhistbins,alpha=0.5,color=plt.cm.Paired(12/12.),zorder=8)
		ax2histx.hist(imz[imzdays <= 2],bins=imzhistbins,alpha=0.5,color=plt.cm.Paired(2/12.),zorder=8)
		ax2histx.hist(imz[imzdays <= 1],bins=imzhistbins,alpha=0.5,color=plt.cm.Paired(1/12.),zorder=10)
		

		#gmrbinspace = (gmrhistbins[0]+gmrhistbins[1])/2.
		#hist1,bin_edges1 = np.histogram(gmr[gmrdays <= 1],bins=gmrhistbins)
		#hist2,bin_edges2 = np.histogram(gmr[gmrdays <= 2],bins=gmrhistbins)
		#fill_between_steps(ax1histy, hist1, [0]*len(hist1), bin_edges1[0:-1]+gmrbinspace,color=plt.cm.Paired(1/12.),alpha=0.5)
		#fill_between_steps(ax1histy, hist2, [0]*len(hist2), bin_edges2[0:-1]+gmrbinspace,color=plt.cm.Paired(2/12.),alpha=0.5)
		for ax in [ax1,ax2]:
			ax.set_xlim([-2,2])
			ax.set_ylim([-2,2])
			
		import pdb; pdb.set_trace()
		
		return

	def findYoung(self,pklfile,riseerrlim=0.25,ztfriselim=0.15,sntype='all'):
		X, y = read_data(pklfile)
		sntype,days_from_expl,grise,griseerr,rrise,rriseerr,\
			irise,iriseerr,zrise,zriseerr,\
			gXrise,gXriseerr,rYrise,rYriseerr,\
			griselim,rriselim,iriselim,zriselim,\
			gXriselim,rYriselim,gdeterr,rdeterr,ideterr,\
			   zdeterr,gXdeterr,rYdeterr = self.getRiseLists(X,y,sntypetouse=sntype)

		plt.clf()
		plt.ion()
		axg = plt.subplot(141); axr = plt.subplot(142); axgx = plt.subplot(143); axry = plt.subplot(144)
		
		# g rise
		magbins = np.linspace(np.min(grise[grise == grise]),
							  np.max(grise[grise == grise]),20)
		axg.hist(griselim[(days_from_expl >= 3) & (griselim == griselim) & (gdeterr < riseerrlim)],bins=magbins,
				 label='>3 days',color='C0',histtype='step',lw=2,ls='-')
		axg.hist(griselim[(days_from_expl >= 0) & (days_from_expl < 3) & (griselim == griselim) & (gdeterr < riseerrlim)],bins=magbins,
				 label='<3 days',color='C1',histtype='step',lw=2,ls='--')
		axg.hist(griselim[(days_from_expl >= 0) &(days_from_expl < 2) & (griselim == griselim) & (gdeterr < riseerrlim)],bins=magbins,
				 label='<2 days',color='C2',histtype='step',lw=2,ls='-.')
		axg.hist(griselim[(days_from_expl >= 0) & (days_from_expl < 1) & (griselim == griselim) & (gdeterr < riseerrlim)],bins=magbins,
				 label='<1 days',color='C3',histtype='step',lw=2,ls='-')
		
		# r rise
		magbins = np.linspace(np.min(rrise[rrise == rrise]),
							  np.max(rrise[rrise == rrise]),20)
		axr.hist(rriselim[(days_from_expl >= 3) & (rriselim == rriselim) & (rdeterr < riseerrlim)],bins=magbins,
				 label='>3 days',color='C0',histtype='step',lw=2,ls='-')
		axr.hist(rriselim[(days_from_expl >= 0) & (days_from_expl < 3) & (rriselim == rriselim) & (rdeterr < riseerrlim)],bins=magbins,
				 label='<3 days',color='C1',histtype='step',lw=2,ls='--')
		axr.hist(rriselim[(days_from_expl >= 0) & (days_from_expl < 2) & (rriselim == rriselim) & (rdeterr < riseerrlim)],bins=magbins,
				 label='<2 days',color='C2',histtype='step',lw=2,ls='-.')
		axr.hist(rriselim[(days_from_expl >= 0) & (days_from_expl < 1) & (rriselim == rriselim) & (rdeterr < riseerrlim)],bins=magbins,
				 label='<1 days',color='C3',histtype='step',lw=2,ls='-')

		# g-PS1/g-ZTF rise
		magbins = np.linspace(np.min(gXrise[(gXrise == gXrise) & (gXriseerr < ztfriselim)]),
							  np.max(gXrise[(gXrise == gXrise) & (gXriseerr < ztfriselim)]),10)
		axgx.hist(gXrise[(days_from_expl >= 3) & (days_from_expl < 15) & (gXriselim == gXriselim) & (gXriseerr < ztfriselim)],bins=magbins,
				 label='>3 days',color='C0',histtype='step',lw=2,ls='-')
		axgx.hist(gXrise[(days_from_expl >= 0) & (days_from_expl < 3) & (gXriselim == gXriselim) & (gXriseerr < ztfriselim)],bins=magbins,
				 label='<3 days',color='C1',histtype='step',lw=2,ls='--')
		axgx.hist(gXrise[(days_from_expl >= 0) & (days_from_expl < 2) & (gXriselim == gXriselim) & (gXriseerr < ztfriselim)],bins=magbins,
				 label='<2 days',color='C2',histtype='step',lw=2,ls='-.')
		axgx.hist(gXrise[(days_from_expl >= 0) & (days_from_expl < 1) & (gXriselim == gXriselim) & (gXriseerr < ztfriselim)],bins=magbins,
				 label='<1 days',color='C3',histtype='step',lw=2,ls='-')

		# r-PS1/R-ZTF rise
		magbins = np.linspace(np.min(rYrise[(rYrise == rYrise) & (rYriseerr < ztfriselim)]),
							  np.max(rYrise[(rYrise == rYrise) & (rYriseerr < ztfriselim)]),10)
		axry.hist(rYrise[(days_from_expl >= 3) & (days_from_expl < 15) & (rYriselim == rYriselim) & (rYriseerr < ztfriselim)],bins=magbins,
				  label='>3 days',color='C0',histtype='step',lw=2,ls='-')
		axry.hist(rYrise[(days_from_expl >= 0) & (days_from_expl < 3) & (rYriselim == rYriselim) & (rYriseerr < ztfriselim)],bins=magbins,
				  label='<3 days',color='C1',histtype='step',lw=2,ls='--')
		axry.hist(rYrise[(days_from_expl >= 0) & (days_from_expl < 2) & (rYriselim == rYriselim) & (rYriseerr < ztfriselim)],bins=magbins,
				  label='<2 days',color='C2',histtype='step',lw=2,ls='-.')
		axry.hist(rYrise[(days_from_expl >= 0) & (days_from_expl < 1) & (rYriselim == rYriselim) & (rYriseerr < ztfriselim)],bins=magbins,
				  label='<1 days',color='C3',histtype='step',lw=2,ls='-')


		for ax,flt in zip([axg,axr,axgx,axry],['g','r','gX','rY']):
			if flt in 'gr':
				ax.set_xlabel('$%s$ rise above prev. 3$\sigma$ limiting mag (3 days)'%flt)
				ax.set_xlim([-1,4])
			else:
				ax.set_xlabel('$%s$ rise (ZTF to PS1)'%flt[0])
				ax.set_xlim([-1,1])
			ax.set_ylabel('N$_{SNe}$',fontsize=15)
			ax.set_yscale('log')

		axg.legend()

		plt.show()
		import pdb; pdb.set_trace()
		
		# based on these histograms, we make some cuts and estimate completeness
		
	def getRiseLists(self,X,y,sntypetouse='all'):

		grise_arr,griseerr_arr,rrise_arr,rriseerr_arr,irise_arr,iriseerr_arr,zrise_arr,zriseerr_arr = \
			np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),\
			np.array([]),np.array([]),np.array([])
		gXrise_arr,gXriseerr_arr,rYrise_arr,rYriseerr_arr = \
			np.array([]),np.array([]),np.array([]),np.array([])
		griselim_arr,rriselim_arr,iriselim_arr,zriselim_arr,gXriselim_arr,rYriselim_arr = \
			np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
		gdeterr_arr,rdeterr_arr,ideterr_arr,zdeterr_arr,gXdeterr_arr,rYdeterr_arr = \
			np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
		days_from_expl,sntype = np.array([]),np.array([])
		
		for i in range(len(X.snid)):
			sn = X.iloc[i]
			if sn['SIM_NON1a'] == 0: continue
			if sntypetouse != 'all' and noniadict[sn['SIM_NON1a']] != sntypetouse: continue

			gmag = -2.5*np.log10(sn['fluxcal_g'])+27.5; gmagerr = 1.086*sn['fluxcalerr_g']/sn['fluxcal_g']
			rmag = -2.5*np.log10(sn['fluxcal_r'])+27.5; rmagerr = 1.086*sn['fluxcalerr_r']/sn['fluxcal_r']
			imag = -2.5*np.log10(sn['fluxcal_i'])+27.5; imagerr = 1.086*sn['fluxcalerr_i']/sn['fluxcal_i']
			zmag = -2.5*np.log10(sn['fluxcal_z'])+27.5; zmagerr = 1.086*sn['fluxcalerr_z']/sn['fluxcal_z']
			Xmag = -2.5*np.log10(sn['fluxcal_X'])+27.5; Xmagerr = 1.086*sn['fluxcalerr_X']/sn['fluxcal_X']
			Ymag = -2.5*np.log10(sn['fluxcal_Y'])+27.5; Ymagerr = 1.086*sn['fluxcalerr_Y']/sn['fluxcal_Y']
			
			mjdmin = 1e8
			for filt in 'griz':
				mjd = sn['mjd_%s'%filt]
				mjd_goodsnr = mjd[sn['fluxcal_%s'%filt]/sn['fluxcalerr_%s'%filt] > 5]
				if len(mjd_goodsnr) and np.min(mjd_goodsnr) < mjdmin: mjdmin = np.min(mjd_goodsnr)

			if mjdmin != 1e8:
				sntype = np.append(sntype,noniadict[sn['SIM_NON1a']])
				days_from_expl = np.append(days_from_expl,
										   mjdmin - sn['SIM_PEAKMJD'] - explosiondict[sn['SIM_NON1a']])

				# if detected in PS1, find how quickly it rose from previous epoch
				iG = np.where(sn['mjd_g'] == mjdmin)[0]
				iGbefore = np.where(sn['mjd_g'] == mjdmin-3)[0]
				if len(iG) and len(iGbefore):
					grise = -(gmag[iG] - gmag[iGbefore])
					griseerr = np.sqrt(gmagerr[iG]**2. + gmagerr[iGbefore]**2.)
					detmagerr = gmagerr[iG]
					griselim = -(gmag[iG] - (-2.5*np.log10(sn['fluxcal_g'][iGbefore]+3*sn['fluxcalerr_g'][iGbefore])+27.5))
					#if griselim > 2: import pdb; pdb.set_trace()
				else: grise = np.nan; griseerr = np.nan; griselim = np.nan; detmagerr = np.nan
				grise_arr = np.append(grise_arr,grise)
				griseerr_arr = np.append(griseerr_arr,griseerr)
				griselim_arr = np.append(griselim_arr,griselim)
				gdeterr_arr = np.append(gdeterr_arr,detmagerr)
				
				iR = np.where(sn['mjd_r'] == mjdmin)[0]
				iRbefore = np.where(sn['mjd_r'] == mjdmin-3)[0]
				if len(iR) and len(iRbefore):
					rrise = -(rmag[iR] - rmag[iRbefore])
					rriseerr = np.sqrt(rmagerr[iR]**2. + rmagerr[iRbefore]**2.)
					rriselim = -(rmag[iR] - (-2.5*np.log10(sn['fluxcal_r'][iRbefore]+3*sn['fluxcalerr_r'][iRbefore])+27.5))
					detmagerr = rmagerr[iR]
				else: rrise = np.nan; rriseerr = np.nan; rriselim = np.nan; detmagerr = np.nan
				rrise_arr = np.append(rrise_arr,rrise)
				rriseerr_arr = np.append(rriseerr_arr,rriseerr)
				rriselim_arr = np.append(rriselim_arr,rriselim)
				rdeterr_arr = np.append(rdeterr_arr,detmagerr)
				
				iI = np.where(sn['mjd_i'] == mjdmin)[0]
				iIbefore = np.where(sn['mjd_i'] == mjdmin-3)[0]
				if len(iI) and len(iIbefore):
					irise = -(imag[iI] - imag[iIbefore])
					iriseerr = np.sqrt(imagerr[iI]**2. + imagerr[iIbefore]**2.)
					iriselim = -(imag[iI] - (-2.5*np.log10(sn['fluxcal_i'][iIbefore]+3*sn['fluxcalerr_i'][iIbefore])+27.5))
					detmagerr = imagerr[iI]
				else: irise = np.nan; iriseerr = np.nan; iriselim = np.nan; detmagerr = np.nan
				irise_arr = np.append(irise_arr,irise)
				iriseerr_arr = np.append(iriseerr_arr,iriseerr)
				iriselim_arr = np.append(iriselim_arr,iriselim)
				ideterr_arr = np.append(ideterr_arr,detmagerr)
				
				iZ = np.where(sn['mjd_z'] == mjdmin)[0]
				iZbefore = np.where(sn['mjd_z'] == mjdmin-3)[0]
				if len(iZ) and len(iZbefore):
					zrise = -(zmag[iZ] - zmag[iZbefore])
					zriseerr = np.sqrt(zmagerr[iZ]**2. + zmagerr[iZbefore]**2.)
					zriselim = -(zmag[iZ] - (-2.5*np.log10(sn['fluxcal_z'][iZbefore]+3*sn['fluxcalerr_z'][iZbefore])+27.5))
					detmagerr = zmagerr[iZ]
				else: zrise = np.nan; zriseerr = np.nan; zriselim = np.nan; detmagerr = np.nan
				zrise_arr = np.append(zrise_arr,zrise)
				zriseerr_arr = np.append(zriseerr_arr,zriseerr)
				zriselim_arr = np.append(zriselim_arr,zriselim)
				zdeterr_arr = np.append(zdeterr_arr,detmagerr)
				
				# find how quickly it rose from ZTF 3 hours earlier
				iG = np.where(sn['mjd_g'] == mjdmin)[0]
				iXbefore = np.where(sn['mjd_X'] == mjdmin-0.12)[0]
				if len(iG) and len(iXbefore):
					gXrise = -(gmag[iG] - Xmag[iXbefore])
					gXriseerr = np.sqrt(gmagerr[iG]**2. + Xmagerr[iXbefore]**2.)
					gXriselim = -(gmag[iG] - (-2.5*np.log10(sn['fluxcal_X'][iXbefore]+3*sn['fluxcalerr_X'][iXbefore])+27.5))
					detmagerr = gmagerr[iG]
				else: gXrise = np.nan; gXriseerr = np.nan; gXriselim = np.nan; detmagerr = np.nan
				gXrise_arr = np.append(gXrise_arr,gXrise)
				gXriseerr_arr = np.append(gXriseerr_arr,gXriseerr)
				gXriselim_arr = np.append(gXriselim_arr,gXriselim)
				gXdeterr_arr = np.append(gXdeterr_arr,detmagerr)
				#if gXrise > 0.5 and gXriseerr < 0.1: import pdb; pdb.set_trace()
				
				iR = np.where(sn['mjd_r'] == mjdmin)[0]
				iYbefore = np.where(sn['mjd_Y'] == mjdmin-0.12)[0]
				if len(iR) and len(iYbefore):
					rYrise = -(rmag[iR] - Ymag[iYbefore])
					rYriseerr = np.sqrt(rmagerr[iR]**2. + Ymagerr[iYbefore]**2.)
					rYriselim = -(rmag[iR] - (-2.5*np.log10(sn['fluxcal_Y'][iYbefore]+3*sn['fluxcalerr_Y'][iYbefore])+27.5))
					detmagerr = rmagerr[iR]
				else: rYrise = np.nan; rYriseerr = np.nan; rYriselim = np.nan; detmagerr = np.nan
				rYrise_arr = np.append(rYrise_arr,rYrise)
				rYriseerr_arr = np.append(rYriseerr_arr,rYriseerr)
				rYriselim_arr = np.append(rYriselim_arr,rYriselim)
				rYdeterr_arr = np.append(rYdeterr_arr,detmagerr)
				
		return(sntype,days_from_expl,grise_arr,griseerr_arr,rrise_arr,rriseerr_arr,
			   irise_arr,iriseerr_arr,zrise_arr,zriseerr_arr,
			   gXrise_arr,gXriseerr_arr,rYrise_arr,rYriseerr_arr,
			   griselim_arr,rriselim_arr,iriselim_arr,zriselim_arr,
			   gXriselim_arr,rYriselim_arr,gdeterr_arr,rdeterr_arr,ideterr_arr,
			   zdeterr_arr,gXdeterr_arr,rYdeterr_arr)
		
if __name__ == "__main__":
	import os
	import optparse
	
	ys = YoungSNPlots()
	usagestring = 'FindYoungSN.py <options>'
	parser = ys.add_options(usage=usagestring)
	options,  args = parser.parse_args()

	if options.serialize:
		os.system('serialize_lsst_model.py $SNDATA_ROOT/SIM/%s_YOUNG'%options.simname)

	plt.rcParams['figure.figsize'] = (10,5)
	if options.perfect:
		plt.rcParams['figure.figsize'] = (15,3)
		ys.risehist_perfect('%s_YOUNG.pkl.gz'%options.simname,sntype=options.sntype)
		plt.savefig('%s_%s_risehist_perfect.png'%(options.simname,options.sntype))

	if options.empirical:
		plt.rcParams['figure.figsize'] = (15,3)
		print('pklfile: %s_YOUNG.pkl.gz'%options.simname)
		ys.findYoung('%s_YOUNG.pkl.gz'%options.simname,sntype=options.sntype)
		plt.savefig('%s_risehist_empirical.png'%options.simname)
