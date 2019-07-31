#!/usr/bin/env python
# D. Jones - 3/6/17
"""get luminosity function adjustment for 
different subtypes"""

subtypedict = {"IIb":{'SN1993J':601,
					  'SN2008ax':603,
					  'SN2011dh':602},
			   "IIP":{'2013am':604,
					  '2013ej':605,
					  '2012aw':606},
			   "IIn":{'2015bh':607,
					  '2011ht':608,
					  '2009ip':609},
			   "Ib":{'iPTF13bvn':610,
					 '2007Y':611},
			   "Ic":{'2007gr':612,
					 '2013ge':613},
			   "Ia":{'SN2011fe':614,
					 'SN2014J':615,
					 'SN2005hk':616}}

from hstsnpipe.tools import snana
import cosmo
import numpy as np
import pylab as plt

bands = 'BVRI'
Nsim = 3000
dR = 0.25
Rbins = np.arange( -23,-8, dR )
dB = 0.05
Bbins = np.arange( -19.7,-18, dB )

zsim = 2.4350179e-9	 # 10 pc for H0=73
zsim = 0.0015	# min z allowed in SNANA = 0.0001
mu = cosmo.mu( zsim, H0=73 ) 

def plotgauss( ax, mu, sigma, scale=1, **kwarg ) :
	""" plot a gaussian on the given axes """
	dx = 0.01
	x = np.arange( -23, -13, dx)
	gauss = np.exp( -(x-mu)**2/(2*sigma**2) )
	gauss = scale * gauss / (gauss.sum() *dx)
	ax.plot( x, gauss, **kwarg )

def doIbc( simroot='lftest', smear=True, dustmodel='none',clobber=False, verbose=1 ): 
	# simulate BVRI band light curves at z=0.001  with no mag smearing and no dust
	plt.clf()
	# Run the simulation
	if not clobber : 
		try : simIbc = snana.SimTable('youngsn')
		except : clobber = True 
	if clobber : 
		simIbc = snana.SimTable('youngsn')

	# read in results and convert to absolute mags
	iIb = np.where( simIbc.SIM_SUBTYPE=='Ib' )
	RIb = simIbc.SIM_PEAKMAG_R[iIb] - mu

	iIc = np.where( (simIbc.SIM_SUBTYPE=='Ic') ) 
	RIc = simIbc.SIM_PEAKMAG_R[iIc] - mu

	# Plot histograms
	RIbhist, RIbedges = np.histogram( RIb , bins=Rbins )
	RIchist, RIcedges = np.histogram( RIc , bins=Rbins )
	
	# Li:2011 / Smartt:2009 : the fraction of all SNe that are Ib/c
	fIbc = 0.20	 # Ib/c fraction 
	fIb = 0.07	 # Ib fraction
	fIc = 0.13	 # Ic sub-fraction includes some peculiars
	plt.plot( RIbedges[:-1], fIb * RIbhist/(dR*RIbhist.sum()), 'g-', drawstyle='steps-mid', lw=2, label='Ic')
	plt.plot( RIcedges[:-1], fIc * RIchist/(dR*RIchist.sum()), 'b-', drawstyle='steps-mid', lw=2, label='Ib')

	ax = plt.gca()

	plotgauss( ax, -17.01, 0.41, fIb, color='k',ls='--', lw=1.5 )  # Type Ib from Li+ 2011
	plotgauss( ax, -16.04, 1.28, fIc, color='k',ls='--', lw=1.5 )  # Type Ic from Li+ 2011
	ax.text(0.95,0.95, 'Type Ib/c Luminosity Functions\n from Li+2011',
			ha='right',va='top',transform=ax.transAxes)


	ax.text( 0.7, 0.4, 'Ic', color='b', ha='right',va='top', fontweight='bold',fontsize='x-large')
	ax.text( 0.4, 0.3, 'Ib', color='g', ha='right',va='top', fontweight='bold',fontsize='x-large')
	ax.set_xlabel( 'absolute R band magnitude')
	ax.legend(loc='upper left')
	
	ax.set_xlim([-20,-12])
	print('avg Ib, Ic mags')
	print(np.median(RIb),np.median(RIc))
	print('suggested Ib, Ic offsets')
	print(-(np.median(RIb)+17.01),-(np.median(RIc)+16.04))
	#print('suggested Ib individual offsets')
	#print(-(np.unique(RIb)+17.01))
	
	return(ax)

def doII( simroot='lftest', smear=True, dustmodel='none',clobber=False, verbose=1 ): 
	# simulate BVRI band light curves at z=0.001  with no mag smearing and no dust
	plt.clf()
	# Run the simulation
	if not clobber : 
		try : simII = snana.SimTable('youngsn')
		except : clobber = True 
	if clobber : 
		simII = snana.SimTable('youngsn')

	# read in results and convert to absolute mags
	iIIb = np.where( simII.SIM_SUBTYPE=='IIb' )
	RIIb = simII.SIM_PEAKMAG_R[iIIb] - mu
	IIbIDX = simII.SIM_NON1a[iIIb]
	
	iIIP = np.where( simII.SIM_SUBTYPE=='IIP' )
	RIIP = simII.SIM_PEAKMAG_R[iIIP] - mu
	IIPIDX = simII.SIM_NON1a[iIIP]
	
	iIIn = np.where( simII.SIM_SUBTYPE=='IIn' )
	RIIn = simII.SIM_PEAKMAG_R[iIIn] - mu
	IInIDX = simII.SIM_NON1a[iIIn]

	# Plot histograms
	RIIbhist, RIIbedges = np.histogram( RIIb , bins=Rbins )
	RIIPhist, RIIPedges = np.histogram( RIIP , bins=Rbins )
	RIInhist, RIInedges = np.histogram( RIIn , bins=Rbins )
	
	# Li:2011 / Smartt:2009 : the fraction of all SNe that are Ib/c
	fII = 0.57	 # Ib/c fraction 
	fIIb = 0.12	 # Ib fraction
	fIIP = 0.70	 # Ic sub-fraction includes some peculiars
	fIIn = 0.09	 # Ic sub-fraction includes some peculiars
	plt.plot( RIIbedges[:-1], fIIb * RIIbhist/(dR*RIIbhist.sum()),
			  'g-', drawstyle='steps-mid', lw=2, label='IIb')
	plt.plot( RIIPedges[:-1], fIIP * RIIPhist/(dR*RIIPhist.sum()),
			  'r-', drawstyle='steps-mid', lw=2, label='IIP')
	plt.plot( RIInedges[:-1], fIIn * RIInhist/(dR*RIInhist.sum()),
			  'p-', drawstyle='steps-mid', lw=2, label='IIn')

	ax = plt.gca()

	plotgauss( ax, -16.65, 1.30, fIIb, color='k',ls='--', lw=1.5 )  # Type Ib from Li+ 2011
	plotgauss( ax, -15.66, 1.23, fIIP, color='k',ls='--', lw=1.5 )  # Type Ic from Li+ 2011
	plotgauss( ax, -16.86, 1.61, fIIn, color='k',ls='--', lw=1.5 )  # Type Ic from Li+ 2011
	ax.text(0.95,0.95, 'Type II Luminosity Functions\n from Li+2011',
			ha='right',va='top',transform=ax.transAxes)


	ax.set_xlabel( 'absolute R band magnitude')
	ax.legend(loc='upper left')
	
	ax.set_xlim([-20,-12])
	print('avg IIb, IIP, IIn mags')
	print(np.median(RIIb),np.median(RIIP),np.median(RIIn))
	print('suggested change in std for IIb, IIP, IIn mags')
	print(1.30-np.std(RIIb),1.23-np.std(RIIP),1.61-np.std(RIIn))
	print('suggested IIb, IIP, IIn offsets')
	print(-(np.median(RIIb)+16.65),-(np.median(RIIP)+15.66),-(np.median(RIIn)+16.86))
	for IIbi in np.unique(IIbIDX):
	   print('suggested IIb individual offsets')
	   print(IIbi,-(np.median(RIIb[IIbIDX == IIbi])+16.65))
	print('suggested IIP individual offsets')
	for IIpi in np.unique(IIPIDX):
	   print('suggested IIP individual offsets')
	   print(IIpi,-(np.median(RIIP[IIPIDX == IIpi])+15.66))
	print('suggested IIn individual offsets')
	for IIni in np.unique(IInIDX):
	   print('suggested IIn individual offsets')
	   print(IIni,-(np.median(RIIn[IInIDX == IIni])+16.86))

	
	return(ax)

def doIa( simroot='lftest', smear=True, dustmodel='none',clobber=False, verbose=1 ): 
	# simulate BVRI band light curves at z=0.001  with no mag smearing and no dust
	plt.clf()
	# Run the simulation
	if not clobber : 
		try : simIa = snana.SimTable('youngsn')
		except : clobber = True 
	if clobber : 
		simIa = snana.SimTable('youngsn')

	# read in results and convert to absolute mags
	iIa = np.where( simIa.SIM_SUBTYPE=='Ia' )
	RIa = simIa.SIM_PEAKMAG_R[iIa] - mu
	IaIDX = simIa.SIM_NON1a[iIa]
	
	# Plot histograms
	RIahist, RIaedges = np.histogram( RIa , bins=Rbins )
	
	# Li:2011 / Smartt:2009 : the fraction of all SNe that are Ib/c
	fIa = 0.24	 # Ib/c fraction 
	plt.plot( RIaedges[:-1], fIa * RIahist/(dR*RIahist.sum()),
			  'g-', drawstyle='steps-mid', lw=2, label='Ia')

	ax = plt.gca()

	plotgauss( ax, -18.67, 0.51, fIa, color='k',ls='--', lw=1.5 )  # Type Ia from Li+ 2011
	ax.text(0.95,0.95, 'Type Ia Luminosity Function\n from Li+2011',
			ha='right',va='top',transform=ax.transAxes)


	ax.set_xlabel( 'absolute R band magnitude')
	ax.legend(loc='upper left')
	
	ax.set_xlim([-20,-12])
	print('avg Ia mags')
	print(np.median(RIa))
	print('suggested change in std for Ia mags')
	print(0.51-np.std(RIa))
	print('suggested Ia offset')
	print(-(np.median(RIa)+18.67))
	print('suggested Ia individual offsets')
	for Iai in np.unique(IaIDX):

	   print(Iai,-(np.median(RIa[IaIDX == Iai])+18.67))

	
	return(ax)



if __name__ == "__main__":
	doIbc()
