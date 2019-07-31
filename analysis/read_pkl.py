import gzip
import pickle
import pandas as pd
from explosion import explosiondict_plasticc
import numpy as np

def read_data(filename):
	"""Read data from pickled file to a pandas dataframe"""
	with gzip.open(filename, 'rb') as f:
		data = pickle.load(f)
	X = to_dataframe(data)
	y = pd.get_dummies(X.type == 0, prefix='SNIa', drop_first=True)
	X = X.drop(columns=['type'])

	return X, y

def to_dataframe(data,simkey1='SIM_TYPE_NAME',simkey2='SIM_MODEL_NAME',simkey3='SIM_TEMPLATE_INDEX'):
	"""Converts from a python dictionary to a pandas dataframe"""
	goodcount,badcount = 0,0
	for idx in data:
		sn = data[idx]
		for filt in 'grizXY':
			sn['mjd_%s' % filt] = np.array(sn[filt]['mjd'])
			sn['fluxcal_%s' % filt] = np.array(sn[filt]['fluxcal'])
			sn['fluxcalerr_%s' % filt] = np.array(sn[filt]['fluxcalerr'])
			sn['photflag_%s' % filt] = np.array(sn[filt]['photflag'])

			exkey = '%s_%s_%i'%(sn['header'][simkey1].replace(' ',''),sn['header'][simkey2].replace(' ',''),sn['header'][simkey3])

			try:
				sn['days_from_firstlight_%s' % filt] = \
					sn['mjd_%s' % filt] - sn['header']['SIM_PEAKMJD'] - explosiondict_plasticc[exkey]
				goodcount += 1
			except:
				badcount += 1
				sn['days_from_firstlight_%s' % filt] = sn['mjd_%s' % filt] - sn['header']['SIM_PEAKMJD'] + 999

			del sn[filt]
		sn['mjd'] = np.concatenate((sn['mjd_g'],sn['mjd_r'],sn['mjd_i'],sn['mjd_z'],sn['mjd_X'],sn['mjd_Y'],))
		sn['FLT'] = np.concatenate((['g']*len(sn['mjd_g']),['r']*len(sn['mjd_r']),
									['i']*len(sn['mjd_i']),['z']*len(sn['mjd_z']),
									['X']*len(sn['mjd_X']),['Y']*len(sn['mjd_Y']),))

		sn.update(sn['header'])
		sn['SIM_TYPE'] = sn['header']['SIM_MODEL_NAME'].split('.')[-1].replace(' ','')
		if 'WFIRST' in sn['SIM_TYPE']:
			sn['SIM_TYPE'] = 'SNIa'
		elif 'MOSFIT' in sn['SIM_TYPE'] or 'NMF' in sn['SIM_TYPE']:
			sn['SIM_TYPE'] = sn['header']['SIM_MODEL_NAME'].split('.')[-1].replace(' ','').split('-')[0]
		if 'NON1ASED' in sn['SIM_TYPE']:
			sn['SIM_TYPE'] = sn['header']['SIM_TYPE_NAME'].replace(' ','')
		if sn['SIM_TYPE'] in ['Ib','Ic']: sn['SIM_TYPE'] = 'SN%s'%sn['SIM_TYPE']
		if sn['SIM_TYPE'] == 'IIP': sn['SIM_TYPE'] = 'SNII'
		
		sn['FULLTYPE'] = exkey
		
		#if 'UNKNOWN' not in sn['header']['SIM_TY_NAME']:
		#	sn['SIM_TYPE'] = sn['header']['SIM_MODEL_NAME']
		#elif 'NON1ASED' in sn['SIM_TYPE']:
		#	sn['SIM_TYPE'] = sn['header']['SIM_TYPE_NAME'].split('.')[-1]
		if sn['SIM_TYPE'] not in ['SNIa-91bg','SNIax','CART','SNIIn','SNIb','SNIc','SNII','SNIa']: import pdb; pdb.set_trace()
		
		del sn['header']
	print('warning: could not figure out explosion times for %i of %i transients (%.2f%%)'%(badcount,goodcount,badcount/goodcount*100))
	return pd.DataFrame.from_dict(data, orient='index')
