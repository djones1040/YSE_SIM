#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Alexandre Boucaud <aboucaud@lal.in2p3.fr>
# Modified by D. Jones

import argparse
from pathlib import Path
import os
import time

from .validutils.io import save_compressed
from .validutils.table import parse_model

def main(simname,outfile_root=None,dirpath='$SNDATA_ROOT/SIM',
		 verbose=False,save=True,filters=None):

	tstart = time.time()
	if not outfile_root:
		outfile_root = simname
	
	# Use model directory as base name for the output file
	filename = "{}.pkl.gz".format(outfile_root)
	dirpath = os.path.expandvars('%s/%s'%(dirpath,simname))
	
	table = parse_model(dirpath,filters=filters)
	if save: save_compressed(table, filename)

	if verbose:
		print("SN data from {} saved to {}".format(dirpath, filename))
		print('serializer took %.1f minutes'%((time.time()-tstart)/60.))

	return table
		
if __name__ == '__main__':
	main()
	
