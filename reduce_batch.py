# reduction of KAST data "in batch"
# assumes that red reduction script is contained in "input.txt" 
# and blue reduction instructions is contained in "input_blue.txt"
# will skip reduction is these are missing
# default mode is to reduce everything; otherwise you can provide folder numbers

import os, numpy, scipy, pandas, glob, sys
from astropy.io import fits
from scipy.interpolate import interp1d

import kastredux
import numpy,scipy,copy
from astropy.io import fits
import matplotlib.pyplot as plt

base_folder = '/Users/caganze/research/J1624/data/'
red_input = 'kastRED_J1624-3212_20200722.fits'
blue_input = 'kastBLUE_J1624-3212_20200722.fits'

def initializeStandards(spt,sdss=True,reset=False,sd=False,esd=False,usd=False,beta=False,giant=False,gamma=False,verbose=False):
	if isinstance(spt,list) == False: spt = [spt]
	if isinstance(spt[0],float) == True or isinstance(spt[0],int) == True: 
		if len(spt) == 2: spt = numpy.arange(spt[0],spt[1]+1)
		spt = [kastredux.typeToNum(s) for s in spt] 

	if sd==True and 'sd' not in spt[0]: spt = ['sd'+s for s in spt]
	if esd==True and 'esd' not in spt[0]: spt = ['esd'+s for s in spt]
	if usd==True and 'usd' not in spt[0]: spt = ['usd'+s for s in spt]
	if beta==True and 'b' not in spt[0]: spt = [s+'b' for s in spt]
	if gamma==True and 'g' not in spt[0]: spt = [s+'g' for s in spt]
	if giant==True and 'I' not in spt[0]: spt = [s+'I' for s in spt]

	for s in spt:
		if s not in list(kastredux.SPTSTDS.keys()) or reset==True:
			f = numpy.array(glob.glob('{}/{}*.txt'.format(kastredux.SPTSTDFOLDER,s.replace('.',''))))
			if len(f) > 1:
				if sdss==True and len(f[['SDSS' in x for x in f]])>0: f = f[['SDSS' in x for x in f]]
# for subdwarfs
				elif (sd==True or esd==True or usd==True) and len(f[['lepine2007' in x for x in f]])>0: 
#					print(len(f))
					f = f[['lepine2007' in x for x in f]]
				else: pass
			if len(f) > 0:
				kastredux.SPTSTDS[s] = kastredux.readSpectrum(f[0],name='{} STD'.format(s))
			else: 
				if verbose==True: print('Warning: cannot find a spectral standard for type {}'.format(s))
	return

initializeStandards([50,75])
initializeStandards([50,75],sd=True)
initializeStandards([50,75],esd=True)
initializeStandards([50,75],usd=True)

#initializeStandards(['sdK{:.1f}'.format(n) for n in numpy.arange(0,10,0.5)])
#initializeStandards(['sdM{:.1f}'.format(n) for n in numpy.arange(0,10,0.5)])
#initializeStandards(['esdK{:.1f}'.format(n) for n in numpy.arange(0,10,0.5)])
#initializeStandards(['esdM{:.1f}'.format(n) for n in numpy.arange(0,10,0.5)])

def reduction_batch(dates):
	for date in dates:
		folder = '{}/{}/'.format(base_folder,date)
# red		
		instruct_file = '{}/{}'.format(folder,red_input)
		if os.path.exists(instruct_file):
			print('\n\nReducing RED data from {}\n\n'.format(date))
			par = kastredux.readInstructions(instruct_file)
# first run a profileCheck
			tmp = kastredux.profileCheck(instruct_file)
# now reduce
			redux = kastredux.reduce(instructions=instruct_file,bias_file='bias_RED.fits',flat_file='flat_RED.fits',mask_file='mask_RED.fits',cal_flux_file='cal_flux_RED.pkl',reset=True)#,cal_wave_file='cal_wave_BLUE.pkl')

# classify results
			print('\n\nClassifying RED data from {}\n\n'.format(date))
			sources = list(par['SOURCE'].keys())
			for src in sources:
				spec = kastredux.readSpectrum('{}/kastRED_{}_20{}.fits'.format(par['REDUCTION_FOLDER'],src,date))
				spec.name = src
				stats = []
				stds = list(kastredux.SPTSTDS.keys())
				for s in stds:
					st,scl = kastredux.compareSpectra_simple(spec,kastredux.SPTSTDS[s],fit_range=[6800,8800],plot=False)
					stats.append(st)
				spt = stds[numpy.argmin(stats)]
				kastredux.compareSpectra_simple(spec,kastredux.SPTSTDS[spt],fit_range=[6800,8800],plot=True,plot_file='{}/classify{}_{}.pdf'.format(par['REDUCTION_FOLDER'],par['MODE'],src))
#			except:
#				print('Problem reducing RED data from {}'.format(date))

# blue		
		instruct_file = '{}/{}'.format(folder,blue_input)
		if os.path.exists(instruct_file):
			print('Reducing BLUE data from {}'.format(date))
			par = kastredux.readInstructions(instruct_file)
# first run a profileCheck
			tmp = kastredux.profileCheck(instruct_file)
# now reduce
			redux = kastredux.reduce(instructions=instruct_file,bias_file='bias_BLUE.fits',flat_file='flat_BLUE.fits',mask_file='mask_BLUE.fits',cal_flux_file='cal_flux_BLUE.pkl',reset=True)#,cal_wave_file='cal_wave_BLUE.pkl')

# classify results
			# sources = list(par['SOURCE'].keys())
			# for src in sources:
			# 	spec = kastredux.readSpectrum('{}/kastBLUE_{}_20{}.fits'.format(par['REDUCTION_FOLDER'],src,date))
			# 	spec.name = src
			# 	stats = []
			# 	stds = list(kastredux.SPTSTDS.keys())
			# 	for s in stds:
			# 		st,scl = kastredux.compareSpectra_simple(spec,kastredux.SPTSTDS[s],fit_range=[4500,5300],plot=False)
			# 		stats.append(st)
			# 	spt = stds[numpy.argmin(stats)]
			# 	kastredux.compareSpectra_simple(spec,kastredux.SPTSTDS[spt],fit_range=[4500,5300],plot=True,plot_file='{}/classify{}_{}.pdf'.format(par['REDUCTION_FOLDER'],par['MODE'],src))
#			except:
#				print('Problem reducing BLUE data from {}'.format(date))

	return

# external function call
if __name__ == '__main__':
	if len(sys.argv) > 1: 
		dates = sys.argv[1:]
	else: 
		dates = glob.glob('{}/[12]*'.format(base_folder))
		dates = [a.split('/')[-1] for a in dates]

	if len(dates)>0: 
		reduction_batch(dates)
