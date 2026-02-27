#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

from scipy.signal import hilbert, correlate, convolve2d, find_peaks
from scipy.stats import gaussian_kde, pearsonr
from carputils.carpio import igb
from scipy.ndimage import gaussian_filter

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

# _________________________________________________________________________________________________
def read_igb_file(igb_filepath, verbose=False):
  if verbose:
    print("reading {} ...".format(igb_filepath))

  data, _, _ = igb.read(igb_filepath)
  return data

# _________________________________________________________________________________________________
def write_igb(igb_filepath, igb_data, verbose=False):
  if verbose:
    print("writing {} ...".format(igb_filepath))

  Ndat, Nslices = igb_data.shape
  igbobj = igb.IGBFile(igb_filepath, 'w')
  igb_header = {'x': Ndat, 'y': 1, 'z': 1, 't': Nslices,
                'type': 'float', 'systeme': 'little_endian'}
  igbobj.write(igb_data.astype(np.single()).T, igb_header)
  igbobj.close()

# _________________________________________________________________________________________________
def compute_phasemap_from_vm(vm_data, axis=-1, verbose=False):
  if verbose:
    print("computing phasemap ...")

  ht_data = hilbert(vm_data, axis=axis)
  return np.atan2(ht_data.imag, ht_data.real)

# _________________________________________________________________________________________________
def extract_frontface_data(data):
  Npts = int(data.shape[1]/2)
  Mpts = int(np.sqrt(Npts))
  return data[:,:Npts].reshape((data.shape[0], Mpts, Mpts))

# _________________________________________________________________________________________________
def znorm(ts):
  return (ts - np.mean(ts))/np.std(ts)

# _________________________________________________________________________________________________
def preprocess_signal(sig):
  peaks, _ = find_peaks(sig, distance=100)
  normed_sig = np.zeros(len(sig))

  for it in range(len(peaks) - 1):
    normed_sig[peaks[it]:peaks[it+1]] = znorm(sig[peaks[it]:peaks[it+1]])

  return normed_sig

# _________________________________________________________________________________________________
def compute_dtw(sigA, sigB, window=100, dstart=0.0):
  Na, Nb = D_shape = (len(sigA), len(sigB))
  D = np.ones(D_shape)*1e50
  D[0,0] = dstart
  w = max(window, np.abs(Na - Nb))
  
  for it in range(1, Na):
    for jt in range(max(1, it-w), min(Nb, it+w)):
      cost = (sigA[it-1] - sigB[jt-1])**2
      D[it,jt] = cost + min(min(D[it-1,jt], D[it,jt-1]), D[it-1,jt-1])
  
  return np.sqrt(D[it,jt]), D

# _________________________________________________________________________________________________
def compute_dtw_batched(sigA, sigB, window=100, peakdist=200, normalize=True):
  peaksA, _ = find_peaks(sigA, distance=peakdist)
  peaksB, _ = find_peaks(sigB, distance=peakdist)

  Np = min(len(peaksA), len(peaksB))

  peaksA = peaksA[:Np]
  peaksB = peaksB[:Np]
  total_dist = 0.0

  for it in range(Np - 1):
    sigAn = sigA[peaksA[it]:peaksA[it+1]]
    sigBn = sigB[peaksB[it]:peaksB[it+1]]

    if normalize:
      sigAn = znorm(sigAn)
      sigBn = znorm(sigBn)

    dist, Dmat = compute_dtw(sigAn, sigBn, window=window, dstart=0.0)
    total_dist += Dmat[-1,-1]

  return np.sqrt(total_dist)

# _________________________________________________________________________________________________
def dtw_task(it_item):
  #sigA, sigB = phasemap_dtw_data[it]
  sigA, sigB = it_item
  dtw_res = compute_dtw_batched(sigA, sigB, window=100, normalize=False)
  return dtw_res

# _________________________________________________________________________________________________
def main(args):
  simdir = args["simdir"]
  verbose = args["verbose"]
  vm_rd_igb = "{}/rd_250um/sim/vm.igb".format(simdir)
  vm_pie_igb  = "{}/pie_1000um/vm_mv.igb".format(simdir)

  # load Vm
  vm_rd  = read_igb_file(vm_rd_igb, verbose).T
  vm_pie = read_igb_file(vm_pie_igb, verbose).T

  # compute phasemaps for RD and PIE
  phasemap_rd = compute_phasemap_from_vm(vm_rd, 0, verbose)
  phasemap_pie = compute_phasemap_from_vm(vm_pie, 0, verbose)

  #write_igb("{}/rd_250um/sim/phasemap.igb".format(simdir), phasemap_rd.T, verbose)
  #write_igb("{}/pie_1000um/phasemap.igb".format(simdir), phasemap_pie.T, verbose)
  #return

  phasemap_rd = extract_frontface_data(phasemap_rd)
  #vm_pie = gaussian_filter(vm_pie[:,:Npts].reshape((vm_pie.shape[0], Mpts, Mpts)), sigma=1, axes=(1,2)).reshape((vm_pie.shape[0], Npts))
  phasemap_pie = extract_frontface_data(phasemap_pie)
  #print(phasemap_rd.shape)
  #print(phasemap_pie.shape)

  # compute correlation (frontface only!)
  phasemap_rd_matched = phasemap_rd[:,::4,::4]
  pcc_data, _ = pearsonr(phasemap_rd_matched, phasemap_pie, axis=0)
  print("mean", np.mean(pcc_data), "std", np.std(pcc_data), "min", np.amin(pcc_data), "max", np.amax(pcc_data), "med", np.median(pcc_data))

  '''
  Nt, Npie = phasemap_pie.shape[0], np.prod(phasemap_pie.shape[1:])
  phasemap_rd_dtw = phasemap_rd_matched.reshape((Nt, Npie))
  phasemap_pie_dtw = phasemap_pie.reshape((Nt, Npie))
  phasemap_dtw_data = np.moveaxis(np.array((phasemap_rd_dtw, phasemap_pie_dtw)).T, 2, 1)

  # Runs across multiple threads
  print("compute DTW ...")
  with ProcessPoolExecutor(max_workers=32) as executor: # ProcessPoolExecutor ThreadPoolExecutor
    #dtw_data = list(executor.map(dtw_task, [*range(Npie)]))
    dtw_data = list(executor.map(dtw_task, phasemap_dtw_data))

  dtw_data = np.array(dtw_data)
  dtw_data = dtw_data.reshape(pcc_data.shape)'''

  ts = 600

  if False:
    # phase singularity detection (alternative, experimental)
    grad_phase_rd = np.array(np.gradient(phasemap_rd, axis=(1,2), edge_order=1))
    grad_phase_pie = np.array(np.gradient(phasemap_pie, axis=(1,2), edge_order=1))
    #gx = phasemap_rd[:,1:,:] - phasemap_rd[:,:-1,:]
    #gy = phasemap_rd[:,:,1:] - phasemap_rd[:,:,:-1]
    #print(gx.shape, gy.shape)
    #grad_phase_rd = np.array([gx[:,:,:-1], gy[:,:-1,:]])
    #print(phasemap_rd.shape)
    print(grad_phase_rd.shape)
    
    #for it in range():
    #  test = convolve2d(grad_phase_pie[it], phase_kernel)
    kx = np.array([[ 1,-1],
                  [ 0, 0]])
    ky = np.array([[-1, 0],
                  [ 1, 0]])
    #kx = np.array([[-0.5,  0.0,  0.5],
    #               [-1.0,  0.0,  1.0],
    #               [-0.5,  0.0,  0.5]])
    #ky = np.array([[ 0.5,  1.0,  0.5],
    #               [ 0.0,  0.0,  0.0],
    #               [-0.5, -1.0, -0.5]])
    #kx = np.array([[ 1.0,  0.0, -1.0],
    #               [ 1.0,  0.0, -1.0],
    #               [ 0.0,  0.0,  0.0]])
    #ky = np.array([[-1.0, -1.0,  0.0],
    #               [ 0.0,  0.0,  0.0],
    #               [ 1.0,  1.0,  0.0]])
    #kx = -np.array([[ 0.5,  0.0,  0.0, -0.5],
    #               [ 1.0,  0.0,  0.0, -1.0],
    #               [ 1.0,  0.0,  0.0, -1.0],
    #               [ 0.5,  0.0,  0.0, -0.5]])
    #ky = np.array([[-0.5, -1.0, -1.0, -0.5],
    #               [ 0.0,  0.0,  0.0,  0.0],
    #               [ 0.0,  0.0,  0.0,  0.0],
    #               [ 0.5,  1.0,  1.0,  0.5]])
    print(kx)
    print(ky)
    testx = convolve2d(grad_phase_rd[0,ts], ky)
    testy = convolve2d(grad_phase_rd[1,ts], kx)

    '''grad_phase_pie = np.array(np.gradient(phasemap_pie, axis=(1,2)))
    print(grad_phase_pie.shape)
    
    #for it in range():
    #  test = convolve2d(grad_phase_pie[it], phase_kernel)
    phase_kernel = np.ones((2, 2))
    testx = convolve2d(grad_phase_pie[0,1000], phase_kernel)
    testy = convolve2d(grad_phase_pie[1,1000], phase_kernel)'''
    #path = testx[1:,1:] + testy[:-1,:-1]
    #path = testx[:-1,:-1] + testy[1:,1:]
    #path = testx[:-1,:-1] + testy[:-1,:-1]
    #path = testx[1:,1:] + testy[:-1,:-1]
    path = testx + testy
    print(testx.shape)
    print(testy.shape)

    #test = np.array([testx, testy])
    #print(test.shape)

    #path = np.linalg.norm(test, axis=0)
    print(path.shape)
    #asd 

  # plot visualization
  ax = []
  rows, cols = (2, 8)
  #plt.rc('text', usetex=True)
  #plt.rc('font', family='Times New Roman', size=14)
  fig = plt.figure(figsize=(24, 5))
  grd = gs.GridSpec(rows, cols)

  for row in range(rows):
    for col in range(cols-2):
      #it = row*cols + col
      ax.append(fig.add_subplot(grd[row, col]))
      ax[-1].set_xticks([])
      ax[-1].set_yticks([])
  
  ax.append(fig.add_subplot(grd[:, cols-2:]))
  ax[-1].set_xticks([])
  ax[-1].set_yticks([])

  '''
  vr = np.pi

  ax[0].plot(phasemap_rd[:,0,0])
  ax[0].plot(phasemap_pie[:,0,0])

  im1 = ax[1].imshow(phasemap_rd[1000], cmap='hsv', vmin=-vr, vmax=vr, origin='lower')
  im2 = ax[2].imshow(phasemap_pie[1000], cmap='hsv', vmin=-vr, vmax=vr, origin='lower')'''

  vr = np.pi
  im0 = ax[0].imshow(phasemap_rd[350],  vmin=-vr, vmax=vr, cmap='hsv', origin='lower')
  im1 = ax[1].imshow(phasemap_rd[500],  vmin=-vr, vmax=vr, cmap='hsv', origin='lower')
  im2 = ax[2].imshow(phasemap_rd[650],  vmin=-vr, vmax=vr, cmap='hsv', origin='lower')
  im3 = ax[3].imshow(phasemap_rd[1850], vmin=-vr, vmax=vr, cmap='hsv', origin='lower')
  im4 = ax[4].imshow(phasemap_rd[3250], vmin=-vr, vmax=vr, cmap='hsv', origin='lower')
  im5 = ax[5].imshow(phasemap_rd[5000], vmin=-vr, vmax=vr, cmap='hsv', origin='lower')

  im6  = ax[6].imshow(phasemap_pie[350], vmin=-vr, vmax=vr, cmap='hsv', origin='lower', extent=[0,320,0,320])#aspect='auto'
  im7  = ax[7].imshow(phasemap_pie[500], vmin=-vr, vmax=vr, cmap='hsv', origin='lower', extent=[0,320,0,320])#aspect='auto'
  im8  = ax[8].imshow(phasemap_pie[650], vmin=-vr, vmax=vr, cmap='hsv', origin='lower', extent=[0,320,0,320])#aspect='auto'
  im9  = ax[9].imshow(phasemap_pie[1850], vmin=-vr, vmax=vr, cmap='hsv', origin='lower', extent=[0,320,0,320])#aspect='auto'
  im10 = ax[10].imshow(phasemap_pie[3250], vmin=-vr, vmax=vr, cmap='hsv', origin='lower', extent=[0,320,0,320])#aspect='auto'
  im11 = ax[11].imshow(phasemap_pie[5000], vmin=-vr, vmax=vr, cmap='hsv', origin='lower', extent=[0,320,0,320])#aspect='auto'
  
  #im1 = ax[1].imshow(grad_phase_rd[0,1000], vmin=-vr, vmax=vr, origin='lower') #cmap='hsv'
  #im2 = ax[2].imshow(grad_phase_rd[1,1000], vmin=-vr, vmax=vr, origin='lower') #cmap='hsv'
  #im1 = ax[1].imshow(testx, vmin=-vr, vmax=vr, origin='lower') #cmap='hsv'
  #im2 = ax[2].imshow(testy, vmin=-vr, vmax=vr, origin='lower') #cmap='hsv'
  #im3 = ax[3].imshow(path, vmin=-2*vr, vmax=2*vr, origin='lower') #, cmap='hsv'

  #im2 = ax[1].imshow(phasemap_rd[400], cmap='hsv', vmin=-vr, vmax=vr)#, cmap='Greys' np.pi
  #im2 = ax[1].imshow(np.atan2(vm_rd[400], ht_rd[400]), cmap='hsv', vmin=-vr, vmax=vr)#, cmap='Greys' np.pi
  
  #im1 = ax[0].imshow(phasemap_rd[50], cmap='hsv', vmin=-vr, vmax=vr)#, cmap='Greys' np.pi
  #im2 = ax[1].imshow(phasemap_rd[400], cmap='hsv', vmin=-vr, vmax=vr)#, cmap='Greys'
  #print(vm_rd[2].shape)
  #im1 = ax[0].imshow(vm_rd[50], cmap='hsv', vmin=-90, vmax=40)#, cmap='Greys'
  #im2 = ax[1].imshow(vm_rd[400], cmap='hsv', vmin=-90, vmax=40)#, cmap='Greys'

  #ax[2].plot(phasemap_rd[:,0,0])
  #ax[2].plot(phasemap_pie[:,0,0])

  #iml = ax[-1].imshow(pcc_data, vmin=-1.0, vmax=1.0, origin='lower')#, cmap='Greys' np.pi
  #im3 = ax[3].imshow(dtw_data, origin='lower')#, cmap='Greys' np.pi
  #fig.colorbar(iml, orientation='vertical')

  iml = ax[-1].imshow(phasemap_rd[350], vmin=-vr, vmax=vr, cmap='hsv', origin='lower', extent=[0,320,0,320])
  fig.colorbar(iml, orientation='vertical')

  #fig.colorbar(im0, orientation='vertical')
  #fig.colorbar(im1, orientation='vertical')
  #fig.colorbar(im2, orientation='vertical')
  #fig.colorbar(im3, orientation='vertical')
  pdf_filepath = "{}/A2_phase_correlation_v2.pdf".format(args["outdir"])
  plt.tight_layout()
  plt.savefig(pdf_filepath, dpi=200, bbox_inches='tight')
  plt.show()

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Visualization of phasefield correlation.')
  parser.add_argument('--simdir',  help='Simulation directory.', required=True)
  parser.add_argument('--outdir',  help='Output directory.',     required=True)
  parser.add_argument('--verbose', help='Control verbosity.',    default=True)
  args = vars(parser.parse_args())
  main(args)
