#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

from scipy.signal import periodogram
from scipy.stats import gaussian_kde
from carputils.carpio import igb
from matplotlib.ticker import ScalarFormatter

# pie-solver --compare --mshA=msh_rd --latA=lat_rd --lrtA=lrt_rd --mshB=msh_pie --latB=lat_pie --lrtB=lrt_rd --out=cmpdir --nodal --np=32
# Example: python3 scripts/frequency_maps.py --simdir=./results/sim/A2_functional

# _________________________________________________________________________________________________
def read_igb_file(igb_filepath):
  data, _, _ = igb.read(igb_filepath)
  return data

def compute_pdf(data, xmin, xmax, N=100):
  xx = np.linspace(xmin, xmax, N)
  pdf = gaussian_kde(data.ravel())
  return xx, pdf(xx)

# _________________________________________________________________________________________________
def main(args):
  #cases = ["A1_structural", "A2_functional", "A3_wholeheart"]
  simdir = args["simdir"]

  # Extend ScalarFormatter
  class MyScalarFormatter(ScalarFormatter):
      # Override '_set_format' with your own
      def _set_format(self):
          self.format = '%d'  # Show 2 decimals
  custom_formatter = MyScalarFormatter(useMathText=True)

  # dominant frequency from PCL using PDF
  #pcl_rd_igb  = "{}/comp_1000/pclA_ms.igb".format(simdir)
  #pcl_pie_igb = "{}/comp_1000/pclB_ms.igb".format(simdir)
  #pcl_rd  = read_igb_file(pcl_rd_igb)[:,1:]
  #pcl_pie = read_igb_file(pcl_pie_igb)[:,1:]
  
  #fmap_rd  = 1e3/pcl_rd
  #fmap_pie = 1e3/pcl_pie

  #ftrace_rd  = compute_pdf(fmap_rd, 0.0, 8.0)
  #ftrace_pie = compute_pdf(fmap_pie, 0.0, 8.0)

  # dominant frequency from Vm using PSD
  vm_rd_igb = "{}/rd_250um/sim/vm.igb".format(simdir)
  vm_pie_igb  = "{}/pie_1000um/vm_mv.igb".format(simdir)
  vm_rd  = read_igb_file(vm_rd_igb).T
  vm_pie = read_igb_file(vm_pie_igb).T

  f_rd,  PSD_rd  = periodogram(vm_rd,  fs=1e3, axis=0)
  f_pie, PSD_pie = periodogram(vm_pie, fs=1e3, axis=0)

  amax_rd  = np.argmax(PSD_rd,  axis=0)
  amax_pie = np.argmax(PSD_pie, axis=0)

  df_rd =  f_rd[amax_rd]
  df_pie = f_pie[amax_pie]
  
  Npts = int(vm_rd.shape[1]/2)
  Mpts = int(np.sqrt(Npts))
  dfmap_rd = df_rd[:Npts].reshape((Mpts, Mpts))

  Npts = int(vm_pie.shape[1]/2)
  Mpts = int(np.sqrt(Npts))
  dfmap_pie = df_pie[:Npts].reshape((Mpts, Mpts))

  # plot visualization
  ax = []
  rows, cols = (1, 3)
  #plt.rc('text', usetex=True)
  #plt.rc('font', family='Times New Roman', size=14)
  fig = plt.figure(figsize=(15, 5))
  grd = gs.GridSpec(rows, cols)

  for row in range(rows):
    for col in range(cols):
      it = row*cols + col
      ax.append(fig.add_subplot(grd[row, col]))

  #fmin, fmax = (2.5, 4.5)
  fmin, fmax = (0.0, 16.0)

  #ax[0].plot(ftrace_rd[0],  ftrace_rd[1])
  #ax[0].plot(ftrace_pie[0], ftrace_pie[1])
  #ax[0].set(xlim=(fmin, fmax), xlabel='Frequency [Hz]', ylabel='Intensity [-]', title="Dominant Frequency (PCL+PDF)")
  #ax[1].imshow(phasemap_rd[1000], origin='lower') #cmap='hsv', 
  avg_PSD_rd = np.mean(PSD_rd, axis=1)
  avg_PSD_pie = np.mean(PSD_pie, axis=1)
  #std_PSD_rd = np.std(PSD_rd, axis=1)
  #std_PSD_pie = np.std(PSD_pie, axis=1)

  ax[0].semilogy(f_rd,  avg_PSD_rd, label="RD")
  ax[0].semilogy(f_pie, avg_PSD_pie, label="PIE")
  # plotting STD ranges in PSD is problematic when AVG - STD <= 0!
  #lines = ax[0].get_lines()
  #ax[0].fill_between(f_rd, np.maximum(avg_PSD_rd-std_PSD_rd, 1e-4), avg_PSD_rd+std_PSD_rd ,alpha=0.3, facecolor=lines[0].get_color())
  #ax[0].fill_between(f_pie, np.maximum(avg_PSD_pie-std_PSD_pie, 1e-4), avg_PSD_pie+std_PSD_pie ,alpha=0.3, facecolor=lines[1].get_color())
  ax[0].set(xlim=(fmin, fmax+1.0), ylim=(1e-2, 1e4), xlabel='Frequency [Hz]', ylabel='PSD [$V_{m}^{2}$/Hz]', title="Averaged PSD")
  ax[0].xaxis.set_major_formatter(custom_formatter)
  ax[0].set_xticks(np.arange(0, 17, 1))
  ax[0].legend(loc="lower right")

  dfmin, dfmax = (2.5, 7.5)

  im_rd = ax[1].imshow(dfmap_rd, vmin=dfmin, vmax=dfmax, origin='lower', cmap='plasma')
  #fig.colorbar(im_rd, orientation='vertical')
  ax[1].set_title("Dominant Frequency Map (RD)")
  ax[1].set_xticks([])
  ax[1].set_yticks([])

  im_pie = ax[2].imshow(dfmap_pie, vmin=dfmin, vmax=dfmax, origin='lower', cmap='plasma', aspect="auto")#, extent=[0,320,0,320]
  fig.colorbar(im_pie, orientation='vertical')
  ax[2].set_title("Dominant Frequency Map (PIE)")
  ax[2].set_xticks([])
  ax[2].set_yticks([])

  #print(np.amin(dfmap_rd), np.amax(dfmap_rd), np.amin(dfmap_pie), np.amax(dfmap_pie)) # -> 2.5 - 7.5

  pdf_filepath = "{}/A2_dominant_frequencies.pdf".format(args["outdir"])
  plt.tight_layout()
  plt.savefig(pdf_filepath, dpi=200, bbox_inches='tight')
  plt.show()

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Visualization of dominant frequency maps.')
  parser.add_argument('--simdir', help='Simulation directory.', required=True)
  parser.add_argument('--outdir', help='Output directory.',     required=True)
  args = vars(parser.parse_args())
  main(args)