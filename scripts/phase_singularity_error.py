#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

# _________________________________________________________________________________________________
def linear_fit(x, y):
  coef = np.polyfit(x, y, 1)
  poly1d_fn = np.poly1d(coef)
  return poly1d_fn(x)

# _________________________________________________________________________________________________
def main(args):
  simdir = args["simdir"]

  trace_rd  = np.load("{}/trajectory_data_rd.npy".format(simdir))#*0.25
  trace_pie = np.load("{}/trajectory_data_pie.npy".format(simdir))

  trace_t = np.arange(340, 5000, 2)
  trace_dif_mse = np.mean((trace_rd - trace_pie)**2, axis=1)
  trace_dif_rmse = np.sqrt(trace_dif_mse)
  lin_mse = linear_fit(trace_t, trace_dif_mse)
  lin_rmse = linear_fit(trace_t, trace_dif_rmse)

  # plot visualization
  ax = []
  rows, cols = (1, 2)
  #plt.rc('text', usetex=True)
  #plt.rc('font', family='Times New Roman', size=14)
  fig = plt.figure(figsize=(10, 5))
  grd = gs.GridSpec(rows, cols)

  for row in range(rows):
    for col in range(cols):
      it = row*cols + col
      ax.append(fig.add_subplot(grd[row, col]))

  ax[0].plot(trace_t, trace_dif_mse, color='black')
  ax[0].plot(trace_t, lin_mse, ls='-', color='red')
  ax[0].set(xlabel='Time [ms]', ylabel='MSE [$mm^{2}$]', title="Phase Singularity Error over time")

  ax[1].plot(trace_t, trace_dif_rmse, color='black')
  ax[1].plot(trace_t, lin_rmse, ls='-', color='red')
  ax[1].set(xlabel='Time [ms]', ylabel='$L^2$ Distance [mm]', title="Phase Singularity Error over time")

  pdf_filepath = "{}/A2_singularity_error.pdf".format(args["outdir"])
  plt.tight_layout()
  plt.savefig(pdf_filepath, dpi=200, bbox_inches='tight')
  plt.show()

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Visualization of phase singularity trajectory error.')
  parser.add_argument('--simdir', help='Simulation directory.', required=True)
  parser.add_argument('--outdir', help='Output directory.',     required=True)
  args = vars(parser.parse_args())
  main(args)
