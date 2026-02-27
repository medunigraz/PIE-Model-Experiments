#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

from scipy.stats import gaussian_kde
from carputils.carpio import igb

# _________________________________________________________________________________________________
def read_igb_file(filepath, tmax=11):
  data, _, _ = igb.read(filepath)
  return data[:,:tmax].ravel()

# _________________________________________________________________________________________________
def main(args):
  cases = ["B1_restitution", "B2_curvature", "B3_diffusion"]
  case_labels = ["Restitution Conditions (B1)", "Curved Wavefronts (B2)", "Diffusion Conditions (B3)"]

  # plot visualization
  ax = []
  rows, cols = (3, 3)
  #plt.rc('text', usetex=True)
  #plt.rc('font', family='Times New Roman', size=14)
  fig = plt.figure(figsize=(12, 12))
  grd = gs.GridSpec(rows, cols)

  new_ticks = [[0, 1000, 2000, 3000, 4000, 5000],
                [0.3, 0.35, 0.4, 0.45, 0.5, 0.55],
                [0, 1000, 2000, 3000, 4000, 5000],
                [0, 10, 20, 30, 40, 50, 60, 70],
                [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
                [300, 310, 320, 330, 340, 350, 360],
                [0, 20, 40, 60, 80, 100],
                [0.2, 0.4, 0.6, 0.8, 1.0],
                [300, 330, 360, 390, 420, 450]]
  
  for row in range(rows):
    for col in range(cols):
      it = row*cols + col
      ax.append(fig.add_subplot(grd[row, col]))
      #ax[-1].locator_params(nbins=5)
      #ax[-1].locator_params(axis='y', nbins=5)
      #ax[-1].locator_params(axis='x', nbins=5)
      #ax[-1].xaxis.set_major_locator(plt.MaxNLocator(3))
      #ax[-1].yaxis.set_major_locator(plt.MaxNLocator(3))

      compdir = []

      if row == 0:
        compdir = "./results/sim/{}/comp_1000_11".format(cases[0])
      elif row == 1:
        compdir = "./results/sim/{}/comp_1000".format(cases[1])
      elif row == 2:
        compdir = "./results/sim/{}/comp_1000".format(cases[2])
      else:
        print("Invalid testcase!")

      data_rd = []
      data_pie = []

      if col == 0:
        data_rd  = read_igb_file("{}/latA_ms.igb".format(compdir))
        data_pie = read_igb_file("{}/latB_ms.igb".format(compdir))
        ax[-1].set_xlabel("LAT RD  [ms]")
        ax[-1].set_ylabel("LAT PIE [ms]")
      elif col == 1:
        data_rd  = read_igb_file("{}/velA_m_s.igb".format(compdir))
        data_pie = read_igb_file("{}/velB_m_s.igb".format(compdir))
        ax[-1].set_xlabel("CV RD  [m/s]")
        ax[-1].set_ylabel("CV PIE [m/s]")
        #ax[-1].title.set_text(cases[row])
        ax[-1].title.set_text(case_labels[row])
      elif col == 2:
        data_rd  = read_igb_file("{}/lrtA_ms.igb".format(compdir))
        data_pie = read_igb_file("{}/lrtB_ms.igb".format(compdir))
        ax[-1].set_xlabel("LRT RD  [ms]")
        ax[-1].set_ylabel("LRT PIE [ms]")
      else:
        print("Invalid quantity!")

      sorted_idx = np.argsort(data_rd)
      data_rd  = data_rd[sorted_idx]
      data_pie = data_pie[sorted_idx]

      ax[-1].plot([data_rd[0], data_rd[-1]], [data_rd[0], data_rd[-1]], color='black')
      ax[-1].scatter(data_rd, data_pie, s=2, color='red')
      ax[-1].set_xticks(new_ticks[it])
      ax[-1].set_yticks(new_ticks[it])#, labels=[f'\\${x:1.2f}' for x in xticks]
      #ax[-1].plot(data_rd, data_pie, 'o', markersize=2, color='red', markevery=10)
  
  pdf_filepath = "{}/outliers.pdf".format(args["outdir"])
  plt.tight_layout()

  #plt.locator_params(nbins=5)
  #plt.locator_params(axis='y', nbins=5)
  #plt.locator_params(axis='x', nbins=5)

  plt.savefig(pdf_filepath, dpi=300, bbox_inches='tight')
  plt.savefig("{}/outliers.png".format(args["outdir"]), dpi=200, bbox_inches='tight')
  plt.show()

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Visualization of computed error density maps.')
  parser.add_argument('--outdir',    help='Output directory.',  required=True)
  args = vars(parser.parse_args())
  main(args)