#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import matplotlib as mpl
from matplotlib import ticker
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, LogNorm

from scipy.stats import gaussian_kde
from carputils.carpio import igb

# _________________________________________________________________________________________________
def read_igb_file(filepath, tmax=11):
  data, _, _ = igb.read(filepath)
  return data[:,:tmax].ravel()

# _________________________________________________________________________________________________
def main(args):
  cases = ["A_restitution", "B_curvature", "C_diffusion"]
  case_labels = ["(A) Restitution Conditions", "(B) Curved Wavefronts", "(C) Diffusion Conditions"]

  # plot visualization
  ax = []
  rows, cols = (3, 3)
  plt.rc('text', usetex=True)
  plt.rc('font', family='Times New Roman', size=14)
  fig = plt.figure(figsize=(16, 12))
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

      ticks = np.array(new_ticks[it])
      xymin = ticks[0]
      xymax = ticks[-1]

      # line in reference paper
      ax[-1].plot([xymin, xymax], [xymin, xymax], color='lightgray', zorder=0) #, ls=':'
      #ax[-1].scatter(data_rd, data_pie, s=2, color='red')

      gist_heat_r = mpl.colormaps['hot'] # mpl.colormaps['gist_heat_r']
      newcmp = ListedColormap(gist_heat_r(np.linspace(0.2, 0.7, 128)[::-1]))

      if False:
        # (i) 2D KDE plot
        xx, yy = np.mgrid[xymin:xymax:100j, xymin:xymax:100j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([data_rd, data_pie])
        kernel = gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)
        lvls = np.linspace(np.amin(f), 1e-6, 8)
        eps = 1e-6#(np.amax(f) - np.amin(f))/100.0
        cfset = ax[-1].contourf(xx, yy, f+eps, cmap=newcmp, norm=LogNorm(vmin=np.amin(f)+eps, vmax=np.amax(f)+eps))
        cb = fig.colorbar(cfset, ax=ax[-1], label='density')
      else:
        # (ii) hexbin plot
        hb = ax[-1].hexbin(data_rd, data_pie, gridsize=70, cmap=newcmp, mincnt=1, vmin=1)
        bmin, bmax = int(np.amin(hb.get_array())), int(np.amax(hb.get_array()))
        print(bmin, bmax)
        cb = fig.colorbar(hb, ax=ax[-1], label='counts', ticks=[1, bmax])

      ax[-1].set_xticks(ticks)
      ax[-1].set_yticks(ticks)#, labels=[f'\\${x:1.2f}' for x in xticks]
      #ax[-1].plot(data_rd, data_pie, 'o', markersize=2, color='red', markevery=10)
      ax[-1].set_aspect('equal')
      ax[-1].set_xlim([xymin, xymax])
      ax[-1].set_ylim([xymin, xymax])
  
  pdf_filepath = "{}/outliers.pdf".format(args["outdir"])
  plt.tight_layout()

  plt.savefig(pdf_filepath, dpi=300, bbox_inches='tight')
  plt.savefig("{}/outliers.png".format(args["outdir"]), dpi=200, bbox_inches='tight')
  plt.show()

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Visualization of computed error density maps.')
  parser.add_argument('--outdir',    help='Output directory.',  required=True)
  args = vars(parser.parse_args())
  main(args)