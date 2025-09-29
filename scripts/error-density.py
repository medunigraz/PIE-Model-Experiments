#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

from scipy.stats import gaussian_kde
from carputils.carpio import igb


# _________________________________________________________________________________________________
def read_igb_file(filepath):
  data, _, _ = igb.read(filepath)
  return data[:,:11].ravel()


# _________________________________________________________________________________________________
def main(args):
  cmp_dir = args["compdir"]
  out_dir = args["outdir"]
  sim_id = cmp_dir.split("/")[-2]

  # load data
  lat_data = []
  vel_data = []
  lrt_data = []

  lat_data = read_igb_file("{}/lat_dif_ms.igb".format(cmp_dir))
  vel_data = read_igb_file("{}/vel_dif_m_s.igb".format(cmp_dir))
  lrt_data = read_igb_file("{}/lrt_dif_ms.igb".format(cmp_dir))

  all_data = [lat_data[::-1], vel_data[::-1], lrt_data[::-1]]

  # visualization properties
  axis_cfg = {"B1_restitution": [[-4.0, 4.0], [-0.02, 0.02], [ -4.0,  4.0]],
              "B2_curvature":   [[-4.0, 4.0], [-0.3,  0.3 ], [ -8.0,  8.0]],
              "B3_diffusion":   [[-8.0, 8.0], [-0.3,  0.3 ], [-32.0, 32.0]]}

  n_points = 100
  data_ranges = axis_cfg["B1_restitution"]
  data_label  = ["LAT", "CV", "LRT"]
  units       = ["ms", "m/s", "ms"]

  if sim_id in axis_cfg.keys():
    data_ranges = axis_cfg[sim_id]

  # plot visualization
  ax = []
  rows, cols = (1, 3)
  plt.rc('text', usetex=True)
  plt.rc('font', family='Times New Roman', size=14)
  fig = plt.figure(figsize=(12, 2))
  grd = gs.GridSpec(rows, cols)
  latex_str = sim_id
  
  for row in range(rows):
    for col in range(cols):
      it = row*cols + col
      ax.append(fig.add_subplot(grd[row, col]))

      pdata = all_data[col]
      print(pdata.shape)
      xmin, xmax = data_ranges[col]
      xx = np.linspace(xmin, xmax, n_points)
      ax[it].xaxis.set_major_locator(plt.MaxNLocator(5))
      xrange = np.arange(xmin, xmax+1e-6, (xmax-xmin)/4)
      xrange_labels = []

      if col != 1:
        xrange_labels = ["{:d}".format(int(item)) for item in xrange]
      else:
        xrange_labels = ["{:.2f}".format(item) for item in xrange]

      ax[it].set_xticks(xrange)
      ax[it].set_xticklabels(xrange_labels)

      #d = np.clip(pdata, xmin, xmax)
      d = pdata
      pdf = gaussian_kde(d)
      curve = pdf(xx)

      # data stats
      mean = np.mean(d)
      std  = np.std(d)
      rmse = np.sqrt(np.mean(d**2))
      vmin = np.amin(d)
      vmax = np.amax(d)
      print("{}:\tMEAN: {:.3f} STD: {:.3f} RMSE: {:.3f} ... MIN: {:.3f} MAX: {:.3f}".format(data_label[col], mean, std, rmse, vmin, vmax))
      latex_str += " & ${:.3f}$ & ${:.3f}$ & ${:.3f}$".format(mean, std, rmse)

      x_mean = mean
      y_mean = pdf(x_mean)
      polygon = ax[it].fill_between(xx, np.zeros(n_points), curve, zorder=len(pdata)+1, color='none')
      verts = np.vstack([p.vertices for p in polygon.get_paths()])
      gradient = ax[it].imshow(np.linspace(0, 1, 256).reshape(1, -1), cmap='turbo', aspect='auto', extent=[verts[:, 0].min(), verts[:, 0].max(), verts[:, 1].min(), verts[:, 1].max()])
      gradient.set_clip_path(polygon.get_paths()[0], transform=plt.gca().transData)

      ax[it].plot(xx, curve, c='k', zorder=len(pdata)+1)
      ax[it].vlines([x_mean], 0.0, y_mean, colors="black", linestyles=':',lw=2, zorder=len(pdata)+1)
      ax[it].spines['top'].set_visible(False)
      ax[it].spines['right'].set_visible(False)
      ax[it].spines['left'].set_visible(False)
      ax[it].set_xlabel("Signed {} Error [{}]".format(data_label[col], units[col]))
      ax[it].set_yticklabels([])
      ax[it].tick_params(axis='y', which='both', length=0)

      if col == 0:
        ax[it].title.set_text("{} Error Density".format(data_label[col]))
      else:
        ax[it].title.set_text("{} Error Density".format(data_label[col]))

  latex_str += " \\\\"
  print(latex_str)

  pdf_filepath = "{}/{}-densities.pdf".format(out_dir, sim_id)
  plt.tight_layout()
  plt.savefig(pdf_filepath, dpi=200, bbox_inches='tight')
  plt.show()

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Visualization of computed error density maps.')
  parser.add_argument('--compdir',   help='Compare directory.', required=True)
  parser.add_argument('--outdir',    help='Output directory.',  required=True)
  args = vars(parser.parse_args())
  main(args)