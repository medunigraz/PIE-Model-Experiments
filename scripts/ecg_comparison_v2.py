#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import numpy as np
import matplotlib.text as mpltext
import matplotlib.lines as mpllines
import matplotlib.pyplot as mplplot
import matplotlib.gridspec as mplgridspec

from scipy.stats import pearsonr
from scipy.signal import find_peaks, correlate

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
def pcc_dtw_per_lead(ecgA, ecgB, leads, distance=200, window=100, title="", verbose=False):
  ncc_data = np.zeros(len(leads))
  pcc_data = np.zeros(len(leads))
  rmse_data = np.zeros(len(leads))
  dtw_data = np.zeros(len(leads))
  ecg_rd_pie_t = np.linspace(0.0, 5.0, len(ecgB['ecg']["I"]))

  for it, lead in enumerate(leads):
    lead_data_rd = np.array(ecgA['ecg'][lead], dtype=float)
    lead_data_pie = np.array(ecgB['ecg'][lead], dtype=float)
    fac = np.round(len(lead_data_rd)/len(lead_data_pie))
    valid_idx = (np.arange(len(lead_data_pie))*fac).astype(np.int32)
    lead_data_rd = lead_data_rd[valid_idx]

    ncc_data[it] = (correlate(lead_data_rd/lead_data_rd.std(), lead_data_pie/lead_data_pie.std(), mode="valid")/len(lead_data_rd))[0]
    pcc_data[it], _ = pearsonr(lead_data_rd, lead_data_pie)
    rmse_data[it] = np.mean((lead_data_rd - lead_data_pie)**2)
    dtw_data[it] = compute_dtw_batched(lead_data_rd, lead_data_pie, window=window, normalize=False)

    peaks_rd, _ = find_peaks(lead_data_rd, distance=distance)
    peaks_pie, _ = find_peaks(lead_data_pie, distance=distance)

    if verbose:
      mplplot.plot(ecg_rd_pie_t, lead_data_rd)
      mplplot.plot(ecg_rd_pie_t, lead_data_pie)
      mplplot.plot(peaks_rd/1000.0, lead_data_rd[peaks_rd], "x")
      mplplot.plot(peaks_pie/1000.0, lead_data_pie[peaks_pie], "x")
      mplplot.title("{} - {}".format(title, lead))
      mplplot.show()

  dtw_min = np.amin(dtw_data)
  dtw_max = np.amax(dtw_data)
  dtw_mean = np.mean(dtw_data)
  pcc_min = np.amin(pcc_data)
  pcc_max = np.amax(pcc_data)
  pcc_mean = np.mean(pcc_data)

  if verbose:
    print("{}: DTW: min {:.3f} max {:.3f} PCC: min {:.3f} max {:.3f}".format(title, dtw_min, dtw_max, pcc_min, pcc_max))

  # latex table:
  dtw_latex = "DTW"
  pcc_latex = "PCC"

  for it, lead in enumerate(leads):
    dtw_latex += " & ${:.3f}$".format(dtw_data[it])
    pcc_latex += " & ${:.3f}$".format(pcc_data[it])
  
  dtw_latex += " & ${:.3f}$".format(dtw_mean)
  pcc_latex += " & ${:.3f}$".format(pcc_mean)

  print(dtw_latex)
  print(pcc_latex)

# _________________________________________________________________________________________________
class LineCollection:
  LINE_STYLES = ('solid', 'dashed', 'dashdot', 'dotted')

  def __init__(self):
    self._visible = True
    self._style_idx = 0
    self._lines = list()
  
  def append(self, line):
    if isinstance(line, mpllines.Line2D):
      self._lines.append(line)

  def toggle_visibility(self):
    self._visible = not self._visible
    for line in self._lines:
      line.set_visible(self._visible)

  def next_style(self):
    self._style_idx = (self._style_idx + 1) % len(LineCollection.LINE_STYLES)
    for line in self._lines:
      line.set_linestyle(LineCollection.LINE_STYLES[self._style_idx])

  @property
  def line_style(self):
    return LineCollection.LINE_STYLES[self._style_idx]
  
  @property
  def visible(self):
    return self._visible

# _________________________________________________________________________________________________
def dat2dict(filepath):
  header = []
  with open(filepath, 'r') as file:
    header = file.readline()

  leads = header[:-1].split(" ")[1:]
  data = np.loadtxt(filepath, skiprows=1).T*1e-3
  jdata = {"t": data[0]}
  jecgs = {}
  data = data[1:]

  for it in range(len(leads)):
    jecgs[leads[it]] = data[it]

  jdata["ecg"] = jecgs
  
  return jdata

# _________________________________________________________________________________________________
def main():
  if not len(sys.argv) > 1:
    name = os.path.basename(sys.argv[0])
    print('usage: "{} ECG0.json [ECG1.json ...]"\n'.format(name))
    exit(1)

  width = 20
  height = 7.5
  #height = 10.0
  #width = 2.0*height
  # width = 1.618*height
  files = list()

  for arg in sys.argv[1:]:
    arg_data = arg.split(':')
    file, scale, alpha = arg_data.pop(0), 1.0, 1.0

    if len(arg_data) > 0:
      scale = float(arg_data.pop(0))

    if len(arg_data) > 0:
      alpha = float(arg_data.pop(0))

    if file.endswith('.json') and os.path.exists(file):
      files.append((file, scale, alpha))
      continue

    if file.endswith('.dat') and os.path.exists(file):
      files.append((file, scale, alpha))
      continue

    print('unknown argument "{}" !'.format(arg))

  ecgs, leads = dict(), set()

  for file, scale, alpha in files:
    jdata = []

    if file.endswith('.json'):
      fp = open(file, 'r')
      jdata = json.load(fp)
      fp.close()
    
    if file.endswith('.dat'):
      jdata = dat2dict(file)
    
    #jdata['ecg']['aVR'] = -1.0*np.array(jdata['ecg']['aVR'], dtype=float) # we stored maVR in aVR

    for lead in jdata['ecg'].keys():
      leads.add(lead)

    ecgs[file] = (jdata, scale, alpha)

  leads = tuple(["I", "II", "III", "maVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"])

  # compute PCC and DTW
  #ecg_rd, _, _  = ecgs["./results/ecg/ecg_rd_250um.json"]
  #ecg_pie, _, _ = ecgs["./results/ecg/ecg_pie_1000um.json"]
  #pcc_dtw_per_lead(ecg_rd, ecg_pie, leads, distance=200, window=100, title="DTW and PCC")
  
  lcols = ['black', 'red', 'goldenrod', 'mediumvioletred', 'dodgerblue', 'goldenrod', 'mediumvioletred', 'dodgerblue']
  #lcols = ['black', 'red', 'mediumblue', 'royalblue', 'deepskyblue', 'darkgreen', 'limegreen', 'lightgreen']
  lalph = 0.8*np.array([255.0, 255.0, 255.0, 255.0, 255.0, 255.0, 255.0, 255.0])/255.0

  #lcols = ['black', 'red', 'royalblue', 'royalblue', 'royalblue', 'limegreen', 'limegreen', 'limegreen']
  #lalph = np.array([255.0, 255.0, 200.0, 150.0, 100.0, 200.0, 150.0, 100.0])/255.0
  vlins = np.array([125, 500, 1950, 3650, 4890, 5000])*1e-3

  # plot
  na = 3
  nb = 4
  fig = mplplot.figure(figsize=(width, height))
  gs = mplgridspec.GridSpec(na, nb, left=0.01, right=0.99, top=0.99, bottom=0.01, hspace=0.0, wspace=0.0)
  axes = dict()
  ax_idx, ref_axes = [0, 0], None

  for k, lead in enumerate(leads):
    if ref_axes is None:
      ax = fig.add_subplot(gs[ax_idx[0], ax_idx[1]])
      ref_axes = ax
    else:
      ax = fig.add_subplot(gs[ax_idx[0], ax_idx[1]], sharex=ref_axes, sharey=ref_axes)
      
    ax_idx[0] += 1

    if ax_idx[0] >= na:
      ax_idx[0] = 0
      ax_idx[1] += 1

    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])
    #ax.spines[['left', 'right', 'bottom', 'top']].set_visible(False)
    ax.axhline(0.0, color='#000000', alpha=0.5)
    axes[lead] = ax

  lbl2coll_map = dict()
  ecg_it = 0

  for label, (ecg, scale, alpha) in ecgs.items():
    line_coll = LineCollection()

    for lead in ecg['ecg']:
      lead_data = np.array(ecg['ecg'][lead], dtype=float)*scale

      line = []
      #if ecg_it == 0 or (ecg_it > 1 and ecg_it < 5):
      #if ecg_it == 0 or ecg_it > 4:
      line, = axes[lead].plot(ecg['t'], lead_data, color=lcols[ecg_it], alpha=lalph[ecg_it], label=label, linestyle=LineCollection.LINE_STYLES[0])
      #if ecg_it == 0:
      #  line.set_zorder(20)

      line_coll.append(line)

    lbl2coll_map[label] = line_coll
    ecg_it = ecg_it + 1
  
  for lead in leads:
    #axes[lead].text(0.01, 0.95, lead, transform=axes[lead].transAxes, size=12)
    tbox = axes[lead].text(0.007, 0.91, lead, transform=axes[lead].transAxes, size=16)
    tbox.set_bbox(dict(facecolor='white', alpha=1.0, edgecolor='black', linewidth=1.5))
    axes[lead].vlines(vlins, -1.5, 1.5, ls=":", color='gray')
    axes[lead].set(xlim=[0.0,5.0], ylim=[-1.5, 1.5])

  legline2coll_map = dict()

  if True:
    legend = ref_axes.legend()

    for leg_text, leg_line in zip(legend.get_texts(), legend.get_lines()):
      leg_text.set_picker(True)
      leg_line.set_picker(5.0)
      legline2coll_map[leg_line] = lbl2coll_map[leg_text.get_text()]

  def onpick(event):
    artist = event.artist

    if isinstance(artist, mpltext.Text):
      label = artist.get_text()

      if (line_coll := lbl2coll_map.get(label, None)) is not None:
        line_coll.toggle_visibility()
        alpha = 1.0 if line_coll.visible else 0.2
        artist.set_alpha(alpha)
        fig.canvas.draw()

    elif isinstance(artist, mpllines.Line2D):
      if (line_coll := legline2coll_map.get(artist, None)) is not None:
        line_coll.next_style()
        artist.set_linestyle(line_coll.line_style)
        fig.canvas.draw()

  fig.canvas.mpl_connect('pick_event', onpick)
  #pdf_filepath = "./results/ecg_comp_A.pdf"
  #pdf_filepath = "./results/ecg_comp_B.pdf"
  pdf_filepath = "./results/A3_wholeheart_ECG.pdf"
  mplplot.savefig(pdf_filepath, dpi=300, bbox_inches='tight')  
  mplplot.show()

# _________________________________________________________________________________________________
if __name__ == '__main__':
    main()
