#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

import neurokit2 as nk
from carputils.carpio import igb
from scipy.stats import pearsonr, gaussian_kde
from scipy.signal import argrelmax, find_peaks, butter, lfilter, freqz, correlate

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
  peaksA, _ = find_peaks(sigA, height=0.1, distance=peakdist)
  peaksB, _ = find_peaks(sigB, height=0.1, distance=peakdist)

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
def measure_RR_interval(lead_t, lead_f):
  peak_indices, _ = find_peaks(lead_f, height=0.1, distance=250)

  if False:
    print("Number of peaks: {:d}".format(len(peak_indices)))
    fig = plt.figure(figsize=(20, 8))
    plt.hlines([0.0], lead_t[0], lead_t[-1], ls=":", color="gray")
    plt.plot(lead_t, lead_f)
    plt.plot(lead_t[peak_indices], lead_f[peak_indices], "x")
    plt.tight_layout()
    plt.show()

  RR_intervals  = np.diff(lead_t[peak_indices])*1e3 # to ms
  RR_avg = np.mean(RR_intervals)
  RR_std = np.std(RR_intervals)

  return RR_avg, RR_std, peak_indices

# _________________________________________________________________________________________________
def load_json(filepath):
  dict = []

  with open(filepath, 'r') as file:
    dict = json.load(file)

  return dict

# _________________________________________________________________________________________________
def load_ecg(ecg_filepath):
  with open(ecg_filepath, 'r') as ecg_file:
    dict = json.load(ecg_file)

  t = np.array(dict["t"])
  ecg_dict = dict["ecg"]

  if "maVR" in ecg_dict.keys():
    ecg_dict['aVR'] = -np.array(ecg_dict.pop("maVR"))

  if t[1] - t[0] < 5e-3:
    t = t*1e3

  return t, ecg_dict

# _________________________________________________________________________________________________
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

# _________________________________________________________________________________________________
def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

# _________________________________________________________________________________________________
def main(args):
  ecgA_filepath = args["ecgA"]
  ecgB_filepath = args["ecgB"]
  output_dir = args["out"]
  dtw_distance = args["distance"]
  dtw_window   = args["window"]

  # load ECG data
  ecgA  = load_json(ecgA_filepath)
  ecgB = load_json(ecgB_filepath)

  # check if ECG data is complete (TODO: maVR - aVR)
  leads = tuple(["I", "II", "III", "maVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"])

  if not 't' in ecgA.keys():
    print("Error: ecgA is missing time data!")
    return
  if not 'ecg' in ecgA.keys():
    print("Error: ecgA is missing ecg data!")
    return
  if not 't' in ecgB.keys():
    print("Error: ecgB is missing time data!")
    return
  if not 'ecg' in ecgB.keys():
    print("Error: ecgB is missing ecg data!")
    return

  if np.sum(np.array(ecgA['t']) - np.array(ecgB['t'])) > 1e-3:
    print("Error: time data difference between A and B!")
    return

  for lead in leads:
    if not lead in ecgA['ecg'].keys():
      print("Error: ecgA is missing lead data ({})!".format(lead))
      return
    if not lead in ecgB['ecg'].keys():
      print("Error: ecgB is missing lead data ({})!".format(lead))
      return

  print("Successfully loaded ECG data!")
  #t_min = np.amin(ecg_t)
  #t_max = np.amax(ecg_t)
  #cutoff_frequency = 60 # in Hz
  #filter_order     = 3  # number of sequential filters
  #sampling_rate = int((len(ecg_t) - 1)/(t_max-t_min))

  # PCC and DTW
  if True:
    #ncc_data = np.zeros(len(leads))
    pcc_data = np.zeros(len(leads))
    rmse_data = np.zeros(len(leads))
    dtw_data = np.zeros(len(leads))
    #ecg_rd_pie_t = np.linspace(0.0, 5.0, len(ecgB['ecg']["I"]))

    for it, lead in enumerate(leads):
      lead_data_rd = np.array(ecgA['ecg'][lead], dtype=float)
      lead_data_pie = np.array(ecgB['ecg'][lead], dtype=float)
      fac = np.round(len(lead_data_rd)/len(lead_data_pie))
      valid_idx = (np.arange(len(lead_data_pie))*fac).astype(np.int32)
      lead_data_rd = lead_data_rd[valid_idx]

      #ncc_data[it] = (correlate(lead_data_rd/lead_data_rd.std(), lead_data_pie/lead_data_pie.std(), mode="valid")/len(lead_data_rd))[0]
      pcc_data[it], _ = pearsonr(lead_data_rd, lead_data_pie)
      rmse_data[it] = np.mean((lead_data_rd - lead_data_pie)**2)
      dtw_data[it] = compute_dtw_batched(lead_data_rd, lead_data_pie, window=dtw_window, normalize=False)

    dtw_min = np.amin(dtw_data)
    dtw_max = np.amax(dtw_data)
    dtw_mean = np.mean(dtw_data)
    dtw_std = np.std(dtw_data)
    pcc_min = np.amin(pcc_data)
    pcc_max = np.amax(pcc_data)
    pcc_mean = np.mean(pcc_data)
    pcc_std = np.std(pcc_data)

    #if verbose:
    print("DTW: min {:.3f} max {:.3f} PCC: min {:.3f} max {:.3f}".format(dtw_min, dtw_max, pcc_min, pcc_max))

    # latex table:
    dtw_latex = "DTW"
    pcc_latex = "PCC"

    for it, lead in enumerate(leads):
      if np.abs(dtw_data[it] - dtw_min) < 1e-6 or np.abs(dtw_data[it] - dtw_max) < 1e-6:
        dtw_latex += " & $\\mathbf{{{:.3f}}}$".format(dtw_data[it])
      else:
        dtw_latex += " & ${:.3f}$".format(dtw_data[it])

      if np.abs(pcc_data[it] - pcc_min) < 1e-6 or np.abs(pcc_data[it] - pcc_max) < 1e-6:
        pcc_latex += " & $\\mathbf{{{:.3f}}}$".format(pcc_data[it])
      else:
        pcc_latex += " & ${:.3f}$".format(pcc_data[it])
    
    dtw_latex += " & ${:.3f}\\pm{:.3f}$".format(dtw_mean, dtw_std)
    pcc_latex += " & ${:.3f}\\pm{:.3f}$".format(pcc_mean, pcc_std)

    print(dtw_latex)
    print(pcc_latex)

    # plot
    if False:
      ax = []
      rows, cols = (len(leads), 1)
      #plt.rc('text', usetex=True)
      #plt.rc('font', family='Times New Roman', size=14)
      fig = plt.figure(figsize=(20, 8))
      grd = gs.GridSpec(rows, cols)
      ecg_t = np.array(ecgA["t"])

      for it, lead in enumerate(leads):
        ax.append(fig.add_subplot(grd[it, 0]))
        lead_data_A = np.array(ecgA['ecg'][lead], dtype=float)
        lead_data_B = np.array(ecgB['ecg'][lead], dtype=float)
        peaksA, _ = find_peaks(lead_data_A, height=0.1, distance=dtw_distance)
        peaksB, _ = find_peaks(lead_data_B, height=0.1, distance=dtw_distance)
        ax[it].plot(ecg_t, lead_data_A)
        ax[it].plot(ecg_t, lead_data_B)
        ax[it].plot(ecg_t[peaksA], lead_data_A[peaksA], "x")
        ax[it].plot(ecg_t[peaksB], lead_data_B[peaksB], "x")

      plt.tight_layout()
      plt.show()

  # RR intervals
  if True:
    #RR_lead = 'II'
    #RR_lead = 'I' # 'III'
    #ecg_t = np.array(ecgA["t"])
    #ecgA_RR = np.array(ecgA['ecg'][RR_lead])
    #ecgB_RR = np.array(ecgB['ecg'][RR_lead])
    #print(np.amin(ecgA_RR), np.amax(ecgA_RR), np.amin(ecgB_RR), np.amax(ecgB_RR))
    #RR_mean_A, RR_std_A, peaksA = measure_RR_interval(ecg_t, ecgA_RR)
    #RR_mean_B, RR_std_B, peaksB = measure_RR_interval(ecg_t, ecgB_RR)
    #print("RR")
    #print("A", "mean: {:.1f} std: {:.1f}".format(RR_mean_A, RR_std_A))
    #print("B", "mean: {:.1f} std: {:.1f}".format(RR_mean_B, RR_std_B))
    #print("A", "${:d}\\pm{:d}$".format(int(RR_mean_A), int(RR_std_A)))
    #print("B", "${:d}\\pm{:d}$".format(int(RR_mean_B), int(RR_std_B)))

    RR_leads = leads # ['I']
    ecg_t = np.array(ecgA["t"])
    RR_means_A = np.zeros(len(RR_leads))
    RR_means_B = np.zeros(len(RR_leads))
    RR_latex_A = ""
    RR_latex_B = ""

    for it, RR_lead in enumerate(RR_leads):
      ecgA_RR = np.array(ecgA['ecg'][RR_lead])
      ecgB_RR = np.array(ecgB['ecg'][RR_lead])
      RR_mean_A, RR_std_A, peaksA = measure_RR_interval(ecg_t, ecgA_RR)
      RR_mean_B, RR_std_B, peaksB = measure_RR_interval(ecg_t, ecgB_RR)

      RR_means_A[it] = RR_mean_A
      RR_means_B[it] = RR_mean_B
      RR_latex_A += "${:d}$ & ".format(int(RR_mean_A))
      RR_latex_B += "${:d}$ & ".format(int(RR_mean_B))

    RR_latex_A += "${:d}\\pm{:d}$".format(int(np.mean(RR_means_A)), int(np.std(RR_means_A)))
    RR_latex_B += "${:d}\\pm{:d}$".format(int(np.mean(RR_means_B)), int(np.std(RR_means_B)))

    print(RR_latex_A)
    print(RR_latex_B)

    # plot visualization
    if False:
      ax = []
      rows, cols = (3, 1)
      #plt.rc('text', usetex=True)
      #plt.rc('font', family='Times New Roman', size=14)
      fig = plt.figure(figsize=(20, 8))
      grd = gs.GridSpec(rows, cols)
      
      for row in range(rows):
        for col in range(cols):
          it = row*cols + col
          ax.append(fig.add_subplot(grd[row, col]))
          ax[-1].hlines([0.0], ecg_t[0], ecg_t[-1], ls=":", color="gray")
          ax[-1].set_xlabel("time [s]")

      ax[0].plot(ecg_t, ecgA_RR)
      ax[0].plot(ecg_t, ecgB_RR)
      ax[0].plot(ecg_t[peaksA], ecgA_RR[peaksA], "x")
      ax[0].plot(ecg_t[peaksB], ecgB_RR[peaksB], "x")
      ax[0].set_title("Lead {}".format(RR_lead))
      ax[0].set_ylabel("{} [mV]".format(RR_lead))

      #pdf_filepath = "{}/ecg_QRST.pdf".format(args["outdir"])
      plt.tight_layout()
      #plt.savefig(pdf_filepath, dpi=200, bbox_inches='tight')
      plt.show()

  print("----------------")

  # RR, QT, QRS (experimental)
  if False:
    ecg_t = np.array(ecgA["t"])
    #ecg_aVL = np.array(ecg_rd["ecg"]["aVL"])
    ecg_rd_II = np.array(ecgA["ecg"]["II"])
    ecg_rd_V2 = np.array(ecgA["ecg"]["V2"])
    ecg_rd_V3 = np.array(ecgA["ecg"]["V3"])

    ecg_pie_II = np.array(ecgB["ecg"]["II"])
    ecg_pie_V2 = np.array(ecgB["ecg"]["V2"])
    ecg_pie_V3 = np.array(ecgB["ecg"]["V3"])

    t_min = np.amin(ecg_t)
    t_max = np.amax(ecg_t)
    cutoff_frequency = 60 # in Hz
    filter_order     = 3  # number of sequential filters
    sampling_rate = int((len(ecg_t) - 1)/(t_max-t_min))

    #ecg_II = butter_lowpass_filter(ecg_II, 60, sampling_rate)
    #ecg_II = nk.signal_sanitize(ecg_II)
    #ecg_II = nk.ecg_clean(ecg_II, sampling_rate=sampling_rate, method='neurokit')

    #aVL_peaks, _ = find_peaks(ecg_aVL, height=0.01, distance=50)
    #II_ppeaks, _ = find_peaks(ecg_II, height=0.01, distance=50)
    #II_npeaks, _ = find_peaks(-ecg_II, height=0.01, distance=50)
    #II_peaks, _ = find_peaks(np.abs(ecg_II), height=0.01, distance=30)
    #II_p2peaks, _ = find_peaks(ecg_II+0.2, height=0.01, distance=50)
    #II_p2peaks = II_p2peaks[6:]
    #print(aVL_peaks)
    #V2_peaks, _ = find_peaks(ecg_V2+0.6, height=0.01, distance=60)
    #print(ecg_V2[V2_peaks])
    #
    #print(V2_S)

    # RD
    V2_peaks, _ = find_peaks(ecg_rd_V2, height=0.01, distance=50)
    V2_peaks2, _ = find_peaks(ecg_rd_V2+0.6, height=0.01, distance=60)
    V2_S = V2_peaks2[ecg_rd_V2[V2_peaks2] < 0.0]

    V2_zc = np.where(np.diff(np.sign(ecg_rd_V2)))[0]
    V2_zc = V2_zc[2:]
    V2_zc = V2_zc[::2]

    RR_intervals_rd  = np.diff(ecg_t[V2_peaks])
    QT_intervals_rd  = ecg_t[V2_zc[2:]] - ecg_t[V2_S[:-2]] #np.diff(ecg_t[V2_zc])
    QRS_intervals_rd = np.diff(ecg_t[V2_S])

    #PIE
    V2_peaks_pie, _ = find_peaks(ecg_pie_V2, height=0.01, distance=50)
    V2_peaks_pie = V2_peaks_pie[1:]
    V2_peaks2_pie, _ = find_peaks(ecg_pie_V2+0.6, height=0.01, distance=100)
    V2_S_pie = V2_peaks2_pie[ecg_pie_V2[V2_peaks2_pie] < 0.0]

    V2_zc_pie = np.where(np.diff(np.sign(ecg_pie_V2)))[0]
    V2_zc_pie = V2_zc_pie[1:]
    V2_zc_pie = V2_zc_pie[::2]

    RR_intervals_pie  = np.diff(ecg_t[V2_peaks_pie])
    QT_intervals_pie  = ecg_t[V2_zc_pie[1:]] - ecg_t[V2_S_pie[:-1]] #np.diff(ecg_t[V2_zc])
    QRS_intervals_pie = np.diff(ecg_t[V2_S_pie])

    # post process
    RR_A_dat = RR_intervals_rd*1e3
    RR_A_avg = np.mean(RR_A_dat)
    RR_A_std = np.std(RR_A_dat)

    RR_B_dat = RR_intervals_pie*1e3
    RR_B_avg = np.mean(RR_B_dat)
    RR_B_std = np.std(RR_B_dat)

    QT_A_dat = QT_intervals_rd*1e3
    QT_A_avg = np.mean(QT_A_dat)
    QT_A_std = np.std(QT_A_dat)

    QT_B_dat = QT_intervals_pie*1e3
    QT_B_avg = np.mean(QT_B_dat)
    QT_B_std = np.std(QT_B_dat)

    QRS_A_dat = QRS_intervals_rd*1e3
    QRS_A_avg = np.mean(QRS_A_dat)
    QRS_A_std = np.std(QRS_A_dat)

    QRS_B_dat = QRS_intervals_pie*1e3
    QRS_B_avg = np.mean(QRS_B_dat)
    QRS_B_std = np.std(QRS_B_dat)

    print("RR")
    print("A", "mean: {:.1f} std: {:.1f}".format(RR_A_avg, RR_A_std), RR_A_dat)
    print("B", "mean: {:.1f} std: {:.1f}".format(RR_B_avg, RR_B_std), RR_B_dat)
    print("QT")
    print("A", "mean: {:.1f} std: {:.1f}".format(QT_A_avg, QT_A_std), QT_A_dat)
    print("B", "mean: {:.1f} std: {:.1f}".format(QT_B_avg, QT_B_std), QT_B_dat)
    print("QRS")
    print("A", "mean: {:.1f} std: {:.1f}".format(QRS_A_avg, QRS_A_std), QRS_A_dat)
    print("B", "mean: {:.1f} std: {:.1f}".format(QRS_B_avg, QRS_B_std), QRS_B_dat)

    # latex
    print("A:")
    print("RR  & ${:d}\\pm{:d}$ ".format(int(RR_A_avg), int(RR_A_std)))
    print("QT  & ${:d}\\pm{:d}$ ".format(int(QT_A_avg), int(QT_A_std)))
    print("QRS & ${:d}\\pm{:d}$ ".format(int(QRS_A_avg), int(QRS_A_std)))

    print("B:")
    print("RR  & ${:d}\\pm{:d}$ ".format(int(RR_B_avg), int(RR_B_std)))
    print("QT  & ${:d}\\pm{:d}$ ".format(int(QT_B_avg), int(QT_B_std)))
    print("QRS & ${:d}\\pm{:d}$ ".format(int(QRS_B_avg), int(QRS_B_std)))

    # plot visualization
    if False:
      ax = []
      rows, cols = (3, 1)
      #plt.rc('text', usetex=True)
      #plt.rc('font', family='Times New Roman', size=14)
      fig = plt.figure(figsize=(20, 8))
      grd = gs.GridSpec(rows, cols)
      
      for row in range(rows):
        for col in range(cols):
          it = row*cols + col
          ax.append(fig.add_subplot(grd[row, col]))
          ax[-1].hlines([0.0], ecg_t[0], ecg_t[-1], ls=":", color="gray")
          ax[-1].set_xlabel("time [s]")
      
      #ax[-1].plot(ecg_t, ecg_aVL)
      #ax[-1].plot(ecg_t[pos_peaks], ecg_aVL[pos_peaks], "x")
      #ax[-1].plot(ecg_t, ecg_II)
      #ax[0].plot(ecg_t, ecg_V2)
      #ax[-1].plot(ecg_t[II_ppeaks], ecg_II[II_ppeaks], "x")
      #ax[-1].plot(ecg_t[II_npeaks], ecg_II[II_npeaks], "x")
      #ax[-1].plot(ecg_t[II_peaks], ecg_II[II_peaks], "x")
      #ax[-1].plot(ecg_t[II_p2peaks], ecg_II[II_p2peaks], "x")
      #ax[-1].plot(ecg_t[V2_peaks], ecg_V2[V2_peaks], "x")
      #ax[0].plot(ecg_t[V2_zc], ecg_V2[V2_zc], "x")
      #ax[0].plot(ecg_t[V2_S], ecg_V2[V2_S], "x")

      ax[0].plot(ecg_t, ecg_rd_II)
      ax[0].plot(ecg_t, ecg_pie_II)
      ax[0].set_title("Lead II")
      ax[0].set_ylabel("II [mV]")

      ax[1].plot(ecg_t, ecg_rd_V2)
      ax[1].plot(ecg_t, ecg_pie_V2)
      #ax[1].plot(ecg_t[V2_peaks], ecg_rd_V2[V2_peaks], "x")
      #ax[1].plot(ecg_t[V2_zc], ecg_rd_V2[V2_zc], "x")
      #ax[1].plot(ecg_t[V2_S], ecg_rd_V2[V2_S], "x")
      #ax[1].plot(ecg_t[V2_peaks_pie], ecg_pie_V2[V2_peaks_pie], "x")
      ax[1].plot(ecg_t[V2_zc_pie], ecg_pie_V2[V2_zc_pie], "x")
      ax[1].plot(ecg_t[V2_S_pie], ecg_pie_V2[V2_S_pie], "x")
      ax[1].set_title("Lead V2")
      ax[1].set_ylabel("V2 [mV]")

      ax[2].plot(ecg_t, ecg_rd_V3)
      ax[2].plot(ecg_t, ecg_pie_V3)
      ax[2].set_title("Lead V3")
      ax[2].set_ylabel("V3 [mV]")

      pdf_filepath = "{}/ecg_QRST.pdf".format(args["outdir"])
      plt.tight_layout()
      #plt.savefig(pdf_filepath, dpi=200, bbox_inches='tight')
      plt.show()

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Computes quantiative error metrics between two ECGs.')
  parser.add_argument('-a', '--ecgA',      help='Reference ECG.',          required=True)
  parser.add_argument('-b', '--ecgB',      help='Test ECG.',               required=True)
  parser.add_argument('-o', '--out',       help='Output directory.',       required=True)
  parser.add_argument('-d', '--distance',  help='DTW distance [ms].',      default=200)
  parser.add_argument('-w', '--window',    help='DTW window [ms].',        default=100)
  parser.add_argument('-n', '--normalize', help='DTW normalization flag.', default=False)
  parser.add_argument('-v', '--verbose',   help='Verbosity flag.',         default=False)
  args = vars(parser.parse_args())
  main(args)