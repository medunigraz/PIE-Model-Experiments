#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import argparse
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import carputils.carpio as cuio

# _________________________________________________________________________________________________
class fnc_data:
  def __init__(self, label, col, fnc, pln_dir, fnc_id):
    self.label = label
    self.col = col
    self.pln_dir = pln_dir + "/"

    gdata  = fnc["CV"]["G"]
    cvdata = fnc["CV"]["CVref"]
    scdata = fnc["EP"]["space-constant"]

    self.v0 = np.array([cvdata["vf"], cvdata["vs"], cvdata["vn"]])
    self.ge = np.array([gdata["gel"], gdata["get"], gdata["gen"]])
    self.gi = np.array([gdata["gil"], gdata["git"], gdata["gin"]])
    self.sc = np.array([scdata["lf"], scdata["ls"], scdata["ln"]])
    self.beta = gdata["surf2vol"]*1e6
    self.C_m = 0.01
    self.gm = (self.ge*self.gi)/(self.ge+self.gi)
    self.Dm = self.gm/(self.beta*self.C_m)
    
    self.fapd = self.load_APDr(self.pln_dir + fnc["EP"]["APDres-file"])
    self.fcvr = self.load_CVr(self.pln_dir + fnc["CV"]["CVres-file"])
    self.fap  = self.load_AP2D(self.pln_dir + fnc["EP"]["AP-file"])

    self.vm_rest = self.load_IMP(self.pln_dir + fnc["EP"]["init"])
    scm_dir_l = self.pln_dir + "tissue_init/" + fnc_id + "_l"
    scm_dir_t = self.pln_dir + "tissue_init/" + fnc_id + "_t"
    mesh_pts = np.loadtxt(scm_dir_l + "/mesh/cable.pts", skiprows=1, dtype=np.float32) #um
    self.SCV_f = self.load_SCV(mesh_pts, scm_dir_l+"/scm/", self.vm_rest)
    self.SCV_t = self.load_SCV(mesh_pts, scm_dir_t+"/scm/", self.vm_rest)

  def load_AP2D(self, ap2d_path):
    ap2d_data = np.loadtxt(ap2d_path)
    tlen = np.int32(ap2d_data[0,0])
    ap_t = ap2d_data[1:,0]
    ap_di = ap2d_data[0,1:]
    ap_vm = ap2d_data[1:,1:]

    if tlen != len(ap_t):
      print("load_AP2D: tlen mismatch!")

    if tlen*len(ap_di) != np.prod(ap_vm.shape):
      print("load_AP2D: size mismatch!")

    return ap_t, ap_di, ap_vm

  def load_CVr(self, CVr_filepath):
    return np.loadtxt(CVr_filepath).T

  def load_APDr(self, APDr_filepath):
    return np.loadtxt(APDr_filepath).T

  def load_IMP(self, imp_path):
    vm_rest = -90.0

    with open(imp_path) as imp_file:
      vm_rest = float(imp_file.readline().split("#")[0])

    return vm_rest

  def load_SCV(self, mesh_pts, path, Vm_rest_mV):
    mesh_dx_mm = 0.25
    UM2MM = 1.0e-3
    vm_path = "{}/vm.igb".format(path)

    if not os.path.exists(vm_path):
      raise RuntimeError("load_SCV(): vm file not found: {}".format(vm_path))

    vm_igb = cuio.igb.IGBFile(vm_path)
    vm_hdr = vm_igb.header()
    vm_dat = vm_igb.data().reshape((vm_hdr['t'], -1))

    num_time_slices = vm_hdr['t']
    delta_time = vm_hdr['inc_t']

    mesh_xpos = np.unique(mesh_pts[:,0])*UM2MM
    pts_idx = list(list() for _ in range(mesh_xpos.size))
    tolerance = mesh_dx_mm*0.25

    for i, pt in enumerate(mesh_pts):
      xidx = np.where(np.isclose(mesh_xpos, pt[0]*UM2MM, atol=tolerance, rtol=0.0))[0][0]
      pts_idx[xidx].append(i)

    avg_data = np.zeros((num_time_slices, len(pts_idx)), dtype=float)

    for t in range(num_time_slices):
      for i, idx in enumerate(pts_idx):
        data = sum(vm_dat[t][j] for j in idx)
        avg_data[t][i] = data/len(idx)

    mesh_xpos += 10.0
    k = (avg_data[-1][1]-avg_data[-1][0])/(mesh_xpos[1]-mesh_xpos[0])
    ct = []
    cf = []

    if abs(k) > 0.0:
      d = avg_data[-1][0] - k*mesh_xpos[0]
      x_star = (Vm_rest_mV-d)/k
      dx_star = abs(x_star-mesh_xpos[0])

      ct = [mesh_xpos[0], x_star]
      cf = [k*mesh_xpos[0]+d, k*x_star+d]

    #print(ct, cf)
    tdat = mesh_xpos
    fdat = avg_data[i]

    return [tdat, fdat], [ct, cf]

  def harmonic_mean(ge, gi):
    return (ge*gi)/(ge+gi)

# _________________________________________________________________________________________________
def main(args):
  pln_filepath = args["pln"]

  if not os.path.exists(pln_filepath):
    print("Simulation plan not found!")
    return

  tissue_formats = {"healthy-250um": ["Healthy",    np.array([220.0, 50.0, 32.0, 255.0])/255.0], 
                    "bz-250um":      ["Borderzone", np.array([0.0, 90.0, 181.0, 255.0])/255.0]}

  k_range = np.linspace(-2000.0, 2000.0, 21)

  vlim=[-90, 50]
  xlim_ap  = [-10.0, 500.0]
  ylim_ap  = [-90, 50]
  xlim_apd = [0, 310]
  ylim_apd = [0, 380]
  xlim_cv  = xlim_apd
  ylim_cv  = [0, 0.9]
  xlim_cur = [-2000, 2000]
  ylim_cur = [0, 1.25]

  # obtain calibration data from plan
  calib_data = []
  plan_data = []

  with open(pln_filepath) as plan_file:
    plan_data = json.load(plan_file)

  for key, fnc in plan_data["functions"].items():
    label, col = tissue_formats[key]
    calib_data.append(fnc_data(label, col, fnc, os.path.dirname(os.path.realpath(pln_filepath)), key))

  # plot calibration data
  fig = plt.figure(figsize=(24, 12))
  ax = []
  rows, cols = (2, 3)
  grd = gs.GridSpec(rows, cols)
  #plt.rc('text', usetex=True)
  #plt.rc('font', family='Times New Roman', size=24)
    
  for row in range(0, rows):
    for col in range(0, cols):
      it = row*cols + col
      ax.append(fig.add_subplot(grd[row, col]))
      ax[it].grid(True)
        
  for cid, cdata in enumerate(calib_data[::-1]):
    ap_row = cid*cols
    ap_t, ap_apd, ap_vm = cdata.fap

    ap_t[-1] = xlim_ap[1]/0.9
    ap_vm[0] = cdata.vm_rest
    alpha_step = 1.0/len(ap_apd)
    alpha = 1.0

    for it in range(len(ap_apd)):
      ax[ap_row].plot(ap_t*0.9, ap_vm[:, it], color=cdata.col, alpha=alpha)
      alpha -= alpha_step/1.2

    ax[ap_row].set_xlabel("Time [ms]")
    ax[ap_row].set_ylabel("V_m [mV]") #$V_{\\textrm{m}}$ [mV]
    ax[ap_row].set(xlim=xlim_ap, ylim=ylim_ap)
    ax[ap_row].title.set_text("AP Shape - {}".format(cdata.label))

    # APDR data
    apdr_di, apdr_apd = cdata.fapd
    ax[1].plot(apdr_di, apdr_apd, color=cdata.col, linewidth=2, label=cdata.label)
    #ax[1].vlines([di_crit], ylim_apd[0], apdr_apd[0], ls=":", color=cdata.col)
    ax[1].set(xlim=xlim_apd, ylim=ylim_apd)
    ax[1].set_xlabel("$DI^(n-1) [ms]") #$\\textrm{DI}^{n-1}$ [ms]
    ax[1].set_ylabel("APD [ms]") #$\\textrm{APD}$ [ms]
    ax[1].title.set_text("APD Restitution Curves")

    # CVC data
    v0 = cdata.v0
    Dm = cdata.Dm
    ax[2].plot(k_range, v0[0] - k_range*Dm[0], color=cdata.col, linewidth=2, label=cdata.label)
    ax[2].plot(k_range, v0[1] - k_range*Dm[1], color=cdata.col, linestyle='--', linewidth=2, label=cdata.label)
    ax[2].plot(0.0, v0[0], marker='o', color=cdata.col)
    ax[2].plot(0.0, v0[1], marker='o', color=cdata.col)
    ax[2].set(xlim=xlim_cur, ylim=ylim_cur)
    ax[2].set_xlabel("kappa [1/m]") #$\\kappa$ [1/m]
    ax[2].set_ylabel("CV [m/s]") #$\\textrm{CV}$ [m/s]
    ax[2].title.set_text("CV - Curvature Dependance")

    # CVR data
    cvr_di, cvr_f, cvr_t = cdata.fcvr
    ax[4].plot(cvr_di, cvr_f, color=cdata.col, linewidth=2, label=cdata.label)
    ax[4].plot(cvr_di, cvr_t, color=cdata.col, linewidth=2, label=cdata.label, linestyle='--')
    #ax[4].vlines([di_crit], ylim_cv[0], cvr_f[0], ls=":", color=cdata.col)
    ax[4].set(xlim=xlim_cv, ylim=ylim_cv)
    ax[4].set_xlabel("DI^(n-1) [ms]") #$\\textrm{DI}^{n-1}$ [ms]
    ax[4].set_ylabel("CV [m/s]") #$\\textrm{CV}$ [m/s]
    ax[4].title.set_text("CV Restitution Curves")

    # SC data
    ax[5].plot(cdata.SCV_f[0][0], cdata.SCV_f[0][1], color=cdata.col, linewidth=2, label=cdata.label)
    ax[5].plot(cdata.SCV_f[1][0], cdata.SCV_f[1][1], color=cdata.col, linestyle=':', linewidth=2, marker='o')
    ax[5].plot(cdata.SCV_t[0][0], cdata.SCV_t[0][1], color=cdata.col, linestyle='--', linewidth=2, label=cdata.label)
    ax[5].plot(cdata.SCV_t[1][0], cdata.SCV_t[1][1], color=cdata.col, linestyle=':', linewidth=2, marker='o')
    ax[5].set_xlim([0.0, 2.0])
    ax[5].hlines(cdata.vm_rest, 0.0, 5.0, ls=":", color="black")
    ax[5].set_xlabel("lambda_zeta [mm]") #$\\lambda_{\\zeta}$ [mm]
    ax[5].set_ylabel("Vm [mV]") #$V_{\\textrm{m}}$ [mV]
    ax[5].title.set_text("Space Constant")

  plt.tight_layout()
  pdf_filepath = "{}/calibration-results.pdf".format(args["odir"])
  plt.savefig(pdf_filepath, dpi=200, bbox_inches='tight')
  plt.show()
  return

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Script to visualize PIE calibration results.')
  parser.add_argument('--pln',  help='simulation plan',  required=True)
  parser.add_argument('--odir', help='output directory', required=True)
  args = vars(parser.parse_args())
  main(args)
