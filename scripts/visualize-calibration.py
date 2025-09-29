#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

import carputils.carpio as cuio

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
    #SCV_dir_f = "{}/difcal/sc_{}".format(calib_dir, tissue_f.lower())
    #SCV_dir_t = "{}/difcal/sc_{}".format(calib_dir, tissue_t.lower())
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

# -----------------------------------------------------------------------------
def main(args):
  pln_filepath = args["pln"]

  if not os.path.exists(pln_filepath):
    print("Simulation plan not found!")
    return

  tissue_formats = {"healthy-250um": ["Healthy",    np.array([220.0, 50.0, 32.0, 255.0])/255.0], 
                    "bz-250um":      ["Borderzone", np.array([0.0, 90.0, 181.0, 255.0])/255.0]}

  k_range = np.linspace(-2000.0, 2000.0, 21)
  #vm_rest = [-85.11, -85.6708] # get from sv files

  vlim=[-90, 50]
  #xlim_ap  = [-0.05, 1.5]
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
  plt.rc('text', usetex=True)
  plt.rc('font', family='Times New Roman', size=24)
    
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
    ax[ap_row].set_ylabel("$V_{\\textrm{m}}$ [mV]")
    ax[ap_row].set(xlim=xlim_ap, ylim=ylim_ap)
    ax[ap_row].title.set_text("AP Shape - {}".format(cdata.label))

    # APDR data
    apdr_di, apdr_apd = cdata.fapd
    ax[1].plot(apdr_di, apdr_apd, color=cdata.col, linewidth=2, label=cdata.label)
    #ax[1].vlines([di_crit], ylim_apd[0], apdr_apd[0], ls=":", color=cdata.col)
    ax[1].set(xlim=xlim_apd, ylim=ylim_apd)
    ax[1].set_xlabel("$\\textrm{DI}^{n-1}$ [ms]")
    ax[1].set_ylabel("$\\textrm{APD}$ [ms]")
    ax[1].title.set_text("APD Restitution Curves")

    # CVC data
    v0 = cdata.v0
    Dm = cdata.Dm
    ax[2].plot(k_range, v0[0] - k_range*Dm[0], color=cdata.col, linewidth=2, label=cdata.label)
    ax[2].plot(k_range, v0[1] - k_range*Dm[1], color=cdata.col, linestyle='--', linewidth=2, label=cdata.label)
    ax[2].plot(0.0, v0[0], marker='o', color=cdata.col)
    ax[2].plot(0.0, v0[1], marker='o', color=cdata.col)
    ax[2].set(xlim=xlim_cur, ylim=ylim_cur)
    ax[2].set_xlabel("$\\kappa$ [1/m]")
    ax[2].set_ylabel("$\\textrm{CV}$ [m/s]")
    ax[2].title.set_text("CV - Curvature Dependance")

    # CVR data
    cvr_di, cvr_f, cvr_t = cdata.fcvr
    ax[4].plot(cvr_di, cvr_f, color=cdata.col, linewidth=2, label=cdata.label)
    ax[4].plot(cvr_di, cvr_t, color=cdata.col, linewidth=2, label=cdata.label, linestyle='--')
    #ax[4].vlines([di_crit], ylim_cv[0], cvr_f[0], ls=":", color=cdata.col)
    ax[4].set(xlim=xlim_cv, ylim=ylim_cv)
    ax[4].set_xlabel("$\\textrm{DI}^{n-1}$ [ms]")
    ax[4].set_ylabel("$\\textrm{CV}$ [m/s]")
    ax[4].title.set_text("CV Restitution Curves")

    # SC data
    ax[5].plot(cdata.SCV_f[0][0], cdata.SCV_f[0][1], color=cdata.col, linewidth=2, label=cdata.label)
    ax[5].plot(cdata.SCV_f[1][0], cdata.SCV_f[1][1], color=cdata.col, linestyle=':', linewidth=2, marker='o')
    ax[5].plot(cdata.SCV_t[0][0], cdata.SCV_t[0][1], color=cdata.col, linestyle='--', linewidth=2, label=cdata.label)
    ax[5].plot(cdata.SCV_t[1][0], cdata.SCV_t[1][1], color=cdata.col, linestyle=':', linewidth=2, marker='o')
    ax[5].set_xlim([0.0, 2.0])
    ax[5].hlines(cdata.vm_rest, 0.0, 5.0, ls=":", color="black")
    ax[5].set_xlabel("$\\lambda_{\\zeta}$ [mm]")
    ax[5].set_ylabel("$V_{\\textrm{m}}$ [mV]")
    ax[5].title.set_text("Space Constant")


  #pdf_filepath = "{}/calibration-results-2.pdf".format(calib_dir)
  #print(pdf_filepath)
  plt.tight_layout()
  #fig.subplots_adjust(bottom=0.15)
  #plt.savefig(pdf_filepath, dpi=200, bbox_inches='tight') #, bbox_inches='tight'
  plt.show()
  return







  asd

  # load mesh pts
  #SCV_mesh = "{}/difcal/mesh_0.25_10.00/block".format(calib_dir)
  #mesh_pts = mread.read_points('{}.pts'.format(SCV_mesh))

  # load calibration data
  #for tid, tissue_type in enumerate(tissue_types):


  # visualize
  fig = plt.figure(figsize=(24, 12))
  ax = []
  rows, cols = (2, 4)
  grd = gs.GridSpec(rows, cols) # , width_ratios=[1, 2], height_ratios=[4, 1]
  #colmap = cm.get_cmap('turbo')
  #plt.rcParams.update({'font.size': 18})
  plt.rc('text', usetex=True)
  plt.rc('font', family='Times New Roman', size=24) #, weight='bold'
    
  for row in range(0, rows):
    for col in range(0, cols):
      it = row*cols + col

      if col != 1:
        ax.append(fig.add_subplot(grd[row, col])) # , adjustable='box', aspect=1.0
      else:
        ax.append(fig.add_subplot(grd[row, col], projection='3d'))
      #ax[it].xaxis.set_major_locator(ticker.MaxNLocator(5))
      #ax[it].yaxis.set_major_locator(ticker.MaxNLocator(5))
      
      if col != 1:
        ax[it].grid(True)
  
  for tid, tissue_type in enumerate(tissue_types):
    ap_row = tid*cols
    tissue = "{}-{:d}um".format(tissue_type.lower(), res_um)
    tissue_f = "{}_f".format(tissue_type)
    tissue_t = "{}_t".format(tissue_type)

    AP_filepath = "{}/AP2_{}_{}_bcl_{:d}.0.dat".format(calib_dir, tissue, cell_model, bcl_ms)
    APDr_filepath = "{}/APDr2_{}_{}_bcl_{:d}.0.dat".format(calib_dir, tissue, cell_model, bcl_ms)
    CVr_filepath = "{}/CVr3_{}_bcl_{:d}.0.dat".format(calib_dir, tissue, bcl_ms)

    ap_t, ap_di, ap_vm = load_AP2D(AP_filepath)
    apdr_di, apdr_apd = load_APDr(APDr_filepath)
    cvr_di, cvr_f, cvr_t = load_CVr(CVr_filepath)
    SCV_dir_f = "{}/difcal/sc_{}".format(calib_dir, tissue_f.lower())
    SCV_dir_t = "{}/difcal/sc_{}".format(calib_dir, tissue_t.lower())
    SCV_f = load_SCV(mesh_pts, SCV_dir_f, vm_rest[tid])
    SCV_t = load_SCV(mesh_pts, SCV_dir_t, vm_rest[tid])

    gdata = funcs[tissue]["CV"]["G"]
    cvdata = funcs[tissue]["CV"]["CVref"]
    epdata = funcs[tissue]["EP"]["space-constant"]
    v0 = np.array([cvdata["vf"], cvdata["vs"], cvdata["vn"]])
    ge = np.array([gdata["gel"], gdata["get"], gdata["gen"]])
    gi = np.array([gdata["gil"], gdata["git"], gdata["gin"]])
    sc = np.array([epdata["lf"], epdata["ls"], epdata["ln"]])
    beta = gdata["surf2vol"]*1e6
    gm = harmonic_mean(ge, gi)
    Dm = gm/(beta*C_m)

    # use same DI_crit
    di_crit = np.amin(cvr_di)
    valid_range = np.where(apdr_di >=  di_crit)
    apdr_di = (apdr_di[valid_range])[:-1]
    apdr_apd = (apdr_apd[valid_range])[:-1]

    ap_t[-1] = xlim_ap[1]/0.9
    ap_vm[0] = vm_rest[tid]
    alpha_step = 1.0/len(ap_di)
    alpha = 1.0

    for it in range(len(ap_di)):
      #ax[ap_row].plot(ap_t, ap_vm[:, it]) #ax[0].plot(ap_t[:400], ap_vm[:400, it])
      #ax[ap_row].plot(ap_t*apdr_apd[it], ap_vm[:, it])
      ax[ap_row].plot(ap_t*0.9, ap_vm[:, it], color=colors[tid], alpha=alpha)
      alpha -= alpha_step/1.2

    ax[ap_row].set_xlabel("Time [ms]")
    ax[ap_row].set_ylabel("$f_{\\textrm{AP}}$ [mV]")
    ax[ap_row].set(xlim=xlim_ap, ylim=ylim_ap)
    ax[ap_row].title.set_text("AP Shape - {}".format(tissue_type))

    # traces AP
    #ap_aspect = (xlim_ap[1] - xlim_ap[0])/(ap_di[0] - ap_di[-1])
    #ax[ap_row+1].imshow(ap_vm.T, aspect=ap_aspect, extent=[xlim_ap[0], xlim_ap[1], ap_di[-1], ap_di[0]], vmin=vlim[0], vmax=vlim[1], cmap='coolwarm', interpolation="gaussian") #cmap='turbo' aspect=(len(ap_t)/len(ap_di))
    #ax[ap_row+1].set_xlabel("Time [ms]")
    ##ax[ap_row+1].set_ylabel("$\\textrm{DI}^{n-1}$ [ms]")
    #ax[ap_row+1].set_ylabel("$\\textrm{APD}^{n}$ [ms]")
    #ax[ap_row+1].title.set_text("AP Density - {}".format(tissue_type))

    # 3D AP plot
    t_ap, apd_ap = np.meshgrid(ap_t, ap_di)
    print(t_ap.shape, apd_ap.shape, ap_vm.T.shape)
    ls = LightSource(300, 50) # To use a custom hillshading mode, override the built-in shading and pass in the rgb colors of the shaded surface calculated from "shade".
    #col_healthy = ls.shade(ap_vm.T, cmap=cm.coolwarm, vert_exag=0.1, blend_mode='soft')
    #ax[ap_row+1].set(xlim=xlim_ap, ylim=ylim_ap, zlim=[-90.0, 40.0])
    ax[ap_row+1].plot_surface(t_ap, apd_ap, ap_vm.T, linewidth=0, cmap='coolwarm', alpha=None, rstride=100, vmin=-90.0, vmax=40.0) #rstride=1, cstride=100, antialiased=True
    #ax[ap_row+1].plot_surface(t_ap, apd_ap, ap_vm.T, linewidth=0, facecolors=col_healthy, shade=False)

    ax[2].plot(apdr_di, apdr_apd, color=colors[tid], linewidth=2, label=tissue_type)
    ax[2].vlines([di_crit], ylim_apd[0], apdr_apd[0], ls=":", color=colors[tid])

    ax[3].plot(t_cvk, v0[0] - t_cvk*Dm[0], color=colors[tid], linewidth=2, label=tissue_f)
    ax[3].plot(t_cvk, v0[1] - t_cvk*Dm[1], color=colors[tid], linestyle='--', linewidth=2, label=tissue_t)
    ax[3].plot(0.0, v0[0], marker='o', color=colors[tid])
    ax[3].plot(0.0, v0[1], marker='o', color=colors[tid])

    ax[6].plot(cvr_di,     cvr_f, color=colors[tid], linewidth=2, label=tissue_f)
    ax[6].plot(cvr_di,     cvr_t, color=colors[tid], linewidth=2, label=tissue_t, linestyle='--')
    ax[6].vlines([di_crit], ylim_cv[0], cvr_f[0], ls=":", color=colors[tid])

    ax[7].plot(SCV_f[0][0], SCV_f[0][1], color=colors[tid], linewidth=2, label=tissue_f)
    ax[7].plot(SCV_f[1][0], SCV_f[1][1], color=colors[tid], linestyle=':', linewidth=2, marker='o')
    ax[7].plot(SCV_t[0][0], SCV_t[0][1], color=colors[tid], linestyle='--', linewidth=2, label=tissue_t)
    ax[7].plot(SCV_t[1][0], SCV_t[1][1], color=colors[tid], linestyle=':', linewidth=2, marker='o')

    

    
  
  '''
  for tid, tissue_type in enumerate(tissue_types):
    tissue = "{}-{:d}um".format(tissue_type.lower(), res_um)
    tissue_f = "{}_f".format(tissue_type)
    tissue_t = "{}_t".format(tissue_type)
    AP_filepath = "{}/AP_{}_{}_bcl_{:d}.0.dat".format(calib_dir, tissue, cell_model, bcl_ms)
    APDr_filepath = "{}/APDr_{}_{}_bcl_{:d}.0.dat".format(calib_dir, tissue, cell_model, bcl_ms)
    CVr_filepath = "{}/CVr_{}_bcl_{:d}.0.dat".format(calib_dir, tissue, bcl_ms)
    SCV_dir_f = "{}/difcal/sc_{}".format(calib_dir, tissue_f.lower())
    SCV_dir_t = "{}/difcal/sc_{}".format(calib_dir, tissue_t.lower())
    SCV_f = load_SCV(mesh_pts, SCV_dir_f, vm_rest[tid])
    SCV_t = load_SCV(mesh_pts, SCV_dir_t, vm_rest[tid])

    ap_t, ap_di, ap_vm = load_AP2D(AP_filepath)
    apdr_di, apdr_apd = load_APDr(APDr_filepath)
    cvr_di, cvr_f, cvr_t = load_CVr(CVr_filepath)
    ap_vm[0] = vm_rest[tid] # set first data points to Vm_rest for visualization

    gdata = funcs[tissue]["CV"]["G"]
    cvdata = funcs[tissue]["CV"]["CVref"]
    epdata = funcs[tissue]["EP"]["space-constant"]
    v0 = np.array([cvdata["vf"], cvdata["vs"], cvdata["vn"]])
    ge = np.array([gdata["gel"], gdata["get"], gdata["gen"]])
    gi = np.array([gdata["gil"], gdata["git"], gdata["gin"]])
    sc = np.array([epdata["lf"], epdata["ls"], epdata["ln"]])
    beta = gdata["surf2vol"]*1e6
    gm = harmonic_mean(ge, gi)
    Dm = gm/(beta*Cm)


    for col, apd90_ms in enumerate(apdr_apd):
      ap_vm[:,col] = np.interp(ap_t*apdr_apd[-1], ap_t*apd90_ms, ap_vm[:,col])

    ap_vm = ap_vm[:,::-1]
    #print(ap_vm[:,::-1].shape)
    #asdas

    ap_row = tid*cols
    #ap_t_lbl = [str(xl) for xl in ap_t]

    for it in range(len(ap_di)):
      #ax[ap_row].plot(ap_t, ap_vm[:, it]) #ax[0].plot(ap_t[:400], ap_vm[:400, it])
      ax[ap_row].plot(ap_t*apdr_apd[it], ap_vm[:, it])

    ax[ap_row].set_xlabel("Time [ms]")
    ax[ap_row].set_ylabel("Transmembrane Voltage [mV]")
    ax[ap_row].set_ylim(vlim)
    ax[ap_row].title.set_text("$V_\\textrm{m}$ Shape - Healthy")

    #ax[1].imshow(ap_vm, aspect=4.0, extent=[np.amin(ap_t), np.amax(ap_t), np.amin(ap_di), np.amax(ap_di)], vmin=vlim[0], vmax=vlim[1], cmap='turbo')
    ax[ap_row+1].imshow(ap_vm.T, aspect=(len(ap_t)/len(ap_di)), vmin=vlim[0], vmax=vlim[1], cmap='coolwarm') #cmap='turbo'
    #extent=[np.amin(ap_t), np.amax(ap_t), np.amin(ap_di), np.amax(ap_di)]
    #ax[ap_row+1].imshow(ap_vm.T, aspect=(len(ap_t)/len(ap_di)), extent=extent, vmin=vlim[0], vmax=vlim[1], cmap='turbo')
    #ax[ap_row+1].set_xticklabels(ap_t_lbl)
    ax[ap_row+1].set_xlabel("Time [ms]")
    ax[ap_row+1].set_ylabel("DI [ms]")
    ax[ap_row+1].title.set_text("Vm Density - Healthy")

    #Nxticks = 5
    #Nxstep = len(ap_t)/Nxticks

    #xticks = [ap_t[int(it*Nxstep)]* for it in range(Nxticks)]
    #xticks = [ for tick in range(0, )]
    #len(ap_t)
    #ax[ap_row+1].set_xticks(xticks)
    #ax[ap_row+1].set_xticklabels([0,1,2])

    yticks = ax[ap_row+1].get_yticks().tolist()
    #print(yticks)
    #asdas
    Nyticks = len(yticks)
    Nystep = len(ap_di)/Nyticks
    yticks = [ap_di[int(it*Nystep)] for it in range(Nyticks)]
    print(yticks)
    #ax[ap_row+1].set_yticks(yticks)
    ax[ap_row+1].set_yticklabels(yticks)
    #ax[ap_row+1].set_yticklabels([0,1,2])
    #ticks_loc = [ap_t[np.clip(int(xl), 0, 500-1)] for xl in ax[ap_row+1].get_xticks().tolist()]
    #print(ticks_loc)
    #ax[ap_row+1].set_xticklabels(ticks_loc)
    #sdas

    ax[2].plot(apdr_di, apdr_apd, color=colors[tid], linewidth=2, label=tissue_type)
    ax[6].plot(cvr_di, cvr_f, color=colors[tid], linewidth=2, label=tissue_f)
    ax[6].plot(cvr_di, cvr_t, color=colors[tid], linestyle='--', linewidth=2, label=tissue_t)

    ax[3].plot(t_cvk, v0[0] - t_cvk*Dm[0], color=colors[tid], linewidth=2, label=tissue_f)
    ax[3].plot(t_cvk, v0[1] - t_cvk*Dm[1], color=colors[tid], linestyle='--', linewidth=2, label=tissue_t)

    ax[7].plot(SCV_f[0][0], SCV_f[0][1], color=colors[tid], linewidth=2, label=tissue_f)
    ax[7].plot(SCV_f[1][0], SCV_f[1][1], color=colors[tid], linestyle=':', linewidth=2, marker='o')
    ax[7].plot(SCV_t[0][0], SCV_t[0][1], color=colors[tid], linestyle='--', linewidth=2, label=tissue_t)
    ax[7].plot(SCV_t[1][0], SCV_t[1][1], color=colors[tid], linestyle=':', linewidth=2, marker='o')

  
  #ax[2].plot(apd_borderz[0], apd_borderz[1], color=col_borderz, linewidth=2, label="borderzone")
  ax[2].set_xlim([0, 300])
  ax[2].set_ylim([0, 400])
  ax[2].set_xlabel("DI [ms]")
  ax[2].set_ylabel("APD [ms]")
  ax[2].title.set_text("APD Restitution Curves")
  #ax[2].legend()

  #ax[3].plot(cvr_di, cvr_f, color=col_healthy, linewidth=2, label="{}_f".format(tissue_type))
  #ax[3].plot(cvr_di, cvr_t, color=col_healthy, linestyle='--', linewidth=2, label="{}_t".format(tissue_type))
  #ax[3].plot(cv_borderz[0][0], cv_borderz[0][1], color=col_borderz, linewidth=2, label="borderzone_f")
  #ax[3].plot(cv_borderz[1][0], cv_borderz[1][1], color=col_borderz, linestyle='--', linewidth=2, label="borderzone_t")
  #ax[3].set_xlim([200, 600])
  ax[6].set_ylim([0, 1.2])
  ax[6].set_xlabel("DI [ms]")
  ax[6].set_ylabel("CV [m/s]")
  ax[6].title.set_text("CV Restitution Curves")
  #ax[6].legend()

  ax[3].set_xlim([-2000, 2000])
  ax[3].set_ylim([0, 1.5])
  ax[3].vlines([0], 0.0, 1.5, ls=":", color="black")
  ax[3].set_xlabel("k [1/m]")
  ax[3].set_ylabel("CV [m/s]")
  ax[3].title.set_text("CV - Curvature Dependance")
  #ax[3].legend()

  ax[7].set_xlim([0.0, 2.0])
  ax[7].hlines(vm_rest, 0.0, 5.0, ls=":", color="black")
  ax[7].set_xlabel("Distance from Stimulus [mm]")
  ax[7].set_ylabel("Transmembrane Voltage [mV]")
  ax[7].title.set_text("Space Constant")
  #ax[7].legend()

  lines_labels = [axi.get_legend_handles_labels() for axi in [ax[6]]]
  lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
  fig.legend(lines, labels, loc='lower center', ncol=4)
  #fig.legend(lines, labels, loc=(0, -0.1), ncol=4)
  '''

  ax[2].set(xlim=xlim_apd, ylim=ylim_apd)
  ax[2].set_xlabel("$\\textrm{DI}^{n-1}$ [ms]")
  ax[2].set_ylabel("$f_{\\textrm{APD}}$ [ms]")
  ax[2].title.set_text("APD Restitution Curves")

  ax[3].set(xlim=xlim_cur, ylim=ylim_cur)
  ax[3].set_xlabel("$\\kappa$ [1/m]")
  ax[3].set_ylabel("$f_{\\textrm{CVC}}$ [m/s]")
  ax[3].title.set_text("CV - Curvature Dependance")

  ax[6].set(xlim=xlim_cv, ylim=ylim_cv)
  ax[6].set_xlabel("$\\textrm{DI}^{n-1}$ [ms]")
  ax[6].set_ylabel("$f_{\\textrm{CVR}}$ [m/s]")
  ax[6].title.set_text("CV Restitution Curves")

  ax[7].set_xlim([0.0, 2.0])
  ax[7].hlines(vm_rest, 0.0, 5.0, ls=":", color="black")
  ax[7].set_xlabel("$\\lambda_{\\zeta}$ [mm]")
  ax[7].set_ylabel("$V_{\\textrm{m}}$ [mV]")
  ax[7].title.set_text("Space Constant")

  #pdf_filepath = "{}/calibration-results.pdf".format(calib_dir)
  pdf_filepath = "{}/calibration-results-2.pdf".format(calib_dir)
  #print(pdf_filepath)
  plt.tight_layout()
  #fig.subplots_adjust(bottom=0.15)
  plt.savefig(pdf_filepath, dpi=200, bbox_inches='tight') #, bbox_inches='tight'
  plt.show()

# -----------------------------------------------------------------------------
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Script to visualize calibration results.')
  parser.add_argument('--pln', help='Input: simulation plan', required=True)
  #parser.add_argument('--odir', help='Output: pdf directory',        required=True)
  args = vars(parser.parse_args())
  main(args)
