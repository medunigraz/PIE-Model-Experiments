#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import warnings
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import imageio
from skimage import measure
import carputils.carpio as cuio
from scipy.spatial import distance
import matplotlib.animation as animation
from scipy.ndimage import gaussian_filter

def colored_line_between_pts(x, y, c, ax, **lc_kwargs):
    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Check color array size (LineCollection still works, but values are unused)
    if len(c) != len(x) - 1:
        warnings.warn(
            "The c argument should have a length one less than the length of x and y. "
            "If it has the same length, use the colored_line function instead."
        )

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, **lc_kwargs)
    lc.set_array(c)

    return ax.add_collection(lc)

# -----------------------------------------------------------------------------
def main(args):
  simtype = "MD" # MD REK

  base_path = "/home/tom/workspace/eikonal-experiments/results/sim/6_rotor"

  if simtype == "MD":
    vm_path = "{}/rd_250um/sim/vm.igb".format(base_path)
  elif simtype == "REK":
    vm_path = "{}/ek_1000um/vm_mv.igb".format(base_path)
  else:
    print("should not happen!")
    return

  if not os.path.exists(vm_path):
    raise RuntimeError("load_SCV(): vm file not found: {}".format(vm_path))

  vm_igb = cuio.igb.IGBFile(vm_path)
  vm_hdr = vm_igb.header()
  vm_dat = vm_igb.data().reshape((vm_hdr['t'], -1))
  Npts = int(vm_dat.shape[1]/2)
  Mpts = int(np.sqrt(Npts))
  min_len = 20

  t_start = 340
  t_end = 5000
  t_step = 2
  dt = t_step*1e-3

  fig, ax = plt.subplots()
  #last_contour = []
  trajectory = np.zeros((t_end-t_start, 2))
  vm_slices = np.zeros((t_end-t_start, Mpts, Mpts))
  contour_vm_list = []
  contour_gvm_list = []
  max_vmc = -1
  max_vmgc = -1

  for t in range(t_start, t_end, t_step):
    vm_slice = vm_dat[t,:Npts].reshape((Mpts, Mpts))
    #vm_slice = gaussian_filter(vm_slice, sigma=1)
    t_slice = int((t-t_start)/t_step)
    vm_slices[t_slice] = vm_slice
    
    measured_contours_vm = measure.find_contours(vm_slice, -75.0)
    contours_vm = [np.array(item) for item in measured_contours_vm if len(item) > min_len]
    contour_vm = np.vstack(contours_vm)

    if t > t_start:
      grad_vm = (vm_slices[t_slice-1] - vm_slices[t_slice])/dt
      measured_contours_gvm = measure.find_contours(grad_vm, 0.0)
      contours_gvm = [np.array(item) for item in measured_contours_gvm if len(item) > min_len]
      contour_gvm = np.vstack(contours_gvm)

      cdists = distance.cdist(contour_vm, contour_gvm).min(axis=1)
      sdists = np.linalg.norm(trajectory[t_slice-1] - contour_vm, axis=1)

      if t > t_start + 3*t_step:
        tdists = 0.9*cdists+0.1*sdists
        #tdists = 0.8*cdists+0.2*sdists
        #tdists = 0.7*cdists+0.3*sdists
        #tdists = 0.6*cdists+0.4*sdists
      else:
        tdists = cdists

      idx = np.argmin(tdists)
      trajectory[t_slice] = contour_vm[idx]
      contour_gvm_list.append(contours_gvm)
    else:
      contour_gvm_list.append(contours_vm)
      trajectory[t_slice] = [0.0, 0.0]
    
    contour_vm_list.append(contours_vm)

    max_vmc = max(max_vmc, len(contour_vm_list[-1]))
    max_vmgc = max(max_vmgc, len(contour_gvm_list[-1]))

  #im = ax.imshow(vm_slices[0], vmin=-90.0, vmax=40.0, cmap='coolwarm', interpolation="gaussian")
  im = ax.imshow(np.ones(vm_slices[0].shape)*(0.0), cmap='Greys')
  print(max_vmc, max_vmgc)
  pvm = []
  pvmg = []

  for it in range(max_vmc):
    p1, = ax.plot(contour_vm_list[0][0][:,1], contour_vm_list[0][0][:,0], linewidth=2, color='blue')
    pvm.append(p1)

  #for it in range(max_vmgc):
  #  p2, = ax.plot(contour_vm_list[0][0][:,1], contour_vm_list[0][0][:,0], linewidth=2, color='red')
  #  pvmg.append(p2)

  #p1, = ax.plot(contour_vm_list[0][:,1], contour_vm_list[0][:,0], linewidth=2, color='blue')
  #p2, = ax.plot(contour_gvm_list[0][:,1], contour_gvm_list[0][:,0], linewidth=2, color='red')
  p3, = ax.plot(trajectory[:,1], trajectory[:,0], linewidth=2, color='black')
  p4, = ax.plot([trajectory[0,1]], [trajectory[0,0]], marker='o', color='yellow', markeredgecolor='black') 

  #dydx = (trajectory[1:-1,0] - trajectory[2:,0])/(trajectory[1:-1,1] - trajectory[2:,1])
  #col = np.arange(0, len(trajectory[1:,1]-1), 1)/len(trajectory[1:,1]-1)
  #nline = len(trajectory[1:,1]-1)
  #ncol = 5#max(nline-10, nline)
  #col = np.zeros((nline-1,))
  #col[nline-ncol-1:] = np.arange(0, ncol, 1)/ncol
  #line = colored_line_between_pts(trajectory[1:,1], trajectory[1:,0], col, ax, linewidth=2, cmap="Greys")
  #plt.show()
  #asdas

  ax.set_xticks([])
  ax.set_yticks([])
  plt.gca().invert_yaxis()
  #plt.show()

  def update(frame):
    #im.set_array(vm_slices[frame])
    contours_vm = contour_vm_list[frame]
    contours_gvm = contour_gvm_list[frame]

    for it in range(max_vmc):
      if it < len(contours_vm):
        pvm[it].set_xdata(contours_vm[it][:,1])
        pvm[it].set_ydata(contours_vm[it][:,0])
      else:
        pvm[it].set_xdata([])
        pvm[it].set_ydata([])

    #for it in range(max_vmgc):
    #  if it < len(contours_gvm):
    #    pvmg[it].set_xdata(contours_gvm[it][:,1])
    #    pvmg[it].set_ydata(contours_gvm[it][:,0])
    #  else:
    #    pvmg[it].set_xdata([])
    #    pvmg[it].set_ydata([])

    p3.set_xdata(trajectory[1:frame,1])
    p3.set_ydata(trajectory[1:frame,0])

    p4.set_xdata([trajectory[frame,1]])
    p4.set_ydata([trajectory[frame,0]])

    #plt.savefig("scripts/trajectories_{}/traj_{}_{:d}.png".format(simtype, simtype, t_step*frame + 340), bbox_inches='tight')
    plt.savefig("png2/trajectories_{}/traj_{}_{:04d}.png".format(simtype, simtype, t_step*frame + 340), bbox_inches='tight', dpi=200)#, dpi=200
    
    return (pvm, pvmg, p3, p4, im)

  num_frames = int((t_end-t_start)/t_step)
  ani = animation.FuncAnimation(fig=fig, func=update, frames=num_frames, interval=100, repeat=False)#1/10
  plt.show()
  #mplplot.tight_layout()
  #writer = animation.FFMpegWriter(fps=10, metadata=dict(artist='Thomas Schrotter'), bitrate=1800)
  #ani.save("trajectory-{}-vis.mp4".format(simtype), writer=writer)
  #ani.save(filename="scripts/trajectories/test.png", writer="pillow")

# -----------------------------------------------------------------------------
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Script to visualize calibration results.')
  parser.add_argument('--idir', help='Input: calibration directory', required=False)
  #parser.add_argument('--odir', help='Output: pdf directory',        required=True)
  args = vars(parser.parse_args())
  main(args)