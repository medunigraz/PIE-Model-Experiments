#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

import numpy as np
import matplotlib.pyplot as plt
import carputils.carpio as cuio
import matplotlib.animation as animation

from skimage import measure
from scipy.spatial import distance
from scipy.ndimage import gaussian_filter

# _________________________________________________________________________________________________
def main(args):
  sim_model = args["model"]
  input_dir = args["idir"]
  output_dir = args["odir"]
  mod = 1.0

  if sim_model == "rd":
    vm_path = "{}/rd_250um/sim/vm.igb".format(input_dir)
    #msh_path = 
    mod = 0.25
  elif sim_model == "pie":
    vm_path = "{}/pie_1000um/vm_mv.igb".format(input_dir)
    #msh_path = 
  else:
    print("Invalid simulation model type specified!")
    return

  if not os.path.exists(vm_path):
    raise RuntimeError("Vm file not found: {}".format(vm_path))

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
  t_range = np.arange(t_start, t_end, t_step, dtype=int)

  fig, ax = plt.subplots()
  trajectory = np.zeros((len(t_range), 2))
  vm_slices = np.zeros((len(t_range), Mpts, Mpts))
  contour_vm_list = []
  contour_gvm_list = []
  max_vmc = -1
  max_vmgc = -1

  for it, t in enumerate(t_range):
    vm_slice = vm_dat[t,:Npts].reshape((Mpts, Mpts))
    vm_slice = gaussian_filter(vm_slice, sigma=1)
    #t_slice = int((t-t_start)/t_step)
    vm_slices[it] = vm_slice
    
    measured_contours_vm = measure.find_contours(vm_slice, -75.0)
    contours_vm = [np.array(item) for item in measured_contours_vm if len(item) > min_len]
    contour_vm = np.vstack(contours_vm)

    if t > t_start:
      grad_vm = (vm_slices[it-1] - vm_slices[it])/dt
      measured_contours_gvm = measure.find_contours(grad_vm, 0.0)
      contours_gvm = [np.array(item) for item in measured_contours_gvm if len(item) > min_len]
      contour_gvm = np.vstack(contours_gvm)

      cdists = distance.cdist(contour_vm, contour_gvm).min(axis=1)
      sdists = np.linalg.norm(trajectory[it-1] - contour_vm, axis=1)

      if t > t_start + 3*t_step:
        tdists = 0.9*cdists+0.1*sdists
        #tdists = 0.8*cdists+0.2*sdists
        #tdists = 0.7*cdists+0.3*sdists
        #tdists = 0.6*cdists+0.4*sdists
      else:
        tdists = cdists

      idx = np.argmin(tdists)
      trajectory[it] = contour_vm[idx]
      contour_gvm_list.append(contours_gvm)
    else:
      contour_gvm_list.append(contours_vm)
      trajectory[it] = [0.0, 0.0]
    
    contour_vm_list.append(contours_vm)

    max_vmc = max(max_vmc, len(contour_vm_list[-1]))
    max_vmgc = max(max_vmgc, len(contour_gvm_list[-1]))

  im = ax.imshow(np.ones(vm_slices[0].shape)*(0.0), cmap='Greys')
  pvm = []
  pvmg = []

  for it in range(max_vmc):
    p1, = ax.plot(contour_vm_list[0][0][:,1], contour_vm_list[0][0][:,0], linewidth=2, color='blue')
    pvm.append(p1)

  p3, = ax.plot(trajectory[:,1], trajectory[:,0], linewidth=2, color='black')
  p4, = ax.plot([trajectory[0,1]], [trajectory[0,0]], marker='o', color='yellow', markeredgecolor='black') 

  ax.set_xticks([])
  ax.set_yticks([])
  plt.gca().invert_yaxis()

  # output trajectory as file
  np.save("{}/trajectory_data_{}.npy".format(input_dir, sim_model), trajectory*mod, allow_pickle=False)

  def update(frame):
    contours_vm = contour_vm_list[frame]

    for it in range(max_vmc):
      if it < len(contours_vm):
        pvm[it].set_xdata(contours_vm[it][:,1])
        pvm[it].set_ydata(contours_vm[it][:,0])
      else:
        pvm[it].set_xdata([])
        pvm[it].set_ydata([])

    p3.set_xdata(trajectory[1:frame,1])
    p3.set_ydata(trajectory[1:frame,0])

    p4.set_xdata([trajectory[frame,1]])
    p4.set_ydata([trajectory[frame,0]])

    tra_dir = "{}/trajectories_{}".format(output_dir, sim_model)

    if not os.path.exists(tra_dir):
      os.makedirs(tra_dir)

    plt.savefig("{}/traj_{}_{:04d}.png".format(tra_dir, sim_model, t_step*frame + 340), bbox_inches='tight', dpi=200)
    
    return (pvm, pvmg, p3, p4, im)

  num_frames = int((t_end-t_start)/t_step)
  ani = animation.FuncAnimation(fig=fig, func=update, frames=num_frames, interval=100, repeat=False)
  plt.show()

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Phase singularity tracking based on Vm.')
  parser.add_argument('--idir', help='simulation directory', required=True)
  parser.add_argument('--model', help='simulation model',    required=True)
  parser.add_argument('--odir', help='output directory',     required=True)
  args = vars(parser.parse_args())
  main(args)