#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import argparse
import subprocess
import numpy as np
import meshgen as mg

import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

np.set_printoptions(suppress=True)

# _________________________________________________________________________________________________
def main(args):
  setup_ID = "A4_scaling"
  runs = args["runs"]
  modes = args["mode"].split(",")
  setup_dir = "setups/{}".format(setup_ID)
  result_dir = args["outdir"]
  Nproc = args["np"]
  mesh_res_um = args["res"]

  depth_start = np.float32(args["dstart"])*1e-4
  depth_end   = np.float32(args["dend"])*1e-4
  depth_step  = np.float32(args["dstep"])*1e-4
  mesh_width_cm = np.arange(depth_start, depth_end + 1e-3, depth_step)

  init_time_data = np.zeros((len(mesh_width_cm), 8))
  solve_time_data = np.zeros((len(mesh_width_cm), 8))
  labels = ["RD", "PIE[-]", "PIE[A]", "PIE[C]", "PIE[AC]", "PIE[ACK]", "PIE[ACKD]"]
  pln_filepath = "{}/A4_scaling.plan.json".format(setup_dir)

  for it, width_cm in enumerate(mesh_width_cm):
    msh_dir = "{}/mesh_{:.3f}cm".format(setup_dir, width_cm)
    msh_filepath = "{}/mesh_{:.3f}cm".format(msh_dir, width_cm)
    pts_filepath = "{}.bpts".format(msh_filepath)

    # generate mesh if not present
    if (not os.path.isfile(pts_filepath)) or args["remesh"]:
      if not os.path.exists(msh_dir):
        os.makedirs(msh_dir)

      dims_cm = [8, 8, width_cm]
      msh = mg.mesher_gen(dims_cm, mesh_res_um, 0, msh_filepath)
      msh.add_block([0, 0, 0], dims_cm, 1, False)
      msh.generate()

    #pts_data = np.fromfile(pts_filepath, dtype=np.float32, offset=1024) # carp_bin
    Npts = np.fromfile(pts_filepath, dtype=np.uint32, count=1, sep=' ')[0]
    init_time_data[it,0] = Npts
    solve_time_data[it,0] = Npts

    # run CARP simulation
    if "RD" in modes:
      rd_outdir = "{}/RD/{:.3f}cm".format(result_dir, width_cm)
      rd_timefile = "{}/exec_times.txt".format(rd_outdir)

      if (not os.path.isfile(rd_timefile) or args["recompute"]):
        if not os.path.exists(rd_outdir):
          os.makedirs(rd_outdir)

        subprocess.run("pie-solver --pln2par --msh={} --pln={} --out={}".format(msh_filepath, pln_filepath, rd_outdir).split(' '))
        t1 = time.time()
        subprocess.run("mpirun -np {:d} openCARP +F {}/lat.par".format(Nproc, rd_outdir).split(' '))
        t2 = time.time()
        subprocess.run("mpirun -np {:d} openCARP +F {}/prp.par".format(Nproc, rd_outdir).split(' '))
        t3 = time.time()
        subprocess.run("mpirun -np {:d} openCARP +F {}/sim.par".format(Nproc, rd_outdir).split(' '))
        t4 = time.time()

        np.savetxt("{}/exec_times.txt".format(rd_outdir), np.array([t2 - t1, t3 - t2, t4 - t3, t4 - t1], dtype=np.float32))

      rd_data = np.loadtxt(rd_timefile)
      rd_init = np.sum(rd_data[:2])
      rd_solve = rd_data[3]
    else:
      rd_init = 0.0
      rd_solve = 0.0
    
    init_time_data[it,1] = rd_init
    solve_time_data[it,1] = rd_solve

    # run PIE simulation
    if "PIE" in modes:
      pie_it = 0
      pie_configs = {"a": "--no-cvr --no-apdr --no-curv --no-dif --dat=32",
                     "b": "--no-cvr --no-curv --no-dif --dat=32",
                     "c": "--no-apdr --no-curv --no-dif --dat=32",
                     "d": "--no-curv --no-dif --dat=32",
                     "e": "--no-dif --dat=32",
                     "f": "--dat=32"}
      
      for key, flags in pie_configs.items():
        pie_outdir = "{}/PIE/{}_{:.3f}cm".format(result_dir, key, width_cm)
        pie_timefile = "{}/exec_times.txt".format(pie_outdir)

        if (not os.path.isfile(pie_timefile)) or args["recompute"]:
          if not os.path.exists(pie_outdir):
            os.makedirs(pie_outdir)

          pie_cmd = "pie-solver --msh={} --pln={} --out={} --np={:d} {}".format(msh_filepath, pln_filepath, pie_outdir, Nproc, flags)
          subprocess.run(pie_cmd.split(' '))

        pie_init, pie_solve = np.loadtxt(pie_timefile)
        init_time_data[it,2+pie_it] = pie_init
        solve_time_data[it,2+pie_it] = pie_solve
        pie_it += 1

  # plot
  if args["plot"] == True:
    outfile_id = "{}/{}".format(result_dir, setup_ID)
    np.savetxt("{}_init_times.txt".format(outfile_id), init_time_data)
    np.savetxt("{}_solve_times.txt".format(outfile_id), solve_time_data)

    fig = plt.figure(figsize=(10, 5))
    ax = []
    rows, cols = (1, 1)
    grd = gs.GridSpec(rows, cols)

    #plt.rc('text', usetex=True)
    #plt.rc('font', family='Times New Roman', size=18)

    xr = [0, 3.2e6]
    yr = [2e-1, 3e4]

    titles = ["Simulation Time"]
    ylabels = ["Time in seconds"]

    for row in range(0, rows):
      for col in range(0, cols):
        it = row*cols + col
        ax.append(fig.add_subplot(grd[row, col]))
        ax[it].set(xlim=xr, ylim=yr)
        ax[it].set_title(titles[it])
        ax[it].set(xlabel='Number of Nodes', ylabel=ylabels[it])
        ax[it].set_yscale('log')
        ax[it].grid(which='minor')
        ax[it].grid(which='major')
    
    #for it in range(solve_time_data.shape[1]-1):
    #  if init_time_data[:,1+it].sum() > 1e-3:
    #    ax[0].plot(init_time_data[:,0], init_time_data[:,1+it], label=labels[it])

    for it in range(solve_time_data.shape[1]-1):
      if solve_time_data[:,1+it].sum() > 1e-3:
        ax[0].plot(solve_time_data[:,0], solve_time_data[:,1+it], label=labels[it])

    plt.tight_layout()
    plt.savefig("{}_times.pdf".format(outfile_id), dpi=200, bbox_inches='tight')
    plt.show()

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Performance scaling comparison benchmark between RD and PIE solvers.')
  parser.add_argument('--mode',      help='select method(s) to evaluate (RD, PIE) (comma separated).', type=str, default="RD,PIE")
  parser.add_argument('--res',       help='specifies generated mesh resolution in um.',  type=int, default=250)
  parser.add_argument('--dstart',    help='start mesh depth in um.',                     type=int, default=250)
  parser.add_argument('--dend',      help='end mesh depth in um.',                       type=int, default=2500)
  parser.add_argument('--dstep',     help='step mesh depth in um.',                      type=int, default=250)
  parser.add_argument('--np',        help='number of threads to use for the benchmark.', type=int, default=32)
  parser.add_argument('--remesh',    help='flag to force mesh generation.',              action='store_true', default=False)
  parser.add_argument('--recompute', help='flag to force recomputation of all results.', action='store_true', default=False)
  parser.add_argument('--plot',      help='flag to plot the results.',                   action='store_true', default=False)
  parser.add_argument('--runs',      help='flag conduct runtime analysis.',              type=int, default=1)
  parser.add_argument('--outdir',    help='specifies output directory.',                 type=str, default="results")
  args = vars(parser.parse_args())
  main(args)
