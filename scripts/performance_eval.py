#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import argparse
import subprocess
import numpy as np

# _________________________________________________________________________________________________
def read_pie_runtime(pie_dir):
  return np.loadtxt("{}/exec_times.txt".format(pie_dir))[1]

# _________________________________________________________________________________________________
def read_gizmo_runtime(gizmo_dir):
  with open("{}/stdout.txt".format(gizmo_dir), "r") as text_file:
    for line in text_file:
      if "# run simulation" in line:
        return (' '.join(line.rstrip().split())).split(' ')[4]

  print("Error: GIZMO output does not contain simulation time!")

# _________________________________________________________________________________________________
def main(args):
  msh_filepath = args["msh"]
  pln_filepath = args["pln"]
  out_dir      = args["out"]
  num_proc = int(args["np"])
  num_runs = int(args["runs"])

  # PIE simulations
  pie_cmd = "./bin/pie-solver --msh={} --pln={} --out={} --np={:d} --dat=32".format(msh_filepath, pln_filepath, out_dir, num_proc)
  pie_runtimes = np.zeros(num_runs)

  for run in range(num_runs):
    subprocess.call(pie_cmd.split(' ') ,stdout=subprocess.DEVNULL)
    pie_runtimes[run] = read_pie_runtime(out_dir)

  # LF reconstruction using GIZMO
  gizmo_cmd = "./bin/gizmo --model-name {} --model-format carp_bin --sim-plan {} --extern-vm {}/vm_mv.igb --trace-output json --output-dir {} --np {:d} --lin-solver amgcl".format(msh_filepath, pln_filepath, out_dir, out_dir, num_proc)
  gizmo_runtimes = np.zeros(num_runs)
  #print(gizmo_cmd)

  for run in range(num_runs):
    with open("{}/stdout.txt".format(out_dir), "w") as stdout_file:
      p1 = subprocess.Popen(gizmo_cmd.split(' '), stdin=subprocess.PIPE, stdout=stdout_file)
      p1.communicate(input="o".encode())[0]
    gizmo_runtimes[run] = read_gizmo_runtime(out_dir)

  if num_runs > 1:
    print(pie_runtimes, np.mean(pie_runtimes), np.std(pie_runtimes))
    print(gizmo_runtimes, np.mean(gizmo_runtimes), np.std(gizmo_runtimes))
    print("\\SI{{{:d}}}{{\\second}} + \\SI{{{:d}}}{{\\second}}".format(int(np.mean(pie_runtimes)), int(np.mean(gizmo_runtimes))))
  elif num_runs == 1:
    print(pie_runtimes[0])
    print(gizmo_runtimes[0])
  else:
    print("Error: invalid number of runs specified!")

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Performance evaluation of PIE simulation and lead-field reconstruction.')
  parser.add_argument('--msh',  help='Input mesh filename.',    required=True)
  parser.add_argument('--pln',  help='Input simplan filename.', required=True)
  parser.add_argument('--out',  help='Output directory.',       required=True)
  parser.add_argument('--np',   help='Number of processors.',   default=32)
  parser.add_argument('--runs', help='Number of runs.',         default=5)
  args = vars(parser.parse_args())
  main(args)
