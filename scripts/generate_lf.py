#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import argparse
import subprocess
import numpy as np

# Example: python3 ./scripts/generate_sigma.py --mshdir=setups/A3_wholeheart/mesh_3000um/ --pln=setups/A3_wholeheart/A3_wholeheart-3k.plan.json --gmodel=Clerc76 --np=16

# _________________________________________________________________________________________________
class g_set:
  def __init__(self, g_il, g_it, g_el, g_et):#name, 
    #self.name = name
    self.g_i = np.array([g_il, g_it, g_it])
    self.g_e = np.array([g_el, g_et, g_et])

# _________________________________________________________________________________________________
def read_vtx_file(vtx_filepath):
  return np.loadtxt(vtx_filepath, skiprows=2, dtype=int)

# _________________________________________________________________________________________________
def load_plan(pln_filepath):
  plan_dict = []

  with open(pln_filepath, 'r') as plan_file:
    plan_dict = json.load(plan_file)

  return plan_dict

def save_plan(pln_filepath, plan_dict):
  with open(pln_filepath, 'w') as plan_file:
    json.dump(plan_dict, plan_file, indent=2)

# _________________________________________________________________________________________________
def extract_g_m(fnc, g_model):
  g_m = []

  if fnc["EP"] != None:
    #fnc_g = fnc["CV"]["G"]
    #g_i = np.array([fnc_g["gel"], fnc_g["gen"], fnc_g["get"]])
    #g_e = np.array([fnc_g["gil"], fnc_g["gin"], fnc_g["git"]])
    g_m = (g_model.g_e*g_model.g_i)/(g_model.g_e+g_model.g_i)
  else:
    g_m = np.ones(3)*fnc["CV"]["G"]["gbath"]

  return g_m

# _________________________________________________________________________________________________
def write_sigma_file(pln_dict, g_model, sigma_filepath):
  fnc_g_m_map = {}

  for label, fnc in pln_dict["functions"].items():
    fnc_g_m_map[label] = extract_g_m(fnc, g_model)

  content = ""

  for label, cfg in pln_dict["config"].items():
    g_m = fnc_g_m_map[cfg["func"]]

    for tag in cfg["tags"]:
      content += "{} {} {} {}\n".format(tag, g_m[0], g_m[1], g_m[2])

  with open(sigma_filepath, "w") as sigma_file:
    sigma_file.write(content)

# _________________________________________________________________________________________________
def generate_lfdata(msh_path, pln_dict, g_model, g_scale=1e-3, num_proc=8):
  msh_dir = os.path.dirname(msh_path)
  lfd_dir = "{}/leadfields".format(msh_dir) #lfdata
  g_m_file = "{}/sigma-file.greg".format(lfd_dir)
  ref_vtx = read_vtx_file("{}/lf_ref.vtx".format(msh_dir))
  src_vtx_file = "{}/lf_src.vtx".format(msh_dir)

  if not os.path.exists(lfd_dir):
    os.makedirs(lfd_dir)

  if not os.path.exists(g_m_file):
    write_sigma_file(pln_dict, g_model, g_m_file)

  gfsolve_cmd = "gizmo_gfsolve --model-name {} --model-format carp_bin --sigma-file {} --sigma-scale {:f} --ref-vtx {:d} --src-vtx {} --out-format {}/grfnc_ref_%R_src_%S --binary-gf --np {:d} --lin-solver amgcl".format(msh_path, g_m_file, g_scale, ref_vtx, src_vtx_file, lfd_dir, num_proc)
  print(gfsolve_cmd)
  subprocess.run(gfsolve_cmd.split(' '))

# _________________________________________________________________________________________________
def main(args):
  msh_filepath = args["msh"]
  pln_filepath = args["pln"]
  g_model = args["gmodel"]
  g_scale = args["gscale"]
  num_proc = args["np"]

  # check if gizmo is found

  #if not os.path.exists():
  #  print("Error: mesh dir not found: {}".format(msh_filepath)) # also check for vtx files!
  #  return
  
  if not os.path.exists(pln_filepath):
    print("Error: plan file not found: {}".format(pln_filepath))
    return

  pln_dict = load_plan(pln_filepath)
  g_settings = {"Default":      g_set(0.174,  0.019,  0.625, 0.236 ),
                "Clerc76":      g_set(0.17,   0.019,  0.62,  0.24  ),
                "Roberts79":    g_set(0.28,   0.026,  0.22,  0.13  ),
                "Roberts82":    g_set(0.34,   0.060,  0.12,  0.08  ),
                "Courtemanche": g_set(1.0189, 0.2197, 3.716, 2.7749)}
  
  g_vals = g_settings["Default"]

  if g_model in g_settings:
    g_vals = g_settings[g_model]
    generate_lfdata(msh_filepath, pln_dict, g_vals, g_scale, num_proc)
  else:
    print("Error: invalid g_model specified: {}, must be either Clerc76, Roberts79, Roberts82 or Courtemanche.".format(g_model))
    return

# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Generates sigma file for LF.')
  parser.add_argument('--msh',    help='Mesh file.',            required=True)
  parser.add_argument('--pln',    help='Simulation plan.',      required=True)
  parser.add_argument('--gmodel', help='Conductivity model.',   required=True)
  parser.add_argument('--gscale', help='Conductivity scaling.', default=1e-3, type=float)
  parser.add_argument('--np',     help='Number of processors.', default=8, type=int)
  args = vars(parser.parse_args())
  main(args)