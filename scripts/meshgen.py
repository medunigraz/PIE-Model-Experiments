#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess

import numpy as np


# _________________________________________________________________________________________________
class mesher_gen:
  def __init__(self, dimensions_cm, resolution_um, fibdir_deg, msh_filepath):
    dimensions_cm = np.array(dimensions_cm, np.float32)
    size_cm = dimensions_cm.copy()
    size_cm[size_cm < 0.0] = resolution_um/1e4
    center_cm = (size_cm/2.0)
    
    self.region_id = 0
    self.command = ["mesher"]
    self.command += ["-resolution", "{{{} {} {}}}".format(resolution_um, resolution_um, resolution_um)]
    self.command += ["-size", "{{{:f} {:f} {:f}}}".format(size_cm[0], size_cm[1], size_cm[2])]
    self.command += ["-center", "{{{:f} {:f} {:f}}}".format(center_cm[0], center_cm[1], center_cm[2])]
    self.command += ["-fibers.rotEpi", "{}".format(fibdir_deg)]
    self.command += ["-fibers.rotEndo", "{}".format(fibdir_deg)]
    self.command += ["-mesh", "{}".format(msh_filepath)]
    
    self.center_cm = center_cm
    self.regions = []
    self.rids = ""
    self.tags = ""
    self.msh_filepath = msh_filepath
    
  def add_block(self, p0_cm, p1_cm, tag, is_bath):
    rid = self.region_id
    p0_cm = np.array(p0_cm, dtype=np.float32)
    p1_cm = np.array(p1_cm, dtype=np.float32) - self.center_cm

    if not is_bath:
      self.rids += "{},".format(rid+1)
      self.tags += "{},".format(tag)
    
    self.regions += ["-regdef[{}].type".format(rid), "{}".format(0)]
    self.regions += ["-regdef[{}].tag".format(rid),  "{}".format(tag)]
    self.regions += ["-regdef[{}].bath".format(rid), "{:g}".format(is_bath)]
    self.regions += ["-regdef[{}].p0".format(rid),   "{{{:f} {:f} {:f}}}".format(p0_cm[0], p0_cm[1], p0_cm[2])]
    self.regions += ["-regdef[{}].p1".format(rid),   "{{{:f} {:f} {:f}}}".format(p1_cm[0], p1_cm[1], p1_cm[2])]
    self.region_id += 1
  
  def add_sphere(self, p0_cm, rad_cm, tag, is_bath):
    rid = self.region_id
    p0_cm = np.array(p0_cm, dtype=np.float32)

    if not is_bath:
      self.rids += "{},".format(rid+1)
      self.tags += "{},".format(tag)
    
    self.regions += ["-regdef[{}].type".format(rid), "{}".format(1)]
    self.regions += ["-regdef[{}].tag".format(rid),  "{}".format(tag)]
    self.regions += ["-regdef[{}].bath".format(rid), "{:g}".format(is_bath)]
    self.regions += ["-regdef[{}].p0".format(rid),   "{{{:f} {:f} {:f}}}".format(p0_cm[0], p0_cm[1], p0_cm[2])]
    self.regions += ["-regdef[{}].rad".format(rid),  "{:f}".format(rad_cm)]
    self.region_id += 1
  
  def add_cylinder(self, p0_cm, p1_cm, rad_cm, clen, tag, is_bath):
    rid = self.region_id
    p0_cm = np.array(p0_cm, dtype=np.float32)

    if not is_bath:
      self.rids += "{},".format(rid+1)
      self.tags += "{},".format(tag)
    
    self.regions += ["-regdef[{}].type".format(rid), "{}".format(2)]
    self.regions += ["-regdef[{}].tag".format(rid),  "{}".format(tag)]
    self.regions += ["-regdef[{}].bath".format(rid), "{:g}".format(is_bath)]
    self.regions += ["-regdef[{}].p0".format(rid),   "{{{:f} {:f} {:f}}}".format(p0_cm[0], p0_cm[1], p0_cm[2])]
    self.regions += ["-regdef[{}].p1".format(rid),   "{{{:f} {:f} {:f}}}".format(p1_cm[0], p1_cm[1], p1_cm[2])]
    self.regions += ["-regdef[{}].rad".format(rid),  "{:f}".format(rad_cm)]
    self.regions += ["-regdef[{}].cyllen".format(rid),  "{:f}".format(clen)]
    self.region_id += 1
  
  def generate(self):
    if self.region_id > 0:
      self.command += ["-numRegions", "{}".format(self.region_id)]
      self.command += ["-first_reg", "{}".format(0)]
      self.command += self.regions
      
    # generate mesh in carp_txt format
    dry_cmd = " ".join(self.command)
    print(dry_cmd)
    subprocess.run(self.command)

    # convert to binary
    mt_cmd = "meshtool extract mesh -msh={} -tags={} -submsh={} -ifmt=carp_txt -ofmt=carp_bin".format(self.msh_filepath, self.tags[:-1], self.msh_filepath)
    subprocess.run(mt_cmd.split(' '))

    # remove carp_txt files
    rm_cmd = "rm {} {} {} {} {} {} {}".format(self.msh_filepath+".eidx", 
                                              self.msh_filepath+".elem", 
                                              self.msh_filepath+".lon", 
                                              self.msh_filepath+".nod", 
                                              self.msh_filepath+".pts", 
                                              self.msh_filepath+".vec", 
                                              self.msh_filepath+".vpts")
    subprocess.run(rm_cmd.split(' '))


# _________________________________________________________________________________________________
def main(args):
  mesh_id = args["msh"]
  base_dir = os.path.dirname(os.path.abspath(sys.argv[0]))+"/../"
  setup_dir = base_dir + "setups/{}/".format(mesh_id)
  resolutions_um = [int(res_um) for res_um in args["res"].split(',')]
  
  for resolution_um in resolutions_um:
    msh_dir = "{}mesh_{}um/".format(setup_dir, resolution_um)
    msh_filepath = "{}mesh_{}um".format(msh_dir, resolution_um)
    res_cm = resolution_um*1e-4
    
    if os.path.isfile(msh_filepath+".pts") or os.path.isfile(msh_filepath+".bpts") or os.path.isfile(msh_filepath+".vtk"):
      continue
    
    if not os.path.exists(msh_dir):
      os.makedirs(msh_dir)

    if mesh_id == "A1_anatomical":
      if resolution_um == 1000:
        print("Provided by repository.")
      else:
        r_scar = 3.75
        r_bz = 4.0
        w_bz = 0.5
        dims_cm = np.array([10, 10,  res_cm])
        cent_cm = dims_cm/2.0
        p0 = [cent_cm[0] - w_bz/2.0, cent_cm[1] - r_scar, 0]
        p1 = [cent_cm[0] + w_bz/2.0, dims_cm[1] - p0[1], dims_cm[2]]
        msh = mesher_gen(dims_cm, resolution_um, 90, msh_filepath)
        msh.add_block([0, 0, 0], dims_cm, 1, False)
        msh.add_sphere(cent_cm, r_bz, 2, False)
        msh.add_sphere(cent_cm, r_scar, 3, False)
        msh.add_block(p0, p1, 4, False)
        msh.generate()
        mt_cmd = "meshtool extract mesh -msh={} -tags=1,2,4 -submsh={} -ifmt=carp_bin -ofmt=carp_bin".format(msh_filepath, msh_filepath)
        subprocess.run(mt_cmd.split(' '))
    elif mesh_id == "A2_functional":
      dims_cm = [8, 8, res_cm]
      msh = mesher_gen(dims_cm, resolution_um, 0, msh_filepath)
      msh.add_block([0, 0, 0], dims_cm, 1, False)
      msh.generate()
    elif mesh_id == "A3_wholeheart":
      print("Please download A3_wholeheart from Zenodo.")
    elif mesh_id == "B1_restitution":
      dims_cm = [25, 1, res_cm]
      msh = mesher_gen(dims_cm, resolution_um, 90, msh_filepath)
      msh.add_block([0, 0, 0], dims_cm, 1, False)
      msh.generate()
    elif mesh_id == "B2_curvature":
      dims_cm = [4, 2, res_cm]
      msh = mesher_gen(dims_cm, resolution_um, 0, msh_filepath)
      msh.add_block([0,   0, 0], dims_cm, 1, False)
      msh.add_block([1,   0, 0], [3, 0.9, res_cm], 2, True)
      msh.add_block([1, 1.1, 0], [3,   2, res_cm], 2, True)
      msh.generate()
    elif mesh_id == "B3_diffusion":
      dims_cm = [2, 6, res_cm]
      msh = mesher_gen(dims_cm, resolution_um, 90, msh_filepath)
      msh.add_block([0, 0, 0], [1, 6, res_cm], 1, False)
      msh.add_block([1, 0, 0], [2, 6, res_cm], 2, False)
      msh.generate()
    else:
      print("meshgen.py: invalid mesh ID! aborting ...")
      raise


# _________________________________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Mesh generator for EP and EK simulations.')
  parser.add_argument('--msh', help='Mesh ID to generate',                         required=True)
  parser.add_argument('--res', help='Mesh resolution(s), comma separated values.', required=True)
  parser.add_argument('--out', help='Output directory',                            required=True)
  args = vars(parser.parse_args())
  main(args)
