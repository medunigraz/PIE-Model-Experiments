#!/usr/bin/env python3

import os
import argparse
import numpy as np

# example: 
# >> python3 scripts/dat2bdat.py --idir=./setup/lead_field --odir=./setup/lead_field_bin

def main(args):
  idir = args["idir"]
  odir = args["odir"]

  if not os.path.exists(idir):
    print("Error: input directory does not exist: \n{}".format(idir))
    return
  
  if not os.path.exists(odir):
    os.makedirs(odir)
  
  for file in os.listdir(idir):
    ext = file[-4:]
    
    if ext == ".dat":
      data = np.loadtxt("{}/{}".format(idir, file))
      data.astype('float64').tofile("{}/{}.bdat".format(odir, file[:-4]))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Evaluation of the Pharma Study.')
  parser.add_argument('--idir', help='Input LF directory.',  type=str)
  parser.add_argument('--odir', help='Output LF directory.', type=str)
  args = vars(parser.parse_args())
  main(args)