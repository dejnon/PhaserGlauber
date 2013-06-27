#!/usr/bin/python
# Script for running automated execution of ragne of values
# C++ Program that is being used has the following call syntax:
# programname [maxt] [l] [w0] [cmean] [csigma] [cmodename]
# [maxt]      - Maximal time-step threshold (int)
# [l]         - Latice size (int)
# [w0]        - Ordering parameter (float), [0,1]
# [cmean]     - Mean value of c-parameter, (flt), [0,1]
# [csigma]    - C-parameter's standard deviation (or c=[cmean-csigma; cmean+sigma]), (flt), [0,1]
# [cmodename] - C-mode: 0-well / 1-gaussian / 2-uniform / 3-triangle, (int), {0,1,2,3}
# [verbose]   - Display detailed system progression, (int), {0,1}

import os
import argparse
import numpy as np
import subprocess
import multiprocessing

W0_RANGE = np.linspace(0.51, 1.0, 5)
W0_RANGE = [W0_RANGE[0]]
CM_RANGE = np.linspace(0.0, 1.0, 30)
CS_RANGE = np.linspace(0.0, 1.0, 30)
# MODE = [0,1,2,3] # - triangle

maxt = 10000000
l = 100
verbose = 0
averages = 5
cmodename = 3

output_path = "./_w0-0.51_cmean-all_csigma-all_maxt-10-7_l-100_w0grain-1_cgrain-30/"
current_path = os.getcwd()

if not os.path.exists(output_path):
    os.makedirs(output_path)
os.chdir(output_path)

def run_simulation(param):
  os.system(param)

list_of_simulations = []
if os.system("gcc %s/main.cpp -o ../main -lstdc++ -lgsl -lgslcblas -lm" % current_path) == 0:
  for w0 in W0_RANGE:
    for cmean in CM_RANGE:
      for csigma in CS_RANGE:
          filename = "w0-%.4f_cmean-%.4f_csigma-%.4f_maxt-%d_l-%d.csv" % (w0, cmean, csigma, maxt, l)
          simul = "%s/main %d %d %f %f %f %d %d %d %s" % (current_path, maxt, l, w0, cmean, csigma, cmodename, verbose, averages, filename)
          list_of_simulations.append(simul)
  # Run 12 instances at the time
  pool = multiprocessing.Pool(16)
  pool.map(run_simulation, list_of_simulations)



