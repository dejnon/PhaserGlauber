import os
import argparse
import numpy as np
import subprocess

W0_RANGE = np.linspace(0.51, 1.0, 5)
CM_RANGE = np.linspace(0.0, 1.0, 10)
CS_RANGE = np.linspace(0.0, 1.0, 10)
# MODE = [0,1,2,3] # - triangle


maxt = 1000000
l = 100
verbose = 0
averages = 100

if os.system("gcc main.cpp -o main -lstdc++ -lgsl -lgslcblas") == 0:
  cmodename = 3
  for w0 in W0_RANGE:
    for cmean in CM_RANGE:
      for csigma in CS_RANGE:
          filename = "w0-%.4f_cmean-%.4f_csigma-%.4f_maxt-%d_l-%d.csv" % (w0, cmean, csigma, maxt, l)
          simul = "./main %d %d %f %f %f %d %d %d %s" % (maxt, l, w0, cmean, csigma, cmodename, verbose, averages, filename)
          p = os.system(simul)

