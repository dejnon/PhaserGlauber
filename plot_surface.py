import os
import argparse
import numpy as np
import subprocess
import csv
import itertools
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import glob


catalogue = "_w0-0.1.51_cmean-all_csigma-all_maxt-10-7_l-100_w0grain-1_cgrain-30"
x_name = "cmean"
y_name = "csigma"
z_name = "BondDensity"
col_no = 2 # t-2 mcs-3 magnetization-4 bonddens-5
W0_RANGE = np.linspace(0.51, 1.0, 5)
w0 = W0_RANGE[0]
title = "Triangular "+z_name+" W0:%.2f Averages:5 Grain:30x30 MaxT:10^7" % w0

x = CM_RANGE = np.linspace(0.0, 1.0, 30)
y = CS_RANGE = np.linspace(0.0, 1.0, 30)

def fun(x,y):
  filename = [
    glob.glob(catalogue+'/w0-%.4f*%s-%.4f*%s-%.4f*' % (w0, x_name, x, y_name, y)),
    glob.glob(catalogue+'/w0-%.4f*%s-%.4f*%s-%.4f*' % (w0, y_name, y, x_name, x))
  ]
  if filename == [[],[]]:
    print str(x)+" "+str(y)
    raise Exception('No file foud matching criteria')
  filename = filename[0] if filename[0] else filename[1]

  with open(filename[0]) as csv_file:
    reader = csv.reader(csv_file, delimiter=' ', quotechar='"')
    column = []
    for row in reader:
      row = list(filter(("").__ne__, row)) # dubble spaces... ehhh...
      column.append(float(row[col_no]))
    average = sum(column) / len(column)
  return average

X, Y = np.meshgrid(x, y)
zs = np.array([fun(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
Z = zs.reshape(X.shape)


fig = plt.figure(figsize=(8, 6), dpi=100)
fig.text(.4, .95, title)
# bigfigure (comment any other)
ax = fig.add_subplot(111, projection='3d')
ax.view_init(20, -130)
p = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
cb = fig.colorbar(p, shrink=0.5)
ax.set_xlabel(x_name)
ax.set_ylabel(y_name)
ax.set_zlabel(z_name)


# # original
# ax = fig.add_subplot(221, projection='3d')
# ax.view_init(30, 180)
# ax.plot_surface(X, Y, Z, rstride=4, cstride=4, linewidth=0)
# ax.set_xlabel(x_name)
# ax.set_ylabel(y_name)
# ax.set_zlabel(z_name)

# # tails colored
# ax = fig.add_subplot(222, projection='3d')
# ax.view_init(20, -130)
# p = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# cb = fig.colorbar(p, shrink=0.5)
# ax.set_xlabel(x_name)
# ax.set_ylabel(y_name)
# ax.set_zlabel(z_name)

# # transparent
# ax = fig.add_subplot(223, projection='3d')
# ax.plot_surface(X, Y, Z, rstride=4, cstride=4, alpha=0.25)
# ax.view_init(50, -120)
# cset = ax.contour(X, Y, Z, zdir='z', offset=1, cmap=cm.coolwarm)
# cset = ax.contour(X, Y, Z, zdir='x', offset=1, cmap=cm.coolwarm)
# cset = ax.contour(X, Y, Z, zdir='y', offset=1, cmap=cm.coolwarm)
# ax.set_xlabel(x_name)
# ax.set_ylabel(y_name)
# ax.set_zlabel(z_name)



# ax = fig.add_subplot(224, projection='3d')
# ax.plot_surface(X, Y, Z, rstride=4, cstride=4, alpha=0.25)
# cset = ax.contour(X, Y, Z, zdir='z', offset=1, cmap=cm.coolwarm)
# cset = ax.contour(X, Y, Z, zdir='x', offset=1, cmap=cm.coolwarm)
# cset = ax.contour(X, Y, Z, zdir='y', offset=1, cmap=cm.coolwarm)
# ax.view_init(15, 100)
# ax.set_xlabel(x_name)
# ax.set_ylabel(y_name)
# ax.set_zlabel(z_name)

# plt.show()
fig.tight_layout()
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5,10.5)
fig.savefig(title.replace(" ","_").replace(":", "-").replace("^","-")+".png", dpi=100)
# fig.savefig(title.replace(" ","_").replace(":", "-").replace("^","-")+".eps", dpi=100)

