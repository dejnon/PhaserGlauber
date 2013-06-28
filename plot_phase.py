
# main 100000000 100 0.510000 0.2 0.2 3 1 3 przejscie.csv > przejscie.txt
# 2D array of 0 and 1 obtained from the abote command
img = [[0,0,1,0],[0,0,0,0]]

import matplotlib.pyplot as plt
import numpy as np


# im = plt.matshow(img)
# # plt.colorbar(im, orientation='horizontal')
# plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(img, interpolation='nearest')

# plt.show()
fig.tight_layout()
fig.set_size_inches(10,100)
fig.savefig("micro.png", dpi=100)

