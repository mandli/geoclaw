#!/usr/bin/env python

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

outdir = Path('_output')
d = np.loadtxt(outdir / 'fort.FG1.valuemax')
b = np.loadtxt(outdir / 'fort.FG1.aux1')
topo = b[:,4]  # for level 3
plt.figure(401)
plt.clf()
