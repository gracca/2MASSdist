#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# plot2MASSdist.py
#
# Copyright (C) 2014 Germán A. Racca
# E-Mail: <gracca[AT]gmail[DOT]com>
# License: GPLv3+

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import ScalarFormatter, MultipleLocator
from astropy import coordinates as coord

# True = fixed bins / False = variable bins
whatbin = False

# object
obj    = 'BHR001'
fname  = 'AvD-' + obj
if whatbin:
    figpng = fname + '-fixed.png'
    dist = 150
    derr = 50
    ymax = 2.0  # max Av in bottom plot
else:
    figpng = fname + '-variable.png'
    dist = 150
    derr = 50
    ymax = 2.0 # max Av in bottom plot

# galactic latitude
coo = coord.Galactic.from_name(obj)
#coo = coord.Galactic.from_name('BHR062')
lat = coo.b.degree

# function for fixed bins
def fix_bin(start, w, end):
    """Calculate centers and intervals for fixed bins"""
    s = -w / 2
    c = np.arange(start, end, s)
    b = [(i-s, i+s) for i in c]
    return c, b

# function for variable bins
def var_bins(n, start):
    """Calculate centers and intervals for variable bins"""
    c = [start]
    for i in range(0, n):
        c.append(c[i] - c[i]*0.09)
    b = [(c[i] + c[i]*0.09, c[i+1]) for i in range(0, n)]
    return c, b

# function for the obscuration
def extin(b, r):
    """Calculate the extinction for latitude b and distance r
       following the obscuration model of Bahcall & Soneira (1980)
    """
    H = 100.0
    b = b * np.pi / 180.0
    Ainf = 0.165 * (1.192 - np.abs(np.tan(b))) / np.abs(np.sin(b))
    Avis = Ainf * (1.0 - np.exp(-np.abs(np.sin(b)) * r / H))
    return Avis

# main program
if whatbin:
    print("Fixed bins: %s \n" % obj)
    c, b = fix_bin(1000, 30, 90)
else:
    print("Variable bins: %s \n" % obj)
    c, b = var_bins(27, 1000)

for i, j in zip(c, b):
   print("%.1f --> (%.1f, %.1f) --> %.1f" % (i, j[0], j[1], j[0] - j[1]))
print("")

# calculations
Av, sigAv, D, sigD = np.loadtxt(fname, usecols=(2, 3, 4, 5), unpack=True)

meanAvD = np.empty((0, 6), dtype=int)
for ci, bi in zip(c, b):
    res = np.empty((0, 4), dtype=int)
    for av, s_av, d, s_d in zip(Av, sigAv, D, sigD):
        if d < bi[0] and d >= bi[1]:
            res = np.append(res, [[av, s_av, d, s_d]], axis=0)
    n = res.shape[0]
    m_av   = res[:,0].mean()
    m_s_av = res[:,0].std() / np.sqrt(n)
    m_d    = res[:,2].mean()
    m_s_d  = res[:,3].mean()
    meanAvD = np.append(meanAvD, [[n, ci, m_av, m_s_av, m_d, m_s_d]], axis=0)
    print("%3d %.2f %.4f %.4f %.4f %.4f" % (n, ci, m_av, m_s_av, m_d, m_s_d))

# plot Av vs. D
fig = plt.figure(figsize=(9, 8))
#fig.suptitle(obj, fontsize=16)
gs = gridspec.GridSpec(2, 1, height_ratios=[1.5, 1])

ax1 = plt.subplot(gs[0])
ax1.set_xlabel("Distance [pc]", fontsize=14)
ax1.set_ylabel(r"A$\mathrm{_V}$ [mag]", fontsize=14)
ax1.set_xscale('log')
ax1.set_xlim(10, 2000)
ax1.set_ylim(0, 4)
ax1.xaxis.set_major_formatter(ScalarFormatter())  # no scientific notation
ax1.yaxis.set_major_locator(MultipleLocator(0.5))
ax1.yaxis.set_minor_locator(MultipleLocator(0.25))
ax1.errorbar(D, Av, xerr=sigD, yerr=sigAv, fmt='o')
ax1.text(0.1, 0.9, obj, fontsize=16, ha='center', va='center',
        transform=ax1.transAxes)  # label
ax1.plot([dist, dist], [0, 5], 'k--', lw=2)  # mark distance

# plot obscuration
x = np.arange(0, 2010, 10)
y = extin(lat, x)
ax1.plot(x, y, 'r-', lw=2)

# plot mean Av vs. mean D
ax2 = plt.subplot(gs[1])
ax2.set_xlabel("Distance [pc]", fontsize=14)
ax2.set_ylabel(r"mean A$\mathrm{_V}$ [mag]", fontsize=14)
ax2.set_xlim(0, 600)
ax2.set_ylim(0, ymax)
ax2.xaxis.set_minor_locator(MultipleLocator(20))
ax2.yaxis.set_major_locator(MultipleLocator(0.2))
ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
ax2.errorbar(meanAvD[:, 4], meanAvD[:, 2], yerr=meanAvD[:, 3], fmt='s-')
ax2.plot([dist, dist], [0, 3], 'k--', lw=2)  # mark distance
ax2.text(0.825, 0.1, 'D = ' + str(dist) + u" \u00B1 " + str(derr) + ' pc',
         fontsize=16, ha='center', va='center',
         transform=ax2.transAxes)  # label

print("\nPlot for " + obj + " is ready!\n")

# save and show figure
fig.savefig(figpng, bbox_inches='tight')
plt.show()

