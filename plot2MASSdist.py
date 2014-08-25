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

obj    = 'NGC7023'
fname  = 'AvD-' + obj
figpng = fname + '.png'

# latitud galáctica
#coo = coord.Galactic.from_name('bhr ' + obj)
coo = coord.Galactic.from_name(obj)
lat = coo.b.degree

# número de bins
nbins = 22

# centros
c = [1000]
for i in range(0, nbins):
    c.append(c[i] - c[i]*0.09)

# bins
b = [(c[i] + c[i]*0.09, c[i+1]) for i in range(0, nbins)]

for i, j in zip(c, b):
    print("%.1f --> (%.1f, %.1f) --> %.1f" % (i, j[0], j[1], j[0] - j[1]))
print("")

# cálculos
Av, sigAv, D, sigD = np.loadtxt(fname, unpack=True)

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

def extin(b, r):
    """Calcula la extinción para una determinada dirección b y distancia r
       según el modelo 'obscuration' de Bahcall & Soneira (1980)
    """
    H = 100.0
    b = b * np.pi / 180.0
    Ainf = 0.165 * (1.192 - np.abs(np.tan(b))) / np.abs(np.sin(b))
    Avis = Ainf * (1.0 - np.exp(-np.abs(np.sin(b)) * r / H))
    return Avis

# gráfico Av vs. D
fig = plt.figure(figsize=(9, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[1.5, 1])

ax1 = plt.subplot(gs[0])
ax1.set_xlabel("Distance [pc]", fontsize=14)
ax1.set_ylabel(r"A$\mathrm{_V}$ [mag]", fontsize=14)
ax1.set_xscale('log')
ax1.set_xlim(10, 2000)
ax1.set_ylim(0, 4)
ax1.xaxis.set_major_formatter(ScalarFormatter())  # sin notación científica
ax1.yaxis.set_major_locator(MultipleLocator(0.5))
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1.errorbar(D, Av, xerr=sigD, yerr=sigAv, fmt='o')

# gráfico obscuration
x = np.arange(0, 1510, 10)
y = extin(lat, x)
ax1.plot(x, y, 'r-', lw=2)

# gráfico Av mean vs. D mean
ax2 = plt.subplot(gs[1])
ax2.set_xlabel("Distance [pc]", fontsize=14)
ax2.set_ylabel(r"mean A$\mathrm{_V}$ [mag]", fontsize=14)
ax2.text(0.15, 0.9, obj, fontsize=16, ha='center', va='center',
        transform=ax2.transAxes)
ax2.set_xlim(0, 600)
ax2.set_ylim(0, 1)
ax2.xaxis.set_minor_locator(MultipleLocator(20))
ax2.yaxis.set_major_locator(MultipleLocator(0.2))
ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
ax2.errorbar(meanAvD[:, 1], meanAvD[:, 2], yerr=meanAvD[:, 3], fmt='s-')

fig.savefig(figpng, bbox_inches='tight')
plt.show()
