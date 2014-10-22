#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# calc2MASSdist.py
#
# Copyright (C) 2014 Germán A. Racca
# E-Mail: <gracca[AT]gmail[DOT]com>
# License: GPLv3+


"""
Extinction law converted to 2MASS system by Nielbock & Chini (2005):

E(J-H) / E(H-K) = 2.2

Aj / E(J-H) = 2.4

Ak = 0.106 x Av

Av = 19.4834 x E(H-K)
Av =  8.8561 x E(J-H)
"""

import numpy as np
import pynotify

obj = 'BHR001'
file_in = 'out-2mass-' + obj
file_ou = 'AvD-' + obj
figpng  = file_ou + '.png'

# intrinsic data
SpT, Mk, JHin, HKin, sJHin, sHKin = np.loadtxt('intrinsic.dat',
                                               dtype=('3a,f,f,f,f,f'),
                                               unpack=True)

# observed data
a, d, j, s_j, h, s_h, k, s_k = np.loadtxt(file_in, unpack=True)

# colors and their observed errors
jh, hk, jk = j - h, h - k, j - k
s_jh, s_hk = np.sqrt(s_j**2 + s_h**2), np.sqrt(s_h**2 + s_k**2)

# select only (J - K) <= 0.75 mag
k = k[jk<=0.75]
jh, hk = jh[jk<=0.75], hk[jk<=0.75]
s_jh, s_hk = s_jh[jk<=0.75], s_hk[jk<=0.75]
s_j, s_h, s_k = s_j[jk<=0.75], s_h[jk<=0.75], s_k[jk<=0.75]
a, d = a[jk<=0.75], d[jk<=0.75]

# covariance matrix
def covinv(sJ, sH, sK):
    """Calculate the covariance matrix and its inverse"""
    C = np.array([[sJ**2 + sH**2, -sH**2], [-sH**2, sH**2 + sK**2]])
    return np.linalg.inv(C)

# extinction range
ext = np.arange(0, 5.01, 0.01)

# calculation of chi-square
nst = 0
res = np.empty((0, 9), dtype=int)
print(obj + ":\n")
for A, D, K, JH, HK, sJ, sH, sK in zip(a, d, k, jh, hk, s_j, s_h, s_k):
    nst += 1
    print("Star %3d: (%f, %f)" % (nst, A, D))
    res1 = np.empty((0, 6), dtype=int)
    for Av in ext:
        chi = np.array([])
        for JHms, HKms in zip(JHin, HKin):
            # intrinsic colors
            JH0 = JH - 0.1129*Av
            HK0 = HK - 0.0513*Av
            # elements of inverse covariance matrix
            B = covinv(sJ, sH, sK)
            B11, B12, B22 = B[0][0], B[0][1], B[1][1]
            # chi-square
            ChiSq = B11*(JH0 - JHms)**2 + \
                    2*B12*(JH0 - JHms)*(HK0 - HKms) + \
                    B22*(HK0 - HKms)**2
            chi = np.append(chi, ChiSq)
        arg = chi.argmin()
        res1 = np.append(res1, [[chi.min(), Av, SpT[arg], Mk[arg],
                                 JHin[arg], HKin[arg]]], axis=0)
    ChiMin  = res1[:,0].astype(float).min()
    Min     = res1[:,0].astype(float).argmin()
    AvMin   = res1[:,1].astype(float)[Min]
    SpTMin  = res1[:,2].astype(str)[Min]
    MkMin   = res1[:,3].astype(float)[Min]
    JHinMin = res1[:,4].astype(float)[Min]
    HKinMin = res1[:,5].astype(float)[Min]
    print("%.4f %.2f %s %.3f %.2f %.3f %.3f" % (ChiMin, AvMin, SpTMin, K,
                                                MkMin, JHinMin, HKinMin))
    res = np.append(res, [[A, D, ChiMin, AvMin, SpTMin, K, MkMin,
                           JHinMin, HKinMin]], axis=0)

# calculation of distance
# d(pc) = 10^[(K - Mk + 5 - Ak)/5]
ra   = res[:,0].astype(float)
dec  = res[:,1].astype(float)
K    = res[:,5].astype(float)
Mk   = res[:,6].astype(float)
Av   = res[:,3].astype(float)
Ak   = 0.106*Av
dist = np.power(10, (K - Mk + 5 - Ak)/5)

# calculation of uncertainties
sigK     = s_k
sigMk    = 0.4
covJH_HK = np.corrcoef(s_jh, s_hk)[0][1]*s_jh*s_hk  # Pearson corr. coef.
sigAv    = np.sqrt(4.4**2*s_jh**2 + 9.7**2*s_hk**2 + 2*43.1*covJH_HK)
sigD     = 0.46*dist*np.sqrt(sigK**2 + sigMk**2 + sigAv**2)

# select only ChiSq <= 0.1
ChiSq = res[:,2].astype(float)
ra    = ra[ChiSq<=0.1]
dec   = dec[ChiSq<=0.1]
Av    = Av[ChiSq<=0.1]
sigAv = sigAv[ChiSq<=0.1]
D     = dist[ChiSq<=0.1]
sigD  = sigD[ChiSq<=0.1]

AvD = np.vstack((ra, dec, Av, sigAv, D, sigD)).T
np.savetxt(file_ou, AvD, fmt='%.6f %.6f %.2f %.2f %.2f %.2f')

print("\n" + obj + ":")
print("Number of stars after selection criteria: %3d" % a.size)
print("Number of stars with ChiSquare <= 0.1: %3d" % Av.size)
print("\n-- End --")

# send notification
pynotify.init('End')
notice = pynotify.Notification("calc2MASSdist - " + obj,
                               "Fitting finished, see terminal output.")
notice.show()

