#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#-----------------------------------------------------------------------#
# calc2MASSdist.py                                                      #
#                                                                       #
# Copyright (C) 2014 Germán A. Racca - <gracca[AT]gmail[DOT]com>        #
#                                                                       #
# This program is free software: you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# This program is distributed in the hope that it will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this program. If not, see <http://www.gnu.org/licenses/>.  #
#-----------------------------------------------------------------------#


"""
Ley de extinción convertida al sistema 2MASS por Nielbock & Chini (2005):

E(J-H) / E(H-K) = 2.2

Aj / E(J-H) = 2.4

Ak = 0.106 x Av

Av = 19.4834 x E(H-K)
Av =  8.8561 x E(J-H)
"""

import numpy as np

obj = 'NGC7023-F2'
file_in = 'out-2mass-' + obj
file_ou = 'AvD-' + obj
figpng  = file_ou + '.png'

# datos intrínsecos
SpT, Mk, JHin, HKin, sJHin, sHKin = np.loadtxt('intrinsic.dat',
                                               dtype=('3a,f,f,f,f,f'),
                                               unpack=True)

# datos observados
a, d, j, s_j, h, s_h, k, s_k = np.loadtxt(file_in, unpack=True)

# colores y sus errores observados
jh, hk, jk = j - h, h - k, j - k
s_jh, s_hk = np.sqrt(s_j**2 + s_h**2), np.sqrt(s_h**2 + s_k**2)

# restringimos a (J - K) <= 0.75 mag
k = k[jk<=0.75]
jh, hk = jh[jk<=0.75], hk[jk<=0.75]
s_jh, s_hk = s_jh[jk<=0.75], s_hk[jk<=0.75]
s_j, s_h, s_k = s_j[jk<=0.75], s_h[jk<=0.75], s_k[jk<=0.75]
a, d = a[jk<=0.75], d[jk<=0.75]

# función: matríz de covarianza
def covinv(sJ, sH, sK):
    """Calcula la matríz de covarianza y su inversa"""
    C = np.array([[sJ**2 + sH**2, -sH**2], [-sH**2, sH**2 + sK**2]])
    return np.linalg.inv(C)

# intervalo de extinción
ext = np.arange(0, 5.01, 0.01)

# cálculo de chi-cuadrado
nst = 0
res = np.empty((0, 7), dtype=int)
for A, D, K, JH, HK, sJ, sH, sK in zip(a, d, k, jh, hk, s_j, s_h, s_k):
    nst += 1
    print("Star %3d: (%f, %f)" % (nst, A, D))
    res1 = np.empty((0, 6), dtype=int)
    for Av in ext:
        chi = np.array([])
        for JHms, HKms in zip(JHin, HKin):
            # colores intrínsecos
            JH0 = JH - 0.1129*Av
            HK0 = HK - 0.0513*Av
            # elementos matríz de covarianza inversa
            B = covinv(sJ, sH, sK)
            B11, B12, B22 = B[0][0], B[0][1], B[1][1]
            # chi-cuadrado
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
    print("%.4f %.2f %s %.3f %.2f %.3f %.3f" % (ChiMin, AvMin, SpTMin, K, MkMin, JHinMin, HKinMin))
    res = np.append(res, [[ChiMin, AvMin, SpTMin, K, MkMin,
                           JHinMin, HKinMin]], axis=0)

# cálculo de las distancias
# d(pc) = 10^[(K - Mk + 5 - Ak)/5]
K  = res[:,3].astype(float)
Mk = res[:,4].astype(float)
Av = res[:,1].astype(float)
Ak = 0.106*Av
dist = np.power(10, (K - Mk + 5 - Ak)/5)

# cálculo de las incertezas
sigK     = s_k
sigMk    = 0.4
covJH_HK = np.corrcoef(s_jh, s_hk)[0][1]*s_jh*s_hk  # coef. de corr. de Pearson
sigAv    = np.sqrt(4.4**2*s_jh**2 + 9.7**2*s_hk**2 + 2*43.1*covJH_HK)
sigD     = 0.46*dist*np.sqrt(sigK**2 + sigMk**2 + sigAv**2)

# seleccionamos sólo ChiSq <= 0.1
ChiSq = res[:,0].astype(float)
Av    = Av[ChiSq<=0.1]
sigAv = sigAv[ChiSq<=0.1]
D     = dist[ChiSq<=0.1]
sigD  = sigD[ChiSq<=0.1]

AvD = np.vstack((Av, sigAv, D, sigD)).T
np.savetxt(file_ou, AvD, fmt='%.2f')

print("\n-- Fin --")