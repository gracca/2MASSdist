#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# data2MASSdist.py
#
# Copyright (C) 2014 Germán A. Racca
# E-Mail: <gracca[AT]gmail[DOT]com>
# License: GPLv3+

import subprocess

obj   = 'NGC7023-F2'
fname = 'out-2mass-' + obj

# ejecutamos el comando
command = ['vizquery',
           '-source=II/246',
           '-c.bm=60x60',
           #'-c.bm=30x30',
           '-out=RAJ2000 DEJ2000 Jmag Jcmsig Hmag Hcmsig Kmag Kcmsig',
           '-out.form=mini',
           '-c=314.744713 +67.160805',
           #'-c=BHR 142',
           'Jcmsig=<0.035',
           'Hcmsig=<0.035',
           'Kcmsig=<0.035',
           'Qflg=AAA']

# guardamos los datos
with open(fname, 'wb') as out:
    p = subprocess.Popen(command, stdout=out)
    p.wait()

# borramos las primeras 49 líneas
lines = open(fname).readlines()
open(fname, 'w').writelines(lines[49:-1])
