#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# data2MASSdist.py
#
# Copyright (C) 2014 Germán A. Racca
# E-Mail: <gracca[AT]gmail[DOT]com>
# License: GPLv3+

import subprocess

obj   = 'BHR001'
fname = 'out-2mass-' + obj
lobj  = '-c=%s' % obj

# command to execute in bash
command = ['vizquery',
           '-source=II/246',
           '-c.bm=45x45',
           '-out=RAJ2000 DEJ2000 Jmag Jcmsig Hmag Hcmsig Kmag Kcmsig',
           '-out.form=mini',
           lobj,
           'Jcmsig=<0.03',
           'Hcmsig=<0.03',
           'Kcmsig=<0.03',
           'Qflg=AAA']

# save data
with open(fname, 'wb') as out:
    p = subprocess.Popen(command, stdout=out)
    p.wait()

# erease first 49 lines
lines = open(fname).readlines()
open(fname, 'w').writelines(lines[49:-1])

print("\nData for " + obj + " is ready!\n")
