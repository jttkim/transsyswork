#!/usr/bin/python

import sys
sys.path = ['/home/jtkim/lib/python'] + sys.path


import cgi

import matplotlib
import matplotlib.pyplot


fig = matplotlib.pyplot.figure()
matplotlib.pyplot.plot(range(10))

f = sys.stdout

# f.write('HTTP/1.1 200 OK\r\n')
f.write('Content-Type: image/png\r\n')
f.write('\r\n')
fig.savefig(f, format = 'png')

