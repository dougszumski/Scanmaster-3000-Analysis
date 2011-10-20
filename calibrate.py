#!/usr/bin/python
# -*- coding: utf-8 -*-

#     Doug S. Szumski  <d.s.szumski@gmail.com>  08-08-2011
#     Resistor calibration script for quad channel STM module
# 
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program; if not, write to the Free Software
#     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

import numpy as np
import matplotlib.pylab as plt
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from pylab import plot, title, show , legend
import sys

test_res = 50808 #Ohms from Keithly source meter, on high accuracy mode (stable measurement)

def fileinput(filename):
	#Read data from input file
	infile=open(filename,"r")
	ch0 = []
	ch1 = []
	ch2 = []
	ch3 = []
	lines = infile.readlines()
	for line in lines:
		line = line.split()
		ch0.append(float(line[0]))	
		ch1.append(float(line[1]))	
		ch2.append(float(line[2]))	
		ch3.append(float(line[3]))	
	#print float(line[1]), float(line[2]), float(line[3]), float(line[4])
	return ch0, ch1, ch2, ch3
	infile.close

def slicer(data, minval, maxval):
	slice_begin = 0
	for i in range (0,len(data)):
		if ((data[i] > minval) and (slice_begin == 0)):
			slice_begin = i
		if ((data[i] > maxval) and (slice_begin != 0)): 
			slice_end = i
			return slice_begin, slice_end


#Recommended sweep voltage of 1V for 50k resistor to generate data file
#Record the file using the DAQ software
#Be aware that there is about 2mV peak to peak 50Hz interference on the sample voltage
#This was checked on both the Picoscan and the newer Agilent system. Bit of a joke, 
#but could be removed with a filter. 
#Depending on the sweep rate the noise will probably show up on the HS channels.

#This is the data file we are reading in. 
ch0,ch1,ch2,ch3 = fileinput("output_10V_new.txt")

maxloc = ch0.index(max(ch0)) #return the location of the maximum value
minloc = ch0.index(min(ch0))

#TODO: If the data is non-contiguous this won't work and you'll get an error. 
#Basically look through the calibration data and delete the bit after the sweep.

if (minloc > maxloc):
	print "ERROR: Calibration data is non-contiguous. Manually remove 'jump' to continue."
	sys.exit(0)

ch0 = ch0[minloc:maxloc]
ch1 = ch1[minloc:maxloc]
ch2 = ch2[minloc:maxloc]
ch3 = ch3[minloc:maxloc]
x = np.linspace(-1,1,len(ch0))

#On this channel we take the whole slice of data spanning the range of the DAQ
slice_begin, slice_end = slicer(ch1, -10.0, 10.0)
ch1_ss = ch1[slice_begin:slice_end]
x1_ss = x[slice_begin:slice_end]

#On this channel the limter switches on at about 8-9V so look at the section
#when the limiter is off
slice_begin, slice_end = slicer(ch2, -8.0, 8.0)
ch2_ss = ch2[slice_begin:slice_end]
x2_ss = x[slice_begin:slice_end]

#This is the section from the above set with the limiter on
slice_begin, slice_end = slicer(ch2, -10.4, -9.0)
ch2_ss2 = ch2[slice_begin:slice_end]
x2_ss2 = x[slice_begin:slice_end]

#Same again with the limiter
slice_begin, slice_end = slicer(ch3, -8.0, 8.0)
ch3_ss = ch3[slice_begin:slice_end]
x3_ss = x[slice_begin:slice_end]

(ar,br)=polyfit(x,ch0,1)
xr0=polyval([ar,br],x)
Rf = ar*test_res
print "Channel 0 gradient is: ", ar, " and intercept: ", br
print  "Effective resistance (ohms): ", Rf

(ar,br)=polyfit(x1_ss,ch1_ss,1)
xr1=polyval([ar,br],x1_ss)
Rf = ar*test_res
print "Channel 1 gradient is: ", ar, " and intercept: ", br
print  "Effective resistance (ohms): ", Rf

(ar,br)=polyfit(x2_ss,ch2_ss,1)
xr2=polyval([ar,br],x2_ss)
Rf = ar*test_res
print "Channel 2, limiter off gradient is: ", ar, " and intercept: ", br
print  "Effective resistance (ohms): ", Rf

(ar,br)=polyfit(x2_ss2,ch2_ss2,1)
xr22=polyval([ar,br],x2_ss2)
Rf = ar*test_res
print "Channel 2, limiter on gradient is: ", ar, " and intercept: ", br
print  "Effective resistance (ohms): ", Rf

(ar,br)=polyfit(x3_ss,ch3_ss,1)
xr3=polyval([ar,br],x3_ss)
Rf = ar*test_res
print "Channel 3 gradient is: ", ar, " and intercept: ", br
print  "Effective resistance (ohms): ", Rf


xr3_theory=polyval([ar,br],x3_ss)

#matplotlib ploting
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Quad channel STM calibration')
ax.plot(x,ch0,'k.', label = "Ch0: 10uA/V")
ax.plot(x,xr0, 'r-', label = "Ch0: polyfit")
ax.plot(x,ch1,'r.', label = "Ch1: 1uA/V")
ax.plot(x1_ss,xr1, 'b-', label = "Ch1: polyfit")
ax.plot(x,ch2,'g.', label = "Ch2: 10nA/V")
ax.plot(x2_ss,xr2, 'k-', label = "Ch2: polyfit: limiter off")
ax.plot(x2_ss2,xr22, 'c-', label = "Ch2: polyfit: limiter on")
ax.plot(x,ch3,'b.', label = "Ch3: 1nA/V")
ax.plot(x3_ss,xr3, 'g-', label = "Ch3: polyfit")
ax.grid(True)
ax.legend(loc="upper left")
#ax.set_xlim([-0.2,.2])
ax.set_xlabel('Input voltage (V)')
ax.set_ylabel('Output voltage (V)')

plt.show()

