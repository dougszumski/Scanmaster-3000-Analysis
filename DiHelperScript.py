#     Doug S. Szumski  <d.s.szumski@gmail.com>  13-04-2012
#     Quick helper script for DI system data. 
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
import os
import sys

def groupavg(data):

    #Returns the average value of a string of data points
    tmp = 0.00
    for value in data:
        tmp += value
    tmp /= len(data)
    return tmp

def fileoutput(filename,data_1,data_2,data_3,data_4):
   
    #Writes quad channel data to file 
    points = len(data_1)
    with open(filename, 'w') as FILE:
        for i in range(0, points):
            FILE.write('%s' % data_1[i] + '\t %s' % data_2[i] + '\t %s'
                       % data_3[i] + '\t %s \n' % data_4[i])

def chopper(data, minSegLen, maxSegLen, l_th, u_th):

    """Chops the continuous data file into individual scans."""

    data_window = []
    # The number of points to not include from the determined 'end' (prevents problems with other channels)
    offset = 20 
    window_length = 50
    l_plat_len = 0
    u_plat_len = 0
    counter = 0
    start = 0
    stop = 0
    filecounter = 0
    output_dir = '/choppedScans'

    #Check to see if output directory exists
    if os.path.exists(os.getcwd() + output_dir) != True: 
        os.mkdir(os.getcwd() + output_dir)
        print "Created output directory"
    else:
        print "WARNING: Output directory already exists. Files will be overwritten"

    #First of all check the polarity of the STM data and invert if necessary
    dat_avg = np.average(data)
    print "LSx1 channel average is: ", dat_avg, "V"
    if (dat_avg < -0.5):
        print "Assuming negative tip-substrate bias"
        #Invert the data
        data = data * -1.0  
    elif (dat_avg > 0.5):
        print "Assuming positive tip-substrate bias"
    else:
        print "ABORTING!!! Low set point current detected; automatic polarity check failed"
        sys.exit("Failed bias check")

    #Now look through the data stream for I(s) scans
    for value in data:
        #Use a rolling data window
        data_window.append(value)
        if len(data_window) > window_length:
            # Get rid off the first value to make way for the last
            del data_window[0]
        if len(data_window) == window_length:
            # Full window so do some checks for lower threshold
            if groupavg(data_window) < l_th:
                # Found a plateau so increment the length counter
                l_plat_len += 1
            if groupavg(data_window) > u_th:
                # Found a plateau so increment the length counter
                u_plat_len += 1
            if groupavg(data_window) < u_th and u_plat_len > 0 and stop == 0:
                # Hopefully this is the tip retraction point and from here on the current decays
                # Stop must be zero otherwise we might end up with background, then decay!
                start = counter  # could make this "counter - u_plat_len" to get all the plateaus but not of interest
                # print "Found upper plateau: ", u_plat_len, " points long, stopping at: ", start
                u_plat_len = 0
            if groupavg(data_window) > l_th and l_plat_len > 0 and start != 0:
                # We found the end of the background plateau at the previous counter value
                # start must be > 0 to ensure we have already found the initial current decay
                stop = counter - offset
                # start is at the location of the last found upper plateau
                # print "Found lower plateau: ", l_plat_len, " point long, stopping at: ", stop........
                points = stop - start
                lowchan_pass = True
                background_level = groupavg(data[stop - 10:stop])
                if background_level > l_th:
                    lowchan_pass = False
                    # The level will depend on the 'characteristic level' of the experiment.
                    # For bare contacts this should of the order of mV. For in-situ work it will depend on the tip leakage current. 
                    print 'REJECTED: Voltage decayed to:', background_level, '(V). Threshold level:', l_th, '(V)'
                # Save the measurements if they pass some tests:
                # Each scan should have a characteristic length. 
                # In this case the length is controlledscan manually and will depend on the scan duration in seconds.
                # If the scan duration is changed adjust the acceptable range to the characteristic length
                # This acts as a filter to prevent multiple length scans, or unreasonably short I(s) scans getting through
                if points > minSegLen and points < maxSegLen and lowchan_pass:
                    data_slice = data[start:stop]
                    # Generate the filename
                    # Later on this can calculate a sdev and dynamically exclude
                    filename = os.getcwd() + output_dir + "/" + 'slice' + str(filecounter).zfill(4) + '.txt'
                    print 'Reconstructed I(s) scan containing: ', points, 'data points as: ', filename
                    #FIXME: Channels are duplicated to avoid modifiying the scanmaster code 
                    #If this is a useful script add feature to scanmaster
                    fileoutput(filename, data_slice, data_slice, data_slice, data_slice)
                    filecounter += 1
                l_plat_len = 0
                start = 0
                stop = 0
        counter += 1

dataStream = []
filename = ''

#TODO: Deal with sys args properly!
if (len(sys.argv) == 6):
    filename = sys.argv[1] 
    print "Input filename:", filename
    minScanLength = int(sys.argv[2])
    print "Minimum scan length:", minScanLength
    maxScanLength = int(sys.argv[3])
    print "Maximum scan length:", maxScanLength
    scanLowThresh = float(sys.argv[4])
    print "Background voltage level:", scanLowThresh
    scanUpThresh = float(sys.argv[5])
    print "Setpoint voltage level:", scanUpThresh
else:
    print "Invalid arguments\n"
    print "Usage: python DiHelperScript <FILENAME> <MIN SCAN LENGTH> <MAX SCAN LENGTH> <BACKGROUND LEVEL VOLTAGE> <SETPOINT LEVEL VOLTAGE>\n"
    print "Example: python DiHelperScript data.txt 1000 100000 0.01 9.9\n"
    quit()
    

print "Reading data file..."
with open(filename, 'r') as FILE:
    lines = FILE.readlines()
    for line in lines:
        dataStream.append(float(line))

print "Chopping data into scans..."
chopper(dataStream, minScanLength, maxScanLength, scanLowThresh, scanUpThresh)
print "Finished reconstructing scans."
