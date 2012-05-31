#     Doug S. Szumski  <d.s.szumski@gmail.com>  13-04-2012
#     Quick helper script for MI system data. 
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

import os
import sys
import numpy as np

#Threshold level for current average used for inversion of data in the case of
#negative tip-substrate bias. You may need to adjust this if you decrease your 
#setpoint current
scanAvgThreshold = 0.1
    
#Threshold level for an acceptable background leakage. Should normally be much lower unless in echem. environment with large tip leakage.
backgroundThreshold = 1.0

def fileoutput(filename,data_1,data_2,data_3,data_4):
    #Writes quad channel data to file 
    points = len(data_1)
    with open(filename, 'w') as FILE:
        for i in range(0, points):
            FILE.write('%s' % data_1[i] + '\t %s' % data_2[i] + '\t %s'
                       % data_3[i] + '\t %s \n' % data_4[i])

def backgroundPresent(i_dat, bcorfac):
    #Check background level is below a certain threshold
    end = len(i_dat)
    final_int = int(end * bcorfac) 
    begin = int(end - final_int)
    # Calculate the average of the background
    fin_avg = sum(i_dat[begin:end]) / final_int
    # Subtract the average from every current in the I(s) measurement
    if (abs(fin_avg) < backgroundThreshold):
        return True
    else:
        return False

#Deal with the sys args
if (len(sys.argv) == 3):
    inputFolder = str(sys.argv[1]) 
    print "Input folder:", inputFolder
    outputFolder = str(sys.argv[2])
    print "Output folder:", outputFolder
else:
    print "Invalid arguments\n"
    print "Usage: python picohelper.py <INPUT FOLDER> <OUTPUT FOLDER>\n"
    quit()
    
inputPath = os.getcwd() + "/" + inputFolder
outputPath = os.getcwd() +  "/" + outputFolder
scanList = os.listdir(inputPath)

#Check to see if output directory exists
if os.path.exists(outputPath) != True: 
    os.mkdir(outputPath)
    print "Created output directory"
else:
    print "WARNING: Output directory already exists. Files will be overwritten"

#Now process the scans...
filecounter = 0
filename = "scan"
print "Converting scans..."

for scan in scanList:
    scanName = filename + str(filecounter).zfill(4) + ".txt"
    scanLoc = os.path.abspath(inputPath) + "/" + scan
    with open(scanLoc, 'r') as FILE:
        lines = FILE.readlines()
        #Trucate the scan data, this should theoretically work for any PicoScan scan length assuming the preamble doesn't change
        lines = lines[105:-1]
        i_list = []        
        for line in lines: 
            line = line.split()
            i_list.append(float(line[1]))
    #Check the polarity of the data and invert if negative bias used
    dat_avg = np.average(i_list)
    scanLength = len(i_list)
    print "Channel average is: ", dat_avg, "nA"
    if (dat_avg < -scanAvgThreshold):
        print "Negative tip-substrate bias detected; inverting measurement..."
        for i in range(scanLength):
            i_list[i] = i_list[i] * -1.0
    elif (dat_avg > scanAvgThreshold):
        print "Positive tip-substrate bias detected"
    else:
        print "WARNING: Automatic polarity check failed: skipping file."
    #Output four columns of the same data TODO: Convert scanmaster so this isn't necessary
    if backgroundPresent(i_list, 0.20):
        fileoutput((outputPath + "/" + scanName),i_list,i_list,i_list,i_list)
        print "Saved scan of length ", scanLength, " points as: ", scanName
        filecounter += 1
    else:
        print "Scan ignored: Didn't decay below background threshold"

print "Finished converting scans."
