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

import os
import sys

def fileoutput(filename,data_1,data_2,data_3,data_4):
    #Writes quad channel data to file 
    points = len(data_1)
    with open(filename, 'w') as FILE:
        for i in range(0, points):
            FILE.write('%s' % data_1[i] + '\t %s' % data_2[i] + '\t %s'
                       % data_3[i] + '\t %s \n' % data_4[i])

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

filecounter = 0
filename = "scan"
print "Converting scans..."

for scan in scanList:
    scanname = filename + str(filecounter).zfill(4) + ".txt"
    scanLoc = os.path.abspath(inputPath) + "/" + scan
    with open(scanLoc, 'r') as FILE:
        lines = FILE.readlines()
        #Trucate the scan data, this should theoretically work for any PicoScan scan length assuming the preamble doesn't change
        lines = lines[105:-1]
        i_list = []        
        for line in lines: 
            line = line.split()
            i_list.append(float(line[1]))
    #Output four columns of the same data TODO: Convert scanmaster so this isn't necessary
    fileoutput((outputPath + "/" + scanname),i_list,i_list,i_list,i_list)
    print "Saved scan as: ", scanname
    filecounter += 1

#TODO: Parameterise "scale" in scanmaster so it can be set to 1 instead of 1e9 -- then scans should load in scanmaster
#TODO: Check for output directory and create a new one if it doens't exist. 


print "Finished converting scans."
