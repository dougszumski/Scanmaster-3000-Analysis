#!/usr/bin/python
# -*- coding: utf-8 -*-

#     Doug S. Szumski  <d.s.szumski@gmail.com>  14-06-2012
#     with contributions from Richard Brooke
#
#     Processing and analysis script for quad channel STM module
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

#Updates are posted at: https://github.com/dougszumski/Scanmaster-3000
#Git repo: git://github.com/dougszumski/Scanmaster-3000.git
#For push access contact Doug

# For embedded graph:
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
    NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import sys
from matplotlib.colors import LogNorm

#Multithreading support
import threading
from Queue import Queue

# Dialogue
from Tkinter import *
from tkFileDialog import askopenfilename, asksaveasfile
from tkColorChooser import askcolor
from tkMessageBox import askquestion, showerror, showinfo
from tkSimpleDialog import askfloat
from tkMessageBox import askokcancel
from FileDialog import LoadFileDialog
from Dialog import Dialog
import signal

# File IO
import shutil  # High level file from scipy import stats
import os  # import command line for file operations
from subprocess import call
import string
import cPickle
import gzip

# Plotting stuff
import numpy as np
from scipy import stats as sci_stats
import matplotlib.pylab as plt
import statistics
from scipy import linspace, polyval, polyfit, sqrt, stats, randn

#2D corr stuff, FIXME: is some of this pylab stuff redundant?
from pylab import imshow, plt, contour, xlabel, ylabel, title, figure, barh, bar
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.mlab as mlab
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Debugging
import random
import time

# Retro Fortran modules -- recompile using: f2py -c -m plat_seek plat_seek8.f90 etc..
import plat_seek

class Data:
    #Initialise a load of lists / variables for global use
    #TODO Add functions to set/read these attributes rather than addressing the directly (good practice)
    def __init__(self):
        #G0 in nS 
        self.G0 = 77480.917 
        # Sampling distance interval - calculated in scaninput for each file
        self.interval = 0.00
        # Four current channels from the quad amp
        self.i_dat_all_ls = []
        self.i_dat_all_hs = []
        self.i10_dat_all_ls = []
        self.i10_dat_all_hs = []
        # The recombined current data
        self.i_dat_all_combined = []
        # Reconstructed spatial array to go with all current measurements
        self.s_dat_all_ls = []
        # The current working filename
        self.filename = []
        # The current file output number (for numbering selected data saved to new file)
        self.file_newnum = 0
        # Filtered data used in linear regression fitting
        self.i_dat_filtered = []
        self.s_dat_filtered = []
        # Linear regression related stuff for plotting
        self.polyfit_rescaled = []
        self.err = 0.00  # Error on least squares fit for linear reg.
    def clear_current_lists(self):
        # Four current channels from the quad amp
        self.i_dat_all_ls = []
        self.i_dat_all_hs = []
        self.i10_dat_all_ls = []
        self.i10_dat_all_hs = []
        # The recombined current data
        self.i_dat_all_combined = []
        # Reconstructed spatial array to go with all current measurements
        self.s_dat_all_ls = []
        print "Current lists cleared"
        
def scanInput(name):

    # Reads in the file, filename, reads number of measurements and then extracts I(s) data into lists
    # Read data from input file

    with open(name, 'r') as infile:
        s_dat_ls = []
        i_dat_ls = []
        i_dat_hs = []
        i10_dat_ls = []
        i10_dat_hs = []
        lines = infile.readlines()

        # Reconstruct the x-axis using variables from the menu
        sampsec = controller.sampsec.get()
        srange = controller.srange.get()
        sduration = controller.sduration.get()

        # Total distance travelled by tip in data segment =  (retraction rate * points in measurment) / sampling rate
        # This is the distance the tip travels between data points only
        data.interval = srange / sduration / sampsec
        position = 0.00
        for line in lines:
            line = line.split()
            i_dat_ls.append((float(line[0])))  
            i10_dat_ls.append((float(line[1])))
            i_dat_hs.append((float(line[2])))
            i10_dat_hs.append((float(line[3])))
            s_dat_ls.append(position)
            position += data.interval
        return (i_dat_ls, i10_dat_ls, i_dat_hs, i10_dat_hs, s_dat_ls)
 
def contourPlot(savefig=False):
    """ Plots 2d histograms of current distance scans """

    # Generate list of all distances. This list will be the same length as the current list and therefore
    # the indices will be directly related.
    s_list = []
    for s_dat in data.s_dat_all_ls:
        for value in s_dat:
            s_list.append(value)

    # Generate list of all currents in the current distance scans
    offset = float(controller.offset.get())
    i_list = []
    for i_dat in data.i_dat_all_combined:
        for value in i_dat:
            #FIXME: abs shouldn't be used here because of data inversion. The offset
            #should avoid problems if big enough, and abs should only invert wild points
            if value > 0.0:
                i_list.append(np.log10(value+offset))
            else:   
            #FIXME Shift negative data 'under the bed' for now
                i_list.append(np.log10(1e-9+offset))
                
            #i_list.append(np.log10(abs(value+offset)))
            #i_list.append(value) #uncomment for linear plot
    if (len(i_list) < 1):
        error = showerror('Error', 'No data in memory')
        return

    # Plot the 2D histogram
    (H, xedges, yedges) = np.histogram2d(i_list, s_list,
            bins=(controller.ybin_contour.get(),
            controller.xbin_contour.get()),
            range=[[controller.ymin_contour.get(),
            controller.ymax_contour.get()],
            [controller.xmin_contour.get(),
            controller.xmax_contour.get()]], normed=True)
    (H.shape, xedges.shape, yedges.shape)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]  # Don't forget -1 means the last item in the list!
    
    plt.imshow(H, origin='lower', extent=extent, interpolation='nearest', norm=LogNorm(vmin=0.01, vmax=0.5))
    title = data.filename[0:string.rfind(data.filename, '/')]
    plt.title(title, fontsize=8)
    plt.xlabel('Distance (nm)')
    plt.ylabel('log_10[current (nA)]')
    cb = plt.colorbar(ticks=[0.01, 0.1, 0.5])
    cb.set_ticklabels([0.01,0.1,0.5])
    cb.set_label('log_10[counts] (normalised)')
    if savefig:
        #TODO: FIX ISSUE WITH COLORBAR
        plt.savefig(title+"_2d.png", format='png')
        plt.close()
    else:
        plt.show()

def correlationHist(savefig=False):
    """Plots a 2D correlation histogram as per Makk et al."""
    
    #TODO: Turn this into a module?

    #Histogram settings from GUI
    numBins = controller.corrbins.get()
    logscale = controller.corrlogscale.get()
    bias = controller.bias.get()
    if logscale: 
        print "Plotting 2D correlation histogram with log scale"
        histLowerLim = controller.corrcurrentmin.get()
        histUpperLim =  controller.corrcurrentmax.get()
    else: 
        print "Plotting 2D correlation histogram with linear scale"
        histLowerLim = 10**controller.corrcurrentmin.get()
        histUpperLim = 10**controller.corrcurrentmax.get()
    print "Loading data..."

    #Create a matrix containing all the scan data 
    #TODO: This is redundant and a waste of resources, convert scan container to a matrix in the data class?
    numScans = len(data.i_dat_all_combined)
    scanMatrix = np.zeros(numScans*numBins).reshape(numScans,numBins)
    for i in range(numScans):
        if logscale:
            hist, binEdges = np.histogram(np.log10(data.i_dat_all_combined[i]), numBins, (histLowerLim,histUpperLim))   
        else:
            hist, binEdges = np.histogram(data.i_dat_all_combined[i], numBins, (histLowerLim,histUpperLim))   
        scanMatrix[i,:] = hist

    print "Computing 2D correlation histogram..."
    #Scan histogram
    scanHist = np.zeros(numBins)
    #Sum the scans together
    for i in range(numScans):
        scanHist += scanMatrix[i,:]
    #Normalise
    scanHist /= numScans

    #Cross product matrix: Product of individual scans averaged over all scans
    cpMatrix = np.zeros(numBins*numBins).reshape(numBins,numBins)
    for i in range(numScans):
        cpMatrix += scanMatrix[i,:] * scanMatrix[i,:].reshape(numBins,1)
    cpMatrix /= numScans

    #Covarience matrix (normalised by default)
    covMatrix = cpMatrix[:] - scanHist[:] * scanHist[:].reshape(numBins,1)

    #Correlation matrix
    devMeanVal = np.zeros(numBins) 
    for i in range(numScans):
        devMeanVal += (scanMatrix[i,:] - scanHist[:])**2
    normFactor = ( devMeanVal[:] * devMeanVal[:].reshape(numBins,1) )**0.5 
    #Normalise the normalisation factor!
    normFactor /= numScans
    corrMatrix = covMatrix / normFactor

    #Colour map to match ACSNANO paper:
    acsnanoMap = {'red':   ((0.0, 0.0, 0.0),
                            (0.52, 0.0, 1.0),
                            (1.0, 1.0, 1.0)),

                 'green':  ((0.0, 0.0, 0.0),
                            (0.48, 0.0, 0.5),
                            (0.50, 0.5, 0.0),
                            (0.50, 0.0, 0.5),
                            (0.52, 0.5, 0.0),
                            (0.52, 0.0, 1.0),
                            (1.0, 0.0, 0.0)),

                 'blue':   ((0.0, 0.0, 0.0), 
                            (0.48, 1.0, 0.0),
                            (1.0, 0.0, 0.0)) 
                }
    blackMap =   {'red':   ((0.0, 0.0, 0.0),
                            (1.0, 0.0, 0.0)),

                 'green':  ((0.0, 0.0, 0.0),
                            (1.0, 0.0, 0.0)),

                 'blue':   ((0.0, 0.0, 0.0), 
                            (1.0, 0.0, 0.0)) 
                }
    acsnano = LinearSegmentedColormap('acsNanoMap', acsnanoMap)
    black = LinearSegmentedColormap('blackMap', blackMap)

    #Plot the correlation matrix
    print "Plotting 2D correlation histogram..."
    fig = plt.figure(1, figsize=(10.0,10.0))
    corrPlot = plt.subplot(111)
    extent = [histLowerLim,histUpperLim,histLowerLim,histUpperLim] 
    cax = corrPlot.imshow(corrMatrix, origin='lower', extent=extent, cmap=acsnano, vmax = 1.0, vmin = -1.0, interpolation='nearest')
    cb = plt.colorbar(cax, shrink=0.75)
    cb.set_label(r'$H_{i,j}^{corr}$', fontsize=18)
    corrPlot.contour(corrMatrix, origin='lower', extent=extent, cmap=black, vmax = 1.0, vmin = -1.0)
    corrPlot.set_aspect(1.)
    if logscale:
        corrPlot.set_xlabel(r'$log_{10}(G/G_{0}$)', fontsize =18)
        corrPlot.set_ylabel(r'$log_{10}(G/G_{0}$)', fontsize =18)
    else:
        corrPlot.set_xlabel(r'$G/G_{0}$', fontsize =18)
        corrPlot.set_ylabel(r'$G/G_{0}$', fontsize =18)
    corrPlot.grid(True, color='w', linestyle='-', which='major', linewidth=1)

    #Setup plot for histograms on the side
    divider = make_axes_locatable(corrPlot)
    axHisty = divider.append_axes('right', 1.5, pad=0.1, sharey=corrPlot)
    plt.setp(axHisty.get_yticklabels(), visible=False)
    plt.setp(axHisty.get_xticklabels(), visible=False)

    #Generate the current histogram on the side
    i_list = []
    for i in range(len(data.i_dat_all_combined)):
        for value in data.i_dat_all_combined[i]:
            if logscale:
                value = np.log10(value)
                if (value > histLowerLim and value < histUpperLim):
                    i_list.append(value)
            else:
                if (value > histLowerLim and value < histUpperLim):
                    i_list.append(value)
    axHisty.hist(i_list, bins = numBins, orientation='horizontal', normed = 1)
    #Plot the histogram on the top
    axHistx = divider.append_axes('top', 1.5, pad=0.1 , sharex=corrPlot)
    plt.setp(axHistx.get_yticklabels(), visible=False)
    plt.setp(axHistx.get_xticklabels(), visible=False)
    axHistx.hist(i_list, bins = numBins, normed = 1)    

    #TODO Add auto y_max on linear histograms? Add auto rescale if plot partially blank?

    #This is the alternative and more efficient way, since the data is already
    #binned, but unresolved problems with sharing the axes. 
    #
    #for i in range(len(scanHist)):
    #    axHisty = barh(i, scanHist[i], 1.0)

    title = data.filename[0:string.rfind(data.filename, '/')]
    plt.title(title, fontsize=10)
    plt.draw()
    if savefig:
        plt.savefig(title+"_ch.png", format='png')
        plt.close()
    else:
        plt.show()

def export_current_data():
    # Writes current array to file for external analysis
    # TODO: Could be useful to have the option to choose:
    #       log of current, sorted current...
    
    #Start by assembling the current list
    i_list = []
    for i_dat in data.i_dat_all_combined:
        for value in i_dat:
            i_list.append(value)

    #Check it's got data in, otherwise flag an error
    if (len(i_list) < 1):
        error = showerror('Error', 'No data in memory')
        return

    #Write it to file 
    with asksaveasfile(mode='w', title='Save current list as...') as FILE:
        if FILE:
            for value in i_list:
                FILE.write('%s \n' % value)
            info = showinfo('Notice', 'File saved')

def export_scan_data():
    """Exports scans stored in memory"""

    if (len(data.i_dat_all_combined) < 1):
        error = showerror('Error', 'No data in memory')
        return
    #Save current directory (should be at script level)
    original_directory = os.getcwd()
    #Make new directoy to work in:
    working_dir = "exportedScans"
    if os.path.exists(working_dir) != True:
        os.mkdir(working_dir)
    try:
        listLen = len(data.i_dat_all_combined)
        os.chdir(os.path.abspath(working_dir))
        print "Saving: ", listLen, " scans"
        for i in range(listLen):
            with open("combined" + str(i).zfill(4) + ".txt", 'w') as FILE:
                for value in data.i_dat_all_combined[i]:
                    FILE.write('%s \n' % value)
    finally:
        os.chdir(original_directory)
        info = showinfo('Notice', str(listLen) +' scans saved to:\n ../' + working_dir)
     
def plat_seeker(current_trace):

    # Search data in current_trace for plateaus
    num_data_points = len(current_trace)
    plat_max_grad = controller.plat_max_grad.get()
    background_tol = controller.background_tol.get()
    max_plat_cur = controller.max_plat_cur.get()
    max_points_plat = controller.max_points_plat.get()
    fractional_plat_tol = controller.fractional_plat_tol.get()
    min_points_plat = controller.min_points_plat.get()

    # Attempt to fit plateaus of decreasing length to current trace using 'legacy' Fortran code for speed
    (p_dat, p_avg, p_crd) = plat_seek.s_dev(
        max_points_plat,
        min_points_plat,
        plat_max_grad,
        background_tol,
        fractional_plat_tol,
        max_plat_cur,
        data.interval,
        current_trace,
        num_data_points,
        )
    #TODO Add this in the fortran module
    # Did we find a plateau? - more efficient way to do this is in the Fortran module
    if sum(p_avg) > 0.00:
        p_loc = True
    else:
        p_loc = False

    # Return list of length current_trace, but with only the plateau average non-zero (if found)
    return (p_avg, p_loc, p_crd)


def datInput(start, finish, bcorfac):

    # Generates sequential filenames, currently in the format filenameXXXX.ivs where the number XXXX is set in the GUI.
    # For each filename it then calls scaninput (see above), corrects the current data for a background offset
    # by averaging over the final fraction of the data set in the GUI, then appends the I(s) data to a list

    file_prefix = data.filename[0:-8]  # includes file path 4 digits
    file_ext = '.txt'
    file_number = data.filename[-8:-4]  # 4 digits
    print 'Reading files...'
    for i in range(start, finish + 1):
        # Generate filename
        auto_name = file_prefix + str(i).zfill(4) + file_ext  # 4 digits
        # auto_name = file_prefix + str(i).zfill(3) + file_ext #3digits
        print 'Reading...', auto_name
        #File prefix for filtered saved directory, extracted from the current directory
        temp = auto_name[0:string.rfind(auto_name, '/')]
        savedir = temp[string.rfind(temp, '/')+1:len(temp)]
        # Read in data from file
        (i_dat_ls, i10_dat_ls, i_dat_hs, i10_dat_hs, s_dat_ls) = \
            scanInput(auto_name)
        # Work out where to chop the channels whilst the data is stored as a voltage
        # Threshold set from the GUI: Normally these are 0.1V 
        th_1 = controller.th_1.get()
        th_2 = controller.th_2.get()
        th_3 = controller.th_3.get() 
        bias = controller.bias.get()
        gainfactor = controller.gainfactor.get()
        #Now have a look at the background levels on the HS channels to see if they've saturated
        end = len(i_dat_hs)
        final_int = int(end * bcorfac)
        begin = int(end - final_int)
        # Calculate the absolute average of the background 
        fin_avg_hs1 = abs(sum(i_dat_hs[begin:end]) / final_int)
        fin_avg_hs10 = abs(sum(i10_dat_hs[begin:end]) / final_int)
        #Now see if the channels have saturated
        sat_flag_hs1 = False
        sat_flag_hs10 = False     
        #TODO Add the 6.0V here as a parameter in the GUI. It is a function of the limiter circuit.   
        if (fin_avg_hs1 > 6.0):
            print "WARNING: Channels HSx1 and HSx10 are saturated and have been ignored"
            sat_flag_hs1 = True
            #It follows that HSx10 must have also saturated
            sat_flag_hs10 = True
        elif (fin_avg_hs10 > 6.0):
            sat_flag_hs10 = True
        #Get the resistor values
        i_ls_res = int(controller.lowres.get())
        i_hs_res = int(controller.highres.get())
        scale = controller.currentScaleFactor.get()  # convert to nanoamps - normally 1e9
        #Make a copy of the list (note to self: direct assignment creates a pointer to the list which caused minor hair loss)
        #Important: This contains voltages so they can be compared directly to the thresholds
        i_dat_combined = list(i_dat_ls)
        #Set this to true to use only the LSx1 channel
        single_chan = False
        if single_chan:
            #Test for only one channel
            for i in range(len(i_dat_combined)):
                #if i_dat_combined[i] 
                i_dat_combined[i] = i_dat_combined[i] / (i_ls_res) * scale / (bias*data.G0) 
        else:
            #Stitch the channels together, using the least sensitive non-saturated channel (top down approach)
            for i in range(len(i_dat_combined)): 
                if (abs(i_dat_combined[i]) < th_1):
                    #Output below threshold so try and replace it with a higher sensitivity measurement
                    if (abs(i10_dat_ls[i]) > th_2):
                        #LSx10 channel is still in operation so use the measurment from there
                        i_dat_combined[i] = i10_dat_ls[i] / (i_ls_res * gainfactor) * scale / (bias*data.G0) 
                    elif ((abs(i_dat_hs[i]) > th_3) and not sat_flag_hs1):
                        #Limiter is off, HSx1 channel is in operation
                        i_dat_combined[i] = i_dat_hs[i] / i_hs_res * scale / (bias*data.G0) 
                    elif not sat_flag_hs10:
                        #If HSx1 is below threshold then use the HSx10 channel regardless
                        i_dat_combined[i] = i10_dat_hs[i] / (i_hs_res * gainfactor) * scale / (bias*data.G0) 
                    elif not sat_flag_hs1:
                        #If the above fails because HSx10 is saturated fall back to HSx1 as that's the best we have
                        i_dat_combined[i] = i_dat_hs[i] / i_hs_res * scale / (bias*data.G0) 
                    else:
                        #If both HSx1 and HSx10 channels are saturated then use LSx10
                        #NOTE: In this case the measurement is probably useless for anything but a break junction
                        i_dat_combined[i] = i10_dat_ls[i] / (i_ls_res * gainfactor) * scale / (bias*data.G0) 
                else:
                    #LSx1 channel is fine
                    i_dat_combined[i] = i_dat_combined[i] / (i_ls_res) * scale / (bias*data.G0) 

        #Now convert the individual channels, could have convoluted this with the above
        #for efficiency at the expense of clarity
        #This isn't required really, but its nice to have on the graph plots for fine tuning
        for i in range(len(i_dat_ls)):
                i_dat_ls[i] = i_dat_ls[i] / i_ls_res * scale / (bias*data.G0) 
                i10_dat_ls[i] = i10_dat_ls[i] / (i_ls_res * gainfactor) * scale / (bias*data.G0) 
                #TODO (low priority) If these are saturated they will appear as almost flat lines on the plot
                #Notify the plotter not to plot them to avoid confusion?
                if not sat_flag_hs1: i_dat_hs[i] = i_dat_hs[i] / i_hs_res * scale / (bias*data.G0) 
                if not sat_flag_hs10: i10_dat_hs[i] = i10_dat_hs[i] / (i_hs_res * gainfactor) * scale / (bias*data.G0) 
        
        #Correct the background: this is the important one
        i_dat_combined = back_correct(i_dat_combined, bcorfac)
        #This is done so that a negative decay level doesn't get removed from the scatter plot
        #Really, these should all decay to the leakage of the opamp, but not in the case of in-situ STM
        i_dat_ls = back_correct(i_dat_ls, bcorfac)
        i10_dat_ls = back_correct(i10_dat_ls, bcorfac)
        if not sat_flag_hs1: i_dat_hs = back_correct(i_dat_hs, bcorfac)
        if not sat_flag_hs10:  i10_dat_hs = back_correct(i10_dat_hs, bcorfac)

        # Append individual data calculations to list
        data.i_dat_all_ls.append(i_dat_ls)
        data.i10_dat_all_ls.append(i10_dat_ls)
        data.i_dat_all_hs.append(i_dat_hs)
        data.i10_dat_all_hs.append(i10_dat_hs)
        data.s_dat_all_ls.append(s_dat_ls)
        data.i_dat_all_combined.append(i_dat_combined)

        # Fit a linear regression line to data for autofiltering and plotting if required
        #TODO: Each of these tests saves the files in a folder, but if the folder already contains files 
        # They either must be deleted or not saved at all to avoid 'contamination'
        if controller.autocheck_linfit.get():
            data.i_dat_filtered = []
            data.s_dat_filtered = []
            for i in range(0, len(i_dat_combined)):
                if i_dat_combined[i] < 5000 and i_dat_combined[i] \
                    > 0.01:
                    data.i_dat_filtered.append(np.log10(i_dat_combined[i]))
                    data.s_dat_filtered.append(s_dat_ls[i])
            (ar, br) = polyfit(data.s_dat_filtered,
                               data.i_dat_filtered, 1)
            xr = polyval([ar, br], data.s_dat_filtered)
            data.polyfit_rescaled = []
            # A new low in the art of bodgery, because we're plotting on a log scale take the inverse log first!
            # TODO: Improve efficiency by getting rid of this bodge
            for value in xr:
                data.polyfit_rescaled.append(10 ** value)
            data.err = sqrt(sum((xr - data.i_dat_filtered) ** 2)
                            / len(data.i_dat_filtered))
            # Save the data if MSE less than a certain amount, put this in the GUI if it works...
            if data.err < 0.4:
                print 'PASSED: polyfit test', data.err
                autosave(auto_name, savedir + '_fil_linfit')

        if controller.autocheck_avgcur.get():
            avgcur = 0.00
            for i in i_dat_combined:
                avgcur += np.log10(abs(i))
            if avgcur < controller.datfillogi.get():
                print 'PASSED: average log(I) check', data.err
                autosave(auto_name, savedir + '_fil_avgcur')
            else:
                print 'FAILED: average log(I) check with value:', avgcur

        if controller.autocheck_pltfit.get():
            # Have a look for some plateaus to highlight; you never know, there might actually be some!
            (plat_avg, plat_loc, plat_crd) = plat_seeker(i_dat_ls)
            if plat_loc:
                print 'PASSED: plateaus were found'
                autosave(auto_name, savedir + '_fil_pltfit')
            else:
                print 'FAILED: no plateaus found'

        # Plot/check if set in the GUI
        if controller.plot_dat.get() > 0:
            # If we haven't already looked for plateaus then we better go do it now so we can plot them on the graph
            if controller.autocheck_pltfit.get() == False:
                plat_avg, plat_loc, plat_crd = plat_seeker(i_dat_ls)
            egraph.updater(
                s_dat_ls,
                i_dat_ls,
                i10_dat_ls,
                i_dat_hs,
                i10_dat_hs,
                i_dat_combined,
                plat_avg,
                auto_name,
                )
        if controller.check_dat.get() > 0:
            userinput(auto_name)

def plat_sync():
    """Synchronises I(s) scans from G0 plateaus"""
    filecounter = 0
    #FIXME: these are used to convert the data from current back to volts - not ideal 
    #Other options are: call plat_sync from dat input, create copy of voltage arrays etc..    
    i_ls_res = int(controller.lowres.get())
    i_hs_res = int(controller.highres.get())
    gainfactor = controller.gainfactor.get()
    scale = controller.currentScaleFactor.get()  # convert to nanoamps
    #Save current directory (should be at script level)
    original_directory = os.getcwd()
    #Make new directoy to work in:
    working_dir = data.filename[0:string.rfind(data.filename, '/')]+"_psf"
    print "FILES SAVED HERE", working_dir
    if os.path.exists(working_dir) != True:
        os.mkdir(working_dir)
    try:
        os.chdir(os.path.abspath(working_dir))
        for i in range(len(data.i_dat_all_combined)):
            plat_avg, plat_loc, plat_crd = plat_seeker(data.i_dat_all_ls[i])
            n_plats = max(plat_crd)
            print "Found", n_plats, "plateau(s)"
            plat_start = 0
            plat_end = 0
            save_flag = False
            for j in range(len(plat_avg)):
                #TODO ADD the range to the GUI
                if ( (plat_avg[j] > 5000.0) and (plat_avg[j] < 10000) and (plat_start == 0)): 
                    #We've found the start of the plateau
                    plat_num = plat_crd[j]
                    print "Plateau number", plat_num, "is in range, with value", plat_avg[j]
                    plat_start = j
                if ((plat_start > 0) and (plat_crd[j] != plat_num)):
                    #We've found the end of the plateau
                    plat_end = j-1
                    #Ignore others for now (in a perfect world there won't be any if the range is tight)
                    save_flag = True
                    break
            if save_flag == True:
                filename = 'psync' + str(filecounter).zfill(4) + '.txt'
                #FIXME: Converting back from currents to voltages, only to reconvert is poor. If this works, clean it up.
                for j in range(len(data.i_dat_all_ls[i])):
                    data.i_dat_all_ls[i][j] = data.i_dat_all_ls[i][j] * i_ls_res / scale
                    data.i10_dat_all_ls[i][j] = data.i10_dat_all_ls[i][j] * (i_ls_res * gainfactor) / scale
                    data.i_dat_all_hs[i][j] = data.i_dat_all_hs[i][j] * i_hs_res / scale
                    data.i10_dat_all_hs[i][j] = data.i10_dat_all_hs[i][j] * (i_hs_res * gainfactor) / scale 
                #TODO Add the choice of start /end of plateau syncing
                fileoutput(filename, 
                        data.i_dat_all_ls[i][plat_start:],
                        data.i10_dat_all_ls[i][plat_start:],
                        data.i_dat_all_hs[i][plat_start:],
                        data.i10_dat_all_hs[i][plat_start:],
                        )
                print "File saved as:", filename
                filecounter += 1
    finally:
        os.chdir(original_directory)

def back_correct(i_dat, bcorfac):

    # Correct data for background offset
    # Define begin and end for background average
    end = len(i_dat)
    final_int = int(end * bcorfac)  # add correction factor later
    begin = int(end - final_int)

    # Calculate the average of the background
    fin_avg = sum(i_dat[begin:end]) / final_int

    # Subtract the average from every current in the I(s) measurement
    for j in range(end):
        i_dat[j] = i_dat[j] - fin_avg
    return i_dat


def userinput(auto_name):
    # Used to generate save/ignore option at the command line for rejecting or keeping data
    # If the user saves the data it is copied to the folder 'processed' which is created in the path of the script
    targetfile = auto_name  # Local copy for recursion
    var = raw_input('Save(s)/Ignore(d)? : ')
    if var == 's':
        print 'Saving data...'
        data.file_newnum += 1
        file_ext = '.txt'
        new_name = 'man' + str(data.file_newnum).zfill(4) + file_ext
        dirname = 'processed'
        pathvar = './' + dirname + '/'  # Copy data to this directory
        if not os.path.isdir(pathvar):  # Create directory if not exist
            os.mkdir(pathvar)
        shutil.copy(targetfile, pathvar + new_name)  # Copy file to pathvar
    elif var == 'd':
        print 'Data ignored'
    else:
        print 'Invalid input, try again'
        userinput(targetfile)


def autosave(auto_name, dirname):
    # Copy a file auto_name in a directory with no questions asked
    targetfile = auto_name  # Local copy for recursion
    print 'Saving data...'
    data.file_newnum += 1
    file_ext = '.txt'
    new_name = 'man' + str(data.file_newnum).zfill(4) + file_ext
    pathvar = './' + dirname + '/'  # Copy data to this directory
    if not os.path.isdir(pathvar):  # Create directory if not exist
        os.mkdir(pathvar)
    shutil.copy(targetfile, pathvar + new_name)  # Copy file to pathvar

def groupavg(data):
    # Returns the average value of a string of data points
    tmp = 0.00
    for value in data:
        tmp += value
    tmp /= len(data)
    return tmp

def fileoutput(
    filename,
    data_1,
    data_2,
    data_3,
    data_4,
    ):
    """Writes quad channel data to file."""
    
    points = len(data_1)
    with open(filename, 'w') as FILE:
        for i in range(0, points):
            FILE.write('%s' % data_1[i] + '\t %s' % data_2[i] + '\t %s'
                       % data_3[i] + '\t %s \n' % data_4[i])

def chopper(
    data_1,
    data_2,
    data_3,
    data_4,
    l_th,
    u_th,
    filecounter,
    output_dir,
    ):
    """Chops the continuous data file into individual scans."""
   
    data_window = []
    # The number of points to not include from the determined 'end' (prevents problems with other channels)
    int_bodge = 20 
    window_length = 50
    l_plat_len = 0
    u_plat_len = 0
    counter = 0
    start = 0
    stop = 0
    
    # Only data above u_th should be discarded as data below l_th is required for the background correction
    # Detect it the tip substrate bias is positive or negative to allow auto-inversion of data
    # The average of the LSx1 channel is a good, but inefficient measure of this
    # TODO: If the warning appears frequently parameterise this in the GUI with a manual override

    dat_avg = np.average(data_1)
    print "LSx1 channel average is: ", dat_avg, "V"
    if (dat_avg < -0.5):
        print "Assuming negative tip-substrate bias"
        #Invert the data
        #FIXME get rid of the now redundant inversion check when reading the split files back in
        for i in range(len(data_1)):
                data_1[i] = data_1[i] * -1.0
                data_2[i] = data_2[i] * -1.0
                data_3[i] = data_3[i] * -1.0
                data_4[i] = data_4[i] * -1.0
    elif (dat_avg > 0.5):
        print "Assuming positive tip-substrate bias"
    else:
        print "ABORTING!!! Low set point current detected; automatic polarity check failed"
        sys.exit("Failed bias check")

    # filecounter = 0
    # Implement minimum

    for value in data_1:
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
            if groupavg(data_window) < u_th and u_plat_len > 0 and stop \
                == 0:
                # Hopefully this is the tip retraction point and from here on the current decays
                # Stop must be zero otherwise we might end up with background, then decay!
                start = counter  # could make this "counter - u_plat_len" to get all the plateaus but not of interest
                # print "Found upper plateau: ", u_plat_len, " points long, stopping at: ", start
                u_plat_len = 0
            if groupavg(data_window) > l_th and l_plat_len > 0 \
                and start != 0:
                # We found the end of the background plateau at the previous counter value
                # start must be > 0 to ensure we have already found the initial current decay
                stop = counter - int_bodge
                # start is at the location of the last found upper plateau
                # print "Found lower plateau: ", l_plat_len, " point long, stopping at: ", stop........
                points = stop - start
                # A bit of a 'hatchet job', basically check CH2 remains saturated for x number of points,
                # akin to checking the gradient of CH1
                # TODO: Bug here: this assumes the channnel is saturated which is only the case for BJ experiments
                sat_level = data_2[start]
                sat_pass = True
                if len(data_2) >= (start+50): #Get rid of obscure error which happened in T54
                    for i in range(0, 50):
                        if sat_level != data_2[start + i]:
                            sat_pass = False
                            print "REJECTED: No stable set-point current detected"

                        # print sat_level, data_2[i]
                # Dislocations due to moving the scanner during a I(s) measurements need to be removed.
                # For example:
                # 4.01572.... 10.3392.... 10.454.... 10.3357
                # 0.0261655.... 0.252062.... 6.06029.... 10.3357
                # But for now we don't worry about that here and deal with it later in the stitcher
                # Now check to see if the high sens channel has reached zero, otherwise it's a #*$% end point / measurement

                lowchan_pass = True
                background_level = groupavg(data_3[stop - 10:stop])
                if background_level > l_th:
                    lowchan_pass = False
                    # The level will depend on the 'characteristic level' of the experiment.
                    # For bare contacts this should of the order of mV. For in-situ work it will depend on the tip leakage current. 
                    print 'REJECTED: HS x10 voltage decayed to:', background_level, '(V). Threshold level:', l_th, '(V)'
   
                # Save the measurements if they pass some tests:
                # Each scan should have a characteristic length. If you wanted to be flash you could plot a histogram of
                # Scan length. Then choose only scans within a couple of SDs of the mean length.
                # In this case the length is controlled manually and will depend on the scan duration in seconds.
                # If the scan duration is changed adjust the acceptable range to the characteristic length
                # This acts as a filter to prevent multiple length scans, or unreasonably short I(s) scans getting through

                if points > controller.min_seg_len.get() and points < controller.max_seg_len.get() and sat_pass \
                    and lowchan_pass:
                    data_slice1 = data_1[start:stop]
                    data_slice2 = data_2[start:stop]  # Matching one from x10 channel
                    data_slice3 = data_3[start:stop]
                    data_slice4 = data_4[start:stop]

                    # Generate the filename
                    # Later on this can calculate a sdev and dynamically exclude

                    filename = os.getcwd()[:-3] + output_dir + "/" + 'slice' + str(filecounter).zfill(4) \
                        + '.txt'
                    print 'Reconstructed I(s) scan containing: ', points, \
                        'data points as: ', filename
                    fileoutput(
                        filename,
                        data_slice1,
                        data_slice2,
                        data_slice3,
                        data_slice4,
                        )
                    filecounter += 1
                l_plat_len = 0
                start = 0
                stop = 0
        counter += 1
    return filecounter


def tea_break_maker():

    """Deals with reading all raw data files from the ADC and splitting them into I(s) scans using chopper"""
    # Attempt to prevent laziness by asking some questions about the experiment before processing the data
    # After these have been answered they'll be plenty of time for a tea break.

    questions = [
        "Experiment name                 :",
        "Sweep time (seconds)            :",
        "Sweep range (nm)                :",
        "Substrate material              :",
        "Substrate SAM?                  :",
        "Tip material                    :",
        "Tip SAM?                        :",
        "Environment                     :",
        "Magnetic field?                 :",
        "Preamp L sense (kOhm)           :",
        "Preamp H sense (MOhm)           :",
        "Other comments                  :",
        ]

    # Have a poke around in the raw directory so we can ask for comments on each file within
    files = os.listdir(os.path.abspath('raw'))
    for filename in files:
        # For each precious measurement we must have a description
        questions.append(filename
                         + " details (bias/setpoint/other)	:")

    #Generate some dialogue to guide/warn the user
    if (len(files) < 1):
        error = showerror('Disaster', 'No data found in ../raw')
        return
    info = showinfo('Scan reconstructor', 'Refer to the terminal')
    answers = []
   
    # Now ask the questions
    if not os.access('nfo.txt', os.F_OK):
        for question in questions:
            var = raw_input(question)
            answers.append(var)
        # Write the questions and answer to file for later viewing,
        # but if the file exists then don't overwrite it
        with open('nfo.txt', 'w') as FILE:
            for i in range(0, len(answers)):
                info = questions[i] + '\t' + answers[i] + '\n'
                FILE.write(info)
    else:
        print "Skipping questions: nfo.txt already exists"

    # Now for the chopping
    # Assumes you've put the raw data in ../raw/
    # Go into the raw data directory
    original_directory = os.getcwd()
    try:
        os.chdir(os.path.abspath('raw'))
        files = os.listdir(os.getcwd())
        #Start reconstruction of scans
        print 'This is a good time to have tea break...'
        process_files(files)
    finally:
        print 'Finished reconstructing I(s) scans'
        os.chdir(original_directory)

def process_files(files):
    def producer(q, files):
        for filename in files:
            thread = InitChopper(filename)
            thread.start()
            print "------------->   THREAD NAME ----------->", thread.getName()
            #Put the thread in the queue, arg True blocks if necessary until a slot is available 
            #FIXME: Waits for thread to finish -- not parallel. Rewrite chopper in C?
            thread.join()
            q.put(thread, True)

    def consumer(q, total_files):
        counter = 0
        while counter < total_files:
            #Remove and return an item from the queue, arg True blocks if necessary until item available
            thread = q.get(True)
            #Wait for thread to terminate
            thread.join()
            counter +=1

    q = Queue(3)
    prod_thread = threading.Thread(target=producer, args=(q, files))
    cons_thread = threading.Thread(target=consumer, args=(q, len(files)))
    prod_thread.start()
    cons_thread.start()
    #Wait for producer to terminate
    prod_thread.join() 
    cons_thread.join()

class InitChopper(threading.Thread):
    def __init__(self, filename):
        self.raw_data_filename = filename
        threading.Thread.__init__(self)

    def run(self):
        #NEED TO USE absolute paths, none of this CHDIR business
        # Reset the file counter for the chopped scans
        filecounter = 0
        data = True
        #Open the raw data file which should be gzipped 
        with gzip.open(self.raw_data_filename) as gz_data:
            # Create an output folder for the individual scans
            dirname = self.raw_data_filename[0:-7]
            output = 'chopped' + dirname[6:]
            print 'Reconstructing I(s) scans into the directory', output, '...'
            # Go out of the raw data folder........
            raw_directory = os.getcwd()
            try:
                #TODO MAKE THREAD SAFE -- PUT A LOCK HERE
                os.chdir(os.pardir)
                # Make the output folder for the I(s) scans if it doesn't exist already and cd into it
                if os.path.exists(output) != True:
                    os.mkdir(output)
                os.chdir(raw_directory)

                #os.chdir(os.path.abspath(output))
                # Reconstruct the I(s) scans only if the the folder is empty       
                #if os.listdir(os.getcwd()):
                #    print "Skipping scan reconstruction: target folder:", "../" + output, "is not empty"
                #    return
                #TODO REMOVE LOCK  -- change below path to abs 
                while data == True:
                    #Create empty current lists
                    i_list_ls_x1 = []
                    i_list_ls_x10 = []
                    i_list_hs_x1 = []
                    i_list_hs_x10 = []
                    for i in range(controller.chunksize.get()):
                        line = gz_data.readline()
                        if not line:
                            print "End of raw data: last chunk is", i," lines long."
                            #Don't do any more while loops
                            data = False
                            #Leave the for loop with i lines of data        
                            break
                        line = line.split()
                        #Populate the current lists
                        i_list_ls_x1.append(float(line[0]))
                        i_list_ls_x10.append(float(line[1]))
                        i_list_hs_x1.append(float(line[2])) 
                        i_list_hs_x10.append(float(line[3]))
                    #Now reconstruct the scans, but only if there is at least one line of data in the chunk
                    if line:
                        filecounter = chopper(
                            i_list_ls_x1,
                            i_list_ls_x10,
                            i_list_hs_x1,
                            i_list_hs_x10,
                            controller.scan_l_th.get(),
                            controller.scan_u_th.get(),
                            filecounter,
                            output
                            )
            finally:
                os.chdir(raw_directory)

class egraph:

    # Plots embedded matplotlib graphs in main window
    def __init__(self, myParent):
        #Deal with different screen resolutions
        if root.winfo_screenheight() <= 768:
            #Low res so downsize figures, ideal for XGA, WXGA but no smaller (unlikely!)
            self.f = Figure(figsize=(10, 7), dpi=80)
        else:
            #This is fine for SXGA and higher
            self.f = Figure(figsize=(14, 9), dpi=80)
        #Add the subplots
        self.ax = self.f.add_subplot(131)
        self.ax2 = self.f.add_subplot(132)
        self.ax3 = self.f.add_subplot(133)
        self.canvas = FigureCanvasTkAgg(self.f, master=myParent)
        #Add the toolbar
        toolbar = NavigationToolbar2TkAgg( self.canvas, root )
        toolbar.pack(side=BOTTOM, fill=BOTH, expand=1)
        toolbar.update()
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=BOTTOM)

    def confScanPlot(self):
        """Plots histogram with overlaid Gaussian"""
        
        #Reconfigure the axes if they're setup for scan plotting
        try:
            self.f.delaxes(self.ax)
            self.f.delaxes(self.ax2)
            self.f.delaxes(self.ax3)
            self.g_ax = self.f.add_subplot(111)
        except ValueError:
            #Clear the exisiting plot
            self.g_ax.cla() 
                    
        # Fetch KDE parameters from control widgets
        kde_bandwidth = controller.kde_bandwidth.get()
        kde_start = controller.kde_bandwidth  # The lowest the KDE start can go
        kde_stop = controller.kde_stop.get()
        kde_points = controller.kde_points.get()

        #Fetch the plot limits
        self.g_xmin = controller.gF_xmin.get()
        self.g_xmax = controller.gF_xmax.get()
        self.g_ymin = controller.gF_ymin.get()
        self.g_ymax = controller.gF_ymax.get()
        #Set plot limits
        self.g_ax.set_xlim([self.g_xmin, self.g_xmax])
        self.g_ax.set_ylim([self.g_ymin, self.g_ymax])

        # Plot KDE function and read all current data to a single list
        # Initialse current list and append currents from individual files
        i_list = []

        #FIXME CLean this up when polished
        DEBUG = False      
        if DEBUG:
            for i in range(1000):
                i_list.append(random.gauss(1,0.5))    
        else:
            for reading in data.i_dat_all_combined:
                for value in reading:
                    if value > 0:
                        i_list.append(np.log10(value))
       
        print 'Lists appended, found:', len(i_list), 'data points.'

        if (len(i_list) < 1):
            self.error = showerror('Error', 'No data in memory')
            return
        
        #Compute the KDE
        self.x = np.linspace(min(i_list),max(i_list),kde_points)
        self.z = statistics.pdf(i_list, self.x, h=kde_bandwidth, kernel='E')
        
        #Plot KDE and histogram
        self.g_ax.set_xlabel(r'$log_{10}(G/G_{0}$)')
        self.g_ax.set_ylabel(r'$Density$')

        #Set the title
        try:
            self.g_title = data.filename[0:string.rfind(data.filename, '/')]
        except AttributeError:
            self.g_title = "Folder undefined"
        self.g_ax.set_title(self.g_title)

        #Plot the graph
        self.g_ax.plot(self.x, self.z,'b', label='KDE', linewidth=3)
        self.g_ax.text(self.g_xmin, self.g_ymax*0.95, ' ' + str(len(data.i_dat_all_combined)) + ' scans', {'color' : 'k', 'fontsize' : 18})
        self.g_ax.legend()
        self.g_ax.grid()
        self.canvas.show()

    def plotGaussian(self):
        
        #Fetch the Gaussian parameters
        mu1 = controller.gF_mu1.get()
        sigma1 = controller.gF_sigma1.get() 
        scale1 = controller.gF_scale1.get()
        ydisp1 = controller.gF_ydisp1.get()
        mu2 = controller.gF_mu2.get()
        sigma2 = controller.gF_sigma2.get() 
        scale2 = controller.gF_scale2.get()
        ydisp2 = controller.gF_ydisp2.get()
        mu3 = controller.gF_mu3.get()
        sigma3 = controller.gF_sigma3.get() 
        scale3 = controller.gF_scale3.get()
        ydisp3 = controller.gF_ydisp3.get()
        mu4 = controller.gF_mu4.get()
        sigma4 = controller.gF_sigma4.get() 
        scale4 = controller.gF_scale4.get()
        ydisp4 = controller.gF_ydisp4.get()
        
        #Clear axes and replot everything
        self.g_ax.cla() #TODO: Rather than clearing and replotting can you remove last previous plot??
        self.g_ax.grid()
        self.g_ax.plot(self.x, self.z,'b', label='KDE', linewidth=2)
        self.fit1 = mlab.normpdf( self.x, mu1, sigma1)*scale1 + ydisp1    
        self.g_ax.plot(self.x, self.fit1, 'r-', label='Fit 1', linewidth=2)
        self.fit2 = mlab.normpdf( self.x, mu2, sigma2)*scale2 + ydisp2   
        self.g_ax.plot(self.x, self.fit2, 'g-', label='Fit 2', linewidth=2)
        self.fit3 = mlab.normpdf( self.x, mu3, sigma3)*scale3 + ydisp3   
        self.g_ax.plot(self.x, self.fit3, 'c-', label='Fit 3', linewidth=2)
        self.fit4 = mlab.normpdf( self.x, mu4, sigma4)*scale4 + ydisp4   
        self.g_ax.plot(self.x, self.fit4, 'm-', label='Fit 4', linewidth=2)

        combo = self.fit1 + self.fit2 + self.fit3 + self.fit4
        self.g_ax.plot(self.x, combo, 'k-', label='Combined Fit', linewidth=3)
        self.g_ax.legend()

        #Set some labels
        self.g_ax.text(self.g_xmin, self.g_ymax*0.95, ' ' + str(len(data.i_dat_all_combined)) + ' scans', {'color' : 'k', 'fontsize' : 16})
        self.g_ax.text(self.g_xmin, self.g_ymax*0.90, r' $\mu_{1} = $' + str(mu1), {'color' : 'k', 'fontsize' : 16})
        self.g_ax.text(self.g_xmin, self.g_ymax*0.85, r' $\sigma_{1} = $' + str(sigma1), {'color' : 'k', 'fontsize' : 16})
        self.g_ax.text(self.g_xmin, self.g_ymax*0.80, r' $\mu_{2} = $' + str(mu2), {'color' : 'k', 'fontsize' : 16})
        self.g_ax.text(self.g_xmin, self.g_ymax*0.75, r' $\sigma_{2} = $' + str(sigma2), {'color' : 'k', 'fontsize' : 16})
        self.g_ax.text(self.g_xmin, self.g_ymax*0.70, r' $\mu_{3} = $' + str(mu3), {'color' : 'k', 'fontsize' : 16})
        self.g_ax.text(self.g_xmin, self.g_ymax*0.65, r' $\sigma_{3} = $' + str(sigma3), {'color' : 'k', 'fontsize' : 16})
        self.g_ax.text(self.g_xmin, self.g_ymax*0.60, r' $\mu_{4} = $' + str(mu4), {'color' : 'k', 'fontsize' : 16})
        self.g_ax.text(self.g_xmin, self.g_ymax*0.55, r' $\sigma_{4} = $' + str(sigma4), {'color' : 'k', 'fontsize' : 16})
        
        #TODO: A lot of replication between here and confScanplot 
        #Plot KDE and histogram
        self.g_ax.set_xlabel(r'$log_{10}(G/G_{0}$)')
        self.g_ax.set_ylabel(r'$Density$')
        self.g_ax.set_title(self.g_title)
        self.g_ax.set_xlim([self.g_xmin, self.g_xmax])
        self.g_ax.set_ylim([self.g_ymin, self.g_ymax])
        self.canvas.show()
        
    def updater(
        self,
        s_dat_ls,
        i_dat_ls,
        i10_dat_ls,
        i_dat_hs,
        i10_dat_hs,
        i_dat_combined,
        plat_data,
        title,
        ):

        #Configure the axes
        try:
            #Delete the gaussian plot if it exists
            self.f.delaxes(self.g_ax)
            #Then setup the scan axes
            self.ax = self.f.add_subplot(131)
            self.ax2 = self.f.add_subplot(132)
            self.ax3 = self.f.add_subplot(133)
        except ValueError:
            pass
        except AttributeError:
            pass

        #The limit of the x,y axis plotted in various figures
        xlim = float(controller.xfac.get())
        ylim = float(controller.yfac.get())
        #Constant offset to add to each current value to shift the scan from negative region
        offset = float(controller.offset.get())
        
        self.ax.cla()  # Clear current axes
        # PLOT 1
        self.ax.set_yscale('log')
        self.ax.plot(s_dat_ls, i_dat_ls, 'k.', label='LS x1')
        self.ax.plot(s_dat_ls, i10_dat_ls, 'b.', label='LS x10')
        self.ax.plot(s_dat_ls, i_dat_hs, 'g.', label='HS x1')
        self.ax.plot(s_dat_ls, i10_dat_hs, 'r.', label='HS x10')
        self.ax.grid(True)
        self.ax.set_xlim([0, xlim])
        self.ax.set_xlabel(r'$Distance (nm)$')
        self.ax.set_ylabel(r'$G/G_{0}$')
        self.ax.legend()

        # PLOT 2
        self.ax2.cla()  # Clear current axes
        self.ax2.set_yscale('log')
        #Add the current offset (only for cosmetic purposes to stop negative data getting ignored)
        i_dat_shifted = []
        for value in i_dat_combined:
            shifted_val = value + offset
            i_dat_shifted.append(shifted_val)
        self.ax2.plot(s_dat_ls, i_dat_shifted, 'k.',
                      label='Combined data')#, shift (nA):' + str(offset))
        if controller.autocheck_linfit.get():
            self.ax2.plot(data.s_dat_filtered, data.polyfit_rescaled,
                          'r', label='MSE:' + str(data.err))
        self.ax2.grid(True)
        self.ax2.set_xlim([0, xlim])
        self.ax.set_xlabel(r'$Distance (nm)$')
        self.ax2.legend()

        # PLOT 3
        self.ax3.cla()  # Clear current axes
        self.ax3.plot(s_dat_ls, i_dat_ls, 'k.-', label='LS x1')
        self.ax3.plot(s_dat_ls, plat_data, 'b.-', label='Plateau fitting')
        self.ax3.grid(True)
        self.ax3.set_ylim([0, ylim])
        self.ax3.set_xlim([0, xlim])
        self.ax.set_xlabel(r'$Distance (nm)$')
        self.ax3.legend()

        self.canvas.show()

#Callback for dynamic updating of plot
def callback(var):
    #Catch invalid variables
    try: 
        egraph.plotGaussian()
        return True
    except ValueError:
        return False

# Main GUI
class controller:

    # This generates the main GUI and contains all variables the use can modify
    def __init__(self, myParent):

        #Use an interrupt to terminate the script on Ctrl+C
        signal.signal(signal.SIGINT, self.signal_handler)
    
        #Save the root path of the gui
        root_path = os.getcwd() 

        # Constants for controlling layout
        button_width = 15
        button_padx = '2m'
        button_pady = '1m'
        buttons_frame_padx = '3m'
        buttons_frame_pady = '2m'
        buttons_frame_ipadx = '3m'
        buttons_frame_ipady = '1m'

        # Initiate Variables
        self.stavar = IntVar()
        self.finvar = IntVar()
        self.bcorfac = DoubleVar()

        # Data processing variables
        self.bkgnd_tol = IntVar()
        self.xfac = DoubleVar()
        self.yfac = DoubleVar()
        self.offset = DoubleVar()
        self.currentScaleFactor = DoubleVar()
        self.plot_dat = IntVar()
        self.check_dat = IntVar()
        self.check_dat2 = IntVar()
        self.autocheck_linfit = IntVar()
        self.autocheck_avgcur = IntVar()
        self.autocheck_pltfit = IntVar()
        self.auto_read = IntVar()
        self.clear_global_data = IntVar()
        self.bias = DoubleVar()

        # Current Density plots
        self.xmin_cd = DoubleVar()
        self.xmax_cd = DoubleVar()
        self.ymin_cd = DoubleVar()
        self.ymax_cd = DoubleVar()
        self.bins_cd = IntVar()
        self.kde_bandwidth = DoubleVar()
        self.kde_stop = DoubleVar()
        self.kde_points = IntVar()
        self.kdePeakLabel = IntVar()
        self.kdePeakStart = DoubleVar()
        self.kdePeakStop = DoubleVar()
        

        # Resistor division
        self.lowres = DoubleVar()
        self.highres = DoubleVar()
        self.gainfactor = DoubleVar()

        # Stitch volages
        self.th_1 = DoubleVar()
        self.th_2 = DoubleVar()
        self.th_3 = DoubleVar()

         # STM/ADC variables sampsec, srange, sduration
        self.sampsec = IntVar()
        self.srange = DoubleVar()
        self.sduration = DoubleVar()

        # Contour plot parameters
        self.xmin_contour = DoubleVar()
        self.xmax_contour = DoubleVar()
        self.ymin_contour = DoubleVar()
        self.ymax_contour = DoubleVar()
        self.xbin_contour = IntVar()
        self.ybin_contour = IntVar()

        #Gaussian fit parameters
        self.gF_mu1 = DoubleVar()
        self.gF_sigma1 = DoubleVar()
        self.gF_scale1 = DoubleVar()
        self.gF_ydisp1 = DoubleVar()
        self.gF_mu2 = DoubleVar()
        self.gF_sigma2 = DoubleVar() 
        self.gF_scale2 = DoubleVar()
        self.gF_ydisp2 = DoubleVar()
        self.gF_mu3 = DoubleVar()
        self.gF_sigma3 = DoubleVar() 
        self.gF_scale3 = DoubleVar()
        self.gF_ydisp3 = DoubleVar()
        self.gF_mu4 = DoubleVar()
        self.gF_sigma4 = DoubleVar() 
        self.gF_scale4 = DoubleVar()
        self.gF_ydisp4 = DoubleVar()
        self.gF_xmin = DoubleVar()
        self.gF_xmax = DoubleVar()
        self.gF_ymin = DoubleVar()
        self.gF_ymax = DoubleVar()

        #2D Correlation plot parameters
        self.corrcurrentmin = DoubleVar()
        self.corrcurrentmax = DoubleVar()
        self.corrbins = IntVar()
        self.corrlogscale = IntVar()

        # Plateau fitting parameters
        self.plat_max_grad = DoubleVar()
        self.background_tol = IntVar()
        self.max_plat_cur = DoubleVar()
        self.max_points_plat = IntVar()
        self.fractional_plat_tol = DoubleVar()
        self.min_points_plat = IntVar()

        #Chopper parameters
        self.scan_u_th = DoubleVar()
        self.scan_l_th = DoubleVar()
        self.chunksize = IntVar()
        self.max_seg_len = IntVar()
        self.min_seg_len = IntVar()

        # Data filtering
        self.datfillogi = IntVar()

        #Check to see if some select variables have been pickled:

        try: 
            open("settings")
        except:
            print "Settings file not found, using defaults."
            # Set default values for variables
        
            # Resistor defaults
            self.lowres.set(10000)  # Calibrate from resistor measurements, assumes output voltage is correct from STM
            self.highres.set(10000000)  # Calibrate from resistor measurements, assumes output voltage is correct from STM
            self.gainfactor.set(10.0) # Unless the hardware has been modified this shouldn't need to be changed 
            self.th_1.set(1.0)  # Ch2 is simply Ch1 x 10 so at this stage Ch2 is at 90% maximum and can take over
            self.th_2.set(0.08)  # Around here the limiter turns off, plus a bit more to stop dodgy overlap
            self.th_3.set(1.0)  # Ch3 is on and can take over
            
            #Data input
            self.stavar.set(1)
            self.finvar.set(5)
            self.bcorfac.set(0.20)
            self.xfac.set(4.00)
            self.yfac.set(25000)
            self.currentScaleFactor.set(1e9)
            self.offset.set(0.00)
            self.plot_dat.set(0)
            self.check_dat.set(0)
            self.check_dat2.set(0)
            self.auto_read.set(0)
            self.clear_global_data.set(1)
            self.bias.set(0.1)

            # Data filtering defaults:
            self.datfillogi.set(-1000)

            # Current density estimate plot:
            self.xmin_cd.set(-6.0)
            self.xmax_cd.set(1.0)
            self.ymin_cd.set(0.0)
            self.ymax_cd.set(0.4) # Current density estimate plot:
            self.bins_cd.set(1000)
            self.kde_bandwidth.set(0.10)
            self.kde_stop.set(4)
            self.kde_points.set(1000)
            self.kdePeakLabel.set(1)
            self.kdePeakStart.set(0.0)
            self.kdePeakStop.set(1.5)         

            # STM/ADC parameters used for reconstructing distance defaults
            self.sampsec.set(10000)
            self.srange.set(6)
            self.sduration.set(0.3)

            # Contour plot parameter defaults
            self.xmin_contour.set(0.00)
            self.xmax_contour.set(3.0)
            self.ymin_contour.set(-6.0)
            self.ymax_contour.set(1.0)
            self.xbin_contour.set(250)
            self.ybin_contour.set(250)

            #Gaussian fit parameters - defaults to TMB/Au/100mV
            self.gF_mu1.set(0.135)
            self.gF_sigma1.set(0.22)
            self.gF_scale1.set(0.023)
            self.gF_ydisp1.set(0.00)
            self.gF_mu2.set(0.67)
            self.gF_sigma2.set(0.17)
            self.gF_scale2.set(0.0285)
            self.gF_ydisp2.set(0.00)
            self.gF_mu3.set(0.96)
            self.gF_sigma3.set(0.17)
            self.gF_scale3.set(0.0175)
            self.gF_ydisp3.set(0.00)
            self.gF_mu4.set(1.14)
            self.gF_sigma4.set(0.17)
            self.gF_scale4.set(0.0145)
            self.gF_ydisp4.set(0.00)
            self.gF_xmin.set(-4.0)
            self.gF_xmax.set(2.0)
            self.gF_ymin.set(0.0)
            self.gF_ymax.set(0.25)

            # 2D Correlation plot parameter defaults
            self.corrcurrentmin.set(-5.0)
            self.corrcurrentmax.set(1.0)
            self.corrbins.set(300)
            self.corrlogscale.set(1)

            # Plateau fitting parameters
            self.plat_max_grad.set(50000.0)
            self.background_tol.set(5)
            self.max_plat_cur.set(10000.0)
            self.max_points_plat.set(100)
            self.max_points_plat.set(100)
            self.fractional_plat_tol.set(0.20)
            self.min_points_plat.set(30)

            #Chopper parameters
            self.scan_u_th.set(4.0)
            self.scan_l_th.set(0.009)
            self.chunksize.set(1000000)
            self.max_seg_len.set(3200)
            self.min_seg_len.set(2200)

        else:
            print "Loading custom settings from file 'settings'."
            data = self.restoreData()
            try:
                #Now set all the variables. See the defaults section above for what they are
                #TODO: This could be done in a 'cleaner way' by iterating over the variable dicionary returned by saveData()
                self.lowres.set(data[ 'lowres'])  
                self.highres.set(data[ 'highres'])  
                self.gainfactor.set(data[ 'gainfactor'])
                self.th_1.set(data[ 'th_1'])  
                self.th_2.set(data[ 'th_2']) 
                self.th_3.set(data[ 'th_3']) 
                self.stavar.set(data[ 'stavar'])  
                self.finvar.set(data[ 'finvar'])  
                self.bcorfac.set(data[ 'bcorfac'])  
                self.xfac.set(data[ 'xfac'])  
                self.yfac.set(data[ 'yfac'])  
                self.offset.set(data[ 'offset'])  
                self.currentScaleFactor.set(data[ 'currentScaleFactor'])
                self.plot_dat.set(data[ 'plot_dat'])  
                self.check_dat.set(data[ 'check_dat'])  
                self.check_dat2.set(data[ 'check_dat2'])  
                self.auto_read.set(data[ 'auto_read'])  
                self.clear_global_data.set(data[ 'clear_global_data']) 
                self.bias.set(data[ 'bias']) 
                self.datfillogi.set(data[ 'datfillogi'])  
                self.xmin_cd.set(data[ 'xmin_cd'])  
                self.xmax_cd.set(data[ 'xmax_cd'])  
                self.ymin_cd.set(data[ 'ymin_cd'])  
                self.ymax_cd.set(data[ 'ymax_cd'])  
                self.bins_cd.set(data[ 'bins_cd'])  
                self.kde_bandwidth.set(data[ 'kde_bandwidth'])  
                self.kde_stop.set(data[ 'kde_stop'])  
                self.kde_points.set(data[ 'kde_points'])  
                self.kdePeakLabel.set(data[ 'kdePeakLabel']) 
                self.kdePeakStart.set(data[ 'kdePeakStart']) 
                self.kdePeakStop.set(data[ 'kdePeakStop'])  
                self.sampsec.set(data[ 'sampsec'])  
                self.srange.set(data[ 'srange'])  
                self.sduration.set(data[ 'sduration'])  
                self.xmin_contour.set(data[ 'xmin_contour'])  
                self.xmax_contour.set(data[ 'xmax_contour'])  
                self.ymin_contour.set(data[ 'ymin_contour'])  
                self.ymax_contour.set(data[ 'ymax_contour'])  
                self.xbin_contour.set(data[ 'xbin_contour'])  
                self.ybin_contour.set(data[ 'ybin_contour'])  
                self.gF_mu1.set(data['gF_mu1'])
                self.gF_sigma1.set(data['gF_sigma1'])
                self.gF_scale1.set(data['gF_scale1'])
                self.gF_ydisp1.set(data['gF_ydisp1'])
                self.gF_mu2.set(data['gF_mu2'])
                self.gF_sigma2.set(data['gF_sigma2'])
                self.gF_scale2.set(data['gF_scale2'])
                self.gF_ydisp2.set(data['gF_ydisp2'])
                self.gF_mu3.set(data['gF_mu3'])
                self.gF_sigma3.set(data['gF_sigma3'])
                self.gF_scale3.set(data['gF_scale3'])
                self.gF_ydisp3.set(data['gF_ydisp3'])
                self.gF_mu4.set(data['gF_mu4'])
                self.gF_sigma4.set(data['gF_sigma4'])
                self.gF_scale4.set(data['gF_scale4'])
                self.gF_ydisp4.set(data['gF_ydisp4'])
                self.gF_xmin.set(data['gF_xmin'])
                self.gF_xmax.set(data['gF_xmax'])
                self.gF_ymin.set(data['gF_ymin'])
                self.gF_ymax.set(data['gF_ymax'])
                self.corrcurrentmin.set(data[ 'corrcurrentmin'])
                self.corrcurrentmax.set(data[ 'corrcurrentmax'])
                self.corrbins.set(data[ 'corrbins'])
                self.corrlogscale.set(data[ 'corrlogscale'])
                self.plat_max_grad.set(data[ 'plat_max_grad'])  
                self.background_tol.set(data[ 'background_tol'])  
                self.max_plat_cur.set(data[ 'max_plat_cur'])  
                self.max_points_plat.set(data[ 'max_points_plat'])  
                self.fractional_plat_tol.set(data[ 'fractional_plat_tol'])  
                self.min_points_plat.set(data[ 'min_points_plat'])  
                self.scan_u_th.set(data[ 'scan_u_th'])  
                self.scan_l_th.set(data[ 'scan_l_th'])  
                self.chunksize.set(data[ 'chunksize'])  
                self.max_seg_len.set(data[ 'max_seg_len'])  
                self.min_seg_len.set(data[ 'min_seg_len'])  
            except KeyError:
                print "ERROR: Settings file incomplete."
                print "You are probably using a new version of the script with an"
                print "old settings file. Delete and recreate the settings file."

        #Trace the relevant variables for dynamic updating of embedded plot
        #name is name of tk var, index if var i an array, otherwise empty string, 
        self.gF_mu1.trace("w", lambda name, index, mode, gF_mu1=self.gF_mu1: callback(gF_mu1))
        self.gF_mu2.trace("w", lambda name, index, mode, gF_mu2=self.gF_mu2: callback(gF_mu2))
        self.gF_mu3.trace("w", lambda name, index, mode, gF_mu3=self.gF_mu3: callback(gF_mu3))
        self.gF_mu4.trace("w", lambda name, index, mode, gF_mu4=self.gF_mu4: callback(gF_mu4))
        self.gF_sigma1.trace("w", lambda name, index, mode, gF_sigma1=self.gF_sigma1: callback(gF_sigma1))
        self.gF_sigma2.trace("w", lambda name, index, mode, gF_sigma2=self.gF_sigma2: callback(gF_sigma2))
        self.gF_sigma3.trace("w", lambda name, index, mode, gF_sigma3=self.gF_sigma3: callback(gF_sigma3))
        self.gF_sigma4.trace("w", lambda name, index, mode, gF_sigma4=self.gF_sigma4: callback(gF_sigma4))
        self.gF_scale1.trace("w", lambda name, index, mode, gF_scale1=self.gF_scale1: callback(gF_scale1))
        self.gF_scale2.trace("w", lambda name, index, mode, gF_scale2=self.gF_scale2: callback(gF_scale2))
        self.gF_scale3.trace("w", lambda name, index, mode, gF_scale3=self.gF_scale3: callback(gF_scale3))
        self.gF_scale4.trace("w", lambda name, index, mode, gF_scale4=self.gF_scale4: callback(gF_scale4))
        self.gF_ydisp1.trace("w", lambda name, index, mode, gF_ydisp1=self.gF_ydisp1: callback(gF_ydisp1))
        self.gF_ydisp2.trace("w", lambda name, index, mode, gF_ydisp2=self.gF_ydisp2: callback(gF_ydisp2))
        self.gF_ydisp3.trace("w", lambda name, index, mode, gF_ydisp3=self.gF_ydisp3: callback(gF_ydisp3))
        self.gF_ydisp4.trace("w", lambda name, index, mode, gF_ydisp4=self.gF_ydisp4: callback(gF_ydisp4))

        # Menu bar at the top of the main window
        self.mBar = Frame(myParent, relief=RAISED, borderwidth=2)
        self.mBar.pack(fill=X)

        # File menu
        self.FileMenu = Menubutton(self.mBar, text='File', underline=0)
        self.FileMenu.pack(side=LEFT, padx='2m')
        self.FileMenu.menu = Menu(self.FileMenu)
        self.FileMenu.menu.add_command(label='Open...', underline=0,
                command=self.filename_browse)
        self.FileMenu.menu.add('separator')
        self.FileMenu.menu.add_command(label='Quit', underline=0,
                command=self.quit)
        self.FileMenu['menu'] = self.FileMenu.menu

        # Plot menu -- everything related to plotting graphs
        self.PlotMenu = Menubutton(self.mBar, text='Plot', underline=0)
        self.PlotMenu.pack(side=LEFT, padx='2m')
        self.PlotMenu.menu = Menu(self.PlotMenu)
        self.PlotMenu.menu.add_command(label='Logarithmic current density plot'
                , underline=0, command=self.kdePlot)
        self.PlotMenu.menu.add_command(label='Linear current histogram'
                , underline=1, command=self.linearPlot)
        self.PlotMenu.menu.add_command(label='2D current-distance histogram'
                , underline=1, command=contourPlot)
        self.PlotMenu.menu.add_command(label='2D correlation histogram'
                , underline=1, command=correlationHist)
        self.PlotMenu.menu.add_command(label='Embedded logarithmic current density plot (testing)'
                , underline=1, command=egraph.confScanPlot)
        self.PlotMenu['menu'] = self.PlotMenu.menu

        # Scan analysis menu -- everything related to reading in and filtering individual I(s) scans
        self.ScanAnalysis = Menubutton(self.mBar, text='Scan analysis',
                underline=0)
        self.ScanAnalysis.pack(side=LEFT, padx='2m')
        self.ScanAnalysis.menu = Menu(self.ScanAnalysis)
        self.ScanAnalysis.menu.add_command(label='Read scans to memory'
                , underline=0, background='grey', activebackground='red'
                , command=self.readFiles)
        self.ScanAnalysis.menu.add_command(label='Automatically generate current density plots'
                , underline=0, background='grey', activebackground='red'
                , command=self.autoPlot)
        self.ScanAnalysis.menu.add_command(label='Sync scans in memory to plateau'
                , underline=0, background='grey', activebackground='red'
                , command=self.plat_syncer)
        #self.ScanAnalysis.menu.add('separator')
        self.ScanAnalysis.menu.add_command(label='Export current list to file'
                , underline=0, background='grey', activebackground='green'
                , command=export_current_data)
        self.ScanAnalysis.menu.add_command(label='Export scans in memory to file'
                , underline=0, background='grey', activebackground='green'
                , command=export_scan_data)
        self.ScanAnalysis.menu.add('separator')
        self.ScanAnalysis.menu.add_checkbutton(label='Clear global data on read'
                , underline=0, variable=self.clear_global_data)
        self.ScanAnalysis.menu.add_checkbutton(label='Read all scans in folder'
                , underline=0, variable=self.auto_read)
        self.ScanAnalysis.menu.add_checkbutton(label='Plot scans on-the-fly'
                , underline=0, variable=self.plot_dat)
        self.ScanAnalysis.menu.add('separator')
        self.ScanAnalysis.menu.add_checkbutton(label='Filter scans using manual selection'
                , underline=0, variable=self.check_dat)
        self.ScanAnalysis.menu.add_checkbutton(label='Filter scans using linear regression (alpha)'
                , underline=0, variable=self.autocheck_linfit)
        self.ScanAnalysis.menu.add_checkbutton(label='Filter scans using average current'
                , underline=0, variable=self.autocheck_avgcur)
        self.ScanAnalysis.menu.add_checkbutton(label='Filter scans using plateau fitting'
                , underline=0, variable=self.autocheck_pltfit)
        self.ScanAnalysis['menu'] = self.ScanAnalysis.menu

        # Data menu -- everything related to data input
        self.ImportData = Menubutton(self.mBar, text='Process raw data',
                underline=0)
        self.ImportData.pack(side=LEFT, padx='2m')
        self.ImportData.menu = Menu(self.ImportData)
        self.ImportData.menu.add_command(label='Reconstruct I(s) scans', underline=0,
                background='grey', activebackground='red',
                command=self.importdata)
        self.ImportData.menu.add('separator')
        self.ImportData.menu.add_command(label='Info', underline=0,
                command=self.recon_scan_info)
        self.ImportData.menu.add('separator')

        self.ImportData['menu'] = self.ImportData.menu

        # Settings menu -- lots of variables to play with
        self.SettingsMenu = Menubutton(self.mBar, text='Settings',
                underline=0)
        self.SettingsMenu.pack(side=LEFT, padx='2m')
        self.SettingsMenu.menu = Menu(self.SettingsMenu)
        self.SettingsMenu.menu.add_command(label='Data input',
                underline=0, command=self.data_params)
        self.SettingsMenu.menu.add_command(label='Data filtering',
                underline=0, command=self.data_filter)
        self.SettingsMenu.menu.add_command(label='STM / ADC configuration'
                , underline=0, command=self.adc_params)
        self.SettingsMenu.menu.add_command(label='2D histogram'
                , underline=0, command=self.contour_params)
        self.SettingsMenu.menu.add_command(label='2D correlation histogram'
                , underline=0, command=self.correlation_params)
        self.SettingsMenu.menu.add_command(label='Current density histogram',
                underline=0, command=self.kde_params)
        self.SettingsMenu.menu.add_command(label='Gaussian fitting',
                underline=0, command=self.gaussian_params)
        self.SettingsMenu.menu.add_command(label='Quad channel module calibration'
                , underline=0, command=self.quadchannel_params)
        self.SettingsMenu.menu.add_command(label='Plateau fitting',
                underline=0, command=self.plateau_fitting_params)
        self.SettingsMenu.menu.add_command(label='Scan reconstruction',
                underline=0, command=self.chopper_params)
        self.SettingsMenu.menu.add('separator')
        self.SettingsMenu.menu.add_command(label='Save settings'
                , underline=0, background='grey', activebackground='red'
                , command=self.saveData)
        self.SettingsMenu['menu'] = self.SettingsMenu.menu
    
    def signal_handler(self, signal, frame):
        print '\nKeyboard interrupt; terminating immediately!'
        sys.exit(0)
    
    def saveData(self):

        # Pickle a dictionary containing variables to be saved
        # See the default settings above for what they do
        data = {'lowres' : self.lowres.get(),
                'highres' : self.highres.get(),
                'gainfactor' : self.gainfactor.get(),
                'th_1' : self.th_1.get(),
                'th_2' : self.th_2.get(),
                'th_3' : self.th_3.get(),
                'stavar' : self.stavar.get(),
                'finvar' : self.finvar.get(),
                'bcorfac' : self.bcorfac.get(),
                'xfac' : self.xfac.get(),
                'yfac' : self.yfac.get(),
                'offset' : self.offset.get(),
                'currentScaleFactor' : self.currentScaleFactor.get(),
                'plot_dat' : self.plot_dat.get(),
                'check_dat' : self.check_dat.get(),
                'check_dat2' : self.check_dat2.get(),
                'auto_read' : self.auto_read.get(),
                'clear_global_data' : self.clear_global_data.get(),
                'bias' : self.bias.get(),
                'datfillogi' : self.datfillogi.get(),
                'xmin_cd' : self.xmin_cd.get(),
                'xmax_cd' : self.xmax_cd.get(),
                'ymin_cd' : self.ymin_cd.get(),
                'ymax_cd' : self.xmax_cd.get(),
                'bins_cd' : self.bins_cd.get(),
                'kde_bandwidth' : self.kde_bandwidth.get(),
                'kde_stop' : self.kde_stop.get(),
                'kde_points' : self.kde_points.get(),
                'kdePeakLabel' : self.kdePeakLabel.get(),
                'kdePeakStart' : self.kdePeakStart.get(),
                'kdePeakStop' : self.kdePeakStop.get(),
                'sampsec' : self.sampsec.get(),
                'srange' : self.srange.get(),
                'sduration' : self.sduration.get(),
                'xmin_contour' : self.xmin_contour.get(),
                'xmax_contour' : self.xmax_contour.get(),
                'ymin_contour' : self.ymin_contour.get(),
                'ymax_contour' : self.ymax_contour.get(),
                'xbin_contour' : self.xbin_contour.get(),
                'ybin_contour' : self.ybin_contour.get(),
                'gF_mu1' : self.gF_mu1.get(),
                'gF_sigma1': self.gF_sigma1.get(),
                'gF_scale1' : self.gF_scale1.get(),
                'gF_ydisp1' : self.gF_ydisp1.get(),
                'gF_mu2' : self.gF_mu2.get(),
                'gF_sigma2' : self.gF_sigma2.get(),
                'gF_scale2' : self.gF_scale2.get(),
                'gF_ydisp2' : self.gF_ydisp2.get(),
                'gF_mu3' : self.gF_mu3.get(),
                'gF_sigma3' : self.gF_sigma3.get(),
                'gF_scale3' : self.gF_scale3.get(),
                'gF_ydisp3' : self.gF_ydisp3.get(),
                'gF_mu4' : self.gF_mu4.get(),
                'gF_sigma4' : self.gF_sigma4.get(),
                'gF_scale4' : self.gF_scale4.get(),
                'gF_ydisp4' : self.gF_ydisp4.get(),
                'gF_xmin' : self.gF_xmin.get(),
                'gF_xmax' : self.gF_xmax.get(),
                'gF_ymin' : self.gF_ymin.get(),
                'gF_ymax' : self.gF_ymax.get(),
                'corrcurrentmin' : self.corrcurrentmin.get(),
                'corrcurrentmax' : self.corrcurrentmax.get(),
                'corrbins' : self.corrbins.get(),
                'corrlogscale' : self.corrlogscale.get(),
                'plat_max_grad' : self.plat_max_grad.get(),
                'background_tol' : self.background_tol.get(),
                'max_plat_cur' : self.max_plat_cur.get(),
                'max_points_plat' : self.max_points_plat.get(),
                'fractional_plat_tol' : self.fractional_plat_tol.get(),
                'min_points_plat' : self.min_points_plat.get(),
                'scan_u_th' : self.scan_u_th.get(),
                'scan_l_th' : self.scan_l_th.get(),
                'chunksize' : self.chunksize.get(),
                'max_seg_len' : self.max_seg_len.get(),
                'min_seg_len' : self.min_seg_len.get(),
                }

        with open("settings", 'w') as settingsfile:
            cPickle.dump(data, settingsfile)
            print 'Settings saved'

    def restoreData(self):
        # Return an unpickled dictionary with stored variables 
        with open("settings", 'r') as settingsfile:
            data = cPickle.load(settingsfile)
            return data
 
    def filename_browse(self):
        # Browse to select filename
        data.filename = askopenfilename(title='Select I(s) scan...',
                filetypes=[('Text files', '*.txt'), ('Gzipped files', '*.gz'), ('All files', '*')])

    def quit(self):
        self.ans = askokcancel('Verify exit', 'Are you sure?')
        if self.ans:
            quit(self)

    def recon_scan_info(self):
        # Displays the nfo.txt for the dataset
        #TODO: There is a minor bug here, in that the program only looks in the current directory for nfo.txt
        # If you're analysing data in a different folder to location of scanmaster it won't find the file. 
        if os.access('nfo.txt', os.F_OK):
    
            # Configure the data input parameters in a new window
            self.recon_params = Toplevel()
            self.recon_params.title('Experiment notes')
            self.recon_frame = Frame(self.recon_params)
            self.recon_frame.pack(side=TOP, padx=50, pady=5)
            self.text = Text(self.recon_frame, width=100, height=20)
            self.text.grid()
            with open('nfo.txt', 'r') as FILE:
                lines = FILE.readlines()
                number = 1.0
                for line in lines:
                    self.text.insert(number, line)
                    number += 1.0
        else:
            #Display an error if file isn't there
            self.error = showerror('Error', 'nfo.txt not found')

    def chopper_params(self):

        # Configure chopper parameters
        self.chopper_params = Toplevel()
        self.chopper_params.title('Chopper parameters')

        # Put the parameters in a frame
        self.chopper_frame = Frame(self.chopper_params)
        self.chopper_frame.pack(side=TOP, padx=50, pady=5)

        # Draw the controls; see the text for what they do
        Label(self.chopper_frame, text='Minimum segment length (points)').grid(row=0,
                column=0, sticky=W)
        self.minseglength = Spinbox(
            self.chopper_frame,
            from_=1,
            to=10000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.min_seg_len,
            )
        self.minseglength.grid(row=0, column=1)

        Label(self.chopper_frame, text='Maximum segment length (points)').grid(row=1,
                column=0, sticky=W)
        self.maxseglength = Spinbox(
            self.chopper_frame,
            from_=1,
            to=10000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.max_seg_len,
            )
        self.maxseglength.grid(row=1, column=1)

        Label(self.chopper_frame, text='Scan level upper threshold (V)').grid(row=2,
                column=0, sticky=W)
        self.u_th = Spinbox(
            self.chopper_frame,
            from_=-10.0,
            to=10.0,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.scan_u_th,
            )
        self.u_th.grid(row=2, column=1)

        Label(self.chopper_frame, text='Scan level lower threshold (V)').grid(row=3,
                column=0, sticky=W)
        self.l_th = Spinbox(
            self.chopper_frame,
            from_=-10.0,
            to=10.0,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.scan_l_th,
            )
        self.l_th.grid(row=3, column=1)

        Label(self.chopper_frame, text='Raw data chunk input size').grid(row=4,
                column=0, sticky=W)
        self.l_th = Spinbox(
            self.chopper_frame,
            from_=0,
            to=10000000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.chunksize,
            )
        self.l_th.grid(row=4, column=1)


    def quadchannel_params(self):

        # Configure the quad channel system in a new window
        self.qc_params = Toplevel()
        self.qc_params.title('Quad channel module calibration')

        # Put the parameters in frames
        self.resistance_frame = Frame(self.qc_params)
        self.resistance_frame.pack(side=TOP, padx=70, pady=5)
        self.cutoff_frame = Frame(self.qc_params)
        self.cutoff_frame.pack(side=TOP, padx=70, pady=5)

        # Feedback resistor settings
        Label(self.resistance_frame, text='Feedback resistance (L):'
              ).grid(row=0, column=0, sticky=W)

        Spinbox(
            self.resistance_frame,
            from_=0,
            to=100000000000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.lowres,
            ).grid(row=0, column=1)

        Label(self.resistance_frame, text='Feedback resistance (H):'
              ).grid(row=1, column=0, sticky=W)

        Spinbox(
            self.resistance_frame,
            from_=0,
            to=100000000000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.highres,
            ).grid(row=1, column=1)

        Label(self.resistance_frame, text='Gain Factor'
              ).grid(row=2, column=0, sticky=W)

        Spinbox(
            self.resistance_frame,
            from_=1,
            to=100,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gainfactor,
            ).grid(row=2, column=1)

        # Cutoff voltage settings for channel switchover
        Label(self.cutoff_frame, text='Cutoff voltages (V)'
              ).grid(row=0, column=0, sticky=W)

        Label(self.cutoff_frame, text='Channel 1:').grid(row=1,
                column=0, sticky=W)

        Spinbox(
            self.cutoff_frame,
            from_=0,
            to=10,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.th_1,
            ).grid(row=1, column=1)

        Label(self.cutoff_frame, text='Channel 2:').grid(row=2,
                column=0, sticky=W)

        Spinbox(
            self.cutoff_frame,
            from_=0,
            to=10,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.th_2,
            ).grid(row=2, column=1)

        Label(self.cutoff_frame, text='Channel 3:').grid(row=3,
                column=0, sticky=W)

        Spinbox(
            self.cutoff_frame,
            from_=0,
            to=10,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.th_3,
            ).grid(row=3, column=1)

    def kde_params(self):

        # Configure the KDE parameters in a new window
        self.kde_params = Toplevel()
        self.kde_params.title('Current density plot')

        # Put the parameters in a frame
        self.kde_frame = Frame(self.kde_params)
        self.kde_frame.pack(side=TOP, pady=5)

        Label(self.kde_frame, text='Plot dimensions:', pady=10).grid(row=0,
                column=0, sticky=W)

        Label(self.kde_frame, text='x-axis minimum:').grid(row=1,
                column=0, sticky=W)
        self.xmin = Spinbox(
            self.kde_frame,
            from_=-100000,
            to=100000,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.xmin_cd,
            )
        self.xmin.grid(row=1, column=1)

        Label(self.kde_frame, text='x-axis maximum:').grid(row=2,
                column=0, sticky=W)
        self.xmax = Spinbox(
            self.kde_frame,
            from_=0,
            to=100000,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.xmax_cd,
            )
        self.xmax.grid(row=2, column=1)

        Label(self.kde_frame, text='y-axis minimum:').grid(row=3,
                column=0, sticky=W)
        self.ymin = Spinbox(
            self.kde_frame,
            from_=-10,
            to=10,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.ymin_cd,
            )
        self.ymin.grid(row=3, column=1)

        Label(self.kde_frame, text='y-axis maximum:').grid(row=4,
                column=0, sticky=W)
        self.ymax = Spinbox(
            self.kde_frame,
            from_=-10,
            to=10,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.ymax_cd,
            )
        self.ymax.grid(row=4, column=1)

        Label(self.kde_frame, text='Histogram parameters:', pady=10).grid(row=5,
                column=0, sticky=W)

        Label(self.kde_frame, text='Bin count:').grid(row=6,
                column=0, sticky=W)
        self.xbin = Spinbox(
            self.kde_frame,
            from_=0,
            to=1000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.bins_cd,
            )
        self.xbin.grid(row=6, column=1)

        Label(self.kde_frame, text='KDE parameters:', pady=10).grid(row=7,
                column=0, sticky=W)
    
        Label(self.kde_frame, text='Bandwidth:').grid(row=8,
                column=0, sticky=W)
        self.bandwidth = Spinbox(
            self.kde_frame,
            from_=-10,
            to=10,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.kde_bandwidth,
            )
        self.bandwidth.grid(row=8, column=1)

        Label(self.kde_frame, text='Stop position:').grid(row=9,
                column=0, sticky=W)
        self.stop = Spinbox(
            self.kde_frame,
            from_=-10,
            to=10,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.kde_stop,
            )
        self.stop.grid(row=9, column=1)

        Label(self.kde_frame, text='Points:').grid(row=10,
                column=0, sticky=W)
        self.points = Spinbox(
            self.kde_frame,
            from_=0,
            to=1000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.kde_points,
            )
        self.points.grid(row=10, column=1)

        Label(self.kde_frame, text='KDE auto peak maximum labelling:', pady=10).grid(row=11,
                column=0, sticky=W)

        Label(self.kde_frame, text='Attempt to auto label peak maximum:').grid(row=12,
                column=0, sticky=W)
        self.peakMax = Checkbutton(self.kde_frame, variable=self.kdePeakLabel)
        self.peakMax.grid(row=12, column=1)

        Label(self.kde_frame, text='Label minimum (log10[nA]):').grid(row=13,
                column=0, sticky=W)
        self.peakStart = Spinbox(
            self.kde_frame,
            from_=-10,
            to=10,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.kdePeakStart,
            )
        self.peakStart.grid(row=13, column=1)

        Label(self.kde_frame, text='Label maximum (log10[nA]):').grid(row=14,
                column=0, sticky=W)
        self.peakStop = Spinbox(
            self.kde_frame,
            from_=-10,
            to=10,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.kdePeakStop,
            )
        self.peakStop.grid(row=14, column=1)

    def data_filter(self):

        # Configure the data input parameters in a new window
        self.data_filter = Toplevel()
        self.data_filter.title('Data filtering')

        # Put the parameters in a frame
        self.data_frame = Frame(self.data_filter)
        self.data_frame.pack(side=TOP, padx=50, pady=5)

        # Draw the controls; see the text for what they do
        Label(self.data_frame, text='Cumulative log(I) threshold:'
              ).grid(row=0, column=0, sticky=W)
        self.start = Spinbox(
            self.data_frame,
            from_=-10000,
            to=10000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.datfillogi,
            )
        self.start.grid(row=0, column=1)

    def data_params(self):

        # Configure the data input parameters in a new window
        self.data_params = Toplevel()
        self.data_params.title('Data parameters')

        # Put the parameters in a frame
        self.data_frame = Frame(self.data_params)
        self.data_frame.pack(side=TOP, padx=50, pady=5)

        # Draw the controls; see the text for what they do
        Label(self.data_frame, text='Manual file input range start:'
              ).grid(row=0, column=0, sticky=W)
        self.start = Spinbox(
            self.data_frame,
            from_=1,
            to=10000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.stavar,
            )
        self.start.grid(row=0, column=1)

        Label(self.data_frame, text='Manual file input range finish:'
              ).grid(row=1, column=0, sticky=W)
        self.finish = Spinbox(
            self.data_frame,
            from_=1,
            to=10000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.finvar,
            )
        self.finish.grid(row=1, column=1)

        Label(self.data_frame, text='Background correction factor:'
              ).grid(row=2, column=0, sticky=W)
        self.bcorfactor = Spinbox(
            self.data_frame,
            from_=0.05,
            to=0.20,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.bcorfac,
            )
        self.bcorfactor.grid(row=2, column=1)

        Label(self.data_frame, text='Scan distance cutoff (nm):'
              ).grid(row=3, column=0, sticky=W)
        self.xfactor = Spinbox(
            self.data_frame,
            from_=0.1,
            to=6.0,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.xfac,
            )
        self.xfactor.grid(row=3, column=1)

        Label(self.data_frame, text='Current offset for OTF plot (nA):'
              ).grid(row=4, column=0, sticky=W)
        self.goffset = Spinbox(
            self.data_frame,
            from_=0.0,
            to=1000.0,
            increment=1.0,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.offset,
            )
        self.goffset.grid(row=4, column=1)

        Label(self.data_frame, text='Current scale factor:'
              ).grid(row=5, column=0, sticky=W)
        self.scalefac = Spinbox(
            self.data_frame,
            from_=0.0,
            to=10000000000.0,
            increment=1.0,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.currentScaleFactor,
            )
        self.scalefac.grid(row=5, column=1)

        Label(self.data_frame, text='y-axis limit for BJ plot:'
              ).grid(row=6, column=0, sticky=W)
        self.yfactor = Spinbox(
            self.data_frame,
            from_=0.1,
            to=50000.0,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.yfac,
            )
        self.yfactor.grid(row=6, column=1)

        Label(self.data_frame, text='Tip sample bias (V):'
              ).grid(row=7, column=0, sticky=W)
        self.Bias = Spinbox(
            self.data_frame,
            from_=0.1,
            to=10.0,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.bias,
            )
        self.Bias.grid(row=7, column=1)

    def adc_params(self):

        # Configure the ADC / STM parameters
        self.adc_params = Toplevel()
        self.adc_params.title('ADC / STM parameters')

        # Put the parameters in a frame
        self.adc_frame = Frame(self.adc_params)
        self.adc_frame.pack(side=TOP, padx=50, pady=5)

        # Draw the controls; see the text for what they do
        Label(self.adc_frame, text='Samples per second:').grid(row=0,
                column=0, sticky=W)
        self.sampsecond = Spinbox(
            self.adc_frame,
            from_=1,
            to=10000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.sampsec,
            )
        self.sampsecond.grid(row=0, column=1)

        Label(self.adc_frame, text='Scan range (nm):').grid(row=1,
                column=0, sticky=W)
        self.scanrange = Spinbox(
            self.adc_frame,
            from_=0.1,
            to=15.0,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.srange,
            )
        self.scanrange.grid(row=1, column=1)

        Label(self.adc_frame, text='Sweep duration (s):').grid(row=2,
                column=0, sticky=W)
        self.scanduration = Spinbox(
            self.adc_frame,
            from_=0.1,
            to=1.0,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.sduration,
            )
        self.scanduration.grid(row=2, column=1)

    def contour_params(self):

        # Configure the contour parameters
        self.contour_params = Toplevel()
        self.contour_params.title('2D histogram parameters')

        # Put the parameters in a frame
        self.contour_frame = Frame(self.contour_params)
        self.contour_frame.pack(side=TOP, padx=50, pady=5)

        # Draw the controls; see the text for what they do
        Label(self.contour_frame, text='x-axis minimum:').grid(row=0,
                column=0, sticky=W)
        self.xmin = Spinbox(
            self.contour_frame,
            from_=-100000,
            to=100000,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.xmin_contour,
            )
        self.xmin.grid(row=0, column=1)

        Label(self.contour_frame, text='x-axis maximum:').grid(row=1,
                column=0, sticky=W)
        self.xmax = Spinbox(
            self.contour_frame,
            from_=0,
            to=100000,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.xmax_contour,
            )
        self.xmax.grid(row=1, column=1)

        Label(self.contour_frame, text='y-axis minimum:').grid(row=2,
                column=0, sticky=W)
        self.ymin = Spinbox(
            self.contour_frame,
            from_=-10,
            to=10,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.ymin_contour,
            )
        self.ymin.grid(row=2, column=1)

        Label(self.contour_frame, text='y-axis maximum:').grid(row=3,
                column=0, sticky=W)
        self.ymax = Spinbox(
            self.contour_frame,
            from_=-10,
            to=10,
            increment=0.1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.ymax_contour,
            )
        self.ymax.grid(row=3, column=1)

        Label(self.contour_frame, text='x bin count:').grid(row=4,
                column=0, sticky=W)
        self.xbin = Spinbox(
            self.contour_frame,
            from_=0,
            to=1000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.xbin_contour,
            )
        self.xbin.grid(row=4, column=1)

        Label(self.contour_frame, text='y bin count:').grid(row=5,
                column=0, sticky=W)
        self.ybin = Spinbox(
            self.contour_frame,
            from_=0,
            to=1000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.ybin_contour,
            )
        self.ybin.grid(row=5, column=1)

    def correlation_params(self):

        # Configure the contour parameters
        self.correlation_params = Toplevel()
        self.correlation_params.title('2D correlation histogram parameters')

        # Put the parameters in a frame
        self.correlation_frame = Frame(self.correlation_params)
        self.correlation_frame.pack(side=TOP, padx=50, pady=5)

        # Draw the controls; see the text for what they do
        Label(self.correlation_frame, text='Current minimum [Log10(nA)] :').grid(row=0,
                column=0, sticky=W)
        self.currentmin = Spinbox(
            self.correlation_frame,
            from_=--10,
            to=10,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.corrcurrentmin,
            )
        self.currentmin.grid(row=0, column=1)

        Label(self.correlation_frame, text='Current maximum [Log10(nA)] :').grid(row=1,
                column=0, sticky=W)
        self.currentmax = Spinbox(
            self.correlation_frame,
            from_=-10,
            to=10,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.corrcurrentmax,
            )
        self.currentmax.grid(row=1, column=1)

        Label(self.correlation_frame, text='Number of bins:').grid(row=2,
                column=0, sticky=W)    
        self.corrnumbins = Spinbox(
            self.correlation_frame,
            from_=10,
            to=1000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.corrbins,
            )
        self.corrnumbins.grid(row=2, column=1)

        Label(self.correlation_frame, text='Use log scale:').grid(row=3,
                column=0, sticky=W)
        self.logscale = Checkbutton(self.correlation_frame, variable=self.corrlogscale)
        self.logscale.grid(row=3, column=1)

    def gaussian_params(self):

        self.gauss_params = Toplevel()
        self.gauss_params.title('Gaussian fitting parameters')
        self.gauss_frame = Frame(self.gauss_params)
        self.gauss_frame.pack(side=TOP, padx=50, pady=5)

        
            
        Label(self.gauss_frame, text='Gaussian 1 settings:', pady=10).grid(row=0,
                column=0, sticky=W)

        Label(self.gauss_frame, text='Mu 1:').grid(row=1,
                column=0, sticky=W)
        self.mu1 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_mu1,
            )
        self.mu1.grid(row=1, column=1)

        Label(self.gauss_frame, text='Sigma 1:').grid(row=2,
                column=0, sticky=W)
        self.sigma1 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_sigma1,
            )
        self.sigma1.grid(row=2, column=1)

        Label(self.gauss_frame, text='Scale factor 1:').grid(row=3,
                column=0, sticky=W)
        self.scale1 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_scale1,
            )
        self.scale1.grid(row=3, column=1)

        Label(self.gauss_frame, text='y-disp 1:').grid(row=4,
                column=0, sticky=W)
        self.ydisp1 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_ydisp1,
            )
        self.ydisp1.grid(row=4, column=1)

        Label(self.gauss_frame, text='Gaussian 2 settings:', pady=10).grid(row=5,
                column=0, sticky=W)

        Label(self.gauss_frame, text='Mu 2:').grid(row=6,
                column=0, sticky=W)
        self.mu2 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_mu2,
            )
        self.mu2.grid(row=6, column=1)

        Label(self.gauss_frame, text='Sigma 2:').grid(row=7,
                column=0, sticky=W)
        self.sigma2 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_sigma2,
            )
        self.sigma2.grid(row=7, column=1)

        Label(self.gauss_frame, text='Scale factor 2:').grid(row=8,
                column=0, sticky=W)
        self.scale2 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_scale2,
            )
        self.scale2.grid(row=8, column=1)

        Label(self.gauss_frame, text='y-disp 2:').grid(row=9,
                column=0, sticky=W)
        self.ydisp2 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_ydisp2,
            )
        self.ydisp2.grid(row=9, column=1)

        Label(self.gauss_frame, text='Gaussian 3 settings:', pady=10).grid(row=10,
                column=0, sticky=W)

        Label(self.gauss_frame, text='Mu 3:').grid(row=11,
                column=0, sticky=W)
        self.mu3 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_mu3,
            )
        self.mu3.grid(row=11, column=1)

        Label(self.gauss_frame, text='Sigma 3:').grid(row=12,
                column=0, sticky=W)
        self.sigma3 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_sigma3,
            )
        self.sigma3.grid(row=12, column=1)

        Label(self.gauss_frame, text='Scale factor 3:').grid(row=13,
                column=0, sticky=W)
        self.scale3 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_scale3,
            )
        self.scale3.grid(row=13, column=1)

        Label(self.gauss_frame, text='y-disp 3:').grid(row=14,
                column=0, sticky=W)
        self.ydisp3 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_ydisp3,
            )
        self.ydisp3.grid(row=14, column=1)

        Label(self.gauss_frame, text='Gaussian 4 settings:', pady=10).grid(row=15,
                column=0, sticky=W)

        Label(self.gauss_frame, text='Mu 4:').grid(row=16,
                column=0, sticky=W)
        self.mu4 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_mu4,
            )
        self.mu4.grid(row=16, column=1)

        Label(self.gauss_frame, text='Sigma 4:').grid(row=17,
                column=0, sticky=W)
        self.sigma4 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_sigma4,
            )
        self.sigma4.grid(row=17, column=1)

        Label(self.gauss_frame, text='Scale factor 4:').grid(row=18,
                column=0, sticky=W)
        self.scale4 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_scale4,
            )
        self.scale4.grid(row=18, column=1)

        Label(self.gauss_frame, text='y-disp 4:').grid(row=19,
                column=0, sticky=W)
        self.ydisp4 = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_ydisp4,
            )
        self.ydisp4.grid(row=19, column=1)

        Label(self.gauss_frame, text='Axis limits:', pady=10).grid(row=20,
                column=0, sticky=W)
    
        Label(self.gauss_frame, text='x-axis minimum:').grid(row=21,
                column=0, sticky=W)
        self.xmin = Spinbox(
            self.gauss_frame,
            from_=-10,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_xmin,
            )
        self.xmin.grid(row=21, column=1)

        Label(self.gauss_frame, text='x-axis maximum:').grid(row=22,
                column=0, sticky=W)
        self.xmax = Spinbox(
            self.gauss_frame,
            from_=-10,  #FIXME put here xmin?
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_xmax,
            )
        self.xmax.grid(row=22, column=1)

        Label(self.gauss_frame, text='y-axis minimum:').grid(row=23,
                column=0, sticky=W)
        self.ymin = Spinbox(
            self.gauss_frame,
            from_=0,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_ymin,
            )
        self.ymin.grid(row=23, column=1)

        Label(self.gauss_frame, text='y-axis maximum:').grid(row=24,
                column=0, sticky=W)
        self.ymax = Spinbox(
            self.gauss_frame,
            from_=0,
            to=10,
            increment=0.001,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.gF_ymax,
            )
        self.ymax.grid(row=24, column=1)

    def plateau_fitting_params(self):

        # Configure the contour parameters
        self.plateau_fitting_params = Toplevel()
        self.plateau_fitting_params.title('Plateau fitting parameters')

        # Put the parameters in a frame
        self.plateau_fitting_frame = Frame(self.plateau_fitting_params)
        self.plateau_fitting_frame.pack(side=TOP, padx=50, pady=5)

        # Draw the controls; see the text for what they do
        Label(self.plateau_fitting_frame,
              text='Minimum data points per plateau:').grid(row=0,
                column=0, sticky=W)
        self.xmin = Spinbox(
            self.plateau_fitting_frame,
            from_=6,
            to=30,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.min_points_plat,
            )
        self.xmin.grid(row=0, column=1)

        Label(self.plateau_fitting_frame,
              text='Maximum data points per plateau:').grid(row=1,
                column=0, sticky=W)
        self.xmax = Spinbox(
            self.plateau_fitting_frame,
            from_=10,
            to=100,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.max_points_plat,
            )
        self.xmax.grid(row=1, column=1)

        Label(self.plateau_fitting_frame,
              text='Background tolerence level (standard deviations)'
              ).grid(row=2, column=0, sticky=W)
        self.ymin = Spinbox(
            self.plateau_fitting_frame,
            from_=0,
            to=100,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.background_tol,
            )
        self.ymin.grid(row=2, column=1)

        Label(self.plateau_fitting_frame,
              text='Maximum plateau gradient (A/m):').grid(row=3,
                column=0, sticky=W)
        self.ymax = Spinbox(
            self.plateau_fitting_frame,
            from_=-1000000,
            to=1000000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.plat_max_grad,
            )
        self.ymax.grid(row=3, column=1)

        Label(self.plateau_fitting_frame,
              text='Maximum plateau current (nA):').grid(row=4,
                column=0, sticky=W)
        self.xbin = Spinbox(
            self.plateau_fitting_frame,
            from_=0,
            to=100000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.max_plat_cur,
            )
        self.xbin.grid(row=4, column=1)

        Label(self.plateau_fitting_frame,
              text='Maximum plateau deviation from the average:'
              ).grid(row=5, column=0, sticky=W)
        self.ybin = Spinbox(
            self.plateau_fitting_frame,
            from_=0.00,
            to=10000.0,
            increment=0.01,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.fractional_plat_tol,
            )
        self.ybin.grid(row=5, column=1)

    def kdePlot(self, savefig=False):

        # Fetch KDE parameters from control widgets
        xmin = self.xmin_cd.get()
        xmax = self.xmax_cd.get()
        ymin = self.ymin_cd.get()
        ymax = self.ymax_cd.get()
        bins = self.bins_cd.get()
        kde_bandwidth = self.kde_bandwidth.get()
        kde_start = kde_bandwidth  # The lowest the KDE start can go
        kde_stop = self.kde_stop.get()
        kde_points = self.kde_points.get()

        # Plot KDE function and read all current data to a single list
        i_list = []
        for reading in data.i_dat_all_combined:
            for value in reading:
                #FIXME You wouldn't normally expect negative currents after badkground subtraction
                #but overshoot in the op amp can sometimes cause them. Need to chuck them away
                #as opposed to takings the abs. value to avoid creating a false peak!
                #This has the effect of vastly reducing the background peak!
                if value > 0:
                    i_list.append(np.log10(value))
        print 'Lists appended, found:', len(i_list), 'data points.'

        #Check that some data is in memory before attempting to plot
        if (len(i_list) < 1):
            self.error = showerror('Error', 'No data in memory')
            return
        
        #This bit is experimental -- plot 5 KDES as a function of scan number
        #numScans = len(data.i_dat_all_combined)    
        #if numScans > 100:
            #Make some current lists
        #    iListA = data.i_dat_all_combined[scansPerList*(scanLists-5):scansPerList*(scanLists-4)]
        #    iListB = data.i_dat_all_combined[scansPerList*(scanLists-4):scansPerList*(scanLists-3)]
        #    iListC = data.i_dat_all_combined[scansPerList*(scanLists-3):scansPerList*(scanLists-2)]
        #    iListD = data.i_dat_all_combined[scansPerList*(scanLists-2):scansPerList*(scanLists-1)]
        #    iListE = data.i_dat_all_combined[scansPerList*(scanLists-1):scansPerList*(scanLists-0)]
            #Divide up the scans
        #    scanLists = 5
        #    scansPerList = numScans / scanLists

        #Fit KDE and setup plot
        #FIXME: Scrap the optimal bandwidth for now.. does it make sense on log scale?
        #haa = statistics.bandwidth(i_list, kernel='Epanechnikov')
        #bandwidth = haa /2
        #print "Optimal bandwidth:", haa
        x = np.linspace(min(i_list),max(i_list),kde_points)
        z = statistics.pdf(i_list, x, h=kde_bandwidth, kernel='E')
       
        #Autolabelling of the KDE peak maximum in a certain window. 
        #Crude, but works well for some molecules as a first pass estimate
        if self.kdePeakLabel.get():
            flag = False
            flag2 = False
            for i in range(len(x)):
                if (x[i] > self.kdePeakStart.get()) and (flag == False):
                    sliceStart = i
                    flag = True
                if (x[i] > self.kdePeakStop.get()) and (flag2 ==False):
                    sliceEnd = i 
                    flag2 = True
            zSS = z[sliceStart:sliceEnd]
            xSS = x[sliceStart:sliceEnd]
            print zSS.argmax(), xSS[zSS.argmax()]
            HC = xSS[zSS.argmax()]
  
        #Setup plot
        plt.hist(i_list, bins=bins, facecolor='black', normed=1)
        plt.plot(x, z,'b', label='KDE', linewidth=3)
        plt.axis([xmin, xmax, ymin, ymax])
        plt.legend(bbox_to_anchor=(1.05, 1),loc=2,borderaxespad=0.)
        plt.grid()
        plt.ylabel(r'$Density$', fontsize=18)
        plt.xlabel(r'$log_{10}(G/G_{0}$)', fontsize=18)       
        if self.kdePeakLabel.get():
            plt.annotate(str(HC)[0:5], (HC, max(zSS)*1.2),
                xytext=(0.5, 0.5), textcoords='axes fraction',
                arrowprops=dict(facecolor='blue', shrink=0.05),
                fontsize=14,
                horizontalalignment='right', verticalalignment='top')
        title = data.filename[0:string.rfind(data.filename, '/')]
        plt.title(title, fontsize=16)
        if savefig:
            plt.savefig(title+".png", format='png') # can also be pdf etc..
            plt.cla() 
        else:
            plt.show()

    def linearPlot(self):

        # Fetch KDE parameters from control widgets
        i_list = []
        for reading in data.i_dat_all_combined:
            for value in reading:
                # i_list.append(np.log10(float(value)))
                # i_list.append(value)
                if ( (value > 0) ): #and (value < 0.5) ):
                    i_list.append(value)
        
        if (len(i_list) < 1):
            self.error = showerror('Error', 'No data in memory')
            return

        print 'Lists appended, found:', len(i_list), 'data points.'
        plt.hist(i_list, bins=2000, facecolor='black', normed=1)
        plt.grid()
        plt.ylabel(r'$Density$', fontsize=18)
        plt.xlabel(r'$G/G_{0}$', fontsize=18)
        plt.show()

    def inputRange(self):
        """Returns the start/finish numbers for the I(s) scans to read in"""
        if (len(data.filename) < 1):
            self.error = showerror('Error', 'Input folder not defined')
            return
        # Reset the new file counter so autoprocessed / manually processed files are labelled
        # from zero everytime.
        data.file_newnum = 0
        # Set the number of files to read in automatically
        auto_read = controller.auto_read.get()
        if auto_read == 1:
            # Automatically set file input range if so desired by the user
            target_folder = data.filename[0:string.rfind(data.filename, '/')]
            # Have a poke around the user selected data directory for the number of files
            # which must be sequentially numbered
            files = os.listdir(os.path.abspath(target_folder))
            filenumber_list = []
            for filename in files:
                filenumber_list.append(int(filename[-8:-4]))
            start = min(filenumber_list)
            finish = max(filenumber_list)
        else:
            # Use the default or user set input range
            start = int(self.stavar.get())
            finish = int(self.finvar.get())
        return start, finish

    def readFiles(self):
        """ Read in the I(s) scans, stitch together and correct the background """
        start, finish = self.inputRange()
        bcorfac = self.bcorfac.get()
        #Clear the current lists if you don't want to combine data sets
        if controller.clear_global_data.get() == 1:
            data.clear_current_lists()
        datInput(start, finish, bcorfac)
        print 'End of file input'

    def autoPlot(self):
        """ Read in the I(s) scans in all "chopped" folders and plot graphs on fixed scale """ 
        print "Automatically generating current density plots..."
        contents = os.listdir(os.getcwd())
        dirList = []
        for item in contents:
            if (len(item) >= 8) and (item[:7] == "chopped") and (item[-4:] != ".png"):
                dirList.append(item)
        for directory in dirList:
            scans = os.listdir(os.path.abspath(directory))
            if len(scans) == 0: return
            data.filename = os.path.abspath(directory) + '/' + scans[0]
            start, finish = self.inputRange()
            bcorfac = self.bcorfac.get()
            #Clear the current lists if you don't want to combine data sets
            if controller.clear_global_data.get() == 1:
                data.clear_current_lists()
            datInput(start, finish, bcorfac)
            print "Generating current density plot..."
            #self.kdePlot(True)
            #contourPlot(True)
            #print "Syncing scans"
            #plat_sync()
            correlationHist(True)
        print "Automated plotting complete."
        
    def plat_syncer(self):
        """ Sync the I(s) scans to the start or end of a G0 plateau """
        start, finish = self.inputRange()
        plat_sync()

    def importdata(self):
        # Needs to call tea break maker
        tea_break_maker()
   
root = Tk()
root.title('Scanmaster 3000 v0.59')
egraph = egraph(root)
#Create the data container
data = Data()
controller = controller(root)
root.mainloop()
