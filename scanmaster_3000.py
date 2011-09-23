#!/usr/bin/python
# -*- coding: utf-8 -*-

#     Doug S. Szumski  <d.s.szumski@gmail.com>  01-08-2011
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

# Dialogue

from Tkinter import *
from tkFileDialog import askopenfilename, asksaveasfile
from tkColorChooser import askcolor
from tkMessageBox import askquestion, showerror, showinfo
from tkSimpleDialog import askfloat
from tkMessageBox import askokcancel
from FileDialog import LoadFileDialog
from Dialog import Dialog

# File IO

import shutil  # High level file from scipy import stats
import os  # import command line for file operations
from subprocess import call
import string
import cPickle

# KDE stuff

import numpy as np
from scipy import stats as sci_stats
import matplotlib.pylab as plt
import statistics
from scipy import linspace, polyval, polyfit, sqrt, stats, randn

# Retro Fortran modules -- recompile using: f2py -c -m plat_seek plat_seek8.f90 etc..

import plat_seek


class data:
    #Initialise a load of lists / variables for global use
    def __init__(self):
        # Sampling distance interval - calculated in scaninput for each file
        self.interval = 0.00
        # Four current channels from the quad amp
        self.i_dat_all_ls = []
        self.i_dat_all_hs = []
        self.i10_dat_all_ls = []
        self.i10_dat_all_hs = []
        # Plateau data
        self.plat_dat_ls = []
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


def scaninput(name):

    # Reads in the file, filename, reads number of measurements and then extracts I(s) data into two lists
    # Read data from input file

    infile = open(name, 'r')
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
    # This is the distance the tip travels between data points

    data.interval = srange / sduration / sampsec
    position = 0.00
    for line in lines:
        line = line.split()
        i_dat_ls.append((float(line[0])))  # ignore current sign, correct for resistor /v
        i10_dat_ls.append((float(line[1])))
        i_dat_hs.append((float(line[2])))
        i10_dat_hs.append((float(line[3])))
        s_dat_ls.append(position)
        position += data.interval
    return (i_dat_ls, i10_dat_ls, i_dat_hs, i10_dat_hs, s_dat_ls)
    infile.close


def contour_plot():

    # Plots 2d histograms of current distance scans
    # Generate list of all currents in the current distance scans

    i_list = []
    for i_dat in data.i_dat_all_combined:
        for value in i_dat:
            i_list.append(np.log10(abs(value)))
            #i_list.append(abs(value)) #uncomment for linear plot
    if (len(i_list) < 1):
        error = showerror('Error', 'No data in memory')
        return

    # Generate list of all distances. This list will be the same length as the current list and therefore
    # the indices will be directly related.
    s_list = []
    for s_dat in data.s_dat_all_ls:
        for value in s_dat:
            s_list.append(value)

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
    plt.imshow(H, origin='lower', extent=extent, interpolation='nearest')
    cb = plt.colorbar()
    cb.set_label('counts (normalised)')
    plt.title('I(s) scan 2D histogram')
    plt.xlabel('Distance (nm)')
    plt.ylabel('log10[current (nA)]')
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
    FILE = asksaveasfile(mode='w', title='Save current list as...')
    if FILE:
        for value in i_list:
            FILE.write('%s \n' % value)
        info = showinfo('Notice', 'File saved')
        FILE.close()


def plat_seeker(current_trace):

    # Search data in current_trace for plateaus

    # Fetch variables from GUI for plateau seek

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

    # Did we find a plateau? - more efficient way to do this is in the Fortran module

    if sum(p_avg) > 0.00:
        locplat = True
    else:
        locplat = False

    # Return list of length current_trace, but with only the plateau average non-zero (if found)

    return (p_avg, locplat)


def dat_input(start, finish, bcorfac):

    # Generates sequential filenames, currently in the format filenameXXXX.ivs where the number XXXX is set in the GUI.
    # For each filename it then calls scaninput (see above), corrects the current data for a background offset
    # by averaging over the final fraction of the data set in the GUI, then appends the I(s) data to a list

    file_prefix = data.filename[0:-8]  # includes file path 4 digits
    # file_prefix = data.filename[0:-7] #includes file path 3 digits
    file_ext = '.txt'
    file_number = data.filename[-8:-4]  # 4 digitis
    # file_number = data.filename[-7:-4] # 3 digits
    #FIXME: Delete: initiliaised in the data class
    data.i_dat_all_ls = []  # current
    data.i10_dat_all_ls = []
    data.i_dat_all_hs = []
    data.i10_dat_all_hs = []
    data.s_dat_all_ls = []  # distance
    data.i_dat_all_combined = []
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
            scaninput(auto_name)

        # Work out where to chop the channels whilst the data is stored as a voltage
        # Threshold set from the GUI: Normally these are 0.1V 
        th_1 = controller.th_1.get()
        th_2 = controller.th_2.get()
        th_3 = controller.th_3.get() 

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
            print "Channels HSx1 and HSx10 are saturated and have been ignored"
            sat_flag_hs1 = True
            #It follows that HSx10 must have also saturated
            sat_flag_hs10 = True
        elif (fin_avg_hs1 > 6.0):
            print "Channel HSx10 saturated and has been ignored"
            sat_flag_hs10 = True
            
        #Get the resistor values
        i_ls_res = int(controller.lowres.get())
        i_hs_res = int(controller.highres.get())
        scale = 1e9  # convert to nanoamps

        #Figure out if the scan is negative or not, using sum for average and invert if it is
        #TODO Delete this when sure as it should now be redundant following inversion upon initial data
        # stream reading
        #if (sum(i_dat_ls) < 0.0):
        #    #Invert the measurement
        #    for i in range(len(i_dat_ls)):
        #        i_dat_ls[i] = i_dat_ls[i] * -1.0 
        #        i10_dat_ls[i] = i10_dat_ls[i] * -1.0
        #        i_dat_hs[i] = i_dat_hs[i] * -1.0
        #        i10_dat_hs[i] = i10_dat_hs[i] * -1.0

        #Make a copy of the list (note to self: direct assignment creates a pointer to the list which caused minor hair loss)
        #Important: This contains voltages so they can be compared directly to the thresholds
        i_dat_combined = list(i_dat_ls)

        #Set this to true to use only the LSx1 channel
        single_chan = False
        if single_chan:
            #Test for only one channel
            for i in range(len(i_dat_combined)):
                #if i_dat_combined[i] 
                i_dat_combined[i] = i_dat_combined[i] / (i_ls_res) * scale
        else:
            #Stitch the channels together
            for i in range(len(i_dat_combined)): 
                #FIXME Convert to amps from volts, this used to be in an 'else' but 
                #moved here (less efficient) due to obscure bug?
                i_dat_combined[i] = i_dat_combined[i] / (i_ls_res) * scale
                if (i_dat_combined[i] < th_1):
                    #Output below threshold so try and replace it with a higher sensitivity measurement
                    if (i10_dat_ls[i] > th_2):
                        #LSx10 channel is still in operation so use the measurment from there
                        i_dat_combined[i] = i10_dat_ls[i] / (i_ls_res * 10) * scale
                        #i_dat_combined[i] = i_dat_combined[i] / (i_ls_res) * scale
                    elif ( (i_dat_hs[i] > th_3) and not sat_flag_hs1):
                        #Limiter should be off so check the HSx1 channelfirst
                        i_dat_combined[i] = i_dat_hs[i] / i_hs_res * scale
                    elif not sat_flag_hs10:
                        #If HSx1 is below threshold then use the HSx10 channel
                        i_dat_combined[i] = i10_dat_hs[i] / (i_hs_res * 10) * scale  

        #Now convert the individual channels to currents, could have convoluted this with the above
        #for efficiency at the expense of clarity
        #This isn't required really, but its nice to have on the graph plots for fine tuning
        for i in range(len(i_dat_ls)):
                i_dat_ls[i] = i_dat_ls[i] / i_ls_res * scale
                i10_dat_ls[i] = i10_dat_ls[i] / (i_ls_res * 10) * scale
                #TODO (low priority) If these are saturated they will appear as almost flat lines on the plot
                #Notify the plotter not to plot them to avoid confusion?
                if not sat_flag_hs1: i_dat_hs[i] = i_dat_hs[i] / i_hs_res * scale
                if not sat_flag_hs10: i10_dat_hs[i] = i10_dat_hs[i] / (i_hs_res * 10) * scale 
        
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

            (plat_data, locplat) = plat_seeker(i_dat_ls)
            if locplat:
                print 'PASSED: plateaus were found', data.err
                autosave(auto_name, savedir + '_fil_pltfit')
            else:
                print 'FAILED: no plateaus found'

        # Plot/check if set in the GUI

        if controller.plot_dat.get() > 0:

            # If we haven't already looked for plateaus then we better go do it now so we can plot them on the graph

            if controller.autocheck_pltfit.get() == False:
                (plat_data, locplat) = plat_seeker(i_dat_ls)
            egraph.updater(
                s_dat_ls,
                i_dat_ls,
                i10_dat_ls,
                i_dat_hs,
                i10_dat_hs,
                i_dat_combined,
                plat_data,
                auto_name,
                )
        if controller.check_dat.get() > 0:
            userinput(auto_name)



def back_correct(i_dat, bcorfac):

    # Correct data for background offset
    # Define begin and end for background average

    end = len(i_dat)
    final_int = int(end * bcorfac)  # add correction factor later
    begin = int(end - final_int)

    # Calculate the average of the background

    fin_avg = sum(i_dat[begin:end]) / final_int

    # Subtract the average from every current in the I(s) measurement

    for j in range(0, end):
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

    # Save a file auto_name in a directory with no questions asked

    targetfile = auto_name  # Local copy for recursion
    print 'Saving data...'
    data.file_newnum += 1
    file_ext = '.txt'
    new_name = 'man' + str(data.file_newnum).zfill(4) + file_ext
    pathvar = './' + dirname + '/'  # Copy data to this directory
    if not os.path.isdir(pathvar):  # Create directory if not exist
        os.mkdir(pathvar)
    shutil.copy(targetfile, pathvar + new_name)  # Copy file to pathvar


def rawfileInput(filename):
    print 'Opening: ', filename
    infile = open(filename, 'r')
    i_list_ls_x1 = []
    i_list_ls_x10 = []
    i_list_hs_x1 = []
    i_list_hs_x10 = []
    lines = infile.readlines()
    for line in lines:
        a = line.split()
        i_list_ls_x1.append(float(a[0]))
        i_list_ls_x10.append(float(a[1]))
        i_list_hs_x1.append(float(a[2])) 
        i_list_hs_x10.append(float(a[3]))
    print 'Data loaded...'
    return (i_list_ls_x1, i_list_ls_x10, i_list_hs_x1, i_list_hs_x10)
    infile.close

def groupAvg(data):

    # Returns the average value of a string of data points
    tmp = 0.00
    for value in data:
        tmp += value
    tmp /= len(data)
    # return abs(tmp)
    return tmp


def fileOutput(
    filename,
    data,
    data2,
    data3,
    data4,
    points,
    ):
    FILE = open(filename, 'w')
    for i in range(0, points):
        FILE.write('%s' % data[i] + '\t %s' % data2[i] + '\t %s'
                   % data3[i] + '\t %s \n' % data4[i])
    FILE.close()


def chopper(
    data,
    data_2,
    data_3,
    data_4,
    l_th,
    u_th,
    filecounter,
    ):

    # Chop the continuous data file into individual scans
    # Only data above u_th should be discarded as data below l_th is required for the background correction

    data_window = []
    int_bodge = 20  # The number of points to not include from the determined 'end' (prevents problems with other channels)
    window_length = 50
    l_plat_len = 0
    u_plat_len = 0
    counter = 0
    start = 0
    stop = 0

    #Detect it the tip substrate bias is positive or negative to allow auto-inversion of data
    #The average of the LSx1 channel is a good, but inefficient measure of this
    #TODO: If the warning appears frequently parameterise this in the GUI with a manual override

    dat_avg = np.average(data)
    print "LSx1 channel average is: ", dat_avg, "V"
    if (dat_avg < -0.5):
        print "Assuming negative tip-substrate bias"
        #Invert the data
        #FIXME get rid of the now redundant inversion check when reading the split files back in
        for i in range(len(data)):
                data[i]   = data[i]   * -1.0
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

    for value in data:
        data_window.append(value)
        if len(data_window) > window_length:

            # Get rid off the first value to make way for the last

            del data_window[0]
        if len(data_window) == window_length:

            # Full window so do some checks for lower threshold

            if groupAvg(data_window) < l_th:

                # Found a plateau so increment the length counter

                l_plat_len += 1
            if groupAvg(data_window) > u_th:

                # Found a plateau so increment the length counter

                u_plat_len += 1
            if groupAvg(data_window) < u_th and u_plat_len > 0 and stop \
                == 0:

                # Hopefully this is the tip retraction point and from here on the current decays
                # Stop must be zero otherwise we might end up with background, then decay!

                start = counter  # could make this "counter - u_plat_len" to get all the plateaus but not of interest

                # print "Found upper plateau: ", u_plat_len, " points long, stopping at: ", start

                u_plat_len = 0
            if groupAvg(data_window) > l_th and l_plat_len > 0 \
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
                background_level = groupAvg(data_3[stop - 10:stop])
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
                    data_slice = data[start:stop]
                    data_slice2 = data_2[start:stop]  # Matching one from x10 channel
                    data_slice3 = data_3[start:stop]
                    data_slice4 = data_4[start:stop]

                    # Generate the filename
                    # Later on this can calculate a sdev and dynamically exclude

                    filename = 'slice' + str(filecounter).zfill(4) \
                        + '.txt'
                    print 'Reconstructed I(s) scan containing: ', points, \
                        'data points as: ', filename
                    fileOutput(
                        filename,
                        data_slice,
                        data_slice2,
                        data_slice3,
                        data_slice4,
                        points,
                        )
                    filecounter += 1
                l_plat_len = 0
                start = 0
                stop = 0
        counter += 1
    return filecounter


def tea_break_maker():

    # Deals with reading all raw data files from the ADC and splitting them into I(s) scans using chopper
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

    os.chdir(os.path.abspath('raw'))
    for (path, dirs, files) in os.walk(os.getcwd()):
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
    os.chdir(os.pardir)  # 'cd..

    # Now ask the questions
    if not os.access('nfo.txt', os.F_OK):
        for question in questions:
            var = raw_input(question)
            answers.append(var)

        # Write the questions and answer to file for later viewing,
        # but if the file exists then don't overwrite it

        FILE = open('nfo.txt', 'w')
        for i in range(0, len(answers)):
            info = questions[i] + '\t' + answers[i] + '\n'
            FILE.write(info)
        FILE.close()
    else:
        print "Skipping questions: nfo.txt already exists"

    # Now for the splitting:
    # Assumes you've put the raw data in ../raw/

    print 'This is a good time to have tea break...'

    # Go into the raw data directory
    check_flag = False #This is so the check for exiating files is done once only

    os.chdir(os.path.abspath('raw'))
    for (path, dirs, files) in os.walk(os.getcwd()):
        for raw_data_filename in files:

            # Put the split files in a directory named after the file but without the extension

            dirname = raw_data_filename[0:-4]

            # Make a new folder for the split raw data files and cd into it

            os.chdir(os.pardir)  # cd..
            if os.path.exists(dirname) != True:
                os.mkdir(dirname)
            os.chdir(os.path.abspath(dirname))  # returns absolute path of output and changes to that directory

            # This is the location of the data we want to split up
            
            dataloc = '../raw/' + raw_data_filename
            
            # See if split file aa exists and if it doesn't split up the raw data
            # TODO: make this more generic instead of relying on default output 'aa' from split
            
            if not os.access('aa', os.F_OK):
                print 'Splitting up', raw_data_filename, 'into the directory', dirname, '...'
                # NOTE: distro dependent. Decrease number of lines to reduce RAM requirements, at the expense of
                # destroying scans split between files.  
                call(['split', '--verbose', '-l 1000000', dataloc, ''])
            else:
                print 'Using existing data: target folder not empty'

            # Now start chopping the split raw file

            output = 'chopped' + dirname[6:]
            print 'Reconstructing I(s) scans into the directory', \
                output, '...'

            # The filenumber must be reset for every raw data file

            filecounter = 0
            for (path, dirs, files) in os.walk(os.getcwd()):

                # For every split raw data segment chop it up

                for split_raw_data_segment in files:

                    # Read in the ADC channels to individual lists

                    (i_list_ls_x1, i_list_ls_x10, i_list_hs_x1,
                     i_list_hs_x10) = \
                        rawfileInput(split_raw_data_segment)

                    # Go out of the split raw data folder........

                    os.chdir(os.pardir)

                    # Make the output folder for the I(s) scans if it doesn't exist already and cd into it

                    if os.path.exists(output) != True:
                        os.mkdir(output)
                    os.chdir(os.path.abspath(output))  # returns absolute path of output and change to that directory

                    # Reconstruct the I(s) scans only if the the folder is empty
                    # TODO: make this more generic -- get rid of dependence on specified file        
                    #Time period for processing on one split file and estimate remaining time.  
                    # FIXME: Bug here: Current test only works for first folder! Fix

                    if os.access('slice0000.txt', os.F_OK) and not check_flag:
                        print "Skipping scan reconstruction: target folder:", "../" + output, "is not empty"
                        break                    

                    filecounter = chopper(
                        i_list_ls_x1,
                        i_list_ls_x10,
                        i_list_hs_x1,
                        i_list_hs_x10,
                        controller.scan_l_th.get(),
                        controller.scan_u_th.get(),
                        filecounter,
                        )

                    check_flag = True                    

                    # L_th normally 0.009, u_th 4.0
                    # Go back to the split raw data folder

                    os.chdir(os.pardir)  # cd..
                    os.chdir(os.path.abspath(dirname))
    print 'Finished reconstructing I(s) scans'

    # Go back to the root data folder
    # TODO: All this directory changing works fine, but is a bit confusing to the layman and
    # should really be simplified somehow

    os.chdir(os.pardir)


class egraph:

    # Plots embedded matplotlib graphs in main window
    # NOTE: To edit the x,y range of the plots change xmin, xmax, ymin, ymax below
    def __init__(self, myParent):
        #Deal with different screen resolutions
        if root.winfo_screenwidth() <= 1024:
            #Low res so downsize figures, ideal for XGA, but no smaller (unlikely!)
            self.f = Figure(figsize=(10, 7), dpi=80)
        else:
            #This is fine for SXGA and higher
            self.f = Figure(figsize=(14, 9), dpi=80)
        #Add the subplots
        self.ax = self.f.add_subplot(131)
        self.ax2 = self.f.add_subplot(132)
        self.ax3 = self.f.add_subplot(133)
        self.canvas = FigureCanvasTkAgg(self.f, master=myParent)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=BOTTOM)  # , fill=BOTH, expand=1)

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

        xlim = float(controller.xfac.get())
        self.ax.cla()  # Clear current axes

        # PLOT 1
        # ax.set_title('Current-distance measurements', fontsize=22)

        self.ax.set_yscale('log')
        self.ax.plot(s_dat_ls, i_dat_ls, 'k.', label='LS x1')
        self.ax.plot(s_dat_ls, i10_dat_ls, 'b.', label='LS x10')
        self.ax.plot(s_dat_ls, i_dat_hs, 'g.', label='HS x1')
        self.ax.plot(s_dat_ls, i10_dat_hs, 'r.', label='HS x10')
        self.ax.grid(True)

        # self.ax.set_ylim([1e-8,1e4])

        self.ax.set_xlim([0, xlim])
        self.ax.set_xlabel('Distance (nm)')
        self.ax.set_ylabel('Current (nA)')
        self.ax.legend()

        # PLOT 2

        self.ax2.cla()  # Clear current axes
        self.ax2.set_yscale('log')

        # print i_dat_combined

        self.ax2.plot(s_dat_ls, i_dat_combined, 'k.',
                      label='Combined data')
        if controller.autocheck_linfit.get():
            self.ax2.plot(data.s_dat_filtered, data.polyfit_rescaled,
                          'r', label='MSE:' + str(data.err))
        self.ax2.grid(True)

        # self.ax.set_ylim([1e-8,1e4])

        self.ax2.set_xlim([0, xlim])
        self.ax2.set_xlabel('Distance (nm)')

        # self.ax2.set_ylabel('Current nA')

        self.ax2.legend()

        # PLOT 3

        self.ax3.cla()  # Clear current axes

        # self.ax2.set_yscale('log')
        # print i_dat_combined

        self.ax3.plot(s_dat_ls, i_dat_ls, 'k.-', label='LS x1')
        self.ax3.plot(s_dat_ls, plat_data, 'b.-', label='Plateau fitting')
        self.ax3.grid(True)
        self.ax3.set_ylim([0, 25000])
        self.ax3.set_xlim([0, 2])
        self.ax3.set_xlabel('Distance (nm)')

        # self.ax3.set_ylabel('Current nA')

        self.ax3.legend()

        self.canvas.show()

# Main GUI

class controller:

    # This generates the main GUI and contains all variables the use can modify
    def __init__(self, myParent):
    
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
        self.plot_dat = IntVar()
        self.check_dat = IntVar()
        self.check_dat2 = IntVar()
        self.autocheck_linfit = IntVar()
        self.autocheck_avgcur = IntVar()
        self.autocheck_pltfit = IntVar()
        self.auto_read = IntVar()

         # KDE
        self.kde_bandwidth = DoubleVar()
        self.kde_stop = IntVar()
        self.kde_points = IntVar()

         # Resistor division
        self.lowres = DoubleVar()
        self.highres = DoubleVar()

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
            self.lowres.set(99559)  # Calibrated from resistor measurements, assumes output voltage is correct from STM
            self.highres.set(91474235)  # Calibrated from resistor measurements, assumes output voltage is correct from STM
            self.th_1.set(0.1)  # Ch2 is simply Ch1 x 10 so at this stage Ch2 is at 90% maximum and can take over
            self.th_2.set(0.1)  # Around here the limiter turns off, plus a bit more to stop dodgy overlap
            self.th_3.set(0.1)  # Ch3 is on and can take over
            
            #Data input
            self.stavar.set(1)
            self.finvar.set(5)
            self.bcorfac.set(0.20)
            self.xfac.set(2.00)
            self.plot_dat.set(0)
            self.check_dat.set(0)
            self.check_dat2.set(0)
            self.auto_read.set(0)

            # Data filtering defaults:
            self.datfillogi.set(-1000)

            # KDE stuff:
            self.kde_bandwidth.set(0.1)
            self.kde_stop.set(4)
            self.kde_points.set(1000)
            
            # STM/ADC parameters used for reconstructing distance defaults
            self.sampsec.set(10000)
            self.srange.set(4)
            self.sduration.set(0.3)

            # Contour plot parameter defaults
            self.xmin_contour.set(0.00)
            self.xmax_contour.set(5.0)
            self.ymin_contour.set(-4.0)
            self.ymax_contour.set(4.0)
            self.xbin_contour.set(100)
            self.ybin_contour.set(100)

            # Plateau fitting parameters
            self.plat_max_grad.set(50000.0)
            self.background_tol.set(5)
            self.max_plat_cur.set(10000.0)
            self.max_points_plat.set(100)
            self.fractional_plat_tol.set(0.20)
            self.min_points_plat.set(30)

            #Chopper parameters
            self.scan_u_th.set(4.0)
            self.scan_l_th.set(0.009)
            self.max_seg_len.set(3200)
            self.min_seg_len.set(2200)

        else:
            print "Loading custom settings from file 'settings'."
            data = self.restoreData()
            
            #Now set all the variables. See the defaults section above for what they are
            #TODO: This could be done in a 'cleaner way' by iterating over the variable dicionary returned by saveData()
            self.lowres.set(data[ 'lowres'])  
            self.highres.set(data[ 'highres'])  
            self.th_1.set(data[ 'th_1'])  
            self.th_2.set(data[ 'th_2']) 
            self.th_3.set(data[ 'th_3']) 
            self.stavar.set(data[ 'stavar'])  
            self.finvar.set(data[ 'finvar'])  
            self.bcorfac.set(data[ 'bcorfac'])  
            self.xfac.set(data[ 'xfac'])  
            self.plot_dat.set(data[ 'plot_dat'])  
            self.check_dat.set(data[ 'check_dat'])  
            self.check_dat2.set(data[ 'check_dat2'])  
            self.auto_read.set(data[ 'auto_read'])  
            self.datfillogi.set(data[ 'datfillogi'])  
            self.kde_bandwidth.set(data[ 'kde_bandwidth'])  
            self.kde_stop.set(data[ 'kde_stop'])  
            self.kde_points.set(data[ 'kde_points'])  
            self.sampsec.set(data[ 'sampsec'])  
            self.srange.set(data[ 'srange'])  
            self.sduration.set(data[ 'sduration'])  
            self.xmin_contour.set(data[ 'xmin_contour'])  
            self.xmax_contour.set(data[ 'xmax_contour'])  
            self.ymin_contour.set(data[ 'ymin_contour'])  
            self.ymax_contour.set(data[ 'ymax_contour'])  
            self.xbin_contour.set(data[ 'xbin_contour'])  
            self.ybin_contour.set(data[ 'ybin_contour'])  
            self.plat_max_grad.set(data[ 'plat_max_grad'])  
            self.background_tol.set(data[ 'background_tol'])  
            self.max_plat_cur.set(data[ 'max_plat_cur'])  
            self.max_points_plat.set(data[ 'max_points_plat'])  
            self.fractional_plat_tol.set(data[ 'fractional_plat_tol'])  
            self.min_points_plat.set(data[ 'min_points_plat'])  
            self.scan_u_th.set(data[ 'scan_u_th'])  
            self.scan_l_th.set(data[ 'scan_l_th'])  
            self.max_seg_len.set(data[ 'max_seg_len'])  
            self.min_seg_len.set(data[ 'min_seg_len'])  

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
        self.PlotMenu.menu.add_command(label='Logarithmic current histogram'
                , underline=0, command=self.kde_plot)
        self.PlotMenu.menu.add_command(label='Linear current histogram'
                , underline=1, command=self.linear_data_plot)
        self.PlotMenu.menu.add_command(label='Current-distance 2D histogram'
                , underline=1, command=contour_plot)
        self.PlotMenu['menu'] = self.PlotMenu.menu

        # Scan analysis menu -- everything related to reading in and filtering individual I(s) scans

        self.ScanAnalysis = Menubutton(self.mBar, text='Scan analysis',
                underline=0)
        self.ScanAnalysis.pack(side=LEFT, padx='2m')
        self.ScanAnalysis.menu = Menu(self.ScanAnalysis)
        self.ScanAnalysis.menu.add_command(label='Read scans to memory'
                , underline=0, background='grey', activebackground='red'
                , command=self.readfiles)
        #self.ScanAnalysis.menu.add('separator')
        self.ScanAnalysis.menu.add_command(label='Export current list to file'
                , underline=0, background='grey', activebackground='green'
                , command=export_current_data)
        self.ScanAnalysis.menu.add('separator')
        self.ScanAnalysis.menu.add_checkbutton(label='Read all scans in folder'
                , underline=0, variable=self.auto_read)
        self.ScanAnalysis.menu.add('separator')
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
        self.SettingsMenu.menu.add_command(label='2D histogram plotting'
                , underline=0, command=self.contour_params)
        self.SettingsMenu.menu.add_command(label='KDE plotting',
                underline=0, command=self.kde_params)
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


    def saveData(self):

        # Pickle a dictionary containing variables to be saved
        # See the default settings above for what they do

        data = {'lowres' : self.lowres.get(),
                'highres' : self.highres.get(),
                'th_1' : self.th_1.get(),
                'th_2' : self.th_2.get(),
                'th_3' : self.th_3.get(),
                'stavar' : self.stavar.get(),
                'finvar' : self.finvar.get(),
                'bcorfac' : self.bcorfac.get(),
                'xfac' : self.xfac.get(),
                'plot_dat' : self.plot_dat.get(),
                'check_dat' : self.check_dat.get(),
                'check_dat2' : self.check_dat2.get(),
                'auto_read' : self.auto_read.get(),
                'datfillogi' : self.datfillogi.get(),
                'kde_bandwidth' : self.kde_bandwidth.get(),
                'kde_stop' : self.kde_stop.get(),
                'kde_points' : self.kde_points.get(),
                'sampsec' : self.sampsec.get(),
                'srange' : self.srange.get(),
                'sduration' : self.sduration.get(),
                'xmin_contour' : self.xmin_contour.get(),
                'xmax_contour' : self.xmax_contour.get(),
                'ymin_contour' : self.ymin_contour.get(),
                'ymax_contour' : self.ymax_contour.get(),
                'xbin_contour' : self.xbin_contour.get(),
                'ybin_contour' : self.ybin_contour.get(),
                'plat_max_grad' : self.plat_max_grad.get(),
                'background_tol' : self.background_tol.get(),
                'max_plat_cur' : self.max_plat_cur.get(),
                'max_points_plat' : self.max_points_plat.get(),
                'fractional_plat_tol' : self.fractional_plat_tol.get(),
                'min_points_plat' : self.min_points_plat.get(),
                'scan_u_th' : self.scan_u_th.get(),
                'scan_l_th' : self.scan_l_th.get(),
                'max_seg_len' : self.max_seg_len.get(),
                'min_seg_len' : self.min_seg_len.get(),
                }

        file = open("settings", 'w')
        cPickle.dump(data, file)
        file.close()
        print 'Settings saved'

    def restoreData(self):

        # Return an unpickled dictionary with stored variables 
        file = open("settings", 'r')
        data = cPickle.load( file)
        file.close()
        return data
 
    def filename_browse(self):

        # Browse to select filename
        data.filename = askopenfilename(title='Select I(s) scan...',
                filetypes=[('Text files', '*.txt'), ('All files', '*')])

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
            FILE = open('nfo.txt', 'r')
            lines = FILE.readlines()
            number = 1.0
            for line in lines:
                self.text.insert(number, line)
                number += 1.0
            FILE.close()
        else:

            #Display an error if file isn't there
            self.error = showerror('Error', 'nfo.txt not found')

    def chopper_params(self): #MARKER

        # Configure chopper parameters
        self.chopper_params = Toplevel()
        self.chopper_params.title('Chopper parameters')

        # Put the parameters in a frame
        self.chopper_frame = Frame(self.chopper_params)
        self.chopper_frame.pack(side=TOP, padx=50, pady=5)

        # Draw the controls; see the text for what they do
        Label(self.chopper_frame, text='Minimum segment length (points)').grid(row=0,
                column=0)
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
                column=0)
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
                column=0)
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
                column=0)
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
              ).grid(row=0, column=0)

        Spinbox(
            self.resistance_frame,
            from_=0,
            to=100000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.lowres,
            ).grid(row=0, column=1)

        Label(self.resistance_frame, text='Feedback resistance (H):'
              ).grid(row=1, column=0)

        Spinbox(
            self.resistance_frame,
            from_=1000000,
            to=1000000000,
            increment=1,
            width=10,
            wrap=True,
            validate='all',
            textvariable=self.highres,
            ).grid(row=1, column=1)

        # Cutoff voltage settings for channel switchover
        Label(self.cutoff_frame, text='Cutoff voltages (V)'
              ).grid(row=0, column=0)

        Label(self.cutoff_frame, text='Channel 1:').grid(row=1,
                column=0)

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
                column=0)

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
                column=0)

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
        self.kde_params.title('KDE parameters')

        # Put the parameters in a frame
        self.kde_frame = Frame(self.kde_params)
        self.kde_frame.pack(side=TOP, pady=5)

        # KDE parameters
        Scale(
            self.kde_frame,
            label='KDE bandwidth:',
            variable=self.kde_bandwidth,
            from_=0,
            to=1,
            resolution=0.01,
            length=280,
            tickinterval=0,
            showvalue=YES,
            orient='horizontal',
            ).pack()

        Scale(
            self.kde_frame,
            label='KDE fit stop:',
            variable=self.kde_stop,
            from_=1,
            to=50,
            resolution=1,
            length=280,
            tickinterval=0,
            showvalue=YES,
            orient='horizontal',
            ).pack()

        Scale(
            self.kde_frame,
            label='KDE points:',
            variable=self.kde_points,
            from_=0,
            to=1000,
            resolution=50,
            length=280,
            tickinterval=0,
            showvalue=YES,
            orient='horizontal',
            ).pack()

    def data_filter(self):

        # Configure the data input parameters in a new window
        self.data_filter = Toplevel()
        self.data_filter.title('Data filtering')

        # Put the parameters in a frame
        self.data_frame = Frame(self.data_filter)
        self.data_frame.pack(side=TOP, padx=50, pady=5)

        # Draw the controls; see the text for what they do
        Label(self.data_frame, text='Cumulative log(I) threshold:'
              ).grid(row=0, column=0)
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
              ).grid(row=0, column=0)
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
              ).grid(row=1, column=0)
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
              ).grid(row=2, column=0)
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
              ).grid(row=3, column=0)
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

    def adc_params(self):

        # Configure the ADC / STM parameters
        self.adc_params = Toplevel()
        self.adc_params.title('ADC / STM parameters')

        # Put the parameters in a frame
        self.adc_frame = Frame(self.adc_params)
        self.adc_frame.pack(side=TOP, padx=50, pady=5)

        # Draw the controls; see the text for what they do
        Label(self.adc_frame, text='Samples per second:').grid(row=0,
                column=0)
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
                column=0)
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
                column=0)
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
                column=0)
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
                column=0)
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
                column=0)
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
                column=0)
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
                column=0)
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
                column=0)
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
                column=0)
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
                column=0)
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
              ).grid(row=2, column=0)
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
                column=0)
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
                column=0)
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
              ).grid(row=5, column=0)
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

    def kde_plot(self):

        # Fetch KDE parameters from control widgets
        kde_bandwidth = self.kde_bandwidth.get()
        kde_start = kde_bandwidth  # The lowest the KDE start can go
        kde_stop = self.kde_stop.get()
        kde_points = self.kde_points.get()

        # Plot KDE function and read all current data to a single list
        # Initialse current list and append currents from individual files
        i_list = []
        for reading in data.i_dat_all_combined:
            for value in reading:

            # i_list.append(np.log10(float(value)))
            # i_list.append(value)

                if value > 0:
                    i_list.append(np.log10(value))
        print 'Lists appended, found:', len(i_list), 'data points.'

        if (len(i_list) < 1):
            self.error = showerror('Error', 'No data in memory')
            return

        # Fit KDE and setup plot
        # x = np.linspace(min(i_list),max(i_list),kde_points)
        # haa = statistics.bandwidth(i_list, kernel='Epanechnikov')
        # bandwidth = haa /2
        # z = statistics.pdf(i_list, x, h=bandwidth, kernel='E')
        # print "Optimal bandwidth:", haa
        # Plot KDE and histogram
        # plot_y_lim = ( max(z) ) * 1.1
        # plot_y_lim = 0.12
        # plt.xscale('log')
        plt.hist(i_list, bins=500, facecolor='black', normed=1)
        # qc= 7.74
        # oc= 7.08....
        # plt.axvline(np.log10(qc), def __init__(self, myParent):
        # plt.axvline(oc, c='r')
        # plt.axvline(oc*2, c='r')
        # plt.axvline(oc*3, c='r')
        # plt.axvline(oc*4, c='r')
        # plt.axvline(oc*5, c='r')
        # plt.plot(x, z,'b', label='KDE', linewidth=3)
        # plt.axis([0, kde_stop, 0, plot_y_lim])
        # plt.legend()
        plt.grid()
        plt.ylabel('Density', fontsize=14)
        plt.xlabel('log_10[current (nanoamps)]', fontsize=14)
        # plt.title('Molecular conductance: 0 degrees', fontsize=16)
        plt.show()

    def linear_data_plot(self):

        # Fetch KDE parameters from control widgets
        i_list = []
        for reading in data.i_dat_all_combined:
            for value in reading:
                # i_list.append(np.log10(float(value)))
                # i_list.append(value)
                if ( (value > 0) and (value < 0.5) ):
                    i_list.append(value)
        
        if (len(i_list) < 1):
            self.error = showerror('Error', 'No data in memory')
            return

        print 'Lists appended, found:', len(i_list), 'data points.'
        plt.hist(i_list, bins=2000, facecolor='black', normed=1)
        plt.grid()
        plt.ylabel('Density', fontsize=14)
        plt.xlabel('current (nanoamps)', fontsize=14)
        plt.show()

    def readfiles(self):

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

            # Save the current working directory
            cwd_bak = os.getcwd()

            # Have a poke around the user selected data directory for the number of files
            # which must be sequentially numbered
            os.chdir(target_folder)
            filenumber_list = []
            for (path, dirs, files) in os.walk(os.getcwd()):
                for filename in files:
                    filenumber_list.append(int(filename[-8:-4]))
            start = min(filenumber_list)
            finish = max(filenumber_list)

            # Go back to the current working directory
            os.chdir(cwd_bak)
        else:
            # Use the default or user set input range
            start = int(self.stavar.get())
            finish = int(self.finvar.get())

        bcorfac = self.bcorfac.get()
        dat_input(start, finish, bcorfac)
        print 'End of file input'

    def importdata(self):
        
        # Needs to call tea break maker
        tea_break_maker()


root = Tk()
root.title('Scanmaster 3000 v0.42')

egraph = egraph(root)
controller = controller(root)
root.mainloop()

