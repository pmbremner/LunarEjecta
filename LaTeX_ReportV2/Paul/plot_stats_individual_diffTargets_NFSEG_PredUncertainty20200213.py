#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#       Process individual difference targets and determine the median
#       and quartiles Q1 and Q3. Plot the results, along with the Mean
#       on Cumulative Probability and Normal-Score plots.
#
#       Run this script in the directory where the files for the difference 
#       target files are located.
#
#
#       Written 20200213. PMBremner
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


import os
import shutil
#import pandas as pd
import numpy as np
#from numpy.linalg import inv
from scipy import stats

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#import collections  # to use OrderDict
#from copy import deepcopy as dp

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#       FUNCTIONS
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Calculate the MEDIAN and Q1,Q3 quartiles of an 1D array
def calc_median_Q1_Q3(array):
    
    # Sort the array
    array = np.sort(array,kind='mergesort')
    
    # Determine the size of the array only once
    n = tdiff.size
    
    # Check where there are an even or odd number of elements -- 0=even, 1=odd
    evenodd = n % 2
    
    if (evenodd==0):
        # Average the two elements in positions n/2 and (n+2)/2
        n1 = np.int(n/2)
        n2 = np.int((n+2)/2)
        median = 0.5*(array[n1] + array[n2])
        
    elif (evenodd==1):
        # Report the element in position (n+1)/2
        n1 = (n+1)/2
        median = array[n1]
    else:
        raise Exception("The evenness of the array len could not be determined!")
    # END IF
    
    # Find the first and third Quartiles
    intcheck = n % 4
    if (intcheck==0): # Integer found, take mean of k and k+1 positions
        position = np.int(n/4)
        Q1 = 0.5*(array[position]+array[position+1])
        
        position = np.int(n*0.75)
        Q3 = 0.5*(array[position]+array[position+1])
        
    else: # No integer, round up use that position
        position = np.int(np.ceil(n/4))
        Q1 = array[position]
        
        position = np.int(np.ceil(n*0.75))
        Q3 = array[position]
        
    # END IF
    
    
    return median,Q1,Q3



# Read in the Mean and Standard Deviation values in the stats results file
# This function expects the following format:
#       0. TargetName_SP1
#       1. TargetName_SP2
#       2. Group_SP1
#       3. Group_SP2
#       4. Mean_SP1
#       5. Stdv_SP1
#       6. NumSamp_SP1
#       7. Mean_SP2
#       8. Stdv_SP2
#       9. NumSamp_SP2
#       10. Mean_SP2-1
#       11. Stdv_SP2-1
#       12. NumSamp_SP2-1
#       13. Weight_SP1
#       14. Weight_SP2
#       15. ObsClass
#       16. ObsClass_Num
#       17. ObsValue
def read_stats_file(statsfile,meancol,stddvcol):
    precalc_stats = []
    
    with open(statsfile,'r')as sf:
        for i,sline in enumerate(sf):
            # initialize the tmp array
            tmp = 3*[None]
            
            # Skip the header line
            if (i>0):
                stmp = sline.strip().split()
                tmp[0]=('diff_' +
                        ".".join([stmp[1],stmp[3]]) +
                        "-" +
                        ".".join([stmp[0],stmp[2]]) +
                        '.dat')
                tmp[1] = stmp[meancol]
                tmp[2] = stmp[stddvcol]
                
                precalc_stats.append(tmp)
            #
        #
    #
    
    # Return array in the form [0]=created_diff_filename,[1]=meancol,[2]=stddvcol
    return precalc_stats
#

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#       MAIN PROGRAM
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Specify the directory of interest
#DIR = /beohome/a/MultiBasin/NFSEG/MODFLOW/NFSEGv1_1/20191028_predictive_uncertainty/Data/targetfiles_diff
#DIR='../../Data/targetfiles_diff/'
DIR='./'


# Name the target list file
targetlist=(DIR+'target_list_NFRWSPA.txt')


# Read in the pre-calculated Mean and Standard Deviations from stats file
# The Read function expects the following format:
#       0. TargetName_SP1
#       1. TargetName_SP2
#       2. Group_SP1
#       3. Group_SP2
#       4. Mean_SP1
#       5. Stdv_SP1
#       6. NumSamp_SP1
#       7. Mean_SP2
#       8. Stdv_SP2
#       9. NumSamp_SP2
#       10. Mean_SP2-1
#       11. Stdv_SP2-1
#       12. NumSamp_SP2-1
#       13. Weight_SP1
#       14. Weight_SP2
#       15. ObsClass
#       16. ObsClass_Num
#       17. ObsValue
#statsfile = '../../Data/stats_table_case_007h_pred_2009_2035_NFPA_calib_only.dat'
statsfile = '../stats_table_case_007h_pred_2009_2035_NFPA_calib_only.dat'
precalc_stats = read_stats_file(statsfile,10,11)


# Name the output summary file
# Open and write out the header information
#outputlist=('../../Data/stats_table_case_007h_pred_2009_2035_NFPA_calib-and-weighted_only_mean-stddv-med-q1-q3.dat')
outputlist=('../stats_table_case_007h_pred_2009_2035_NFPA_calib-and-weighted_only_mean-stddv-med-q1-q3.dat')
outf = open(outputlist,'w')
outf.write ('DO NOT EDIT. AUTO-GENERATED FILE from plot_stats_individual_diffTargets_NFSEG_PredUncertainty20200213.py\n')
outf.write ('Filename ObsClass mean stddev median Q1 Q3\n')


# Check if the Normal Score figure directory exists
# If True: remove/remake it to start fresh
# If False: make it
FIGDIR1 = ('../Figures/tdiff_normscore')
if (os.path.exists(FIGDIR1)):
    #os.rmdir(FIGDIR1)
    shutil.rmtree(FIGDIR1, ignore_errors=True)
    os.mkdir(FIGDIR1)
else:
    os.mkdir(FIGDIR1)
#

# Check if the Cumulative Frequency figure directory exists
# If True: remove/remake it to start fresh
# If False: make it
FIGDIR2 = ('../Figures/tdiff_cmfreq')
if (os.path.exists(FIGDIR2)):
    #os.rmdir(FIGDIR2)
    shutil.rmtree(FIGDIR2, ignore_errors=True)
    os.mkdir(FIGDIR2)
else:
    os.mkdir(FIGDIR2)
#


# Define plot properties
# -----------------------------------------------------------------------------
linewidth_means = 2.5 #4.0 #2.5
linewidth_medians = 1.5 #2.5 #4.0
linestyle_mean = '-'
linestyle_stdplus = '-.'
linestyle_stdneg = ':'
linestyle_median = '--'
linestyle_Q1 = ':'
linestyle_Q3 = '-.'
#
color_sample = 'tab:blue'
color_ref = 'tab:orange'
color_mean = 'tab:red'
color_stdplus = 'tab:pink'
color_stdneg = 'tab:pink'
color_median = 'k'
color_Q1 = 'tab:gray'
color_Q3 = 'tab:gray'
#
label_sample = 'Target Diff'
label_ref = 'Ref Norm Distribution'
label_mean = 'Mean'
label_stdplus = 'Std dev +'
label_stdneg = 'Std dev -'
label_median = 'Median'
label_Q1 = 'Q1'
label_Q3 = 'Q3'
#
textfontsize = 12
legendfontsize = 11.5
textstyle = 'normal' #'bold'
labelfontsize = 18
ticklabelsize = 14
#
# Set whether to actually make the plots, or just make the output data file
plotflag = True
#
# Set whetherto simply show the plots (False), or save them (True)
saveplotflag = True
# -----------------------------------------------------------------------------



# March through all the targets in the list
with open(targetlist,'r') as tarf:
    tarinfo = []
    
    for j,tline in enumerate(tarf):
        
        # Skip 2 header lines
        if (j>1):
            tline = tline.rstrip()
            Filename, ObsClass, NonZeroWeight_SP1, NonZeroWeight_SP2 = tline.split()
            
            # Check for Zero Weighted targets
            if (NonZeroWeight_SP1=='y' and NonZeroWeight_SP2=='y'):
                #tarinfo.append([Filename, ObsClass, NonZeroWeight_SP1, NonZeroWeight_SP2])
                print (Filename)
                
                # Find the current target in the list of precalc_stats
                # precalc_stats FORMAT: filename, mean, standard deviation
                curstat = []
                for stat in precalc_stats:
                    if (stat[0]==Filename):
                        curstat = stat
                        break
                    # END IF
                # END FOR LOOP OVER stat
                
                # Capture the two stats needed
                curmean = curstat[1]
                curstddv = curstat[2]
                
                #print (curstat,tarinfo[0][0])
                
                # Read the current target's values from all the runs
                # NOTE: Only use the Calibrated Models
                with open((DIR+Filename),'r') as fdiff:
                    
                    tdiff = []
                    
                    for i,dline in enumerate(fdiff):
                        # Skip the header lines
                        if (i>5):
                            #0 RunNum
                            #1 ResFileName
                            #2 ModeledVal_SP1
                            #3 ModeledVal_SP2
                            #4 ModDIFF_SP2-SP1
                            #5 CalibrationStatus_SP1
                            #6 CalibrationStatus_SP2
                            
                            if (dline.strip().split()[5]=='Y' and dline.strip().split()[6]=='Y'):
                                tdiff.append(dline.strip().split()[4])
                            #
                    #  END FOR LOOP dline
                #  END WITH -- CLOSE FILE fdiff
                
                # Make tdiff a numpy array
                tdiff = np.asarray(tdiff,dtype=np.float)
                
                
                # Calculate the median
                # Calculate the first, second, and 3rd quartiles, Q1,Q2,Q3
                median,Q1,Q3 = calc_median_Q1_Q3(tdiff)
                
                # Put all relavent values together in one array and
                # write line to file
                tarinfo.append([Filename, ObsClass, curmean, curstddv, median, Q1, Q3])
                outf.write (Filename + ' ' +
                            ObsClass + ' ' +
                            str(curmean) + ' ' +
                            str(curstddv) + ' ' +
                            str(median) + ' ' +
                            str(Q1) + ' ' +
                            str(Q3) + '\n')
                
                
                # Create the the float versions of the calculated statistics
                mean_val = np.float(curmean)
                stddev_val = np.float(curstddv)
                sdplus_val = mean_val + stddev_val
                sdminus_val = mean_val - stddev_val
                median_val = np.float(median)
                Q1_val = np.float(Q1)
                Q3_val = np.float(Q3)
                
                
                # Plot the Cumulative Probability with the Stats and
                # Plot the tdiff values as function of Normal Score
                # *************************************************************
                
                # Get the ordered list of plotting positions for the Residual data
                pos = np.zeros(tdiff.size,dtype=np.float)
                tdiff_ordered = np.sort(tdiff,kind='mergesort')
                for ip in range(pos.size): pos[ip] = ((ip+1)-0.4) / (pos.size+0.2)
                
                # Obtain Zscores pretending that the distribution is Normal (used for comparison)
                pos = stats.zscore(pos)
                
                # Construct the reference lines for the Normal Score Plots
                # The lines will show the following properties:
                #       1. mean
                #       2. plus and minus standard deviations
                #       3. median
                #       4. Q1, Q3 quartiles
                # -----------------------------------------------
                # Construct the reference line
                # The line will be y = residuals_mean + (residuals_stddv * x)
                # Find only the extreme two points to be able to draw the line
                refy_min = mean_val + (stddev_val * pos.min())
                refy_max = mean_val + (stddev_val * pos.max())
                refline_x = np.array([pos.min(),pos.max()],dtype=np.float)
                refline_y = np.array([refy_min,refy_max],dtype=np.float)
                
                norm_line_y_1 = np.array([tdiff_ordered.min(), tdiff_ordered.max()],dtype=np.float)
                reduct = (tdiff_ordered.max()-tdiff_ordered.min())*0.25 # 50% total reduction, an quarter from each end
                norm_line_y_2 = np.array([(tdiff_ordered.min()+reduct), (tdiff_ordered.max()-reduct)],dtype=np.float)
                
                norm_y_min = min(refy_min,tdiff_ordered.min())
                norm_y_max = max(refy_max,tdiff_ordered.max())
                norm_y_range = norm_y_max - norm_y_min
                norm_y_min = norm_y_min - norm_y_range*0.025
                norm_y_max = norm_y_max + norm_y_range*0.025
                
                idx = (np.abs(tdiff_ordered - mean_val)).argmin()
                norm_mean_x = pos[idx]
                norm_mean_line_x = np.array([norm_mean_x,norm_mean_x],dtype=np.float)
                
                idx = (np.abs(tdiff_ordered - sdplus_val)).argmin()
                norm_sdplus_x = pos[idx]
                norm_sdplus_line_x = np.array([norm_sdplus_x,norm_sdplus_x],dtype=np.float)
                
                idx = (np.abs(tdiff_ordered - sdminus_val)).argmin()
                norm_sdminus_x = pos[idx]
                norm_sdminus_line_x = np.array([norm_sdminus_x,norm_sdminus_x],dtype=np.float)
                
                idx = (np.abs(tdiff_ordered - median_val)).argmin()
                norm_median_x = pos[idx]
                norm_median_line_x = np.array([norm_median_x,norm_median_x],dtype=np.float)
                
                idx = (np.abs(tdiff_ordered - Q1_val)).argmin()
                norm_Q1_x = pos[idx]
                norm_Q1_line_x = np.array([norm_Q1_x,norm_Q1_x],dtype=np.float)
                
                idx = (np.abs(tdiff_ordered - Q3_val)).argmin()
                norm_Q3_x = pos[idx]
                norm_Q3_line_x = np.array([norm_Q3_x,norm_Q3_x],dtype=np.float)
                # -----------------------------------------------
                
                
                # Obtain the cumulative frequency of the tdiff values
                numbins = np.int(200)
                cmfreq = stats.cumfreq(tdiff_ordered, numbins=numbins)
                cmfreq_x = cmfreq.lowerlimit + np.linspace(0,cmfreq.binsize*cmfreq.cumcount.size,cmfreq.cumcount.size)
                
                # Create a random Normal Distribution as reference
                # np.random.normal (mean, standard_deviation, array_size)
                samples = np.random.normal (mean_val,stddev_val,tdiff.size)
                cmfreq_norm = stats.cumfreq(samples, numbins=numbins)
                cmfreq_norm_x = cmfreq_norm.lowerlimit + np.linspace(0,cmfreq_norm.binsize*cmfreq_norm.cumcount.size,cmfreq_norm.cumcount.size)
                
                
                cmfreq_x_min = min(cmfreq_x.min(),cmfreq_norm_x.min(),tdiff_ordered.min())
                cmfreq_x_max = max(cmfreq_x.max(),cmfreq_norm_x.max(),tdiff_ordered.max())
                cmfreq_x_range = cmfreq_x_max - cmfreq_x_min
                cmfreq_x_min = cmfreq_x_min - cmfreq_x_range*0.025
                cmfreq_x_max = cmfreq_x_max + cmfreq_x_range*0.025
                
                
                # Construct the reference lines for the Cumulative Freq Plots
                # The lines will show the following properties:
                #       1. mean
                #       2. plus and minus standard deviations
                #       3. median
                #       4. Q1, Q3 quartiles
                # -----------------------------------------------
                # Make an array for the x-coordinates of the lines
                # The same x-coordinates are applied to all lines
                cmfreq_line_x_1 = np.array([tdiff_ordered.min(), tdiff_ordered.max()],dtype=np.float)
                reduct = (tdiff_ordered.max()-tdiff_ordered.min())*0.25 # 50% total reduction, an quarter from each end
                cmfreq_line_x_2 = np.array([(tdiff_ordered.min()+reduct), (tdiff_ordered.max()-reduct)],dtype=np.float)
                
                idx = (np.abs(cmfreq_x - mean_val)).argmin()
                cmfreq_mean_y = cmfreq.cumcount[idx]
                cmfreq_mean_line_y = np.array([cmfreq_mean_y,cmfreq_mean_y],dtype=np.float)
                
                idx = (np.abs(cmfreq_x - sdplus_val)).argmin()
                cmfreq_sdplus_y = cmfreq.cumcount[idx]
                cmfreq_sdplus_line_y = np.array([cmfreq_sdplus_y,cmfreq_sdplus_y],dtype=np.float)
                
                idx = (np.abs(cmfreq_x - sdminus_val)).argmin()
                cmfreq_sdminus_y = cmfreq.cumcount[idx]
                cmfreq_sdminus_line_y = np.array([cmfreq_sdminus_y,cmfreq_sdminus_y],dtype=np.float)
                
                idx = (np.abs(cmfreq_x - median_val)).argmin()
                cmfreq_median_y = cmfreq.cumcount[idx]
                cmfreq_median_line_y = np.array([cmfreq_median_y,cmfreq_median_y],dtype=np.float)
                
                idx = (np.abs(cmfreq_x - Q1_val)).argmin()
                cmfreq_Q1_y = cmfreq.cumcount[idx]
                cmfreq_Q1_line_y = np.array([cmfreq_Q1_y,cmfreq_Q1_y],dtype=np.float)
                
                idx = (np.abs(cmfreq_x - Q3_val)).argmin()
                cmfreq_Q3_y = cmfreq.cumcount[idx]
                cmfreq_Q3_line_y = np.array([cmfreq_Q3_y,cmfreq_Q3_y],dtype=np.float)
                # -----------------------------------------------
                
                # -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
                if plotflag:
                    
                                    
                    # Make and Save Plots
                    # -----------------------------------------------
                    fig1 = plt.figure(figsize=(10,6))
                    # plot = fig.add_axes([left corner,bottom corner,width,height], label='name')
                    normplot = fig1.add_axes([0.0,0.0,0.70,1.0], label='axes1')
                    fig2 = plt.figure(figsize=(10,6))
                    cmfreqplot = fig2.add_axes([0.0,0.0,0.70,1.0], label='axes3')
                    
                    # Construct the plot label
                    part1='.'.join(Filename.split('_')[1:])
                    part2='.'.join(part1.split('.')[:-1])
                    TName1 = part2.split('-')[0]
                    TName2 = part2.split('-')[1]
                    #TName1 = 'abcdefghijklmnopqrstuvwxyzabcdef'
                    #TName1 = 'qr01_crystal_sprgrp.qs_spring01.NonZeroWeight.dat'
                    labelstr = (TName1 + '\n' + TName2 + '\n' + ObsClass)
                    
                    # Plot the data and the reference norm distribution
                    normplot.scatter(pos, tdiff_ordered, marker='D', s=35, label=label_sample
                                     , color=color_sample, edgecolor=color_sample, linewidths=1.00, zorder=3)
                    normplot.plot(refline_x, refline_y, label=label_ref
                                  , color=color_ref, linewidth=3.50, zorder=4)
                    
                    
                    # Place the calculated statistics on the plot
                    zorder = 2 #5
#                    normplot.plot(norm_mean_line_x, norm_line_y_1
#                                  , label=label_mean, color=color_mean, zorder=zorder
#                                  , linestyle=linestyle_mean, linewidth=linewidth_means)
                    normplot.plot([pos.min(),norm_mean_x], [mean_val,mean_val]
                                  , label=label_mean, color=color_mean, zorder=zorder
                                  , linestyle=linestyle_mean, linewidth=linewidth_means)
                    normplot.plot([norm_mean_x,norm_mean_x], [tdiff_ordered.min(),mean_val]
                                  , color=color_mean, zorder=zorder
                                  , linestyle=linestyle_mean, linewidth=linewidth_means)
                                        
#                    normplot.plot(norm_sdplus_line_x, norm_line_y_1
#                                  , label=label_stdplus, color=color_stdplus, zorder=zorder
#                                  , linestyle=linestyle_stdplus, linewidth=linewidth_means)
                    normplot.plot([pos.min(),norm_sdplus_x], [sdplus_val,sdplus_val]
                                  , label=label_stdplus, color=color_stdplus, zorder=zorder
                                  , linestyle=linestyle_stdplus, linewidth=linewidth_means)
                    normplot.plot([norm_sdplus_x,norm_sdplus_x], [tdiff_ordered.min(),sdplus_val]
                                  , color=color_stdplus, zorder=zorder
                                  , linestyle=linestyle_stdplus, linewidth=linewidth_means)
                                        
#                    normplot.plot(norm_sdminus_line_x, norm_line_y_1
#                                  , label=label_stdneg, color=color_stdneg, zorder=zorder
#                                  , linestyle=linestyle_stdneg, linewidth=linewidth_means)
                    normplot.plot([pos.min(),norm_sdminus_x], [sdminus_val,sdminus_val]
                                  , label=label_stdneg, color=color_stdneg, zorder=zorder
                                  , linestyle=linestyle_stdneg, linewidth=linewidth_means)
                    normplot.plot([norm_sdminus_x,norm_sdminus_x], [tdiff_ordered.min(),sdminus_val]
                                  , color=color_stdneg, zorder=zorder
                                  , linestyle=linestyle_stdneg, linewidth=linewidth_means)
                                        
#                    normplot.plot(norm_median_line_x, norm_line_y_2
#                                  , label=label_median, color=color_median, zorder=zorder
#                                  , linestyle=linestyle_median, linewidth=linewidth_medians)
                    normplot.plot([pos.min(),norm_median_x], [median_val,median_val]
                                  , label=label_median, color=color_median, zorder=zorder
                                  , linestyle=linestyle_median, linewidth=linewidth_medians)
                    normplot.plot([norm_median_x,norm_median_x], [tdiff_ordered.min(),median_val]
                                  , color=color_median, zorder=zorder
                                  , linestyle=linestyle_median, linewidth=linewidth_medians)
                                        
#                    normplot.plot(norm_Q3_line_x, norm_line_y_2
#                                  , label=label_Q3, color=color_Q3, zorder=zorder
#                                  , linestyle=linestyle_Q3, linewidth=linewidth_medians)
                    normplot.plot([pos.min(),norm_Q3_x], [Q3_val,Q3_val]
                                  , label=label_Q3, color=color_Q3, zorder=zorder
                                  , linestyle=linestyle_Q3, linewidth=linewidth_medians)
                    normplot.plot([norm_Q3_x,norm_Q3_x], [tdiff_ordered.min(),Q3_val]
                                  , color=color_Q3, zorder=zorder
                                  , linestyle=linestyle_Q3, linewidth=linewidth_medians)
                                        
#                    normplot.plot(norm_Q1_line_x, norm_line_y_2
#                                  , label=label_Q1, color=color_Q1, zorder=zorder
#                                  , linestyle=linestyle_Q1, linewidth=linewidth_medians)
                    normplot.plot([pos.min(),norm_Q1_x], [Q1_val,Q1_val]
                                  , label=label_Q1, color=color_Q1, zorder=zorder
                                  , linestyle=linestyle_Q1, linewidth=linewidth_medians)
                    normplot.plot([norm_Q1_x,norm_Q1_x], [tdiff_ordered.min(),Q1_val]
                                  , color=color_Q1, zorder=zorder
                                  , linestyle=linestyle_Q1, linewidth=linewidth_medians)
                    
                    
                    zorder = 5; s1 = 110; s2 = 160
                    normplot.scatter(norm_mean_x, mean_val
                                  , marker='o', color=color_mean, zorder=zorder
                                  , s=s2, linewidth=linewidth_means)
                    normplot.scatter(norm_sdplus_x, sdplus_val
                                  , marker='o', color=color_stdplus, zorder=zorder
                                  , s=s2, linewidth=linewidth_means)
                    normplot.scatter(norm_sdminus_x, sdminus_val
                                  , marker='o', color=color_stdneg, zorder=zorder
                                  , s=s2, linewidth=linewidth_means)
                    normplot.scatter(norm_median_x, median_val
                                  , marker='o', color=color_median, zorder=zorder
                                  , s=s1, linewidth=linewidth_medians)
                    normplot.scatter(norm_Q3_x, Q3_val
                                  , marker='o', color=color_Q3, zorder=zorder
                                  , s=s1, linewidth=linewidth_medians)
                    normplot.scatter(norm_Q1_x, Q1_val
                                  , marker='o', color=color_Q1, zorder=zorder
                                  , s=s1, linewidth=linewidth_medians)
                    
                    
                    # bbox_to_anchor = (x, y, width, height)
                    legend1 = normplot.legend(fontsize=legendfontsize, loc='upper left'
                                              , bbox_to_anchor=(1.01,1.0), borderaxespad=0.0
                                              , handlelength=3.0, framealpha=1.0)
                    
                    # Label the plot
                    text_y = 0.0 #0.54 #0.97 # in axis coordinates
                    text_x = 1.02 #0.05 # in axis coordinates
                    normplot.text(text_x,text_y, labelstr, fontsize=textfontsize, fontweight=textstyle
                                  , color='k', verticalalignment='bottom', horizontalalignment='left'
                                  , rotation=90, transform=normplot.transAxes
                                  , bbox=dict(edgecolor='none',facecolor='w',alpha=1.0), zorder=6)
                    
                    
                    # Adjust the look and size of the axes
                    normplot.set_ylim([norm_y_min,norm_y_max])
                    normplot.set_ylabel('Delta T', fontsize=labelfontsize, fontweight='bold')
                    normplot.set_xlabel('Normal score', fontsize=labelfontsize, fontweight='bold')
                    normplot.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
                    normplot.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
                    normplot.tick_params(axis='both', direction='out', left='on'
                                         , labelleft='on', labelsize=ticklabelsize
                                         , which='major', length=10, width=1.5
                                         , gridOn='on')
                    normplot.tick_params(axis='both', direction='out', left='on'
                                         , labelleft='on', labelsize=ticklabelsize
                                         , which='minor', length=6, width=1.5)
                    
                    
                    # Plot the data and the reference norm distribution
                    #cmfreqplot.bar(cmfreq_x,cmfreq.cumcount, width=cmfreq.binsize)
                    cmfreqplot.scatter(cmfreq_x, cmfreq.cumcount
                                       , marker='D', s=40, label=label_sample
                                       , color=color_sample, edgecolor=color_sample
                                       , linewidths=0.25, zorder=3)
                    cmfreqplot.plot(cmfreq_norm_x, cmfreq_norm.cumcount
                                    , label=label_ref, color=color_ref, linewidth=3.50, zorder=4)
                    
                    
                    # Place the calculated statistics on the plot
                    zorder = 2 #5
#                    cmfreqplot.plot(cmfreq_line_x_1, cmfreq_mean_line_y
#                                    , label=label_mean, color=color_mean, zorder=zorder
#                                    , linestyle=linestyle_mean, linewidth=linewidth_means)
                    cmfreqplot.plot([tdiff_ordered.min(),mean_val], [cmfreq_mean_y,cmfreq_mean_y]
                                    , label=label_mean, color=color_mean, zorder=zorder
                                    , linestyle=linestyle_mean, linewidth=linewidth_means)
                    cmfreqplot.plot([mean_val,mean_val], [0,cmfreq_mean_y]
                                    , color=color_mean, zorder=zorder
                                    , linestyle=linestyle_mean, linewidth=linewidth_means)
                    
#                    cmfreqplot.plot(cmfreq_line_x_1, cmfreq_sdplus_line_y
#                                    , label=label_stdplus, color=color_stdplus, zorder=zorder
#                                    , linestyle=linestyle_stdplus, linewidth=linewidth_means)
                    cmfreqplot.plot([tdiff_ordered.min(),sdplus_val], [cmfreq_sdplus_y,cmfreq_sdplus_y]
                                    , label=label_stdplus, color=color_stdplus, zorder=zorder
                                    , linestyle=linestyle_stdplus, linewidth=linewidth_means)
                    cmfreqplot.plot([sdplus_val,sdplus_val], [0,cmfreq_sdplus_y]
                                    , color=color_stdplus, zorder=zorder
                                    , linestyle=linestyle_stdplus, linewidth=linewidth_means)
                                        
#                    cmfreqplot.plot(cmfreq_line_x_1, cmfreq_sdminus_line_y
#                                    , label=label_stdneg, color=color_stdneg, zorder=zorder
#                                    , linestyle=linestyle_stdneg, linewidth=linewidth_means)
                    cmfreqplot.plot([tdiff_ordered.min(),sdminus_val], [cmfreq_sdminus_y,cmfreq_sdminus_y]
                                    , label=label_stdneg, color=color_stdneg, zorder=zorder
                                    , linestyle=linestyle_stdneg, linewidth=linewidth_means)
                    cmfreqplot.plot([sdminus_val,sdminus_val], [0,cmfreq_sdminus_y]
                                    , color=color_stdneg, zorder=zorder
                                    , linestyle=linestyle_stdneg, linewidth=linewidth_means)
                                        
#                    cmfreqplot.plot(cmfreq_line_x_2, cmfreq_median_line_y
#                                    , label=label_median, color=color_median, zorder=zorder
#                                    , linestyle=linestyle_median, linewidth=linewidth_medians)
                    cmfreqplot.plot([tdiff_ordered.min(),median_val], [cmfreq_median_y,cmfreq_median_y]
                                    , label=label_median, color=color_median, zorder=zorder
                                    , linestyle=linestyle_median, linewidth=linewidth_medians)
                    cmfreqplot.plot([median_val,median_val], [0,cmfreq_median_y]
                                    , color=color_median, zorder=zorder
                                    , linestyle=linestyle_median, linewidth=linewidth_medians)
                                        
#                    cmfreqplot.plot(cmfreq_line_x_2, cmfreq_Q1_line_y
#                                    , label=label_Q1, color=color_Q1, zorder=zorder
#                                    , linestyle=linestyle_Q1, linewidth=linewidth_medians)
                    cmfreqplot.plot([tdiff_ordered.min(),Q1_val], [cmfreq_Q1_y,cmfreq_Q1_y]
                                    , label=label_Q1, color=color_Q1, zorder=zorder
                                    , linestyle=linestyle_Q1, linewidth=linewidth_medians)
                    cmfreqplot.plot([Q1_val,Q1_val], [0,cmfreq_Q1_y]
                                    , color=color_Q1, zorder=zorder
                                    , linestyle=linestyle_Q1, linewidth=linewidth_medians)
                                        
#                    cmfreqplot.plot(cmfreq_line_x_2, cmfreq_Q3_line_y
#                                    , label=label_Q3, color=color_Q3, zorder=zorder
#                                    , linestyle=linestyle_Q3, linewidth=linewidth_medians)
                    cmfreqplot.plot([tdiff_ordered.min(),Q3_val], [cmfreq_Q3_y,cmfreq_Q3_y]
                                    , label=label_Q3, color=color_Q3, zorder=zorder
                                    , linestyle=linestyle_Q3, linewidth=linewidth_medians)
                    cmfreqplot.plot([Q3_val,Q3_val], [0,cmfreq_Q3_y]
                                    , color=color_Q3, zorder=zorder
                                    , linestyle=linestyle_Q3, linewidth=linewidth_medians)
                    
                    zorder = 5; s1 = 110; s2 = 160
                    cmfreqplot.scatter(mean_val, cmfreq_mean_y
                                    , marker='o', color=color_mean, zorder=zorder
                                    , s=s2, linewidth=linewidth_means)
                    cmfreqplot.scatter(sdplus_val, cmfreq_sdplus_y
                                    , marker='o', color=color_stdplus, zorder=zorder
                                    , s=s2, linewidth=linewidth_means)
                    cmfreqplot.scatter(sdminus_val, cmfreq_sdminus_y
                                    , marker='o', color=color_stdneg, zorder=zorder
                                    , s=s2, linewidth=linewidth_means)
                    cmfreqplot.scatter(median_val, cmfreq_median_y
                                    , marker='o', color=color_median, zorder=zorder
                                    , s=s1, linewidth=linewidth_medians)
                    cmfreqplot.scatter(Q1_val, cmfreq_Q1_y
                                    , marker='o', color=color_Q1, zorder=zorder
                                    , s=s1, linewidth=linewidth_medians)
                    cmfreqplot.scatter(Q3_val, cmfreq_Q3_y
                                    , marker='o', color=color_Q3, zorder=zorder
                                    , s=s1, linewidth=linewidth_medians)
                    
                    # bbox_to_anchor = (x, y, width, height)
                    legend2 = cmfreqplot.legend(fontsize=legendfontsize, loc='upper left'
                                                , bbox_to_anchor=(1.01,1.0), borderaxespad=0.0
                                                , handlelength=3.0, framealpha=1.0)
                    
                    
                    # Label the plot
                    text_y = 0.0 #0.98 # in axis coordinates
                    text_x = 1.02 #0.02 # in axis coordinates
                    cmfreqplot.text(text_x,text_y, labelstr, fontsize=textfontsize, fontweight=textstyle
                                    , color='k', verticalalignment='bottom', horizontalalignment='left'
                                    , rotation=90, transform=cmfreqplot.transAxes
                                    , bbox=dict(edgecolor='none', facecolor='w', alpha=1.0), zorder=6)
                    
                    
                    # Adjust the look and size of the axes
                    cmfreqplot.set_xlim([cmfreq_x_min,cmfreq_x_max])
                    cmfreqplot.set_ylabel('Cumulative Frequency', fontsize=labelfontsize, fontweight='bold')
                    cmfreqplot.set_xlabel('Delta T', fontsize=labelfontsize, fontweight='bold')
                    cmfreqplot.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
                    cmfreqplot.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
                    cmfreqplot.tick_params(axis='both', direction='out', left='on'
                                         , labelleft='on', labelsize=ticklabelsize
                                         , which='major', length=10, width=1.5
                                         , gridOn='on')
                    cmfreqplot.tick_params(axis='both', direction='out', left='on'
                                         , labelleft='on', labelsize=ticklabelsize
                                         , which='minor', length=6, width=1.5)
                    
                    # SAVE OR SHOW FIGURES
                    if saveplotflag:
                        figurename_prefix = '.'.join(Filename.split('.')[:-1])
                        fig1.savefig((FIGDIR1 + '/' + figurename_prefix + 'normplot.eps'), dpi=None
                                         , facecolor='w', edgecolor='w'
                                         , orientation='portrait', papertype='letter'
                                         , format='eps', transparent=False
                                         , bbox_inches='tight', pad_inches=0.1
                                         , frameon=False)
                        fig1.savefig((FIGDIR1 + '/' + figurename_prefix + 'normplot.png'), dpi=None
                                         , facecolor='w', edgecolor='w'
                                         , orientation='portrait', papertype='letter'
                                         , format='png', transparent=False
                                         , bbox_inches='tight', pad_inches=0.1
                                         , frameon=False)
                        plt.close(fig1)
                        #
                        fig2.savefig((FIGDIR2 + '/' + figurename_prefix + 'cmfreqplot.eps'), dpi=None
                                         , facecolor='w', edgecolor='w'
                                         , orientation='portrait', papertype='letter'
                                         , format='eps', transparent=False
                                         , bbox_inches='tight', pad_inches=0.1
                                         , frameon=False)
                        fig2.savefig((FIGDIR2 + '/' + figurename_prefix + 'cmfreqplot.png'), dpi=None
                                         , facecolor='w', edgecolor='w'
                                         , orientation='portrait', papertype='letter'
                                         , format='png', transparent=False
                                         , bbox_inches='tight', pad_inches=0.1
                                         , frameon=False)
                        plt.close(fig2)
                    else:
                        plt.show()
                    # -----------------------------------------------
                    
                    #plotflag = False  # tmp option !!! PMB
                # -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
                
                # *************************************************************
            # END IF
        #
    #
#

outf.close()

#plt.show()