#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#       Fit a model to the means and standard deviations of
#       individual difference targets. Calculate the Sum of
#       Squared Error Residuals and the R^2 coefficient of determination.
#
#       Plot the models with the data points, Normal-Score plots of the
#       Residuals, and a scatter plot of the Residuals vs the Predicted
#       Values from the model fit to the data.
#
#       Run this script in the directory where the stats file is located.
#
#
#       Written 20200116. PMBremner
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import pandas as pd
import numpy as np
from numpy.linalg import inv
from scipy import stats

# qt5 error fix from Tim
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
import collections  # to use OrderDict
from copy import deepcopy as dp
from decimal import *


# Allow the use of latex
#rc('font', **{'family':'serif','serif':['Palatino']})
#rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font', **{'family':'sans-serif','sans-serif':['Lato']})
#rc('text',usetex=True)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#       FUNCTIONS
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#  A Polynomial Regression

# Setup and calculate the Least Squares solution of d = Gm
# for an over-determined system

# d = the data generated from a model

# m = the model used to describe the data

# G = the kernel matrix that transforms the model space into the data space
# G contains the independent variables, and is constructed from the
# partial derivatives of the data with respect to the model variables

# GT is the transpose of G

# The over-determined system solution is:
# m = [GT*G]^-1 *GT * d


def polyregression(independentV,dependentV,polyorder=1):
    
    # Get the size of each array so that it only needs to be done once
    DEPsize = dependentV.size
    #INDEPsize = independentV.size
    
    
    # Initialize the model estimate array
    # Add 1 to the polyorder to get the correct Model size
    # EXAMPLE: Linear model is polyorder = 1. So Model size = polyorder + 1 = 2
    Linsize = polyorder+1
    MODsize = Linsize
    mest = np.zeros(MODsize)
    
    
    #print (independentV,dependentV)  # just a quick check to ensure correct values obtained
    
    # Initialize the kernel matrix, G, and its transpose, GT
    G = np.zeros((DEPsize,MODsize),dtype=np.float)
    GT = np.zeros((MODsize,DEPsize),dtype=np.float)
    
    
    # Populate G and GT
    for i in range(DEPsize):
        for j in range(MODsize):
            G[i][j] = np.power(independentV[i],j)
            GT[j][i] = np.power(independentV[i],j)
    #
    
    # Convert G and GT into matrix objects
    G = np.matrix(G)
    GT = np.matrix(GT)
    
    
    # Calculate the model from the data
    GTG = np.dot(GT,G)
    GTGinv = inv(GTG)
    GTGinvGT = np.dot(GTGinv,GT)
    mest = np.dot(GTGinvGT,dependentV)
    
    
    #print (mest)
    
    
    # Calculate the estimated data, dest, from the model estimate, mest
    dest = np.zeros(DEPsize,dtype=np.float)
    #for i in range(dest.size):
    #    dest[i] = np.dot(G[i],np.matrix.transpose(mest))
    ##
    dest = np.dot(G,np.matrix.transpose(mest))
    
    
    # Calculate the arithmetic mean of dest
    #dest_bar = (np.sum(dest))/dest.size
    dependentV_bar = (np.sum(dependentV))/DEPsize
    
    
    # Calculate the Sum of Squares of dest
    SS_dest = 0.0
    for d in dependentV:
        #SS_dest = (d-dest_bar)*(d-dest_bar)
        SS_dest = SS_dest + (d-dependentV_bar)*(d-dependentV_bar)
    #
    
    
    # Calculate the error residuals
    residuals = np.zeros(DEPsize,dtype=np.float)
    for i in range(DEPsize):
        residuals[i] = dependentV[i] - dest[i]
    #
    
    
    # Calculate the Sum of Squared Error Residuals
    SSE = 0
    SSE = np.sum(residuals*residuals)
    
    
    # Determine the R2 (R^2), or coefficient of determination
    R2 = 1 - (SSE/SS_dest)
    
    
    #print (R2)
    
    return mest,dest,residuals,R2


# Linearized form of the Power Law regression
#
# y = C_0 x^b  =>  log(y) = C_1 + C_2*log(x)
# where C_0 = 10^C_1  and  b = C_2
#
# Since taking the logs of the vectors results in a
# linear form of the equation, this function simply
# converts the values then calls polyregression to
# do the inversion, polyorder=1.
# Values are reverted back to real values at the
# end before returning values to the main program.
#
def powerlaw_regression(independentV,dependentV,polyorder=1):

    # Convert values to log base 10 values
    log_indV = np.log10(independentV, dtype=np.float)
    log_depV = np.log10(dependentV, dtype=np.float)

    # Call polyregression to calculate the linear model of log values
    mest,dest,residuals,R2 = polyregression(log_indV, log_depV, 1)

    # Take the base 10 exponential to revert back to real values
    mest[0,0] = np.power(10,mest[0,0])
    dest = np.power(10,dest)

    print ('Shape of power law model: mest={}  dest={}\n'.format(mest.shape,dest.shape))
    print ('Components of mest: mest[0]={}  mest[0][0]={}\n'.format(mest[0],mest[0,0]))
    print ('Squeezed model matrix: mest={}\n'.format(np.squeeze(np.asarray(mest))))

    return mest,dest,residuals,R2
# -----------------------------------------------------------------------------

#  T-Tests

# Perform one of two T-Tests:
# 1) Equal Variance (Pooled) T-Test
# 2) Unequal Variance T-Test
    
# When to use 1) or 2):
# Question 1) Are the two sample sets the same or related?
#       Yes -- use a Paired (Dependent) T-Test instead
#       No -- Use test 1) or 2)
# Question 2) Are the two sample sets of the same size?
#       Yes -- use test 1)
#       No -- Move to next question
# Question 3) Do the two sample sets have the same variances?
#       Yes -- use test 1)
#       No -- use test 2)

def t_test_equalvar (mean1,mean2,factor1,factor2,n1,n2,type='stddev'):
    
    # Calculate:
    # degrees of freedom = dof = n1 + n2 - 2
    #
    # t = (mean1 - mean2) / [s * sqrt((1/n1) + (1/n2))]
    #
    # where s = pooled sample standard deviation
    # s = sqrt( ((n1-1)*var1 + (n2-1)*var2) / dof )
    #
    # where n1,n2 = number of samples in datasets 1 and 2
    # var1,var2 = the variances of datasets 1 and 2
    # and mean1,mean2 = the mean values of datasets 1 and 2
    
    # Test if type = stddev, or variance
    if (type=='stddev'):
        var1 = factor1*factor1
        var2 = factor2*factor2
    elif (type=='variance'):
        var1 = factor1
        var2 = factor2
    else:
        print ('Invalid type input for t-test. Use "stddev" or "variance"')
        raise NameError(type)
    # END IF
    
    
    # Calculate the Degrees of Freedom (dof)
    dof = n1 + n2 - 2
    
    # Calculate the t_value
    numerator = mean1 - mean2
    
    # Calculate the pooled sample standard deviation, pssd
    pssd = ((n1-1)*var1)+((n2-1)*var2)
    pssd = pssd/dof
    pssd = np.sqrt(pssd)
    
    # Calculate the other factor in the denominator
    rootfactor = (1.0/n1) + (1.0/n2)
    rootfactor = np.sqrt(rootfactor)
    
    t_value = numerator / (pssd * rootfactor)
    
    return t_value,dof

def t_test_unequalvar (mean1,mean2,factor1,factor2,n1,n2,type='stddev'):
    
    # Calculate:
    # degrees of freedom = dof = [numerator] / [denominator]
    # where numerator = [ (var1/n1) + (var2/n2) ]^2
    # and denominator = [ ((var1/n1)^2/(n1-1)) + ((var2/n2)^2/(n2-1)) ]
    #
    # t = (mean1 - mean2) / sqrt( (var1/n1) + (var2/n2) )
    #
    # where n1,n2 = number of samples in datasets 1 and 2
    # var1,var2 = the variances of datasets 1 and 2
    # and mean1,mean2 = the mean values of datasets 1 and 2
    
    # Test if type = stddev, or variance
    if (type=='stddev'):
        var1 = factor1*factor1
        var2 = factor2*factor2
    elif (type=='variance'):
        var1 = factor1
        var2 = factor2
    else:
        print ('Invalid type input for t-test. Use "stddev" or "variance"')
        raise NameError(type)
    # END IF
    
    
    # Calculate some needed parts first
    part1 = (var1 / n1)
    part2 = (var2 / n2)
    part3 = part1 + part2
    
    # Calculate the Degrees of Freedom (dof)
    frac1 = (part1*part1) / (n1-1)
    frac2 = (part2*part2) / (n2-1)
    dof = (part3*part3) / (frac1 + frac2)
    
    # Calculate the t_value
    numerator = mean1 - mean2
    
    t_value = numerator / np.sqrt(part3)
    
    return t_value,dof

# -----------------------------------------------------------------------------

#   Rank-Sum

def rank_arrays(xi,yi):
    
    # Concatenate the arrays
    xiyi=np.concatenate((xi,yi),axis=0)
    
    # Determine if there are ties
    # Start by assuming 'no', then
    # count the unique array values
    # and any that do not have only
    # a single occurance
    tie = 'no'
    a,b=np.unique(xiyi,return_counts=True)
    for i in b:
        if (i!=1): tie='yes'
        if (tie=='yes'): break
    # END for loop
    
    # Rank the new concatenated array
    xiyi_ranks = stats.rankdata(xiyi,method='average')
    
    # Retrieve the ranked arrays, separated
    xi_ranks = np.asarray(xiyi_ranks[:xi.size],dtype=np.float)
    yi_ranks = np.asarray(xiyi_ranks[xi.size:],dtype=np.float)
    
    return xi_ranks, yi_ranks, tie


def sum_ranks(ranks_x,ranks_y):
    
    Wrs_x = np.sum(ranks_x)
    Wrs_y = np.sum(ranks_y)
    
    return Wrs_x,Wrs_y


def compute_joint_ranks(Wrs_x,Wrs_y):
    
    Wrs_joint = np.sum([Wrs_x,Wrs_y])
    
    return Wrs_joint


def rank_sum_exact_test(xi,yi):
    
    # Determine the ranks of each data point (jointly)
    ranks_x,ranks_y,tie = rank_arrays(xi,yi)
    
    
    # Sum the ranks
    Wrs_x, Wrs_y = sum_ranks(ranks_x,ranks_y)
    #if (tie == 'yes'):
    #    Wrs_joint = compute_joint_ranks(Wrs_x,Wrs_y)
    
    
    # Determine the size of the ranks arrays
    nx = ranks_x.size
    ny = ranks_y.size
    
    # Determine the joint size of the ranks arrays
    #N_joint = nx + ny
    
    # Determine which ranks set is smaller: ranks_x or ranks_y
    if (nx < ny):
        #n = nx
        Wrs = Wrs_x
    elif (nx >= ny):
        #n = ny
        Wrs = Wrs_y
    # End IF
    
    return Wrs


def rank_sum_large_sample_approximation(xi,yi):
    
    # Determine the ranks of each data point (jointly)
    ranks_x,ranks_y,tie = rank_arrays(xi,yi)
    
    
    # Sum the ranks
    Wrs_x, Wrs_y = sum_ranks(ranks_x,ranks_y)
    #if (tie == 'yes'):
    #    Wrs_joint = compute_joint_ranks(Wrs_x,Wrs_y)
    
    
    # Determine the size of the ranks arrays
    nx = ranks_x.size
    ny = ranks_y.size
    
    # Determine the joint size of the ranks arrays
    N_joint = nx + ny
    
    # Determine which ranks set is smaller: ranks_x or ranks_y
    if (nx < ny):
        n = nx
        Wrs = Wrs_x
    elif (nx >= ny):
        n = ny
        Wrs = Wrs_y
    # End IF
    
    
    # Calculate the mean and standard deviation of the ranks distribution
    mu_w = n*(N_joint+1)/2.0
    if (tie == 'no'):
        sigma_w = np.sqrt((nx*ny*(N_joint+1))/12.0)
    elif (tie == 'yes'): # Correction for ties
        # Determine the sum of squares of each ranks element
        ranks2 = np.sum([(ranks_x*ranks_x),(ranks_y*ranks_y)],dtype=np.float)
        
        # Calculate the corrected standard deviation
        sigma_w = np.sqrt(
                        ((nx*ny)/(N_joint*(N_joint-1.0)))*
                        ranks2-
                        ((nx*ny*(N_joint+1.0)*(N_joint+1.0))/(4.0*(N_joint-1.0)))
                        )
    # End if
    
    
    # Determine which is larger: Wrs or mu_w
    # of the smaller set of ranks: ranks_x or ranks_y
    # Calculate Zrs, the standardized form of the test statistic
    # accordingly
    if (Wrs > mu_w):
        Zrs = (Wrs - 0.5 - mu_w)/sigma_w
    elif (Wrs == mu_w): Zrs = 0.0
    elif (Wrs < mu_w):
        Zrs = (Wrs + 0.5 - mu_w)/sigma_w
    # End IF
    
    return Zrs
    
    


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#       MAIN PROGRAM
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

filename = '../../Data/stats_table_case_007h_pred_2009_2035_NFPA_calib_only.dat.xlsx'

print (filename)


# Define the PATH to the Figure Directories
# FIGDIR1 = '../../Data/Figures/RegressionModel_PolyRegrOrder'
# FIGDIR2 = '../../Data/Figures/RegressionModel_ResidVsPred'
# FIGDIR3 = '../../Data/Figures/RegressionModel_ResidVsNorm'
FIGDIR1 = '../../Data/Figures_tmp/RegressionModel_PolyRegrOrder'
FIGDIR2 = '../../Data/Figures_tmp/RegressionModel_ResidVsPred'
FIGDIR3 = '../../Data/Figures_tmp/RegressionModel_ResidVsNorm'


# Specify the sheets that need to be processed
#
# NOTE: * The sheetnames MUST match the names in the Excel file
#       ** The Observation Class names MUST be in the same order as the sheetnames
#
# ----------------------------------------------
sheetnames = ['DeltaT_NonZero_baseflow',
              'DeltaT_NonZero_springflow',
              'DeltaT_NonZero_SAS',
              'DeltaT_NonZero_UFA',
              'DeltaT_NonZero_LFA']

obsnames = ['Baseflow',
            'Spring_flow',
            'SAS',
            'UFA',
            'LFA']

# Set the bsfl value to the element of the Baseflow defined above
# NOTE: if no baseflow included, set to any unused value
bsfl = 0

# Set the springs value to the element of the Spring_flow defined above
# NOTE: if no springs included, set to any unused value
springs = 1
# ----------------------------------------------


print ('Begin reading Excel file and assigning values . . .')



# Get the means and standard deviations
# and assign them to numpy arrays
column1 = 'Mean_AdjustSigWhereNeed'
column2 = 'Stdv_AdjustSign'
column3 = 'NumSamp_SP2-1'
column_powerX = 'abs(Mean_SP2-1)'
column_powerY = 'Stdv/abs(Mean)_SP2-1'

# Construct a dictionary to hold the values from each observation class
placeholder=np.zeros(2)
masterlist = collections.OrderedDict([('obsname','name'),
                                      ('mean',placeholder),
                                      ('stddv',placeholder),
                                      ('numsamp',placeholder),
                                      ('nsize',0),
                                      ('variance',placeholder),
                                      ('absMean',placeholder),
                                      ('StddvOverAbsMean',placeholder)])
obsclasses = []
for i in range(len(sheetnames)): obsclasses.append(dp(masterlist))


for i in range(len(sheetnames)):
    dataf = pd.read_excel(filename,sheet_name=sheetnames[i],header=0)
    obsclasses[i]['obsname'] = obsnames[i]
    obsclasses[i]['mean'] = np.array(dataf[column1],dtype=np.float)
    obsclasses[i]['stddv'] = np.array(dataf[column2],dtype=np.float)
    obsclasses[i]['numsamp'] = np.array(dataf[column3],dtype=np.int)
    obsclasses[i]['nsize'] = obsclasses[i]['mean'].size
    obsclasses[i]['variance'] = np.zeros(obsclasses[i]['nsize'],dtype=np.float)
    obsclasses[i]['variance'] = [j*j for j in obsclasses[i]['stddv']]
    obsclasses[i]['absMean'] = np.array(dataf[column_powerX], dtype=np.float)
    obsclasses[i]['StddvOverAbsMean'] = np.array(dataf[column_powerY], dtype=np.float)
# End FOR LOOP


# Define plot properties
# -----------------------------------------------------------------------------
linewidth_mod = 3.50 #2.5 #4.0 #2.5
linewidth_ref = 3.50 #1.5 #2.5 #4.0
#
linestyle_mod = '--'
linestyle_ref = '--'
#
color_sample = 'tab:blue'
color_ref = 'tab:orange'
color_drawdown = 'tab:red'
color_rebound = 'tab:pink'
#
label_sample = 'mean(Target Diff)'
label_ref = 'Ref Norm Distribution'
label_mod = 'Est Model'
label_modpt = 'Est Model Value'
label_drawdown = 'Draw Down'
label_rebound = 'Rebound'
label_reduc = 'Reduction'
label_incr = 'Increase'
#
textfontsize = 14 #14 #12
legendfontsize = 11.5
textstyle = 'normal' #'bold'
labelfontsize = 18 #20 #18 #24
ticklabelsize = 14 #16 #14
#
# Set whether to actually make the plots, or just make the output data file
plotflag = True
#
# Set whetherto simply show the plots (False), or save them (True)
saveplotflag = True
# -----------------------------------------------------------------------------

# Designate an integer to indicate when to use the Power Law
PowLaw_flag = np.int(999)

#tmpobsclasslist = [2]
#for i_obs in tmpobsclasslist:
for i_obs in range(len(sheetnames)):
    # i_obs = 4

    # Select the polynomial order -- order 1 = line, 2 = parabola, etc.
    polyorder_range = np.asarray([1,2,3,4,5,PowLaw_flag],dtype=np.int)
    #polyorder_range = np.asarray([6],dtype=np.int)
    for polyorder in polyorder_range:
        if (polyorder==PowLaw_flag):
            # Assign the columns of interest to numpy arrays
            independentV = obsclasses[i_obs]['absMean']  # absolute value of the mean
            dependentV = obsclasses[i_obs]['StddvOverAbsMean']  # stddv over the absolute value of the mean
            numsamp = obsclasses[i_obs]['numsamp']  # number of samples
            n = obsclasses[i_obs]['nsize']  # number of means
            var = obsclasses[i_obs]['variance']  # variance

            # Call the Power Law regression instead
            mest, dest, residuals, R2 = powerlaw_regression(independentV, dependentV, polyorder)
        else:
            # Assign the columns of interest to numpy arrays
            independentV = obsclasses[i_obs]['mean']  # mean
            dependentV = obsclasses[i_obs]['stddv']  # stddv
            numsamp = obsclasses[i_obs]['numsamp']  # number of samples
            n = obsclasses[i_obs]['nsize']  # number of means
            var = obsclasses[i_obs]['variance']  # variance

            # Call the polyregression
            mest,dest,residuals,R2 = polyregression(independentV,dependentV,polyorder)
        #
        mest=np.asarray(mest,dtype=np.float)
        dest=np.asarray(dest,dtype=np.float)
        print ('Polyorder={}  mest={}\n'.format(polyorder,mest))
        print ('R2={}\n'.format(R2))
        
        
        # Calculate some points from the model in order to plot the line
        # ------------------------------------------------
        numpoints = 300 #polyorder * 20
        independentV_embellished = np.linspace(independentV.min(),independentV.max(),numpoints, dtype=np.float)
        
        # Specify the size of the regression model
        if not polyorder==PowLaw_flag:
            modcomps = np.linspace(0, polyorder, polyorder+1, dtype=np.int)
        
        ordered_indpdV_dest = np.zeros([independentV_embellished.size,2],dtype=np.float)
        for pnt in range(independentV_embellished.size):
            est_data = np.float(0.0)
            if (polyorder==PowLaw_flag):
                est_data = mest[0,0]*np.power(independentV_embellished[pnt],mest[0,1])
            else:
                for comp in modcomps:
                    est_data = est_data + mest[0,comp]*np.power(independentV_embellished[pnt],comp) #mest[0][comp]
                #
            #
            # for imod,mod in enumerate(mest[0]):
            #     est_data = est_data + mod*np.power(independentV_embellished[pnt],(imod))
            # #
            ordered_indpdV_dest[pnt] = [independentV_embellished[pnt],est_data]
        #
        #ordered_indpdV_dest = np.sort(ordered_indpdV_dest,axis=0,kind='mergesort')
        # ________________________________________________
        # endmod endmod endmod endmod endmod endmod endmod
        
        
        # Get the ordered list of plotting positions for the Residual data
        pos = np.zeros(residuals.size,dtype=np.float)
        residuals_ordered = np.sort(residuals,kind='mergesort')
        for ip in range(pos.size): pos[ip] = ((ip+1)-0.4) / (pos.size+0.2)
        # NEED Z-SCORES !!! PMB !!!
        pos = stats.zscore(pos)
        
        
        # Calculate the mean and standard deviation of the residuals
        residuals_mean = np.sum(residuals) / residuals.size
        residuals_stddv = np.sqrt( np.sum((residuals-residuals_mean)*(residuals-residuals_mean)) / residuals.size )
        
        # Construct the reference line
        # The line will be y = residuals_mean + (residuals_stddv * x)
        # Find only the extreme two points to be able to draw the line
        refy_1 = residuals_mean + (residuals_stddv * pos.min())
        refy_2 = residuals_mean + (residuals_stddv * pos.max())
        refline_x = np.array([pos.min(),pos.max()],dtype=np.float)
        refline_y = np.array([refy_1,refy_2],dtype=np.float)
        
        ## Calculate the mean of the independentV and dependentV values
        #independentV_mean = np.sum(independentV) / independentV.size
        #dependentV_mean = np.sum(dependentV) / dependentV.size
        #
        ## Calculate the standard deviation of the independentV and dependentV values
        #independentV_stddv = np.sqrt( np.sum((independentV-independentV_mean)*(independentV-independentV_mean)) / independentV.size )
        #dependentV_stddv = np.sqrt( np.sum((dependentV-dependentV_mean)*(dependentV-dependentV_mean)) / dependentV.size )
        
        
#        for i in range(n-1):
#            t_value,dof = t_test_equalvar(independentV[i],independentV[i+1],var[i],var[i+1],numsamp[i],numsamp[i+1],'variance')
#            #t_value,dof = t_test_equalvar(mean_bf[i],mean_bf[i+1],var_bf[i],var_bf[i+1],522,522,'variance')
#            #print (t_value, dof)
#        #
        
        
        # Use the rank-sum approach
        # ----------------------------------------------
        # TMP Assignment !!! PMB !!!
        xi = np.asarray([1.7,0.59,4.0,1.1,0.87,1.2,1.1,1.6,1.3,3.2],dtype=np.float)
        yi = np.asarray([9.7,0.7,0.9,0.3,1.3,0.7,0.92,0.36,1.0,0.5],dtype=np.float)
        
        # The Exact Test Statistic
        Wrs = rank_sum_exact_test(xi,yi)
        print (Wrs)
        
        # The Large Sample Approximation
        Zrs = rank_sum_large_sample_approximation(xi,yi)
        print (Zrs)
        # ----------------------------------------------
        
        
        # Construct the plot label
        quantizeprec='0.00000000001'
        if (polyorder==1):
            polylabel = (r'$y = c_0 + c_{1}x$')
            #polylabel = (r'$y = a + bx$')
            modlabel = (r'$[c_{0} = $' + str(Decimal(mest[0,0]).quantize(Decimal(quantizeprec)))
                        + ', $c_{1} = $' + str(Decimal(mest[0,1]).quantize(Decimal(quantizeprec)))
                        + '$]$')
        elif (polyorder==2):
            polylabel = (r'$y = c_0 + c_{1}x + c_{2}x^2$')
            #polylabel = (r'$y = a + bx + cx^2$')
            modlabel = (r'$[c_{0} = $' + str(Decimal(mest[0,0]).quantize(Decimal(quantizeprec)))
                        + ', $c_{1} = $' + str(Decimal(mest[0,1]).quantize(Decimal(quantizeprec)))
                        + '\n..., $c_{2} = $' + str(Decimal(mest[0,2]).quantize(Decimal(quantizeprec)))
                        + '$]$')
        elif (polyorder==3):
            polylabel = (r'$y = c_0 + c_{1}x + c_{2}x^2 + c_{3}x^3$')
            #polylabel = (r'$y = a + bx + cx^2 + dx^3$')
            modlabel = (r'$[c_{0} = $' + str(Decimal(mest[0,0]).quantize(Decimal(quantizeprec)))
                        + ', $c_{1} = $' + str(Decimal(mest[0,1]).quantize(Decimal(quantizeprec)))
                        + '\n..., $c_{2} = $' + str(Decimal(mest[0,2]).quantize(Decimal(quantizeprec)))
                        + ', $c_{3} = $' + str(Decimal(mest[0,3]).quantize(Decimal(quantizeprec)))
                        + '$]$')
        elif (polyorder==4):
            polylabel = (r'$y = c_0 + c_{1}x + c_{2}x^2 + c_{3}x^3 + c_{4}x^4$')
            #polylabel = (r'$y = a + bx + cx^2 + dx^3 + ex^4$')
            modlabel = (r'$[c_{0} = $' + str(Decimal(mest[0,0]).quantize(Decimal(quantizeprec)))
                        + ', $c_{1} = $' + str(Decimal(mest[0,1]).quantize(Decimal(quantizeprec)))
                        + '\n..., $c_{2} = $' + str(Decimal(mest[0,2]).quantize(Decimal(quantizeprec)))
                        + ', $c_{3} = $' + str(Decimal(mest[0,3]).quantize(Decimal(quantizeprec)))
                        + '\n..., $c_{4} = $' + str(Decimal(mest[0,4]).quantize(Decimal(quantizeprec)))
                        + '$]$')
        elif (polyorder==5):
            polylabel = (r'$y = c_0 + c_{1}x + c_{2}x^2 + c_{3}x^3 + c_{4}x^4 + c_{5}x^5$')
            #polylabel = (r'$y = a + bx + cx^2 + dx^3 + ex^4 + fx^5$')
            modlabel = (r'$[c_{0} = $' + str(Decimal(mest[0,0]).quantize(Decimal(quantizeprec)))
                        + ', $c_{1} = $' + str(Decimal(mest[0,1]).quantize(Decimal(quantizeprec)))
                        + '\n..., $c_{2} = $' + str(Decimal(mest[0,2]).quantize(Decimal(quantizeprec)))
                        + ', $c_{3} = $' + str(Decimal(mest[0,3]).quantize(Decimal(quantizeprec)))
                        + '\n..., $c_{4} = $' + str(Decimal(mest[0,4]).quantize(Decimal(quantizeprec)))
                        + ', $c_{5} = $' + str(Decimal(mest[0,5]).quantize(Decimal(quantizeprec)))
                        + '$]$')
        elif (polyorder==6):
            polylabel = (r'$y = c_0 + c_{1}x + c_{2}x^2 + c_{3}x^3 + c_{4}x^4 + c_{5}x^5 + c_{6}x^6$')
            #polylabel = (r'$y = a + bx + cx^2 + dx^3 + ex^4 + fx^5 + gx^6$')
            modlabel = (r'$[c_{0} = $' + str(Decimal(mest[0,0]).quantize(Decimal(quantizeprec)))
                        + ', $c_{1} = $' + str(Decimal(mest[0,1]).quantize(Decimal(quantizeprec)))
                        + '\n..., $c_{2} = $' + str(Decimal(mest[0,2]).quantize(Decimal(quantizeprec)))
                        + ', $c_{3} = $' + str(Decimal(mest[0,3]).quantize(Decimal(quantizeprec)))
                        + '\n..., $c_{4} = $' + str(Decimal(mest[0,4]).quantize(Decimal(quantizeprec)))
                        + ', $c_{5} = $' + str(Decimal(mest[0,5]).quantize(Decimal(quantizeprec)))
                        + '\n..., $c_{6} = $' + str(Decimal(mest[0,6]).quantize(Decimal(quantizeprec)))
                        + '$]$')
        elif (polyorder==PowLaw_flag):
            polylabel = (r'$y = c x^b$')
            # polylabel = (r'$y = c_0 * x^b$')
            modlabel = (r'$[c = $' + str(Decimal(mest[0, 0]).quantize(Decimal(quantizeprec)))
                        + ', $b = $' + str(Decimal(mest[0, 1]).quantize(Decimal(quantizeprec)))
                        + '$]$')
        #
        #labelstr = (obsnames[i_obs] + '\nPoly Regression Order ' + str(polyorder))
        #labelstr = (obsnames[i_obs] + '\nPoly Regression Order: ' + str(polyorder)
        #            + ' , R2: ' + str(R2) + '\nRank Sum Test Result: ' + str(Zrs))
        if (polyorder == PowLaw_flag):
            labelstr = (obsnames[i_obs] + '\nPower Law Regression'
                        + ' , $\mathrm{R}^{2}$: ' + str(R2) + '\nPower Law Model: ' + str(polylabel)
                        + '\n' + modlabel)
        else:
            labelstr = (obsnames[i_obs] + '\nPoly Regression Order: ' + str(polyorder)
                        + ' , $\mathrm{R}^{2}$: ' + str(R2) + '\nPolynomial Model:\n' + str(polylabel)
                        + '\n' + modlabel)
        
        if (i_obs==bsfl):
            reftype_neg_label =dp(label_reduc)
            reftype_pos_label =dp(label_incr)
        elif (i_obs==springs):
            reftype_neg_label =dp(label_reduc)
            reftype_pos_label =dp(label_rebound)
        else:
            reftype_neg_label =dp(label_drawdown)
            reftype_pos_label =dp(label_rebound)
        
        
        
        # -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
        if plotflag:
            
            # Make plots
            #fig1,(resplot,normplot) = plt.subplots(2,1,figsize=(5,10))
            
            # Make and Save Plots
            # -----------------------------------------------
            # plot = fig.add_axes([left corner,bottom corner,width,height], label='name')
            fig1 = plt.figure(figsize=(10,6))
            regrplot = fig1.add_axes([0.0,0.0,0.70,1.0], label='axes1')
            fig2 = plt.figure(figsize=(10,6))
            resplot = fig2.add_axes([0.0,0.0,0.70,1.0], label='axes2')
            fig3 = plt.figure(figsize=(10,6))
            normplot = fig3.add_axes([0.0,0.0,0.70,1.0], label='axes3')
            
            # Turn off the top and right axes
            fig1.gca().spines['top'].set_visible(False)
            fig1.gca().spines['right'].set_visible(False)
            fig2.gca().spines['top'].set_visible(False)
            fig2.gca().spines['right'].set_visible(False)
            fig3.gca().spines['top'].set_visible(False)
            fig3.gca().spines['right'].set_visible(False)
            
            # Plot the data points with the model fit
            # -----------------------------------------------
            # Plot the data and the reference norm distribution
            
            regrplot.plot(ordered_indpdV_dest[:,0], ordered_indpdV_dest[:,1], label=label_mod
                          , color=color_ref, linestyle=linestyle_mod, linewidth=linewidth_mod, zorder=3)
            regrplot.scatter(independentV,dest, marker='o', s=100, label=label_modpt
                             , color=color_ref, edgecolor=color_ref, linewidths=1.00, zorder=3)
            regrplot.scatter(independentV,dependentV, marker='D', s=35, label=label_sample
                             , color=color_sample, edgecolor=color_sample, linewidths=1.00, zorder=4, alpha=0.75)
            
            # bbox_to_anchor = (x, y, width, height)
            legend1 = regrplot.legend(fontsize=legendfontsize, loc='upper left'
                                      , bbox_to_anchor=(1.01,1.0), borderaxespad=0.0
                                      , handlelength=3.0, framealpha=1.0)
            
            ymin,ymax = regrplot.get_ylim()
            xmin,xmax = regrplot.get_xlim()
            xspan_neg = [xmin,0]
            xspan_pos = [0,xmax]
            #regrplot.axhspan(tickmin,0,color='lightgray')
            #fig1.gca().fill_between(independentV,ymin,0,color='lightgray') #,where=independentV<=0
            if (ymin<0): fig1.gca().fill_between(xspan_neg,ymin,0,color='lightgray',alpha=0.75,zorder=1,label=reftype_neg_label)
            if (ymax>0): fig1.gca().fill_between(xspan_pos,0,ymax,color='mistyrose',alpha=0.75,zorder=1,label=reftype_pos_label)
            
            
            # Label the plot
            text_y = 0.0 #0.54 #0.97 # in axis coordinates
            text_x = 1.02 #0.05 # in axis coordinates
            regrplot.text(text_x,text_y, labelstr, fontsize=textfontsize, fontweight=textstyle
                          , color='k', verticalalignment='bottom', horizontalalignment='left'
                          , rotation=90, transform=regrplot.transAxes
                          , bbox=dict(edgecolor='none',facecolor='w',alpha=1.0), zorder=6)
            
            
            # Adjust the look and size of the axes
            #normplot.set_ylim([norm_y_min,norm_y_max])
            #deltachar = u'\u0394'
            #regrplot.set_ylabel('Delta T Stddev (sign adjusted)', fontsize=labelfontsize, fontweight='bold')
            # regrplot.set_xlabel('Delta T Mean', fontsize=labelfontsize, fontweight='bold')
            if (polyorder == PowLaw_flag):
                regrplot.set_ylabel(r'$s_{\Delta \mathrm{T}} / \left\| m_{\Delta \mathrm{T}} \right\|$', fontsize=labelfontsize, fontweight='bold')
                regrplot.set_xlabel(r'$\left\| m_{\Delta \mathrm{T}} \right\|$', fontsize=labelfontsize, fontweight='bold')
            else:
                regrplot.set_ylabel(r'$s_{\Delta \mathrm{T}}$ (sign adjusted)', fontsize=labelfontsize, fontweight='bold')
                regrplot.set_xlabel(r'$m_{\Delta \mathrm{T}}$', fontsize=labelfontsize, fontweight='bold')
            #

            regrplot.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
            regrplot.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
            regrplot.tick_params(axis='both', direction='out', left='on'
                                 , labelleft='on', labelsize=ticklabelsize
                                 , which='major', length=10, width=1.5
                                 , gridOn='on')
            regrplot.tick_params(axis='both', direction='out', left='on'
                                 , labelleft='on', labelsize=ticklabelsize
                                 , which='minor', length=6, width=1.5)
            # _______________________________________________
            # endf1 endf1 endf1 endf1 endf1 endf1 endf1 endf1
            
            
            # Plot up residuals
            # -----------------------------------------------
            
            # Plot a line along y=zero
            zeroline_x = np.array([np.min(dest),np.max(dest)],dtype=np.float)
            zeroline_y = np.array([0,0],dtype=np.float)
            
            # Plot the data and the reference norm distribution
            resplot.plot(zeroline_x,zeroline_y, label=label_ref, color=color_ref
                         , linestyle=linestyle_ref, linewidth=linewidth_ref, zorder=3)
            #resplot.plot(zeroline_x,zeroline_y,color='k',linewidth=0.75)
            resplot.scatter(dest,residuals, marker='D', s=35, label=label_sample
                             , color=color_sample, edgecolor=color_sample, linewidths=1.00, zorder=4, alpha=0.75)
            #resplot.scatter(dest,residuals,color='seagreen',edgecolor='k',linewidths=0.5)
            
            
            
            # bbox_to_anchor = (x, y, width, height)
            legend1 = resplot.legend(fontsize=legendfontsize, loc='upper left'
                                      , bbox_to_anchor=(1.01,1.0), borderaxespad=0.0
                                      , handlelength=3.0, framealpha=1.0)
            
            # Label the plot
            text_y = 0.0 #0.54 #0.97 # in axis coordinates
            text_x = 1.02 #0.05 # in axis coordinates
            resplot.text(text_x,text_y, labelstr, fontsize=textfontsize, fontweight=textstyle
                          , color='k', verticalalignment='bottom', horizontalalignment='left'
                          , rotation=90, transform=resplot.transAxes
                          , bbox=dict(edgecolor='none',facecolor='w',alpha=1.0), zorder=6)
            
            
            # Adjust the look and size of the axes
            #normplot.set_ylim([norm_y_min,norm_y_max])
            resplot.set_ylabel('residuals', fontsize=labelfontsize, fontweight='bold')
            resplot.set_xlabel('predicted values', fontsize=labelfontsize, fontweight='bold')
            resplot.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
            resplot.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
            resplot.tick_params(axis='both', direction='out', left='on'
                                 , labelleft='on', labelsize=ticklabelsize
                                 , which='major', length=10, width=1.5
                                 , gridOn='on')
            resplot.tick_params(axis='both', direction='out', left='on'
                                 , labelleft='on', labelsize=ticklabelsize
                                 , which='minor', length=6, width=1.5)
            # _______________________________________________
            # endf2 endf2 endf2 endf2 endf2 endf2 endf2 endf2
            
            
            
            # Plot the Residuals in a Normal-Scores Plot
            # -----------------------------------------------
            # Plot the data and the reference norm distribution
            normplot.plot(refline_x, refline_y, label=label_ref, color=color_ref
                          , linestyle=linestyle_ref, linewidth=linewidth_ref, zorder=3)
            #normplot.plot(refline_x,refline_y,color='k',linewidth=0.75)
            normplot.scatter(pos,residuals_ordered, marker='D', s=35, label=label_sample
                             , color=color_sample, edgecolor=color_sample, linewidths=1.00, zorder=4, alpha=0.75)
            #normplot.scatter(pos,residuals_ordered,color='seagreen',edgecolor='k',linewidths=0.5)
            
            
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
            #normplot.set_ylim([norm_y_min,norm_y_max])
            normplot.set_ylabel('residuals', fontsize=labelfontsize, fontweight='bold')
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
            # _______________________________________________
            # endf3 endf3 endf3 endf3 endf3 endf3 endf3 endf3
            
            
            
            # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            
            
            # SAVE OR SHOW FIGURES
            if saveplotflag:
                if (polyorder == PowLaw_flag):
                    figurename_prefix = ('RegressionModel_PowerLawRegression_' + obsnames[i_obs])
                else:
                    figurename_prefix = ('RegressionModel_PolyRegrOrder_' + str(polyorder) + '_' + obsnames[i_obs])
                #
                fig1.savefig((FIGDIR1 + '/' + figurename_prefix + '.eps'), dpi=None
                                 , facecolor='w', edgecolor='w'
                                 , orientation='portrait', papertype='letter'
                                 , format='eps', transparent=False
                                 , bbox_inches='tight', pad_inches=0.1
                                 ) #frameon=False)
                fig1.savefig((FIGDIR1 + '/' + figurename_prefix + '.png'), dpi=None
                                 , facecolor='w', edgecolor='w'
                                 , orientation='portrait', papertype='letter'
                                 , format='png', transparent=False
                                 , bbox_inches='tight', pad_inches=0.1
                                 )
                plt.close(fig1)
                #
                if (polyorder == PowLaw_flag):
                    figurename_prefix = ('RegressionModel_ResidVsPred_PowerLawRegression_' + obsnames[i_obs])
                else:
                    figurename_prefix = ('RegressionModel_ResidVsPred_PolyRegrOrder_' + str(polyorder) + '_' + obsnames[i_obs])
                fig2.savefig((FIGDIR2 + '/' + figurename_prefix + '.eps'), dpi=None
                                 , facecolor='w', edgecolor='w'
                                 , orientation='portrait', papertype='letter'
                                 , format='eps', transparent=False
                                 , bbox_inches='tight', pad_inches=0.1
                                 )
                fig2.savefig((FIGDIR2 + '/' + figurename_prefix + '.png'), dpi=None
                                 , facecolor='w', edgecolor='w'
                                 , orientation='portrait', papertype='letter'
                                 , format='png', transparent=False
                                 , bbox_inches='tight', pad_inches=0.1
                                 )
                plt.close(fig2)
                #
                if (polyorder == PowLaw_flag):
                    figurename_prefix = ('RegressionModel_ResidVsNorm_PowerLawRegression_' + obsnames[i_obs])
                else:
                    figurename_prefix = ('RegressionModel_ResidVsNorm_PolyRegrOrder_' + str(polyorder) + '_' + obsnames[i_obs])
                fig3.savefig((FIGDIR3 + '/' + figurename_prefix + '.eps'), dpi=None
                                 , facecolor='w', edgecolor='w'
                                 , orientation='portrait', papertype='letter'
                                 , format='eps', transparent=False
                                 , bbox_inches='tight', pad_inches=0.1
                                 )
                fig3.savefig((FIGDIR3 + '/' + figurename_prefix + '.png'), dpi=None
                                 , facecolor='w', edgecolor='w'
                                 , orientation='portrait', papertype='letter'
                                 , format='png', transparent=False
                                 , bbox_inches='tight', pad_inches=0.1
                                 )
                plt.close(fig3)
            else:
                plt.show()
            # -----------------------------------------------
            
            #plotflag = False  # tmp option !!! PMB
        # -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
        
    #  END FOR LOOP OVER POLYORDER
#  END FOR LOOP OVER I_OBS