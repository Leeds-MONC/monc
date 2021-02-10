###
###
### SCRIPT TO READ IN ASCOS RADIOSONDE DATA AND OUTPUT FOR MONC
###
###

from __future__ import print_function
import time
import datetime
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import os
import seaborn as sns

#### import python functions
import sys
sys.path.insert(1, '../py_functions/')
from time_functions import calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct, readMatlabData
from physFuncts import calcThetaE, calcThetaVL
from pyFixes import py3_FixNPLoad

def quicklooksSonde(data, sondenumber):

    '''
    Notes from Tom's readme:
    Each sounding is in the Xdate_time_EDT variables. The data is organized such that column 7
    is the altitude, 3 is temperature, 4 is relative humidity, 10 is specific humidity, 12 is
    wind speed, and the rest you can probably figure out by looking at the variables. Some of
    them makes no sense at all; this is a quick monitiring file generated by the sounding
    system.

    There is also a variable named "sond" which is interpolated to same heights (in height). In
    that is also 13=potential temperature and 14=equivalent potential temperature. Time is the
    nominal times (00, 06, 12 etc UTC). The actual times are in the log file.
    '''

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.figure(figsize=(12,5))
    plt.rc('legend',fontsize=MED_SIZE)
    plt.subplots_adjust(top = 0.9, bottom = 0.12, right = 0.92, left = 0.1,
            hspace = 0.22, wspace = 0.4)

    yylim = 2.4e3

    plt.subplot(131)
    plt.plot(data['sonde'][:,2], data['sonde'][:,6])
    plt.ylabel('Z [m]')
    plt.xlabel('Temperature [K]')
    plt.ylim([0,yylim])
    plt.xlim([265,276])

    plt.subplot(132)
    plt.plot(data['sonde'][:,9], data['sonde'][:,6])
    plt.xlabel('Spec. Hum. [g/kg]')
    plt.ylim([0,yylim])

    plt.subplot(133)
    plt.plot(data['sonde'][:,3], data['sonde'][:,6])
    plt.xlabel('Rel. Hum. [%]')
    plt.ylim([0,yylim])

    plt.savefig('../FIGS/Quicklooks_' + sondenumber + '.png')
    plt.show()


def LEM_Load(data, sondenumber):

    '''
    Data copied from gillian/LEM/tom_arc1/morr2712/NAMELIST/nmlFORCE_ASCOS1
    '''
    #### dz is 5m between 0 and 1157.5 m, then linearly interpolated to 10m between 1285m and 2395m
    data['ascos1'] = {}
    data['ascos1']['z'] = np.zeros(361)
    data['ascos1']['z'][1:233] = np.arange(2.5,1157.6,5.0)
    data['ascos1']['z'][233:249] = np.array([1162.5002, 1167.5026, 1172.5162, 1177.5726, 1182.7501,
                            1188.1971, 1194.1285, 1200.77, 1208.27, 1216.6285, 1225.6971, 1235.2501, 1245.0726,
                            1255.0162, 1265.0026, 1275.0002])
    data['ascos1']['z'][249:] = np.arange(1285., 2396., 10.0)

    data['ascos1']['th'] = np.array([269.2,269.2029,269.2059,269.209,269.2124,269.216,
                            269.22,269.2244,269.2293,269.2348,269.241,269.2479,269.2556,
                            269.2642,269.2736,269.2841,269.2956,269.3083,269.322,269.3369,
                            269.3531,269.3704,269.3889,269.4086,269.4295,269.4514,269.4745,
                            269.4986,269.5235,269.5493,269.5759,269.603,269.6307,269.6587,
                            269.687,269.7154,269.7437,269.7719,269.7998,269.8273,269.8543,
                            269.8807,269.9063,269.9311,269.9551,269.9781,270.0002,270.0212,
                            270.0413,270.0603,270.0784,270.0955,270.1117,270.127,270.1414,
                            270.1551,270.1681,270.1805,270.1923,270.2036,270.2145,270.225,
                            270.2353,270.2454,270.2553,270.2652,270.275,270.2848,270.2947,
                            270.3047,270.3148,270.325,270.3354,270.346,270.3569,270.3679,
                            270.3792,270.3907,270.4025,270.4145,270.4267,270.4393,270.4521,
                            270.4652,270.4787,270.4924,270.5065,270.521,270.5359,270.5513,
                            270.5672,270.5837,270.6007,270.6185,270.6371,270.6565,270.6768,
                            270.6982,270.7208,270.7446,270.7697,270.7963,270.8246,270.8545,
                            270.8863,270.9199,270.9557,270.9936,271.0337,271.0762,271.1211,
                            271.1684,271.2183,271.2708,271.3259,271.3835,271.4438,271.5067,
                            271.5721,271.64,271.7104,271.7831,271.858,271.9351,272.0143,
                            272.0953,272.1782,272.2628,272.3489,272.4364,272.5251,272.6151,
                            272.706,272.7978,272.8905,272.9837,273.0776,273.1719,273.2666,
                            273.3617,273.4569,273.5524,273.6479,273.7436,273.8393,273.935,
                            274.0306,274.1262,274.2218,274.3173,274.4127,274.508,274.6032,
                            274.6983,274.7933,274.8881,274.9829,275.0776,275.1722,275.2667,
                            275.3611,275.4554,275.5496,275.6437,275.7377,275.8317,275.9255,
                            276.0193,276.113,276.2067,276.3002,276.3937,276.4871,276.5804,
                            276.6737,276.7668,276.8599,276.9529,277.0459,277.1388,277.2316,
                            277.3243,277.417,277.5096,277.6021,277.6945,277.7869,277.8792,
                            277.9714,278.0636,278.1556,278.2477,278.3396,278.4315,278.5232,
                            278.615,278.7066,278.7982,278.8897,278.9811,279.0724,279.1637,
                            279.2549,279.346,279.437,279.5279,279.6187,279.7095,279.8001,
                            279.8907,279.9811,280.0714,280.1616,280.2517,280.3417,280.4315,
                            280.5211,280.6106,280.7,280.7891,280.8781,280.9668,281.0554,281.1437,
                            281.2318,281.3197,281.4073,281.4947,281.5817,281.6685,281.755,
                            281.8412,281.9272,282.0128,282.0981,282.1831,282.2678,282.3522,
                            282.4363,282.5201,282.6036,282.6868,282.7698,282.8524,282.9349,
                            283.0171,283.099,283.1808,283.2623,283.3436,283.4248,283.5058,
                            283.5866,283.6673,283.7478,283.8282,283.9085,283.9886,284.0687,
                            284.1486,284.2285,284.3083,284.3879,284.4675,284.547,284.6265,
                            284.7058,284.7851,284.8643,284.9434,285.0224,285.1014,285.1803,
                            285.2591,285.3378,285.4164,285.495,285.5735,285.6519,285.7302,
                            285.8084,285.8865,285.9645,286.0425,286.1204,286.1982,286.2758,
                            286.3534,286.431,286.5084,286.5857,286.663,286.7401,286.8172,
                            286.8942,286.9711,287.0479,287.1246,287.2012,287.2777,287.3542,
                            287.4305,287.5068,287.583,287.6591,287.7351,287.811,287.8868,
                            287.9626,288.0382,288.1138,288.1893,288.2646,288.34,288.4152,
                            288.4903,288.5653,288.6403,288.7152,288.7899,288.8646,288.9392,
                            289.0138,289.0882,289.1626,289.2368,289.311,289.3851,289.4591,
                            289.5331,289.6069,289.6807,289.7544,289.828,289.9015,289.9749,
                            290.0483,290.1216,290.1948,290.2679,290.341,290.414,290.4869,
                            290.5598,290.6326,290.7053,290.778,290.8506,290.9232,290.9957,
                            291.0682,291.1406,291.213,291.2853,291.3576,291.4299,291.5022,
                            291.5744,291.6466,291.7189,291.7911])

    data['pressure'] = data['sonde'][:,7]*1e2
    data['temperature'] = data['sonde'][:,2]
    data['q'] = data['sonde'][:,9]
    data['z'] = data['sonde'][:,6]

    data['theta'], data['thetaE'] = calcThetaE(data['temperature'], data['pressure'], data['q'])


    ####    --------------- FIGURE

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.figure(figsize=(4,5))
    plt.rc('legend',fontsize=MED_SIZE)
    plt.subplots_adjust(top = 0.9, bottom = 0.12, right = 0.95, left = 0.25,
            hspace = 0.22, wspace = 0.4)

    yylim = 2.4e3

    plt.plot(data['ascos1']['th'], data['ascos1']['z'], label = 'ASCOS1')
    plt.plot(data['theta'], data['z'], label = 'SONDE')
    plt.ylabel('Z [m]')
    plt.xlabel('$\Theta$ [K]')
    plt.ylim([0,yylim])
    plt.xlim([265,295])

    plt.savefig('../FIGS/Quicklooks_LEM-ASCOS1_' + sondenumber + '.png')
    plt.show()


def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    '''
    Python script to build initialisation data for MONC from radiosondes
    '''

    print ('Import ASCOS radiosonde data:')
    print ('...')

    print ('Load radiosonde data from Jutta...')
    sondes = readMatlabData('../DATA/radiosondes.mat')

    print ('')
    # print (sondes.keys())

    ## -------------------------------------------------------------
    ## Load radiosonde from 20180827 1200UTC (as in Tom's work)
    ## -------------------------------------------------------------
    data = {}
    sondenumber = 'X080827_12_EDT'
    data['sonde'] = sondes[sondenumber]

    ## -------------------------------------------------------------
    ## Quicklook plots of chosen sonde
    ## -------------------------------------------------------------
    figure = quicklooksSonde(data, sondenumber)

    ## -------------------------------------------------------------
    ## Read in data from LEM namelists
    ## -------------------------------------------------------------
    figure = LEM_Load(data, sondenumber)

    ## -------------------------------------------------------------
    ## save out working data for testing
    ## -------------------------------------------------------------
    np.save('working_data', data)

    # -------------------------------------------------------------
    # FIN.
    # -------------------------------------------------------------
    END_TIME = time.time()
    print ('******')
    print ('')
    print ('End: ' + time.strftime("%c"))
    print ('')


if __name__ == '__main__':

    main()
