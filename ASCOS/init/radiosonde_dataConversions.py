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
from scipy.interpolate import interp1d

#### import python functions
import sys
sys.path.insert(1, '../py_functions/')
from time_functions import calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct, readMatlabData
from physFuncts import calcThetaE, calcThetaVL, adiabatic_lwc, calcTemperature
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

def LEM_LoadTHREF(data, sondenumber):

    '''
    Load initialisation reference potential temperature profile
        -- Data copied from gillian/LEM/tom_arc1/morr2712/NAMELIST/nmlFORCE_ASCOS1
    '''
    #### dz is 5m between 0 and 1157.5 m, then linearly interpolated to 10m between 1285m and 2395m
    data['ascos1'] = {}
    data['ascos1']['z'] = np.zeros(361)
    data['ascos1']['z'][1:233] = np.arange(2.5,1157.6,5.0)
    data['ascos1']['z'][233:249] = np.array([1162.5002, 1167.5026, 1172.5162, 1177.5726, 1182.7501,
                            1188.1971, 1194.1285, 1200.77, 1208.27, 1216.6285, 1225.6971, 1235.2501, 1245.0726,
                            1255.0162, 1265.0026, 1275.0002])
    data['ascos1']['z'][249:] = np.arange(1285., 2396., 10.0)

    data['ascos1']['thref'] = np.array([269.2,269.2029,269.2059,269.209,269.2124,269.216,
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
    data['q'] = data['sonde'][:,9]/1e3
    data['z'] = data['sonde'][:,6]

    data['theta'], data['thetaE'] = calcThetaE(data['temperature'], data['pressure'], data['q'])

    ####    --------------- FIGURE
    #
    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=MED_SIZE)
    # plt.rc('axes',labelsize=MED_SIZE)
    # plt.rc('xtick',labelsize=MED_SIZE)
    # plt.rc('ytick',labelsize=MED_SIZE)
    # plt.figure(figsize=(4,5))
    # plt.rc('legend',fontsize=MED_SIZE)
    # plt.subplots_adjust(top = 0.9, bottom = 0.12, right = 0.95, left = 0.25,
    #         hspace = 0.22, wspace = 0.4)
    #
    # yylim = 2.4e3
    #
    # plt.plot(data['ascos1']['thref'], data['ascos1']['z'], label = 'ASCOS1')
    # plt.plot(data['theta'], data['z'], label = 'SONDE')
    # plt.ylabel('Z [m]')
    # plt.xlabel('$\Theta$ [K]')
    # plt.ylim([0,yylim])
    # plt.xlim([265,295])
    #
    # plt.savefig('../FIGS/Quicklooks_LEM-ASCOS1_' + sondenumber + '.png')
    # plt.show()

    return data

def LEM_LoadTHINIT_QINIT1(data,sondenumber):


    '''
    Calculate initialisation potential temperature and moisture profiles
        -- LEM data in gillian/LEM/tom_arc1/morr2712/UPDATES/setprofileASCOS.f
    '''

    data['ascos1']['thinit'] = data['ascos1']['thref']
    data['ascos1']['qinit1'] = np.zeros(np.size(data['ascos1']['thinit']))

    sondei1 = np.where(data['z'] > 1.5e3)
    ascosi1 = np.where(data['ascos1']['z'] > 1.5e3)
    data['ascos1']['qinit1'][ascosi1] = data['q'][sondei1[0][0]]

    sondei2 = np.where(np.logical_and(data['z'] > 500, data['z'] < 1000))
    maxindex = np.where(data['q'][sondei2] == np.nanmax(data['q'][sondei2]))
    # print (data['q'][sondei2][maxindex])
    # print (data['z'][sondei2][maxindex])
        ### inversion of 0.00257 kg/kg at 654 m

    ascosi2 = np.where(np.logical_and(data['ascos1']['z'] > 650, data['ascos1']['z'] < 660))
    # print (data['ascos1']['z'][ascosi2[0][0]])

    data['ascos1']['qinit1'][ascosi2[0][0]] = data['q'][sondei2][maxindex]

    sondei3 = np.where(np.logical_and(data['z'] > data['z'][sondei2][maxindex], data['z'] < data['z'][sondei1[0][0]]))
    # print (data['z'][ascosi2[0][0]])
    # print (data['ascos1']['z'][ascosi1[0][0]])
    ascosi3 = np.where(np.logical_and(data['ascos1']['z'] > data['ascos1']['z'][ascosi2[0][0]], data['ascos1']['z'] < data['ascos1']['z'][ascosi1[0][0]]))

    # print ([data['z'][sondei3[0][0]],data['z'][sondei1[0][0]]])
    # print ([data['q'][sondei3[0][0]],data['q'][sondei1[0][0]]])
    # print (data['ascos1']['z'][int(ascosi3[0][0])+2:ascosi3[0][-1]])

    interp_qinit1 = interp1d([data['z'][sondei3[0][0]],data['z'][sondei1[0][0]]],[data['q'][sondei3[0][0]],data['q'][sondei1[0][0]]])
    q1temp = interp_qinit1(data['ascos1']['z'][int(ascosi3[0][0])+2:ascosi3[0][-1]])
    data['ascos1']['qinit1'][int(ascosi3[0][0])+2:ascosi3[0][-1]] = q1temp
    data['ascos1']['qinit1'][int(ascosi1[0][0]-1)] = data['q'][sondei1[0][0]]

    data['ascos1']['qinit1'][0:int(ascosi3[0][0])+2] = data['q'][sondei2][maxindex]

    ### build new height array for namelist initialisation
    nml_Z = np.arange(0., 700., 50.)
    nml_Z = np.append(nml_Z, np.arange(700., 2301., 100.))
    print (nml_Z)

    ### build qinit1 array
    interp_qinit1 = interp1d(data['ascos1']['z'],data['ascos1']['qinit1'])
    nml_qinit1 = interp_qinit1(nml_Z)

    ### build thref array
    interp_thref = interp1d(data['ascos1']['z'],data['ascos1']['thref'])
    nml_thref = interp_thref(nml_Z)

    ### manually append last value to 2400m (since last Z in ASCOS1 is 2395m and above interpolation range)
    nml_Z = np.append(nml_Z, 2400.)
    nml_qinit1 = np.append(nml_qinit1, data['ascos1']['qinit1'][-1])
    nml_thref = np.append(nml_thref, data['ascos1']['thref'][-1])

    ### save to dictionary so data can be easily passed to next function
    data['monc'] = {}
    data['monc']['z'] = nml_Z
    data['monc']['thref'] = nml_thref
    data['monc']['thinit'] = nml_thref
    data['monc']['qinit1'] = nml_qinit1


    ####    --------------- FIGURE

    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=MED_SIZE)
    # plt.rc('axes',labelsize=MED_SIZE)
    # plt.rc('xtick',labelsize=MED_SIZE)
    # plt.rc('ytick',labelsize=MED_SIZE)
    # plt.figure(figsize=(8,5))
    # plt.rc('legend',fontsize=MED_SIZE)
    # plt.subplots_adjust(top = 0.9, bottom = 0.12, right = 0.95, left = 0.1,
    #         hspace = 0.22, wspace = 0.4)
    #
    # yylim = 2.4e3

    # plt.subplot(121)
    # plt.plot(data['ascos1']['thinit'], data['ascos1']['z'], label = 'ASCOS1')
    # plt.plot(data['theta'], data['z'], label = 'SONDE')
    # plt.plot(nml_thref, nml_Z, '.', label = 'monc-namelist')
    # plt.ylabel('Z [m]')
    # plt.xlabel('$\Theta$ [K]')
    # plt.ylim([0,yylim])
    # plt.xlim([265,295])
    #
    # plt.subplot(122)
    # plt.plot(data['ascos1']['qinit1'][data['ascos1']['qinit1'] > 0], data['ascos1']['z'][data['ascos1']['qinit1'] > 0], label = 'ASCOS1')
    # plt.plot(data['q'], data['z'], label = 'SONDE')
    # plt.plot(nml_qinit1, nml_Z, '.', label = 'monc-namelist')
    # plt.xlabel('q [kg/kg]')
    # plt.grid('on')
    # plt.ylim([0,yylim])
    # plt.legend()
    #
    # plt.savefig('../FIGS/Quicklooks_LEM-ASCOS1-MONCnmlist_thinit-qinit1_' + sondenumber + '.png')
    # plt.show()

    ### print out to terminal in format for monc namelists
    print ('z_init_pl_theta = ')
    for line in nml_Z: sys.stdout.write('' + str(line).strip() + ',')
    print ('')

                # z_init_pl_theta = 0.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0,550.0,600.0,650.0,700.0,800.0,900.0,1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0,2100.0,2200.0,2300.0,2400.0

    print ('f_init_pl_theta = ')
    for line in nml_thref: sys.stdout.write('' + str(np.round(line,2)).strip() + ',')
    print ('')

                # f_init_pl_theta = 269.2,269.24,269.36,269.59,269.87,270.09,270.22,270.32,270.43,270.58,270.78,271.14,271.75,272.57,273.5,275.41,277.28,279.12,280.92,282.59,283.47,284.27,285.06,285.85,286.62,287.39,288.15,288.9,289.64,290.38,291.1,291.79

    print ('f_init_pl_q = ')
    for line in nml_qinit1: sys.stdout.write('' + str(np.round(line,5)).strip() + ',')
    print ('')

                # f_init_pl_q = 0.00257,0.00257,0.00257,0.00257,0.00257,0.00257,0.00257,0.00257,0.00257,0.00257,0.00257,0.00257,0.00257,0.00257,0.00246,0.00228,0.00211,0.00193,0.00176,0.00159,0.00141,0.00124,0.00105,0.00105,0.00105,0.00105,0.00105,0.00105,0.00105,0.00105,0.00105,0.00105


    ### namelist entries, saved for later:
    #names_init_pl_q=vapour
    #z_init_pl_q=0.,2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5,82.5,87.5,92.597.5,102.5,107.5,112.5,117.5,122.5,127.5,132.5,137.5,142.5,147.5,152.5,157.5,162.5,167.5,172.5,177.5,182.5,187.5,192.5,197.5,202.5,207.5,212.5,217.5,222.5,227.5,232.5,237.5,242.5,247.5,252.5,257.5,262.5,267.5,272.5,277.5,282.5,287.5,292.5,297.5,302.5,307.5,312.5,317.5,322.5,327.5,332.5,337.5,342.5,347.5,352.5,357.5,362.5,367.5,372.5,377.5,382.5,387.5,392.5,397.5,402.5,407.5,412.5,417.5,422.5,427.5,432.5,437.5,442.5,447.5,452.5,457.5,462.5,467.5,472.5,477.5,482.5,487.5,492.5,497.5,502.5,507.5,512.5,517.5,522.5,527.5,532.5,537.5,542.5,547.5,552.5,557.5,562.5,567.5,572.5,577.5,582.5,587.5,592.5,597.5,602.5,607.5,612.5,617.5,622.5,627.5,632.5,637.5,642.5,647.5,652.5,657.5,662.5,667.5,672.5,677.5,682.5,687.5,692.5,697.5,702.5,707.5,712.5,717.5,722.5,727.5,732.5,737.5,742.5,747.5,752.5,757.5,762.5,767.5,772.5,777.5,782.5,787.5,792.5,797.5,802.5,807.5,812.5,817.5,822.5,827.5,832.5,837.5,842.5,847.5,852.5,857.5,862.5,867.5,872.5,877.5,882.5,887.5,892.5,897.5,902.5,907.5,912.5,917.5,922.5,927.5,932.5,937.5,942.5,947.5,952.5,957.5,962.5,967.5,972.5,977.5,982.5,987.5,992.5,997.5,1002.5,1007.5,1012.5,1017.5,1022.5,1027.5,1032.5,1037.5,1042.5,1047.5,1052.5,1057.5,1062.5,1067.5,1072.5,1077.5,1082.5,1087.5,1092.5,1097.5,1102.5,1107.5,1112.5,1117.5,1122.5,1127.5,1132.5,1137.5,1142.5,1147.5,1152.5,1157.5,1162.5002,1167.5026,1172.5162,1177.5726,1182.7501,1188.1971,1194.1285,1200.77,1208.27,1216.6285,1225.6971,1235.2501,1245.0726,1255.0162,1265.0026,1275.0002,1285.,1295.,1305.,1315.,1325.,1335.,1345.,1355.,1365.,1375.,1385.,1395.,1405.,1415.,1425.,1435.,1445.,1455.,1465.,1475.,1485.,1495.,1505.,1515.,1525.,1535.,1545.,1555.,1565.,1575.,1585.,1595.,1605.,1615.,1625.,1635.,1645.,1655.,1665.,1675.,1685.,1695.,1705.,1715.,1725.,1735.,1745.,1755.,1765.,1775.,1785.,1795.,1805.,1815.,1825.,1835.,1845.,1855.,1865.,1875.,1885.,1895.,1905.,1915.,1925.,1935.,1945.,1955.,1965.,1975.,1985.,1995.,2005.,2015.,2025.,2035.,2045.,2055.,2065.,2075.,2085.,2095.,2105.,2115.,2125.,2135.,2145.,2155.,2165.,2175.,2185.,2195.,2205.,2215.,2225.,2235.,2245.,2255.,2265.,2275.,2285.,2295.,2305.,2315.,2325.,2335.,2345.,2355.,2365.,2375.,2385.,2395.

    #z_init_pl_theta=0.,2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5,82.5,87.5,92.597.5,102.5,107.5,112.5,117.5,122.5,127.5,132.5,137.5,142.5,147.5,152.5,157.5,162.5,167.5,172.5,177.5,182.5,187.5,192.5,197.5,202.5,207.5,212.5,217.5,222.5,227.5,232.5,237.5,242.5,247.5,252.5,257.5,262.5,267.5,272.5,277.5,282.5,287.5,292.5,297.5,302.5,307.5,312.5,317.5,322.5,327.5,332.5,337.5,342.5,347.5,352.5,357.5,362.5,367.5,372.5,377.5,382.5,387.5,392.5,397.5,402.5,407.5,412.5,417.5,422.5,427.5,432.5,437.5,442.5,447.5,452.5,457.5,462.5,467.5,472.5,477.5,482.5,487.5,492.5,497.5,502.5,507.5,512.5,517.5,522.5,527.5,532.5,537.5,542.5,547.5,552.5,557.5,562.5,567.5,572.5,577.5,582.5,587.5,592.5,597.5,602.5,607.5,612.5,617.5,622.5,627.5,632.5,637.5,642.5,647.5,652.5,657.5,662.5,667.5,672.5,677.5,682.5,687.5,692.5,697.5,702.5,707.5,712.5,717.5,722.5,727.5,732.5,737.5,742.5,747.5,752.5,757.5,762.5,767.5,772.5,777.5,782.5,787.5,792.5,797.5,802.5,807.5,812.5,817.5,822.5,827.5,832.5,837.5,842.5,847.5,852.5,857.5,862.5,867.5,872.5,877.5,882.5,887.5,892.5,897.5,902.5,907.5,912.5,917.5,922.5,927.5,932.5,937.5,942.5,947.5,952.5,957.5,962.5,967.5,972.5,977.5,982.5,987.5,992.5,997.5,1002.5,1007.5,1012.5,1017.5,1022.5,1027.5,1032.5,1037.5,1042.5,1047.5,1052.5,1057.5,1062.5,1067.5,1072.5,1077.5,1082.5,1087.5,1092.5,1097.5,1102.5,1107.5,1112.5,1117.5,1122.5,1127.5,1132.5,1137.5,1142.5,1147.5,1152.5,1157.5,1162.5002,1167.5026,1172.5162,1177.5726,1182.7501,1188.1971,1194.1285,1200.77,1208.27,1216.6285,1225.6971,1235.2501,1245.0726,1255.0162,1265.0026,1275.0002,1285.,1295.,1305.,1315.,1325.,1335.,1345.,1355.,1365.,1375.,1385.,1395.,1405.,1415.,1425.,1435.,1445.,1455.,1465.,1475.,1485.,1495.,1505.,1515.,1525.,1535.,1545.,1555.,1565.,1575.,1585.,1595.,1605.,1615.,1625.,1635.,1645.,1655.,1665.,1675.,1685.,1695.,1705.,1715.,1725.,1735.,1745.,1755.,1765.,1775.,1785.,1795.,1805.,1815.,1825.,1835.,1845.,1855.,1865.,1875.,1885.,1895.,1905.,1915.,1925.,1935.,1945.,1955.,1965.,1975.,1985.,1995.,2005.,2015.,2025.,2035.,2045.,2055.,2065.,2075.,2085.,2095.,2105.,2115.,2125.,2135.,2145.,2155.,2165.,2175.,2185.,2195.,2205.,2215.,2225.,2235.,2245.,2255.,2265.,2275.,2285.,2295.,2305.,2315.,2325.,2335.,2345.,2355.,2365.,2375.,2385.,2395.
    #f_init_pl_theta=269.2,269.2029,269.2059,269.209,269.2124,269.216,269.22,269.2244,269.2293,269.2348,269.241,269.2479,269.2556,269.2642,269.2736,269.2841,269.2956,269.3083,269.322,269.3369,269.3531,269.3704,269.3889,269.4086,269.4295,269.4514,269.4745,269.4986,269.5235,269.5493,269.5759,269.603,269.6307,269.6587,269.687,269.7154,269.7437,269.7719,269.7998,269.8273,269.8543,269.8807,269.9063,269.9311,269.9551,269.9781,270.0002,270.0212,270.0413,270.0603,270.0784,270.0955,270.1117,270.127,270.1414,270.1551,270.1681,270.1805,270.1923,270.2036,270.2145,270.225,270.2353,270.2454,270.2553,270.2652,270.275,270.2848,270.2947,270.3047,270.3148,270.325,270.3354,270.346,270.3569,270.3679,270.3792,270.3907,270.4025,270.4145,270.4267,270.4393,270.4521,270.4652,270.4787,270.4924,270.5065,270.521,270.5359,270.5513,270.5672,270.5837,270.6007,270.6185,270.6371,270.6565,270.6768,270.6982,270.7208,270.7446,270.7697,270.7963,270.8246,270.8545,270.8863,270.9199,270.9557,270.9936,271.0337,271.0762,271.1211,271.1684,271.2183,271.2708,271.3259,271.3835,271.4438,271.5067,271.5721,271.64,271.7104,271.7831,271.858,271.9351,272.0143,272.0953,272.1782,272.2628,272.3489,272.4364,272.5251,272.6151,272.706,272.7978,272.8905,272.9837,273.0776,273.1719,273.2666,273.3617,273.4569,273.5524,273.6479,273.7436,273.8393,273.935,274.0306,274.1262,274.2218,274.3173,274.4127,274.508,274.6032,274.6983,274.7933,274.8881,274.9829,275.0776,275.1722,275.2667,275.3611,275.4554,275.5496,275.6437,275.7377,275.8317,275.9255,276.0193,276.113,276.2067,276.3002,276.3937,276.4871,276.5804,276.6737,276.7668,276.8599,276.9529,277.0459,277.1388,277.2316,277.3243,277.417,277.5096,277.6021,277.6945,277.7869,277.8792,277.9714,278.0636,278.1556,278.2477,278.3396,278.4315,278.5232,278.615,278.7066,278.7982,278.8897,278.9811,279.0724,279.1637,279.2549,279.346,279.437,279.5279,279.6187,279.7095,279.8001,279.8907,279.9811,280.0714,280.1616,280.2517,280.3417,280.4315,280.5211,280.6106,280.7,280.7891,280.8781,280.9668,281.0554,281.1437,281.2318,281.3197,281.4073,281.4947,281.5817,281.6685,281.755,281.8412,281.9272,282.0128,282.0981,282.1831,282.2678,282.3522,282.4363,282.5201,282.6036,282.6868,282.7698,282.8524,282.9349,283.0171,283.099,283.1808,283.2623,283.3436,283.4248,283.5058,283.5866,283.6673,283.7478,283.8282,283.9085,283.9886,284.0687,284.1486,284.2285,284.3083,284.3879,284.4675,284.547,284.6265,284.7058,284.7851,284.8643,284.9434,285.0224,285.1014,285.1803,285.2591,285.3378,285.4164,285.495,285.5735,285.6519,285.7302,285.8084,285.8865,285.9645,286.0425,286.1204,286.1982,286.2758,286.3534,286.431,286.5084,286.5857,286.663,286.7401,286.8172,286.8942,286.9711,287.0479,287.1246,287.2012,287.2777,287.3542,287.4305,287.5068,287.583,287.6591,287.7351,287.811,287.8868,287.9626,288.0382,288.1138,288.1893,288.2646,288.34,288.4152,288.4903,288.5653,288.6403,288.7152,288.7899,288.8646,288.9392,289.0138,289.0882,289.1626,289.2368,289.311,289.3851,289.4591,289.5331,289.6069,289.6807,289.7544,289.828,289.9015,289.9749,290.0483,290.1216,290.1948,290.2679,290.341,290.414,290.4869,290.5598,290.6326,290.7053,290.778,290.8506,290.9232,290.9957,291.0682,291.1406,291.213,291.2853,291.3576,291.4299,291.5022,291.5744,291.6466,291.7189,291.7911

    return data

def LEM_LoadQINIT2(data, sondenumber):

    '''
    Calculate adiabatic lwc up to 650m (main inversion):
        -- take thref and pressure, & calculate temperature
        # lwc_adiabatic(ii,liquid_bases(jj):liquid_tops(jj)) = dlwc_dz.*[1:(liquid_tops(jj)-liquid_bases(jj)+1)].*dheight;
    '''

    print (data['z'])
    print (data['monc']['z'][1:])

    data['monc']['pressure'] = np.zeros(np.size(data['monc']['z']))
    interp_pres = interp1d(data['z'],data['pressure'])
    data['monc']['pressure'][1:] = interp_pres(data['monc']['z'][1:])
    data['monc']['pressure'][0] = 102050. ## reference surface pressure from mcf

    ### calculate temperature from thref and pressure
    temp_T = calcTemperature(data['monc']['thref'], data['monc']['pressure'])

    ### adapt temperature array
    data['monc']['temperature'] = np.zeros(np.size(data['monc']['z']))
    data['monc']['temperature'][0] = temp_T[0]
    data['monc']['temperature'][1] = temp_T[1] - 0.3
    data['monc']['temperature'][2] = temp_T[2] - 0.5
    data['monc']['temperature'][3] = temp_T[3] - 0.7
    data['monc']['temperature'][4] = temp_T[4] - 0.8
    data['monc']['temperature'][5] = temp_T[5] - 0.3
    data['monc']['temperature'][6:11] = temp_T[6:11]
    data['monc']['temperature'][11] = temp_T[11] - 0.25
    data['monc']['temperature'][12] = temp_T[12] - 0.6

    ### calculate qinit2
    ### interpolate free troposphere temperatures from radiosonde onto monc namelist gridding
    interp_temp = interp1d(data['z'], data['temperature'])
    data['monc']['temperature'][13:] = interp_temp(data['monc']['z'][13:])

    ### calculate adiabatic lwc rate of change
    dlwcdz, dqldz, dqdp = adiabatic_lwc(data['monc']['temperature'], data['monc']['pressure'])
    dheight = data['monc']['z'][1:] - data['monc']['z'][:-1]

    print (dheight)
    print (dlwcdz)

    freetrop_index = np.where(data['monc']['z'] > 600.0)
    dheight[int(freetrop_index[0][0]):] = 0.0   ## ignore points in the free troposphere
    blcloud_index = np.where(data['monc']['z'] < 300.0)
    dheight[blcloud_index] = 0.0   ## ignore points in the free troposphere

    data['monc']['qinit2'] = dlwcdz[:-1] * dheight
    data['monc']['qinit2'] = np.append(data['monc']['qinit2'], 0.)

    ### adapt theta (thinit and thref) based on revised temperature profile
    data['monc']['thinit'], thetaE = calcThetaE(data['monc']['temperature'], data['monc']['pressure'], data['monc']['qinit1'])
    data['monc']['thref'] = data['monc']['thinit']

    ####    --------------- FIGURE

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.figure(figsize=(13,5))
    plt.rc('legend',fontsize=MED_SIZE)
    plt.subplots_adjust(top = 0.9, bottom = 0.12, right = 0.95, left = 0.1,
            hspace = 0.22, wspace = 0.4)

    yylim = 2.4e3

    plt.subplot(151)
    plt.plot(data['ascos1']['thinit'], data['ascos1']['z'], color = 'steelblue', label = 'ASCOS1')
    plt.plot(data['theta'], data['z'], color = 'darkorange', label = 'SONDE')
    plt.plot(data['monc']['thref'], data['monc']['z'], 'k.', label = 'monc-namelist')
    plt.ylabel('Z [m]')
    plt.xlabel('$\Theta$ [K]')
    plt.grid('on')
    plt.ylim([0,yylim])
    plt.xlim([265,295])

    plt.subplot(152)
    plt.plot(data['ascos1']['qinit1'][data['ascos1']['qinit1'] > 0]*1e3, data['ascos1']['z'][data['ascos1']['qinit1'] > 0], color = 'steelblue', label = 'LEM-ASCOS1')
    plt.plot(data['q']*1e3, data['z'], color = 'darkorange', label = 'SONDE')
    plt.plot(data['monc']['qinit1']*1e3, data['monc']['z'], 'k.', label = 'monc-namelist')
    plt.xlabel('qinit1 [g/kg]')
    plt.grid('on')
    plt.ylim([0,yylim])
    plt.legend(bbox_to_anchor=(0.25, 1.01, 1., .102), loc=3, ncol=3)

    plt.subplot(153)
    # plt.plot(data['ascos1']['qinit1'][data['ascos1']['qinit1'] > 0], data['ascos1']['z'][data['ascos1']['qinit1'] > 0], color = 'steelblue', label = 'ASCOS1')
    plt.plot(data['pressure'], data['z'], color = 'darkorange', label = 'SONDE')
    plt.plot(data['monc']['pressure'], data['monc']['z'], 'k.', label = 'monc-namelist')
    plt.xlabel('pressure [Pa]')
    plt.grid('on')
    plt.ylim([0,yylim])
    plt.xlim([7e4, 10.5e4])

    plt.subplot(154)
    # plt.plot(data['ascos1']['qinit1'][data['ascos1']['qinit1'] > 0], data['ascos1']['z'][data['ascos1']['qinit1'] > 0], color = 'steelblue', label = 'ASCOS1')
    plt.plot(data['temperature'], data['z'], color = 'darkorange', label = 'SONDE')
    plt.plot(data['monc']['temperature'], data['monc']['z'], 'k.', label = 'monc-namelist')
    plt.xlabel('temperature [K]')
    plt.grid('on')
    plt.ylim([0,yylim])
    plt.xlim([265,275])

    plt.subplot(155)
    plt.plot(data['monc']['qinit2']*1e3, data['monc']['z'][:], 'k.', label = 'monc-namelist')
    plt.xlabel('qinit2 [g/kg]')
    plt.grid('on')
    plt.ylim([0,yylim])
    # plt.xlim([265,275])

    # plt.savefig('../FIGS/Quicklooks_LEM-ASCOS1-MONCnmlist_thinit-qinit1_' + sondenumber + '.png')
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
    # figure = quicklooksSonde(data, sondenumber)

    ## -------------------------------------------------------------
    ## Read in data from LEM namelists
    ## -------------------------------------------------------------
    data = LEM_LoadTHREF(data, sondenumber)
    data = LEM_LoadTHINIT_QINIT1(data, sondenumber)
    data = LEM_LoadQINIT2(data, sondenumber)

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
