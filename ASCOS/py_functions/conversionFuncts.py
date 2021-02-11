"""
Functions for data conversion
==============================

"""

import numpy as np
# from __future__ import print_function

def reGrid_Sondes(data1, data2, data3, obs, doy, var):

    """
    Function to convert radiosonde and IFS data onto UM vertical grid
    """

    from scipy.interpolate import interp1d

    ### 6-hourly time binning for model
    ### um['time'][:24:6].data
    ###     BUT there is a problem since we have 25 timesteps (i.e. [24] == [25])
    ###     need to pick out where we have a repeated time value, then remove it so
    ###     that the time array can be indexed easily

    #### ---------------------------------------------------------------
    ### build list of variables names wrt input data [OBS, UM, CASIM, IFS]
    #### ---------------------------------------------------------------
    if var == 'temp':
        varlist = ['temperature','temperature','temperature','temperature']
    elif var == 'thetaE':
        # varlist = ['epottemp','thetaE','thetaE','thetaE']     # use sonde file's epottemp
        varlist = ['thetaE','thetaE','thetaE','thetaE']         # use sonde calculated thetaE
    elif var == 'theta':
        # varlist = ['pottemp','theta','theta','theta']     # use sonde file's pottemp
        varlist = ['theta','theta','theta','theta']         # use sonde calculated theta
    elif var == 'q':
        varlist = ['mr','q','q','q']

    ### stop double counting of 0000 and 2400 from model data
    temp = np.zeros([len(data1['time'])])
    for i in range(0, len(temp)-1):
        if data1['time'][i] == data1['time'][i+1]:
            continue
        else:
            temp[i] = data1['time'][i]
    ii = np.where(temp != 0.0)      ### picks out where data are non-zero

    #### ---------------------------------------------------------------
    #### save hourly temperature model profiles (using the ii index defined by the time indices)
    #### ---------------------------------------------------------------
    data1[var + '_hrly'] = np.squeeze(data1[varlist[1]][ii,:])
    data2[var + '_hrly'] = np.squeeze(data2[varlist[2]][ii,:])
    data3[var + '_hrly'] = np.squeeze(data3[varlist[3]][ii,:])

    #### ---------------------------------------------------------------
    #### explicitly save 6-hourly temperature model profiles and time binning for ease
    #### ---------------------------------------------------------------
    ### can use temp for all model data since they are on the same (hourly) time binning
    data1['time_6hrly'] = data1['time_hrly'][::6]
    data2['time_6hrly'] = data2['time_hrly'][::6]
    data3['time_6hrly'] = data3['time_hrly'][::6]
    data1[var + '_6hrly'] = data1[var + '_hrly'][::6]
    data2[var + '_6hrly'] = data2[var + '_hrly'][::6]
    data3[var + '_6hrly'] = data3[var + '_hrly'][::6]

    #### ---------------------------------------------------------------
    #### index to only look at altitudes <10km
    #### ---------------------------------------------------------------
    iTim = 0        ### initialised
    iObs = np.where(obs['sondes']['gpsaltitude'][:,iTim] <= 11000)
    iUM = np.where(data1['height'] <= 11000)
    iIFS = np.where(data3['height'][iTim,:] <= 11000)

    #### ---------------------------------------------------------------
    #### remove flagged IFS heights
    #### ---------------------------------------------------------------
    data3['height'][data3['height'] == -9999] = 0.0
            #### set all heights to zero if flagged. setting to nan caused problems
            ####        further on
    data3['height_hrly'] = np.squeeze(data3['height'][ii,:])  ### need to explicitly save since height coord changes at each timedump

    #### ---------------------------------------------------------------
    #### START INTERPOLATION
    #### ---------------------------------------------------------------
    print ('')
    print ('Defining IFS temperature profile as a function:')
    print ('using ifs.height[i,:] to define temperature profiles...')
    data3[var + '_hrly_UM'] = np.zeros([np.size(data3['time_hrly'],0),len(data1['height'][iUM[0][3:]])])
    for iTim in range(0,np.size(data3['time_hrly'],0)):
        # print (iTim)
        iIFSind = np.where(data3['height_hrly'][iTim,:] <= 11000)
        if np.all(data3['height_hrly'][iTim,:] == 0.0):
            data3[var + '_hrly_UM'][iTim,:] = np.nan
        else:
            fnct_IFS = interp1d(np.squeeze(data3['height_hrly'][iTim,iIFSind]), np.squeeze(data3[var + '_hrly'][iTim,iIFSind]))
            data3[var + '_hrly_UM'][iTim,:] = fnct_IFS(data1['height'][iUM[0][3:]].data)
    print ('...')
    print ('IFS(UM Grid) function worked!')
    print (var + ' IFS data now on UM vertical grid')
    print ('*****')
    ### assign for easier indexing later
    data3[var + '_6hrly_UM'] = data3[var + '_hrly_UM'][::6,:]
    data3[var + '_6hrly'] = data3[var + '_hrly'][::6,:]
    data3['height_6hrly'] = data3['height_hrly'][::6,:]  ### need to explicitly save since height coord changes at each timedump

    #### INTERPOLATION TESTING:
    # print (data3['temp_hrly_UM'].shape)
    # print (data3['time_hrly'][::6].shape)
    # print (data1['temp_hrly'][:,iUM[0][3:]].shape)
    # print (data1['time_hrly'][::6].shape)
    # for i in range(0, np.size(data3['temp_6hrly_UM'],0)):
    #     fig = plt.figure()
    #     plt.plot(data3['temp_6hrly_UM'][i,:],data1['height'][iUM[0][3:]], label = 'interpd')
    #     plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height_6hrly'][i,iIFS]), label = 'height indexed')
    #     plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height'][0,iIFS]), label = 'height0')
    #     plt.title('IFS test ' + str(data3['time_6hrly'][i]))
    #     plt.legend()
    #     plt.savefig('../FIGS/regrid/IFS_test_doy' + str(data3['time_6hrly'][i]) + '.png')
    #     if i == 0:
    #         plt.show()
    #     else:
    #         plt.close()

    print ('')
    print ('Defining Sonde temperature profile as a function for the UM:')
    obs['sondes'][var + '_allSondes_UM'] = np.zeros([np.size(obs['sondes']['doy'],0),len(data1['height'][iUM[0][3:]])])
    for iTim in range(0,np.size(obs['sondes']['doy'],0)):
        # print 'iTim = ', str(iTim)
        fnct_Obs = interp1d(np.squeeze(obs['sondes']['gpsaltitude'][iObs,iTim]), np.squeeze(obs['sondes'][varlist[0]][iObs,iTim]))
        obs['sondes'][var + '_allSondes_UM'][iTim,:] = fnct_Obs(data1['height'][iUM[0][3:]].data)
    print ('...')
    print ('Sonde(UM Grid) function worked!')
    print ('All ' + var + ' sonde data now on UM vertical grid.')
    print ('*****')
    #
    # print ('')
    # print ('Defining Sonde temperature profile as a function for the IFS:')
    # obs['sondes'][var + '_allSondes_IFS'] = np.zeros([np.size(obs['sondes']['doy'],0),len(data1['height'][0,iIFS])])
    # for iTim in range(0,np.size(obs['sondes']['doy'],0)):
    #     # print 'iTim = ', str(iTim)
    #     iIFS = np.where(data3['height'][iTim,:] <= 11000)
    #     fnct_ObsIFS = interp1d(np.squeeze(obs['sondes']['gpsaltitude'][iObs,iTim]), np.squeeze(obs['sondes'][varlist[0]][iObs,iTim]))
    #     obs['sondes'][var + '_allSondes_UM'][iTim,:] = fnct_ObsIFS(data3['height'][iTim,iIFS])
    # print ('...')
    # print ('Sonde(IFS Grid) function worked!')
    # print ('All ' + var + ' sonde data now on IFS_DATA vertical grid.')
    # print ('*****')

    #### ---------------------------------------------------------------
    #### ONLY LOOK AT SONDES FROM THE DRIFT
    #### ---------------------------------------------------------------
    drift = np.where(np.logical_and(obs['sondes']['doy'] >= 225.9, obs['sondes']['doy'] <= 258.0))

    ### save in dict for ease
    obs['sondes']['doy_drift'] = obs['sondes']['doy'][drift]
    obs['sondes']['drift'] = drift
    obs['sondes'][var + '_driftSondes_UM'] = obs['sondes'][var + '_allSondes_UM'][drift[0],:]

    #### INTERPOLATION TESTING - IFS + SONDE + UM_RA2M:
    # print (obs['sondes']['doy_drift'].shape)
    # print (obs['sondes']['temp_allSondes_UM'][drift[0],:].shape)
    # if var == 'temp':
    #     for i in range(0, np.size(obs['sondes']['doy_drift'])):
    #         plt.plot(np.squeeze(obs['sondes']['temperature'][iObs,drift[0][i]]) + 273.15,np.squeeze(obs['sondes']['gpsaltitude'][iObs,drift[0][i]]), '--', color = 'k', label = 'sonde-original')
    #         plt.plot(obs['sondes']['temp_driftSondes_UM'][i,:] + 273.15,data1['height'][iUM[0][3:]], color = 'k', label = 'sonde-interpd')
    #         plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height_6hrly'][i,iIFS]), '--', color = 'darkorange', label = 'ifs-Zindexed')
    #         plt.plot(data3['temp_6hrly_UM'][i,:],data1['height'][iUM[0][3:]], color = 'darkorange', label = 'ifs-interpd')
    #         plt.plot(data1['temp_6hrly'][i,iUM[0][3:]], data1['height'][iUM[0][3:]], color = 'steelblue', label = 'um_ra2m')
    #         plt.plot(data2['temp_6hrly'][i,iUM[0][3:]], data2['height'][iUM[0][3:]], color = 'forestgreen', label = 'um_casim-100')
    #         plt.title('REGRID test ' + str(np.round(obs['sondes']['doy_drift'][i],2)))
    #         plt.legend()
    #         plt.savefig('../FIGS/regrid/REGRID_Ttest_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
    #         if i == 0:
    #             plt.show()
    #         else:
    #             plt.close()
    # elif var == 'q':
    #     for i in range(0, np.size(obs['sondes']['doy_drift'])):
    #         plt.plot(np.squeeze(obs['sondes']['mr'][iObs,drift[0][i]]), np.squeeze(obs['sondes']['gpsaltitude'][iObs,drift[0][i]]), '--', color = 'k', label = 'sonde-original')
    #         plt.plot(obs['sondes'][var + '_driftSondes_UM'][i,:], data1['height'][iUM[0][3:]], color = 'k', label = 'sonde-interpd')
    #         plt.plot(np.squeeze(data3[var + '_6hrly'][i,iIFS])*1e3,np.squeeze(data3['height_6hrly'][i,iIFS]), '--', color = 'darkorange', label = 'ifs-Zindexed')
    #         plt.plot(data3[var + '_6hrly_UM'][i,:]*1e3,data1['height'][iUM[0][3:]], color = 'darkorange', label = 'ifs-interpd')
    #         plt.plot(data1[var + '_6hrly'][i,iUM[0][3:]]*1e3, data1['height'][iUM[0][3:]], color = 'steelblue', label = 'um_ra2m')
    #         plt.plot(data2[var + '_6hrly'][i,iUM[0][3:]]*1e3, data2['height'][iUM[0][3:]], color = 'forestgreen', label = 'um_casim-100')
    #         plt.title('REGRID test ' + str(np.round(obs['sondes']['doy_drift'][i],2)))
    #         plt.legend()
    #         plt.savefig('../FIGS/regrid/REGRID_Qtest_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
    #         if i == 0:
    #             plt.show()
    #         else:
    #             plt.close()

    #### ---------------------------------------------------------------
    #### make some dictionary assignments for use later
    #### ---------------------------------------------------------------
    data1['universal_height'] = data1['height'][iUM[0][3:]]
    data1['universal_height_UMindex'] = iUM[0][3:]

    return data1, data2, data3, obs, drift
