"""
Functions for loading specific datasets
==============================

"""

import numpy as np
# from __future__ import print_function

def load_radar(proj, day):

    """
    Function to load in radar data

    Written to allow different load in procedure for different projects

    Takes:
     'proj' argument, options:
        1. moccha
         'day' argument, options:
            1. date in format YYYYMMDD (string)
            2. all (string), loads all data (20180814 to 20180914)

    Use example:
    data = load_radar('moccha', 'all')

    """

    from netCDF4 import Dataset
    from time_functions import calcTime_Date2DOY

    if proj == 'moccha':

        ### choose variables to load
        var_list = ['time','range','Zh']

        ### make empty data dictionary
        data = {}
        temp = {}

        if day != 'all':

            ### define filename
            data_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
            filename = data_dir + 'mmcr/' + day + '_oden_mira.nc'

            ### load file
            nc = Dataset(filename,'r')

            ### populate dictionary
            for var in var_list:
                data[var] = nc.variables[var][:]

                ### transpose 2D data so that time is 0th axis
                if np.ndim(data[var]) == 2:
                    data[var] = np.transpose(data[var])

            ### find date in DOY format
            date = calcTime_Date2DOY(day)

            ### update time array to reference base date
            data['time'] = data['time']/24.0 + date

            ### remove flagged data points
            data['Zh'][data['Zh'] == -999.0] = np.nan

            nc.close()

        else:

            ### define filename
            data_dir = '/home/gillian/MOCCHA/ODEN/DATA/'

            moccha_names = ['20180814_oden_','20180815_oden_','20180816_oden_',
            '20180817_oden_','20180818_oden_','20180819_oden_','20180820_oden_',
            '20180821_oden_','20180822_oden_','20180823_oden_','20180824_oden_',
            '20180825_oden_','20180826_oden_','20180827_oden_','20180828_oden_',
            '20180829_oden_','20180830_oden_','20180831_oden_','20180901_oden_',
            '20180902_oden_','20180903_oden_','20180904_oden_','20180905_oden_',
            '20180906_oden_','20180907_oden_','20180908_oden_','20180909_oden_',
            '20180911_oden_','20180912_oden_','20180913_oden_','20180914_oden_']

            for name in moccha_names:

                filename = data_dir + 'mmcr/' + name + 'mira.nc'

                ### load file
                nc = Dataset(filename,'r')

                ### initialise array with first file
                if name == '20180814_oden_':
                    ### populate dictionary
                    for var in var_list:
                        data[var] = nc.variables[var][:]

                        ### transpose 2D data so that time is 0th axis
                        if np.ndim(data[var]) == 2:
                            data[var] = np.transpose(data[var])

                    ### find date in DOY format
                    date = calcTime_Date2DOY(name[:8])

                    ### update time array to reference base date
                    data['time'] = data['time']/24.0 + date

                ### subsequent dates are appended
                else:
                    ### populate dictionary
                    for var in var_list:
                        temp[var] = nc.variables[var][:]

                        ### transpose 2D data so that time is 0th axis
                        if np.ndim(temp[var]) == 2:
                            temp[var] = np.transpose(temp[var])

                    ### find date in DOY format
                    date = calcTime_Date2DOY(name[:8])

                    ### update time array to reference base date
                    temp['time'] = temp['time']/24.0 + date

                    ### append loaded data to existing data dictionary
                    for var in var_list:
                        if var == 'range':
                            continue
                        elif np.ndim(data[var]) == 1:
                            data[var] = np.append(data[var], temp[var])
                        elif np.ndim(data[var]) == 2:
                            data[var] = np.append(data[var], temp[var], 0)

                nc.close()

            ### remove flagged data points
            data['Zh'][data['Zh'] == -999.0] = np.nan

    else:
        print('*** invalid project name ***')

    return data
