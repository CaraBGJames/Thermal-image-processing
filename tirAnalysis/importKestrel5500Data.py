#!/usr/bin/python

from datetime import datetime, date
from pytz import timezone, country_timezones

import pandas


names = ("datetime",
         "alt",
         "T_wb",
         "dir_mag",
         "dir_true",
         "p",
         "Crosswind",
         "T_dp",
         "T",
         "HSI",
         "alt_density",
         "Headwind",
         "RH",
         "wind",
         "P_bar",
         "T_wc")

cols2names = {'FORMATTED DATE_TIME':           'datetime',
              'Altitude':                      'alt',
              'Psychro Wet Bulb Temperature':  'T_wb',
              'Direction - True':              'dir_true',
              'Direction - Mag':               'dir_mag',
              'Station Pressure':              'p',
              'Crosswind':                     'Crosswind',
              'Dew Point':                     'T_dp',
              'Temperature':                   'T',
              'Heat Stress Index':             'HSI',
              'Density Altitude':              'alt_density',
              'Headwind':                      'Headwind',
              'Relative Humidity':             'RH',
              'Wind Speed':                    'wind',
              'Barometric Pressure':           'P_bar',
              'Wind Chill':                    'T_wc'}


def importKestrel5500Data(filename,
                          daterange=(None, None),
                          tzinfo={'localtz': 'utc', 'toUTC': False},
                          exportNames=False):
    '''
    A utility to automate the import of data logs from the Kestrel 5500 weather
    station.

    Parameters:
    -----------
    filename : str
        The name of the data log to be read
    daterange : tuple of datetime objects (optional)
        Start and stop times of the dataframe
    tzinfo: dict containing the fields
        localtz : str
            Two-letter code denoting the local timezone (default is 'utc'),
            though this is not recognised by 'country_timezones'
        toUTC : bool
            Convert timezone to UTC or not (default is False)
    exportNames : bool
        export the column names used as a tuple

    Returns:
    --------
    df : pandas.DataFrame
        A data frame containing the formatted contents of <filename>

    Structure & namespace:
    ----------------------
    The data logs have the following structure, abbreviated as:
             FORMATTED DATE_TIME    datetime    (used as index)
                        Altitude    alt
    Psychro Wet Bulb Temperature    T_wb
                 Direction - Mag    dir_mag
                Direction - True    dir_true
                Station Pressure    p
                       Crosswind    Crosswind
                       Dew Point    T_dp
                     Temperature    T
               Heat Stress Index    HSI
                Density Altitude    alt_density
                        Headwind    Headwind
               Relative Humidity    RH
                      Wind Speed    wind
             Barometric Pressure    P_bar
                      Wind Chill    T_wc
    '''
    # Try to parse the timezone information.  Defaults to UTC if the
    # 'localtz' parameter passed in tzinfo doesn't work out
    try:
        tz = timezone(country_timezones(tzinfo['localtz']))
    except KeyError:
        tz = timezone('UTC')

    df = pandas.read_csv(filename,
                         sep=',',                   # csv separated by commas
                         parse_dates=[0],           # FORMATTED DATE_TIME
                         skiprows=[0, 1, 2, 4],     # header is 12 rows long
                         na_values='***',           # ignore these
                         thousands=',')             # thousands separator

    df.rename(cols2names, axis='columns', inplace=True)
    df.index = df['datetime']
    df.drop(['datetime'], axis='columns', inplace=True)

    # Localise the dataframe according to the timezone
    df = df.tz_localize(tz)
    # Convert to UTC if requested
    if tzinfo['toUTC'] is not False:
        df = df.tz_localize('UTC')

    # Filter the dataframe by date, if necessary.  Note that the logical
    # slicing is done on the index
    if len(daterange) is not 2:
        raise ValueError('daterange must be of length 2')
    if daterange[0] is not None:
        df = df[df.index > daterange[0]]
    if daterange[1] is not None:
        df = df[df.index < daterange[1]]

    if exportNames:
        return df, names
    else:
        return df
