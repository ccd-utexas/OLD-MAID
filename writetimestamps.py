#!/usr/bin/env python

# Imports.
# Standard libraries.
from __future__ import absolute_import, division, print_function
import os
import csv
import sys
import datetime as dt
import dateutil.parser
# Installed packages.
import pandas as pd
from bs4 import BeautifulSoup
# Local modules.
import read_spe

def main(fpath_spe):    
    # Create a dataframe of timestamps.
    # Check the trigger type, timestamp for beginning data acquisition, and ticks per second.
    # View the data frame.
    fpath_spe  = os.path.abspath(fpath_spe)
    spe = read_spe.File(fpath_spe)
    
    if hasattr(spe, 'footer_metadata'):
        footer_metadata = BeautifulSoup(spe.footer_metadata, "xml")
        trigger_response = footer_metadata.find(name='TriggerResponse').text
        ts_begin = footer_metadata.find(name='TimeStamp', event='ExposureStarted').attrs['absoluteTime']
        dt_begin = dateutil.parser.parse(ts_begin)
        ticks_per_second = int(footer_metadata.find(name='TimeStamp', event='ExposureStarted').attrs['resolution'])
    else:
        print(("WARNING: No XML footer metadata.\n" +
               "Unknown trigger response.\n" +
               "Using file creation time as absolute timestamp.\n" +
               "Assuming 1E6 ticks per seconds."), file=sys.stderr)
        trigger_response = ""
        dt_begin = dt.datetime.utcfromtimestamp(os.path.getctime(fpath_spe))
        ticks_per_second = 1E6
    idx_metadata_map = {}
    for idx in xrange(spe.get_num_frames()):
        (frame, metadata) = spe.get_frame(idx)
        idx_metadata_map[idx] = metadata
    spe.close()
    df_metadata = pd.DataFrame.from_dict(idx_metadata_map, orient='index')
    df_metadata = df_metadata.set_index(keys='frame_tracking_number')
    df_metadata = df_metadata[['time_stamp_exposure_started', 'time_stamp_exposure_ended']].applymap(lambda x: x / ticks_per_second)
    df_metadata = df_metadata[['time_stamp_exposure_started', 'time_stamp_exposure_ended']].applymap(lambda x : dt_begin + dt.timedelta(seconds=x))
    df_metadata[['diff_time_stamp_exposure_started', 'diff_time_stamp_exposure_ended']] = df_metadata - df_metadata.shift()
    print("Trigger response = {tr}".format(tr=trigger_response))
    print("Absolute timestamp = {dt_begin}".format(dt_begin=dt_begin))
    print("Ticks per second = {tps}".format(tps=ticks_per_second))
    df_metadata.head()
    
    # Write out as CSV to source directory of SPE file.
    fpath_csv = os.path.splitext(fpath_spe)[0]+'_timestamps.csv'
    df_metadata.to_csv(fpath_csv, quoting=csv.QUOTE_NONNUMERIC)

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print('ERROR: must provide filename of SPE file.')
    else:
        main(sys.argv[1])
