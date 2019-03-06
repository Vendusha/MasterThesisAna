#!/usr/bin/env python
import logging
import log_helper
import csv
import numpy as np
import os
from scipy.optimize import curve_fit
######Only for the purpose of debugging 
import argparse
import sys
import os

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

class MCASpectrum:
    """A class to hold our spectrum measurement data and meta data (such as live_time)"""
    def __init__(self, filename):
        self.filename = filename
        self.x = np.array(np.zeros(0))   ## creates a new empty array; will later store our x values
        self.y = np.array(np.zeros(0))    ## creates a new empty array; will later store our y values
        self.name = os.path.splitext(     ## a more descriptive name, can be used e.g. in legends
                os.path.basename(filename))[0] ## init with file name without extension
        self.live_time = 0
        self.time = 0  ## start time of measurement; TODO: use a well-defined type to store this info
        self.int=0
        self.error_total=0
        self.real_time=0
        self.error_y = np.array(np.zeros(0))    ## creates a new empty array; will later store our y_error values
        self.mu_0=0
        self.gauss_int=0
        self.p=0
        self.sigma=0
    def subtract(self, m):
        self.y = self.y - m.y
        ## these are spectra: cannot have counts below 0, so remove these here and set them to 0 instead:
        ## self.y[self.y < 0] = 0;
    def scale(self, scale):
        self.y *= scale
    def calibrate(self, slope, intercept):
        self.x = self.x*slope + intercept
    def Gaussian_fit(self,a,x0,sigma):
        popt, pcov = curve_fit(gauss_function, self.x, self.y, p0 = [a, x0, sigma]) #popt[1] the middle of the Gaussian, popt[2] is the sigma, popt[0] is the multiplication
        self.mu_0=popt[1];
        self.p=popt[0]
        self.sigma=popt[2]
        self.gauss_int=np.trapz(gauss_function(self.x,self.p,self.mu_0,self.sigma),x=self.x)
    def rebin_factor2(self, n): #TODO n as a number of binning
        self.y = self.y[1::2] + self.y[::2]
        self.x = self.x[::2]




def read_mca(filename):
    """Reads in a data file (csv format) stored by the Amptek or Maestra MCA DAQ software packages and returns a 'MCASpectrum' object. """
    log = log_helper.init_logging('read_mca') ## set up logging
    m = MCASpectrum(filename) ## create a new Spectrum measurement object; this is what we return in the end
    log.info("Reading data from file '" + filename + "'")
    m = MCASpectrum(filename) ## create a new Spectrum measurement object; this is what we return in the end
    m4=[MCASpectrum(filename),MCASpectrum(filename),MCASpectrum(filename),MCASpectrum(filename)] ## for the 4 channel MCA, define a python tuple of 4 objects
    if str(filename).endswith('.mca'):
        # Amptek file format
        log.debug("Identified as 'Amptek' format")

        with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
            reader = csv.reader(f) ## use the python csv module to parse the file
            interval = [0, 0] ## start/stop channel numbers used to assign correct x values to the data points
            ## first parse the "header" of the data file (until the '$DATA:' line) containing all the meta data
            index=0;
            for row in reader:
                if 'REAL_TIME' in row[0]:
                    ##real time of the measurement
                    real_time = float(row[0].split(' - ')[1])
                    m.real_time = real_time
                if 'LIVE_TIME' in row[0]:
                    ## this item gives the live_time of the measurement
                    live_time = float(row[0].split(' - ')[1])
                    m.live_time = live_time
                if 'START_TIME' in row[0]:
                    log.info('Parsing start time info')
                    m.time = row[0].split(' - ')[1]
                if row[0] == '<<DATA>>':
                    row = next(reader)
                    while (row[0] != '<<END>>'):
                        index +=1
                        if index>2:
                            m.y = np.append(m.y, int(row[0]))
                            m.int=m.int+int(row[0])
                        else:
                            m.y = np.append(m.y, 0)
                        row = next(reader)
                # determine bin boundaries from upper/lower thresholds
                # NOTE: check that software actually removes these channels from the data and not just sets them to 0!
                if 'MCSL=' in row[0]:
                    log.info("Parsing DATA header info")
                    words = row[0].split('=')
                    interval[0] = int(words[1].split(';')[0])
                if 'MCAC=' in row[0]:
                    words = row[0].split('=')
                    interval[1] =int(words[1].split(';')[0])-1
                    log.info("Done with header parsing")
                    ## finds the maximum number of bins in the spectra and creates a range of the channels spaced by 1
            m.x = np.arange(interval[0], interval[1]+1, 1)
            log.info("Loaded all data from file")
    elif str(filename).endswith('.Spe'):
        # Maestro ASCII file
        log.debug("Identified as 'Maestro' format")
        with open(filename, newline='') as f:
            reader = csv.reader(f) ## use the python csv module to parse the file
            interval = []          ## start/stop channel numbers used to assign correct x values to the data points
            ## first parse the "header" of the data file (until the '$DATA:' line) containing all the meta data
            for row in reader:
                if row[0] == '$MEAS_TIM:':
                    ## this item gives the live_time of the measurement
                    log.debug("Parsing MEAS_TIM header info")
                    row = next(reader)
                    live_time = [int(s) for s in row[0].split(' ')]
                    m.live_time = duration[1] ## two parts: real time/live time; take the second
                if row[0] == '$DATA:':
                    ## this is the last part of the header and contains the start/stop channel numbers
                    log.debug("Parsing DATA header info")
                    row = next(reader)
                    interval = [int(s) for s in row[0].split(' ')]
                    ## "DATA" is the last item: stop with the header processing
                    break
            ## TODO: make sure that the file does not end before we have parsed the header!
            log.debug("Done with header parsing")
            nchannel = int(interval[1]-interval[0])+1
            m.y = np.array(np.zeros(nchannel))
            ## continue, now reading data
            for idx, row in enumerate(reader):
                if idx >= nchannel:
                    break
                m.y[idx] = int(row[0])
            m.x = np.arange(interval[0], interval[1]+1,1)
            log.debug("Loaded all data from file")
    elif str(filename).endswith('asc'):
        # log.debug("Identified as HZB format")
        # m.x,m.y=np.loadtxt(filename,unpack=True)
        # print(m.y)
        enable_loading=0;
        enable_livetime_readout=1
        with open(filename, newline='') as f:
            reader = f.readlines()
            i=0;
            for row in reader:
                if enable_livetime_readout==1 and "events /" in row:
                    FirstSplit=row.split("events /")
                    SecondSplit=FirstSplit[1].split("s, cpsMax")
                    m.live_time=int(SecondSplit[0])
                    enable_livetime_readout=0;
                if (enable_loading==0 and row[0:5]!='TrigA' and row[0:2]!="20"):
                    enable_loading=1
                if (enable_loading==1 and i<3):
                    i+=1
                if (i==3): # skip two lines
                    m.y = np.append(m.y, int(row.split()[1]))
                    m.int+=int(row.split()[1])
            m.int=m.int/m.live_time
            m.error_total=np.sqrt(m.int)/m.live_time
        return(m)
    elif str(filename).endswith('.mpa'):
        enable_loading=0 # first we are processing bunch of different headers
        # FastCom MPA file
        log.debug("Identified as 'FastCom MCA' format")
        with open(filename, newline='') as f:
            reader = f.readlines()
            # reader = csv.reader(f) ## use the python csv module to parse the file
            interval = []          ## start/stop channel numbers used to assign correct x values to the data points
            k=0;
            ## first parse the "header" of the data file (until the '$DATA:' line) containing all the meta data
            for row in reader:
                if (row[0:5]=='[CDAT'):
                    break;
                # if (enable_loading==1 and row[0:5]=='[DATA'):
                #     enable_loading=0;
                #     k+=1

                if (enable_loading==1):
                    m4[k].y = np.append(m4[k].y, int(row.split()[1]))
                    m4[k].int+= int(row.split()[1])                ###########################FOR THE HELIUM TUBES IT IS ENOUGH TO DO THE LOADING OF THE 6000 channels######################
                if (enable_loading==1 and len(m4[k].y)==8192):
                    enable_loading=0;
                    k+=1
                if row[0:4] == '[ADC': ###for extracting of the information one has to look inside of ADC settings section
                    i=int(row[4])-1 # extract the number of current channel and renaming the channels from 1-4 to 0-3
                    ## this item gives the live_time of the measurement
                if row[0:9]=="livetime=":
                    m4[i].live_time=float(row[9:])
                    ## this item gives the real_time of the measurement
                    log.debug("Parsing MEAS_TIM header info")
                if row[0:9]=="realtime=":
                    m4[i].real_time=float(row[9:])
                if row[0:6] == '[DATA]':
                    enable_loading=1;
            for j in range(0,4):
                    m4[j].x=np.arange(0,len(m4[j].y),1)
                    # m4[j].int=m4[j].int/m4[j].live_time
                    # m4[j].error_total=np.sqrt(m4[j].int)/m4[j].live_time
        return(m4)
    else:
        log.error("Unknown extension/format: {}".format(os.path.splitext(os.path.basename(filename))[1]))
        return None
    return m
# if __name__ == "__main__":
#     log = log_helper.init_logging('read_mca')

#     # command line argument parsing
#     parser = argparse.ArgumentParser(
#         prog=os.path.basename(sys.argv.pop(0)),
#         description="Plotter for MCA data files recorded with Maestro/Amptek MCA DAQ"
#     )
#     parser.add_argument(
#         "-l",
#         "--log-level",
#         default="info",
#         help=
#         "Sets the verbosity of log messages where LEVEL is either debug, info, warning or error",
#         metavar="LEVEL")
#     parser.add_argument(
#         "-i",
#         "--interactive",
#         action='store_true',
#         help=
#         "Drop into an interactive IPython shell instead of showing default plots"
#     )
#     parser.add_argument(
#         "file",
#         nargs='+',
#         help=
#         "The MCA data file(s) to be processed; additional info STRING to be included in the plot legend can be added by specifiying FILE:STRING"
#     )
#     args = parser.parse_args(sys.argv)
#     if (args.file!=None):
#         filename=(args.file)
#         print (filename)
#         read_mca(TEST023.mca)
