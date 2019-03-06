#!/usr/bin/env python3
"""

Plot tool for MCA files

"""
import logging
import copy
import pdb; 
import sys
import os
import argparse
import matplotlib
matplotlib.use('Qt5Agg')  # nice, but issue with interactive use e.g. in
                          # Jupyter; see
                          # http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
import matplotlib.pyplot as plt
import numpy as np
import glob
import log_helper
import read_mca
import re
from scipy.optimize import curve_fit
from scipy import stats
import scipy.signal as signal
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
def exp_fit(x, a, b, c):
    return a*np.exp(-b*x) + c
def log_fit(x, a, b):
    return np.log(a) -b*x
def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
#plots time vs number of counts: results in the exponential count
def SumArray(y):
    y_all=y[0]
    for item in range (0,len(y)):
        y_all=y_all+y[item]
    return(y_all)
def CalibrationFit(x,y,s):
    if (s=="He"):    
        Constants=[0.191,0.764] # the self calibrating points for He tube in MeV
    else:
        Constants=[1.330,1.17,0.6617] # the calibration of the Co-60 and Cs-137
        if (s=="CeBr"):
            ParameterFit,BoundaryMin,BoundaryMax,BoundaryMaxWide=[1096,966,533],[1072,924,500],[1118,980,570],[30,20,20]
        elif(s=="LaBr"):
            ParameterFit,BoundaryMin,BoundaryMax,BoundaryMaxWide=[555,485,276],[542,473,250],[573,500,980],[20,20,25]
    Gaussian_Fit=[]
    y_all=SumArray(y)
    plt.figure()
    plt.plot(x[0],y_all);
    f = open("datacalib", "w")
    f.write("# x y \n ")        # column names
    np.savetxt(f, np.array([x[0],y_all]).T)
    x_base=x[0];
    if (s=="He"):
        Gaussian_Fit.append(720+np.argmax(y_all[800:1400]))
        Gaussian_Fit.append(3990+np.argmax(y_all[4000:6000]))
        print ("Gaussian Fit")
        print (Gaussian_Fit);
    else:
        for i in range (0,len(Constants)):
            popt, pcov = curve_fit(gauss_function,x_base[BoundaryMin[i]:BoundaryMax[i]],y_all[BoundaryMin[i]:BoundaryMax[i]], p0 = [6000, ParameterFit[i], 7]) #popt[1] the middle of the Gaussian, popt[2] is the sigma, popt[0] is the multiplication He-Tube
        plt.plot(x_base, gauss_function(x[0],popt[0],popt[1],popt[2]));
        Gaussian_Fit.append(popt[1])
    scale, intercept, r_value, p_value, std_err = stats.linregress(Gaussian_Fit,Constants)
    Gaussian_Fit=np.asarray(Gaussian_Fit)
    Constants=np.asarray(Constants)
    # plotRoutine(np.array([Gaussian_Fit,Gaussian_Fit]),np.array([Constants,scale*Gaussian_Fit+intercept]),"Channel no. [-]","Energy",["Data","Fit"],np.array([0,0]),"Calibration")
    plt.figure()
    plt.plot(Gaussian_Fit,Constants,'x',Gaussian_Fit,scale*Gaussian_Fit+intercept)
    return(scale,intercept)

def plotRoutine(xPlotArray,yPlotArray,xlabel,ylabel,legendPlotArray,errorPlotArray,title):
    plt.figure()
    for i in range(len(xPlotArray)):
        # print (xPlotArray[i])
        plt.errorbar(xPlotArray[i],yPlotArray[i],yerr=errorPlotArray[i]);               #plot the data
    plt.title(title)
    plt.legend(legendPlotArray)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.yscale("log")
def load_folder_basic(folder):
    t_array_real=[[],[],[],[]] # here we store the time after each files
    t_livetime=[0,0,0,0] # initializing the array of the livetimes of each of the measurements
    integral=[0,0,0,0] # here we store the total number of events in each file
    error=[0,0,0,0] # here we store the total number of events in each file
    integralPeak=[]
    x=[[],[],[],[]] # array of arrays for the channel values
    y=[[],[],[],[]] # array of arrays for the count values
    t=[0,0,0,0]; # initializing time 0 of the beginning of the first measurment
    mu_0=[[],[],[],[]]; # array of arrays for the mu_0 values
    for id, file in  enumerate(sorted(glob.glob("./"+folder+"/*.mpa"), key=numericalSort)): #looping through the folder
                a=read_mca.read_mca(file)
                for i in range(0,4):
                    for j in range(0,5):
                        a[i].rebin_factor2(5)
                    t[i]=t[i]+a[i].real_time;
                    t_array_real[i].append(t[i]);
                    integral[i]+=a[i].int;
                    if (id==0):
                        x[i]=np.append(x[i],a[i].x)
                        y[i]=np.append(y[i],a[i].y)
                    else:
                        y[i]=y[i]+a[i].y
                    t_livetime[i]=t_livetime[i]+a[i].live_time
    for i in range(0,4):
        if (t_livetime[i]!=0):
            integral[i]=integral[i]/t_livetime[i]
            error[i]=np.sqrt(integral[i])/t_livetime[i]
    return(t_array_real,t_livetime,x,y,t,mu_0,integral,a,error)
def FWHM_plot(folder):
    FWHM=[]
    for file in (sorted(glob.glob("./"+folder+"/*.asc"))):
        a=read_mca.read_mca(file)
        maximum=np.max(a.y)
        index=np.where(a.y==maximum)
        # maximum=argmax(a.y)
        popt, pcov = curve_fit(gauss_function,a.x,a.y, p0 = [1000, maximum, 50]) #popt[1] the middle of the Gaussian, popt[2] is the sigma, popt[0] is the multiplication He-Tube
        plt.plot(a.x,a.y)
        plt.plot(a.x,gauss_function(a.x,popt[0],popt[1],popt[2]) )
        FWHM.append(popt[1])
    return(FWHM)
    
    plt.plot()
def load_folder(folder,s): # goes through the folder, loads all the files, extracts the live times and plots the number of events over the live times
    integral_HZB_BM=[]
    error_HZB_BM=[]
    folder_Name=[]
    integral_VND_BM=[]
    error_VND_BM=[]
    t_array_real=[[],[],[],[]] # here we store the time after each files
    t_livetime=[0,0,0,0] # initializing the array of the livetimes of each of the measurements
    integral=[0,0,0,0] # here we store the total number of events in each file
    integralPeak=[]
    x=[[],[],[],[]] # array of arrays for the channel values
    y=[[],[],[],[]] # array of arrays for the count values
    t=[0,0,0,0]; # initializing time 0 of the beginning of the first measurment
    mu_0=[[],[],[],[]]; # array of arrays for the mu_0 values
    if (s=="He"):
        p=[6000,7230,80]
    elif (s=="LaBr"):
        p=[50,604,14]
    elif (s=="CeBr"):
        p=[10,1180,35]
    Flux_graph_length_data=0
    for idx,folder_file in enumerate(sorted(glob.glob("./"+folder+"/*"))):
        # if os.path.splitext(file)[1] == ".mp3":
        if (os.path.splitext(folder_file)[1]!=".asc"):
            print("I loop here")
            print(folder_file)
            folder_Name.append(folder_file)
            for LookFolder in (sorted(glob.glob("./"+folder+"/*"))):
                if (LookFolder==folder_file+".asc"):
                    HZB_BM_data=read_mca.read_mca(LookFolder)
                    integral_HZB_BM.append(HZB_BM_data.int)
                    error_HZB_BM.append(HZB_BM_data.error_total)
                    print ("found match")
                    Flux_graph_length_data+=1
                    break
            (t_array_real,t_livetime,x,y,t,mu_0,integral,a,error_VND)=load_folder_basic(folder_file)
            # plt.figure()
            # for i in range (0,4):
            #     plt.plot(x[i],y[i])
            integral_VND_BM.append(integral)
            error_VND_BM.append(error_VND)
        # if (s=="He"):
        #     mu_0.append(np.argmax(4000+a.y[4000:6000]))
        # else:
        #     try:
        #         a.Gaussian_fit(p[0],p[1],p[2])
        #         mu_0.append(a.mu_0)
        #         integralPeak.append(a.gauss_int)
        #     except RuntimeError:
        #         break
     # np.concatenate(np.asarray(integral_VND_BM),np.asarray(integral_Total))
    integral_VND_BM=np.asarray(integral_VND_BM)
    error_VND_BM=np.asarray(error_VND_BM)
    error_HZB_BM=np.asarray(error_HZB_BM)
    return(t_array_real,integral,x,y,t_livetime,mu_0,a,integral_VND_BM,integral_HZB_BM,folder_Name,error_HZB_BM,error_VND_BM)
def GainMatch(mu_ref,x,y,mu):
    for i in range(len(y)):
        # print (mu_ref)
        # print (mu[i])
        bin_shift=int(round((mu_ref-mu[i])/(x[1]-x[0])))
        y[i]=np.roll(y[i],bin_shift)
    return(x,y)
def SumRebinScaleCalibrate(rebin_factor,a,slope,intercept,total_livetime):
    # a.calibrate(slope,intercept)
    for i in range(0,rebin_factor):
        a.rebin_factor2(5)
    # a.error_y=np.sqrt(a.y/total_livetime**2)
    # a.scale(1/total_livetime)
    return(a)
def plot_all(spectra):
    """ plots all MCASpectrum objects in the spectra list in one plot """
    fig, ax = plt.subplots()
    for s in spectra:
        # if we have bin "centers" on x:
        plt.plot(s.x, s.y, label=s.name, ls="steps")
        # if we have bin boundaries on x:
        #plot_binned_data(ax, s.x, s.y)
    plt.legend()
    plt.ylim(ymin=1)
    # plt.yscale('log')
    # plot ratio if more than one file loaded
    fig, ax = plt.subplots()
    if len(spectra)>1:
        for idx, s in enumerate(spectra):
            if idx==0:
                continue
            ratio = s.y/spectra[0].y
            plt.plot(s.x, ratio, label=s.name, ls="steps")
    plt.title("Ratio to {}".format(spectra[0].name))
    plt.axhline(y=1, color = 'black')
    plt.legend()

# to plot already-binned data: 
# from: https://stackoverflow.com/a/19361027
def plot_binned_data(axes, binedges, data, *args, **kwargs):
    #The dataset values are the bin centres
    x = (binedges[1:] + binedges[:-1]) / 2.0
    #The weights are the y-values of the input binned data
    weights = data
    return axes.hist(x, bins=binedges, weights=weights, *args, **kwargs)


if __name__ == "__main__":
    log = log_helper.init_logging('plot_mca')

    # command line argument parsing
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv.pop(0)),
        description="Plotter for MCA data files recorded with Maestro/Amptek MCA DAQ"
    )
    parser.add_argument(
        "-l",
        "--log-level",
        default="info",
        help=
        "Sets the verbosity of log messages where LEVEL is either debug, info, warning or error",
        metavar="LEVEL")
    parser.add_argument(
        "-i",
        "--interactive",
        action='store_true',
        help=
        "Drop into an interactive IPython shell instead of showing default plots"
    )
    parser.add_argument(
        "file",
        nargs='+',
        help=
        "The MCA data file(s) to be processed; additional info STRING to be included in the plot legend can be added by specifiying FILE:STRING"
    )

    parser.add_argument(
        "-f",
        "--folder",
        help=
        "Pass the folder and it will sum all the files in the folder"
    )
    parser.add_argument(
        "-b",
        "--background",
        help=
        "Pass the background folder"
    )
    
    parser.add_argument(
        "-a",
        "--activation",
        help=
        "Pass the background folder"
    )
    parser.add_argument(
        "-c",
        "--calibration",
        help=
        "Pass the calibration folder"
    )
    parser.add_argument(
        "-other",
        "--other",
        help=
        "Pass the asc folder"
    )


    
    # parse the arguments

    args = parser.parse_args(sys.argv)
    # if (args.activation!=None):
        # activation_plot(args.activation)
    # FWHM=FWHM_plot(args.fwhm)
    # print(FWHM)
    # pdb.set_trace()
    if (args.folder!=None):
        detector="He"
        (t_array_real,integral,x,y,t_livetime,mu_0,a,VND_int,HZB_BM_int,folder_Name,error_HZB,error_VND)=load_folder(args.folder,detector)
        (t_array_real_noV,integral_noV,x_noV,y_noV,t_livetime_noV,mu_0_noV,a_noV,VND_int_noV,HZB_BM_int_noV,folder_Name_noV,error_HZB_noV,error_VND_noV)=load_folder(args.background,detector)
        # (t_array_real_full,integral_full,x_full,y_full,t_livetime_full,mu_0_full,a_full,VND_int_full,HZB_BM_int_full,folder_Name_full,error_HZB_full,error_VND_full)=load_folder(args.other,detector)
        plt.figure()
        for i in range(0, 4):
            plt.plot(x[i],y[i]/t_livetime[i])
        plt.legend(["Channel1","Channel2","Channel3","Channel4"] )
        plt.title(folder_Name)
        plt.ylabel("n/s")
        plt.xlabel("Channels[-]")
        # VND_int-=VND_int_noV
        ######################################Kelly's graph##################
        # plt.figure()
        # # z=[5,10,13,15,18,25,30] ##for the data with varying the slit
        # z=[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23] ##for the other data
        # z1=[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23] ##for the other data
        l=[2,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
        z=[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
        # l=[5,8,11,14,17,20,22]
        # l=[5,7,9,10,11,13,14,15,16,17,19,21,22]
        fig,ax1=plt.subplots()
        ax1.set_ylabel("Vanadium BM uncorrected [n/s]")
        for i in range(0, 1):
            plt.errorbar(l,VND_int[:,i],yerr=error_VND[:,i],fmt='x')
            plt.errorbar(z,VND_int_noV[:,i],yerr=error_VND_noV[:,i], fmt='o')
        # for i in range(0, 4):
            # plt.errorbar(l,VND_int_full[:,i],yerr=error_VND_full[:,i],fmt='o')
        plt.xlabel("Chopper opening [degree] (flux regulation)")
        # plt.legend(["Counts on the Helium tubes (V monitor)"],loc="upper left")
        plt.legend(["Counts when the Vanadium is present","Counts when the Vanadium is not present"])
        # plt.title("Counts on Helium Tubes as dependent on the ")
        # plt.savefig('V-monitor_foreground/background.png')

        ################################BASIC DEPENDENCY GRAPH corrected ######################
        # ax2=ax1.twinx()
        # ax2.set_ylabel("HZB BM uncorrected [n/s]")
        # plt.plot(z,HZB_BM_int_noV,"x",color="red")
        # plt.legend(["Counts on HZB_referrence_monitor"],loc="lower right" )
        # plt.title("Scan"+args.folder+"the referrence monitor with the Vanadium Out, relative counts")

        # #############################################Richards graph_corrected###################
        # plt.figure()
        # y_error_numerator_noV=np.sqrt(error_VND_noV[:,i]**2+error_HZB_noV**2)
        # y_error_noV=np.sqrt((y_error_numerator_noV/(HZB_BM_int_noV-VND_int[:,i]))**2+(error_HZB_noV/HZB_BM_int_noV)**2)
        # for i in range(0, 4):
        #     plt.errorbar(HZB_BM_int_noV,(HZB_BM_int_noV-VND_int[:,i])/HZB_BM_int_noV,yerr=y_error_noV,xerr=error_HZB, fmt='x')
        # plt.legend(["Channel1","Channel2","Channel3","Channel4"] )
        # plt.title("Scan"+args.folder+"the referrence monitor with the Vanadium Out, relative counts")
        # plt.ylabel("(Vanadium BM- HZB BM)/HZB BM n/s")
        # plt.xlabel("HZB BM[n/s]")
        #######################################CORRELATIN FOR THE ABSOLUTE FLUX##############################
        #######################################CORRELATIN FOR THE ABSOLUTE FLUX##############################
        #######################################CORRELATIN FOR THE ABSOLUTE FLUX##############################
 ################################BASIC DEPENDENCY GRAPH corrected for absolute flux ######################
        # plt.figure()
        # pdb.set_trace()
        # HZB_BM_int=np.asarray(HZB_BM_int)/(3.9*10**(-3))
        # HZB_BM_int_noV=np.asarray(HZB_BM_int_noV)/(3.9*10**(-3))
        # error_HZB_noV=np.asarray(error_HZB_noV)/(3.9*10**(-3))
        # ##########################percentage of the neutrons scattered
        # for i in range(0,4):
        #     if (i==0):
        #         VND_int[i]=VND_int[i]=np.asarray(VND_int[i])*(4*3.14*5**2)/25
        #     else:
        #         VND_int[i]=VND_int[i]=np.asarray(VND_int[i])*(4*3.14*15**2)/25

        # VND_int=np.asarray(VND_int)/(np.exp(-7.09*10**(22)*5.08*10**(-24)*3.15*0.1)) ###the neutrons that are scattered
        # # EXP(-7.09*10^(22)*5.08*10^(-24)*2*A3*10^(-1)))
        # for i in range(0, 4):
        #     plt.errorbar(HZB_BM_int_noV,VND_int[:,i],yerr=error_VND[:,i],xerr=error_HZB_noV, fmt='x')
        # plt.legend(["Channel1","Channel2","Channel3","Channel4"] )
        # plt.title("Scan"+args.folder+"the referrence monitor with the Vanadium Out, absolute dependency")
        # plt.ylabel("Vanadium BM n/s")
        # plt.xlabel("HZB BM absolute flux [n/s]")

        # # #############################################Richards graph_corrected for absolute flux ###################
        # plt.figure()
        # y_error_numerator_noV=np.sqrt(error_VND_noV[:,i]**2+error_HZB_noV**2)
        # y_error_noV=np.sqrt((y_error_numerator_noV/(HZB_BM_int_noV-VND_int[:,i]))**2+(error_HZB_noV/HZB_BM_int_noV)**2)
        # for i in range(0, 4):
        #     plt.errorbar(HZB_BM_int_noV,(HZB_BM_int_noV-VND_int[:,i])/HZB_BM_int_noV,yerr=y_error_noV,xerr=error_HZB, fmt='x')
        # plt.legend(["Channel1","Channel2","Channel3","Channel4"] )
        # plt.title("Scan"+args.folder+"the referrence monitor with the Vanadium Out, absolute dependency")
        # plt.ylabel("(Vanadium BM- HZB BM)/HZB BM [n/s")
        # plt.xlabel("HZB BM absolute flux,[n/s]")


        plt.show()
        # (t_array_real_bg,integral_bg,x_bg,y_bg,t_livetime_bg, mu_0_bg,a)=load_folder(args.background,detector)
        # (_t_array_real_calib,_integral_calib,x_calib,y_calib,_t_livetime_calib,_mu_0_calib,a)=load_folder(args.calibration,detector)
        # (scale,intercept)=CalibrationFit(x_calib,y_calib,detector);
        # (beta_decay)=plot_sum(args.folder,args.background,scale,intercept,detector)
        #integral_bg=np.sum(np.asarray(integral_bg))/np.sum(np.asarray(t_livetime_bg))*np.asarray(t_livetime)
        #integral_bg=np.sum(np.asarray(integral_bg))/np.sum(np.asarray(t_livetime_bg))*np.asarray(t_livetime)
        #integral_new=np.asarray(integral)-np.asarray(integral_bg)
        # activation_plot_exponential(t_array_real,integralPeak,'Integration over the 1.43Peak')
    # parse file names and extract additionally provided info 
    fileNames = []
    fileDescr = []
    for thisFile in args.file:
        s = thisFile.strip().split(':', 1)  # try to split the string
        if (len(s) == 1):
            # didn't work, only have one entry; use file name instead
            fileNames.append(s[0])
            fileDescr.append(s[0])
        else:
            fileNames.append(s[0])
            fileDescr.append(s[1])

    if args.log_level:
        # Convert log level to upper case to allow the user to specify --log-level=DEBUG or --log-level=debug
        numeric_level = getattr(logging, args.log_level.upper(), None)
        if not isinstance(numeric_level, int):
            log.error('Invalid log level: %s' % args.log_level)
            sys.exit(2)
        log.setLevel(numeric_level)
        # apply setting also to our data reader
        log_reader = log_helper.init_logging('read_mca') ## set up logging
        log_reader.setLevel(numeric_level)

    log.debug("Command line arguments used: %s ", args)
    log.debug("Libraries loaded:")
    log.debug("   - Matplotlib version {}".format(matplotlib.__version__))
    log.debug("   - Numpy version {}".format(np.__version__))

    # now load measurements from file(s)
    msrmts = [] # list to hold all measurements
    for idx, f in enumerate(fileNames):
        s = read_mca.read_mca(f)
        if s is None:
            log.error("Error encountered, skipping file {}".format(f))
            continue
        s.name = fileDescr[idx]
        msrmts.append(s)

    if log.error.counter>0:
        log.warning("There were "+str(log.error.counter)+" error messages reported")

    if (args.interactive):
        print(" Interactive IPython shell ")
        print(" ========================= ")
        print(" Quick command usage:")
        print("  - 'who' or 'whos' to see all (locally) defined variables")
        print(
            "  - list of loaded data files: 'fileNames' with labels stored in 'fileDescr'"
        )
        print("  - data is stored in numpy arrays")
        print("  - access to data via list of MCASpectrum objects called 'msrmts'")
        print("    e.g.")
        print("       for m in msrmts:")
        print("          print(m.y)")
        print("  - if the plots are shown only as black area, run '%gui qt'")
        print("  - to make cmd prompt usable while plots are shown use 'plt.ion()' for interactive mode")
        import IPython
        IPython.embed()
    else:
        # put matplotlib into non-interactive mode
        plt.ioff()
        print(msrmts)
        plot_all(msrmts)
        log.info("Drawing plots....")
        ## final commands to show resulting plots (not needed in an interactive session)
        plt.show(block=False)  # block to allow "hit enter to close"
        plt.pause(0.001)  # <-------
        input("<Hit Enter To Close>")
        plt.close('all')
    

    if log.error.counter>0:
        log.warning("There were "+str(log.error.counter)+" error messages reported")

