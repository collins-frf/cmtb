# -*- coding: utf-8 -*-
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os, getopt, sys, shutil, glob, logging, yaml, re, pickle
import datetime as DT
from subprocess import check_output
import numpy as np
from getdatatestbed.getDataFRF import getObs, getDataTestBed
from testbedutils import fileHandling
from oct2py import octave
import matlab
import math
import netCDF4 as nc

Q = 0
def assim_currents(X, currents_object, obs, i):
    awac_value, awac_indice = find_nearest(X, currents_object['xFRF'])
    obs['v']['d'] = np.append(obs['v']['d'], (matlab.float64([currents_object['aveV'][i]])))
    obs['v']['ind'] = np.append(obs['v']['ind'], (matlab.float64([awac_indice])))
    return obs

def assim_waves(X, waves_object, obs, i):
    hs_value, hs_indice = find_nearest(X, waves_object['xFRF'])
    obs['H']['d'] = np.append(obs['H']['d'], (matlab.float64([waves_object['Hs'][i]])))
    obs['H']['ind'] = np.append(obs['H']['ind'], (matlab.float64([hs_indice])))
    return obs

def calculate_Ch(prior, spec8m, X, delta_t, Q):
    # calculate Q(x,t) (measured process error) according to holman et al 2013
    # Q = Cq * Hmo ^ (-[(x-x0)/sigma_x]^2)
    # deltat is difference from survey to assimilation step;
    # add q each time-step to increase uncertainty as the Ch decreases in posterior posterior.ch + new S *N *S with just the tiny delta_t
    Cq = 0.067  # from sandy duck experiment
    Hmo = np.mean(spec8m['Hs'])  # significant wave height of highest 1/3 of waves
    x0 = 50  # x0 and sigma_x reflect the typical location of breaking waves at this particular beach
    sigma_x = 150
    delta_t = delta_t/(60*60*24)
    print("Delta t:", delta_t)
    xx = np.meshgrid(X)
    Lx = 25  # decorrelation length scale, larger Lx leads to smoother results
    N = np.exp((-(xx - np.transpose(xx)) ** 2) / (2 * Lx))
    Q = Q + (Cq * Hmo) ** (np.power(-((X - x0) / sigma_x), 2))*delta_t
    S = np.diag(np.sqrt(Q))
    Ch = S * N * S

    prior['Ch'] = Ch  # populate elevation covariance
    return prior, Q

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def Master_1DVar_run(inputDict):
    """This function will run CMS with any version prefix given start, end, and timestep.

    Args:
      inputDict: a dictionary that is read from the input yaml

    Returns:
      None

    """
    ## unpack input Dictionary
    version_prefix = inputDict['modelSettings']['version_prefix']
    endTime = inputDict['endTime']
    startTime = inputDict['startTime']
    simulationDuration = inputDict['simulationDuration']
    workingDir = os.path.join(inputDict['workingDirectory'], 'waveModels')
    generateFlag = inputDict['generateFlag']
    runFlag = inputDict['runFlag']
    analyzeFlag = inputDict['analyzeFlag']
    pFlag = inputDict['pFlag']
    model = inputDict.get('model', '1DVar')
    modeldir = inputDict['modelExecutable']
    log = inputDict.get('logging', True)

    # __________________pre-processing checks________________________________
    fileHandling.checkVersionPrefix(model, inputDict)
    # __________________input directories________________________________
    baseDir = os.getcwd()  # location of working directory
    # check executable
    if inputDict['modelExecutable'].startswith(baseDir):  # change to relative path
        inputDict['modelExecutable'] = re.sub(baseDir, '', inputDict['modelExecutable'])
    workingDirectory = os.path.join(workingDir, model.lower(), version_prefix)
    inputDict['netCDFdir'] = os.path.join(inputDict['netCDFdir'], 'waveModels')
    inputDict['path_prefix'] = workingDirectory
    # ______________________ Logging  ____________________________
    # auto generated Log file using start_end time?
    # LOG_FILENAME = fileHandling.logFileLogic(workingDirectory, version_prefix, startTime, endTime, log=log)
    # __________________get time list to loop over________________________________
    try:
        projectEnd = DT.datetime.strptime(endTime, '%Y-%m-%dT%H:%M:%SZ')
        projectStart = DT.datetime.strptime(startTime, '%Y-%m-%dT%H:%M:%SZ')
    except TypeError:  # if input date was parsed as datetime
        projectEnd = endTime
        projectStart = startTime
    # This is the portion that creates a list of simulation end times
    dt_DT = DT.timedelta(0, simulationDuration * 60 * 60)  # timestep in datetime
    # make List of Datestring items, for simulations
    dateStartList = [projectStart]
    dateStringList = [dateStartList[0].strftime("%Y-%m-%dT%H:%M:%SZ")]
    for i in range(int(np.ceil((projectEnd - projectStart).total_seconds() / dt_DT.total_seconds())) - 1):
        dateStartList.append(dateStartList[-1] + dt_DT)
        dateStringList.append(dateStartList[-1].strftime("%Y-%m-%dT%H:%M:%SZ"))
    # fileHandling.displayStartInfo(projectStart, projectEnd, version_prefix, LOG_FILENAME, model)
    # ______________________________gather all data _____________________________
    if generateFlag == True:
        go = getObs(projectStart, projectEnd)  # initialize get observation
        bathygo = getDataTestBed(projectStart, projectEnd)
        #bathyTransect = bathygo.getBathyTransectFromNC(method=1)  # grab bathymetry transects
        bathyTransect = bathygo.getBathyIntegratedTransect(method=1)  # grab bathymetry transects

        pier_wind = go.getWind()

        adop_35m_currents = go.getCurrents(gaugenumber='adop-3.5m', roundto=1)
        awac_45m_currents = go.getCurrents(gaugenumber='awac-4.5m', roundto=1)
        awac_6m_currents = go.getCurrents(gaugenumber='awac-6m', roundto=1)
        #awac_8m_currents = go.getCurrents(gaugenumber='awac-8m', roundto=1)
        #awac_11m_currents = go.getCurrents(gaugenumber='awac-11m', roundto=1)

        spec100m = go.getWaveSpec(gaugenumber='xp100m', specOnly=False)  # grab 100m array observations
        spec125m = go.getWaveSpec(gaugenumber='xp125m', specOnly=False)  # grab 125m pressure array obs
        #spec150m = go.getWaveSpec(gaugenumber='xp150m', specOnly=False)  # grab 150m pressure array obs
        spec150m = None
        spec200m = go.getWaveSpec(gaugenumber='xp200m', specOnly=False)  # grab 200m pressure array obs
        spec8m = go.getWaveSpec(gaugenumber='8m-array', specOnly=False)  # grab 8m array observations


        speed = pier_wind['windspeed']
        direc = np.deg2rad(pier_wind['winddir']-71.8)
        cd = 0.002 #reniers et al. 2004
        tauw = cd*(speed**2)*np.cos(direc)
        wind_indice_skip = len(tauw)/len(dateStringList)

        h = bathyTransect['elevation'][200,35:]
        X = np.linspace(35, np.max(bathyTransect['xFRF']), len(h))

        delta_t = spec8m['epochtime'][0] - bathyTransect['time'].timestamp() #in days

        # load example mat file to populate prior and observations with data for model run
        input = octave.load('./data/waveModels/1DVar/input_oct.mat')
        prior = input['prior']

        X = np.reshape(X, (len(X), 1))
        h = np.reshape(h, (len(h), 1))
        Q=0
        prior, Q = calculate_Ch(prior, spec8m, X, delta_t, Q)
        # populate prior with data
        # Prior.ka_drag is constant
        # Prior.Ctheta0 is constant
        # Prior.CH0 is constant
        # Prior.Cka is constant
        prior['x'] = X  # populate with bathy cross-shore distance data
        prior['h'] = -h  # populate with elevation data
        prior['hErr'] = .1 # populate with elevation error data .1

        obs_indices_h = []
        obs_indices_v = []

    # ________________________________________________ RUN LOOP ________________________________________________
    i = 0
    print(dateStringList)
    for time in dateStringList:
        obs = input['obs']
        obs['tauw'] = tauw[0]
        print(obs['tauw'])
        try:
            print('-------------------------------Beginning Simulation {}-------------------------------'.format(
                DT.datetime.now()))

            if runFlag == True:  # run model
                os.chdir(modeldir)  # changing locations to where the model will be ran from for octave
                print(os.getcwd())
                print('Running {} Simulation'.format(model.upper()))
                dt = DT.datetime.now()

                prior['theta0'] = np.deg2rad(spec8m['waveDp'][i]-71.8)  # populate with theta0 data
                prior['H0'] = spec8m['Hs'][i]  # populate with offshore wave height data
                prior['sigma'] = 2 * np.pi * spec8m['peakf'][i]  # calculate offshore sigma value"""
                print("theta: ", prior['theta0'])
                # select indices of bathymetry where we have Hs measurements
                obs['H']['d'] = []
                obs['H']['ind'] = []
                if spec100m != None:
                    try:
                        obs = assim_waves(X, spec100m, obs, i)
                    except:
                        pass
                if spec125m != None:
                    try:
                        obs = assim_waves(X, spec125m, obs, i)
                    except:
                        pass
                if spec150m != None:
                    try:
                        obs = assim_waves(X, spec150m, obs, i)
                    except:
                        pass
                if spec200m != None:
                    try:
                        obs = assim_waves(X, spec200m, obs, i)
                    except:
                        pass
                if spec8m != None:
                    try:
                        obs = assim_waves(X, spec8m, obs, i)
                    except:
                        pass

                for obs_timestep_ind in obs['H']['ind']:
                    obs_indices_h.append(obs_timestep_ind)
                try:
                    obs['H']['e'] = matlab.float64(obs['H']['e'][0][:len(obs['H']['d'])])
                except:
                    obs['H']['e'] = matlab.float64(obs['H']['e'][:len(obs['H']['d'])])
                print("Wave Heights to Assimilate: ", obs['H']['d'])
                print("Wave Indices to Assimilate: ", obs['H']['ind'])
                print("Wave Errors to Assimilate: ", obs['H']['e'])

                # select indices of bathymetry where we have v measurements
                obs['v']['d'] = []
                obs['v']['ind'] = []
                if adop_35m_currents != None:
                    try:
                        obs = assim_currents(X, adop_35m_currents, obs, i)
                    except:
                        pass
                if awac_45m_currents != None:
                    try:
                        obs = assim_currents(X, awac_45m_currents, obs, i)
                    except:
                        pass
                if awac_6m_currents != None:
                    try:
                        obs = assim_currents(X, awac_6m_currents, obs, i)
                    except:
                        pass
                """if awac_8m_currents != None:
                    try:
                        obs = assim_currents(X, awac_8m_currents, obs, i)
                    except:
                        pass"""
                #if awac_11m_currents != None:
                    #obs = assim_currents(X, awac_11m_currents, obs, i)

                for obs_timestep_ind in obs['v']['ind']:
                    obs_indices_v.append(obs_timestep_ind)
                try:
                    obs['v']['e'] = matlab.float64(obs['v']['e'][0][:len(obs['v']['d'])])
                except:
                    obs['v']['e'] = matlab.float64(obs['v']['e'][:len(obs['v']['d'])])

                print("Currents to Assimilate: ", obs['v']['d'])
                print("Current Indices to Assimilate: ", obs['v']['ind'])
                print("Current Errors to Assimilate: ", obs['v']['e'])

                obs['tauw'] = tauw[int((i+1)*wind_indice_skip)]

                posterior, diagnostics = octave.feval('C:/cmtb/bin/1DVar/assim_1dh.m', prior, obs, 1, nout=2)#, representers = \
                print('Simulation took %s ' % (DT.datetime.now() - dt))

            prior = posterior
            try:
                delta_t = spec8m['time'][i+1] - spec8m['time'][i]
                delta_t = delta_t.seconds
            except:
                pass
            prior, Q = calculate_Ch(prior, spec8m, X, delta_t, Q)

            i += 1
            print('-------------------------------SUCCESS-----------------------------------------')

        except Exception as e:
            print('<< ERROR >> HAPPENED IN THIS TIME STEP\n{}'.format(e))
            logging.exception('\nERROR FOUND @ {}\n'.format(time), exc_info=True)
            os.chdir(modeldir)

    obs_indices_h = np.unique(obs_indices_h)
    obs_indices_v = np.unique(obs_indices_v)
    print("Data assimilated at crosshore locations:")
    print(obs_indices_h*(X[-1]/len(h)))
    print(obs_indices_v*(X[-1]/len(h)))
    Q = np.squeeze(np.transpose(Q))
    fig = plt.figure(figsize=(6, 9))
    grid = gridspec.GridSpec(2, 2, figure=fig)
    ax0 = fig.add_subplot(grid[0, 0])
    ax0.set_title("Prior")
    ax0.set_ylim(ymax=0, ymin=-11)
    ax0.set_xlabel("Crosshore distance (m)")
    ax0.set_ylabel("Depth")
    ax0.set_xlim(xmax=1000)
    ax0.plot(X, h, color='black')
    for obs_loc in obs_indices_h:
        ax0.axvline(x=obs_loc*(X[-1]/len(h)), color='r', linestyle='--')
    for obs_loc in obs_indices_v:
        ax0.axvline(x=obs_loc*(X[-1]/len(h)), color='orange', linestyle='--')

    ax1 = fig.add_subplot(grid[0, 1])
    ax1.set_title("Posterior")
    ax1.set_xlabel("Crosshore distance (m)")
    ax1.set_ylabel("Depth")
    ax1.set_ylim(ymax=0, ymin=-11)
    ax1.set_xlim(xmax=1000)

    ax1.errorbar(X, -posterior['h'], yerr=np.squeeze(posterior['hErr']), color='black', errorevery=5, capsize=2, ecolor='pink')
    #ax1.plot(X, -posterior['h'], color='black')

    for obs_loc in obs_indices_h:
        ax1.axvline(x=obs_loc*(X[-1]/len(h)), color='r', linestyle='--')
    for obs_loc in obs_indices_v:
        ax1.axvline(x=obs_loc*(X[-1]/len(h)), color='orange', linestyle='--')

    ax2 = fig.add_subplot(grid[1, 0])
    ax2.set_title("prior - posterior")
    ax2.scatter(X, h+posterior['h'], color="black")
    ax2.set_xlabel("Crosshore distance (m)")
    ax2.set_ylabel("Elevation change (m)")
    ax2.set_xlim(xmax=1000)

    for obs_loc in obs_indices_h:
        ax2.axvline(x=obs_loc*(X[-1]/len(h)), color='r', linestyle='--')
        temp = obs_loc
    ax2.axvline(x=temp * (X[-1] / len(h)), color='r', linestyle='--', label="wave heights")
    for obs_loc in obs_indices_v:
        ax2.axvline(x=obs_loc*(X[-1]/len(h)), color='orange', linestyle='--')
        temp = obs_loc
    ax2.axvline(x=temp * (X[-1] / len(h)), color='orange', linestyle='--', label="currents")

    ax2.legend()

    ax3 = fig.add_subplot(grid[1, 1])
    ax3.plot(X, Q)
    ax3.set_title("Q")
    #plt.subplots_adjust(hspace=.25)
    plt.show()

if __name__ == "__main__":
    model = '1DVar'
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    print('___________________________________\n___________________________________\n___________________'
          '________________\n')
    print('USACE FRF Coastal Model Test Bed : {}'.format(model))

    args = ['./yaml_files/TestBedExampleInputs/1DVar_Input_example.yml']
    try:
        # assume the user gave the path
        yamlLoc = args[0]
        if os.path.exists('.cmtbSettings'):
            with open('.cmtbSettings', 'r') as fid:
                a = yaml.safe_load(fid)
        with open(os.path.join(yamlLoc), 'r') as f:
            inputDict = yaml.safe_load(f)
        inputDict.update(a)

    except:
        raise IOError(
            'Input YAML file required. See yaml_files/TestBedExampleInputs/{}_Input_example for example yaml file.'.format(
                model))

    Master_1DVar_run(inputDict=inputDict)
