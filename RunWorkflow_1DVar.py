# -*- coding: utf-8 -*-
import os, getopt, sys, shutil, glob, logging, yaml, re, pickle
import datetime as DT
from subprocess import check_output
import numpy as np
from getdatatestbed.getDataFRF import getObs, getDataTestBed
from testbedutils import fileHandling
from oct2py import octave
from scipy.io import savemat
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
Q = 0


def assim_currents(X, currents_object, obs, i):
    """Assimilates currents from the FRF-pulled current object into the format to be added to matlab observations matfile
    possible improvements include utilizing averaging instead of just pulling the measurement closest in time to each
    time in dateStringList
    Args:
        X: gridded crosshore distance data
        currents_object: object returned from getWaveSpec in FRFgetData
        obs: matfile containing observation data from example input
        i: simulation number (to full different observation from waves_object at each time step; waves are taken
        every hour, while we might want to populate the observations with measurements taken every three hours.)

    Returns:
        obs with new wave height data and indices written into it from current simulation timestep

    """
    awac_value, awac_indice = find_nearest(X, currents_object['xFRF'])
    obs['v']['d'] = np.append(obs['v']['d'], float(currents_object['aveV'][i]))
    obs['v']['ind'] = np.append(obs['v']['ind'], int(awac_indice))
    obs['v']['e'] = np.append(obs['v']['e'], .1)
    return obs


def assim_waves(X, waves_object, obs, i):
    """Assimilates waves from the FRF-pulled waves objecst into the format to be added to matlab observations matfile.
    possible improvements include utilizing averaging instead of just pulling the measurement closest in time to each
    time in dateStringList

    Args:
        X: gridded crosshore distance data
        waves_object: object returned from getWaveSpec in FRFgetData
        obs: matfile containing observation data from example input
        i: simulation number (to pull different observation from waves_object at each time step; waves are taken
        every hour, while we might want to populate the observations with measurements taken every three hours.)


    Returns:
        obs with new wave height data and indices written into it from current simulation timestep

    """
    hs_value, hs_indice = find_nearest(X, waves_object['xFRF'])
    obs['H']['d'] = np.append(obs['H']['d'], float(waves_object['Hs'][i]))
    obs['H']['ind'] = np.append(obs['H']['ind'], int(hs_indice))
    obs['H']['e'] = np.append(obs['H']['e'], .2)
    return obs


def calculate_Ch(prior, spec8m, X, delta_t, Q):
    """calculate Q(x,t) (measured process error) according to holman et al 2013
    Q = Cq * Hmo ^ (-[(x-x0)/sigma_x]^2)
    deltat is difference from survey to assimilation step;
    add q each time-step to increase uncertainty as the Ch decreases in posterior
    posterior.ch + new S *N *S with just the tiny delta_t

    Args:
        prior: the prior matfile to write Ch to
        waves_obj[-1]: offshore wave spectra object
        X: gridded distance data in crosshore
        delta_t: time difference in seconds between last simulation (or last measured bathy) and current timestep
        Q: defined above

    Returns:
        prior: with newly written Ch = (old Ch) + Q
        Q: process error
    """

    Cq = 0.067  # from sandy duck experiment
    Hmo = np.mean(spec8m['Hs'])  # significant wave height of highest 1/3 of waves
    x0 = 50  # x0 and sigma_x reflect the typical location of breaking waves at this particular beach
    sigma_x = 150
    delta_t = delta_t/(60*60*24)
    print("Delta t in hours:", delta_t*24)
    xx = np.meshgrid(X)
    Lx = 25  # decorrelation length scale, larger Lx leads to smoother results
    N = np.exp((-(xx - np.transpose(xx)) ** 2) / (2 * Lx**2))
    Q = (Cq * Hmo) ** (np.power(-((X - x0) / sigma_x), 2))*delta_t
    S = np.diag(np.sqrt(Q))
    Ch = S * N * S
    if prior['Ch'] == []:
        prior['Ch'] = np.zeros(shape=Ch.shape)
    prior['Ch'] += Ch  # add process error to prior Ch
    return prior, Q


def set_offshore_conditions(prior, waves_obj, element):
    # population the new "prior" during each timestep with FRF data
    prior['theta0'] = np.deg2rad(waves_obj[-1]['waveDm'][element] - 71.8)  # populate with theta0 data from 8m array
    print("Offshore Wave Direction: ", prior['theta0'])
    prior['H0'] = waves_obj[-1]['Hs'][element]  # populate with offshore wave height data from 8m array
    prior['sigma'] = 2 * np.pi * waves_obj[-1]['peakf'][element]  # calculate offshore sigma value from 8m array
    return prior


def create_prior():
    prior = {}
    prior['x'] = []
    prior['h'] = []
    prior['theta0'] = []
    prior['H0'] = []
    prior['ka_drag'] = .015
    prior['hErr'] = []
    prior['Ctheta0'] = .0305
    prior['CH0'] = .01
    prior['Cka'] = 2.5*10**-5
    prior['Ch'] = []
    prior['sigma'] = []
    return prior


def create_obs():
    obs = {}
    obs['H'] = {}
    obs['v'] = {}
    obs['H']['d'] = []
    obs['H']['ind'] = []
    obs['H']['e'] = []
    obs['v']['d'] = []
    obs['v']['ind'] = []
    obs['v']['e'] = []
    obs['tauw'] = []
    obs['tide'] = []
    return obs


def preprocess_data(waves_obj, current_obj, tide_obj, pier_wind, bathyTransect, dateStringList, simulationDuration):
    prior = create_prior()

    # calculate tauw
    speed = pier_wind['windspeed']
    direc = np.deg2rad(pier_wind['winddir'] - 71.8)
    cd = 0.002  # reniers et al. 2004
    tauw = cd * (speed ** 2) * np.cos(direc)
    wind_indice_skip = len(tauw) / len(dateStringList)
    tide_indice_skip = len(tide_obj) / len(dateStringList)

    # restrict bathy to only include points after indice 35,
    zero_elev = np.argmax(bathyTransect['elevation'][200, :] < 0)
    h = bathyTransect['elevation'][200, zero_elev:]
    X = np.linspace(zero_elev, np.max(bathyTransect['xFRF']), len(h))

    # calculate initial delta_t where it is elapsed time since the bathy was measured to first simulation timestep
    delta_t = waves_obj[-1]['time'][0] - bathyTransect['time']
    #delta_t = bathyTransect['time'].timestamp() - waves_obj[-1]['epochtime'][0]
    delta_t = delta_t.seconds
    delta_t = np.abs(delta_t)
    X = np.reshape(X, (len(X), 1))
    h = np.reshape(h, (len(h), 1))
    Q = 0

    # populate prior with data appropriate for this time period
    prior, Q = calculate_Ch(prior, waves_obj[-1], X, delta_t, Q)
    prior['x'] = X  # populate with bathy cross-shore distance data
    prior['h'] = -h  # populate with elevation data
    prior['hErr'] = .1  # populate with elevation error data currently using .1 m

    # populate obs with initial timestep data

    # create lists for storing locations of observations for plotting
    obs_indices_h = []
    obs_indices_v = []
    obs_list = []

    for element, time in enumerate(dateStringList):
        obs = create_obs()
        obs['tauw'] = tauw[0]
        obs['tide'] = tide_obj[0]
        # populate each obs with tauw and tide data
        if element > 0:
            try:
                obs['tauw'] = tauw[int(element * wind_indice_skip)]
                obs['tide'] = tide_obj[int(element * tide_indice_skip)]
            except:
                obs['tauw'] = tauw[-1]
                obs['tide'] = tide_obj[-1]

        for i in range(len(waves_obj)):
            try:
                # grab measurements exactly on times in dateStringList with element*simulationDuration,
                # since waves are taken every hour
                # no averaging done
                obs = assim_waves(X, waves_obj[i], obs, element*simulationDuration)
            except:
                pass

        # add indices with assimilated wave heights to list for plotting
        for obs_timestep_ind in obs['H']['ind']:
            obs_indices_h.append(obs_timestep_ind)

        for i in range(len(current_obj)):
            try:
                # grab measurements exactly on times in dateStringList with element*simulationDuration/3
                # since currents are measured every 3 hours
                # no averaging done
                obs = assim_currents(X, current_obj[i], obs, int(element*simulationDuration/2))
            except:
                pass

        # add indices with assimilated currents to list for plotting
        for obs_timestep_ind in obs['v']['ind']:
            obs_indices_v.append(obs_timestep_ind)

        obs_list = np.append(obs_list, obs)

    return X, h, Q, tauw, wind_indice_skip, tide_indice_skip, delta_t, prior, obs_list, obs_indices_h, obs_indices_v, zero_elev


def plot_1DVar(X, h, Q, posterior, obs_indices_h, obs_indices_v, projectEnd, zero_elev):
    # plot prior, posterior, difference, and Q
    # plot hErr on the prior and posterior
    # show areas with assimilated data on all plots

    #grab bathy nearest enddate
    bathygo = getDataTestBed(projectEnd, projectEnd)
    bathyTransect = bathygo.getBathyIntegratedTransect(method=1)  # grab bathymetry transects
    print("Final Bathy Survey Date: ", bathyTransect['time'])
    survey_h = bathyTransect['elevation'][200, zero_elev:]

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
    ax2.scatter(X, -(h+posterior['h']), color="black")
    ax2.set_xlabel("Crosshore distance (m)")
    ax2.set_ylabel("Elevation change (m)")

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
    ax3.plot(X, survey_h, color='cyan', label="Survey")
    ax3.plot(X, -posterior['h'], color='red', label="Posterior")
    ax3.set_title("Posterior and Survey\n" + str(bathyTransect['time']))
    for obs_loc in obs_indices_h:
        ax3.axvline(x=obs_loc*(X[-1]/len(h)), color='r', linestyle='--')
        temp = obs_loc
    ax3.axvline(x=temp * (X[-1] / len(h)), color='r', linestyle='--', label="wave heights")
    for obs_loc in obs_indices_v:
        ax3.axvline(x=obs_loc*(X[-1]/len(h)), color='orange', linestyle='--')
        temp = obs_loc
    ax3.axvline(x=temp * (X[-1] / len(h)), color='orange', linestyle='--', label="currents")
    ax3.legend()
    #plt.subplots_adjust(hspace=.25)
    plt.show()


def find_nearest(array, value):
    """This function will run CMS with any version prefix given start, end, and timestep.

    Args:
      inputDict: a dictionary that is read from the input yaml

    Returns:
      None

    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


def Master_1DVar_run(inputDict):
    """This function will run 1DVar given start, end, and timestep found in input yaml

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
    matlabfiledir = inputDict['mainDirectory']
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

        # initiliaze get observation for bathymetry
        bathygo = getDataTestBed(projectStart, projectEnd)
        bathyTransect = bathygo.getBathyIntegratedTransect(method=1)  # grab bathymetry transects
        print("Initial Bathy Survey Date: ", bathyTransect['time'])

        # initialize get observation for measurements
        go = getObs(projectStart, projectEnd)

        # pull wind
        pier_wind = None
        while pier_wind == None:
            pier_wind = go.getWind()

        # pull water level
        tides_list = []
        tides_list = np.append(tides_list, go.getWL()['WL'])

        # pull currents
        current_gauges = ['adop-3.5m', 'awac-4.5m', 'awac-6m', 'awac-8m', 'awac-11m']
        current_obj = []
        for i in range(len(current_gauges)):
            try:
                current_obj = np.append(current_obj, go.getCurrents(gaugenumber=current_gauges[i], roundto=1))
            except:
                print("NetCDF DAP Error while grabbing current data from gauge: ", current_gauges[i])

        # pull waves
        #wave_gauges = ['adop-3.5m', 'awac-4.5m', 'awac-6m', 'awac-8m', 'awac-11m', 'xp100m', 'xp125m', 'xp150m', 'xp200m', '8m-array']
        wave_gauges = ['xp100m', 'xp125m', 'xp150m', 'xp200m', '8m-array']
        waves_obj = []
        for i in range(len(wave_gauges)):
            waves_obj = np.append(waves_obj, go.getWaveData(gaugenumber=wave_gauges[i], spec=False))

        #preprocess data
        X, h, Q, tauw, wind_indice_skip, tide_indice_skip, delta_t, prior, obs_list, obs_indices_h, obs_indices_v, zero_elev = \
            preprocess_data(waves_obj, current_obj, tides_list, pier_wind, bathyTransect, dateStringList, simulationDuration)

        obs_dict = {}
        obs_dict["struct"] = obs_list
        savemat("./data/matlab_files/initialprior.mat", prior)
        savemat("./data/matlab_files/obslist.mat", obs_dict)

    # ________________________________________________ RUN LOOP ________________________________________________
    for element, time in enumerate(dateStringList):
        try:
            print('-------------------------------Beginning Simulation {}-------------------------------'.format(
                DT.datetime.now()))
            
            if runFlag == True:  # run model
                os.chdir(modeldir)  # changing locations to where the model will be ran from for octave
                print('Running {} Simulation'.format(model.upper()))
                dt = DT.datetime.now()

                # set offshore wave conditions used in prior based on new timestop of obs
                # grab measurements exactly on times in dateStringList with element*simulationDuration,
                # since waves are taken every hour
                prior = set_offshore_conditions(prior, waves_obj, element*simulationDuration)

                #grab obs data for current timestep
                obs = obs_list[element]
                print("Assimilation no. " + str(element) + "/" + str(len(dateStringList)))
                print("Assimilating Date: ", dateStringList[element])
                print("Wave heights (m) to assimilate: ", obs['H']['d'])
                print("Wave height indices: ", obs['H']['ind'])
                print("Currents (m/s) to assimilate: ", obs['v']['d'])
                print("Current indices: ", obs['v']['ind'])

                offshore_indice = int(obs['v']['ind'][-1])
                print(offshore_indice)
                if offshore_indice < 242:
                    print("No offshore forcing present in current measurements for this timestep -- skipping assimilation")
                    delta_t = waves_obj[-1]['time'][int(element + 2*simulationDuration)] - waves_obj[-1]['time'][element]
                    delta_t = np.abs(delta_t.seconds)
                    prior, Q = calculate_Ch(prior, waves_obj[-1], X, delta_t, Q)
                    continue
                print("Tauw to assimilate: ", obs['tauw'])
                print("Tide: ", obs['tide'])


                # run 1DVar model with prior and observation input
                # set nout if desired to obtain diagnostics and representer matrices.
                # save matfile of this particular prior and obs
                os.chdir(matlabfiledir)
                savemat("./data/matlab_files/prior_" + str(time[:13]) + ".mat", prior)
                savemat("./data/matlab_files/obs_" + str(time[:13]) + ".mat", obs)
                posterior, diagnostics = octave.feval(modeldir + '/assim_1dh.m', prior, obs, 1, nout=2)#, representers = \
                print('Simulation took %s ' % (DT.datetime.now() - dt))
            os.chdir(matlabfiledir)
            # save matfile of the posterior result
            savemat("./data/matlab_files/posterior_" + str(time[:13]) + ".mat", posterior)
            # make the output of the model the prior for the next timestep
            prior = posterior

            # find delta_t for the next timestep
            try:
                delta_t = waves_obj[-1]['time'][int(element+simulationDuration)] - waves_obj[-1]['time'][element]
                delta_t = np.abs(delta_t.seconds)
            except:
                pass

            #calculate Ch for the next timestep
            prior, Q = calculate_Ch(prior, waves_obj[-1], X, delta_t, Q)

            print('-------------------------------SUCCESS-----------------------------------------')

        except Exception as e:
            print('<< ERROR >> HAPPENED IN THIS TIME STEP\n{}'.format(e))
            logging.exception('\nERROR FOUND @ {}\n'.format(time), exc_info=True)
            os.chdir(modeldir)

    plot_1DVar(X, h, Q, posterior, obs_indices_h, obs_indices_v, projectEnd, zero_elev)

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
