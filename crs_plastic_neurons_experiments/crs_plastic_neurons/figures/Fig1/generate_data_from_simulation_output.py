##################################
# written by: Justus Kromer
##################################

import numpy as np
import sys
import os 
from scipy.interpolate import interp1d
import scipy.sparse

def loadListOfBackupTimesteps( Directory ):
    
    lines=[]
    backupSWithDuplicates=[]

    with open( Directory+'/listOfBackupTimeSteps.txt','r') as f:
        lines.append( f.read() )
        lines=(lines[0].split('\n'))
        for line in lines:
            backupSWithDuplicates.append(line.split(' '))

        backupS = []
        for i in backupSWithDuplicates:
            if i not in backupS:
                backupS.append(i)
    return backupS

# loads complete spike train and weight sequence for order of listed directories
def load_Complete_SpikeTrain_And_Weigth_Trajectories( sortedListOfDirectories ):
    
    import os

    # initialize spikeTrain and weightData
    spikeTrain=[]
    weightData=[]
                
    for Directory in sortedListOfDirectories:
        # load first Trajectory
        if os.path.isfile( Directory+'/listOfBackupTimeSteps.txt' ): 
            # load short data
            # load long data
            backupS =loadListOfBackupTimesteps( Directory )
    
            # load data
            for fileEndings in range(len(backupS)-1):
                
                # get actual file ending
                fileEndingString=backupS[ fileEndings ][1]
                if fileEndingString!='0_sec':
                    
                    #print fileEndingString
                    spikeTimesFile=Directory+'/spikeTimes_'+str(fileEndingString)+'.npy'
                        
                    if os.path.isfile( spikeTimesFile ):
                        if len(spikeTrain)==0:
                            spikeTrain=np.load( spikeTimesFile )
                        else:

                            # ensure that spike trains dont overlap
                            newSpikeTrain = np.load( spikeTimesFile )
                            if len(newSpikeTrain) != 0:
                                #print len(newSpikeTrain), len(spikeTrain), spikeTrain[-1,1], newSpikeTrain[0,1]
                                oldSpikeTrain = spikeTrain[ spikeTrain[:,1]<newSpikeTrain[0,1] ]

                                spikeTrain=np.concatenate( ( oldSpikeTrain, newSpikeTrain ), axis=0 )
                    
                    weightFile=Directory+'/meanWeightTimeSeries_'+str(fileEndingString)+'.npy'
                    if os.path.isfile( weightFile ):
                        if len(weightData)==0:
                            weightData=np.load( weightFile )
                        else:
                            weightData=np.concatenate( ( weightData,np.load(weightFile ) ), axis=0 )
    
        else:
            print( 'Error: no data found in', Directory )
        

    # returns spike train and weight trajectory. 
    # times are in simulation time steps
    return spikeTrain, weightData        

#### calculates Kuramoto order parameter
# resolution ... temporal distance between  (ms)
def piece_wise_calcKuramotoOrderParameter( spikeTimes, tmin, tmax, resolution, arrayOfNeuronIndixes, outputFilename ):

    # delete empty entries
    spikeTimes=spikeTimes[ spikeTimes[:,1]!= 0 ]

    ################################################################################
    ######## the following two lines were added later ##############################
    populationSize=len(arrayOfNeuronIndixes)
    
    # number of grid points for which Kuramoto order parameter is evaluated at once
    NinterPolSteps = int( 1000.0*( tmax-tmin )/float(resolution) )
    processAtOnce_NinterPolSteps = 100000

    # if too many gridpoints are given a piece-wise calculation is performed
    if NinterPolSteps > processAtOnce_NinterPolSteps:
        
        KuramotoOutArray = []
        
        currentInterPolSteps = 0
        lengthOfTimeIntervals = resolution * processAtOnce_NinterPolSteps # ms
        # in order to exclude boundary effects when combining arrays, we consider an overlap of 2000 ms
        TimeStepsOfOverlap = int( 2000.0/resolution ) # time steps

        if 2*TimeStepsOfOverlap > processAtOnce_NinterPolSteps:
            print( 'ERROR: overlap for piece-wise Kuramoto order parameter calculation too long compared to processAtOnce_NinterPolSteps!' )
            return 0
        
        current_Tmin = 1000.0*tmin # ms
        current_Tmax = current_Tmin + lengthOfTimeIntervals # ms
        
        while current_Tmax < 1000*tmax:

            # initialize phases at fixed points
            phases=np.zeros( (populationSize, processAtOnce_NinterPolSteps) )
            arrayOfGridPoints=np.arange( current_Tmin, current_Tmax , resolution)

            # consider only spikes that are between tmin and tmax
            processedSpikeTimes=spikeTimes[ np.logical_and( spikeTimes[:,1]>=current_Tmin , spikeTimes[:,1]<= current_Tmax )  ]        

            # calculate phases in current interval
            krecPhases=0
            for kNeuron in arrayOfNeuronIndixes:

                # get spike train of corresponding neuron
                spikeTrainTemp = processedSpikeTimes[processedSpikeTimes[:,0].astype(int)== kNeuron ][:,1]

                # calc phase function
                if len(spikeTrainTemp) != 0:
                    phaseNPiCrossings=np.concatenate( ( np.full( 1, 1000.0*(2*tmin-tmax) ) , spikeTrainTemp, np.full( 1, 1000.0*(2*tmax-tmin) ) ), axis=0 )
                    PhaseValues=np.linspace(0,len(phaseNPiCrossings)-1, len(phaseNPiCrossings))
                else:
                    phaseNPiCrossings=np.array([-np.inf, np.inf])
                    PhaseValues=np.array([0, 1])

                # linear interpolate phaseNPiCrossings
                phaseFunctionKNeuron=interp1d(phaseNPiCrossings,2*np.pi*PhaseValues)

                phases[krecPhases,:]=phaseFunctionKNeuron(arrayOfGridPoints) 
                krecPhases+=1

            # calc Kuramoto order parameter
            TotalArrayOfKuramotoOrderParameterAtGridPoints=1/float(populationSize)*np.absolute(np.sum( np.exp( 1j*phases ), axis=0 ))
            current_KuramotoOutArray=np.array( [arrayOfGridPoints, TotalArrayOfKuramotoOrderParameterAtGridPoints] )

            current_KuramotoOutArray=np.transpose( current_KuramotoOutArray )

            if len(KuramotoOutArray) == 0:
                KuramotoOutArray = current_KuramotoOutArray[:-TimeStepsOfOverlap]
            else:
                # add to previous Kuramoto array
                KuramotoOutArray = np.concatenate( ( KuramotoOutArray, current_KuramotoOutArray[ TimeStepsOfOverlap:-TimeStepsOfOverlap ] ), axis = 0 )

            # prepare boundaries of next interval 
            current_Tmin = current_Tmax - resolution * 2*TimeStepsOfOverlap # ms
            current_Tmax = current_Tmin + lengthOfTimeIntervals # ms
     
            print( 'current interval (left boundary)' , current_Tmin * 0.001 )
                
                
    # calculate in one single run    
    else:
        phases=np.zeros( (populationSize, NinterPolSteps) )
        arrayOfGridPoints=np.linspace(1000.0*tmin,1000.0*tmax,NinterPolSteps)

        processedSpikeTimes=spikeTimes[ np.logical_and( spikeTimes[:,1]>=1000.0*tmin , spikeTimes[:,1]<= 1000.0*tmax )  ]

            
        krecPhases=0
        for kNeuron in arrayOfNeuronIndixes:

            # get spike train of corresponding neuron
            spikeTrainTemp=processedSpikeTimes[processedSpikeTimes[:,0].astype(int)== kNeuron ][:,1]
            # calc phase function
            if len(spikeTrainTemp) != 0:
                phaseNPiCrossings=np.concatenate( ( np.full( 1, 1000.0*(2*tmin-tmax) ) , spikeTrainTemp, np.full( 1, 1000.0*(2*tmax-tmin) ) ), axis=0 )
                PhaseValues=np.linspace(0,len(phaseNPiCrossings)-1, len(phaseNPiCrossings))
            else:
                phaseNPiCrossings=np.array([-np.inf, np.inf])
                PhaseValues=np.array([0, 1])

            # linear interpolate phaseNPiCrossings
            phaseFunctionKNeuron=interp1d(phaseNPiCrossings,2*np.pi*PhaseValues)

            phases[krecPhases,:]=phaseFunctionKNeuron(arrayOfGridPoints) 
            krecPhases+=1

        # calc Kuramoto order parameter
        TotalArrayOfKuramotoOrderParameterAtGridPoints=1/float(populationSize)*np.absolute(np.sum( np.exp( 1j*phases ), axis=0 ))
        KuramotoOutArray=np.array( [arrayOfGridPoints, TotalArrayOfKuramotoOrderParameterAtGridPoints] )
        
        KuramotoOutArray=np.transpose( KuramotoOutArray )

        
    if outputFilename!='':
       np.save( outputFilename+'.npy' , KuramotoOutArray )
    
    return KuramotoOutArray


###############################################################
# use the following to load 20sec long spike train and mean weight files 
# from simulation output and generate "homogeneous_mwTrajectory_w_MW_seed_SEED.npy"
# and "homogeneous_spikeTrain_w_MW_seed_SEED.npy"
###############################################################
### for sys.argv[2] use
#   homogeneous_network
#   inhomogeneous_network
if sys.argv[1] == "get_traj_meanWeight_spikeTrain_multistability":

    seed = 12
    # specify list of mean initial synaptic weights for which trajectories are loaded
    mwInitValues = [0.0, 0.2, 0.4, 0.6]# np.round( np.arange(0.0,1.1,0.1) , 1 )

    for mw in mwInitValues:
        
        print( mw )

        # The following is the default file path to the simulation output. 
        # Adjust this path if another outputDirectory was used in 2_run.py.
        directory = "../../run_sim/data/"+sys.argv[2]+"/multistability/seed_"+str(seed)+"/mw_"+str(mw)

        # load spike train and weight data
        spikeTrain, weightData = load_Complete_SpikeTrain_And_Weigth_Trajectories( [directory] )

        # load adjaceny matrix
        filename = directory + "/0_sec/synConnections.npz"
        adj = scipy.sparse.load_npz( filename )
        
        # The following line rescales the second column of weightData (mean weight for all NxN entries in weight matrix)
        # by the mean of the adjacency matrix such that the second column of mwTrajectory contains 
        # the mean weight of synaptic connecitons.  
        mwTrajectory = weightData[:,:2]*[[1.0,1.0/np.mean( adj[:1000,:1000] )]]
        
        # the following is the directory in which data are saved
        np.save( "data/"+sys.argv[2]+"_spikeTrain_w_"+str(mw)+"_seed_"+str(seed)+".npy" , spikeTrain )
        np.save( "data/"+sys.argv[2]+"_mwTrajectory_w_"+str(mw)+"_seed_"+str(seed)+".npy" , mwTrajectory )


###############################################################
###############################################################
if sys.argv[1] == "calculate_time_trace_of_Kuramoto_parameter_multistability_mw":

    # time interval of interest
    tmin =  0.0 # sec
    tmax = 13000.0 # sec

    resolution = 1000.0 # ms (one data point every "resolution" ms)

    arrayOfNeuronIndixes = np.arange( 1000 ).astype( int )

    # generate data file
    seed = 12
    mwInitValues = np.round( np.arange(0.0,1.1,0.1) , 1 )
    dataFiles = []

    # mwInitValues = [float(sys.argv[2])]
    mwInitValues = [0.0, 0.2, 0.4, 0.6]#

    # loop over all results
    for mw in mwInitValues:
        for networkType in ["inhomogeneous_network","homogeneous_network"]:
            print( mw, networkType )

            # load spike train
            full_spike_train = "data/"+networkType+"_spikeTrain_w_"+str(mw)+"_seed_"+str(seed)+".npy"
            
            if os.path.isfile( full_spike_train ):
                spikeTimes = np.load( full_spike_train )

                # output folder Kuramoto
                outputFilename = "data/KuramotoOrderParameter_"+networkType+"_mw_"+str(mw)+"_seed_"+str(seed)

                # calculate Kuramoto order parameter and save it to "outputFilename"
                piece_wise_calcKuramotoOrderParameter( spikeTimes*[1,0.1], tmin, tmax, resolution, arrayOfNeuronIndixes, outputFilename )

            else:
                print( "data not found", full_spike_train )

        