##################################
# written by: Justus Kromer
##################################
import sys
import numpy as np
import os

# specify directory in which data will be saved
dataDirectory = "data"

# list of seeds used to generate synaptic connectivity and noise
seed_array = [10, 12, 14, 16, 18]



######################################################################
###  Simulations for CRS for intermediate networks
######################################################################
## A) use the following to simulate CRS with different Tshuffle
#
# 			python 4_run.py CRS_Tshuffle_intermediate_networks
#
if sys.argv[1] == 'CRS_Tshuffle_intermediate_networks':

	NETWORK_TYPE = "intermediate_network"

	# path to simulation script
	pathSherlockHome           = '../CRS'
	# output directory
	pathSherlockScratch        = dataDirectory + '/'+str(NETWORK_TYPE)+'/CRS_Tshuffle'

	# the following output file was used for submission to computing cluster
	# listOfJobsThatNeedToBeDone = 'CRS_Tshuffle_intermediate_network.txt'

	# stimulation time
	Tstim = 7220 # sec
	# time at which stimulation starts
	TstartSim = 3000 # sec
	# stimulation frequency
	fCRarray = [10.0] # Hz
	# e pulse scale (currently not implemented)
	e_pulse_scale = 1.0 # ms
	# number of pulses per stimulus
	number_of_pulses_per_burst = [1]
	# stimulation amplitude
	AstimArray = [1.0]
	# f intra
	fintra = 130.0 # Hz
	# number of separately stimulated subpopulations
	nElectrodeArray = [4]
	# distance between stimulation sites in units of system length scale L=5.0 mm
	dsites = 0.25
	# effective rectangular width of stimulus profile in units of dsites
	erwArray = [0.5] 
	# values of the heterogeneity parameter (lambPar = H)
	lambPar_array = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

	# array of seeds for CR sequence realizations
	seedSequence_array = np.arange(10,300,30).astype(int)
	# Tshuffle periods in seconds
	Tshuffle_sec_array = [0.1, 10.0, 1800.0]


	jobcounter = 0

	# run through seeds
	for seed in seed_array:
		for seedSequence in seedSequence_array:
			for lambPar in lambPar_array:
				for Tshuffle_sec in Tshuffle_sec_array:
					for number_of_pulse in number_of_pulses_per_burst:
						for fCR in fCRarray:
							for Astim in AstimArray:
								for M in nElectrodeArray:
									for erw in erwArray:

										# specify backup directory from which simulation is started
										initalConditionsDirectory = 'data/'+str(NETWORK_TYPE)+'/initial_network/lambPar_'+str(lambPar)+'_seed_'+str(seed)+'/'+str(TstartSim)+'_sec'

										# output directory in which results are saved
										outputdirectory = pathSherlockScratch+'/lambPar_'+str(lambPar)+'_seed_'+str(seed)

										## check whether simulation was already performed
										# by looking for the file "10200_sec/cMatrix.npz"
										finalBackupFilename = outputdirectory + "/SVS_burst_sequence_seed/fCR_"+str(fCR)+"_M_"+str(M)+"_fintra_"+str(fintra)+"_npb_"+str(number_of_pulse)+"_Tshuffle_sec_"+str(Tshuffle_sec)+"_seedSequence_"+str(seedSequence)+"_Astim_"+str(Astim)+"_Tstim_"+str(float(Tstim))+"_dsites_"+str(dsites)+"_erw_"+str(erw)+"/10200_sec/cMatrix.npz"
										
										# if simulation was already performed
										if os.path.isfile( finalBackupFilename ) :               
											# simulation has already been performed 
											print("simulation results found")
										else:
											# simulation needs to be performed
											print("create submission string")
											jobcounter+=1

											# submission string
											shellCommand = 'python '+pathSherlockHome+'/CRS_Tshuffle.py '+ initalConditionsDirectory + ' ' + str(Astim) + ' ' + str(fCR) + ' ' + str(number_of_pulse) + ' ' + str(Tstim) + ' ' + str(M) + ' ' + str(e_pulse_scale) + ' ' + str(fintra) + ' ' + outputdirectory + ' ' + str(seedSequence) +' '+ str(Tshuffle_sec) + ' ' + str(dsites) + ' ' + str(erw) + '\n'

											# print shell command
											print( shellCommand )

											# the following code was used for submission to computing cluster
											# submissionString = 'srun ' + shellCommand

											# # adds job to end of list
											# with open(listOfJobsThatNeedToBeDone, "a") as myfile:
																			
											# 	myfile.write( submissionString )
											# 	myfile.write( '##\n' )
	print(jobcounter)


#################################################################################################
###  Simulations for relaxation after CRS for intermediate networks
#################################################################################################
## C) use the following to simulate relaxation after CRS with different Tshuffle
#
# 			python 4_run.py relaxation_after_shuffled_CRS_intermediate_networks
#
if sys.argv[1] == 'relaxation_after_shuffled_CRS_intermediate_networks':

	NETWORK_TYPE = "intermediate_network"

	# path to simulation script
	pathSherlockCRstim        = dataDirectory + '/'+str(NETWORK_TYPE)+'/CRS_Tshuffle'
	# path to less_relaxation.py
	pathSherlockHome_relax     = "../relaxation_after_stimulation/"
	# the following output file was used for submission to computing cluster
	# listOfJobsThatNeedToBeDone = 'relaxation_shuffled_CRS_network.txt'

	# stimulation time
	Tstim = 7220 # sec
	# time at which stimulation starts
	TstartSim = 3000 # sec
	# stimulation frequency
	fCRarray = [10.0] # Hz
	# e pulse scale (currently not implemented)
	e_pulse_scale = 1.0 # ms
	# number of pulses
	number_of_pulses_per_burst = [1]
	# stimulation amplitude
	AstimArray = [1.0]
	# f intra
	fintra = 130.0 # Hz
	# number of separately stimulated subpopulations
	nElectrodeArray = [4]
	# distance between stimulation sites in units of L
	dsites = 0.25
	# effective rectangular width of stimulus profile in units of dsites
	erwArray = [0.5] 
	# values of the heterogeneity parameter (lambPar = H)
	lambPar_array = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

	# time after which relaxation is started (this is the time at which stimulation ceases)
	TstartRelax = 10200 # sec
	# relaxation time (as a default value this is set to 3 hours and 20 sec). We added 220 sec to make sure that a backup is generated.
	Trelax = 10820+220 # sec
	# Time of backup 3 hours after cessation of stimulation relaxation. 
	TafterRelax = TstartRelax+Trelax # sec


	# array of seeds for CR sequence realizations
	seedSequence_array = np.arange(0,300,10).astype(int)
	# list of shuffle periods for which CRS was simulated
	Tshuffle_sec_array = [0.1, 10.0, 1800.0] # sec

	jobcounter = 0

	# run through seeds
	for seed in seed_array:
		for seedSequence in seedSequence_array:
			for lambPar in lambPar_array:
				for Tshuffle_sec in Tshuffle_sec_array:
					for number_of_pulse in number_of_pulses_per_burst:
						for fCR in fCRarray:
							for Astim in AstimArray:
								for M in nElectrodeArray:
									for erw in erwArray:

										# build path to backup directory from which simulation is started
										outputdirectory     = pathSherlockCRstim +'/lambPar_'+str(lambPar)+'_seed_'+str(seed)
										# path to output of CRS simulations
										pathOutputCRS = outputdirectory + "/SVS_burst_sequence_seed/fCR_"+str(fCR)+"_M_"+str(M)+"_fintra_"+str(fintra)+"_npb_"+str(number_of_pulse)+"_Tshuffle_sec_"+str(Tshuffle_sec)+"_seedSequence_"+str(seedSequence)+"_Astim_"+str(Astim)+"_Tstim_"+str(float(Tstim))+"_seedSeq_"+str(seedSequence)+"_dsites_"+str(dsites)+"_erw_"+str(erw)
										# this is the backup folder fromm which the relaxation simulation is started
										backupDirectory_results_stim    = pathOutputCRS + '/'+str(int(TstartRelax))+'_sec'
										
										# output directory
										outputDirectory_rel_after_stim =  pathOutputCRS + '/relax_after_'+str(int(TstartRelax))+'_sec'

										## check whether CRS simulation is already done
										# this is done by looking for the cMatrix file at time 10200 sec
										finalBackupFilename = pathOutputCRS + "/"+str(TstartRelax)+"_sec/cMatrix.npz"
										
										# print(finalBackupFilename)
										
										# check whether backup directory exists
										if os.path.exists( finalBackupFilename ):

											print( 'backup found' )

											# build submission string
											shellCommand =     'python '+str(pathSherlockHome_relax)+'less_relaxation.py '+ backupDirectory_results_stim + ' ' + str(int(Trelax)) + ' ' + outputDirectory_rel_after_stim + '\n'

											# print shell command
											print( shellCommand )

											# the following code was used for submission to computing cluster
											# submissionString = 'srun ' + shellCommand

											# # adds job to end of list
											# with open(listOfJobsThatNeedToBeDone, "a") as myfile:
																			
											# 	myfile.write( submissionString )
											# 	myfile.write( "##\n" )
												
										else:
											print( 'backup not found' )


