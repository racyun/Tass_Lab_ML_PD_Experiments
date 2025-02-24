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
###  Simulations for CRS for homogeneous and inhomogeneous networks
######################################################################
## A) use the following to simulate CRS with different Tshuffle
#
# 			python 3_run.py CRS_Tshuffle
#
if sys.argv[1] == 'CRS_Tshuffle':

	NETWORK_TYPES = ["homogeneous_network","inhomogeneous_network"]

	for NETWORK_TYPE in NETWORK_TYPES:

		# path to simulation script
		pathSherlockHome           = '../CRS'
		# output directory
		pathSherlockScratch        = dataDirectory + '/'+str(NETWORK_TYPE)+'/CRS_Tshuffle'

		# the following output file was used for submission to computing cluster
		# listOfJobsThatNeedToBeDone = 'CRS_network.txt'

		# stimulation time
		Tstim = 7220 # sec
		# time at which stimulation starts
		TstartSim = 3000 # sec
		# stimulation frequency
		fCRarray = [10.0] # Hz
		# pulse width (currently not implemented)
		e_pulse_scale = 1.0 # ms
		# number of pulses per stimulus
		number_of_pulses_per_burst = [1]
		# stimulation amplitude
		AstimArray = [1.0]
		# intraburst frequency
		fintra = 130.0 # Hz
		# number of separately stimulated subpopulations (only 4 is implemented)
		nElectrodeArray = [4]
		# distance between stimulation sites in units of system length scale L=5.0 mm
		dsites = 0.25
		# effective rectangular width of stimulus profile in units of dsites
		erwArray = [0.5] 

		# array of seeds for CR sequence realizations
		seedSequence_array = np.arange(0,300,10).astype(int)
		# Tshuffle periods in seconds
		Tshuffle_sec_array = [0.1, 10.0, 1800.0] # sec

		jobcounter = 0

		# run through seeds
		for seed in seed_array:
			for seedSequence in seedSequence_array:
				for Tshuffle_sec in Tshuffle_sec_array:
					for number_of_pulse in number_of_pulses_per_burst:
						for fCR in fCRarray:
							for Astim in AstimArray:
								for M in nElectrodeArray:
									for erw in erwArray:

										# specify backup directory from which simulation is started
										initalConditionsDirectory = 'data/'+str(NETWORK_TYPE)+'/initial_network/seed_'+str(seed)+'/'+str(TstartSim)+'_sec'

										# output directory in which results are saved
										outputdirectory = pathSherlockScratch+'/seed_'+str(seed) 

										## Check whether simulation has already been performed
										# by looking for the file "10200_sec/cMatrix.npz" and "meanWeightTimeSeries_3020_sec.npy"
										finalBackupFilename = outputdirectory + "/SVS_burst_sequence_seed/fCR_"+str(fCR)+"_M_"+str(M)+"_fintra_"+str(fintra)+"_npb_"+str(number_of_pulse)+"_Tshuffle_sec_"+str(Tshuffle_sec)+"_seedSequence_"+str(seedSequence)+"_Astim_"+str(Astim)+"_Tstim_"+str(float(Tstim))+"_seedSeq_"+str(seedSequence)+"_dsites_"+str(dsites)+"_erw_"+str(erw)+"/10200_sec/cMatrix.npz"
										firstBackupFilename = outputdirectory + "/SVS_burst_sequence_seed/fCR_"+str(fCR)+"_M_"+str(M)+"_fintra_"+str(fintra)+"_npb_"+str(number_of_pulse)+"_Tshuffle_sec_"+str(Tshuffle_sec)+"_seedSequence_"+str(seedSequence)+"_Astim_"+str(Astim)+"_Tstim_"+str(float(Tstim))+"_seedSeq_"+str(seedSequence)+"_dsites_"+str(dsites)+"_erw_"+str(erw)+"/meanWeightTimeSeries_3020_sec.npy"
		
										# when simulations were already done
										if os.path.isfile( finalBackupFilename ) and os.path.isfile( firstBackupFilename ):               
											print(seed, seedSequence, Tshuffle_sec)
											# simulation has already been perfomed earlier
											print("simulation results found")
											print( firstBackupFilename )
											
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
									
		print( jobcounter )


## B) use the following to simulate non-shuffled CRS with different CR sequences
#
# 			python 3_run.py CRS_non_shuffled
#
if sys.argv[1] == 'CRS_non_shuffled':

	NETWORK_TYPES = ["homogeneous_network","inhomogeneous_network"]

	for NETWORK_TYPE in NETWORK_TYPES:

		# path to simulation script
		pathSherlockHome           = '../CRS'
		# output directory
		pathSherlockScratch        = dataDirectory + '/'+str(NETWORK_TYPE)+'/non_shuffled_CRS'

		# the following output file was used for submission to computing cluster
		# listOfJobsThatNeedToBeDone = 'non_shuffled_CRS_network.txt'

		# stimulation time
		Tstim = 7220 # sec
		# time at which stimulation starts
		TstartSim = 3000 # sec
		# stimulation frequency
		fCRarray = [10.0] # Hz
		# stimulation amplitude
		AstimArray = [1.0]
		# distance between stimulation sites in units of system length scale L=5.0 mm
		dsites = 0.25
		# number of stimulation sites (only 4 is implemented)
		nElectrodeArray = [4]
		# effective rectangular width of stimulus profile in units of dsites
		erwArray = [0.5] 
				
		# stimulation sequences
		seedSequenceArray = ['0_1_2_3','0_1_3_2','0_2_1_3','0_2_3_1','0_3_1_2','0_3_2_1']

		# run through parameter combinations
		for seed in seed_array:
			for sequence in seedSequenceArray:
				for fCR in fCRarray:
					for Astim in AstimArray:		
						for M in nElectrodeArray:
							for erw in erwArray:

								# This is the backup directory from which the simulation is started.
								initalConditionsDirectory = 'data/'+str(NETWORK_TYPE)+'/initial_network/seed_'+str(seed)+'/'+str(TstartSim)+'_sec'

								# output directory in which results are saved
								outputdirectory = pathSherlockScratch+'/seed_'+str(seed) 

								# submission string
								shellCommand = 'python '+ pathSherlockHome +'/CRS_non_shuffled.py '+ initalConditionsDirectory+' '+    str(Astim)+' '+    str(fCR)+' '+    str(1)+' '+                  str(Tstim)+' '+    str(M)+' '+sequence+' '+outputdirectory+' '+str(dsites) + ' ' + str(erw) + '\n'

								# print shell command
								print( shellCommand )

								# the following code was used for submission to computing cluster
								# submissionString = 'srun ' + shellCommand

								# # adds job to end of list
								# with open(listOfJobsThatNeedToBeDone, "a") as myfile:
																
								# 	myfile.write( submissionString )
								# 	myfile.write( '##\n' )

#################################################################################################
###  Simulations for relaxation after CRS for homogeneous and inhomogeneous networks
#################################################################################################
## C) use the following to simulate relaxation after CRS with different Tshuffle
#
# 			python 3_run.py relaxation_after_shuffled_CRS
#
if sys.argv[1] == 'relaxation_after_shuffled_CRS':

	NETWORK_TYPES = ["homogeneous_network","inhomogeneous_network"]

	for NETWORK_TYPE in NETWORK_TYPES:

		# Directory in which results of simulations of shuffled CRS are saved.
		pathSherlockScratch        = dataDirectory + '/'+str(NETWORK_TYPE)+'/CRS_Tshuffle'
		# This is the directory which "less_relaxation.py" is located.
		pathSherlockHome_relax     = "../relaxation_after_stimulation/"
		# the following output file was used for submission to computing cluster
		# listOfJobsThatNeedToBeDone = 'relaxation_shuffled_CRS_network.txt'

		# stimulation time
		Tstim = 7220 # sec
		# time at which stimulation starts
		TstartSim = 3000 # sec
		# stimulation frequency
		fCRarray = [10.0] # Hz
		# e pulse scale (not implemented)
		e_pulse_scale = 1.0 # ms
		# number of pulses
		number_of_pulses_per_burst = [1]
		# stimulation amplitude
		AstimArray = [1.0]
		# intraburst frequency
		fintra = 130.0 # Hz
		# distance between stimulation sites in units of system length scale L=5.0 mm 
		dsites = 0.25
		# number of separately stimulated subpopulations
		nElectrodeArray = [4]
		# effective rectangular width of stimulus profile in units of dsites
		erwArray = [0.5] 


		# time after which relaxation is started (this is the time at which stimulation ceases)
		TstartRelax = 10200 # sec
		# relaxation time (as a default value this is set to 3 hours and 20 sec). We added 220 sec to make sure that a backup is generated.
		Trelax = 10820+220 # sec
		# Time of backup 3 hours after cessation of stimulation relaxation. 
		TafterRelax = TstartRelax+Trelax # sec

		# array of seeds for CR sequence realizations
		seedSequence_array = np.arange(0,300,10).astype(int)
		# Tshuffle periods in seconds
		Tshuffle_sec_array = [0.1, 10.0, 1800.0] 

		# run through seeds
		for seed in seed_array:
			for seedSequence in seedSequence_array:
				for Tshuffle_sec in Tshuffle_sec_array:
					for number_of_pulse in number_of_pulses_per_burst:
						for fCR in fCRarray:
							for Astim in AstimArray:
								for M in nElectrodeArray:
									for erw in erwArray:

										# directory in which backup from stimulation is stored
										parentDirectory = pathSherlockScratch + '/seed_'+str(seed)+'/SVS_burst_sequence_seed/fCR_'+str(fCR)+'_M_'+str(M)+'_fintra_'+str(fintra)+'_npb_'+str(number_of_pulse)+'_Tshuffle_sec_'+str(Tshuffle_sec)+'_seedSequence_'+str(seedSequence)+'_Astim_'+str(Astim)+'_Tstim_'+str(float(Tstim))+'_seedSeq_'+str(seedSequence)+'_dsites_'+str(dsites)+'_erw_'+str(erw)	

										## This is the backup directory from which simulations are started. 
										backupDirectory_results_stim   = parentDirectory+'/'+str(int(TstartRelax))+'_sec'
										## This is the output directory under which simulation output is saved.
										outputDirectory_rel_after_stim = parentDirectory+'/relax_after_'+str(int(TstartRelax))+'_sec'

										# check whether backup directory from which simulation is supposed to start exists
										if os.path.exists( backupDirectory_results_stim ):
											print( 'backup found' )

											# Check whether simulation has already been performed by looking
											# for output (Adjust time as needed. Note that backups are only generated at certain times.)
											testFileExistences = outputDirectory_rel_after_stim + "/relaxation_Tsim_"+str(Trelax)+"/10220_sec/cMatrix.npz"

											# if it does not exist
											if os.path.isfile( testFileExistences ) == False:
												print( 'submit simulation' )
												# build submission string
												shellCommand = 'python '+str(pathSherlockHome_relax)+'less_relaxation.py '+ backupDirectory_results_stim + ' ' + str(int(Trelax)) + ' ' + outputDirectory_rel_after_stim + '\n'

												# print shell command
												print( shellCommand )

												# the following code was used for submission to computing cluster
												# submissionString = 'srun ' + shellCommand

												# # adds job to end of list
												# with open(listOfJobsThatNeedToBeDone, "a") as myfile:
																				
												# 	myfile.write( submissionString )
												# 	myfile.write( "##\n" )
											else:
												print( 'simulation has already been done' )  
												
										else:
											print( 'backup not found' )

## D) use the following to simulate relaxation after non-shuffled CRS 
#
# 			python 3_run.py relaxation_after_non_shuffled_CR
#
if sys.argv[1] == 'relaxation_after_non_shuffled_CR':


	NETWORK_TYPES = ["homogeneous_network","inhomogeneous_network"]

	for NETWORK_TYPE in NETWORK_TYPES:

		# Directory in which results of simulations of non-shuffled CRS are saved.
		pathSherlockScratch        = dataDirectory + '/'+str(NETWORK_TYPE)+'/non_shuffled_CRS'
		# This is the directory which "less_relaxation.py" is located.
		pathSherlockHome_relax     = "../relaxation_after_stimulation/"
		# the following output file was used for submission to computing cluster
		# listOfJobsThatNeedToBeDone = 'relaxation_non_shuffled_CRS_network.txt'

		# parameters
		# list of CR sequences
		sequenceArray = ['0_1_2_3','0_1_3_2','0_2_1_3','0_2_3_1','0_3_1_2','0_3_2_1']
		# stimulation frequency
		fCRarray = [10.0] # Hz
		# stimulation amplitude
		AstimArray = [1.0]
		# stimulation time
		# Tstim = 10200.0 # sec
		Tstim = 7220.0 # sec
		# intraburst frequence (not implemented) (only needed for multiple stimulus pulses per stimulus burst)
		fintra = 130.0 # Hz
		# number of stimulus pulses per stimulus burst (not implemented) 
		number_of_pulses = 1
		# distance between stimulation sites in units of system length scale L=5 mm
		dsites = 0.25
		# effective rectangular width of stimulus profile in units of dsites
		erw = 0.5 
		# number of separately stimulated subpopulations (only implemented for 4)
		nElectrode = 4

		# time after which relaxation is started (this is the time at which stimulation ceases)
		TstartRelax = 10200 # sec
		# relaxation time (as a default value this is set to 3 hours and 20 sec). We added 220 sec to make sure that a backup is generated.
		Trelax = 10820+220 # sec
		# Time of backup 3 hours after cessation of stimulation relaxation. 
		TafterRelax = TstartRelax+Trelax # sec


		# run through seeds
		for seed in seed_array:
			for fCR in fCRarray:
				for Astim in AstimArray:
					for sequence in sequenceArray:

						# directory in which backup from stimulation is stored
						parentDirectory = pathSherlockScratch + "/seed_"+str(seed)+"/fixed_Sequence_CR/"+sequence+"/fCR_"+str(fCR)+"_M_"+str(nElectrode)+"/Astim_"+str(Astim)+"_Tstim_"+str(Tstim)+"_dsites_"+str(dsites)+"_erw_"+str(erw)
						
						# path backup from stimulation from which simulation is started
						backupDirectory_results_stim   = parentDirectory+'/'+str(int(TstartRelax))+'_sec'
						# this is the directory in which output is saved
						outputDirectory_rel_after_stim = parentDirectory+'/relax_after_'+str(int(TstartRelax))+'_sec'

						# check whether backup directory exists
						if os.path.exists( backupDirectory_results_stim ):

							print( 'backup found' )

							# Check whether simulation has already been performed by looking
							# for output (Adjust time as needed. Note that backups are only generated at certain times.)
							testFileExistences = outputDirectory_rel_after_stim + "/relaxation_Tsim_"+str(Trelax)+"/10220_sec/cMatrix.npz"

							# if it does not exist
							if os.path.isfile( testFileExistences ) == False:

								print( 'submit simulation' )
								
								# build submission string
								shellCommand = 'python '+str(pathSherlockHome_relax)+'less_relaxation.py '+ backupDirectory_results_stim + ' ' + str(int(Trelax)) + ' ' + outputDirectory_rel_after_stim + '\n' 

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





