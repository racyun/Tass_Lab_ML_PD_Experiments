##################################
# written by: Justus Kromer
##################################
import sys
import numpy as np
import os

# specify directory in which data will be saved
dataDirectory = "data"

################################################################################################
###  Simulations to analyze coexistence between synchronized and desynchronized state
################################################################################################
# The folowing code can be used to generate shell commands to run 
# simulations for different iniital mean synaptic weights. Results are shown in Figure 1.

# A) use the following for inhomogeneous networks
#
# 			python 1_run.py multistability_inhomogeneous
#
if sys.argv[1] == 'multistability_inhomogeneous':

	# path to simulation script
	pathSherlockHome           = '../multistability/'
	# output directory
	pathSherlockScratch        = dataDirectory + '/inhomogeneous_network/multistability'

	# the following output file was used for submission to computing cluster
	# listOfJobsThatNeedToBeDone = 'multistability_inhomogeneous_network.txt'

	# seed values (Seed=12 was used in Figure 1)
	seed_array = [12]

	# list of mean synaptic weights at t=0 for which simulations are performed
	initialMeanWeights = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

	for seed in seed_array:
		for mw in initialMeanWeights:

			# set output directory
			outputDirectory = pathSherlockScratch + '/seed_'+str(seed)+'/mw_'+str(mw)
			# output directory
			shellCommand = 'python '+str(pathSherlockHome)+'get_multistability_arange_x_inhomogeneous_network.py '+ outputDirectory + ' ' + str(seed) + ' ' + str(mw) + ' ' + '\n'
			
			# print shell command
			print( shellCommand )

			# the following code was used for submission to computing cluster
			# submissionString = 'srun ' + shellCommand

			# # adds job to end of list
			# with open(listOfJobsThatNeedToBeDone, "a") as myfile:
											
			# 	myfile.write( submissionString )
			# 	myfile.write( '##\n' )

	exit()


# B) use the following for homogeneous networks
#
# 			python 2_run.py multistability_homogeneous
#
if sys.argv[1] == 'multistability_homogeneous':

	# path to simulation script
	pathSherlockHome           = '../multistability/'
	# output directory
	pathSherlockScratch        = dataDirectory + '/homogeneous_network/multistability'

	# the following output file was used for submission to computing cluster
	# listOfJobsThatNeedToBeDone = 'multistability_homogeneous_network.txt'

	# seed values
	seed_array = [12]

	# list of mean synaptic weights at t=0 for which simulations are performed
	initialMeanWeights = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

	for seed in seed_array:
		for mw in initialMeanWeights:

			# set output directory
			outputDirectory = pathSherlockScratch + '/seed_'+str(seed)+'/mw_'+str(mw)

			shellCommand = 'python '+str(pathSherlockHome)+'get_multistability_arange_x_homogeneous_network.py '+ outputDirectory + ' ' + str(seed) + ' ' + str(mw) + ' ' + '\n' 
			
			# print shell command
			print( shellCommand )

			# the following code was used for submission to computing cluster
			# submissionString = 'srun ' + shellCommand

			# # adds job to end of list
			# with open(listOfJobsThatNeedToBeDone, "a") as myfile:
											
			# 	myfile.write( submissionString )
			# 	myfile.write( '##\n' )

	exit()
  

