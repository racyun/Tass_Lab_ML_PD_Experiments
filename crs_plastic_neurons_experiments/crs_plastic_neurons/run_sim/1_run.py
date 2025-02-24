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
###  Simulations for synchronized state
######################################################################
# The folowing code can be used to generate shell commands to run 
# simulations for synchronized states.

# A) use the following for inhomogeneous networks
#
# 			python 1_run.py inhomogeneous_network
#
if sys.argv[1] == 'inhomogeneous_network':

	# path to simulation script
	pathSherlockHome           = '../synch_states/'
	# output directory
	pathSherlockScratch        = dataDirectory + '/inhomogeneous_network/initial_network'

	# the following output file was used for submission to computing cluster
	# listOfJobsThatNeedToBeDone = 'synchState_inhomogeneous_network.txt'

	for seed in seed_array:

		# set output directory
		outputDirectory = pathSherlockScratch + '/seed_'+str(seed)

		# generate shell command
		shellCommand = 'python '+str(pathSherlockHome)+'get_stationary_states_arange_x_inhomogeneous_network.py '+ outputDirectory + ' ' + str(seed) +'\n' 
		
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
# 			python 1_run.py homogeneous_network
#
if sys.argv[1] == 'homogeneous_network':

	# path to simulation script
	pathSherlockHome           = '../synch_states/'
	# output directory
	pathSherlockScratch        = dataDirectory + '/homogeneous_network/initial_network'

	# the following output file was used for submission to computing cluster
	# listOfJobsThatNeedToBeDone = 'synchState_homogeneous_network.txt'

	for seed in seed_array:

		# set output directory
		outputDirectory = pathSherlockScratch + '/seed_'+str(seed)

		# generate shell command
		shellCommand = 'python '+str(pathSherlockHome)+'get_stationary_states_arange_x_homogeneous_network.py '+ outputDirectory + ' ' + str(seed) +'\n'

		# print shell command
		print( shellCommand )

		# the following code was used for submission to computing cluster
		# submissionString = 'srun ' + shellCommand

		# # adds job to end of list
		# with open(listOfJobsThatNeedToBeDone, "a") as myfile:
										
		# 	myfile.write( submissionString )
		# 	myfile.write( '##\n' )

	exit()


# C) use the following for intermediate networks 
#
# 			python 1_run.py intermediate_network
#
if sys.argv[1] == 'intermediate_network':

	# path to simulation script
	pathSherlockHome           = '../synch_states/'
	# output directory
	pathSherlockScratch        = dataDirectory + '/intermediate_network/initial_network'

	# the following output file was used for submission to computing cluster
	# listOfJobsThatNeedToBeDone = 'synchState_intermediate_network.txt'

	# The parameter lambPar scales between inhomogeneous network (lambPar=1) 
	# and homogeneous network (lamPar=0). 
	lambPar_array = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

	for lambPar in lambPar_array:
		for seed in seed_array:

			# set output directory
			outputDirectory = pathSherlockScratch + '/lambPar_'+str(lambPar)+'_seed_'+str(seed)

			# generate shell command
			shellCommand = 'python '+str(pathSherlockHome)+'get_stationary_states_arange_x_intermediate_network.py '+ outputDirectory + ' ' + str(seed) + ' ' + str(lambPar) +'\n' 
		
			# print shell command
			print( shellCommand )

			# the following code was used for submission to computing cluster
			# submissionString = 'srun ' + shellCommand

			# # adds job to end of list
			# with open(listOfJobsThatNeedToBeDone, "a") as myfile:
											
			# 	myfile.write( submissionString )
			# 	myfile.write( '##\n' )


	exit()













