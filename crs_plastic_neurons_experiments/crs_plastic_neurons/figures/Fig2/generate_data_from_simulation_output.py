##################################
# written by: Justus Kromer
##################################

import numpy as np
import sys 
import os
import pickle
import scipy.sparse

# specify directory in which simulation data are saved
simDataDirectory = "../../run_sim/data/"


#########################################################
# 1) copy weight and adjacency matrix to output_directory
#########################################################
if sys.argv[1] == '1_copy_weight_and_adjacency_matrix':

	output_directory = "data/cMatrix_conMatrix"

	# set true to load results initial networks
	loadInit = True
	# set true to load results for non-shuffled CR
	loadNonShuffled = False
	# set true to load results for shuffled CR
	loadShuffled = False

	# Generate output_directory if it does not already exist.
	if os.path.isdir(output_directory)== False:
		os.mkdir(output_directory)
		
	# run over seeds
	for seed in [10,12,14,16,18]:
		
		for network_type in ["inhomogeneous","homogenous"]:
			print(seed,network_type)

			##################################################
			## 1) initial network configurations
			if loadInit == True:
				print("loading data for initial networks")
				
				if network_type == "inhomogeneous":
					# inhomogeneous networks
					directory_init = simDataDirectory + "/inhomogeneous_network/initial_network/seed_"+str(seed)+"/3000_sec/"
				elif network_type == "homogenous":
					# homogeneous networks
					directory_init = simDataDirectory + "/homogeneous_network/initial_network/seed_"+str(seed)+"/3000_sec/"

				# construct filename for cMatrix and adjacency matrix
				filename_cMatrix_init             = directory_init+"cMatrix.npz"
				filename_connectivity_matrix_init = directory_init+"synConnections.npz"

				# If both files exist, copy them to output_directory.
				if os.path.isfile( filename_cMatrix_init ) and os.path.isfile( filename_connectivity_matrix_init ):
					cMatrix_init = scipy.sparse.load_npz( filename_cMatrix_init )
					connectivity_matrix_init = scipy.sparse.load_npz( filename_connectivity_matrix_init )

					# save matrices
					scipy.sparse.save_npz( output_directory+"/cMatrix_init_seed_"+str(seed)+"_"+network_type+".npz", cMatrix_init )
					# save matrices
					scipy.sparse.save_npz( output_directory+"/connectivity_matrix_init_seed_"+str(seed)+"_"+network_type+".npz", connectivity_matrix_init )

			##################################################
			## 2) non-shuffled CR
			if loadNonShuffled == True:
				print("loading data for non-shuffled CR")
				for sequence in ['0_1_2_3','0_1_3_2','0_2_1_3','0_2_3_1','0_3_1_2','0_3_2_1']:         
					# non-shuffled CR
					print(seed,network_type,sequence)
					if network_type == "inhomogeneous":
						# inhomogeneous networks
						directory_nonShuffled_cMatrix_acute =        simDataDirectory + "/inhomogeneous_network/non_shuffled_CRS/seed_"+str(seed)+"/fixed_Sequence_CR/"+str(sequence)+"/fCR_10.0_M_4/Astim_1.0_Tstim_7220.0_dsites_0.25_erw_0.5/10200_sec/"
						directory_nonShuffled_cMatrix_long_lasting = simDataDirectory + "/inhomogeneous_network/non_shuffled_CRS/seed_"+str(seed)+"/fixed_Sequence_CR/"+str(sequence)+"/fCR_10.0_M_4/Astim_1.0_Tstim_7220.0_dsites_0.25_erw_0.5/relax_after_10200_sec/relaxation_Tsim_11040/21220_sec/"
					elif network_type == "homogenous":
						# homogeneous networks
						directory_nonShuffled_cMatrix_acute =          simDataDirectory + "/homogeneous_network/non_shuffled_CRS/seed_"+str(seed)+"/fixed_Sequence_CR/"+str(sequence)+"/fCR_10.0_M_4/Astim_1.0_Tstim_7220.0_dsites_0.25_erw_0.5/10200_sec/"
						directory_nonShuffled_cMatrix_long_lasting =   simDataDirectory + "/homogeneous_network/non_shuffled_CRS/seed_"+str(seed)+"/fixed_Sequence_CR/"+str(sequence)+"/fCR_10.0_M_4/Astim_1.0_Tstim_7220.0_dsites_0.25_erw_0.5/relax_after_10200_sec/relaxation_Tsim_11040/21220_sec/"

					filename_nonShuffled_cMatrix_acute =                    directory_nonShuffled_cMatrix_acute+ "cMatrix.npz"
					filename_nonShuffled_connectivity_matrix_acute =        directory_nonShuffled_cMatrix_acute+"synConnections.npz"

					filename_nonShuffled_cMatrix_long_lasting =            directory_nonShuffled_cMatrix_long_lasting+  "cMatrix.npz"
					filename_nonShuffled_connectivity_matrix_long_lasting = directory_nonShuffled_cMatrix_long_lasting+"synConnections.npz"

					# If both files for the acute case exist, copy them to output_directory.
					if os.path.isfile( filename_nonShuffled_cMatrix_acute ) and os.path.isfile( filename_nonShuffled_connectivity_matrix_acute ):
						nonShuffled_cMatrix_acute = scipy.sparse.load_npz( filename_nonShuffled_cMatrix_acute )
						nonShuffled_connectivity_matrix_acute = scipy.sparse.load_npz( filename_nonShuffled_connectivity_matrix_acute )

						# save matrices
						scipy.sparse.save_npz( output_directory+"/nonShuffled_cMatrix_acute_seed_"+str(seed)+"_"+sequence+"_"+network_type+".npz", nonShuffled_cMatrix_acute )
						# save matrices
						scipy.sparse.save_npz( output_directory+"/nonShuffled_connectivity_matrix_acute_seed_"+str(seed)+"_"+sequence+"_"+network_type+".npz", nonShuffled_connectivity_matrix_acute )
					else:
						print("ERROR:",filename_nonShuffled_cMatrix_acute,filename_nonShuffled_connectivity_matrix_acute,"not found")

					# If both files for the long-lasting case exist, copy them to output_directory.
					if os.path.isfile( filename_nonShuffled_cMatrix_long_lasting ) and os.path.isfile( filename_nonShuffled_connectivity_matrix_long_lasting ):
						nonShuffled_cMatrix_long_lasting = scipy.sparse.load_npz( filename_nonShuffled_cMatrix_long_lasting )
						nonShuffled_connectivity_matrix_long_lasting = scipy.sparse.load_npz( filename_nonShuffled_connectivity_matrix_long_lasting )

						# save matrices
						scipy.sparse.save_npz( output_directory+"/nonShuffled_cMatrix_longLasting_seed_"+str(seed)+"_"+sequence+"_"+network_type+".npz", nonShuffled_cMatrix_long_lasting )
						# save matrices
						scipy.sparse.save_npz( output_directory+"/nonShuffled_connectivity_matrix_longLasting_seed_"+str(seed)+"_"+sequence+"_"+network_type+".npz", nonShuffled_connectivity_matrix_long_lasting )
					else:
						print("ERROR:",filename_nonShuffled_cMatrix_long_lasting,filename_nonShuffled_connectivity_matrix_long_lasting,"not found")

			##################################################
			## 3) shuffled CR
			if loadShuffled == True:            
				print("loading data for shuffled CR")
				for Tshuffle in [0.1, 10.0, 1800.0]:
					for seedSeq in np.arange(0,300,10).astype(int):
						print(seed,network_type,Tshuffle,seedSeq)
						# non-shuffled CR
						if network_type == "inhomogeneous":
							# inhomogeneous networks
							directory_shuffled_cMatrix_acute =          simDataDirectory + "/inhomogeneous_network/CRS_Tshuffle/seed_"+str(seed)+"/SVS_burst_sequence_seed/fCR_10.0_M_4_fintra_130.0_npb_1_Tshuffle_sec_"+str(Tshuffle)+"_seedSequence_"+str(seedSeq)+"_Astim_1.0_Tstim_7220.0_seedSeq_"+str(seedSeq)+"_dsites_0.25_erw_0.5/10200_sec/"
							directory_shuffled_cMatrix_long_lasting =   simDataDirectory + "/inhomogeneous_network/CRS_Tshuffle/seed_"+str(seed)+"/SVS_burst_sequence_seed/fCR_10.0_M_4_fintra_130.0_npb_1_Tshuffle_sec_"+str(Tshuffle)+"_seedSequence_"+str(seedSeq)+"_Astim_1.0_Tstim_7220.0_seedSeq_"+str(seedSeq)+"_dsites_0.25_erw_0.5/relax_after_10200_sec/relaxation_Tsim_11040/21220_sec/"
						elif network_type == "homogenous":
							# homogeneous networks
							directory_shuffled_cMatrix_acute =          simDataDirectory +   "/homogeneous_network/CRS_Tshuffle/seed_"+str(seed)+"/SVS_burst_sequence_seed/fCR_10.0_M_4_fintra_130.0_npb_1_Tshuffle_sec_"+str(Tshuffle)+"_seedSequence_"+str(seedSeq)+"_Astim_1.0_Tstim_7220.0_seedSeq_"+str(seedSeq)+"_dsites_0.25_erw_0.5/10200_sec/"
							directory_shuffled_cMatrix_long_lasting =   simDataDirectory +   "/homogeneous_network/CRS_Tshuffle/seed_"+str(seed)+"/SVS_burst_sequence_seed/fCR_10.0_M_4_fintra_130.0_npb_1_Tshuffle_sec_"+str(Tshuffle)+"_seedSequence_"+str(seedSeq)+"_Astim_1.0_Tstim_7220.0_seedSeq_"+str(seedSeq)+"_dsites_0.25_erw_0.5/relax_after_10200_sec/relaxation_Tsim_11040/21220_sec/"
						
						filename_shuffled_cMatrix_acute =                    directory_shuffled_cMatrix_acute+ "cMatrix.npz"
						filename_shuffled_connectivity_matrix_acute =        directory_shuffled_cMatrix_acute+"synConnections.npz"

						filename_shuffled_cMatrix_long_lasting =             directory_shuffled_cMatrix_long_lasting+  "cMatrix.npz"
						filename_shuffled_connectivity_matrix_long_lasting = directory_shuffled_cMatrix_long_lasting+"synConnections.npz"

						# If both files for the acute case exist, copy them to output_directory.
						if os.path.isfile( filename_shuffled_cMatrix_acute ) and os.path.isfile( filename_shuffled_connectivity_matrix_acute ):
							shuffled_cMatrix_acute = scipy.sparse.load_npz( filename_shuffled_cMatrix_acute )
							shuffled_connectivity_matrix_acute = scipy.sparse.load_npz( filename_shuffled_connectivity_matrix_acute )

							# save matrices
							scipy.sparse.save_npz( output_directory+ "/shuffled_cMatrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+network_type+".npz", shuffled_cMatrix_acute )
							# save matrices
							scipy.sparse.save_npz( output_directory+"/shuffled_connectivity_matrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+network_type+".npz", shuffled_connectivity_matrix_acute )
						else:
							print("ERROR:",filename_shuffled_cMatrix_acute,filename_shuffled_connectivity_matrix_acute,"not found")

						# If both files for the long-lasting case exist, copy them to output_directory.
						if os.path.isfile( filename_shuffled_cMatrix_long_lasting ) and os.path.isfile( filename_shuffled_connectivity_matrix_long_lasting ):
							shuffled_cMatrix_long_lasting = scipy.sparse.load_npz( filename_shuffled_cMatrix_long_lasting )
							shuffled_connectivity_matrix_long_lasting = scipy.sparse.load_npz( filename_shuffled_connectivity_matrix_long_lasting )

							# save matrices
							scipy.sparse.save_npz( output_directory+ "/shuffled_cMatrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+network_type+".npz", shuffled_cMatrix_long_lasting )
							# save matrices
							scipy.sparse.save_npz( output_directory+"/shuffled_connectivity_matrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+network_type+".npz", shuffled_connectivity_matrix_long_lasting )
						else:
							print("ERROR:",filename_shuffled_cMatrix_long_lasting,filename_shuffled_connectivity_matrix_long_lasting,"not found")


def getMeanWeight_and_append_to_List( List, filename_cMatrix, filename_connections ):
	# number of STN neurons
	NSTN = 1000
	if os.path.isfile(filename_cMatrix) and os.path.isfile(filename_connections):
		List.append( np.sum( scipy.sparse.load_npz( filename_cMatrix )[:NSTN,:NSTN] )/np.sum( scipy.sparse.load_npz( filename_connections )[:NSTN,:NSTN] )   )
	return List

def dic_save( dic ):
    with open( dic['filename'] + '.pickle', 'wb') as f:
        pickle.dump(dic, f)

############################################################################
# 2) calculate average mean synaptic weight and save results as dictionary
############################################################################
if sys.argv[1] == '2_calculate_average_mean_weight':

	# path to directory in which matrices are stored (see "output_directory" in under "if sys.argv[1] == '1_copy_weight_and_adjacency_matrix':" )
	directory = "data/cMatrix_conMatrix"
	output_directory = "data/"

	# This are the file endings used for "cMatrix_..." and "connectivity_matrix...".
	# In an earlier version of the code inhomogeneous networks were referred to as "circular" that is why some related variables are still called "circular".
	fileEnding_inhomogeneous =  "inhomogeneous" # "circular"
	fileEnding_homogeneous = "homogenous"

	####################################################################################
	# 1a) initial inhomogeneous networks
	mean_Weights_init_circular = {}
	# get data
	mws = []
	for seed in [10,12,14,16,18]:
		# get file names
		cMatrix_filename= directory + "/cMatrix_init_seed_"+str(seed)+"_"+fileEnding_inhomogeneous+".npz"
		conMat_filename = directory + "/connectivity_matrix_init_seed_"+str(seed)+"_"+fileEnding_inhomogeneous+".npz"
		
		mws = getMeanWeight_and_append_to_List( mws, cMatrix_filename, conMat_filename )

	mean_Weights_init_circular["data"] = mws
	mean_Weights_init_circular["mean"] = np.mean( mws )
	mean_Weights_init_circular["std"] = np.std( mws )

	mean_Weights_init_circular['filename']  = output_directory + "/dic_mean_Weights_init_inhomogeneous"
	dic_save( mean_Weights_init_circular )

	####################################################################################
	# 1b) initial homogeneous networks
	mean_Weights_init_homogeneous = {}
	# get data
	mws = []
	for seed in [10,12,14,16,18]:
		# get file names
		cMatrix_filename= directory + "/cMatrix_init_seed_"+str(seed)+"_"+fileEnding_homogeneous+".npz"
		conMat_filename = directory + "/connectivity_matrix_init_seed_"+str(seed)+"_"+fileEnding_homogeneous+".npz"
		mws = getMeanWeight_and_append_to_List( mws, cMatrix_filename, conMat_filename )

	mean_Weights_init_homogeneous["data"] = mws
	mean_Weights_init_homogeneous["mean"] = np.mean( mws )
	mean_Weights_init_homogeneous["std"] = np.std( mws )

	mean_Weights_init_homogeneous['filename']  = output_directory + "/dic_mean_Weights_init_homogeneous"
	dic_save( mean_Weights_init_homogeneous )

	####################################################################################
	# 2a) non-shuffled CR circular networks
	mean_Weights_nonShuffled_circular = {}
	# get data
	mws_acute = []
	mws_long_lasting = []
	for seed in [10,12,14,16,18]:
		for sequence in ['0_1_2_3','0_1_3_2','0_2_1_3','0_2_3_1','0_3_1_2','0_3_2_1']: 
			# get file names
			# a) acute
			cMatrix_filename= directory + "/nonShuffled_cMatrix_acute_seed_"+str(seed)+"_"+sequence+"_"+fileEnding_inhomogeneous+".npz"
			conMat_filename = directory + "/nonShuffled_connectivity_matrix_acute_seed_"+str(seed)+"_"+sequence+"_"+fileEnding_inhomogeneous+".npz"
			mws_acute = getMeanWeight_and_append_to_List( mws_acute, cMatrix_filename, conMat_filename )
			# b) long-lasting
			cMatrix_filename= directory + "/nonShuffled_cMatrix_longLasting_seed_"+str(seed)+"_"+sequence+"_"+fileEnding_inhomogeneous+".npz"
			conMat_filename = directory + "/nonShuffled_connectivity_matrix_longLasting_seed_"+str(seed)+"_"+sequence+"_"+fileEnding_inhomogeneous+".npz"
			mws_long_lasting = getMeanWeight_and_append_to_List( mws_long_lasting, cMatrix_filename, conMat_filename )
			
	mean_Weights_nonShuffled_circular["data_acute"] = mws_acute
	mean_Weights_nonShuffled_circular["mean_acute"] = np.mean( mws_acute )
	mean_Weights_nonShuffled_circular["std_acute"] = np.std( mws_acute )

	mean_Weights_nonShuffled_circular["data_long_lasting"] = mws_long_lasting
	mean_Weights_nonShuffled_circular["mean_long_lasting"] = np.mean( mws_long_lasting )
	mean_Weights_nonShuffled_circular["std_long_lasting"] = np.std( mws_long_lasting )

	mean_Weights_nonShuffled_circular['filename']  = output_directory + "/dic_mean_Weights_nonShuffled_inhomogeneous"
	dic_save( mean_Weights_nonShuffled_circular )

	####################################################################################
	# 2b) non-shuffled CR homogeneous networks
	mean_Weights_nonShuffled_homogeneous = {}
	# get data
	mws_acute = []
	mws_long_lasting = []
	for seed in [10,12,14,16,18]:
		for sequence in ['0_1_2_3','0_1_3_2','0_2_1_3','0_2_3_1','0_3_1_2','0_3_2_1']: 
			# get file names
			# a) acute
			cMatrix_filename= directory + "/nonShuffled_cMatrix_acute_seed_"+str(seed)+"_"+sequence+"_"+fileEnding_homogeneous+".npz"
			conMat_filename = directory + "/nonShuffled_connectivity_matrix_acute_seed_"+str(seed)+"_"+sequence+"_"+fileEnding_homogeneous+".npz"
			mws_acute = getMeanWeight_and_append_to_List( mws_acute, cMatrix_filename, conMat_filename )
			# b) long-lasting
			cMatrix_filename= directory + "/nonShuffled_cMatrix_longLasting_seed_"+str(seed)+"_"+sequence+"_"+fileEnding_homogeneous+".npz"
			conMat_filename = directory + "/nonShuffled_connectivity_matrix_longLasting_seed_"+str(seed)+"_"+sequence+"_"+fileEnding_homogeneous+".npz"
			mws_long_lasting = getMeanWeight_and_append_to_List( mws_long_lasting, cMatrix_filename, conMat_filename )
		
	mean_Weights_nonShuffled_homogeneous["data_acute"] = mws_acute
	mean_Weights_nonShuffled_homogeneous["mean_acute"] = np.mean( mws_acute )
	mean_Weights_nonShuffled_homogeneous["std_acute"] = np.std( mws_acute )

	mean_Weights_nonShuffled_homogeneous["data_long_lasting"] = mws_long_lasting
	mean_Weights_nonShuffled_homogeneous["mean_long_lasting"] = np.mean( mws_long_lasting )
	mean_Weights_nonShuffled_homogeneous["std_long_lasting"] = np.std( mws_long_lasting )

	mean_Weights_nonShuffled_homogeneous['filename']  = output_directory + "/dic_mean_Weights_nonShuffled_homogeneous"
	dic_save( mean_Weights_nonShuffled_homogeneous )

	####################################################################################
	# 3a) shuffled CR inhomogeneous networks Tshuf
	Tshuffle = 0.1 # sec
	mean_Weights_Tshuffle_01_circular = {}
	# get data
	mws_acute = []
	mws_long_lasting = []
	for seed in [10,12,14,16,18]:
		for seedSeq in np.arange(0,300,10).astype(int):
			# get file names
			# a) acute
			cMatrix_filename= directory + "/shuffled_cMatrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			mws_acute = getMeanWeight_and_append_to_List( mws_acute, cMatrix_filename, conMat_filename )
			# b) long-lasting
			cMatrix_filename= directory + "/shuffled_cMatrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			mws_long_lasting = getMeanWeight_and_append_to_List( mws_long_lasting, cMatrix_filename, conMat_filename )

	mean_Weights_Tshuffle_01_circular["data_acute"] = mws_acute
	mean_Weights_Tshuffle_01_circular["mean_acute"] = np.mean( mws_acute )
	mean_Weights_Tshuffle_01_circular["std_acute"] = np.std( mws_acute )

	mean_Weights_Tshuffle_01_circular["data_long_lasting"] = mws_long_lasting
	mean_Weights_Tshuffle_01_circular["mean_long_lasting"] = np.mean( mws_long_lasting )
	mean_Weights_Tshuffle_01_circular["std_long_lasting"] = np.std( mws_long_lasting )

	mean_Weights_Tshuffle_01_circular['filename']  = output_directory + "/dic_mean_Weights_Tshuffle_01_inhomogeneous"
	dic_save( mean_Weights_Tshuffle_01_circular )

	####################################################################################
	# 3b) shuffled CR homogeneous networks Tshuffle = 0.1 sec
	Tshuffle = 0.1 # sec
	mean_Weights_Tshuffle_01_homogeneous = {}
	# get data
	mws_acute = []
	mws_long_lasting = []
	for seed in [10,12,14,16,18]:
		for seedSeq in np.arange(0,300,10).astype(int):
			# get file names
			# a) acute
			cMatrix_filename= directory + "/shuffled_cMatrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			mws_acute = getMeanWeight_and_append_to_List( mws_acute, cMatrix_filename, conMat_filename )
			# b) long-lasting
			cMatrix_filename= directory + "/shuffled_cMatrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			mws_long_lasting = getMeanWeight_and_append_to_List( mws_long_lasting, cMatrix_filename, conMat_filename )

	mean_Weights_Tshuffle_01_homogeneous["data_acute"] = mws_acute
	mean_Weights_Tshuffle_01_homogeneous["mean_acute"] = np.mean( mws_acute )
	mean_Weights_Tshuffle_01_homogeneous["std_acute"] = np.std( mws_acute )

	mean_Weights_Tshuffle_01_homogeneous["data_long_lasting"] = mws_long_lasting
	mean_Weights_Tshuffle_01_homogeneous["mean_long_lasting"] = np.mean( mws_long_lasting )
	mean_Weights_Tshuffle_01_homogeneous["std_long_lasting"] = np.std( mws_long_lasting )

	mean_Weights_Tshuffle_01_homogeneous['filename']  = output_directory + "/dic_mean_Weights_Tshuffle_01_homogeneous"
	dic_save( mean_Weights_Tshuffle_01_homogeneous )

	####################################################################################
	# 3c) shuffled CR inhomogeneous networks Tshuffle = 10.0 sec
	Tshuffle = 10.0 # sec
	mean_Weights_Tshuffle_10_circular = {}
	# get data
	mws_acute = []
	mws_long_lasting = []
	for seed in [10,12,14,16,18]:
		for seedSeq in np.arange(0,300,10).astype(int):
			# get file names
			# a) acute
			cMatrix_filename= directory + "/shuffled_cMatrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			mws_acute = getMeanWeight_and_append_to_List( mws_acute, cMatrix_filename, conMat_filename )
			# b) long-lasting
			cMatrix_filename= directory + "/shuffled_cMatrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			mws_long_lasting = getMeanWeight_and_append_to_List( mws_long_lasting, cMatrix_filename, conMat_filename )

	mean_Weights_Tshuffle_10_circular["data_acute"] = mws_acute
	mean_Weights_Tshuffle_10_circular["mean_acute"] = np.mean( mws_acute )
	mean_Weights_Tshuffle_10_circular["std_acute"] = np.std( mws_acute )

	mean_Weights_Tshuffle_10_circular["data_long_lasting"] = mws_long_lasting
	mean_Weights_Tshuffle_10_circular["mean_long_lasting"] = np.mean( mws_long_lasting )
	mean_Weights_Tshuffle_10_circular["std_long_lasting"] = np.std( mws_long_lasting )

	mean_Weights_Tshuffle_10_circular['filename']  = output_directory + "/dic_mean_Weights_Tshuffle_10_inhomogeneous"
	dic_save( mean_Weights_Tshuffle_10_circular )

	####################################################################################
	# 3d) shuffled CR homogeneous networks Tshuffle = 10.0 sec
	Tshuffle = 10.0 # sec
	mean_Weights_Tshuffle_10_homogeneous = {}
	# get data
	mws_acute = []
	mws_long_lasting = []
	for seed in [10,12,14,16,18]:
		for seedSeq in np.arange(0,300,10).astype(int):
			# get file names
			# a) acute
			cMatrix_filename= directory + "/shuffled_cMatrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			mws_acute = getMeanWeight_and_append_to_List( mws_acute, cMatrix_filename, conMat_filename )
			# b) long-lasting
			cMatrix_filename= directory + "/shuffled_cMatrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			mws_long_lasting = getMeanWeight_and_append_to_List( mws_long_lasting, cMatrix_filename, conMat_filename )

	mean_Weights_Tshuffle_10_homogeneous["data_acute"] = mws_acute
	mean_Weights_Tshuffle_10_homogeneous["mean_acute"] = np.mean( mws_acute )
	mean_Weights_Tshuffle_10_homogeneous["std_acute"] = np.std( mws_acute )

	mean_Weights_Tshuffle_10_homogeneous["data_long_lasting"] = mws_long_lasting
	mean_Weights_Tshuffle_10_homogeneous["mean_long_lasting"] = np.mean( mws_long_lasting )
	mean_Weights_Tshuffle_10_homogeneous["std_long_lasting"] = np.std( mws_long_lasting )

	mean_Weights_Tshuffle_10_homogeneous['filename']  = output_directory + "/dic_mean_Weights_Tshuffle_10_homogeneous"
	dic_save( mean_Weights_Tshuffle_10_homogeneous )

	####################################################################################
	# 3e) shuffled CR inhomogeneous networks Tshuffle = 1800.0 sec
	Tshuffle = 1800.0 # sec
	mean_Weights_Tshuffle_1800_circular = {}
	# get data
	mws_acute = []
	mws_long_lasting = []
	for seed in [10,12,14,16,18]:
		for seedSeq in np.arange(0,300,10).astype(int):
			# get file names
			# a) acute
			cMatrix_filename= directory + "/shuffled_cMatrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			mws_acute = getMeanWeight_and_append_to_List( mws_acute, cMatrix_filename, conMat_filename )
			# b) long-lasting
			cMatrix_filename= directory + "/shuffled_cMatrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_inhomogeneous+".npz"
			mws_long_lasting = getMeanWeight_and_append_to_List( mws_long_lasting, cMatrix_filename, conMat_filename )


	mean_Weights_Tshuffle_1800_circular["data_acute"] = mws_acute
	mean_Weights_Tshuffle_1800_circular["mean_acute"] = np.mean( mws_acute )
	mean_Weights_Tshuffle_1800_circular["std_acute"] = np.std( mws_acute )

	mean_Weights_Tshuffle_1800_circular["data_long_lasting"] = mws_long_lasting
	mean_Weights_Tshuffle_1800_circular["mean_long_lasting"] = np.mean( mws_long_lasting )
	mean_Weights_Tshuffle_1800_circular["std_long_lasting"] = np.std( mws_long_lasting )

	mean_Weights_Tshuffle_1800_circular['filename']  = output_directory + "/dic_mean_Weights_Tshuffle_1800_inhomogeneous"
	dic_save( mean_Weights_Tshuffle_1800_circular )

	####################################################################################
	# 3f) shuffled CR homogeneous networks Tshuffle = 1800.0 sec
	Tshuffle = 1800.0 # sec
	mean_Weights_Tshuffle_1800_homogeneous = {}
	# get data
	mws_acute = []
	mws_long_lasting = []
	for seed in [10,12,14,16,18]:
		for seedSeq in np.arange(0,300,10).astype(int):
			# get file names
			# a) acute
			cMatrix_filename= directory + "/shuffled_cMatrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_acute_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			mws_acute = getMeanWeight_and_append_to_List( mws_acute, cMatrix_filename, conMat_filename )
			# b) long-lasting
			cMatrix_filename= directory + "/shuffled_cMatrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			conMat_filename = directory + "/shuffled_connectivity_matrix_long_lasting_seed_"+str(seed)+"_Tshuffle_"+str(Tshuffle)+"_seedSeq_"+str(seedSeq)+"_"+fileEnding_homogeneous+".npz"
			mws_long_lasting = getMeanWeight_and_append_to_List( mws_long_lasting, cMatrix_filename, conMat_filename )

	mean_Weights_Tshuffle_1800_homogeneous["data_acute"] = mws_acute
	mean_Weights_Tshuffle_1800_homogeneous["mean_acute"] = np.mean( mws_acute )
	mean_Weights_Tshuffle_1800_homogeneous["std_acute"] = np.std( mws_acute )

	mean_Weights_Tshuffle_1800_homogeneous["data_long_lasting"] = mws_long_lasting
	mean_Weights_Tshuffle_1800_homogeneous["mean_long_lasting"] = np.mean( mws_long_lasting )
	mean_Weights_Tshuffle_1800_homogeneous["std_long_lasting"] = np.std( mws_long_lasting )

	mean_Weights_Tshuffle_1800_homogeneous['filename']  = output_directory + "/dic_mean_Weights_Tshuffle_1800_homogeneous"
	dic_save( mean_Weights_Tshuffle_1800_homogeneous )







