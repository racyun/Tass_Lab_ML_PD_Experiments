##################################
# written by: Justus Kromer
##################################
# written and tested for Python 2.7.13
##################################

# imports
from scipy.interpolate import interp1d
import numpy as np 
import itertools
import scipy

##################################
# function: conProbability
#
#	Returns the cartesian product of a and b (list of all possible combinations of elements of 'a' an 'b')
#
#       input: a, b
#			a ... list of elements
#			b ... list of elements
#		return:  
#			list of all possible combinations of elements of a and b
def cartesianProduct(a,b):
	return np.array([x for x in itertools.product(a, b)])

#################################
# places STN neurons along a 1-dim axis in the range of x_STN_min, x_STN_max 
# and GPe neurons between x_GPe_min, x_GPe_max
def placeNeurons_1D( NSTN , NGPe, x_STN_min, x_STN_max, x_GPe_min, x_GPe_max ):
	#
	#       Places NSTN neurons in cuboid associated with subthalamic nucleus (STN) and 
	#		NGPe in cuboid associated with globus pallidus externus (GPe).
	#	    axes are aligned with coordinate system.	
	#
	#       input: NSTN, NGPe
	#           NSTN ... number of STN neurons that need to be placed
	#			NGPe ... number of GPe neurons that need to be placed
	#			rx ... max distance from center in x-direction
	#			ry ... max distance from center in y-direction
	#			rz ... max distance from center in z-direction
	#       return: STNNeuronPositons, GPeNeuronPositons	
	# 			STNNeuronPositons ... numpy array of STN neuron centers in 3d (mm) 
	# 			GPeNeuronPositons ... numpy array of GPe neuron centers in 3d (mm) 

	# 1) place STN neurons
	# (i) start with uniformly distributed positions in 1d .. 
	STNNeuronPositons=np.random.uniform( x_STN_min, x_STN_max, NSTN )
	
	# 2) place GPe neurons
	# (i) start with uniformly distributed positions in 1d .. 
	GPeNeuronPositons=np.random.uniform( x_GPe_min, x_GPe_max, NGPe )

	# return lists of neuron positions in 3d that were placed in STN and GPe volume, respectively 
	return STNNeuronPositons, GPeNeuronPositons

################################################################
### functions for inhomogeneous network
################################################################
### function: circular_network_synReversalsCmatrix_1D
# parameters:
# 		positionsSTNNeurons ... connection probability for purely inhomogeneous networks
#   	positionsGPeNeurons   ... connection probability for purely homogeneous networks
#		NSTN ... total number of STN neurons
#		NGPe ... total number of GPe neurons
#  		Pcon ... connection probability
#		M    ... number of subpopulations (only implemented for M=4)
def circular_network_synReversalsCmatrix_1D(positionsSTNNeurons, positionsGPeNeurons, NSTN, NGPe, Pcon, M ):

	# initiallize return matrix containing synReversals[i,j]
	# 1  ... presynaptic neuron j is connected to postsynapti neuron i by excitatory synapse
	# -1 ... presynaptic neuron j is connected to postsynapti neuron i  by inhibitory synapse
	# 0 ... no connections from j to i
	synReversals=np.zeros( (NSTN+NGPe, NSTN+NGPe) ) 

	# sort neurons according to x-coordinate
	positionsSTNNeurons = np.sort( positionsSTNNeurons )
	positionsGPeNeurons = np.sort( positionsGPeNeurons )

	# randomly fill matrix
	PopulationSize = int( float(NSTN)/float(M) )
	synReversals=np.zeros( (NSTN+NGPe, NSTN+NGPe) ) 
	synReversals[:NSTN,:NSTN] = np.random.choice( [0,1], (NSTN,NSTN), p=[1-Pcon,Pcon] )
	synReversals[PopulationSize:3*PopulationSize,:PopulationSize] = 0
	synReversals[2*PopulationSize:4*PopulationSize,PopulationSize:2*PopulationSize] = 0
	synReversals[3*PopulationSize:4*PopulationSize,2*PopulationSize:3*PopulationSize] = 0
	synReversals[:PopulationSize,3*PopulationSize:4*PopulationSize] = 0
	# # the total number of connections
	# totNumberOfConnection=int( np.round( Pcon * NSTN * NSTN ) )

	# set diagonal to zero
	np.fill_diagonal(synReversals, 0)

	# return result
	return synReversals

################################################################
#  function: sequence_paper_generate_circular_network 
#
#		input:  system_parameters , rnd_state_for_network_generation 
#				system_parameters  .. parameter set used for simulations
#				rnd_state_for_network_generation .. state of random number generator to be used for network generation
#				Pcon .. 
#		output: 
#				synConnections , cMatrix , neuronLoc , sim_objects
#   			synConnections  ... connectivity matrix    entries 1, -1 , 0 for exc., inh. and no connections, respectively
#				cMatrix         ... scipy.sparse matrix synaptic weight matrix  entries [0,1]
#				neuronLoc       ... struct containing  locations of STN and GPe neurons
#				sim_objects		... struct containing objects related to network structure that are generated to 
#									speed up simulations
def sequence_paper_generate_circular_network( system_parameters , rnd_state_for_network_generation ):

	# load needed parameters from system_parameters
	# number of STN neurons
	N_STN = system_parameters['N_STN']
	# number of GPe neurons
	N_GPe = system_parameters['N_GPe']
	# total number of neurons
	N = N_STN+ N_GPe
	# probabilityh for STN -> STN connection
	P_STN_STN = system_parameters['P_STN_STN']
	# probabilityh for STN -> GPe connection
	P_STN_GPe = system_parameters['P_STN_GPe']
	# probabilityh for GPe -> GPe connection
	P_GPe_GPe = system_parameters['P_GPe_GPe']
	# probabilityh for GPe -> STN connection
	P_GPe_STN = system_parameters['P_GPe_STN']

	# synaptic transmission delay in time steps
	StepsTauSynDelaySTNSTN=int(system_parameters['tauSynDelaySTNSTN']/system_parameters['dt']) # time steps
	StepsTauSynDelayGPeGPe=int(system_parameters['tauSynDelayGPeGPe']/system_parameters['dt']) # time steps
	StepsTauSynDelayGPeSTN=int(system_parameters['tauSynDelayGPeSTN']/system_parameters['dt']) # time steps
	StepsTauSynDelaySTNGPe=int(system_parameters['tauSynDelaySTNGPe']/system_parameters['dt']) # time steps

	# max strengths exc coupling
	cMaxExc=system_parameters['cMaxExc']	
	# mean initial strengths exc coupling
	cExcInit=system_parameters['cExcInit']

	# max inh coupling
	cMaxInh=system_parameters['cMaxInh']
	# mean initial strengths inh coupling
	cInhInit=system_parameters['cInhInit']

	# set state of random number generator
	np.random.set_state( rnd_state_for_network_generation  )


	x_STN_min =system_parameters['x_STN_min']
	x_STN_max =system_parameters['x_STN_max']
	
	x_GPe_min =system_parameters['x_GPe_min']
	x_GPe_max =system_parameters['x_GPe_max']

	# get connectivity matrix
	STNCenter, GPeCenter= placeNeurons_1D( N_STN , N_GPe, x_STN_min, x_STN_max, x_GPe_min, x_GPe_max )

	# sort neurons according to x-coordinate
	STNCenter = np.sort( STNCenter )
	GPeCenter = np.sort( GPeCenter )

	# synConnections= variable_distance_synReversalsCmatrix_1D(STNCenter, GPeCenter, N_STN, N_GPe, P_STN_STN, P_STN_GPe, P_GPe_GPe, P_GPe_STN, d_synaptic_length_scale )
	synConnections= circular_network_synReversalsCmatrix_1D(STNCenter, GPeCenter, N_STN, N_GPe, 2*P_STN_STN, 4 )


	# set diagonal to zero 
	diaZero=np.ones( (N,N) )-np.diag( np.ones( N ) )
	synConnections=synConnections*diaZero

	# decouple GPe  ( this is done since we only used STN neurons in our simulations, uncomment if not needed)
	for kNeuron in range(N_STN,N):
		synConnections[:,kNeuron]=np.zeros(N)
		synConnections[kNeuron,:]=np.zeros(N)

	#########################################################################################
	#    in the following additional arrays are introduced to speed up simulations
	# get indicec post and presynaptic neurons to speed up STDP
	PostSynNeurons = {}
	PreSynNeurons = {}

	# max numbers of corresponding synapses
	maxNumberOfPostSynapticNeurons=0
	maxNumberOfPreSynapticNeurons=0

	for kNeuron in range(N_STN):

		PostSynNeurons[kNeuron]=np.nonzero( ( synConnections[:,kNeuron].astype(int) ).tolist() )[0].tolist()
		PreSynNeurons[kNeuron]=np.nonzero( ( synConnections[kNeuron,:].astype(int) ).tolist() )[0].tolist()

		# add random intra-network connections in case of no post/pre neurons
		# this is to help getting fully connected networks
		if len(PostSynNeurons[kNeuron])==0:

			# add random connection
			if kNeuron < N_STN:
				kPost=np.random.choice( range(kNeuron)+range(kNeuron+1,N_STN) )
				synConnections[kPost,kNeuron]=1


		if len(PreSynNeurons[kNeuron])==0:
			# add random connection
			# this is to help getting fully connected networks
			if kNeuron < N_STN:
				kPre=np.random.choice( range(kNeuron)+range(kNeuron+1,N_STN) )
				synConnections[kNeuron,kPre]=1

		# update max numbers of connections
		if maxNumberOfPostSynapticNeurons<len(PostSynNeurons[kNeuron]):
			maxNumberOfPostSynapticNeurons=len(PostSynNeurons[kNeuron])
		if maxNumberOfPreSynapticNeurons<len(PreSynNeurons[kNeuron]):
			maxNumberOfPreSynapticNeurons=len(PreSynNeurons[kNeuron])

	# generate numpy array with post synaptic neurons to speed up simulations ...
	numpyPostSynapticNeurons=np.full((N,maxNumberOfPostSynapticNeurons),N+1)
	numpyPreSynapticNeurons=np.full((N,maxNumberOfPreSynapticNeurons),N+1)

	# ... and corresponding matrix containing transimission delays in time steps
	transmissionDelaysPostSynNeurons=np.full((N,maxNumberOfPostSynapticNeurons),-1.0)
	transmissionDelaysPreSynNeurons=np.full((N,maxNumberOfPreSynapticNeurons),-1.0)

	# gen numpy array with post synaptic neurons
	for kPreSyn in range(N_STN):
		postSynNeuronskPre=PostSynNeurons[kPreSyn]
		for kPostSyn in range(len(postSynNeuronskPre)):
			numpyPostSynapticNeurons[kPreSyn, kPostSyn]=postSynNeuronskPre[kPostSyn]

			kPostNeuron=postSynNeuronskPre[kPostSyn]
			if (kPreSyn < N_STN):
				if (kPostNeuron < N_STN):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelaySTNSTN

				if (kPostNeuron >= N_STN) and (kPostNeuron < N):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelaySTNGPe


			if (kPreSyn >= N_STN) and (kPreSyn < N):
				if (kPostNeuron < N_STN):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelayGPeSTN

				if (kPostNeuron >= N_STN) and (kPostNeuron < N):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelayGPeGPe

	# gen numpy array with post synaptic neurons
	for kPostSyn in range(N_STN):
		preSynNeuronskPre=PreSynNeurons[kPostSyn]
		for kPreSyn in range(len(preSynNeuronskPre)):
			numpyPreSynapticNeurons[kPostSyn, kPreSyn]=preSynNeuronskPre[kPreSyn]

			kPreNeuron=preSynNeuronskPre[kPreSyn]
			if (kPostSyn < N_STN):
				if (kPreNeuron < N_STN):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelaySTNSTN

				if (kPreNeuron >= N_STN) and (kPreNeuron < N):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelayGPeSTN

			if (kPostSyn >= N_STN) and (kPostSyn < N):
				if (kPreNeuron < N_STN):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelaySTNGPe

				if (kPreNeuron >= N_STN) and (kPreNeuron < N):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelayGPeGPe

	# synaptic weight matrix
	cMatrix=np.zeros( (N , N) )

	# initialize synaptic weights by setting random weights to zero so that the mean initial
	# synaptic weights are cExcInit/cMaxExc and cInhInit/cMaxInh for excitatory and inhibitory connections, 
	# respectively
	# mean inital weights
	if cMaxExc != 0:
		meanInitalExcWeight=cExcInit/cMaxExc
	else:
		meanInitalExcWeight=0

	if cMaxInh != 0:
		meanInitalInhWeight=cInhInit/cMaxInh
	else:
		meanInitalInhWeight=0

	# initialize excitatory connections
	P1e=meanInitalExcWeight
	P0e=1-P1e
	cMatrix[:,:N_STN]=np.random.choice([0.0,1.0],(N,N_STN),p=[P0e,P1e])

	# initialize inhibitory connections
	P1i=meanInitalInhWeight
	P0i=1-P1i
	cMatrix[:,N_STN:]=np.random.choice([0.0,1.0],(N,N_GPe),p=[P0i,P1i])

	# filter weights with actual connections according to connectivity matrix
	cMatrix=cMatrix*synConnections
	cMatrix=scipy.sparse.csc_matrix(cMatrix)
	csc_Zero=scipy.sparse.csc_matrix(np.zeros( ( N,N ) ))
	csc_Ones=scipy.sparse.csc_matrix(np.ones( ( N,N ) ))

	# output struct containing neuron positions in mm
	neuronLoc = { 'STN_center_mm' : STNCenter , 'GPe_center_mm' : GPeCenter }

	# output struct containing objects that are related to network structure but only needed during simulation
	sim_objects = { 'max_N_pre' : maxNumberOfPreSynapticNeurons ,'max_N_post' : maxNumberOfPostSynapticNeurons , 'numpyPostSynapticNeurons' : numpyPostSynapticNeurons , 'numpyPreSynapticNeurons' : numpyPreSynapticNeurons , 'td_PostSynNeurons' : transmissionDelaysPostSynNeurons , 'td_PreSynNeurons' : transmissionDelaysPreSynNeurons , 'csc_Zero' : csc_Zero , 'csc_Ones' : csc_Ones	}
	

	# return output
	return synConnections , cMatrix , neuronLoc , sim_objects

################################################################
### functions for homogeneous network
################################################################
# function: synReversalsCmatrix_homogeneous
#
#		input:   NSTN ... total number of STN neurons
#				 NGPe ... total number of GPe neurons
#				 P_STN_STN ... Probability for STN -> STN connections (total number of connections is P_STN_STN * ( NSTN * NSTN ) )
#				 P_STN_GPe ... Probability for STN -> GPe connections (total number of connections is P_STN_GPe * ( NSTN * NGPe ) )
#				 P_GPe_GPe ... Probability for GPe -> GPe connections (total number of connections is P_GPe_GPe * ( NGPe * NGPe ) )
#				 P_GPe_STN ... Probability for GPe -> STN connections (total number of connections is P_GPe_STN * ( NGPe * BSTN ) )
#
#		return:
#			synReversals ... block matrix of integers and dimension (NSTN+NGPe, NSTN+NGPe). synReversals[i,j] contains information about the 
#							 connection between presynatpic neuron j and postsynatpic neuron i
#							 synReversals[i,j] = 1 -> exc. connections from j to i
#							 synReversals[i,j] = -1 -> inh. connections from j to i
#							 synReversals[i,j] = 0 -> no connections from j to i
def synReversalsCmatrix_homogeneous(NSTN, NGPe, P_STN_STN, P_STN_GPe, P_GPe_GPe, P_GPe_STN ):

	# initiallize return matrix containing synReversals[i,j]
	# 1  ... presynaptic neuron j is connected to postsynapti neuron i by excitatory synapse
	# -1 ... presynaptic neuron j is connected to postsynapti neuron i  by inhibitory synapse
	# 0 ... no connections from j to i
	synReversals=np.zeros( (NSTN+NGPe, NSTN+NGPe) ) 

	##################################
	# CONNECTIONS FOR STN -> STN
	##################################

	# total number of STN -> STN connections (round is actually not necessary, just to avoid non integer input )
	totNumberOfConnection=int( np.round( P_STN_STN * NSTN * NSTN ) )

	# implement array of all possible STN -> STN connections
	# first index for pre, second for post synaptic neuron
	allPossibleSTNtoSTNconnections=cartesianProduct( np.arange(NSTN),np.arange(NSTN) )

	# same probabilty for all connections
	probs = np.ones( len(allPossibleSTNtoSTNconnections) )

	# exclude self connections
	indizesOfNonSelfconnections=allPossibleSTNtoSTNconnections[:,0]!=allPossibleSTNtoSTNconnections[:,1]
	probs=probs[indizesOfNonSelfconnections]	# 1/mm
	allPossibleSTNtoSTNconnections=allPossibleSTNtoSTNconnections[indizesOfNonSelfconnections]

	# normalize probabiltiy to one
	probs=1/np.sum(probs)*probs   # probs contains the probability for each connection to be selected if only a single connection was implemented
	
	# implement synaptic connections according to probabilities
	STNSTNconnectionsIndizes=np.random.choice( len(allPossibleSTNtoSTNconnections) , totNumberOfConnection, p=probs, replace=False)
	
	# STNSTNconnections is the array of all connections that were selected. First entry is index of presynaptic neuron, second entry index of the postsynaptic neuron
	STNSTNconnections=allPossibleSTNtoSTNconnections[ STNSTNconnectionsIndizes ]

	# add these connections to the return matrix 'synReversals'
	for connection in STNSTNconnections:
		# add excitatory connection
		# note that synReversals[i,j] is refers to connections from presynaptic neuron j to presynaptic neuron i
		synReversals[ connection[1],  connection[0] ]=1

	##################################	
	# CONNECTIONS FOR STN-> GPe
	##################################
	#  the connections probability for STN -> GPe connections does not depend on the distance

	# total number of STN-> GPe connections (round is actually not necessary, just to avoid non integer input )
	totNumberOfConnection=int(np.round(P_STN_GPe*NSTN*NGPe))

	# implement array of all possible STN -> GPe connections
	# first index for pre, second for post synaptic neuron      neuron indices for STN neurons are 0 - NSTN-1 and for GPe neurons NSTN - NSTN+NGPe-1
	allPossibleSTNtoGPeconnections=cartesianProduct( np.arange(NSTN),np.arange(NSTN, NSTN+NGPe ) )

	# same probabilty for all connections
	probs = np.ones( len(allPossibleSTNtoGPeconnections) )
	# normalize probabiltiy to one
	probs=1/np.sum(probs)*probs   # probs contains the probability for each connection to be selected if only a single connection was implemented
	
	# implement synaptic connections according to probabilities
	STNGPeConnectionsIndizes=np.random.choice(len(allPossibleSTNtoGPeconnections), totNumberOfConnection, p=probs, replace=False)
	
	# STNGPeconnections is the array of all connections that were selected. First entry is index of presynaptic neuron, second entry index of the postsynaptic neuron
	STNGPeconnections=allPossibleSTNtoGPeconnections[ STNGPeConnectionsIndizes ]

	# add these connections to the return matrix 'synReversals'
	for connection in STNGPeconnections:
		# add excitatory connection
		# note that synReversals[i,j] is refers to connections from presynaptic neuron j to presynaptic neuron i
		synReversals[ connection[1],  connection[0] ]=1



	##################################	
	# CONNECTIONS FOR GPe-> GPe
	##################################

	# total number of GPe-> GPe connections (round is actually not necessary, just to avoid non integer input )
	totNumberOfConnection=int(np.round(P_GPe_GPe*NGPe*NGPe))

	# implement array of all possible GPe -> GPe connections
	# first index for pre, second for post synaptic neuron      neuron indices for STN neurons are 0 - NSTN-1 and for GPe neurons NSTN - NSTN+NGPe-1
	allPossibleGPetoGPeconnections=cartesianProduct( np.arange( NSTN, NSTN+NGPe ) ,np. arange( NSTN, NSTN+NGPe )  )

	# same probabilty for all connections
	probs = np.ones( len(allPossibleGPetoGPeconnections) )

	# exclude self connections
	indizesOfNonSelfconnections=allPossibleGPetoGPeconnections[:,0]!=allPossibleGPetoGPeconnections[:,1]
	probs=probs[indizesOfNonSelfconnections]
	allPossibleGPetoGPeconnections=allPossibleGPetoGPeconnections[indizesOfNonSelfconnections]

	# normalize to probs one
	probs=1/np.sum(probs)*probs # probs contains the probability for each connection to be selected if only a single connection was implemented
	
	# implement synaptic connections according to probabilities
	GPeGPeconnectionsIndizes=np.random.choice(len(allPossibleGPetoGPeconnections), totNumberOfConnection, p=probs, replace=False)
	
	# GPeGPeconnections is the array of all connections that were selected. First entry is index of presynaptic neuron, second entry index of the postsynaptic neuron
	GPeGPeconnections=allPossibleGPetoGPeconnections[ GPeGPeconnectionsIndizes ]

	# add these connections to the return matrix 'synReversals'
	for connection in GPeGPeconnections:
		# add inhibitory connection
		# note that synReversals[i,j] is refers to connections from presynaptic neuron j to presynaptic neuron i
		synReversals[ connection[1],  connection[0] ]=-1


	##################################	
	# CONNECTIONS FOR GPe-> STN
	##################################
	#  the connections probability for STN -> GPe connections does not depend on the distance

	# total number of GPe-> STN connections (round is actually not necessary, just to avoid non integer input )
	totNumberOfConnection=int(np.round(P_GPe_STN*NGPe*NSTN))

	# implement array of all possible GPe -> STN connections
	# first index for pre, second for post synaptic neuron      neuron indices for STN neurons are 0 - NSTN-1 and for GPe neurons NSTN - NSTN+NGPe-1
	allPossibleGPetoSTNconnections=cartesianProduct( np.arange(NSTN, NSTN+NGPe ) ,np.arange(NSTN )  )

	# same probabilty for all connections
	probs = np.ones( len(allPossibleGPetoSTNconnections) )

	# normalize to probs one
	probs=1/np.sum(probs)*probs # probs contains the probability for each connection to be selected if only a single connection was implemented

	# GPeSTNconnections is the array of all connections that were selected. First entry is index of presynaptic neuron, second entry index of the postsynaptic neuron
	GPeSTNconnectionsIndizes=np.random.choice(len(allPossibleGPetoSTNconnections), totNumberOfConnection, p=probs, replace=False)
	
	# GPeSTNconnections is the array of all connections that were selected. First entry is index of presynaptic neuron, second entry index of the postsynaptic neuron
	GPeSTNconnections=allPossibleGPetoSTNconnections[ GPeSTNconnectionsIndizes ]

	# add these connections to the return matrix 'synReversals'
	for connection in GPeSTNconnections:
		# add inhibitory connection
		# note that synReversals[i,j] is refers to connections from presynaptic neuron j to presynaptic neuron i
		synReversals[ connection[1],  connection[0] ]=-1
	
	# return result
	return synReversals

################################################################
#  function: generate_connectivity_and_weight_matrix_homogeneous 
#
#		input:  system_parameters , rnd_state_for_network_generation 
#				system_parameters  .. parameter set used for simulations
#				rnd_state_for_network_generation .. state of random number generator to be used for network generation
#
#		output: 
#				synConnections , cMatrix , neuronLoc , sim_objects
#   			synConnections  ... connectivity matrix    entries 1, -1 , 0 for exc., inh. and no connections, respectively
#				cMatrix         ... scipy.sparse matrix synaptic weight matrix  entries [0,1]
#				neuronLoc       ... struct containing  locations of STN and GPe neurons
#				sim_objects		... struct containing objects related to network structure that are generated to 
#									speed up simulations
def generate_connectivity_and_weight_matrix_homogeneous( system_parameters , rnd_state_for_network_generation ):

	# load needed parameters from system_parameters
	# number of STN neurons
	N_STN = system_parameters['N_STN']
	# number of GPe neurons
	N_GPe = system_parameters['N_GPe']
	# total number of neurons
	N = N_STN+ N_GPe
	# probabilityh for STN -> STN connection
	P_STN_STN = system_parameters['P_STN_STN']
	# probabilityh for STN -> GPe connection
	P_STN_GPe = system_parameters['P_STN_GPe']
	# probabilityh for GPe -> GPe connection
	P_GPe_GPe = system_parameters['P_GPe_GPe']
	# probabilityh for GPe -> STN connection
	P_GPe_STN = system_parameters['P_GPe_STN']

	# synaptic transmission delay in time steps
	StepsTauSynDelaySTNSTN=int(system_parameters['tauSynDelaySTNSTN']/system_parameters['dt']) # time steps
	StepsTauSynDelayGPeGPe=int(system_parameters['tauSynDelayGPeGPe']/system_parameters['dt']) # time steps
	StepsTauSynDelayGPeSTN=int(system_parameters['tauSynDelayGPeSTN']/system_parameters['dt']) # time steps
	StepsTauSynDelaySTNGPe=int(system_parameters['tauSynDelaySTNGPe']/system_parameters['dt']) # time steps

	# max strengths exc coupling
	cMaxExc=system_parameters['cMaxExc']	
	# mean initial strengths exc coupling
	cExcInit=system_parameters['cExcInit']

	# max inh coupling
	cMaxInh=system_parameters['cMaxInh']
	# mean initial strengths inh coupling
	cInhInit=system_parameters['cInhInit']

	# set state of random number generator
	np.random.set_state( rnd_state_for_network_generation  )

	x_STN_min =system_parameters['x_STN_min']
	x_STN_max =system_parameters['x_STN_max']
	
	x_GPe_min =system_parameters['x_GPe_min']
	x_GPe_max =system_parameters['x_GPe_max']

	# get connectivity matrix
	STNCenter, GPeCenter= placeNeurons_1D( N_STN , N_GPe, x_STN_min, x_STN_max, x_GPe_min, x_GPe_max )

	# sort neurons according to x-coordinate
	STNCenter = np.sort( STNCenter )
	GPeCenter = np.sort( GPeCenter )

	synConnections= synReversalsCmatrix_homogeneous(N_STN, N_GPe, P_STN_STN, P_STN_GPe, P_GPe_GPe, P_GPe_STN )


	# set diagonal to zero 
	diaZero=np.ones( (N,N) )-np.diag( np.ones( N ) )
	synConnections=synConnections*diaZero

	# decouple GPe  ( this is done since we only used STN neurons in our simulations, uncomment if not needed)
	for kNeuron in range(N_STN,N):
		synConnections[:,kNeuron]=np.zeros(N)
		synConnections[kNeuron,:]=np.zeros(N)

	#########################################################################################
	#    in the following additional arrays are introduced to speed up simulations
	# get indicec post and presynaptic neurons to speed up STDP
	PostSynNeurons = {}
	PreSynNeurons = {}

	# max numbers of corresponding synapses
	maxNumberOfPostSynapticNeurons=0
	maxNumberOfPreSynapticNeurons=0

	for kNeuron in range(N_STN):

		PostSynNeurons[kNeuron]=np.nonzero( ( synConnections[:,kNeuron].astype(int) ).tolist() )[0].tolist()
		PreSynNeurons[kNeuron]=np.nonzero( ( synConnections[kNeuron,:].astype(int) ).tolist() )[0].tolist()

		# add random intra-network connections in case of no post/pre neurons
		# this is to help getting fully connected networks
		if len(PostSynNeurons[kNeuron])==0:

			# add random connection
			if kNeuron < N_STN:
				kPost=np.random.choice( range(kNeuron)+range(kNeuron+1,N_STN) )
				synConnections[kPost,kNeuron]=1


		if len(PreSynNeurons[kNeuron])==0:
			# add random connection
			# this is to help getting fully connected networks
			if kNeuron < N_STN:
				kPre=np.random.choice( range(kNeuron)+range(kNeuron+1,N_STN) )
				synConnections[kNeuron,kPre]=1

		# update max numbers of connections
		if maxNumberOfPostSynapticNeurons<len(PostSynNeurons[kNeuron]):
			maxNumberOfPostSynapticNeurons=len(PostSynNeurons[kNeuron])
		if maxNumberOfPreSynapticNeurons<len(PreSynNeurons[kNeuron]):
			maxNumberOfPreSynapticNeurons=len(PreSynNeurons[kNeuron])

	# generate numpy array with post synaptic neurons to speed up simulations ...
	numpyPostSynapticNeurons=np.full((N,maxNumberOfPostSynapticNeurons),N+1)
	numpyPreSynapticNeurons=np.full((N,maxNumberOfPreSynapticNeurons),N+1)

	# ... and corresponding matrix containing transimission delays in time steps
	transmissionDelaysPostSynNeurons=np.full((N,maxNumberOfPostSynapticNeurons),-1.0)
	transmissionDelaysPreSynNeurons=np.full((N,maxNumberOfPreSynapticNeurons),-1.0)

	# gen numpy array with post synaptic neurons
	for kPreSyn in range(N_STN):
		postSynNeuronskPre=PostSynNeurons[kPreSyn]
		for kPostSyn in range(len(postSynNeuronskPre)):
			numpyPostSynapticNeurons[kPreSyn, kPostSyn]=postSynNeuronskPre[kPostSyn]

			kPostNeuron=postSynNeuronskPre[kPostSyn]
			if (kPreSyn < N_STN):
				if (kPostNeuron < N_STN):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelaySTNSTN

				if (kPostNeuron >= N_STN) and (kPostNeuron < N):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelaySTNGPe


			if (kPreSyn >= N_STN) and (kPreSyn < N):
				if (kPostNeuron < N_STN):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelayGPeSTN

				if (kPostNeuron >= N_STN) and (kPostNeuron < N):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelayGPeGPe

	# gen numpy array with post synaptic neurons
	for kPostSyn in range(N_STN):
		preSynNeuronskPre=PreSynNeurons[kPostSyn]
		for kPreSyn in range(len(preSynNeuronskPre)):
			numpyPreSynapticNeurons[kPostSyn, kPreSyn]=preSynNeuronskPre[kPreSyn]

			kPreNeuron=preSynNeuronskPre[kPreSyn]
			if (kPostSyn < N_STN):
				if (kPreNeuron < N_STN):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelaySTNSTN

				if (kPreNeuron >= N_STN) and (kPreNeuron < N):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelayGPeSTN

			if (kPostSyn >= N_STN) and (kPostSyn < N):
				if (kPreNeuron < N_STN):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelaySTNGPe

				if (kPreNeuron >= N_STN) and (kPreNeuron < N):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelayGPeGPe

	# synaptic weight matrix
	cMatrix=np.zeros( (N , N) )

	# initialize synaptic weights by setting random weights to zero so that the mean initial
	# synaptic weights are cExcInit/cMaxExc and cInhInit/cMaxInh for excitatory and inhibitory connections, 
	# respectively
	# mean inital weights
	if cMaxExc != 0:
		meanInitalExcWeight=cExcInit/cMaxExc
	else:
		meanInitalExcWeight=0

	if cMaxInh != 0:
		meanInitalInhWeight=cInhInit/cMaxInh
	else:
		meanInitalInhWeight=0

	# initialize excitatory connections
	P1e=meanInitalExcWeight
	P0e=1-P1e
	cMatrix[:,:N_STN]=np.random.choice([0.0,1.0],(N,N_STN),p=[P0e,P1e])

	# initialize inhibitory connections
	P1i=meanInitalInhWeight
	P0i=1-P1i
	cMatrix[:,N_STN:]=np.random.choice([0.0,1.0],(N,N_GPe),p=[P0i,P1i])

	# filter weights with actual connections according to connectivity matrix
	cMatrix=cMatrix*synConnections
	cMatrix=scipy.sparse.csc_matrix(cMatrix)
	csc_Zero=scipy.sparse.csc_matrix(np.zeros( ( N,N ) ))
	csc_Ones=scipy.sparse.csc_matrix(np.ones( ( N,N ) ))

	# output struct containing neuron positions in mm
	neuronLoc = { 'STN_center_mm' : STNCenter , 'GPe_center_mm' : GPeCenter }

	# output struct containing objects that are related to network structure but only needed during simulation
	sim_objects = { 'max_N_pre' : maxNumberOfPreSynapticNeurons ,'max_N_post' : maxNumberOfPostSynapticNeurons , 'numpyPostSynapticNeurons' : numpyPostSynapticNeurons , 'numpyPreSynapticNeurons' : numpyPreSynapticNeurons , 'td_PostSynNeurons' : transmissionDelaysPostSynNeurons , 'td_PreSynNeurons' : transmissionDelaysPreSynNeurons , 'csc_Zero' : csc_Zero , 'csc_Ones' : csc_Ones	}
	

	# return output
	return synConnections , cMatrix , neuronLoc , sim_objects

################################################################
### functions for intermediate network
################################################################
### intermediate network case between inhomogeneous and homogeneous network
# parameters:
# 	PconInhomogeneous ... connection probability for purely inhomogeneous networks
#   PconHomogeneous   ... connection probability for purely homogeneous networks
#   lambPar           ... parameter that scales between homogeneous and inhomogeneous networks
# 							for lambPar = 0 the network is homogeneous and for lambPar = 1 the network is inhomogeneous
def transition_network_synReversalsCmatrix_1D(positionsSTNNeurons, positionsGPeNeurons, NSTN, NGPe, PconInhomogeneous, PconHomogeneous, lambPar, M ):

	# initiallize return matrix containing synReversals[i,j]
	# 1  ... presynaptic neuron j is connected to postsynapti neuron i by excitatory synapse
	# -1 ... presynaptic neuron j is connected to postsynapti neuron i  by inhibitory synapse
	# 0 ... no connections from j to i
	synReversals=np.zeros( (NSTN+NGPe, NSTN+NGPe) ) 

	# sort neurons according to x-coordinate
	positionsSTNNeurons = np.sort( positionsSTNNeurons )
	positionsGPeNeurons = np.sort( positionsGPeNeurons )

	# randomly fill matrix
	PopulationSize = int( float(NSTN)/float(M) )
	synReversals=np.zeros( (NSTN+NGPe, NSTN+NGPe) ) 

	# connection probabilities   (y axes is post, x axes is presyn population)
	# homogeneous network:         inhomogeneous network
	#          P P P P 					P 0 0 P
	#          P P P P                  0 0 P P
	#          P P P P                  0 P P P
	#          P P P P	                P P P 0
	#
	# combination of both:
	#		   PB PH PH PB 
	#          PH PH PB PB
	#          PH PB PB PB
	#          PB PB PB PH
	PB  = lambPar * PconInhomogeneous + (1-lambPar) * PconHomogeneous
	PH  = (1-lambPar) * PconHomogeneous

	# top row
	synReversals[3*PopulationSize:4*PopulationSize,:PopulationSize] = np.random.choice( [0,1], (PopulationSize,PopulationSize), p=[1-PB,PB] )
	synReversals[3*PopulationSize:4*PopulationSize,PopulationSize:3*PopulationSize]= np.random.choice( [0,1], (PopulationSize,2*PopulationSize), p=[1-PH,PH] )
	synReversals[3*PopulationSize:4*PopulationSize,3*PopulationSize:4*PopulationSize] = np.random.choice( [0,1], (PopulationSize,PopulationSize), p=[1-PB,PB] )

	# third row
	synReversals[2*PopulationSize:3*PopulationSize,:2*PopulationSize]= np.random.choice( [0,1], (PopulationSize,2*PopulationSize), p=[1-PH,PH] )
	synReversals[2*PopulationSize:3*PopulationSize,2*PopulationSize:4*PopulationSize]= np.random.choice( [0,1], (PopulationSize,2*PopulationSize), p=[1-PB,PB] )

	# second row
	synReversals[1*PopulationSize:2*PopulationSize,:PopulationSize]= np.random.choice( [0,1], (PopulationSize,PopulationSize), p=[1-PH,PH] )
	synReversals[1*PopulationSize:2*PopulationSize,PopulationSize:4*PopulationSize]= np.random.choice( [0,1], (PopulationSize,3*PopulationSize), p=[1-PB,PB] )

	# first row
	synReversals[:PopulationSize,:3*PopulationSize]= np.random.choice( [0,1], (PopulationSize,3*PopulationSize), p=[1-PB,PB] )
	synReversals[:PopulationSize,3*PopulationSize:4*PopulationSize]= np.random.choice( [0,1], (PopulationSize,PopulationSize), p=[1-PH,PH] )

	# # the total number of connections
	# totNumberOfConnection=int( np.round( Pcon * NSTN * NSTN ) )

	# set diagonal to zero
	np.fill_diagonal(synReversals, 0)

	# return result
	return synReversals

### this function generates an intermediate network
# The parameter "lambPar" scales between an inhomogeneous network as in Kromer, J. A. and Tass, P. A. (2024) (lambPar=1) 
#  and homogeneous network (lamPar=0)
def sequence_paper_generate_transition_network( system_parameters , rnd_state_for_network_generation, PconInhomogeneous, PconHomogeneous, lambPar ):

	# load needed parameters from system_parameters
	# number of STN neurons
	N_STN = system_parameters['N_STN']
	# number of GPe neurons
	N_GPe = system_parameters['N_GPe']
	# total number of neurons
	N = N_STN+ N_GPe

	########## THIS PART NEEDS TO BE CHANGED
	# probability for STN -> STN connection
	# P_STN_STN = system_parameters['P_STN_STN']

	# probability for STN -> GPe connection
	P_STN_GPe = system_parameters['P_STN_GPe']
	# probability for GPe -> GPe connection
	P_GPe_GPe = system_parameters['P_GPe_GPe']
	# probability for GPe -> STN connection
	P_GPe_STN = system_parameters['P_GPe_STN']

	# synaptic transmission delay in time steps
	StepsTauSynDelaySTNSTN=int(system_parameters['tauSynDelaySTNSTN']/system_parameters['dt']) # time steps
	StepsTauSynDelayGPeGPe=int(system_parameters['tauSynDelayGPeGPe']/system_parameters['dt']) # time steps
	StepsTauSynDelayGPeSTN=int(system_parameters['tauSynDelayGPeSTN']/system_parameters['dt']) # time steps
	StepsTauSynDelaySTNGPe=int(system_parameters['tauSynDelaySTNGPe']/system_parameters['dt']) # time steps

	# max strengths exc coupling
	cMaxExc=system_parameters['cMaxExc']	
	# mean initial strengths exc coupling
	cExcInit=system_parameters['cExcInit']

	# max inh coupling
	cMaxInh=system_parameters['cMaxInh']
	# mean initial strengths inh coupling
	cInhInit=system_parameters['cInhInit']

	# set state of random number generator
	np.random.set_state( rnd_state_for_network_generation  )


	x_STN_min =system_parameters['x_STN_min']
	x_STN_max =system_parameters['x_STN_max']
	
	x_GPe_min =system_parameters['x_GPe_min']
	x_GPe_max =system_parameters['x_GPe_max']

	# get connectivity matrix
	STNCenter, GPeCenter= placeNeurons_1D( N_STN , N_GPe, x_STN_min, x_STN_max, x_GPe_min, x_GPe_max )

	# sort neurons according to x-coordinate
	STNCenter = np.sort( STNCenter )
	GPeCenter = np.sort( GPeCenter )

	# synConnections= variable_distance_synReversalsCmatrix_1D(STNCenter, GPeCenter, N_STN, N_GPe, P_STN_STN, P_STN_GPe, P_GPe_GPe, P_GPe_STN, d_synaptic_length_scale )
	# synConnections= circular_network_synReversalsCmatrix_1D(STNCenter, GPeCenter, N_STN, N_GPe, 2*P_STN_STN, 4 )
	synConnections= transition_network_synReversalsCmatrix_1D(STNCenter, GPeCenter, N_STN , N_GPe, PconInhomogeneous, PconHomogeneous, lambPar, 4 )

	# set diagonal to zero 
	diaZero=np.ones( (N,N) )-np.diag( np.ones( N ) )
	synConnections=synConnections*diaZero

	# decouple GPe  ( this is done since we only used STN neurons in our simulations, uncomment if not needed)
	for kNeuron in range(N_STN,N):
		synConnections[:,kNeuron]=np.zeros(N)
		synConnections[kNeuron,:]=np.zeros(N)

	#########################################################################################
	#    in the following additional arrays are introduced to speed up simulations
	# get indicec post and presynaptic neurons to speed up STDP
	PostSynNeurons = {}
	PreSynNeurons = {}

	# max numbers of corresponding synapses
	maxNumberOfPostSynapticNeurons=0
	maxNumberOfPreSynapticNeurons=0

	for kNeuron in range(N_STN):

		PostSynNeurons[kNeuron]=np.nonzero( ( synConnections[:,kNeuron].astype(int) ).tolist() )[0].tolist()
		PreSynNeurons[kNeuron]=np.nonzero( ( synConnections[kNeuron,:].astype(int) ).tolist() )[0].tolist()

		# add random intra-network connections in case of no post/pre neurons
		# this is to help getting fully connected networks
		if len(PostSynNeurons[kNeuron])==0:

			# add random connection
			if kNeuron < N_STN:
				kPost=np.random.choice( range(kNeuron)+range(kNeuron+1,N_STN) )
				synConnections[kPost,kNeuron]=1


		if len(PreSynNeurons[kNeuron])==0:
			# add random connection
			# this is to help getting fully connected networks
			if kNeuron < N_STN:
				kPre=np.random.choice( range(kNeuron)+range(kNeuron+1,N_STN) )
				synConnections[kNeuron,kPre]=1

		# update max numbers of connections
		if maxNumberOfPostSynapticNeurons<len(PostSynNeurons[kNeuron]):
			maxNumberOfPostSynapticNeurons=len(PostSynNeurons[kNeuron])
		if maxNumberOfPreSynapticNeurons<len(PreSynNeurons[kNeuron]):
			maxNumberOfPreSynapticNeurons=len(PreSynNeurons[kNeuron])

	# generate numpy array with post synaptic neurons to speed up simulations ...
	numpyPostSynapticNeurons=np.full((N,maxNumberOfPostSynapticNeurons),N+1)
	numpyPreSynapticNeurons=np.full((N,maxNumberOfPreSynapticNeurons),N+1)

	# ... and corresponding matrix containing transimission delays in time steps
	transmissionDelaysPostSynNeurons=np.full((N,maxNumberOfPostSynapticNeurons),-1.0)
	transmissionDelaysPreSynNeurons=np.full((N,maxNumberOfPreSynapticNeurons),-1.0)

	# gen numpy array with post synaptic neurons
	for kPreSyn in range(N_STN):
		postSynNeuronskPre=PostSynNeurons[kPreSyn]
		for kPostSyn in range(len(postSynNeuronskPre)):
			numpyPostSynapticNeurons[kPreSyn, kPostSyn]=postSynNeuronskPre[kPostSyn]

			kPostNeuron=postSynNeuronskPre[kPostSyn]
			if (kPreSyn < N_STN):
				if (kPostNeuron < N_STN):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelaySTNSTN

				if (kPostNeuron >= N_STN) and (kPostNeuron < N):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelaySTNGPe


			if (kPreSyn >= N_STN) and (kPreSyn < N):
				if (kPostNeuron < N_STN):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelayGPeSTN

				if (kPostNeuron >= N_STN) and (kPostNeuron < N):
					transmissionDelaysPostSynNeurons[kPreSyn, kPostSyn]=StepsTauSynDelayGPeGPe

	# gen numpy array with post synaptic neurons
	for kPostSyn in range(N_STN):
		preSynNeuronskPre=PreSynNeurons[kPostSyn]
		for kPreSyn in range(len(preSynNeuronskPre)):
			numpyPreSynapticNeurons[kPostSyn, kPreSyn]=preSynNeuronskPre[kPreSyn]

			kPreNeuron=preSynNeuronskPre[kPreSyn]
			if (kPostSyn < N_STN):
				if (kPreNeuron < N_STN):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelaySTNSTN

				if (kPreNeuron >= N_STN) and (kPreNeuron < N):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelayGPeSTN

			if (kPostSyn >= N_STN) and (kPostSyn < N):
				if (kPreNeuron < N_STN):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelaySTNGPe

				if (kPreNeuron >= N_STN) and (kPreNeuron < N):
					transmissionDelaysPreSynNeurons[kPostSyn, kPreSyn]=StepsTauSynDelayGPeGPe

	# synaptic weight matrix
	cMatrix=np.zeros( (N , N) )

	# initialize synaptic weights by setting random weights to zero so that the mean initial
	# synaptic weights are cExcInit/cMaxExc and cInhInit/cMaxInh for excitatory and inhibitory connections, 
	# respectively
	# mean inital weights
	if cMaxExc != 0:
		meanInitalExcWeight=cExcInit/cMaxExc
	else:
		meanInitalExcWeight=0

	if cMaxInh != 0:
		meanInitalInhWeight=cInhInit/cMaxInh
	else:
		meanInitalInhWeight=0

	# initialize excitatory connections
	P1e=meanInitalExcWeight
	P0e=1-P1e
	cMatrix[:,:N_STN]=np.random.choice([0.0,1.0],(N,N_STN),p=[P0e,P1e])

	# initialize inhibitory connections
	P1i=meanInitalInhWeight
	P0i=1-P1i
	cMatrix[:,N_STN:]=np.random.choice([0.0,1.0],(N,N_GPe),p=[P0i,P1i])

	# filter weits with actual connections according to connectivity matrix
	cMatrix=cMatrix*synConnections
	cMatrix=scipy.sparse.csc_matrix(cMatrix)
	csc_Zero=scipy.sparse.csc_matrix(np.zeros( ( N,N ) ))
	csc_Ones=scipy.sparse.csc_matrix(np.ones( ( N,N ) ))

	# output struct containing neuron positions in mm
	neuronLoc = { 'STN_center_mm' : STNCenter , 'GPe_center_mm' : GPeCenter }

	# output struct containing objects that are related to network structure but only needed during simulation
	sim_objects = { 'max_N_pre' : maxNumberOfPreSynapticNeurons ,'max_N_post' : maxNumberOfPostSynapticNeurons , 'numpyPostSynapticNeurons' : numpyPostSynapticNeurons , 'numpyPreSynapticNeurons' : numpyPreSynapticNeurons , 'td_PostSynNeurons' : transmissionDelaysPostSynNeurons , 'td_PreSynNeurons' : transmissionDelaysPreSynNeurons , 'csc_Zero' : csc_Zero , 'csc_Ones' : csc_Ones	}
	

	# return output
	return synConnections , cMatrix , neuronLoc , sim_objects



