##################################
# written by: Justus Kromer
##################################
# written and tested for Python 2.7.13
##################################

# imports
import sys
import numpy as np
import itertools 

##################################
# the following class is used by CRS_non_shuffled
class fixed_Sequence_CR_overlapping:


	def __init__(self, fCR, totalStimulationTime, M, dt, sequence ):

		# time interval between end of burst and beginning of next burst in ms
		self.fCR=fCR    # Hz

		# total stimulation time in sec
		self.Tstim=totalStimulationTime*1000.0 # ms

		# number of stimulation sites
		self.M=M

		# integration time step used in simulations
		self.dt=dt # ms

		# inverse time between first spikes of subsequent bursts
		self.burstFrequency = 0        # Hz

		# array that contains pre-calculated stimulus
		# kth element is value of stimulation current delivered at k*dt after stimulus onset
		self.signalOnTrain=np.zeros(1)

		# number of time steps in 'signalOnTrain'
		self.lengthSignalOneElectrode=0

		# integration time step at which current stimulus is delivered (beginning of next stimulus)
		# is used to get end of current stimulus and beginning of next one
		self.CurrentStimulusOnset=np.arange(self.M) # time steps

		# directory name for output from simulations using this stimulation protocol
		self.signalName=''

		# start of next cycle
		self.startNextCycle = 0
 
		self.sequence = sequence
		self.Sequence = np.array( sequence.split('_') ).astype( int )

		# print 'sequence'
		# print self.Sequence

	##################################
	#   function: initialize_Fixed_Sequence_CR_overlapping_Chaos_Paper(self )
	def initialize_Fixed_Sequence_CR_overlapping(self ):
	#
	#       initializes signal calculated the full signal of one electrode for one/Nelectrodes signal periods

		# number of integration time steps in a CR cycle
		self.time_steps_per_cycle_period = int(1000./(self.fCR * self.dt))

		# time steps of stimulus onset
		timeStepsToNextStimulusOnset = int( float( self.time_steps_per_cycle_period )/float(self.M) )
		self.stimOnsets = timeStepsToNextStimulusOnset * np.arange( self.M )

		########## generate pulse shape ##########
		### single pulse characteristics
		# pos rectangular pulses of unit amplitude, duration 0.2 ms followed by a 
		# negative counterpart of length 3 ms and amplitude 1/15, interpulse interval 1/130 s
		# positive rectangular pulse
		tStartPosPuls=0.2 # ms
		tStepStartPosPuls= int(tStartPosPuls/self.dt)
		lengthsPosRect=0.4*1.0 # ms
		tSteplengthsPosRect= int(lengthsPosRect/self.dt)

		# normalized such that integral over time in ms yields one
		AmpPosPuls=1.0/(0.4*1.0)

		# negative pulse
		tStartNegPuls=lengthsPosRect+0.2 # ms
		tStepStartNegPuls= int(tStartNegPuls/self.dt)
		lengthsNegRect=0.8*1.0 # ms         # motivated by Tass et al. 2012 (monkey study)
		tSteplengthsNegRect= int(lengthsNegRect/self.dt)
		AmpNegPuls= -(AmpPosPuls*lengthsPosRect)/lengthsNegRect # ensures charge balance by scaling the amplitude


		### minimal interval between subsequent pulses is adjusted to 130 Hz DBS pulsFrequency
		pulsFrequency=130.0 #  Hz
		pulsPeriod=1.0/float(pulsFrequency*0.001) # ms  // approx 7.69 ms

		self.pulsLength = int( (pulsPeriod)/self.dt ) # (number of timesteps until next pulse starts 

		#### in case of burst stimuli one stimulus consistes of 
		# create signal for one electrode
		self.lengthSignalOneElectrode=self.pulsLength
		self.signalOnTrain=np.zeros(self.lengthSignalOneElectrode)


		# directory in which output for this stimulation protocol is saved
		self.signalName='fixed_Sequence_CR/'+self.sequence+'/fCR_'+str(self.fCR)+'_M_'+str(self.M)

		# construct a single stimulus
		for ktimeSteps in range( self.pulsLength ):

			kStep=ktimeSteps

			# add pos pulse
			if ((ktimeSteps)<tSteplengthsPosRect+tStepStartPosPuls) and (tStepStartPosPuls) <= (ktimeSteps):

				self.signalOnTrain[kStep]=AmpPosPuls

			# add neg pulse
			if ((ktimeSteps)<tSteplengthsNegRect+tStepStartNegPuls) and (tStepStartNegPuls) <= (ktimeSteps):

				self.signalOnTrain[kStep]=AmpNegPuls


		# original currents to contact
		# rows contain currents for all time steps during one CR cycle 
		# m th row contains the currents for the m th stimulus activation in the CR cycle
		self.originalCurrentsToContacts = np.zeros( ( self.M , self.time_steps_per_cycle_period ) )
		# run through delivered stimuli
		if self.sequence != "0_0_0_0":
			for k in range( len( self.Sequence  ) ):
				# get stimulation site for stimulus k
				m = self.Sequence[k]

				# run through time step during one stimulus
				for kSignal in range( len(self.signalOnTrain) ):

					if kSignal+self.stimOnsets[k] < len( self.originalCurrentsToContacts[ m ] ):
						self.originalCurrentsToContacts[ m, kSignal+self.stimOnsets[k] ] = self.signalOnTrain[kSignal]
					else:
						print( 'Warning: '+str(m)+'th stimulus would reaches into next cycle and is cut off' )

		else:
			## implement periodic stimulation
			## activate all stimulation sites at 0
			# run through time step during one stimulus
			for kSignal in range( len(self.signalOnTrain) ):
				self.originalCurrentsToContacts[ :, kSignal ] = self.signalOnTrain[kSignal]*np.ones( len( self.Sequence  ) )
	


		# import matplotlib.pyplot as plt 
		# fig = plt.figure()
		# ax = fig.add_subplot(111)
		# ax.imshow (self.originalCurrentsToContacts )

		# ax.set_aspect( 10*len(self.originalCurrentsToContacts)/float(M) )
		# plt.show()

		return 0

	##################################
	#   function: enterNewCycle(self)
	def enterNewCycle(self, timeStep):

		# currents to contacts
		self.current_Currents = np.copy( self.originalCurrentsToContacts )

		#print self.current_Currents
		self.startNextCycle += self.time_steps_per_cycle_period

	##################################
	#   function: getCurrent( self, timeStep )
	def getCurrent( self, timeStep ):
		#print timeStep, self.startNextCycle
		if timeStep == self.startNextCycle:
			currentOutputTimeStep = 0
			self.enterNewCycle( timeStep )
		if timeStep > 0:
			currentOutputTimeStep = timeStep % self.time_steps_per_cycle_period
			#print currentOutputTimeStep
		else:
			currentOutputTimeStep = 0
		if timeStep < 0:
			#print 'exit 1'
			return 0*self.originalCurrentsToContacts[:,0]
		# return current
		#print 'exit 2'
		#print timeStep
		return self.current_Currents[ :,currentOutputTimeStep ]




## the following class is used by CRS_Tshuffle
class SVS_CR_bursts_less_folders:

	## input:           
	#   timeBetweenBursts      ... time between end of burst and beginning of next burst
	#   totalStimulationTime   ... total stimulation time in sec
	#   Nelectrodes            ... number of stimulation sites  (4 for standard CR stimultion)
	#  	fintra 				   ... intraburst frequency
	#   npb                    ... number of pulses per burst
	#   T_shuffle              ... shuffle time after which a new CR sequence is drawn from all possible sequences in sec
	#   seed_sequence          ... seed for creating the random sequence
	def __init__(self, fCR, totalStimulationTime, M, dt, fintra, npb, Tshuffle, seed_sequence ):
		# print("SVS_CR_bursts")
		# time interval between end of burst and beginning of next burst in ms
		self.fCR=fCR    # Hz

		# total stimulation time in sec
		self.Tstim=totalStimulationTime*1000.0 # ms

		# number of stimulation sites
		self.M=M

		# integration time step used in simulations
		self.dt=dt # ms

		# inverse time between first spikes of subsequent bursts
		self.burstFrequency = 0        # Hz

		# array that contains pre-calculated stimulus
		# kth element is value of stimulation current delivered at k*dt after stimulus onset
		self.signalOnTrain=np.zeros(1)

		# number of time steps in 'signalOnTrain'
		self.lengthSignalOneElectrode=0

		# integration time step at which current stimulus is delivered (beginning of next stimulus)
		# is used to get end of current stimulus and beginning of next one
		self.CurrentStimulusOnset=np.arange(self.M) # time steps

		# directory name for output from simulations using this stimulation protocol
		self.signalName=''

		# start of next cycle
		self.startNextCycle = 0
 
		# intraburst frequency in Hz
		self.fintra = fintra 

		# number of pulses per burst
		self.npb = npb

		# shuffling time after which a new sequence is selected
		self.Tshuffle = Tshuffle # sec
		# get shuffling time in time steps
		self.Tshuffle_steps = int( self.Tshuffle * 1000.0/self.dt)

		## start with a random sequence
		self.nextSeedSequence = seed_sequence # integer that is used as seed to create the next CR sequence when Tshuffle has passed
		
		# directory in which output for this stimulation protocol is saved
		self.signalName='SVS_burst_sequence_seed/fCR_'+str(self.fCR)+'_M_'+str(self.M)+'_fintra_'+str(self.fintra)+'_npb_'+str(self.npb)+'_Tshuffle_sec_'+str(self.Tshuffle)+'_seedSequence_'+str( seed_sequence )
		# print( '#################' )
		# print( self.signalName )

	##################################
	#   function: initialize_Fixed_Sequence_CR_overlapping_Chaos_Paper(self )
	def initialize_SVS_CR_bursts(self ):
	#
	#       initializes signal calculated the full signal of one electrode for one/Nelectrodes signal periods

		# number of integration time steps in a CR cycle
		self.time_steps_per_cycle_period = int(1000./(self.fCR * self.dt))

		# time steps of stimulus onset
		timeStepsToNextStimulusOnset = int( float( self.time_steps_per_cycle_period )/float(self.M) )
		self.stimOnsets = timeStepsToNextStimulusOnset * np.arange( self.M )

		##########################################
		########## generate pulse shape ##########
		##########################################
		### single pulse characteristics
		# pos rectangular pulses of unit amplitude, duration 0.2 ms followed by a 
		# negative counterpart of length 3 ms and amplitude 1/15, interpulse interval 1/130 s
		# positive rectangular pulse
		tStartPosPuls=0.2 # ms
		tStepStartPosPuls= int(tStartPosPuls/self.dt)
		lengthsPosRect=0.4*1.0 # ms
		tSteplengthsPosRect= int(lengthsPosRect/self.dt)

		# normalized such that integral over time in ms yields one
		AmpPosPuls=1.0/(0.4*1.0)

		# negative pulse
		tStartNegPuls=lengthsPosRect+0.2 # ms
		tStepStartNegPuls= int(tStartNegPuls/self.dt)
		lengthsNegRect=0.8*1.0 # ms         # motivated by Tass et al. 2012 (monkey study)
		tSteplengthsNegRect= int(lengthsNegRect/self.dt)
		AmpNegPuls= -(AmpPosPuls*lengthsPosRect)/lengthsNegRect # ensures charge balance by scaling the amplitude

		### minimal interval between subsequent pulses is adjusted to 130 Hz DBS pulsFrequency
		# number of pulses per burst
		
		pulsFrequency=self.fintra #  Hz
		pulsPeriod=1.0/float(pulsFrequency*0.001) # ms  // approx 7.69 ms

		self.pulsLength = int( (pulsPeriod)/self.dt ) # (number of timesteps until next pulse starts 

		#### in case of burst stimuli one stimulus consistes of 
		# create signal for one electrode
		self.lengthSignalOneElectrode=self.pulsLength
		self.signalOnTrain=np.zeros(self.lengthSignalOneElectrode)



		#############################################
		########## construct a single stimulus ######
		#############################################
		for ktimeSteps in range( self.pulsLength ):

			kStep = ktimeSteps

			# add pos pulse
			if ((ktimeSteps)<tSteplengthsPosRect+tStepStartPosPuls) and (tStepStartPosPuls) <= (ktimeSteps):

				self.signalOnTrain[kStep]=AmpPosPuls

			# add neg pulse
			if ((ktimeSteps)<tSteplengthsNegRect+tStepStartNegPuls) and (tStepStartNegPuls) <= (ktimeSteps):

				self.signalOnTrain[kStep]=AmpNegPuls

		# original currents to contact
		# rows contain currents for all time steps during one CR cycle 
		# m th row contains the currents for the m th stimulus activation in the CR cycle
		self.originalCurrentsToContacts = np.zeros( ( self.M , self.time_steps_per_cycle_period ) )
		
		##########################################
		########## get next random sequence ######
		##########################################
		# set seed of random number generator
		np.random.seed( self.nextSeedSequence )
		# create random sequence
		self.Sequence = np.random.choice( self.M, self.M, replace=False)
		# print( self.Sequence )
		# get next random seedSequence
		self.nextSeedSequence = np.random.randint(0,2147483647) 


		####################################################################################
		########## get stimulus intensity for individual time bins and subpopulations ######
		####################################################################################
		# check which pattern is delivered
		# self.Sequence = [0,0,0,0] corresponds to periodic stimulation
		# run through stimulation contacts
		for k in range( len( self.Sequence  ) ):

			m = self.Sequence[k]

			### add currents according to individual stimuli
			# run through pulses per burst
			for kpulse in range( self.npb ):

				koffset = kpulse * self.pulsLength
				# print( koffset ,self.time_steps_per_cycle_period )
				# run through indivudual pulses
				for kSignal in range( len(self.signalOnTrain) ):

					kcurrent = ( kSignal + koffset + self.stimOnsets[k] ) % self.time_steps_per_cycle_period

					self.originalCurrentsToContacts[ m, kcurrent ] = self.signalOnTrain[kSignal]

		return 0


	##################################
	#   function: enterNewCycle(self)
	def enterNewCycle(self, timeStep):

		# currents to contacts
		self.current_Currents = np.copy( self.originalCurrentsToContacts )

		#print self.current_Currents
		self.startNextCycle += self.time_steps_per_cycle_period

	##################################
	#   function: getCurrent( self, timeStep )
	#	timeStep ... in units of 0.1 ms
	def getCurrent( self, timeStep ):

		# shuffle sequence
		if ( timeStep % self.Tshuffle_steps == 0 ):
			if timeStep > 0:
				# print('shuffle', timeStep)
				# generate array of stimulus intensities for each time bin and contact for new sequence
				self.initialize_SVS_CR_bursts()
				self.current_Currents = np.copy( self.originalCurrentsToContacts )

		#print timeStep, self.startNextCycle
		if timeStep == self.startNextCycle:
			# print( self.Sequence )
			currentOutputTimeStep = 0
			self.enterNewCycle( timeStep )
		if timeStep > 0:
			currentOutputTimeStep = timeStep % self.time_steps_per_cycle_period
			#print currentOutputTimeStep
		else:
			currentOutputTimeStep = 0
		if timeStep < 0:
			#print 'exit 1'
			return 0*self.originalCurrentsToContacts[:,0]
		# return current
		#print 'exit 2'
		#print timeStep
		return self.current_Currents[ :,currentOutputTimeStep ]









