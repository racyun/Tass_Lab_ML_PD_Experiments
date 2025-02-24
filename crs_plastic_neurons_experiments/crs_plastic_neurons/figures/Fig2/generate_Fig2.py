##################################
# written by: Justus Kromer
##################################

import numpy as np
import matplotlib.pyplot as plt 
import pickle

def dic_save( dic ):
    with open( dic['filename'] + '.pickle', 'wb') as f:
        pickle.dump(dic, f)

def dic_load( filename ):
    with open( filename + '.pickle', 'rb') as f:
        return pickle.load(f)


# directory in which dictionaries for mean weights are saved
directory = "data"

# This are the file endings used for "cMatrix_..." and "connectivity_matrix...".
# In an earlier version of the code inhomogeneous networks were referred to as "circular" that is why some related variables are still called "circular".
fileEnding_inhomogeneous =  "inhomogeneous" # "circular"
fileEnding_homogeneous = "homogeneous"


label_fontsize = 15
ticks_fontsize = 12

size = 80
# width of mean bars
lw_scatter = 0.3
lw_mean = 2
color_mean = "black"
color_scatter = "gray"

fig = plt.figure( figsize = (9,5))

ax_circular_acute = fig.add_subplot(221)
ax_circular_longlasting = fig.add_subplot(222)

ax_homogeneous_acute = fig.add_subplot(223)
ax_homogeneous_longlasting = fig.add_subplot(224)

xInit = 0
xSVS_01   = 1
xSVS_10   = 2
xSVS_1800 = 3
xNotShuffled = 4

#####################################################
###  initial network
#####################################################
### plot mean weights for initial networks
dic_mean_Weights_init_circular    = dic_load( directory + "/dic_mean_Weights_init_"+fileEnding_inhomogeneous )
dic_mean_Weights_init_homogeneous = dic_load( directory + "/dic_mean_Weights_init_"+fileEnding_homogeneous )

## circular network
ax_circular_acute.scatter(       [xInit*np.ones(len(dic_mean_Weights_init_circular["data"]))], dic_mean_Weights_init_circular["data"], color=color_scatter, marker='x', s = size, zorder = -1  )
ax_circular_longlasting.scatter( [xInit*np.ones(len(dic_mean_Weights_init_circular["data"]))], dic_mean_Weights_init_circular["data"], color=color_scatter, marker='x', s = size, zorder = -1  )

## homogeneous network
ax_homogeneous_acute.scatter( [xInit*np.ones(len(dic_mean_Weights_init_homogeneous["data"]))],       dic_mean_Weights_init_homogeneous["data"], color=color_scatter, marker='x', s = size, zorder = -1 )
ax_homogeneous_longlasting.scatter( [xInit*np.ones(len(dic_mean_Weights_init_homogeneous["data"]))], dic_mean_Weights_init_homogeneous["data"], color=color_scatter, marker='x', s = size, zorder = -1 )


for ax in [ax_circular_acute, ax_circular_longlasting]:    
    ax.plot( [xInit-0.5,xInit+0.5], [dic_mean_Weights_init_circular['mean'],   dic_mean_Weights_init_circular['mean']], color=color_mean, lw = lw_mean)
 
for ax in [ax_homogeneous_acute, ax_homogeneous_longlasting]:    
    ax.plot( [xInit-0.5,xInit+0.5], [dic_mean_Weights_init_homogeneous["mean"],dic_mean_Weights_init_homogeneous["mean"]], color=color_mean, lw = lw_mean)
 
#####################################################
### fixed sequence CR
#####################################################
### plot mean weights for CR with fixed sequences

# load results for fixed sequence
dic_mean_Weights_nonShuffled_circular    = dic_load( directory + "/dic_mean_Weights_nonShuffled_"+fileEnding_inhomogeneous )
dic_mean_Weights_nonShuffled_homogeneous = dic_load( directory + "/dic_mean_Weights_nonShuffled_"+fileEnding_homogeneous )
                                                        
#### fixed sequences
## acute
# weights
ax_circular_acute.scatter( xNotShuffled*np.ones( len(dic_mean_Weights_nonShuffled_circular["data_acute"]) ),        dic_mean_Weights_nonShuffled_circular["data_acute"],    color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1 )
ax_homogeneous_acute.scatter( xNotShuffled*np.ones( len(dic_mean_Weights_nonShuffled_homogeneous["data_acute"] ) ), dic_mean_Weights_nonShuffled_homogeneous["data_acute"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1 )
# mean weights
ax_circular_acute.plot(            [xNotShuffled-0.5,xNotShuffled+0.5], [dic_mean_Weights_nonShuffled_circular["mean_acute"],   dic_mean_Weights_nonShuffled_circular["mean_acute"]    ], color=color_mean, lw = lw_mean)
ax_homogeneous_acute.plot(         [xNotShuffled-0.5,xNotShuffled+0.5], [dic_mean_Weights_nonShuffled_homogeneous["mean_acute"] ,dic_mean_Weights_nonShuffled_homogeneous["mean_acute"]], color=color_mean, lw = lw_mean)
## long-lasting
# weights
ax_circular_longlasting.scatter(    xNotShuffled*np.ones( len(dic_mean_Weights_nonShuffled_circular["data_long_lasting"]) ),   dic_mean_Weights_nonShuffled_circular["data_long_lasting"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1 )
ax_homogeneous_longlasting.scatter( xNotShuffled*np.ones( len(dic_mean_Weights_nonShuffled_homogeneous["data_long_lasting"]) ),dic_mean_Weights_nonShuffled_homogeneous["data_long_lasting"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1 )
# mean weights
ax_circular_longlasting.plot(      [xNotShuffled-0.5,xNotShuffled+0.5], [dic_mean_Weights_nonShuffled_circular["mean_long_lasting"],   dic_mean_Weights_nonShuffled_circular["mean_long_lasting"]], color=color_mean, lw = lw_mean)
ax_homogeneous_longlasting.plot(   [xNotShuffled-0.5,xNotShuffled+0.5], [dic_mean_Weights_nonShuffled_homogeneous["mean_long_lasting"],dic_mean_Weights_nonShuffled_homogeneous["mean_long_lasting"]], color=color_mean, lw = lw_mean)





#####################################################
### SVS sequence
#####################################################
### plot mean weights for CR with fixed sequences

# load results for fixed sequence
#dic_Results_SVS = handle_dictionaries.dic_load( "dic_Results_SVS_rerun_Jul_13" )
dic_mean_Weights_Tshuffle_01_circular      = dic_load( directory + "/dic_mean_Weights_Tshuffle_01_"+fileEnding_inhomogeneous  )
dic_mean_Weights_Tshuffle_01_homogeneous   = dic_load( directory + "/dic_mean_Weights_Tshuffle_01_"+fileEnding_homogeneous )
dic_mean_Weights_Tshuffle_10_circular      = dic_load( directory + "/dic_mean_Weights_Tshuffle_10_"+fileEnding_inhomogeneous  )
dic_mean_Weights_Tshuffle_10_homogeneous   = dic_load( directory + "/dic_mean_Weights_Tshuffle_10_"+fileEnding_homogeneous  )
dic_mean_Weights_Tshuffle_1800_circular    = dic_load( directory + "/dic_mean_Weights_Tshuffle_1800_"+fileEnding_inhomogeneous )
dic_mean_Weights_Tshuffle_1800_homogeneous = dic_load( directory + "/dic_mean_Weights_Tshuffle_1800_"+fileEnding_homogeneous )

                                                                         
#############################################                                                                            
### circular network
#############################################
# CR RVS
ax_circular_acute.scatter( xSVS_01*np.ones( len(dic_mean_Weights_Tshuffle_01_circular["data_acute"]) ), dic_mean_Weights_Tshuffle_01_circular["data_acute"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_circular_acute.plot( [xSVS_01-0.5,xSVS_01+0.5], [dic_mean_Weights_Tshuffle_01_circular["mean_acute"],dic_mean_Weights_Tshuffle_01_circular["mean_acute"]], color=color_mean, lw = lw_mean)
# 10 sec
ax_circular_acute.scatter( xSVS_10*np.ones( len(dic_mean_Weights_Tshuffle_10_circular["data_acute"] ) ), dic_mean_Weights_Tshuffle_10_circular["data_acute"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_circular_acute.plot( [xSVS_10-0.5,xSVS_10+0.5], [dic_mean_Weights_Tshuffle_10_circular["mean_acute"],dic_mean_Weights_Tshuffle_10_circular["mean_acute"]], color=color_mean, lw = lw_mean)
# 1800 sec
ax_circular_acute.scatter( xSVS_1800*np.ones( len(dic_mean_Weights_Tshuffle_1800_circular["data_acute"]) ), dic_mean_Weights_Tshuffle_1800_circular["data_acute"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_circular_acute.plot( [xSVS_1800-0.5,xSVS_1800+0.5], [dic_mean_Weights_Tshuffle_1800_circular["mean_acute"],dic_mean_Weights_Tshuffle_1800_circular["mean_acute"]], color=color_mean, lw = lw_mean)

# CR RVS
ax_circular_longlasting.scatter( xSVS_01*np.ones( len(dic_mean_Weights_Tshuffle_01_circular["data_long_lasting"] ) ), dic_mean_Weights_Tshuffle_01_circular["data_long_lasting"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_circular_longlasting.plot( [xSVS_01-0.5,xSVS_01+0.5], [dic_mean_Weights_Tshuffle_01_circular["mean_long_lasting"],dic_mean_Weights_Tshuffle_01_circular["mean_long_lasting"]], color=color_mean, lw = lw_mean)
# 10 sec
ax_circular_longlasting.scatter( xSVS_10*np.ones( len(dic_mean_Weights_Tshuffle_10_circular["data_long_lasting"] ) ), dic_mean_Weights_Tshuffle_10_circular["data_long_lasting"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_circular_longlasting.plot( [xSVS_10-0.5,xSVS_10+0.5], [dic_mean_Weights_Tshuffle_10_circular["mean_long_lasting"],dic_mean_Weights_Tshuffle_10_circular["mean_long_lasting"]], color=color_mean, lw = lw_mean)
# 1800 sec
ax_circular_longlasting.scatter( xSVS_1800*np.ones( len(dic_mean_Weights_Tshuffle_1800_circular["data_long_lasting"] ) ), dic_mean_Weights_Tshuffle_1800_circular["data_long_lasting"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_circular_longlasting.plot( [xSVS_1800-0.5,xSVS_1800+0.5], [dic_mean_Weights_Tshuffle_1800_circular["mean_long_lasting"],dic_mean_Weights_Tshuffle_1800_circular["mean_long_lasting"]], color=color_mean, lw = lw_mean)

#############################################                                                                            
### homogeneous network
#############################################
# CR RVS
ax_homogeneous_acute.scatter( xSVS_01*np.ones( len(dic_mean_Weights_Tshuffle_01_homogeneous["data_acute"] ) ), dic_mean_Weights_Tshuffle_01_homogeneous["data_acute"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_homogeneous_acute.plot( [xSVS_01-0.5,xSVS_01+0.5], [dic_mean_Weights_Tshuffle_01_homogeneous["mean_acute"],dic_mean_Weights_Tshuffle_01_homogeneous["mean_acute"]], color=color_mean, lw = lw_mean)
# 10 sec
ax_homogeneous_acute.scatter( xSVS_10*np.ones( len(dic_mean_Weights_Tshuffle_10_homogeneous["data_acute"] ) ),dic_mean_Weights_Tshuffle_10_homogeneous["data_acute"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_homogeneous_acute.plot( [xSVS_10-0.5,xSVS_10+0.5], [dic_mean_Weights_Tshuffle_10_homogeneous["mean_acute"],dic_mean_Weights_Tshuffle_10_homogeneous["mean_acute"]], color=color_mean, lw = lw_mean)
# 1800 sec
ax_homogeneous_acute.scatter( xSVS_1800*np.ones( len(dic_mean_Weights_Tshuffle_1800_homogeneous["data_acute"] ) ), dic_mean_Weights_Tshuffle_1800_homogeneous["data_acute"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_homogeneous_acute.plot( [xSVS_1800-0.5,xSVS_1800+0.5], [dic_mean_Weights_Tshuffle_1800_homogeneous["mean_acute"],dic_mean_Weights_Tshuffle_1800_homogeneous["mean_acute"]], color=color_mean, lw = lw_mean)

# CR RVS
ax_homogeneous_longlasting.scatter( xSVS_01*np.ones( len(dic_mean_Weights_Tshuffle_01_homogeneous["data_long_lasting"] ) ), dic_mean_Weights_Tshuffle_01_homogeneous["data_long_lasting"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_homogeneous_longlasting.plot( [xSVS_01-0.5,xSVS_01+0.5], [dic_mean_Weights_Tshuffle_01_homogeneous["mean_long_lasting"],dic_mean_Weights_Tshuffle_01_homogeneous["mean_long_lasting"]], color=color_mean, lw = lw_mean)
# 10 sec
ax_homogeneous_longlasting.scatter( xSVS_10*np.ones( len(dic_mean_Weights_Tshuffle_10_homogeneous["data_long_lasting"] ) ), dic_mean_Weights_Tshuffle_10_homogeneous["data_long_lasting"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_homogeneous_longlasting.plot( [xSVS_10-0.5,xSVS_10+0.5], [dic_mean_Weights_Tshuffle_10_homogeneous["mean_long_lasting"],dic_mean_Weights_Tshuffle_10_homogeneous["mean_long_lasting"]], color=color_mean, lw = lw_mean)
# 1800 sec
ax_homogeneous_longlasting.scatter( xSVS_1800*np.ones( len(dic_mean_Weights_Tshuffle_1800_homogeneous["data_long_lasting"] ) ), dic_mean_Weights_Tshuffle_1800_homogeneous["data_long_lasting"], color=color_scatter, marker='x', s = size, lw = lw_scatter, zorder = -1)
ax_homogeneous_longlasting.plot( [xSVS_1800-0.5,xSVS_1800+0.5], [dic_mean_Weights_Tshuffle_1800_homogeneous["mean_long_lasting"],dic_mean_Weights_Tshuffle_1800_homogeneous["mean_long_lasting"]], color=color_mean, lw = lw_mean)

                                                                             
for ax in [ax_circular_acute, ax_homogeneous_acute, ax_circular_longlasting, ax_homogeneous_longlasting]:    
    
    ax.set_xticks([0,1,2,3,4])
    ax.set_xticklabels(["","","","",""], fontsize =  ticks_fontsize)
 
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
    ax.set_yticklabels(["$0$","","$0.4$","","$0.8$",""], fontsize =  ticks_fontsize)
   
    ax.axvline(0.5, color = 'gray', ls = "--", lw=0.5)
    ax.axvline(3.5, color = 'gray', ls = "--", lw=0.5)

    ax.set_ylim(0,0.65)
    ax.set_xlim(-0.7,4.6)
    
for ax_acute in [ax_circular_acute, ax_homogeneous_acute]:
    ax_acute.set_ylabel(r"$\langle w \rangle$", fontsize =  label_fontsize )

for ax in [ ax_homogeneous_acute , ax_homogeneous_longlasting ]:
    ax.set_xticklabels(["init.","$0.1$s","$10$s\n shuffled","$1800$s","           non-shuffled"], fontsize =  ticks_fontsize)
    

       
ax_circular_acute.text(-0.4,0.55, 'acute inh.', fontsize =  ticks_fontsize)    
ax_homogeneous_acute.text(-0.4,0.55, 'acute homog.', fontsize =  ticks_fontsize)  
ax_circular_longlasting.text(-0.4,0.55, 'long-lasting inh.', fontsize =  ticks_fontsize)    
ax_homogeneous_longlasting.text(-0.4,0.55, 'long-lasting homog.', fontsize =  ticks_fontsize)    
       
ax_circular_acute.text(-1.5,0.65,'A', fontsize =  1.5*label_fontsize )
ax_circular_longlasting.text(-1.2,0.65,'B', fontsize =  1.5*label_fontsize )
ax_homogeneous_acute.text(-1.5,0.65,'C', fontsize =  1.5*label_fontsize )
ax_homogeneous_longlasting.text(-1.2,0.65,'D', fontsize =  1.5*label_fontsize )

plt.savefig('Figure_2.png', bbox_inches='tight')
# plt.savefig('Figure_2.pdf', bbox_inches='tight')
# plt.savefig('Figure_2.eps', bbox_inches='tight')
