import numpy as np
import matplotlib.pyplot as plt


tickFontsize = 12
labelFontsize = 15

fig = plt.figure()

ax_Kuramoto_circular  = fig.add_subplot(221)
ax_weight_circular = fig.add_subplot(223)
ax_Kuramoto_homogeneous  = fig.add_subplot(222)
ax_weight_homogeneous = fig.add_subplot(224)

mwValues = np.round( np.arange(0.0,1.1,0.1) , 1 )
colors =   ["0.2","0.4","0.6","0.8","red"]
mwValues = [0.0,   0.2,  0.4,  0.6] # , 0.8]
stepsToSec = 0.001 # sec
seed = 12

### inhomogeneous network
for kmw in range( len( mwValues ) ):
    
    mw = mwValues[kmw]
    print("inhomogeneous network", mw)
    
    # plot trajectory of Kuramoto order parameter
    try:
        KuramotoOrderParameter_scaled = np.load( "data/KuramotoOrderParameter_inhomogeneous_network_mw_"+str(mw)+"_seed_"+str(seed)+".npy" )
        ax_Kuramoto_circular.plot( KuramotoOrderParameter_scaled[:,0] , KuramotoOrderParameter_scaled[:,1], color = colors[kmw] )
    except: pass

    # plot trajectory of weight 
    try:

        weightFile = np.load( "data/inhomogeneous_network_mwTrajectory_w_"+str(mw)+"_seed_"+str(seed)+".npy" )

        # load adj       
        # adjust time to second
        x = 0.0001*weightFile[:,0]
        y = weightFile[:,1]

        ax_weight_circular.plot( x ,y, color = colors[kmw] )

    except: pass

    
### homogeneous network
for kmw in range( len( mwValues ) ):
    
    mw = mwValues[kmw]
    print("homogeneous network", mw)
    
    # plot trajectory of Kuramoto order parameter
    try:                                         
        KuramotoOrderParameter_scaled = np.load( "data/KuramotoOrderParameter_homogeneous_network_mw_"+str(mw)+"_seed_"+str(seed)+".npy" )
        ax_Kuramoto_homogeneous.plot( KuramotoOrderParameter_scaled[:,0] , KuramotoOrderParameter_scaled[:,1], color = colors[kmw] )
    except: pass

    # plot trajectory of weight 
    try:
        weightFile = np.load( "data/homogeneous_network_mwTrajectory_w_"+str(mw)+"_seed_"+str(seed)+".npy" )

        # load adj
        # adjust time to second
        x = 0.0001*weightFile[:,0]
        y = weightFile[:,1]

        ax_weight_homogeneous.plot( x ,y, color = colors[kmw] )
    except: pass    
      
    
for ax in [ax_Kuramoto_circular, ax_weight_circular, ax_Kuramoto_homogeneous, ax_weight_homogeneous]:
#     ax.set_xticks([1200,3000,4800,6600,8400,10200,12000, 13800])
#     ax.set_xticklabels(["","","","","","","",""], fontsize = tickFontsize)
    
#     ax.set_xlim(0,12800)
    ax.set_xticks([100,1000,10000])
    # ax.set_xticklabels(["","",""], fontsize = tickFontsize)
    ax.set_xticklabels(["$100$","$1000$","$10000$"], fontsize = tickFontsize)
    ax.set_xlim(10,12800)
    ax.set_ylim(0,1.05)
    ax.set_xscale('log')
    

for ax in [ax_Kuramoto_circular, ax_Kuramoto_homogeneous]:
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    ax.set_yticklabels(["","","","","",""], fontsize = tickFontsize)
    ax.set_xticklabels(["","",""], fontsize = tickFontsize)
    
for ax in [ax_weight_circular, ax_weight_homogeneous]:
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6])
    ax.set_yticklabels(["","","","","","",""], fontsize = tickFontsize)
    ax.set_ylim(0,0.6)   
    
        
for ax in [ax_weight_circular,ax_weight_homogeneous]:
    # ax.set_xticklabels(["","$0$","","$1$","","$2$","","$3$"], fontsize = tickFontsize)
    # ax.set_xticklabels(["$100$","$1000$","$10000$"], fontsize = tickFontsize)
    ax.set_xlabel("$t$ (sec)", fontsize = labelFontsize )
    

ax_Kuramoto_circular.set_yticklabels(["$0$","","","","","$1$"], fontsize = tickFontsize)
ax_weight_circular.set_yticklabels(["$0$","","","","","$0.5$",""], fontsize = tickFontsize)
    
ax_Kuramoto_circular.set_ylabel(r"$\rho$", fontsize = labelFontsize )
ax_weight_circular.set_ylabel(r"$\langle w \rangle$", fontsize = labelFontsize )


ax_Kuramoto_circular.text( 3, 1.1, "A", fontsize = 1.3*labelFontsize )
ax_Kuramoto_homogeneous.text( 3, 1.1, "B", fontsize = 1.3*labelFontsize )
ax_weight_circular.text( 3, 0.6, "C", fontsize = 1.3*labelFontsize )
ax_weight_homogeneous.text( 3, 0.6, "D", fontsize = 1.3*labelFontsize )

ax_Kuramoto_circular.set_title("inhomogeneous", fontsize = labelFontsize)
ax_Kuramoto_homogeneous.set_title("homogeneous", fontsize = labelFontsize)

# plt.savefig('Figure_1.svg', bbox_inches='tight')
plt.savefig('Figure_1.png', bbox_inches='tight')
# plt.savefig('Figure_1.pdf', bbox_inches='tight')
# plt.savefig('Figure_1.eps', bbox_inches='tight')
