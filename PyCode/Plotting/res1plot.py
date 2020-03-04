import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

os.chdir("../../Data/time_perfomance/")

mbc_folder = "k_MBC/"
ffp_folder = "k_FFP/"

k_aver = ['k_1/', 'k_2/', 'k_4/', 'k_8/']
net_sizes = [64, 128, 256, 512, 1024, 2048, 4096, 8192]

fig, AX = plt.subplots(2,2, figsize=(9.5,8), sharey=True)

for ind, degs in enumerate(k_aver):
    mbc_points = []
    ffp_points = []
    for N in net_sizes:
        mbc_points.append(np.loadtxt(mbc_folder+degs+"N"+str(N)+".dat").mean())
        ffp_points.append(np.loadtxt(ffp_folder+degs+"N"+str(N)+".dat").mean())
    row = int(ind/2)
    col = ind%2
    ### INSET ####
    diff = np.array(mbc_points)-np.array(ffp_points)
    diff[diff<0.000] = 0.00
    #INS = inset_axes(AX[row, col], width=1.2, height=1.0, loc=2, borderpad=2)
    #INS.plot(net_sizes, diff, marker="o", mfc='mediumseagreen', color="mediumseagreen")
    #INS.set_xscale('log')
    #INS.set_yscale('log')
    #for i in ['top', 'bottom', 'left', 'right']:
    #    INS.spines[i].set_linewidth(1.5)
    #INS.tick_params(which='both', direction='in', width=2, top=True, right=True)
    #INS.set_ylabel(r'$\Delta_{diff}$', rotation='vertical', fontsize=14)
    #INS.set_xticklabels()
    ###############################################################################

    AX[row, col].plot(net_sizes, mbc_points, marker="s", color="r", linewidth=2, ms=9, label=r"$MBC$")
    AX[row, col].plot(net_sizes, ffp_points, marker="o", color="dodgerblue", linewidth=2, ms=9, label=r"$FFP$")

for row in range(2):
    for col in range(2):
        #AX[row, col].set_ylim([10**-2, 10**3.0])
        AX[row, col].spines['top'].set_linewidth(2)
        AX[row, col].spines['bottom'].set_linewidth(2)
        AX[row, col].spines['left'].set_linewidth(2)
        AX[row, col].spines['right'].set_linewidth(2)
        AX[row, col].set_yscale('log', basey=10)
        AX[row, col].set_xscale('log', basex=10)
        #AX[row, col].set_aspect('equal', 'box')
        AX[row, col].tick_params(which='both', direction='in', width=2, top=True, right=True, labelsize=16)

# Axis labels #
fig.text(0.0, 0.55, r'$runtime$ (s)', va='center', rotation='vertical', fontsize=24)
AX[1,0].set_xlabel(r'$N$', fontsize=20)
AX[1,1].set_xlabel(r'$N$', fontsize=20)
AX[0,0].set_xticklabels([])
AX[0,1].set_xticklabels([])


# Text for average degree #
# if not in log-log scale #
#AX[0,0].text(5.0, 590.0, r'$\langle k \rangle = 1$', fontsize=24)
#AX[0,1].text(5.0, 470.0, r'$\langle k \rangle = 2$', fontsize=24)
#AX[1,0].text(5.0, 270.0, r'$\langle k \rangle = 4$', fontsize=24)
#AX[1,1].text(5.0, 200.0, r'$\langle k \rangle = 8$', fontsize=24)
# if log-log scale#
AX[0,0].text(70.0, 150.0, r'$\langle k \rangle = 1$', fontsize=28)
AX[0,1].text(70.0, 150.0, r'$\langle k \rangle = 2$', fontsize=28)
AX[1,0].text(70.0, 150.0, r'$\langle k \rangle = 4$', fontsize=28)
AX[1,1].text(70.0, 150.0, r'$\langle k \rangle = 8$', fontsize=28)



AX[1,1].legend(prop={'size':'21'}, bbox_to_anchor=(0.99, 0.34), framealpha=1.)


plt.tight_layout()
plt.subplots_adjust(left=0.1)
plt.show()

    


