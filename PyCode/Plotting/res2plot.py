import os
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt

os.chdir("../../Data/time_perfomance/")

mbc_folder = "MBC/"
ffp_folder = "FFP/"

EDGETYPES = ['INFOTYPE_1/', 'INFOTYPE_2/', 'INFOTYPE_4/', 'INFOTYPE_8/']
k_aver = ['k_1/', 'k_2/', 'k_4/', 'k_8/']
net_sizes = [64, 128, 256, 512, 1024, 2048]

markers = ['o', 's', '^', 'p']
colors = ['red', 'mediumseagreen', 'dodgerblue', 'hotpink']
labels = [r"$K_{type} = 1$", r"$K_{type} = 2$", r"$K_{type} = 4$", r"$K_{type} = 8$"]

fig, AX = plt.subplots(4,2, figsize=(9.2,10), sharey=True)

for ind, info in enumerate(EDGETYPES):
    for ind_deg, degs in enumerate(k_aver):
        mbc_points = []
        ffp_points = []
        for N in net_sizes:
            mbc_points.append(np.loadtxt(mbc_folder+info+degs+"N"+str(N)+".dat").mean())
            ffp_points.append(np.loadtxt(ffp_folder+info+degs+"N"+str(N)+".dat").mean())
        AX[ind_deg, 0].plot(net_sizes, mbc_points, marker=markers[ind], color=colors[ind], linewidth=2, ms=8, label=labels[ind])
        AX[ind_deg, 1].plot(net_sizes, ffp_points, marker=markers[ind], color=colors[ind], linewidth=2, ms=8, label=labels[ind])

for row in range(len(k_aver)-1):
    AX[row,0].set_xticklabels([])
    AX[row,1].set_xticklabels([])

for row in range(len(k_aver)):
    for col in range(2):
        AX[row, col].spines['top'].set_linewidth(2)
        AX[row, col].spines['bottom'].set_linewidth(2)
        AX[row, col].spines['left'].set_linewidth(2)
        AX[row, col].spines['right'].set_linewidth(2)
        AX[row, col].set_yscale('log', basey=10)
        AX[row, col].set_xscale('log', basex=10)
        AX[row, col].tick_params(which='both', direction='in', width=2, top=True, right=True, labelsize=13)

for row in range(len(k_aver)):
    AX[row, 0].set_ylabel(r'$runtime (s)$', fontsize=15)

AX[0,0].set_title(r"$MBC$", fontsize=16)
AX[0,1].set_title(r"$FFP$", fontsize=16)
fig.text(0.51, 0.005, r'$N$', ha='center', fontsize=24)
#fig.text(0.0, 0.5, r'$runtime (s)$', va='center', rotation='vertical', fontsize=24)

#AX[0,0].text(0.0, 15.0, r'$\langle k \rangle = 1$', fontsize=18)
#AX[1,0].text(0.0, 15.0, r'$\langle k \rangle = 2$', fontsize=18)
#AX[2,0].text(0.0, 15.0, r'$\langle k \rangle = 4$', fontsize=18)
#AX[3,0].text(0.0, 15.0, r'$\langle k \rangle = 8$', fontsize=18)

fig.text(0.0, 0.155, r'$\langle k \rangle = 8$', va='center', fontsize=20)
fig.text(0.0, 0.39, r'$\langle k \rangle = 4$', va='center', fontsize=20)
fig.text(0.0, 0.622, r'$\langle k \rangle = 2$', va='center', fontsize=20)
fig.text(0.0, 0.86, r'$\langle k \rangle = 1$', va='center', fontsize=20)

#AX[0,0].text(65.0, 15.0, r'$\langle k \rangle = 1$', fontsize=18)
#AX[1,0].text(65.0, 15.0, r'$\langle k \rangle = 2$', fontsize=18)
#AX[2,0].text(65.0, 15.0, r'$\langle k \rangle = 4$', fontsize=18)
#AX[3,0].text(65.0, 15.0, r'$\langle k \rangle = 8$', fontsize=18)

AX[0,1].legend(prop={'size':'14'}, bbox_to_anchor=(1.02, 1.035), markerscale=1.5)
#mat.rcParams['legend.markerscale'] = 0.0

plt.tight_layout()
plt.subplots_adjust(left=0.20, bottom=0.06)
plt.show()