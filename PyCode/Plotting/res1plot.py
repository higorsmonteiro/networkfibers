import os
import numpy as np
import matplotlib.pyplot as plt

os.chdir("../../Data/time_perfomance/")

mbc_folder = "k_MBC/"
ffp_folder = "k_FFP/"

k_aver = ['k_1/', 'k_2/', 'k_4/', 'k_8/']
net_sizes = [64, 128, 256, 512, 1024, 2048, 4096, 8192]

fig, AX = plt.subplots(2,2, figsize=(10,8))

for ind, degs in enumerate(k_aver):
    mbc_points = []
    ffp_points = []
    for N in net_sizes:
        mbc_points.append(np.loadtxt(mbc_folder+degs+"N"+str(N)+".dat").mean())
        ffp_points.append(np.loadtxt(ffp_folder+degs+"N"+str(N)+".dat").mean())
    row = int(ind/2)
    col = ind%2
    AX[row, col].plot(net_sizes, mbc_points, marker="s", color="r", linewidth=2, ms=8, label="MBC")
    AX[row, col].plot(net_sizes, ffp_points, marker="o", color="dodgerblue", linewidth=2, ms=8, label="FFP")

for row in range(2):
    for col in range(2):
        AX[row, col].spines['top'].set_linewidth(2)
        AX[row, col].spines['bottom'].set_linewidth(2)
        AX[row, col].spines['left'].set_linewidth(2)
        AX[row, col].spines['right'].set_linewidth(2)
        #AX[row, col].set_yscale('log', basey=10)
        #AX[row, col].set_xscale('log', basex=10)
        AX[row, col].tick_params(which='both', direction='in', width=2, top=True, right=True, labelsize=13)

# Axis labels #
AX[0,0].set_ylabel(r'$runtime$ (s)', fontsize=15)
AX[1,0].set_ylabel(r'$runtime$ (s)', fontsize=15)
AX[1,0].set_xlabel(r'$N$', fontsize=15)
AX[1,1].set_xlabel(r'$N$', fontsize=15)
AX[0,0].set_xticklabels([])
AX[0,1].set_xticklabels([])

# Text for average degree #
# if not in log-log scale #
AX[0,0].text(5.0, 590.0, r'$\langle k \rangle = 1$', fontsize=24)
AX[0,1].text(5.0, 470.0, r'$\langle k \rangle = 2$', fontsize=24)
AX[1,0].text(5.0, 270.0, r'$\langle k \rangle = 4$', fontsize=24)
AX[1,1].text(5.0, 200.0, r'$\langle k \rangle = 8$', fontsize=24)
# if log-log scale#



AX[0,1].legend(prop={'size':'16'}, bbox_to_anchor=(1.02, 1.02))


plt.tight_layout()
plt.show()

    


