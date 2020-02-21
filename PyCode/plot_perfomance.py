import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set()

k_folders = ['k_4/', 'k_8/', 'k_16/', 'k_32/']
sizes = [64, 128, 256, 512, 1024, 2048, 4096]
ffp_path = "../Data/time_perfomance/k_FFP/"
mbc_path = "../Data/time_perfomance/k_MBC/"

ffp_time = [[], [], [], []]
mbc_time = [[], [], [], []]

for f_index, ff in enumerate(k_folders):
    for N in sizes:
        mean_ffp = np.loadtxt(ffp_path+ff+"N"+str(N)+".dat").mean()
        mean_mbc = np.loadtxt(mbc_path+ff+"N"+str(N)+".dat").mean()
        ffp_time[f_index].append(mean_ffp)
        mbc_time[f_index].append(mean_mbc)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(10,5))
ax1.plot(sizes, ffp_time[0], "m-o", ms=8, label="ffp")
ax1.plot(sizes, mbc_time[0], "g-s", ms=8, label="mbc")
ax2.plot(sizes, ffp_time[1], "m-o", ms=8, label="ffp")
ax2.plot(sizes, mbc_time[1], "g-s", ms=8, label="mbc")
ax3.plot(sizes, ffp_time[2], "m-o", ms=8, label="ffp")
ax3.plot(sizes, mbc_time[2], "g-s", ms=8, label="mbc")
ax4.plot(sizes, ffp_time[3], "m-o", ms=8, label="ffp")
ax4.plot(sizes, mbc_time[3], "g-s", ms=8, label="mbc")
#ax5.plot(sizes, ffp_time[4], "m-o", ms=8, label="ffp")
#ax5.plot(sizes, mbc_time[4], "g-s", ms=8, label="mbc")

ax1.set_ylim([-10,max(mbc_time[0])+5])
ax2.set_ylim([-10,max(mbc_time[0])+5])
ax3.set_ylim([-10,max(mbc_time[0])+5])
ax4.set_ylim([-10,max(mbc_time[0])+5])
#ax5.set_ylim([-10,max(mbc_time[0])+5])
ax1.set_xticklabels(sizes)
ax2.set_xticklabels(sizes)
ax3.set_xticklabels(sizes)
ax4.set_xticklabels(sizes)
#ax5.set_xticklabels(sizes)

ax1.set_ylabel("time(s)", fontsize=18)
ax1.set_xlabel("N", fontsize=18)
ax2.set_xlabel("N", fontsize=18)
ax3.set_xlabel("N", fontsize=18)
ax1.set_title(r"$\langle k \rangle = 4$", fontsize=18)
ax2.set_title(r"$\langle k \rangle = 8$", fontsize=18)
ax3.set_title(r"$\langle k \rangle = 16$", fontsize=18)
ax4.set_title(r"$\langle k \rangle = 32$", fontsize=18)
#ax5.set_title(r"$\langle k \rangle = 16$", fontsize=18)

#ax1.set_xscale('log')
#ax1.set_yscale('log')
#ax2.set_xscale('log')
#ax2.set_yscale('log')
#ax3.set_xscale('log')
#ax3.set_yscale('log')
#ax4.set_xscale('log')
#ax4.set_yscale('log')
#ax5.set_xscale('log')
#ax5.set_yscale('log')
plt.legend()
plt.show()