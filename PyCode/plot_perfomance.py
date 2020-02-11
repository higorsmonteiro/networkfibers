import numpy as np
import matplotlib.pyplot as plt

ffp_path = "../Data/time_perfomance/k_FFP/k_2/N"
mbc_path = "../Data/time_perfomance/k_MBC/k_2/N"

sizes = [64, 128, 256, 512, 1024]
ffp_time = []
mbc_time = []

for N in sizes:
    mean_ffp = np.loadtxt(ffp_path+str(N)+".dat").mean()
    mean_mbc = np.loadtxt(mbc_path+str(N)+".dat").mean()
    ffp_time.append(mean_ffp)
    mbc_time.append(mean_mbc)

plt.plot(sizes, ffp_time, "b-o", ms=8)
plt.plot(sizes, mbc_time, "m-s", ms=8)
plt.show()