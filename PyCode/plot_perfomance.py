import numpy as np
import matplotlib.pyplot as plt

ffp_path = "../Data/time_perfomance/FFP/ffpN"
mbc_path = "../Data/time_perfomance/MBC/mbcN"

sizes = [32, 64, 128, 256, 512, 1024, 2048, 4096]
ffp_time = []
mbc_time = []

for N in sizes:
    mean_ffp = np.loadtxt(ffp_path+str(N)+".dat").mean()
    mean_mbc = np.loadtxt(mbc_path+str(N)+".dat").mean()
    ffp_time.append(mean_ffp)
    mbc_time.append(mean_mbc)

plt.plot(sizes, ffp_time, "r-o", ms=8)
plt.plot(sizes, mbc_time, "g-s", ms=8)
plt.show()