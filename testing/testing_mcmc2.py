import numpy as np
import corner
import matplotlib.pyplot as plt

data = np.genfromtxt("data/mcmc_pv2_chain.csv",delimiter=",")
data = data[10000:]
for x in data:
    #x[3] = np.arccos(x[3])
    #x[5] = np.exp(x[5])
    #x[6] = np.exp(x[6])
    x[3] = np.arccos(x[3])
    x[4] = np.exp(x[4])
    x[5] = np.exp(x[5])
#labels = [r"RA",r"DEC","psi","iota","phiref",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$","chip","phip"]
labels = [r"RA",r"DEC","psi","iota",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing2.pdf")
#plt.show()
plt.close()