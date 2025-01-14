import matplotlib as mpl
#mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np
import gwatpy.mcmc_routines as gmcmc
import h5py
import scipy
import emcee
injections = np.loadtxt("data/injections.csv",delimiter=',',unpack=True)
data = gmcmc.trim_thin_file("data/injection_output.csv",trim=None,ac=None)
dim = len(data[0])
print("Final samples: ",len(data))
for x in data:
    x[1] = np.arcsin(x[1])
    #x[3] = np.arccos(x[3])
    #x[1] = np.arcsin(x[1])*180./np.pi
    x[6] = np.exp(x[6])
    x[7] = np.exp(x[7])
    #x[0] = np.exp(x[0])
injections[1] = np.arcsin(injections[1])
#injections[3] = np.arccos(injections[3])
injections[6] = np.exp(injections[6])
injections[7] = np.exp(injections[7])
#injections[0] = np.exp(injections[0])
data_thinned = []
for x in np.arange(len(data)):
    if x%1 ==0:
        data_thinned.append(data[x]) 
for i in np.arange(len(data[0])):
    parameter = [x[i] for x in data]
    plt.plot(parameter)
    #plt.show()
    plt.savefig("plots/{}.pdf".format(i))
    plt.close()
ndim, nsamples = 12, len(data) 
#labels = [r"$\alpha$",r"$\sin(\delta)$",r"$\psi$",r"$\cos(\iota)$","$\phi_{ref}$","$t_c$",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$a_{1}$",r"$a_2$",r"$\cos \theta_1$",r"$\cos \theta_2$",r"$\phi_1$",r"$\phi_2$"]
labels = [r"$\alpha$",r"$\sin(\delta)$",r"$\psi$",r"$\cos \iota$","$\phi_{ref}$","$t_c$",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$m_1$",r"$m_2$"]
#labels = [r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$\beta$"]
#labels = [r"$\alpha$",r"$\sin(\delta)$",r"$\psi$",r"$\cos \iota$","$\phi_{ref}$","$t_c$",r"$D_L$",r"$m_1$",r"$m_2$",r"$\chi_{1}$",r"$\chi_2$"]
data_plot=[]
for x in data_thinned:
    #chi1 = x[2]*(x[4])
    #chi2 = x[3]*(x[5])
    #chi1 = x[2]
    #chi2 = x[3]
    #m1 = gpu.calculate_mass1_py(x[0],x[1])
    #m2 = gpu.calculate_mass2_py(x[0],x[1])
    ###chi1 = x[9]*(x[11])
    ###chi2 = x[10]*(x[12])
    ###chi1 = x[9]
    ###chi2 = x[10]
    #m1 = gpu.calculate_mass1_py(x[7],x[8])
    #m2 = gpu.calculate_mass2_py(x[7],x[8])
    #x=np.append(x,(chi1*m1 + chi2*m2 ) /(m1+m2))
    data_plot.append(x)
    #data_plot[-1] = np.insert(data_plot[-1],len(data_plot[-1]),x[7] / (1+x[8]) * np.power(x[8]/(1+x[8])**2, -3./5))
    #data_plot[-1] = np.insert(data_plot[-1],len(data_plot[-1]),x[7]*x[8] / (1+x[8]) * np.power(x[8]/(1+x[8])**2, -3./5))
data_plot = np.asarray(data_plot)
print(np.shape(data_plot))
figure = corner.corner(data_plot, labels=labels,quantiles=[.1,.5,.9], show_titles=True)
axes = np.array(figure.axes).reshape(dim,dim)
for i in np.arange(dim):
    ax = axes[i,i]
    ax.axvline(injections[i])

for yi in np.arange(dim):
    for xi in np.arange(yi):
        ax = axes[yi,xi]
        ax.axvline(injections[xi])
        ax.axhline(injections[yi])
        ax.plot(injections[xi],injections[yi])

plt.savefig("plots/mcmc_injection.pdf")
#plt.savefig("plots/mcmc_experiment.pdf")
plt.close()
