import matplotlib as mpl
#mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np
from phenompy.utilities import calculate_mass1, calculate_mass2, mpc
from gw_analysis_tools_py import mcmc_routines_ext as mcmc
burn = True
labels = [r"$cos\iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$\sqrt{\alpha}$"]
data = np.loadtxt("data/mcmc_output_dCS.csv",delimiter=',')
data2 = np.loadtxt("data/mcmc_output_dCS_2.csv",delimiter=',')
burnlength = 30000
if burn:
    data = data[burnlength:]
    data2 = data2[burnlength:]
data2= data2[:(len(data2)-5000)]
##############################################################
#autocorr = np.loadtxt("data/auto_corr_mcmc_dCS.csv",delimiter=',')
threads = 10
segs = 50
length = len(data2)
acfile=b"data/auto_corr_mcmc_dCS_burned.csv"
mcmc.write_auto_corr_file_from_data_py(acfile, data2, length, 9, segs, .01, threads)

autocorr = np.loadtxt("data/auto_corr_mcmc_dCS_burned.csv",delimiter=',')
lengths = autocorr[0]
autocorr = autocorr[1:]

for i in np.arange(len(autocorr)):
    plt.plot(lengths,autocorr[i], label=labels[i])
plt.legend()
plt.savefig("autocorr_testing_dCS.pdf")
plt.close()
for x in data2:
    x[3] = np.exp(x[3])
    x[4] = np.exp(x[4])
for i in np.arange(9):
    parameter = [x[i] for x in data2]
    plt.plot(parameter)
    #if i == 8:
    #    plt.yscale("log")
    plt.show()
    plt.close()
ndim, nsamples = 9, len(data) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
dataplot = []
for x in data:
    dataplot.append(x)
for x in data2:
    dataplot.append(x)
    #dataplot[-1][-1]=(dataplot[-1][-1])**(1./4.)*(3e5)
dataplot = np.asarray(dataplot)
figure = corner.corner(dataplot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_dCS.pdf")
plt.close()
alphahist = []
for x in dataplot:
    alphahist.append(x[-1])
alphahist = np.asarray(alphahist)
print("Minimum root alpha: ",alphahist.min())
plt.hist(alphahist,bins=100,density=True)
plt.savefig("alpha_hist_dCS.pdf")
plt.close()
##############################################################
#data = np.loadtxt("data/mcmc_output_dCS_hot.csv",delimiter=',')
#if burn:
#    data = data[burnlength:]
##for i in np.arange(9):
##    parameter = [x[i] for x in data]
##    plt.plot(parameter)
##    if i == 8:
##        plt.yscale("log")
##    plt.show()
##    plt.close()
#
##data = data[:30000]
##data = data[:-100]
##chirpmasses = [x[1] for x in data]
##plt.plot(chirpmasses)
##plt.show()
##plt.close()
#ndim, nsamples = 9, len(data) 
##labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
#
#dataplot = []
#for x in data:
#    dataplot.append(x)
#    #dataplot[-1][-1]=(dataplot[-1][-1])**(1./4.)*(3e5)
#figure = corner.corner(dataplot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
#plt.savefig("mcmc_testing_dCS_hot.pdf")
#plt.close()

##############################################################
#autocorr = np.loadtxt("data/auto_corr_mcmc_dCS.csv",delimiter=',')
#threads = 10
#segs = 50
#length = len(data2)
#acfile=b"data/auto_corr_mcmc_dCS_burned.csv"
#mcmc.write_auto_corr_file_from_data_py(acfile, data2, length, 9, segs, .01, threads)
#
#autocorr = np.loadtxt("data/auto_corr_mcmc_dCS_burned.csv",delimiter=',')
#lengths = autocorr[0]
#autocorr = autocorr[1:]
#
#for i in np.arange(len(autocorr)):
#    plt.plot(lengths,autocorr[i], label=labels[i])
#plt.legend()
#plt.savefig("autocorr_testing_dCS.pdf")
#plt.close()

#autocorr = np.loadtxt("data/auto_corr_mcmc_dCS_GPU.csv",delimiter=',')
#lengths = autocorr[0]
#autocorr = autocorr[1:]
#
#for i in np.arange(len(autocorr)):
#    plt.plot(lengths,autocorr[i], label=labels[i])
#plt.legend()
#plt.savefig("autocorr_testing_dCS_GPU.pdf")
#plt.close()

##############################################################
data = np.loadtxt("data/mcmc_output_dCS.csv",delimiter=',')
if burn:
    data = data[10000:]
def chieff(m1, m2, spin1, spin2,):
    return (m1*spin1 + m2*spin2)/(m1+m2)
datatransform= []
for x in data:
    chirpm = x[4]
    symmratio = x[5]
    m1 = calculate_mass1(chirpm, symmratio)
    m2 = calculate_mass2(chirpm, symmratio)
    datatransform.append([m1,m2, chieff(m1,m2,x[6],x[7])])
datatransform = np.asarray(datatransform)

ndim, nsamples = 3, len(datatransform) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
labels = [r"$M_1$",r"$M_2$",r"$\chi_{eff}$"]

figure = corner.corner(datatransform, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_transform_dCS.pdf")
plt.close()
