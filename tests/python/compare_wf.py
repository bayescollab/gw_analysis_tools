import numpy as np
import matplotlib.pyplot as plt
from time import time
# import astropy.cosmology as cosmology
# from astropy.coordinates import Distance
# from astropy import units as u
# import astropy.constants as consts
# from pycbc.types import frequencyseries
# from pycbc.waveform import get_fd_waveform_sequence, get_fd_waveform, get_td_waveform
# from pycbc.waveform import fd_approximants, utils
# from pycbc.types import Array

'''Code for comparing phases with different ppE modifications'''
# data = np.loadtxt("data/error/b-5.csv", delimiter = ',').T
# b = "minus5"
# label1 = "NRT"
# label2 = "ppE_NRT"
#
# phase1 = data[4]
# phase2 = data[5]
#
# plt.plot(phase2 - phase1, label = "phase difference")
# plt.legend()
# plt.savefig("plots/ppE_mod/"+b+"_phasediff.pdf")
# plt.close()
#
# plt.plot(phase1, label=label1)
# plt.plot(phase2, label=label2)
# plt.legend()
# plt.savefig("plots/ppE_mod/"+b+"_phase.pdf")
# plt.close()

'''Code for comparing waveform models for the error calculation'''
# data = np.loadtxt("data/error/betazero/IMRPhenomD_comp_IMRPhenomD_NRT.csv", delimiter=',').T
# models = "NRT_comp_D"
# label1 = "D"
# label2 = "NRT"
data = np.loadtxt("data/error/b5/IMRPhenomD_NRT_comp_ppE_IMRPhenomD_NRT_IMR.csv", delimiter=',').T
models = "b5/NRT_comp_ppE_IMR"
label1 = "NRT"
label2 = "ppE_NRT"
# data = np.loadtxt("data/error/b7/IMRPhenomD_comp_ppE_IMRPhenomD_IMR.csv", delimiter=',').T
# models = "b7/D_comp_ppE"
# label1 = "D"
# label2 = "ppE_D"

freqs = data[6]
wf1 = data[0] + 1j*data[1]
wf2 = data[2] + 1j*data[3]
amp1 = data[0]*data[0] + data[1]*data[1]
amp2 = data[2]*data[2] + data[3]*data[3]
phase1 = data[4]
phase2 = data[5]

plt.loglog(freqs, abs(wf1), label=label1)
plt.loglog(freqs, abs(wf2), label=label2)
plt.axvline(x=2608.32, color = 'black', linestyle = '--')
plt.legend()
plt.title("Waveform")
plt.savefig("plots/wf_comparison/"+models+"_wf.pdf")
plt.close()

plt.plot(freqs, phase1, label=label1)
plt.plot(freqs, phase2, label=label2)
plt.axvline(x=2608.32, color = 'black', linestyle = '--')
plt.legend()
plt.title("Phase comparison")
plt.savefig("plots/wf_comparison/"+models+"_phase.pdf")
plt.close()

plt.plot(freqs, amp1, label=label1)
plt.plot(freqs, amp2, label=label2)
plt.axvline(x=2608.32, color = 'black', linestyle = '--')
plt.legend()
plt.title("Amplitude comparison")
plt.savefig("plots/wf_comparison/"+models+"_amp.pdf")
plt.close()

plt.plot(freqs, phase2 - phase1, label = "phase difference")
plt.axvline(x=2608.32, color = 'black', linestyle = '--')
plt.legend()
plt.title("Phase difference")
plt.savefig("plots/wf_comparison/"+models+"_phasediff.pdf")
plt.close()

plt.plot(freqs, (amp1 - amp2)/amp2, label = "amp difference")
plt.axvline(x=2608.32, color = 'black', linestyle = '--')
plt.legend()
plt.title("Amplitude difference")
plt.savefig("plots/wf_comparison/"+models+"_ampdiff.pdf")
plt.close()

''' Code for comparing waveform models for the error calculation (deprecated file structure)'''
# NRTc = np.loadtxt("data/error/waveform_data.csv", delimiter=",", unpack=True)
# ppENRTc = np.loadtxt("data/error/bgr_waveform_data.csv", delimiter=",", unpack=True)
# plt.plot(NRTc[0],NRTc[1],label="NRT c")
# plt.plot(NRTc[0],ppENRTc[1],label="ppE NRT c", linestyle='-.')
# plt.legend()
# plt.savefig("plots/wfd.pdf")
# plt.close()
#
# ampNRTc = NRTc[1]*NRTc[1] + NRTc[2]*NRTc[2]
# ampppEc = ppENRTc[1]*ppENRTc[1] + ppENRTc[2]*ppENRTc[2]
#
# plt.plot(NRTc[0], ampNRTc, label ="NRT")
# plt.plot(NRTc[0], ampppEc, label ="ppE")
# plt.legend()
# plt.savefig("plots/ampd.pdf")
# plt.close()
#
# plt.plot(NRTc[0],(ampNRTc-ampppEc)/ampppEc,label="cross")
# plt.legend()
# plt.savefig("plots/ampddiff.pdf")
# plt.close()
#
# phaseNRTc = NRTc[3]
# phaseppEc = ppENRTc[3]
#
# plt.plot(NRTc[0],phaseNRTc,label="NRT c")
# plt.plot(NRTc[0],phaseppEc,label="D c")
# plt.xlim(0,250)
# plt.ylim(-100, 100)
# plt.legend()
# plt.savefig("plots/phased.pdf")
# plt.close()
#
# plt.plot(NRTc[0],(phaseNRTc-phaseppEc), label="cross")
# plt.legend()
# plt.savefig("plots/phaseddiff.pdf")
# plt.close()


''' Code for comparing gwat waveforms to LAL waveforms '''
# gwatp = np.loadtxt("data/hpgwat.csv",delimiter=',',unpack=True)
# lalp = np.loadtxt("data/hpLAL.csv",delimiter=',',unpack=True)
# freq = np.loadtxt("data/freqs.csv",delimiter=',',unpack=True)
# plt.plot(freq,gwatp[0],label="gwat p")
# plt.plot(freq,lalp[0],label="LAL p",linestyle='-.')
# plt.legend()
# plt.savefig("plots/wfd.pdf")
# #plt.show()
# plt.close()

# ampgwatp = gwatp[0]*gwatp[0] + gwatp[1]*gwatp[1]
# ampLALp = lalp[0]*lalp[0] + lalp[1]*lalp[1]
# plt.semilogy(freq,ampgwatp,label="gwat p")
# plt.semilogy(freq,ampLALp,label="LAL p")
# plt.legend()
# plt.savefig("plots/ampd.pdf")
# #plt.show()
# plt.close()
#
# plt.plot(freq,(ampgwatp-ampLALp)/ampLALp,label="plus")
# #plt.plot(freq,(ampgwatp-ampLALp),label="plus")
# #plt.xlim([0,400])
# plt.legend()
# plt.savefig("plots/ampddiff.pdf")
# #plt.show()
# plt.close()

# delta_f = freq[1]-freq[0]
# gwatseriesplus = frequencyseries.FrequencySeries(gwatp[0]+ 1j*gwatp[1],delta_f=delta_f)
# phasegwatp = utils.phase_from_frequencyseries(gwatseriesplus)
# LALseriesplus = frequencyseries.FrequencySeries(lalp[0]+ 1j*lalp[1],delta_f=delta_f)
# phaseLALp = utils.phase_from_frequencyseries(LALseriesplus)
# gwatpp = np.loadtxt("data/ppgwat.csv",delimiter=',',unpack=True)
# #phasegwatp = np.arctan(gwatp[1]/gwatp[0] )
# #phaseLALp = np.arctan(lalp[1]/lalp[0] )
# #phasegwatc = np.arctan(gwatc[1]/gwatc[0] )
# #phaseLALc = np.arctan(lalc[1]/lalc[0] )
# plt.plot(freq,phasegwatp,label="gwat p")
# plt.plot(freq,phaseLALp,label="LAL p")
# #plt.plot(freq,gwatpp+(phasegwatp[0]-gwatpp[0]),label="gwat p -- native")
# plt.legend()
# plt.savefig("plots/phased.pdf")
# #plt.show()
# plt.close()
#
# plt.plot(freq,(phasegwatp-phaseLALp)/phaseLALp,label="plus")
# plt.plot(freq,(phasegwatp-phaseLALp),label="plus")
# plt.plot(freq,(phasegwatc-phaseLALc),label="cross")
# plt.plot(freq,(phasegwatp-phaseLALp),label="plus")
# plt.legend()
# plt.savefig("plots/phaseddiff.pdf")
# #plt.show()
# plt.close()
