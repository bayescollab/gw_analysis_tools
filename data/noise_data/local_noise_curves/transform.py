import numpy as np

#data = np.loadtxt("spectrum_DRSE.txt")
####KAGRA SPECIFIC######
#data = np.transpose(data)
#data = [data[0],data[-1]]
#data = np.transpose(data)
#np.savetxt("kagra_DRSE.csv",data,delimiter=',')

data = np.loadtxt("ET-0000A-18_ETDSensitivityCurveTxtFile.txt")
data = np.transpose(data)
data = [data[0],data[-1]]
data = np.transpose(data)
np.savetxt("ET-0000A-18_ETDSensitivityCurveTxtFile.csv",data,delimiter=',')
