import numpy as np
import matplotlib.pyplot as plt
import h5py
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import quad


## Comparison tests
#import tensorflow_probability as tfp



class MCMCOutput:
    filename = ""
    outputFile = None
    betas = None
    betaSchedule = None
    temps = None
    ACVals = None
    trimLengths = None
    dataSets = None
    ensembleSize = 0
    ensembleN = 0
    RJ = False
    evidence = None
    def __init__(self, MCMCOutputFile):
        self.filename = MCMCOutputFile
        self.outputFile = h5py.File(self.filename)
        self.betas = np.array(self.outputFile["MCMC_METADATA"]["CHAIN BETAS"])
        with np.errstate(divide='ignore'):
            self.temps = 1./self.betas
        if "STATUS" in self.outputFile["MCMC_OUTPUT"].keys():
            self.RJ=True
        if not self.RJ:
            self.ACVals = np.array(self.outputFile["MCMC_METADATA"]["AC VALUES"])
        self.trimLengths = np.array(self.outputFile["MCMC_METADATA"]["SUGGESTED TRIM LENGTHS"])


        self.betaSchedule = np.flip(np.unique(np.array(self.outputFile["MCMC_METADATA"]["CHAIN BETAS"])))
        self.ensembleSize = int(len(self.betaSchedule))
        self.ensembleN = int(len(self.betas)/self.ensembleSize)
        if "EVIDENCE" in self.outputFile["MCMC_METADATA"].keys():
            self.evidence = self.outputFile["MCMC_METADATA"]["EVIDENCE"][0]
        return

    #def unpackMCMCData(self,betaID=0, trim=None, thin=None):
    #
    #    if betaID > self.ensembleN:
    #        print("Supplied a betaID larger than the number of ensembles!")
    #        return None, None, None
    #
    #    chainIDs = np.arange(self.chainIndex(0,betaID) , self.chainIndex(self.ensembleN,betaID))
    #
    #
    #    trim_local=0
    #    thin_local=1

    #    if trim is None  and betaID ==0:
    #        trim_local = self.trimLengths[chainIDs[0]]
    #    if thin is None and betaID ==0:
    #        thin_local = np.amax(self.ACVals[chainIDs[0]][:])
    #    if thin_local == 0 :
    #        thin_local=1
    #    data = self.outputFile["MCMC_OUTPUT"]["CHAIN {}".format(chainIDs[0])][trim_local::thin_local]
    #    for x in chainIDs[1:]:
    #        if trim is None and betaID ==0:
    #            trim_local = self.trimLengths[x]
    #        if thin is None and betaID ==0:
    #            thin_local = np.amax(self.ACVals[x][:])
    #        if thin_local == 0:
    #            thin_local=1
    #        data = np.insert(data,-1, self.outputFile["MCMC_OUTPUT"]["CHAIN {}".format(x)][trim_local::thin_local],axis=0)
    #    return data
    def calculateEvidence(integrationSizeCap=None):
        integratedLikelihoods = np.ones(len(betaSchedule))
        for i in np.arange(len(betaSchedule)):
            d = unpackMCMCData(betaID=i,trim=0,thin=1,sizeCap=integrationSizeCap)
            integratedLikelihoods[i] = np.sum(d["logL"])/len(d["logL"])
        f_ = interp1d(betaSchedule, integratedLikelihoods,kind='cubic')
        result = quad(f_,0,1)
        return result[0]

    def unpackMCMCData(self,betaID=0, trim=None, thin=None,sizeCap=None):
        if betaID > self.ensembleSize:
            print("Supplied a betaID larger than the number of ensembles!")
            return None, None, None

        #chainIDs = np.arange(ensembleSize*betaID , ensembleSize*betaID + ensembleN)
        chainIDs = np.arange(self.chainIndex(0,betaID) , self.chainIndex(self.ensembleN,betaID))

        trim_local = 0
        thin_local = 1
        if trim is not None:
            trim_local = trim
        elif betaID ==0 and not self.RJ:
            trim_local = self.trimLengths[chainIDs[0]]

        if thin is not None:
            thin_local = thin
        elif betaID ==0 and not self.RJ:
            thin_local = np.amax(self.ACVals[chainIDs[0]][:])

        self.selectedData = None
        data = self.outputFile["MCMC_OUTPUT"]["CHAIN {}".format(chainIDs[0])][trim_local::thin_local]
        logl = self.outputFile["MCMC_OUTPUT/LOGL_LOGP"]["CHAIN {}".format(chainIDs[0])][trim_local::thin_local,0]
        logp = self.outputFile["MCMC_OUTPUT/LOGL_LOGP"]["CHAIN {}".format(chainIDs[0])][trim_local::thin_local,1]
        status = None
        model_status = None
        if self.RJ:
            status = self.outputFile["MCMC_OUTPUT/STATUS"]["CHAIN {}".format(chainIDs[0])][trim_local::thin_local]
            if "MCMC_OUTPUT/MODEL_STATUS" in self.outputFile.keys():
                model_status = self.outputFile["MCMC_OUTPUT/MODEL_STATUS"]["CHAIN {}".format(chainIDs[0])][trim_local::thin_local]
        for x in chainIDs[1:]:
            if trim is None and betaID ==0 and not self.RJ:
                trim_local = self.trimLengths[x]
            if thin is None and betaID ==0 and not self.RJ:
                thin_local = np.amax(self.ACVals[x][:])
                #print(thin_local)
            if thin_local == 0:
                thin_local=1
            data = np.insert(data,-1, self.outputFile["MCMC_OUTPUT"]["CHAIN {}".format(x)][trim_local::thin_local],axis=0)
            logl = np.insert(logl,-1, self.outputFile["MCMC_OUTPUT/LOGL_LOGP"]["CHAIN {}".format(x)][trim_local::thin_local,0],axis=0)
            logp = np.insert(logp,-1, self.outputFile["MCMC_OUTPUT/LOGL_LOGP"]["CHAIN {}".format(x)][trim_local::thin_local,1],axis=0)
            if self.RJ:
                status = np.insert(status,-1, self.outputFile["MCMC_OUTPUT/STATUS"]["CHAIN {}".format(x)][trim_local::thin_local],axis=0)
                if "MCMC_OUTPUT/MODEL_STATUS" in self.outputFile.keys():
                    model_status = np.insert(model_status,-1, self.outputFile["MCMC_OUTPUT/MODEL_STATUS"]["CHAIN {}".format(x)][trim_local::thin_local],axis=0)

        if sizeCap is not None:
            if data.shape[0] > sizeCap:
                local_trim = int( data.shape[0]/sizeCap)
                data =  data[::local_trim ]
                logl =  logl[::local_trim ]
                logp =  logp[::local_trim ]
                if status is not None:
                    status =  status[::local_trim ]
                if model_status is not None:
                    model_status =  model_status[::local_trim ]
        d = {"logL":logl,"logP":logp}

        labels = list(["Parameter {}".format(x) for x in np.arange(data.shape[1])])
        d["data"] = pd.DataFrame(data,columns=labels)
        if status is not None:
            d["status"] = pd.DataFrame(status,columns=labels)
        if model_status is not None:
            d["model_status"] = model_status
        d["beta"] =self.betaSchedule[betaID]
        if(self.evidence is not None):
            d["Evidence"] = self.evidence
        return d

    def chainIndex(self,ensemble, betaN):
        return ensemble+betaN*self.ensembleN


# Takes the dataObjs for MCMC runs and returns a N-dimensional array for the average R value for each dimension
def gelmanRubinStatistic(dataObjs):
    minEnsembleN = np.amin( np.array( [d.ensembleN for d in dataObjs] ) )
    chainIDs = []
    for d in dataObjs:
        chainIDs.append(np.arange(d.chainIndex(0,0) , d.chainIndex(minEnsembleN,0)))

    chains = []
    lengths = np.ones( ( len(dataObjs), minEnsembleN ))
    for i, d in enumerate(dataObjs):
        chains.append([])
        for j, c in enumerate(chainIDs[i]):
            trim, thin = int(d.trimLengths[c]), int(np.amax(d.ACVals[c]))
            chains[-1].append(d.outputFile["MCMC_OUTPUT"]["CHAIN {}".format(c)][trim::thin])
            lengths[i,j] = len(chains[-1][-1])
    minLength = int(np.amin(lengths))
    dim = len(chains[0][0][0])
    R = []
    #temp = np.ones( (minLength, len(dataObjs), dim))
    #tempv2 = np.zeros( dim)
    for i in np.arange(minEnsembleN):
        chainMeans = np.ones( ( len(dataObjs),dim ) )
        chainVars = np.ones((len(dataObjs), dim) )
        for j in np.arange( len(dataObjs)):
            for k in np.arange(dim):
                chainMeans[j,k] = np.mean(chains[j][i][:minLength, k],axis=0)
                chainVars[j,k] = np.var(chains[j][i][:minLength, k],axis=0,ddof=1)
                #temp[:,j,k] = chains[j][i][:minLength, k]
        totalMean = np.mean(chainMeans,axis=0)
        totalVar = minLength*np.var(chainMeans,axis=0,ddof=1)
        W = np.mean(chainVars, axis=0)
        #R.append( np.sqrt(( (minLength - 1.)/ minLength * W + 1./minLength * totalVar )/W ))
        R.append( ( (minLength - 1.)/ minLength * W + 1./minLength * totalVar )/W )
        #tfprhat = tfp.mcmc.diagnostic.potential_scale_reduction(temp,independent_chain_ndims=1)
        #print(R[-1])
        #print(tfprhat)
        #tempv2 +=tfprhat
    #tempv2/=minEnsembleN
    #print("tensorflow final: ",tempv2)
    R = np.mean(R, axis=0)

    return R


def chainIndex(ensemble, betaN, ensembleN):
    return ensemble+betaN*ensembleN

def RJMCMC_unpack_file(filename,betaID=0):
    f = h5py.File(filename,'r')
    betas = np.array(f["MCMC_METADATA"]["CHAIN BETAS"])
    betaSchedule = np.flip(np.unique(np.array(f["MCMC_METADATA"]["CHAIN BETAS"])))
    ensembleSize = int(len(betaSchedule))
    ensembleN = int(len(betas)/ensembleSize)

    if betaID > ensembleSize:
        print("Supplied a betaID larger than the number of ensembles!")
        return None, None, None

    #chainIDs = np.arange(ensembleSize*betaID , ensembleSize*betaID + ensembleN)
    chainIDs = np.arange(chainIndex(0,betaID,ensembleN) , chainIndex(ensembleN,betaID,ensembleN))

    chains = list(f["MCMC_OUTPUT"].keys())
    # if len(chains) and betaID !=0:
    #     print("This file doesn't have chains hotter than Beta=1!")
    #     return None, None, None
    chains_N = len(chains)
    data = f["MCMC_OUTPUT"]["CHAIN {}".format(chainIDs[0])]
    status = f["MCMC_OUTPUT/STATUS"]["CHAIN {}".format(chainIDs[0])]
    for x in chainIDs[1:]:
        data = np.insert(data,-1, f["MCMC_OUTPUT"]["CHAIN {}".format(x)],axis=0)
        status = np.insert(status,-1, f["MCMC_OUTPUT/STATUS"]["CHAIN {}".format(x)],axis=0)
    model_status = []
    if "MCMC_OUTPUT/MODEL_STATUS" in f.keys():
        model_status = f["MCMC_OUTPUT/MODEL_STATUS"]["CHAIN {}".format(chainIDs[0])]
        for x in chainIDs[1:]:
            model_status = np.insert(model_status,-1, f["MCMC_OUTPUT/MODEL_STATUS"]["CHAIN {}".format(x)],axis=0)
    return data, status,model_status

def MCMC_unpack_file(filename,betaID=0, trim=None, thin=None):

    f = h5py.File(filename,'r')
    betas = np.array(f["MCMC_METADATA"]["CHAIN BETAS"])
    betaSchedule = np.flip(np.unique(np.array(f["MCMC_METADATA"]["CHAIN BETAS"])))
    ensembleSize = int(len(betaSchedule))
    ensembleN = int(len(betas)/ensembleSize)

    if betaID > ensembleN:
        print("Supplied a betaID larger than the number of ensembles!")
        return None, None, None

    chainIDs = np.arange(chainIndex(0,betaID,ensembleN) , chainIndex(ensembleN,betaID,ensembleN))

    chains = list(f["MCMC_OUTPUT"].keys())

    trim_local=0
    thin_local=1
    #if trim is None :
    #    trim_local = f["MCMC_METADATA"]["SUGGESTED TRIM LENGTHS"][chainIDs[0]]
    #if ac is None:
    #    aclist = []
    #    for x in np.arange(len(f["MCMC_METADATA"]["AC VALUES"])):
    #        aclist.append(np.amax(f["MCMC_METADATA"]["AC VALUES"][x][:]))
    #    ac_local = np.mean(aclist)


    # if len(chains) and betaID !=0:
    #     print("This file doesn't have chains hotter than Beta=1!")
    #     return None, None, None
    chains_N = len(chains)
    if trim is None  and betaID ==0:
        trim_local = f["MCMC_METADATA"]["SUGGESTED TRIM LENGTHS"][chainIDs[0]]
    if thin is None and betaID ==0:
        thin_local = np.amax(f["MCMC_METADATA"]["AC VALUES"][chainIDs[0]][:])
    if thin_local == 0 :
        thin_local=1
    data = f["MCMC_OUTPUT"]["CHAIN {}".format(chainIDs[0])][trim_local::thin_local]
    for x in chainIDs[1:]:
        if trim is None and betaID ==0:
            trim_local = f["MCMC_METADATA"]["SUGGESTED TRIM LENGTHS"][x]
        if thin is None and betaID ==0:
            thin_local = np.amax(f["MCMC_METADATA"]["AC VALUES"][x][:])
        if thin_local == 0:
            thin_local=1
        data = np.insert(data,-1, f["MCMC_OUTPUT"]["CHAIN {}".format(x)][trim_local::thin_local],axis=0)
    return data


def RJ_corner(data,status,figsize=None,marginal_bins=20,cov_bins=20,show_quantiles=False,titles=None,marginal_color='black',cov_color='gray',alpha=1):
    data_shape = np.shape(data)
    dim = data_shape[1]

    fig, axes = plt.subplots(nrows=dim,ncols=dim,figsize=figsize,sharex='col')

    for x in np.arange(dim):
        for y in np.arange(x+1):
            axes[x,y].grid(False)
            if y == x:
                #marginalized
                mask = status[:,x] == 1
                if(np.sum(mask) !=0):
                    axes[x,y].hist(data[mask,x],bins=marginal_bins,density=True,edgecolor=marginal_color,histtype='step', color=marginal_color,alpha=alpha)
                axes[x,y].get_yaxis().set_ticks([])
                if titles is not None:
                    if show_quantiles and np.sum(mask) != 0:
                        fifty = np.quantile(data[mask,x],.50)
                        upper = np.quantile(data[mask,x],.84)#1 sigma, for symmetric, gaussian like distribution
                        lower = np.quantile(data[mask,x],.16)#1 sigma, for symmetric, gaussian like distribution
                        axes[x,y].set_title(r"${0:1.2e}^{{+ {1:1.2e} }}_{{- {2:1.2e} }}$".format(fifty,upper-fifty,fifty-lower))
                    else:
                        axes[x,y].set_title(str(titles[x]))
                else:
                    if show_quantiles and np.sum(mask) != 0:
                        fifty = np.quantile(data[mask,x],.50)
                        upper = np.quantile(data[mask,x],.84)#1 sigma, for symmetric, gaussian like distribution
                        lower = np.quantile(data[mask,x],.16)#1 sigma, for symmetric, gaussian like distribution
                        axes[x,y].set_title(r"${0:1.2e}^{{+ {1:1.2e} }}_{{- {2:1.2e} }}$".format(fifty,upper-fifty,fifty-lower))
            else:
                #covariance
                mask = (status[:,x] == 1) & (status[:,y] == 1)
                if(np.sum(mask) !=0):
                    axes[x,y].hexbin(data[mask,y],data[mask,x],bins=cov_bins,mincnt=1,cmap=cov_color,alpha=alpha)
                axes[x,y].get_shared_y_axes().join(axes[0,y])
                if y != 0 :
                    axes[x,y].get_yaxis().set_ticks([])
                else:
                    if titles is not None:
                        axes[x,y].set_ylabel(titles[x])
    for x in np.arange(dim-1):
        for y in np.arange(x+1,dim):
            axes[x,y].axis('off')
    if titles is not None:
        for x in np.arange(dim):
            axes[dim-1,x].set_xlabel(titles[x])
    fig.subplots_adjust(hspace=0.03,wspace=0.03)
    return fig
