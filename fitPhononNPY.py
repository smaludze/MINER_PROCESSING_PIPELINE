#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import qetpy as qp
from qetpy import Noise
from qetpy.sim import TESnoise
from qetpy.plotting import compare_noise, plot_noise_sim
from pprint import pprint
import math
import os
import glob
import sys
import linecache
from array import array
import ROOT
from scipy.signal import butter, lfilter
from scipy import interpolate, optimize
from nptdms import TdmsFile
from ROOT import kRed, kBlue
from time import process_time
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

def butter_bandpass(lowcut, fs, order=5):
    low = lowcut / (0.5 * fs)
    b, a = butter(order, low, btype='lowpass')
    return b, a


def butter_bandpass_filter(data, lowcut, fs, order=5):
    b, a = butter_bandpass(lowcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def fileList(directory):
    flist1 = [glob.glob(directory+'sapphire_ch{}_{}.txt'.format(str(j).zfill(2),str(i).zfill(4))) for i in range(LastFileN-FirstFileN+1)
              for j in Channels]
    flist = [flist1[i][0] for i in range(len(flist1))]
    return flist

#things to change:
# template string for each channel in the beginning !!!
# dictionary with templates for each channel
#read events from the first 5 files
global FirstFileN
FirstFileN = 0
global LastFileNumber
LastFileN = 4
start_time = process_time()
#write down n channel numbers you want to read!!!!!!!!!!!!!!
nDetectors = 4
nChannels = 16
Channels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#include only the directory with "/" at the end. File naming certain format is assumed: "ASCII_set1_ch01_.txt".
#if file name has different format it needs to be changed in the 'fileList' function above.
#files = fileList("/data/Run34_NSC/3420210708_gain20_SiHV_ACD_125V_pulse/")
CutOff = 10000000
#dataset parameters
trigger_position = 0.25
sample_interval = 2000 #us
gain = 20
ADCconv = (3*(2./4096))/(gain*1200.*10.*2) #1000 factor difference with the old digitizer
nSamples = 1024

#dataset directory (noise and pulse)
fd = open("dataset1.txt", "r")
directory = fd.readline().strip() 
noise = fd.readline().strip()
pulse = fd.readline().strip()

print(directory+"/"+pulse+"/")
files = glob.glob(directory+"/"+pulse+"/*/*/*.npy")
#files = glob.glob("/data/Run36_NSC/3620210929_gain50_trgChC_Sapphire17_Shielded_ptOFF*pulse/*.tdms")

#tF = ROOT.TFile(pulse+"_pulse"+".root")
tF = ROOT.TFile("39211203_4mm2_trgBCD_prom10000_calib_extracted_pulses_pulse.root")
tF1 = ROOT.TFile("39211207_sapphire7_4DG20_Digi64Ch_noTR_BIN_1_extracted_pulses_pulse.root")

nF = ROOT.TFile(noise+"_noise.root")
#nF = ROOT.TFile("39211203_901_4mm2_903_trgBCD_prom10000_DRU_extracted_noise_noise.root")
saveF = ROOT.TFile(pulse+"_pulse_fitPhonon"+".root","RECREATE")

#pulse templates for each channel (set input for each channel!!!)
templates = {}
templates[0] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[1] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[2] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[3] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[4] = np.array(np.frombuffer(tF1.Get("pulseCh0_2").GetY(), dtype=np.double))
templates[5] = np.array(np.frombuffer(tF1.Get("pulseCh1_2").GetY(), dtype=np.double))
templates[6] = np.array(np.frombuffer(tF1.Get("pulseCh2_2").GetY(), dtype=np.double))
templates[7] = np.array(np.frombuffer(tF1.Get("pulseCh3_3").GetY(), dtype=np.double))
templates[8] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[9] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[10] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[11] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[12] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[13] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[14] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[15] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))


#noise templates for each channel (set input for each channels!!!)
NoiseTemplates = {}
NoiseTemplates[0] = np.array(np.frombuffer(nF.Get("noise0_template").GetY(), dtype=np.double))
NoiseTemplates[1] = np.array(np.frombuffer(nF.Get("noise1_template").GetY(), dtype=np.double))
NoiseTemplates[2] = np.array(np.frombuffer(nF.Get("noise2_template").GetY(), dtype=np.double))
NoiseTemplates[3] = np.array(np.frombuffer(nF.Get("noise3_template").GetY(), dtype=np.double))
NoiseTemplates[4] = np.array(np.frombuffer(nF.Get("noise4_template").GetY(), dtype=np.double))
NoiseTemplates[5] = np.array(np.frombuffer(nF.Get("noise5_template").GetY(), dtype=np.double))
NoiseTemplates[6] = np.array(np.frombuffer(nF.Get("noise6_template").GetY(), dtype=np.double))
NoiseTemplates[7] = np.array(np.frombuffer(nF.Get("noise7_template").GetY(), dtype=np.double))
NoiseTemplates[8] = np.array(np.frombuffer(nF.Get("noise8_template").GetY(), dtype=np.double))
NoiseTemplates[9] = np.array(np.frombuffer(nF.Get("noise9_template").GetY(), dtype=np.double))
NoiseTemplates[10] = np.array(np.frombuffer(nF.Get("noise10_template").GetY(), dtype=np.double))
NoiseTemplates[11] = np.array(np.frombuffer(nF.Get("noise11_template").GetY(), dtype=np.double))
NoiseTemplates[12] = np.array(np.frombuffer(nF.Get("noise12_template").GetY(), dtype=np.double))
NoiseTemplates[13] = np.array(np.frombuffer(nF.Get("noise13_template").GetY(), dtype=np.double))
NoiseTemplates[14] = np.array(np.frombuffer(nF.Get("noise14_template").GetY(), dtype=np.double))
NoiseTemplates[15] = np.array(np.frombuffer(nF.Get("noise15_template").GetY(), dtype=np.double))


ChannelEvents = [0 for i in range(len(Channels))]
print("Pulse length is: "+str(nSamples))
print("channels: "+str(Channels))

time = np.arange(nSamples,dtype=np.float64)
sl_sec = sample_interval*1e-9
time2 = time*sample_interval
triggerTime = trigger_position*time2[nSamples-1]
sf = (1/(sample_interval*1e-9))

t = {}
nEvents = len(files)
eventsInRoot = 100
print(int(math.floor(nEvents/eventsInRoot)))
for i in range(int(math.floor(nEvents/eventsInRoot))):
	t[i] = ROOT.TTree('data'+str(i), 'main tree'+str(i))

ampOF,chi2OF,chi2OFnop,t0OF,rtft20,mean,IntT,stdpre,stdtail,minn,maxbin,maxtail,riseT,delayT,riseOF = \
dict(),dict(),dict(),dict(),dict(), dict(),dict(),dict(),dict(),dict(),dict(),dict(),dict(),dict(),dict()

#eventTT = np.zeros(nEvents)
#eventNN = np.zeros(nEvents)

for i in range(nChannels):
    ampOF[i] = np.zeros(1, dtype=np.float64)
    chi2OF[i] = np.zeros(1, dtype=np.float64)
    chi2OFnop[i] = np.zeros(1, dtype=np.float64)
    t0OF[i] = np.zeros(1, dtype=np.float64)
    rtft20[i] = np.zeros(1, dtype=np.float64)
    mean[i] = np.zeros(1, dtype=np.float64)
    IntT[i] = np.zeros(1, dtype=np.float64)
    stdpre[i] = np.zeros(1, dtype=np.float64)
    stdtail[i] = np.zeros(1, dtype=np.float64)
    minn[i] = np.zeros(1, dtype=np.float64)
    maxbin[i] = np.zeros(1, dtype=np.float64)
    maxtail[i] = np.zeros(1, dtype=np.float64)
    riseT[i] = np.zeros(1, dtype=np.float64)
    delayT[i] = np.zeros(1, dtype=np.float64)
    riseOF[i] = np.zeros(1, dtype=np.float64)

    t[0].Branch('ampOF'+str(i), ampOF[i], 'ampOF'+str(i)+'/D')
    t[0].Branch('chi2OF' + str(i), chi2OF[i], 'chi2OF' + str(i) + '/D')
    t[0].Branch('t0OF' + str(i), t0OF[i], 't0OF' + str(i) + '/D')
    t[0].Branch('rtft20' + str(i), rtft20[i], 'rtft20' + str(i) + '/D')
    t[0].Branch('mean' + str(i), mean[i], 'mean' + str(i) + '/D')
    t[0].Branch('IntT' + str(i), IntT[i], 'IntT' + str(i) + '/D')
    t[0].Branch('stdpre' + str(i), stdpre[i], 'stdpre' + str(i) + '/D')
    t[0].Branch('stdtail' + str(i), stdtail[i], 'stdtail' + str(i) + '/D')
    t[0].Branch('minn' + str(i), minn[i], 'minn' + str(i) + '/D')
    t[0].Branch('maxbin' + str(i), maxbin[i], 'maxbin' + str(i) + '/D')
    t[0].Branch('maxtail' + str(i), maxtail[i], 'maxtail' + str(i) + '/D')
    t[0].Branch('riseT' + str(i), riseT[i], 'riseT' + str(i) + '/D')
    t[0].Branch('delayT' + str(i), delayT[i], 'delayT' + str(i) + '/D')
    t[0].Branch('riseOF' + str(i), riseOF[i], 'riseOF' + str(i) + '/D')

partXOF,partYOF,partROF,partRXYOF,partthetaOF,partX,partY,partR,partRXY,parttheta,delayXOF,delayYOF,delayX, \
delayY,delayR = dict(),dict(),dict(),dict(),dict(), dict(),dict(),dict(),dict(),dict(),dict(),dict(),dict(),dict(),dict()

for i in range(nDetectors):
    partXOF[i] = np.zeros(1, dtype=np.float64)
    partYOF[i] = np.zeros(1, dtype=np.float64)
    partROF[i] = np.zeros(1, dtype=np.float64)
    partRXYOF[i] = np.zeros(1, dtype=np.float64)
    partthetaOF[i] = np.zeros(1, dtype=np.float64)
    partX[i] = np.zeros(1, dtype=np.float64)
    partY[i] = np.zeros(1, dtype=np.float64)
    partR[i] = np.zeros(1, dtype=np.float64)
    partRXY[i] = np.zeros(1, dtype=np.float64)
    parttheta[i] = np.zeros(1, dtype=np.float64)
    delayXOF[i] = np.zeros(1, dtype=np.float64)
    delayYOF[i] = np.zeros(1, dtype=np.float64)
    delayX[i] = np.zeros(1, dtype=np.float64)
    delayY[i] = np.zeros(1, dtype=np.float64)
    delayR[i] = np.zeros(1, dtype=np.float64)

    t[0].Branch('partXOF'+str(i), partXOF[i], 'partXOF'+str(i)+'/D')
    t[0].Branch('partYOF'+str(i), partYOF[i], 'partYOF'+str(i)+'/D')
    t[0].Branch('partROF'+str(i), partROF[i], 'partROF'+str(i)+'/D')
    t[0].Branch('partRXYOF'+str(i), partRXYOF[i], 'partRXYOF'+str(i)+'/D')
    t[0].Branch('partthetaOF'+str(i), partthetaOF[i], 'partthetaOF'+str(i)+'/D')
    t[0].Branch('partY'+str(i), partY[i], 'partY'+str(i)+'/D')
    t[0].Branch('partX'+str(i), partX[i], 'partX'+str(i)+'/D')
    t[0].Branch('partR'+str(i), partR[i], 'partR'+str(i)+'/D')
    t[0].Branch('parttheta'+str(i), parttheta[i], 'parttheta'+str(i)+'/D')
    t[0].Branch('partRXY'+str(i), partRXY[i], 'partRXY'+str(i)+'/D')
    t[0].Branch('delayXOF'+str(i), delayXOF[i], 'delayXOF'+str(i)+'/D')
    t[0].Branch('delayYOF'+str(i), delayYOF[i], 'delayYOF'+str(i)+'/D')
    t[0].Branch('delayX'+str(i), delayX[i], 'delayX'+str(i)+'/D')
    t[0].Branch('delayY'+str(i), delayY[i], 'delayY'+str(i)+'/D')
    t[0].Branch('delayR'+str(i), delayR[i], 'delayR'+str(i)+'/D')

eventTime = np.zeros(1, dtype=np.float64)
eventNum = np.zeros(1, dtype=np.int64)
trigger = np.zeros(1, dtype=np.int64)
sliceN = np.zeros(1, dtype=np.int64)
fileN = np.zeros(1, dtype=np.int64)
pulseN = np.zeros(1, dtype=np.int64)

t[0].Branch( 'eventTime', eventTime, 'eventTime/D')
t[0].Branch( 'eventNum', eventNum, 'eventNum/I')
t[0].Branch( 'trigger',trigger,'trigger/I')
t[0].Branch( 'sliceN',sliceN,'sliceN/I')
t[0].Branch( 'fileN',fileN,'fileN/I')
t[0].Branch( 'pulseN',pulseN,'pulseN/I')

preBins = int(.2*nSamples)
postBins = int(.2*nSamples)

p = np.zeros(nSamples)
print(len(t))
for nE, filen in enumerate(files):
    #if nE<1000000:
    	#continue
    #if nE == CutOff:
    	#break
    if nE % eventsInRoot == 0 and int(nE/eventsInRoot) > 0 and  False:
        saveF.Write() 
        saveF.Close()
        saveF = ROOT.TFile(pulse+"_pulse_fitPhonon_"+str(int(nE/eventsInRoot))+".root","RECREATE")
        #t[int(nE/eventsInRoot)] = ROOT.TTree('data'+str(int(nE/eventsInRoot)), 'main tree'+str(int(nE/eventsInRoot)))
        #initialize new tree with branches --------------------------------------------------------
        
        for i in range(nChannels):
            '''
            ampOF[i] = np.zeros(1, dtype=np.float64)
            chi2OF[i] = np.zeros(1, dtype=np.float64)
            chi2OFnop[i] = np.zeros(1, dtype=np.float64)
            t0OF[i] = np.zeros(1, dtype=np.float64)
            rtft20[i] = np.zeros(1, dtype=np.float64)
            mean[i] = np.zeros(1, dtype=np.float64)
            IntT[i] = np.zeros(1, dtype=np.float64)
            stdpre[i] = np.zeros(1, dtype=np.float64)
            stdtail[i] = np.zeros(1, dtype=np.float64)
            minn[i] = np.zeros(1, dtype=np.float64)
            maxbin[i] = np.zeros(1, dtype=np.float64)
            maxtail[i] = np.zeros(1, dtype=np.float64)
            riseT[i] = np.zeros(1, dtype=np.float64)
            delayT[i] = np.zeros(1, dtype=np.float64)
            riseOF[i] = np.zeros(1, dtype=np.float64)
            '''
            t[1].Branch('ampOF'+str(i), ampOF[i], 'ampOF'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('ampOF'+str(i), ampOF[i], 'ampOF'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('chi2OF' + str(i), chi2OF[i], 'chi2OF' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('t0OF' + str(i), t0OF[i], 't0OF' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('rtft20' + str(i), rtft20[i], 'rtft20' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('mean' + str(i), mean[i], 'mean' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('IntT' + str(i), IntT[i], 'IntT' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('stdpre' + str(i), stdpre[i], 'stdpre' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('stdtail' + str(i), stdtail[i], 'stdtail' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('minn' + str(i), minn[i], 'minn' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('maxbin' + str(i), maxbin[i], 'maxbin' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('maxtail' + str(i), maxtail[i], 'maxtail' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('riseT' + str(i), riseT[i], 'riseT' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('delayT' + str(i), delayT[i], 'delayT' + str(i) + '/D')
            t[int(nE/eventsInRoot)].Branch('riseOF' + str(i), riseOF[i], 'riseOF' + str(i) + '/D')
        for i in range(nDetectors):
            '''
            partXOF[i] = np.zeros(1, dtype=np.float64)
            partYOF[i] = np.zeros(1, dtype=np.float64)
            partROF[i] = np.zeros(1, dtype=np.float64)
            partRXYOF[i] = np.zeros(1, dtype=np.float64)
            partthetaOF[i] = np.zeros(1, dtype=np.float64)
            partX[i] = np.zeros(1, dtype=np.float64)
            partY[i] = np.zeros(1, dtype=np.float64)
            partR[i] = np.zeros(1, dtype=np.float64)
            partRXY[i] = np.zeros(1, dtype=np.float64)
            parttheta[i] = np.zeros(1, dtype=np.float64)
            delayXOF[i] = np.zeros(1, dtype=np.float64)
            delayYOF[i] = np.zeros(1, dtype=np.float64)
            delayX[i] = np.zeros(1, dtype=np.float64)
            delayY[i] = np.zeros(1, dtype=np.float64)
            delayR[i] = np.zeros(1, dtype=np.float64)
            '''
            t[int(nE/eventsInRoot)].Branch('partXOF'+str(i), partXOF[i], 'partXOF'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('partYOF'+str(i), partYOF[i], 'partYOF'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('partROF'+str(i), partROF[i], 'partROF'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('partRXYOF'+str(i), partRXYOF[i], 'partRXYOF'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('partthetaOF'+str(i), partthetaOF[i], 'partthetaOF'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('partY'+str(i), partY[i], 'partY'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('partX'+str(i), partX[i], 'partX'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('partR'+str(i), partR[i], 'partR'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('parttheta'+str(i), parttheta[i], 'parttheta'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('partRXY'+str(i), partRXY[i], 'partRXY'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('delayXOF'+str(i), delayXOF[i], 'delayXOF'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('delayYOF'+str(i), delayYOF[i], 'delayYOF'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('delayX'+str(i), delayX[i], 'delayX'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('delayY'+str(i), delayY[i], 'delayY'+str(i)+'/D')
            t[int(nE/eventsInRoot)].Branch('delayR'+str(i), delayR[i], 'delayR'+str(i)+'/D')
        eventTime = np.zeros(1, dtype=np.float64)
        eventNum = np.zeros(1, dtype=np.int64)
        trigger = np.zeros(1, dtype=np.int64)
        sliceN = np.zeros(1, dtype=np.int64)
        fileN = np.zeros(1, dtype=np.int64)
        pulseN = np.zeros(1, dtype=np.int64)

        t[int(nE/eventsInRoot)].Branch('eventTime', eventTime, 'eventTime/D')
        t[int(nE/eventsInRoot)].Branch('eventNum', eventNum, 'eventNum/I')
        t[int(nE/eventsInRoot)].Branch('trigger',trigger,'trigger/I')
        t[int(nE/eventsInRoot)].Branch('sliceN',sliceN,'sliceN/I')
        t[int(nE/eventsInRoot)].Branch('fileN',fileN,'fileN/I')
        t[int(nE/eventsInRoot)].Branch('pulseN',pulseN,'pulseN/I')
        
        #-----------------------------------------------------------------------------------------
    #read one pulse------------------------------------------------------------------------------- 
    data = np.loadtxt(filen)
    if (nE+1)%100 == 0:
        print("Events processed: "+str(nE+1))
    #if nE == 1000:
        #break
    for i in range(nChannels):
        p = data[:, i] #i+1
        p*=ADCconv
        p = p-p[:preBins].mean()

        stdpre[i][0] = p[0:preBins].std()
        stdtail[i][0] = p[nSamples - postBins:nSamples].std()
        minn[i][0] = p.min()
        maxtail[i][0] = p[nSamples - postBins:nSamples].max()

        OF = qp.OptimumFilter(p, templates[i], NoiseTemplates[i], sf)  # initialize the OptimumFilter class
        amp, t0, chi2 = OF.ofamp_withdelay()

        chi2OFnop[i][0] = OF.chi2_nopulse()
        ampOF[i][0] = amp
        Integ = amp * templates[i]
        IntT[i][0] = Integ.sum()
        chi2OF[i][0] = chi2
        t0OF[i][0] = t0
        pfilt = butter_bandpass_filter(p, 50000, sf, order=2)
        maxbin = np.argmax(pfilt)
        maxbin1 = np.argmax(pfilt)
        pmax = pfilt.max()
        lowerbin = 0
        upperbin = 0
        for j in range(nSamples):
            if pfilt[maxbin - j] < .2 * pmax:
                lowerbin = maxbin - j
                break

        for j in range(nSamples):
            if maxbin + j > nSamples - 1:
                upperbin = nSamples - 1
                break
            elif pfilt[maxbin + j] < .2 * pmax:
                upperbinA = maxbin + j
                break

        rtft20[i][0] = pfilt[lowerbin:upperbin].sum()
        maxbinT = np.argmax(pfilt[200:(nSamples - 200)])
        maxbinT += 200
        pmaxT = pfilt.max()
        riseT[i][0] = 0
        for j in range(int(nSamples / 4)):
            if pfilt[maxbinT - j] < .2 * pmaxT:
                m = (pfilt[maxbinT - j + 1]) / sl_sec
                riseT[i][0] = (.2 * pmaxT) / m + (maxbinT - j) * sl_sec
                break
        # Calculate 10% and 40% rise timestamps to find risetime of the pulse
        rise10 = 0
        rise20 = 0
        rise40 = 0
        riseOF[i][0] = 0
        for j in range(int(nSamples / 4)):
            if pfilt[maxbinT - j] < .4 * pmaxT and maxbinT - j > 0:
                y = np.array(
                    [pfilt[maxbinT - j - 1], pfilt[maxbinT - j], pfilt[maxbinT - j + 1], pfilt[maxbinT - j + 2]])
                x = np.array(
                    [time2[maxbinT - j - 1], time2[maxbinT - j], time2[maxbinT - j + 1], time2[maxbinT - j + 2]])
                # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA)).roots()
                f = interpolate.interp1d(x, y, kind='linear')
                g = lambda x: f(x) - .4 * pmaxT
                # plt.plot(time2,pfilt,'r.',x,f(x),'g--')
                try:
                    root = optimize.newton(g, time2[maxbinT - j])
                    rise40 = root
                except (ValueError, RuntimeError):
                    rise40 = 0
                break
        for j in range(int(nSamples / 4)):
            if pfilt[maxbinT - j] < .2 * pmaxT and maxbinT - j > 0:
                y = np.array(
                    [pfilt[maxbinT - j - 1], pfilt[maxbinT - j], pfilt[maxbinT - j + 1], pfilt[maxbinT - j + 2]])
                x = np.array(
                    [time2[maxbinT - j - 1], time2[maxbinT - j], time2[maxbinT - j + 1], time2[maxbinT - j + 2]])
                # spl = interpolate.InterpolatedUnivariateSpline(x,y-(.4*pmaxTA)).roots()
                f = interpolate.interp1d(x, y, kind='linear')
                g = lambda x: f(x) - .2 * pmaxT
                # plt.plot(time2,pfilt,'r.')
                try:
                    root = optimize.newton(g, time2[maxbinT - j])
                    rise20 = root
                    delayT[i][0] = (root - triggerTime) * 1e-3
                except (ValueError, RuntimeError):
                    rise20 = 0
                    delayT[i][0] = 0
                break

        for j in range(int(nSamples / 4)):
            if pfilt[maxbinT - j] < .1 * pmaxT and maxbinT - j > 0:
                y = np.array(
                    [pfilt[maxbinT - j - 1], pfilt[maxbinT - j], pfilt[maxbinT - j + 1], pfilt[maxbinT - j + 2]])
                x = np.array(
                    [time2[maxbinT - j - 1], time2[maxbinT - j], time2[maxbinT - j + 1], time2[maxbinT - j + 2]])
                f = interpolate.interp1d(x, y, kind='linear')
                g = lambda x: f(x) - .1 * pmaxT
                try:
                    root = optimize.newton(g, time2[maxbinT - j])
                    rise10 = root
                except (ValueError, RuntimeError):
                    rise10 = 0
                break
        if rise10 != 0 and rise40 != 0:
            riseOF[i][0] = (rise40 - rise10) * 1e-3

    for z in range(nDetectors):
        k = 4 * z
        delayXOF[z][0] = -1 * (.866 * riseOF[k + 3][0] - .866 * riseOF[k + 1][0]) * 1e-3
        delayYOF[z][0] = -1 * (.5 * riseOF[k + 3][0] + .5 * riseOF[k + 1][0] - riseOF[k + 2][0]) * 1e-3
        '''
        delayX[0] = -1 * (.866 * delayTD[0] - .866 * delayTB[0])
        delayY[0] = -1 * (.5 * delayTD[0] + .5 * delayTB[0] - delayTC[0])
        ampOFarray = np.array([ampOFB[0], ampOFC[0], ampOFD[0]])
        delayOFarray = np.array([delayTB[0], delayTC[0], delayTD[0]])
        maxampOFindex = np.argmax(ampOFarray)
        # delayR[0] = delayOFarray[maxampOFindex] - delayTA[0]
        '''
        partXOF[z][0] = (.866 * ampOF[k + 3][0] - .866 * ampOF[k + 1][0]) / (
                    ampOF[k + 3][0] + ampOF[k + 1][0] + ampOF[k + 2][0])
        partYOF[z][0] = (.5 * ampOF[k + 3][0] + .5 * ampOF[k + 1][0] - ampOF[k + 2][0]) / (
                    ampOF[k + 3][0] + ampOF[k + 2][0] + ampOF[k + 1][0])
        '''
        partROF[0] = (ampOFA[0]) / (ampOFA[0] + ampOFB[0] + ampOFC[0] + ampOFD[0])
        partRXYOF[0] = np.sqrt(partYOF[0] * partYOF[0] + partXOF[0] * partXOF[0])
        partthetaOF[0] = np.arctan2(partYOF[0], partXOF[0])
        partX[0] = (.866 * rtft20D[0] - .866 * rtft20B[0]) / (rtft20D[0] + rtft20B[0] + rtft20C[0])
        partY[0] = (.5 * rtft20D[0] + .5 * rtft20B[0] - rtft20C[0]) / (rtft20D[0] + rtft20B[0] + rtft20C[0])
        partR[0] = (rtft20A[0]) / (rtft20A[0] + rtft20B[0] + rtft20C[0] + rtft20D[0])
        partRXY[0] = np.sqrt(partY[0] * partY[0] + partX[0] * partX[0])
        parttheta[0] = np.arctan2(partY[0], partX[0])
        '''

    fileN[0] = int(filen.split('/')[-3].split('_')[1])
    sliceN[0] = int(filen.split('/')[-2].split('_')[1])
    pulseNpy = filen.split('/')[-1]
    pulseN[0] = int(pulseNpy[:len(pulseNpy)-4].split('_')[1])
    #eventTime[0] = eventTT[nE]
    #eventNum[0] = eventNN[nE]

    t[int(nE/eventsInRoot)].Fill()

saveF.Write()
saveF.Close()
stop_time = process_time()
print("Total time needed to process: "+str(round((stop_time-start_time)/60))+" min")
print("Total number of events processed: "+str(len(files)))
