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

std_percent_cut = .3
max_percent_cut = .3

#write down n channel numbers you want to read!!!!!!!!!!!!!!
nChannels = 16
Channels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

#dataset parameters
trigger_position = 0.25
sample_interval = 2000 #us
gain = 20
ADCconv = (3*(2./4096))/(gain*1200.*10.*2) #1000 factor difference with the old digitizer
nSamples = 1024

fd = open("dataset1.txt", "r")
directory = fd.readline().strip()
noise = fd.readline().strip()

print(directory+"/"+noise+"/")
files = glob.glob(directory+"/"+noise+"/*.npy")
#files = glob.glob("/data/Run34_NSC/3420210713_gain20_SiHV_ACD_70V_noise/*.txt")

saveF = ROOT.TFile(noise+"_noise.root","RECREATE")

nEvents = len(files)
ChannelEvents = [0 for i in range(len(Channels))]
print("There are "+str(len(Channels))+" Channels. Total number of Noise Events: "+str(nEvents))
print("Pulse length is: "+str(nSamples))
print("channels: "+str(Channels))

std_cut, max_cut,stdN, maxN,stdL,maxL = dict(), dict(), dict(), dict(), dict(), dict()

for i in range(nChannels):
    stdN[i] = ROOT.TH1D("std_"+str(i), " ; ", 1000, 0, .5e-3)
    maxN[i]= ROOT.TH1D("max_"+str(i), " ; ", 1000, 0, 1e-3)

    stdL[i] = list()
    maxL[i] = list()

time = np.arange(nSamples,dtype=np.float64)
time2 = time*sample_interval
sf = (1/(sample_interval*1e-9))

s = (nChannels, nEvents,nSamples)
p = np.zeros(nSamples)
noise = np.zeros(s)

for nN,filen in enumerate(files):
    data = np.loadtxt(filen)
    for i in range(nChannels):
        p = data[:, i] #i+1
        meanA = p.mean()
        p -= meanA
        p *= ADCconv
        noise[i][nN][:] = p[:]
        stdN[i].Fill(p.std())
        maxN[i].Fill(p.max())
        stdL[i].append(p.std())
        maxL[i].append(p.max())

print("noise: "+str(nEvents))

'''
#manual cuts
stdAmincut = 5e-6
stdAmaxcut = 0.037e-3
maxAmaxcut = 0.08e-3

stdBmincut = 7e-6
stdBmaxcut = 0.03e-3
maxBmaxcut = 0.014e-3

stdCmincut = 2e-6
stdCmaxcut = 0.036e-3
maxCmaxcut = 0.1e-3

stdDmincut = 5e-6
stdDmaxcut = .03e-3
maxDmaxcut = 0.07e-3
'''
nGood, avgNoise, avgPSD, passA,n,PSD = dict(), dict(), dict(), dict(), dict(), dict()
for i in range(nChannels):
    stdL[i].sort()
    maxL[i].sort()
    std_cut[i] = stdL[i][int(std_percent_cut * len(stdL[i]))]
    max_cut[i] = maxL[i][int(max_percent_cut * len(maxL[i]))]
    print("STD_"+str(i)+":" + str(std_cut[i]))
    print("MAX_"+str(i)+":" + str(max_cut[i]))
    nGood[i] = 0

    avgNoise[i] = np.zeros(nSamples)
    avgPSD[i] = np.zeros(int(nSamples / 2))

freq1, throw1 = qp.calc_psd(avgNoise[0],sf, folded_over=False)
freq2 = freq1[0:nSamples-1]

for nE in range(nEvents):
    if nE%100 == 0:
        print(nE)
    for i in range(nChannels):
        passA[i] = False

        p = noise[i][nE]
        if 0.0 not in p and p.std() < std_cut[i] and p.max() < max_cut[i]: passA[i] = True

        if passA[i]:
           chAThrowaway,chNoise = qp.calc_psd(p,sf, folded_over=False)
           chAThrow,chPSD = qp.calc_psd(p,sf, folded_over=True)
           avgPSD[i] += chPSD[0:int(nSamples/2)]
           avgNoise[i] += chNoise
           nGood[i] += 1

for i in range(nChannels):
    print(float(nGood[i]))
    avgNoise[i] = avgNoise[i] / float(nGood[i])
    avgPSD[i] = avgPSD[i] / float(nGood[i])
    avgPSD[i] = np.sqrt(avgPSD[i])

    n[i] = ROOT.TGraph(nSamples,freq1,avgNoise[i])
    n[i].SetName("noise"+str(i)+"_template")
    n[i].Write()

    PSD[i] = ROOT.TGraph(int(nSamples/2),freq2,avgPSD[i])
    PSD[i].SetName("PSD"+str(i))
    PSD[i].Write()

saveF.Write()
print(nEvents)






