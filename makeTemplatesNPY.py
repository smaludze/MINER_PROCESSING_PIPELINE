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
from time import process_time

start_time = process_time()
#write down n channel numbers you want to read!!!!!!!!!!!!!!
nChannels = 12
Channels = [0,1,2,3,4,5,6,7,8,9,10,11]

#dataset parameters
trigger_position = 0.25
sample_interval = 2000 #us
gain = 50
ADCconv = (3*(2./4096))/(gain*1200.*10.*2) #1000 factor difference with the old digitizer
nSamples = 1024

#How many events do you want to read in order to construct template?
CutOff = 20000

fd = open("dataset1.txt", "r")
directory = fd.readline().strip() 
fd.readline()
pulse = fd.readline().strip()
print(directory+"/"+pulse+"/")

#files = glob.glob("/data/Run34_NSC/3420210709_gain50_SiHV_ACD_0V_pulse/*.txt")
files = glob.glob(directory+"/"+pulse+"/*/*/*.npy")
saveF = ROOT.TFile(pulse+"_pulse"+".root","RECREATE")

nEvents = len(files)
CutOff = min(CutOff,nEvents)
ChannelEvents = [0 for i in range(len(Channels))]
print("There are "+str(len(Channels))+" Channels. Total number of Events: "+str(nEvents))
print("Total number of Events to read: "+str(CutOff))
print("Pulse length is: "+str(nSamples))
print("channels: "+str(Channels))

time = np.arange(nSamples,dtype=np.float64)
time2 = time*sample_interval
sf = (1/(sample_interval*1e-9))

p = np.zeros(nSamples)
Pulse = np.zeros((nChannels, CutOff, nSamples))

postBin = int(0.2*nSamples)
preBin = int(0.2*nSamples)

stdA, maxA, maxTail, minA, stdTail = dict(), dict(), dict(), dict(), dict()
for i in range(nChannels):
    stdA[i] = ROOT.TH1D("std_pre"+str(i), " ; ", 1000, 0, .1e-5)
    maxA[i] = ROOT.TH1D("max"+str(i), " ; ", 1000, 0, 3e-5)
    maxTail[i] = ROOT.TH1D("max_tail"+str(i), " ; ", 1000, 0, 1e-6)
    minA[i] = ROOT.TH1D("min"+str(i), " ; ", 1000, -.5e-6, 0)
    stdTail[i] = ROOT.TH1D("std_tail"+str(i), " ; ", 1000, 0, .1e-5)

for nE,filen in enumerate(files):
    data = np.loadtxt(filen)
    for i in range(nChannels):
        p = data[:, i] #i+1
        meanA = p[0:preBin].mean()
        p -= meanA
        p *= ADCconv
        Pulse[i][nE][:] = p[:]
        stdA[i].Fill(p[0:preBin].std())
        maxA[i].Fill(p.max())
        maxTail[i].Fill(p[nSamples - postBin:nSamples].max())
        minA[i].Fill(p.min())
        stdTail[i].Fill(p[nSamples - postBin:nSamples].std())
    if (nE+1)%100 == 0:
        print("events read: "+str(nE+1))
    if nE+1 == CutOff:
        break

#5 magical templates will be generated. From all those 5 chosen one only one will deserve the attention
templates = np.zeros((5, nChannels, nSamples))

freq1, throw1 = qp.calc_psd(p,sf, folded_over=False)
freq2 = freq1[0:nSamples-1]

int1 = [.01e-6, 0.2e-6]
int2 = [0.2e-6,0.5e-6]
int3 = [0.5e-6, 2e-6]
int4 = [2e-6, 3e-6]
int5 = [3e-6, 4e-6]

#cuts defined for each channel
# stdmaxpre, stdmaxpost, mincut, maxtailcut
cuts = {0: [0.05e-6, 0.05e-6, -0.1e-6, 1.5e-6],
        1: [0.03e-6, 0.03e-6, -0.14e-6, 0.5e-6],
        2: [0.04e-6, 0.04e-6, -0.15e-6, 0.5e-6],
        3: [0.03e-6, 0.04e-6, -0.14e-6, 0.5e-6],
        4: [0.04e-6, 0.04e-6, -0.14e-6, 1.5e-6],
        5: [0.04e-3, 0.04e-3, -0.14e-6, 1.5e-6],
        6: [0.04e-3, 0.04e-3, -0.14e-6, 1.5e-6],
        7: [0.04e-3, 0.04e-3, -0.14e-6, 1.5e-6],
        8: [0.04e-3, 0.04e-3, -0.14e-6, 1.5e-6],
        9: [0.04e-3, 0.04e-3, -0.14e-6, 1.5e-6],
        10: [0.04e-3, 0.04e-3, -0.14e-6, 1.5e-6],
        11: [0.04e-3, 0.04e-3, -0.14e-6, 1.5e-6] }
for nE in range(CutOff):
    if nE % 100 == 0:
        print(nE)

        for i in range(nChannels):
            p = Pulse[i][nE][:]

            if p[0:preBin].std() < cuts[i][0] and p[0:preBin].std() > 1e-8 and cuts[i][1] > p[nSamples - postBin:nSamples].std() > 1e-8 and p.min() > cuts[i][2] \
                    and p[nSamples - postBin:nSamples].max() < cuts[i][3]:

                if p.max() > int1[0] and p.max() < int1[1]:
                    templates[0][i][:] += p
                if p.max() > int2[0] and p.max() < int2[1]:
                    templates[1][i][:] += p
                if p.max() > int3[0] and p.max() < int3[1]:
                    templates[2][i][:] += p
                if p.max() > int4[0] and p.max() < int4[1]:
                    templates[3][i][:] += p
                if p.max() > int5[0] and p.max() < int5[1]:
                    templates[4][i][:] += p

for n in range(5):
    for i in range(nChannels):
        if templates[n][i][:].max() != 0:
            templates[n][i][:] /= templates[n][i][:].max()
        else:
            print("No template "+str(n)+" was not generated for Ch "+str(i))
        pulseGraph = ROOT.TGraph(nSamples, time2, templates[n][i][:])
        pulseGraph.SetName("pulseCh"+str(i)+"_"+str(n))
        pulseGraph.Write()

saveF.Write()
stop_time = process_time()
print("Total time needed to process: "+str(round((stop_time-start_time)/60))+" min")

