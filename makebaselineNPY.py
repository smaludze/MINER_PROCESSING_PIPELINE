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


def butter_bandpass(lowcut, fs, order=5):
    low = lowcut / (0.5 * fs)
    b, a = butter(order, low, btype='lowpass')
    return b, a


def butter_bandpass_filter(data, lowcut, fs, order=5):
    b, a = butter_bandpass(lowcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def passCut():
  return True

#write down n channel numbers you want to read!!!!!!!!!!!!!!
nDetectors = 4
nChannels = 16
Channels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

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

files = glob.glob(directory+"/"+noise+"/*.npy")

#tF = ROOT.TFile(pulse+"_pulse"+".root")
tF = ROOT.TFile("39211203_4mm2_trgBCD_prom10000_calib_extracted_pulses_pulse.root")
tF1 = ROOT.TFile("39211207_sapphire7_4DG20_Digi64Ch_noTR_BIN_1_extracted_pulses_pulse.root")

nF = ROOT.TFile(noise+"_noise.root")
#nF = ROOT.TFile("3820211105_sapphire903_new_digitizer_noise_1_noise.root")
saveF = ROOT.TFile(noise+"_noise_baseline"+".root","RECREATE")
#saveF = ROOT.TFile("test.root","RECREATE")

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

nEvents = len(files)
ChannelEvents = [0 for i in range(len(Channels))]
print("There are "+str(len(Channels))+" Channels. Total number of noise Events: "+str(nEvents))
#print("Pulse length is: "+str(nSamples))
#print("channels: "+str(Channels))

time = np.arange(nSamples,dtype=np.float64)
time2 = time*sample_interval
triggerTime = trigger_position*time2[nSamples-1]
sf = (1/(sample_interval*1e-9))

t = ROOT.TTree('data', 'tree with events')

ampOF,chi2OF = dict(), dict()

p = np.zeros(nSamples)
noise = np.zeros((nChannels, nEvents, nSamples))
for i in range(nChannels):
    ampOF[i] = np.zeros(1, dtype=np.float64)
    chi2OF[i] = np.zeros(1, dtype=np.float64)

    t.Branch('ampOF'+str(i), ampOF[i], 'ampOF'+str(i)+'/D')
    t.Branch('chi2OF' + str(i), chi2OF[i], 'chi2OF' + str(i) + '/D')

for n,filen in enumerate(files):
    data = np.loadtxt(filen)
    for i in range(nChannels):
        p = data[:, i] #i+1
        meanA = p.mean()
        p -= meanA
        p *= ADCconv
        noise[i][n] = p[:]

print("noise: "+str(nEvents))
for nE in range(nEvents):
    pass_cut = passCut()
    if pass_cut:
        for i in range(nChannels):
            OF = qp.OptimumFilter(noise[i][nE]+templates[i], templates[i], NoiseTemplates[i], sf)
            amp, chi2 = OF.ofamp_nodelay()
            ampOF[i][0] = amp-1
            chi2OF[i][0] = chi2

    t.Fill()

saveF.Write()
saveF.Close()

