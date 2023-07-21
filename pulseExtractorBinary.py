import struct
from random import randrange
from pathlib import Path
import shutil
import pandas as pd
import numpy as np
import os
import scipy.io
from scipy import signal
from time import process_time
import ROOT
import glob
from time import process_time
import warnings
def rotateMatrix(m):
    return m

#Creates array of sapphire_ch##_###N.txt files
def fileList(directory, fileN):
    flist1 = [glob.glob(directory+'sapphire_ch{}_{}.txt'.format(str(j).zfill(2),str(fileN).zfill(4)))
              for j in Channels]
    flist = [flist1[i][0] for i in range(len(flist1))]
    return flist

#Created directory if it does not exist
#Clears the directory if it exists
def createPath(path):
    if (not os.path.exists(path)):
        os.mkdir(path)
    else:
        pulselist = glob.glob(path+"/*")
        for i in pulselist:
            path1 = Path(i)
            if path1.is_file():
                os.remove(i)
            if path1.is_dir():
                shutil.rmtree(path1, ignore_errors=True)
            #os.rmdir(i)

start_time = process_time()
#warnings.filterwarnings("ignore", category = DeprecationWarning)
nEvents = 0
eventsToRead = 1000000000
totalCrossTime = 0
totalReadingTime = 0
#write down n channel numbers you want to read!!!!!!!!!!!!!!
global nChannels
nChannels = 16
#nChannels = 4
nNoise = 5000
global Channels
Channels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#Channels = [4,5,6,7]

#12 Channels cross-correlation channel parameters
prominace = {0:7000, 1:700, 2:700, 3:700, 4:1000, 5:1000, 6:1000,7:1000, 8:4500, 9:4500, 10:4500,11:4500,12:1500, 13:1500, 14:1500,15:1500}
triggered = {0:False, 1:True, 2:True, 3:True,4:False, 5:True, 6:True, 7:True,8:False, 9:True, 10:True, 11:True,12:False, 13:True, 14:True, 15:True}

#prominace = {0:700, 1:1000, 2:1000, 3:1000}
#triggered = {0:False, 1:True, 2:True, 3:True}


#PulseFolder = os.getcwd()
CoincidenceWindow = 500  #us
SamplingRate = 2 #us
CoincidenceIndex = int(CoincidenceWindow/SamplingRate)
pulseLength = 1024

#event number read for each channel
ChannelEvents = [0 for i in range(len(Channels))]

#reading dataset directory from the file
fd = open("dataset.txt", "r")
directory = fd.readline().strip()
directory1 = "/data1/Run39_NSC" 
pulseData = fd.readline().strip()
pulseData1 = "39211208_sapphire7_4_2_14_DG20_Digi64Ch_noTR_gndFoil_BIN2_SWTall_30min"

#Extracted pulses saved in the same directory
PulseFolder = directory+"/"+pulseData1+"_extracted_pulses"
NoiseFolder = directory+"/"+pulseData1+"_extracted_noise"

#Manual input for save directories
#PulseFolder = '/data/sandro/extractedPulsesTest/Pulse'
#NoiseFolder = '/data/sandro/extractedPulsesTest/Noise'

#Save in the same directory
#PulseFolder = os.getcwd()+"/Pulse"
#NoiseFolder = os.getcwd()+"/Noise"

#create path for noise and pulse
createPath(PulseFolder)
createPath(NoiseFolder)
print(directory+"/"+pulseData+"/")

#tF = ROOT.TFile(pulse+"_pulse"+".root")
tF = ROOT.TFile("39211203_4mm2_trgBCD_prom10000_calib_extracted_pulses_pulse.root")

#pulse templates for each channel (set input for each channel!!!)
templates = {}
#templates[0] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[1] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[2] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[3] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[4] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[5] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[6] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[7] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[8] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[9] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[10] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))
#templates[11] = np.array(np.frombuffer(tF.Get("pulse4A_template").GetY(), dtype=np.double))

templates[0] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[1] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[2] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[3] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[4] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[5] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[6] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[7] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[8] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[9] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[10] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[11] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[12] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[13] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[14] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))
templates[15] = np.array(np.frombuffer(tF.Get("pulseCh0_3").GetY(), dtype=np.double))

#-----------------------------------------------------------------------------------------------------------------------
for fileN in range(0,4):
    time0 = process_time()
    createPath(PulseFolder+"/sapphire_"+str(fileN))
    files = fileList(directory+"/"+pulseData+"/", int(fileN))
    #files = glob.glob(directory+"/"+pulse+"/*.txt")
    nEvent = 0
    headersize = 6 * 4
    minEventN = 1000000
    event_number_list = []
    # Read one channel data to determine settings. It assumes that ch 0 exists.
    for filen in files:
        txtfile = open(filen, "rb")
        isData = True
        while isData:
            headerb = txtfile.read(headersize)
            try:
                header = struct.unpack("IIIIII", headerb)
            except:
                isData = False
                break
            size = header[0]
            channel = header[3]
            eventN = header[4]
            eventT = header[5]
            if channel != Channels[0]:
                break
            if nEvent == 0:
                minEventN = eventN
            if minEventN > eventN:
                minEventN = eventN
            event_number_list.append(eventN)
            datab = txtfile.read(size - headersize)
            data = struct.unpack("H" * int(len(datab) / 2), datab)
            if len(data) <= 0:
                isData = False
                break
            nEvent += 1
            global nSamples
            nSamples = len(data)
        txtfile.close()

    event_number_list.sort() #just in case
    nEvents = nEvents + nEvent

    if nEvents >=eventsToRead:
        break
    if fileN ==0:
        print("There are " + str(len(Channels)) + " Channels.")
        print("Pulse length is: " + str(nSamples))
        print("channels: " + str(Channels))
    print("Number of Events in file: "+str(fileN)+" is " + str(nEvent))

    p = np.zeros(nSamples)
    s = (nChannels, nEvent, nSamples)
    pulse = np.zeros(s)
    pulseExtracted = np.zeros((nChannels, pulseLength))
    eventTT = np.zeros(nEvent)

    for filen in files:
        txtfile = open(filen, "rb")
        isData = True
        while isData:
            headerb = txtfile.read(headersize)
            try:
                header = struct.unpack("IIIIII", headerb)
            except:
                isData = False
                break
            size = header[0]
            channel = header[3]
            eventN = header[4]
            eventT = header[5]
            if channel not in Channels:
                break

            datab = txtfile.read(size - headersize)
            data = struct.unpack("H" * int(len(datab) / 2), datab)
            if len(data) <= 0:
                isData = False
                break
            if (eventN % 10 == 0):
                print(
                    "Reading Channel " + str(channel) + " Event number " + str(eventN) + " from file number " + str(
                        fileN))
            for i in range(nChannels):
                if channel == Channels[i]:
                    try:
                        pulse[i][event_number_list.index(eventN)] = data
                        ChannelEvents[i] += 1
                    except ValueError:
                        print("mismatch in event numbers")
        txtfile.close()

    #print("Event number in each file: "+str(nEvent))
    print("Total Number of events read so far" + str(nEvents))
    print("Event number for each channel:"+str(ChannelEvents))

    time1 = process_time()
    foldercounter = 0

    Xcorr,Xcorrhalf,peaks = dict(),dict(),dict()
    #nEvents = eventsToRead
    for nE in range(nEvent):
        if nE % 1 == 0:
            print(nE)

        #cross correlation and peaks for each channel
        for i in range(nChannels):
            p = pulse[i][nE][:]
            if triggered[i]:
                Xcorr[i] = signal.correlate(p,templates[i], method='fft')
                Xcorrhalf[i] = Xcorr[i][pulseLength:]
                peaks[i] = signal.find_peaks(Xcorrhalf[i], prominence=prominace[i])[0]
            else:
                peaks[i] = np.zeros(nChannels, dtype=int)

        #Identifing peaks not in coincidence
        peaksall1 = np.concatenate([peaks[i] for i in range(nChannels)], axis=0)
        peaksall1.sort()
        peaksall = np.unique(peaksall1)
        #print(peaksall)
        #print("Number of peaks found in"+str(nE)+"th slice:  "+str(len(peaksall)))

        realpeaks = np.zeros(0)
        peakscount = 0  # only for for loop below
        for j, peak in enumerate(peaksall):
            if peaksall[j]!=0:
                if peakscount == 0:
                    realpeaks = np.append(realpeaks, peaksall[j])
                    peakscount += 1
                elif peaksall[j] > realpeaks[peakscount - 1] + CoincidenceIndex:
                    realpeaks = np.append(realpeaks, peaksall[j])
                    peakscount += 1
        #print(realpeaks)
        print("Number of peaks not in coincidence found in" + str(nE) + "th slice:  " + str(len(realpeaks)))

        #Write each pulse into the file
        for k, peak in enumerate(realpeaks):
            foldercounter += 1
            if peak + pulseLength > nSamples or k == len(realpeaks) - 1 or k == len(realpeaks) - 2:
                break
            for i in range(nChannels):
                    if peak+pulseLength< nSamples:
                        pulseExtracted[i] = pulse[i][nE][int(peak):int(peak+pulseLength)]

            pulseExtracted1 = np.rot90(pulseExtracted)
            pulseExtracted2 = np.copy(pulseExtracted1)
            for i in range(nChannels):
                pulseExtracted2[:, i] = pulseExtracted1[:, i][::-1]

            if not os.path.isdir(PulseFolder+"/sapphire_"+str(fileN) + '/slice_' + str(nE)):
                os.mkdir(PulseFolder+"/sapphire_"+str(fileN) + '/slice_' + str(nE))
            np.savetxt(PulseFolder+"/sapphire_"+str(fileN) + '/slice_' + str(nE) + '/pulse_' + str(k) + '.npy',
                       pulseExtracted2)
    time2 = process_time()
    totalCrossTime += time2 - time1
    totalReadingTime += time1 - time0
    if fileN == 0:
        for nN in range(nNoise):
            nE = randrange(nEvent)
            c = randrange(nSamples-pulseLength)
            for i in range(nChannels):
                pulseExtracted[i] = pulse[i][nE][int(c):int(c + pulseLength)]
            pulseExtracted1 = np.rot90(pulseExtracted)
            #pulseExtracted2 = np.copy(pulseExtracted1)
            for i in range(nChannels):
                pulseExtracted2[:, i] = pulseExtracted1[:, i][::-1]
            np.savetxt(NoiseFolder + '/noise_' + str(nN) + '.npy',
                       pulseExtracted2)

stop_time = process_time()
print("Time to read: "+str(round((totalReadingTime)/60))+" min")
print("Cross Correlation time: "+str(round((totalCrossTime)/60))+" min")
print("Total time to process: "+str(round((stop_time-start_time)/60))+" min")
print("Total Events processed: "+int(nEvents))
