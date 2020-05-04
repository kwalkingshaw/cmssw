#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import math
import numpy
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from ROOT import TFile, TTree
from array import array
# import mplhep as hep

ptLSB = 0.25;
etaLSB = 0.0043633231;
phiLSB = 0.0043633231;
hwJets = []
hwSum = [0, 0, 0] #ht, met, phi
hwSums = []
emJets = []
emSum = [0, 0, 0] #ht, met, phi
emSums = []
hwData = []
emData = []
hwDataNoZ = []
emDataNoZ = []
nHw = 0
nEm = 0
nEvNoZ = 0
nPayloadFrames = 15
frameStart = 91
nTxFiles = 1
jetLinkIndices = range(0, 9)
sumLinkIndex = 18
dR = 0.01
eventsPerFrame = 50
nEvents = 50

# iterate through each file
for hwTx in range(0,nTxFiles):
  with open("datFullHwn/" + str(hwTx) + "/tx_summary.txt", "r") as inFile:
    frameIt = -1
    # find a line with 1v in it
    for line in inFile:
      if('1v' in line):
        # count as a frame
        frameIt += 1
        # skip frames once we have analysed enough events
        if nHw >= nEvents: continue
        # skip frames that are not produced by cms data
        if frameIt >= frameStart + nPayloadFrames * (eventsPerFrame + 1):
          continue
        # skip frames before the actual beginning of data
        if frameIt < frameStart: continue
        # splitting links, and excluding the frame header
        linkData = line.split('1v')[1:]
        # obtaning a list of link containing the jets without spaces and new lines
        jetLinks = [linkData[linkIdx].replace(' ','').replace('\n','') for linkIdx in jetLinkIndices]
        sumLink = linkData[sumLinkIndex].replace(' ','').replace('\n','')
        for jetLink in jetLinks:
          # lower jet unpacking if pt not null
          if int(jetLink, 16) & 0xffff:
            jet = [(int(jetLink[8:],16)&0xffff)*ptLSB,
                     ((((int(jetLink[8:],16)>>24)&0xff)*19)+9)*etaLSB,
                     ((((int(jetLink[8:],16)>>16)&0xff)*20)+10)*phiLSB]
            hwJets.append(jet)
          # upper jet unpacking if pt not null
          if (int(jetLink, 16)>>32) & 0xffff:
            jet = [(int(jetLink[:8],16)&0xffff)*ptLSB,
                 ((((int(jetLink[:8],16)>>24)&0xff)*19)+9)*etaLSB,
                 ((((int(jetLink[:8],16)>>16)&0xff)*20)+10)*phiLSB]
            hwJets.append(jet)
        if int(sumLink, 16) != 0:
          hwSum = [(int(sumLink,16)&0xffff)*ptLSB, # ht
            ((int(sumLink,16) >> 16)&0xffff)*ptLSB, # met
            ((int(sumLink,16) >> 32)&0xffff)*ptLSB #mht
          ]
        # handles the end of an event
        if (frameIt - frameStart) % nPayloadFrames == nPayloadFrames - 1:
          # counts events in hw
          nHw+=1
          # add a null jet if no jet has been found, just as a place holder
          hwData.append(hwJets)
          del hwJets
          hwJets = []
          
          hwSums.append(hwSum)
          hwSum = [0, 0, 0]

hardwareOutputFile = TFile("hardwareoutput.root", "RECREATE")
hardwareJetSumsTree = TTree("hardwareOutput", "Tree with output from hardware")

hardware_jetPt = array("f", 3*[0])
hardware_jetEta = array("f", 3*[0]) 
hardware_jetPhi = array("f", 3*[0]) 
hardware_jetMHT = array("f", [0]) 
hardware_jetMET = array("f", [0]) 
hardware_jetHT = array("f", [0]) 
hardware_length = array("i", [0])

hardwareJetSumsTree.Branch("length", hardware_length, "length/i")
hardwareJetSumsTree.Branch("pt", hardware_jetPt, "pt[length]/F")
hardwareJetSumsTree.Branch("eta", hardware_jetEta, "eta[length]/F")
hardwareJetSumsTree.Branch("phi", hardware_jetPhi, "phi[length]/F")
hardwareJetSumsTree.Branch("mht", hardware_jetMHT, "mht/F")
hardwareJetSumsTree.Branch("met", hardware_jetMET, "met/F")
hardwareJetSumsTree.Branch("ht", hardware_jetHT, "ht/F")

for iEv in range(0, nEvents):
  hardware_length[0] = len(hwData[iEv])
  
  for iJet in range(0, 3):
    hardware_jetPt[iJet]  = hwData[iEv][iJet][0] if iJet < hardware_length[0] else 0
    hardware_jetEta[iJet] = hwData[iEv][iJet][1] if iJet < hardware_length[0] else 0
    hardware_jetPhi[iJet] = hwData[iEv][iJet][2] if iJet < hardware_length[0] else 0
  
  hardware_jetHT[0] = hwSums[iEv][0] 
  hardware_jetMET[0] = hwSums[iEv][1]
  hardware_jetMHT[0] = hwSums[iEv][2]

  hardwareJetSumsTree.Fill()

hardwareJetSumsTree.Write()
hardwareOutputFile.Close()

emulatorOutputFile = TFile("CMSSWSums.root")

emulatorJetTree = emulatorOutputFile.Get("SaveJets/EmulatorJets")

event_jetPt = array("f", 3*[0])
event_jetEta = array("f", 3*[0]) 
event_jetPhi = array("f", 3*[0]) 
length = array("i", [0])

emulatorJetTree.SetBranchAddress("pt", event_jetPt)
emulatorJetTree.SetBranchAddress("eta", event_jetEta)
emulatorJetTree.SetBranchAddress("phi", event_jetPhi)
emulatorJetTree.SetBranchAddress("length", length)

for iEv in range(0, nEvents):
  # extracting jets
  emulatorJetTree.GetEntry(iEv)
  # interting a place holder if 
  if length[0] == 0 : 
    emJets = []
  else:
    emJets = [[event_jetPt[iJet], event_jetEta[iJet], event_jetPhi[iJet]] 
      for iJet in range(0, length[0])]
  nEm+=1
  emData.append(emJets)

nDiff = 0

nInterestingEvents = 0
for evIt in range(0,nHw):
  if len(hwData[evIt]) != len(emData[evIt]):
    nInterestingEvents += 1
    nDiff+=1
    continue
  if len(hwData[evIt]) == 0: continue
  nInterestingEvents += 1
  goodJet=0
  for hwJet in hwData[evIt]:
    for emJet in emData[evIt]:
      if hwJet[0] == emJet[0]:
        if (hwJet[1]-emJet[1])<dR:
          if (hwJet[2]-emJet[2])<dR:
            goodJet+=1
  if goodJet < len(hwData[evIt]):
    nDiff+=1



# print("\n\n=====================================================================================")
# print("\t\tFirmware Events: " + str(nHw) + "\t\t" + "Emulator Events: " + str(nEm))
# print("=====================================================================================")
# print("\t\tpT\t" + "eta\t" + "phi\t\t" + "pT\t" + "eta\t" + "phi\t")
# print("=====================================================================================")


# for evIt in range(0,nHw):
#   if hwData[evIt][0][0] > 0:
#     hwDataNoZ.append(hwData[evIt])
#   if emData[evIt][0][0] > 0:
#     emDataNoZ.append(emData[evIt])
#   nEvNoZ+=1


# for evIt in range(0,nHw):
#   if hwData[evIt][0][0] ==0 and emData[evIt][0][0] == 0:
#     continue
#   jetCount=0
#   jetDiff = len(hwData[evIt]) - len(emData[evIt])
#   print("")
#   if jetDiff==0:
#     for jetIt in range(len(hwData[evIt])):
#       print(str(evIt) + "\t\t" + str(hwData[evIt][jetIt][0]) + "\t" + str(hwData[evIt][jetIt][1])[:4] + "\t" + str(hwData[evIt][jetIt][2])[:4] + "\t\t" +
#           str(emData[evIt][jetIt][0]) + "\t" + str(emData[evIt][jetIt][1])[:4] + "\t" + str(emData[evIt][jetIt][2])[:4])
#   if jetDiff>0:
#     for jetIt in range(len(hwData[evIt])):
#       jetCount+=1
#       if jetCount > len(emData[evIt]):
#         emData[evIt].append([0,0,0])
#       print(str(evIt) + "\t\t" + str(hwData[evIt][jetIt][0]) + "\t" + str(hwData[evIt][jetIt][1])[:4] + "\t" + str(hwData[evIt][jetIt][2])[:4]  + "\t\t" +
#           str(emData[evIt][jetIt][0]) + "\t" + str(emData[evIt][jetIt][1])[:4] + "\t" + str(emData[evIt][jetIt][2])[:4])
#   if jetDiff<0:
#     for jetIt in range(len(emData[evIt])):
#       jetCount+=1
#       if jetCount > len(hwData[evIt]):
#         hwData[evIt].append([0,0,0])
#       print(str(evIt) + "\t\t" + str(hwData[evIt][jetIt][0]) + "\t" + str(hwData[evIt][jetIt][1])[:4] + "\t" + str(hwData[evIt][jetIt][2])[:4]  + "\t\t" +
#           str(emData[evIt][jetIt][0]) + "\t" + str(emData[evIt][jetIt][1])[:4] + "\t" + str(emData[evIt][jetIt][2])[:4])

print("\n\nnEvent = " + str(nInterestingEvents) + "\nnDiff = " + str(nDiff) + "\nGood events = " + str((1-float(nDiff)/float(nInterestingEvents))*100) + "%")


# plt.style.use(hep.cms.style.ROOT)
# fig, axs =   plt.subplots(2,3, figsize=(20, 12), gridspec_kw={'height_ratios': [3, 1]})

# fig.patch.set_facecolor( '#ffffff')


# nPtHw,  bPtHw  = np.histogram([jet[0] for event in hwDataNoZ for jet in event], bins=40, range=(0,200))
# nEtaHw, bEtaHw = np.histogram([jet[1] for event in hwDataNoZ for jet in event], bins=18, range=(0,1.5))
# nPhiHw, bPhiHw = np.histogram([jet[2] for event in hwDataNoZ for jet in event], bins=8,  range=(0,0.7))

# meansPt  = [0.5*(bPtHw[i]  + bPtHw[i+1])  for i in range(len(nPtHw))]
# meansEta = [0.5*(bEtaHw[i] + bEtaHw[i+1]) for i in range(len(nEtaHw))]
# meansPhi = [0.5*(bPhiHw[i] + bPhiHw[i+1]) for i in range(len(nPhiHw))]

# nPtEm  = axs[0,0].hist([jet[0] for event in emDataNoZ for jet in event], bins=40, range=(0,200), histtype='step', linewidth=1.5, label='Software', color='#0485d1', zorder=0)[0]
# nEtaEm = axs[0,1].hist([jet[1] for event in emDataNoZ for jet in event], bins=18, range=(0,1.5), histtype='step', linewidth=1.5, label='Software', color='#0485d1', zorder=0)[0]
# nPhiEm = axs[0,2].hist([jet[2] for event in emDataNoZ for jet in event], bins=8,  range=(0,0.7), histtype='step', linewidth=1.5, label='Software', color='#0485d1', zorder=0)[0]

# axs[0,0].scatter(meansPt,  nPtHw,  label='Hardware', c='#000000', linewidths=0.5, s=25, marker='o')
# axs[0,1].scatter(meansEta, nEtaHw, label='Hardware', c='#000000', linewidths=0.5, s=25, marker='o')
# axs[0,2].scatter(meansPhi, nPhiHw, label='Hardware', c='#000000', linewidths=0.5, s=25, marker='o')


# axs[1,0].scatter(meansPt,  [(hw/em) if em>0 else 0 for hw,em in zip(nPtHw,nPtEm)] ,  c='#000000', linewidths=0.5, s=25, zorder=1)
# axs[1,1].scatter(meansEta, [(hw/em) if em>0 else 0 for hw,em in zip(nEtaHw,nEtaEm)], c='#000000', linewidths=0.5, s=25, zorder=1)
# axs[1,2].scatter(meansPhi, [(hw/em) if em>0 else 0 for hw,em in zip(nPhiHw,nPhiEm)], c='#000000', linewidths=0.5, s=25, zorder=1)

# axs[1,0].axhline(y=1, linewidth=1, linestyle='--', c='#929591')
# axs[1,1].axhline(y=1, linewidth=1, linestyle='--', c='#929591')
# axs[1,2].axhline(y=1, linewidth=1, linestyle='--', c='#929591')

# axs[1,0].set(ylim=(0.5,1.5))
# axs[1,1].set(ylim=(0.5,1.5))
# axs[1,2].set(ylim=(0.5,1.5))

# axs[0,0].set(ylabel="Events")
# axs[1,0].set(ylabel="Ratio")

# axs[0,0].legend(prop={'size': 20})
# axs[0,1].legend(prop={'size': 20})
# axs[0,2].legend(prop={'size': 20})

# ymaxPt  = max(np.concatenate([nPtHw,nPtEm]))
# ymaxEta = max(np.concatenate([nEtaHw,nEtaEm]))
# ymaxPhi = max(np.concatenate([nPhiHw,nPhiEm]))

# axs[0,0].set(xlim=(0,200))
# axs[0,1].set(xlim=(0,1.5))
# axs[0,2].set(xlim=(0,0.7))

# axs[0,0].set(ylim=(0,ymaxPt +(0.05*ymaxPt)))
# axs[0,1].set(ylim=(0,ymaxEta+(0.05*ymaxEta)))
# axs[0,2].set(ylim=(0,ymaxPhi+(0.05*ymaxPhi)))

# #axs[0,0].xaxis.labelpad =
# #axs[0,1].xaxis.labelpad = -5
# #axs[0,2].xaxis.labelpad = -5

# #axs[0,0].set_title("Histogrammed PF Jet FW vs EMU: pT, ttbar, 3900 events")  
# #axs[0,1].set_title("Histogrammed PF Jet FW vs EMU: Eta, ttbar, 3900 events")  
# #axs[0,2].set_title("Histogrammed PF Jet FW vs EMU: Phi, ttbar, 3900 events") 

# axs[0,0].set(xlabel="Jet $p_T$ (GeV)")
# axs[0,1].set(xlabel="Jet $\eta$")
# axs[0,2].set(xlabel="Jet $\phi$")

# #plt.xlabel(r'xlabel', ha='right', x=1.)
# #plt.ylabel(r'ylabel', ha='right', y=1.)

# #hep.cms.cmslabel(plt.gca(), data=False, paper=False, llabel='Phase-2 Simulation', rlabel=r'14 TeV  200 PU')

# plt.savefig('hwEmu.pdf')#, bbox_inches='tight')
# plt.show()
