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

ptLSB = 0.25 # 1 in hw is 0.25 GeV
etaLSB = 0.0043633231 # 1 in hw is 0.00436 (in PF units)
phiLSB = 0.0043633231
nPayloadFrames = 15 # after how many frames a new event begins
frameStart = 89 # MAKE SURE THIS POINTS TO THE FIRST FRAME WHERE EVENTS CAN BE FOUND
jetLinkIndices = range(0, 9) # iterable containing links carrying jets
sumLinkIndex = 18 # link carring jets
dR = 0.01 # matching radius for jets
sumsTolerance = 0.05 # error margin for met and mht, eq. to 5% margin
eventsPerFile = 50 # events per file
nEvents = 50000 # number of events to analyse
hardwareOutputFile = TFile("DataEmuComparison.root", "RECREATE") # root file where data-emu stuff will be saved
emulatorOutputFile = TFile("CMSSWSums.root") # root file where to read emulator CMSSW output
nTxFiles = 1000

# DO NOT EDIT, data containers for later
hwJets = [] 
hwSum = [0, 0, 0] #ht, met, mht
hwSums = []
emJets = []
emSum = [0, 0, 0] #ht, met, mht
emSums = []
hwData = []
emData = []
hwDataNoZ = []
emDataNoZ = []
nHw = 0
nEm = 0
nEvNoZ = 0

# FRAME UNPACKING

# iterate through each file
for iTx in range(0, nTxFiles):
  with open("output_ttbar200/ttbar_pu200_pattern_" + str(iTx) + "_output.txt", "r") as inFile:
    frameIt = -1
    # find a line with 1v in it
    for line in inFile:
      if('1v' in line):
        # count as a frame
        frameIt += 1
        # skip frames once we have analysed enough events
        if nHw >= nEvents: continue
        # skip frames that are not produced by cms data
        if frameIt >= frameStart + nPayloadFrames * eventsPerFile:
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

# writing to tree
hardwareOutputFile.cd()
hardwareJetSumsTree = TTree("HardwareEvents", "Tree with output from hardware in a event-like structure")
emulatorJetSumsTree = TTree("EmulatorEvents", "Tree with output from emulator in a event-like structure")
# these trees can be used for analysis in FAST
hardwareJetTree = TTree("HardwareJets", "Flat tree with jets in output from hardware")
emulatorJetTree = TTree("EmulatorJets", "Flat tree with jets in output from emulator")
sumsTree = TTree("Sums", "Flat tree with sums in output from hardware and emulator")

# preparing data holders

hardware_jetPt = array("f", 3*[0])
hardware_jetEta = array("f", 3*[0]) 
hardware_jetPhi = array("f", 3*[0]) 
hardware_MHT = array("f", [0]) 
hardware_MET = array("f", [0]) 
hardware_HT = array("f", [0]) 
hardware_length = array("i", [0])

emulator_jetPt = array("f", 3*[0])
emulator_jetEta = array("f", 3*[0]) 
emulator_jetPhi = array("f", 3*[0]) 
emulator_MHT = array("f", [0]) 
emulator_MET = array("f", [0]) 
emulator_HT = array("f", [0]) 
emulator_length = array("i", [0])

emulator_jetPt_flat = array("f", [0]) 
emulator_jetEta_flat = array("f", [0]) 
emulator_jetPhi_flat = array("f", [0]) 

hardware_jetPt_flat = array("f", [0]) 
hardware_jetEta_flat = array("f", [0]) 
hardware_jetPhi_flat = array("f", [0]) 

# preparing tree branches to write

hardwareJetSumsTree.Branch("length", hardware_length, "length/i")
hardwareJetSumsTree.Branch("pt", hardware_jetPt, "pt[length]/F")
hardwareJetSumsTree.Branch("eta", hardware_jetEta, "eta[length]/F")
hardwareJetSumsTree.Branch("phi", hardware_jetPhi, "phi[length]/F")
hardwareJetSumsTree.Branch("mht", hardware_MHT, "mht/F")
hardwareJetSumsTree.Branch("met", hardware_MET, "met/F")
hardwareJetSumsTree.Branch("ht", hardware_HT, "ht/F")

emulatorJetSumsTree.Branch("length", emulator_length, "length/i")
emulatorJetSumsTree.Branch("pt", emulator_jetPt, "pt[length]/F")
emulatorJetSumsTree.Branch("eta", emulator_jetEta, "eta[length]/F")
emulatorJetSumsTree.Branch("phi", emulator_jetPhi, "phi[length]/F")
emulatorJetSumsTree.Branch("mht", emulator_MHT, "mht/F")
emulatorJetSumsTree.Branch("met", emulator_MET, "met/F")
emulatorJetSumsTree.Branch("ht", emulator_HT, "ht/F")

hardwareJetTree.Branch("pt", hardware_jetPt_flat, "pt/F")
hardwareJetTree.Branch("eta", hardware_jetEta_flat, "eta/F")
hardwareJetTree.Branch("phi", hardware_jetPhi_flat, "phi/F")

emulatorJetTree.Branch("pt", emulator_jetPt_flat, "pt/F")
emulatorJetTree.Branch("eta", emulator_jetEta_flat, "eta/F")
emulatorJetTree.Branch("phi", emulator_jetPhi_flat, "phi/F")

sumsTree.Branch("hwMHT", hardware_MHT, "hwMHT/F")
sumsTree.Branch("hwMET", hardware_MET, "hwMET/F")
sumsTree.Branch("hwHT", hardware_HT, "hwHT/F")
sumsTree.Branch("emuMHT", emulator_MHT, "emuMHT/F")
sumsTree.Branch("emuMET", emulator_MET, "emuMET/F")
sumsTree.Branch("emuHT", emulator_HT, "emuHT/F")

# preparing CMSSW trees to read from



emulatorJetTree_CMSSW = emulatorOutputFile.Get("SaveJets/EmulatorJets")
emulatorMETTree_CMSSW = emulatorOutputFile.Get("SaveSums/genMETL1TMETTree")
emulatorHTTree_CMSSW = emulatorOutputFile.Get("SaveSums/genHTL1THTTree")
# UNCOMMENT FOR MHT, this loads the tree containing the MHT from CMSSW
# emulatorMHTTree_CMSSW = emulatorOutputFile.Get("SaveGenSumsAndL1Sums/genMHTL1TMHTTree")

emulatorJetTree_CMSSW.SetBranchAddress("pt", emulator_jetPt)
emulatorJetTree_CMSSW.SetBranchAddress("eta", emulator_jetEta)
emulatorJetTree_CMSSW.SetBranchAddress("phi", emulator_jetPhi)
emulatorJetTree_CMSSW.SetBranchAddress("length", emulator_length)
emulatorMETTree_CMSSW.SetBranchAddress("l1tMET", emulator_MET)
emulatorHTTree_CMSSW.SetBranchAddress("l1tHT", emulator_HT)
# UNCOMMENT FOR MHT, this loads the branch containing the MHT from CMSSW
# emulatorMHTTree_CMSSW.SetBranchAddress("l1tMHT", emulator_MHT)

#event loop 
for iEv in range(0, nEvents):

  # write hardware event
  hardware_length[0] = len(hwData[iEv])
  
  for iJet in range(0, 3):
    # write hw event in array structure
    hardware_jetPt[iJet]  = hwData[iEv][iJet][0] if iJet < hardware_length[0] else 0
    hardware_jetEta[iJet] = hwData[iEv][iJet][1] if iJet < hardware_length[0] else 0
    hardware_jetPhi[iJet] = hwData[iEv][iJet][2] if iJet < hardware_length[0] else 0
    if iJet < hardware_length[0]:
      # write hw event in flat tree
      hardware_jetPt_flat[0] = hardware_jetPt[iJet]
      hardware_jetEta_flat[0] = hardware_jetEta[iJet]
      hardware_jetPhi_flat[0] = hardware_jetPhi[iJet]
      hardwareJetTree.Fill()

  hardware_HT[0] = hwSums[iEv][0] 
  hardware_MET[0] = hwSums[iEv][1]
  hardware_MHT[0] = hwSums[iEv][2]
  hardwareJetSumsTree.Fill()

  # get emulator event
  emulatorJetTree_CMSSW.GetEntry(iEv)
  emulatorMETTree_CMSSW.GetEntry(iEv)
  emulatorHTTree_CMSSW.GetEntry(iEv)
  # UNCOMMENT FOR MHT, this retrieves the CMSSW MHT for that event
  # emulatorMHTTree_CMSSW.GetEntry(iEv)

  # write emulator event
  
  emJets = []
  for iJet in range(0, emulator_length[0]):
    # writing emu event to flat tree
    emulator_jetPt_flat[0] = emulator_jetPt[iJet]
    emulator_jetEta_flat[0] = emulator_jetEta[iJet]
    emulator_jetPhi_flat[0] = emulator_jetPhi[iJet]
    emulatorJetTree.Fill()
  
  # preparing data structure for jet comparison
    emJets.append([emulator_jetPt[iJet], emulator_jetEta[iJet], emulator_jetPhi[iJet]])
  
  # saving emu event to tree
  emulatorJetSumsTree.Fill()

  sumsTree.Fill()
  nEm+=1
  # preparing data struct for emu-hw agreement
  emSums.append([emulator_HT[0], emulator_MET[0], emulator_MHT[0]])
  emData.append(emJets)

hardwareOutputFile.cd()
hardwareJetSumsTree.Write()
emulatorJetSumsTree.Write()
hardwareJetTree.Write()
emulatorJetTree.Write()
sumsTree.Write()
hardwareOutputFile.Close()


# writing emulator events and obtaining data-emu agreement for jet and sums

nDiff_Jet = 0
nDiff_HT = 0
nDiff_MET = 0
nDiff_MHT = 0

nInterestingEvents_Jet = 0
nInterestingEvents_HT = 0
nInterestingEvents_MET = 0
nInterestingEvents_MHT = 0
for evIt in range(0,nHw):

  # Computing jet agreement %
  # marking event with different number of jets as bad
  if len(hwData[evIt]) != len(emData[evIt]):
    nInterestingEvents_Jet += 1
    nDiff_Jet+=1
  else:
    # skipping empty events
    if len(hwData[evIt]) == 0: continue
    nInterestingEvents_Jet += 1
    goodJet=0
    # jet matching
    for hwJet in hwData[evIt]:
      for emJet in emData[evIt]:
        if hwJet[0] == emJet[0]:
          if (hwJet[1]-emJet[1])<dR:
            if (hwJet[2]-emJet[2])<dR:
              goodJet+=1
    # mark event as bad if not all jets are matched
    if goodJet < len(hwData[evIt]):
      nDiff_Jet+=1

  # Computing sum agreement %
  # skip null sums
  if not (emSums[evIt][0] == 0 and hwSums[evIt][0] == 0):
    nInterestingEvents_HT += 1
    # ht is bad if em!=hw
    if (emSums[evIt][0] != hwSums[evIt][0]):
      nDiff_HT += 1
  
  if not (emSums[evIt][1] == 0 and hwSums[evIt][1] == 0):
    nInterestingEvents_MET += 1
    # met is bad if em-hw/em > sumsTolerance% (to take account for sqrt disagreements)
    # met is bad if em = 0 and hw != 0 and viceversa
    if emSums[evIt][1] > 0 and hwSums[evIt][1] > 0 and abs((emSums[evIt][1] - hwSums[evIt][1]) / emSums[evIt][1]) > sumsTolerance:
      nDiff_MET += 1
  
  if not (emSums[evIt][2] == 0 and hwSums[evIt][2] == 0):
    nInterestingEvents_MHT += 1
    # mht is good is em-hw/em < sumsTolerance% (to take account for sqrt disagreements)
    # mht is bad if em = 0 and hw != 0 and viceversa
    if emSums[evIt][2] > 0 and hwSums[evIt][1] > 0 and abs((emSums[evIt][2] - hwSums[evIt][2]) / emSums[evIt][2]) > sumsTolerance:
      nDiff_MHT += 1

if (nInterestingEvents_Jet > 0): print("\n\nJET:\nnEvent = " + str(nInterestingEvents_Jet) + "\nnDiff = " + str(nDiff_Jet) + "\nGood events = " + str((1-float(nDiff_Jet)/float(nInterestingEvents_Jet))*100) + "%")
if (nInterestingEvents_HT > 0): print("\n\nHT:\nnEvent = " + str(nInterestingEvents_HT) + "\nnDiff = " + str(nDiff_HT) + "\nGood events = " + str((1-float(nDiff_HT)/float(nInterestingEvents_HT))*100) + "%")
if (nInterestingEvents_MET > 0): print("\n\nMET:\nnEvent = " + str(nInterestingEvents_MET) + "\nnDiff = " + str(nDiff_MET) + "\nGood events = " + str((1-float(nDiff_MET)/float(nInterestingEvents_MET))*100) + "%")
if (nInterestingEvents_MHT > 0): print("\n\nMHT:\nnEvent = " + str(nInterestingEvents_MHT) + "\nnDiff = " + str(nDiff_MHT) + "\nGood events = " + str((1-float(nDiff_MHT)/float(nInterestingEvents_MHT))*100) + "%")


# plt.style.use(hep.cms.style.ROOT)
# fig, axs =   plt.subplots(2,3, figsize=s(20, 12), gridspec_kw={'height_ratios': [3, 1]})

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
