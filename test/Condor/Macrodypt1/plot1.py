#! /usr/bin/env python
import sys
import os
import subprocess
from array import array
from ROOT import TH1D,TH2D,TFile,TMath,TCanvas,THStack,TLegend,TPave,TLine,TLatex
from ROOT import gROOT,gStyle,gPad,gStyle
from ROOT import Double,kBlue,kRed,kOrange,kMagenta,kYellow,kCyan,kGreen,kGray,kBlack,kTRUE

gROOT.Macro("~/rootlogon.C")
gStyle.SetOptStat(0)

# ===============
# options
# ===============
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--Lumi', metavar='D', type='float', action='store',
                  default= 2263.,   # 2200., 538.134,
                  dest='Lumi',
                  help='Data Luminosity in pb-1')

parser.add_option('--globalSF', metavar='SF', type='float',
                  default = 1.0,
                  dest='globalSF',
                  help='Global trigger SF (%default default)')

parser.add_option('--plotDir', metavar='P', type='string', action='store',
                  default='test',
                  dest='plotDir',
                  help='output directory of plots')

parser.add_option('--skimType', metavar='S', type='string', action='store',
                  default='CR_Zelel',
                  dest='skimType',
                  help='Skim type: CR_Zelel, CR_Zmumu, SR_Zelel, SR_Zmumu')

parser.add_option('--processDir', metavar='pD', type='string', action='store',
                  default='ana/cat',
                  dest='processDir',
                  help='directory to read histograms from: ana/pre, ana/cnt, ana/sig')

parser.add_option('--var', metavar='T', type='string', action='store',
                  default='ht_zsel',#cutflow, st
                  dest='var',
                  help='variable to plot')

parser.add_option('--sys', metavar='T', type='string', action='store',
                  default='nominal',
                  dest='sys',
                  help='nominal, BTagSFup, BTagSFdown, ScaleUp, ScaleDown, MatchUp, MatchDown')

parser.add_option('--verbose',action='store_true',
                  default=True,
                  dest='verbose',
                  help='verbose switch')

parser.add_option('--rebin', metavar='T', type='int', action='store',
                  default='1',
                  dest='rebin',
                  help='rebin the histograms')

(options,args) = parser.parse_args()
# ==========end: options =============
var = options.var
lumi = options.Lumi
gSF  = options.globalSF
rebinS = options.rebin
pDir = options.processDir
plotDir = options.plotDir
skimType = options.skimType
verbose = options.verbose

if 'elel' in skimType: title = 'e^{#pm}e^{#mp}+jets'
elif 'mumu' in skimType: title = '#mu^{#pm}#mu^{#mp}+jets'
else: title = ''
# ==========add the input ============

execfile("input.py")

# === prepare the input labels and legends ===========
dataLabel     = 'Data_'
topLabel      = 'Top_'
dyLabel       = 'DY_'
dibosonLabel  = 'Diboson_'
wjLabel       = 'WJets_'
sTLabel       = 'ST_'
vvLabel       = 'VV_'

dataLeg       = 'Data'
topLeg        = 't#bar{t}'
dyLeg         = 'Drell-Yan'
dibosonLeg    = 'Diboson'
wjLeg         = 'W+Jets'
sTLeg         = 'Single top'
vvLeg         = 'Diboson'
# === create structure ============
data     = [
            #[f_Data_Oct2015, 1., 1., 1.], # this corresponds to .53 fb-1
            [f_Data_PromptReco, 1., 1., 1.],  #addign this should give 2.2 fb-1  
           ]

top      = [[f_ttbar,         Top_xs,            Top_num,            lumi]]

#dy       = [[f_DYmcnlo,       DY_xs,             DYmcnlo_num,        lumi]]

#dy       = [[f_DYmad,       DY_xs,             DYmad_num,        lumi]] 

dy       = [
            [f_DY100to200,    DY100to200_xs,      DY100to200_num,    lumi],
            [f_DY200to400,    DY200to400_xs,      DY200to400_num,    lumi],
            [f_DY400to600,    DY400to600_xs,      DY400to600_num,    lumi],
            [f_DY600toInf,    DY600toInf_xs,      DY600toInf_num,    lumi],    
            ]
diboson  =[
             [      f_WW,             WW_xs,              WW_num,    lumi],
             [      f_WZto2,           WZto2_xs,          WZto2_num,    lumi],
             [      f_WZto3,           WZto3_xs,          WZto3_num,    lumi],
             [      f_ZZto2,           ZZto2_xs,          ZZto2_num,    lumi],
             [      f_ZZto4,           ZZto4_xs,          ZZto4_num,    lumi],   
             ] 
#wjets    = [
#            [f_WJ100to200,    WJ100to200_xs,      WJ100to200_num,    lumi],
#            [f_WJ200to400,    WJ200to400_xs,      WJ200to400_num,    lumi],
#            [f_WJ400to600,    WJ400to600_xs,      WJ400to600_num,    lumi],
#            [f_WJ600to800,    WJ600to800_xs,      WJ600to800_num,    lumi],     
#            [f_WJ800to1200,   WJ800to1200_xs,     WJ800to1200_num,   lumi],
#            [f_WJ1200to2500,  WJ1200to2500_xs,    WJ1200to2500_num,  lumi],
#            [f_WJ2500toInf,   WJ2500toInf_xs,     WJ2500toInf_num,   lumi], 
#            ]

#st       = [
#            [f_ST_tW_top,     ST_tW_top_xs,       ST_tW_top_num,     lumi],
#            [f_ST_tW_antitop, ST_tW_antitop_xs,   ST_tW_antitop_num, lumi],
#            [f_ST_t,          ST_t_xs,            ST_t_num,          lumi],
#            [f_ST_t_ex1,      ST_t_xs,            ST_t_ex1_num,      lumi],
#            [f_ST_s,          ST_s_xs,            ST_s_num,          lumi], 
#           ]

#vv       = [
#            [f_ZZTo2L2Nu,     ZZTo2L2Nu_xs,       ZZTo2L2Nu_num,    lumi],
#            [f_WZTo2L2Q,      WZTo2L2Q_xs,        WZTo2L2Q_num,     lumi],
#            [f_WWTo2L2Nu,     WWTo2L2Nu_xs,       WWTo2L2Nu_num,    lumi],
#           ]
#tZtZ_700 = [[f_BpBp_tZtZ_700, TpTp700_xs,         TpTp700_num,       lumi]]
#tZbW_700 = [[f_BpTp_tZbW_700, TpTp700_xs,         TpTp700_num,       lumi]]
#tZtH_700 = [[f_BpTp_tZtH_700, TpTp700_xs,         TpTp700_num,       lumi]]
tZtZ_800 = [[f_TpTp_tZtZ_800, TpTp800_xs,         TpTp800tZtZ_num,       lumi]]
#tZbW_800 = [[f_TpTp_tZbW_800, TpTp800_xs,         TpTp800tZbW_num,       lumi]]
tZtH_800 = [[f_TpTp_tZtH_800, TpTp800_xs,         TpTp800tZtH_num,       lumi]]
tZtZ_1000 = [[f_TpTp_tZtZ_1000, TpTp1000_xs,         TpTp1000tZtZ_num,       lumi]]

#tZbW_1000 = [[f_TpTp_tZbW_1000, TpTp1000_xs,         TpTp1000tZbW_num,       lumi]]
tZtH_1000 = [[f_TpTp_tZtH_1000, TpTp1000_xs,         TpTp1000tZtH_num,       lumi]]

tZtZ_1200 = [[f_TpTp_tZtZ_1200, TpTp1200_xs,         TpTp1200tZtZ_num,       lumi]]
#tZbW_1200 = [[f_TpTp_tZbW_1200, TpTp1200_xs,         TpTp1200tZbW_num,       lumi]]
tZtH_1200 = [[f_TpTp_tZtH_1200, TpTp1200_xs,         TpTp1200tZtH_num,       lumi]]

#bZbZ_700 = [[f_BpBp_bZbZ_700, BpBp700_xs,         BpBp700_num,       lumi]]
#bZtW_700 = [[f_BpBp_bZtW_700, BpBp700_xs,         BpBp700_num,       lumi]]
#bZbH_700 = [[f_BpBp_bZbH_700, BpBp700_xs,         BpBp700_num,       lumi]]
#bZbZ_800 = [[f_BpBp_bZbZ_800, BpBp800_xs,         BpBp800_num,       lumi]]
#bZtW_800 = [[f_BpBp_bZtW_800, BpBp800_xs,         BpBp800_num,       lumi]]
#bZbH_800 = [[f_BpBp_bZbH_800, BpBp800_xs,         BpBp800_num,       lumi]]
#bZbZ_1000 = [[f_BpBp_bZbZ_1000, BpBp1000_xs,         BpBp1000_num,       lumi]]
#bZtW_1000 = [[f_BpBp_bZtW_1000, BpBp1000_xs,         BpBp1000_num,       lumi]]
#bZbH_1000 = [[f_BpBp_bZbH_1000, BpBp1000_xs,         BpBp1000_num,       lumi]]


#bZbZ_1200 = [[f_BpBp_bZbZ_1200, BpBp1200_xs,         BpBp1200_num,       lumi]]
#bZtW_1200 = [[f_BpBp_bZtW_1200, BpBp1200_xs,         BpBp1200_num,       lumi]]
#bZbH_1200 = [[f_BpBp_bZbH_1200, BpBp1200_xs,         BpBp1200_num,       lumi]]

#print "----> step 1\n"
h_data     = getHisto(dataLabel,       dataLeg,        pDir, var,  data,     kBlack,     verbose)
#print "----> step 2\n"
#print " -----> ", h_data.GetName()
h_top      = getHisto(topLabel,        topLeg,         pDir, var,  top,      8,          verbose)
h_dy       = getHisto(dyLabel,         dyLeg,          pDir, var,  dy,       90,         verbose)
h_diboson  = getHisto(dibosonLabel,  dibosonLeg,       pDir, var,  diboson,   kBlue,         verbose)
#h_wjets    = getHisto(wjLabel,         wjLeg,          pDir, var,  wjets,    kBlue,      verbose)
#h_st       = getHisto(sTLabel,         sTLeg,          pDir, var,  st,       kCyan,      verbose)
#h_vv       = getHisto(vvLabel,         vvLeg,          pDir, var,  vv,       kRed,       verbose)

#h_tZtZ_700 = getHisto('TT_tZtZ_M700_', 'TT_tZtZ_M700', pDir, var,  tZtZ_700, kRed,    verbose)
#h_tZbW_700 = getHisto('TT_tZbW_M700_', 'TT_tZbW_M700', pDir, var,  tZbW_700, kCyan+4,  verbose)
#h_tZtH_700 = getHisto('TT_tZtH_M700_', 'TT_tZtH_M700', pDir, var,  tZtH_700, kBlue+4, verbose)
h_tZtZ_800 = getHisto('TT_tZtZ_M800_', 'TT_tZtZ_M800', pDir, var,  tZtZ_800, kGreen+4,    verbose)
#h_tZbW_800 = getHisto('TT_tZbW_M800_', 'TT_tZbW_M800', pDir, var,  tZbW_800, kGreen+3,  verbose)
h_tZtH_800 = getHisto('TT_tZtH_M800_', 'TT_tZtH_M800', pDir, var,  tZtH_800, kGreen+2, verbose)
h_tZtZ_1000 = getHisto('TT_tZtZ_M1000_', 'TT_tZtZ_M1000', pDir, var,  tZtZ_1000, kBlue-9,    verbose)
#h_tZbW_1000 = getHisto('TT_tZbW_M1000_', 'TT_tZbW_M1000', pDir, var,  tZbW_1000, kOrange-9,  verbose)
h_tZtH_1000 = getHisto('TT_tZtH_M1000_', 'TT_tZtH_M1000', pDir, var,  tZtH_1000, kMagenta+1, verbose)
h_tZtZ_1200 = getHisto('TT_tZtZ_M1200_', 'TT_tZtZ_M1200', pDir, var,  tZtZ_1200, kBlue+4,    verbose)
#h_tZbW_1200 = getHisto('TT_tZbW_M1200_', 'TT_tZbW_M1200', pDir, var,  tZbW_1200, kBlue+3,    verbose)
h_tZtH_1200 = getHisto('TT_tZtH_M1200_', 'TT_tZtH_M1200', pDir, var,  tZtH_1200, kBlue+2,    verbose)





#h_bZbZ_700 = getHisto('BB_bZbZ_M700_', 'BB_bZbZ_M700', pDir, var,  bZbZ_700, kRed,    verbose)
#h_bZtW_700 = getHisto('BB_bZtW_M700_', 'BB_bZtW_M700', pDir, var,  tZbW_700, kCyan+4,  verbose)
#h_bZbH_700 = getHisto('BB_bZbH_M700_', 'BB_bZbH_M700', pDir, var,  bZbH_700, kBlue+4, verbose)
#h_bZbZ_800 = getHisto('BB_bZbZ_M800_', 'BB_bZbZ_M800', pDir, var,  bZbZ_800, kRed+4,    verbose)
#h_bZtW_800 = getHisto('BB_bZtW_M800_', 'BB_bZtW_M800', pDir, var,  bZtW_800, kRed+3,  verbose)
#h_bZbH_800 = getHisto('BB_bZbH_M800_', 'BB_bZbH_M800', pDir, var,  bZbH_800, kRed+2, verbose)
#h_bZbZ_1000 = getHisto('BB_bZbZ_M1000_', 'BB_bZbZ_M1000', pDir, var,  bZbZ_1000, kMagenta+4,    verbose)
#h_bZtW_1000 = getHisto('BB_bZtW_M1000_', 'BB_bZtW_M1000', pDir, var,  bZtW_1000, kOrange+11,    verbose)
#h_bZbH_1000 = getHisto('BB_bZbH_M1000_', 'BB_bZbH_M1000', pDir, var,  bZbH_1000, kCyan+1,    verbose)
#h_bZbZ_1200 = getHisto('BB_bZbZ_M1200_', 'BB_bZbZ_M1200', pDir, var,  bZbZ_1200, kMagenta+4,    verbose)
#h_bZtW_1200 = getHisto('BB_bZtW_M1200_', 'BB_bZtW_M1200', pDir, var,  bZtW_1200, kMagenta+3,    verbose)
#h_bZbH_1200 = getHisto('BB_bZbH_M1200_', 'BB_bZbH_M1200', pDir, var,  bZbH_1200, kMagenta+2,    verbose)




templates = []
#templates.append(h_data)
templates.append(h_dy)
templates.append(h_top)
templates.append(h_diboson)
#templates.append(h_vv)
#templates.append(h_st)
#templates.append(h_wjets)
#templates.append(h_tZtZ_700)
#templates.append(h_tZbW_700)
#templates.append(h_tZtH_700)
templates.append(h_tZtZ_800)
#templates.append(h_tZbW_800)
templates.append(h_tZtH_800)
templates.append(h_tZtZ_1000)
#templates.append(h_tZbW_1000)
templates.append(h_tZtH_1000)
templates.append(h_tZtZ_1200)
#templates.append(h_tZbW_1200)
templates.append(h_tZtH_1200)

#templates.append(h_bZbZ_700)
#templates.append(h_tZbW_700)
#templates.append(h_bZbH_700)
#templates.append(h_bZbZ_800)
#templates.append(h_bZtW_800)
#templates.append(h_bZbH_800)
#templates.append(h_bZbZ_1000)
#templates.append(h_bZtW_1000)
#templates.append(h_bZbH_1000)
#templates.append(h_bZbZ_1200)
#templates.append(h_bZtW_1200)
#templates.append(h_bZbH_1200)

#f = TFile(plotDir+"/"+skimType+"/"+var+".root", "RECREATE")
#for ihist in templates:
#    ihist.Write()
#f.Close()

f = TFile(plotDir+"/"+skimType+"/"+var+".root", "RECREATE")                                                                                               

#if h_data in templates:

 #   h_data.SetName("mumu1b__DATA")
  #  h_data.Write("")
                                                 

if h_dy in templates:
    
    h_dy.SetName("eeT0Z0b2__DY")
    h_dy.Write("")
    
if h_diboson in templates:
        
    h_diboson.SetName("eeT0Z0b2__VV")
    h_diboson.Write("")

if h_top in templates:                                                                                                                                                                        
                                                                                                                                                                                              
    h_top.SetName("eeT0Z0b2__Top")                                                                                                                                                               
    h_top.Write("")    


if h_tZtZ_800 in templates:

    h_tZtZ_800.SetName("eeT0Z0b2__TT_tZtZ_M800")
    h_tZtZ_800.Write("")

if h_tZtH_800 in templates:

    h_tZtH_800.SetName("eeT0Z0b2__TT_tZtH_M800")
    h_tZtH_800.Write("")

#if h_tZbW_800 in templates:

 #   h_tZbW_800.SetName("TT_tZbW_M800")
  #  h_tZbW_800.Write("")

if h_tZtZ_1000 in templates:

    h_tZtZ_1000.SetName("eeT0Z0b2__TT_tZtZ_M1000")
    h_tZtZ_1000.Write("")

if h_tZtH_1000 in templates:

    h_tZtH_1000.SetName("eeT0Z0b2__TT_tZtH_M1000")
    h_tZtH_1000.Write("")

#if h_tZbW_1000 in templates:

 #   h_tZbW_1000.SetName("TT_tZbW_M1000")
  #  h_tZbW_1000.Write("")

if h_tZtZ_1200 in templates:

    h_tZtZ_1200.SetName("eeT0Z0b2__TT_tZtZ_M1200")
    h_tZtZ_1200.Write("")

if h_tZtH_1200 in templates:

    h_tZtH_1200.SetName("eeT0Z0b2__TT_tZtH_M1200")
    h_tZtH_1200.Write("")

#if h_tZbW_1200 in templates:

 #   h_tZbW_1200.SetName("TT_tZbW_M1200")
 #   h_tZbW_1200.Write("")


f.Close()

#get background uncertainty
h_bkg = h_top.Clone()
h_bkg.Reset()
h_bkg.SetName("total bkg")
h_bkg.Add(h_diboson)
h_bkg.Add(h_dy)
h_bkg.Add(h_top)

#h_bkg.Add(h_vv)
#h_bkg.Add(h_wjets)
#h_bkg.Add(h_st)

#histo properties
nBins = h_bkg.GetNbinsX()
bMin = h_bkg.GetBinLowEdge(1)
bMax = h_bkg.GetBinLowEdge(nBins+1)
bin1 = h_bkg.GetXaxis().FindBin(bMin)
bin2 = h_bkg.GetXaxis().FindBin(bMax)

for ibin in range(0,nBins+1):    
    iTop     = h_top.GetBinContent(ibin)
    iDY      = h_dy.GetBinContent(ibin)
    iDiboson = h_diboson.GetBinContent(ibin)
 #   iWJ      = h_wjets.GetBinContent(ibin)
  #  iVV      = h_vv.GetBinContent(ibin)
    # stat error
    stat_err = (h_bkg.GetBinError(ibin))**2 
    # add approximate systematic uncertainty to each bin
    lumi_err = 0.04**2
    btag_err = 0.017**2
    ID_err   = 0.03**2
    JES_err  = 0.05*0.05
    dy_err   = (0.1*iDY)**2
    top_err  = (0.2*iTop)**2
    diboson_err= (0.3*iDiboson)**2
   # st_err   = (0.3*iTop)**2
   # wjet_err = (0.1*iWJ)**2
   # vv_err   = (0.3*iVV)**2

    new_err = stat_err + lumi_err + btag_err + ID_err + JES_err + dy_err + top_err   + diboson_err # + wjet_err +st_err + vv_err
    if h_bkg.GetBinError(ibin) != 0: h_bkg.SetBinError(ibin, TMath.Sqrt(new_err))

h_bkg.SetMarkerSize(0)
h_bkg.SetLineWidth(2)
h_bkg.SetFillColor(14)
h_bkg.SetLineColor(0)
h_bkg.SetFillStyle(3244)

#histogram to print the total background with stat uncertainty
h_tot = h_top.Clone()
h_tot.Reset()
h_tot.SetName("Total_"+h_tot.GetName().split('_',1)[1])
h_tot.Add(h_top)
h_tot.Add(h_dy)
h_tot.Add(h_diboson)
#h_tot.Add(h_wjets)
#h_tot.Add(h_st)

print h_tot.GetName().split('_',1)[1]

## =========Drawing==============
integralError = Double(5)
# print the latex table:
print '\\begin{tabular}{|c|c| }'
print '\hline'
print 'Sample     & Events  \\\\ '
print '\hline'
count = 0
for ihist in templates :
    #if var != 'cutFlow':overUnderFlow(ihist)
    count = count+1
    if count == 10:
        print '\hline'
    if count == 14:
        print '\hline'
    ihist.IntegralAndError(bin1,bin2,integralError)
    if 'TT' in ihist.GetName() or 'BB' in ihist.GetName():
        print '{0:<5} & {1:<5.2f} $\pm$ {2:<5.2f} \\\\ '.format(ihist.GetName().split('_')[1], ihist.Integral(bin1,bin2), integralError)
    else:      
        print '{0:<5} & {1:<5.2f} $\pm$ {2:<5.2f} \\\\ '.format(ihist.GetName().split('_')[0], ihist.Integral(bin1,bin2), integralError)

print '\hline'
print '{0:<5} & {1:<5.0f} \\\\ '.format('Tot Bkg', h_tot.Integral(bin1,bin2), integralError)
print '\hline'
print '{0:<5} & {1:<5.0f} \\\\ '.format(h_data.GetName().split('_')[0], h_data.Integral())
print '\end{tabular}'
print 'bkg : ', h_bkg.Integral(ibin,bin2), 'tot : ', h_tot.Integral(ibin,bin2)

hs = THStack("","")

for ihist in reversed(templates[0:3]):# 5 was here
    hs.Add(ihist)
    print 'histo added', ihist.GetName()

# Canvas
c1 = TCanvas('c1', 'c1', 800, 600)
c1.Divide(1,2)
scale = (1.0 - 0.3)/0.35

# prepare top pad for original plot
pad = c1.cd(1)
pad.SetPad(0, 0.3, 1, 1)
pad.SetTopMargin(0.1)
pad.SetBottomMargin(0.005)
t = pad.GetTopMargin()

# prepare the 2nd pad
pad = c1.cd(2)
pad.SetPad(0, 0.0, 1, 0.3)
pad.SetTopMargin(0.06)
pad.SetBottomMargin(0.4)
pad.SetTickx(1)
pad.SetTicky(1)
c1.cd(1)

hs.SetMaximum(hs.GetMaximum()*5)
hs.SetMinimum(0.1)
gPad.SetLogy()

hs.Draw("Hist")
h_bkg.Draw("e2 same")
h_data.Draw("same")

for ihist in reversed(templates[3:13]): #5:8
    print 'overlaying, ', ihist.GetName() 
    ihist.Draw("ehist same")

xTitle= h_top.GetXaxis().GetTitle()
yTitle= h_top.GetYaxis().GetTitle()

hs.GetXaxis().SetTitle(xTitle)

setTitle(hs, xTitle)
gPad.RedrawAxis()
ll = TLatex()
ll.SetNDC(kTRUE)
ll.SetTextSize(0.05)
ll.DrawLatex(0.78,0.92, "2.2 fb^{-1} (13 TeV)");#2.2 0r 0.53

cms = TLatex()
cms.SetNDC(kTRUE)
cms.SetTextFont(61)
cms.SetTextSize(0.08)
cms.DrawLatex(0.12, 1-t+0.2*t,"CMS")

sel = TLatex()
sel.SetNDC(kTRUE)
sel.SetTextSize(0.065)

chan = TLatex()
chan.SetNDC(kTRUE)
chan.SetTextSize(0.065)
chan.DrawLatex(0.50, 0.76, title)

prel = TLatex()
prel.SetNDC(kTRUE)
prel.SetTextFont(52)
prel.SetTextSize(0.75*t*0.76)
prel.DrawLatex(0.22,0.92,"Preliminary")

leg.Draw()
gPad.RedrawAxis()



c1.cd(2)
# add the systematic band
h_ratio = h_data.Clone()
h_ratio_bkg = h_bkg.Clone()
h_ratio_bkg.SetDirectory(0)
h_ratio.SetDirectory(0)
h_ratio.Divide(h_data, h_tot)
h_ratio_bkg.Divide(h_bkg, h_tot)

for ibin in range(1, nBins+1):
    if h_bkg.GetBinContent(ibin) == 0: h_ratio_bkg.SetBinContent(ibin,1)

prepareRatio(h_ratio, h_ratio_bkg, scale, xTitle)

line = TLine(bMin, 1, bMax, 1)
line.SetLineColor(kBlack)
h_ratio.Draw("")
h_ratio_bkg.Draw("e2same")
h_ratio.Draw("same")
line.Draw()

gPad.RedrawAxis()

#create a directory if it doesn't exist
m_1 = 'mkdir '+plotDir
m_2 = 'mkdir '+plotDir+"/"+skimType
if not os.path.isdir(plotDir):
    subprocess.call( [m_1], shell=True )
if not os.path.isdir(plotDir+"/"+skimType):
    subprocess.call( [m_2], shell=True )    
    
c1.SaveAs(plotDir+"/"+skimType+"/"+var+"_.pdf")
#c1.SaveAs(plotDir+"/"+skimType+"/"+var+"_.gif")
#raw_input("hold on")
