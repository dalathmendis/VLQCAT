#!/bin/python

import subprocess

dir = 'cnt'

if dir == 'pre':
    options = [

       # ['Iso04_mu_pre'],
       # ['Dz_mu_pre'],
       # ['Dxy_mu_pre'],
       # ['IsPFMuon_mu_pre'],
       # ['IsGlobalMuon_mu_pre'],
       # ['GlbTrkNormChi2_mu_pre'],
       # ['NumberValidMuonHits_mu_pre'],
       # ['NumberMatchedStations_mu_pre'],
       # ['NumberValidPixelHits_mu_pre'],
       # ['NumberTrackerLayers_mu_pre'],



        ['massak4jet1_pre'],
        ['massak4jet2_pre'],       
        ['cutflow'],
        ['nob_st'],

        ['npv_pre'],
        ['nak4_pre'],
        ['ht_pre'],
        ['st_pre'],
        ['met_pre'],
        ['met1_pre'],
        ['metPhi_pre'],
        ['ptak4jet1_pre'],
        ['etaak4jet1_pre'],
      #  ['cvsak4jet1_pre'],
        ['ptak4jet2_pre'],
        ['etaak4jet2_pre'],
       # ['cvsak4jet2_pre'],
      #  ['ptak4jet3_pre'],
      #  ['etaak4jet3_pre'],
      #  ['cvsak4jet3_pre'],
        ['phi_jet1MET_pre'],
       
        ['mass_zmumu_pre'],
        ['mass_Zmumu_pre'],
        ['dr_mumu_pre'],
        ['pt_zmumu_pre'],
        ['pt_mu1_pre'],
        ['eta_mu1_pre'],
        ['pt_mu2_pre'],
        ['eta_mu2_pre'],

        ['mass_zelel_pre'],
        ['mass_Zelel_pre'],
        ['dr_elel_pre'],
        ['pt_zelel_pre'],
        ['pt_el1_pre'],
        ['eta_el1_pre'],
        ['pt_el2_pre'],
        ['eta_el2_pre'],
       
        ['Eta_EB_el_pre'],
        ['Eta_EE_el_pre'],
        ['scEta_EB_el_pre'],
        ['scEta_EE_el_pre'],
        ['Iso03_EB_el_pre'],
        ['Iso03_EE_el_pre'],
        ['scEta_EB_el_pre'],
        ['scEta_EE_el_pre'],
        ['dEtaInSeed_EB_el_pre'],
        ['dEtaInSeed_EE_el_pre'],
        ['dPhiIn_EB_el_pre'],
        ['dPhiIn_EE_el_pre'],
        ['Dz_EB_el_pre'],
        ['Dz_EE_el_pre'],
        ['Dxy_EB_el_pre'],
        ['Dxy_EE_el_pre'],
        ['Full5x5siee_EB_el_pre'],
        ['Full5x5siee_EE_el_pre'],
        ['HoE_EB_el_pre'],
        ['HoE_EE_el_pre'],
        ['ooEmooP_EB_el_pre'],
        ['ooEmooP_EE_el_pre'],
        ['missHits_EB_el_pre'],
        ['missHits_EE_el_pre'],
        ['conveto_EB_El_pre'],
        ['conveto_EE_el_pre'],

        ['Iso04_mu_pre'],
        ['Dz_mu_pre'],
        ['Dxy_mu_pre'],
        ['IsPFMuon_mu_pre'],
        ['IsGlobalMuon_mu_pre'],
        ['GlbTrkNormChi2_mu_pre'],
        ['NumberValidMuonHits_mu_pre'],
        ['NumberMatchedStations_mu_pre'],
        ['NumberValidPixelHits_mu_pre'],
        ['NumberTrackerLayers_mu_pre'],

        ]

elif dir == 'cnt':
    options = [      
        ['massak4jet1_cnt'],
        ['massak4jet2_cnt'], 
        ['nob_ht'],
        ['npv_cnt'],
        ['nak4_cnt'],
        ['ht_cnt'],
        ['st_cnt'],
        ['met_cnt'],
        ['met1_cnt'],
        ['metPhi_cnt'],
        ['ptak4jet1_cnt'],
        ['etaak4jet1_cnt'],
        #['cvsak4jet1_cnt'],
        ['ptak4jet2_cnt'],
        ['etaak4jet2_cnt'],
        #['cvsak4jet2_cnt'],
        #['ptak4jet3_cnt'],
        #['etaak4jet3_cnt'],
        #['cvsak4jet3_cnt'],
        ['phi_jet1MET_cnt'],
       
        ['mass_zmumu_cnt'],
        ['mass_Zmumu_cnt'],
        ['dr_mumu_cnt'],
        ['pt_zmumu_cnt'],
        ['pt_mu1_cnt'],
        ['eta_mu1_cnt'],
        ['pt_mu2_cnt'],
        ['eta_mu2_cnt'],
        ['ht1_cnt'],
        ['st1_cnt'],
       # ['nob_ht'],
        ['b_pt_zmumu'],        
        ['b_st'],



        ['mass_zelel_cnt'],
        ['mass_Zelel_cnt'],
        ['dr_elel_cnt'],
        ['pt_zelel_cnt'],
        ['pt_el1_cnt'],
        ['eta_el1_cnt'],
        ['pt_el2_cnt'],
        ['eta_el2_cnt'],            
        ['b_pt_elel'],
        ]

elif dir == 'sig1':
    options = [
        ['1b_ht'],
        ['1b_st'],
        ['nbjets_sig'],
        ['nbjets_cnt'],
        ['nbjets_0btagcnt'],
        ['nbjets_1btagcnt'],
        ['ht_0btagcnt'],
        #['1b_ht'],
        ['ht_1btagcnt'],
        ['lowmet_ht'],
        ['st_0btagcnt'],
        #['1b_st'],
        ['st_1btagcnt'],
        ['lowmet_st'],
        ]

elif dir == 'sig':
    options = [
        ['massak4jet1'],
        ['massak4jet2'],
        ['mass_Zmumu'],
        ['mass_Zelel'],
        ['npv'],
        ['nak4'],
        ['ht'],
        ['st'],
        ['met'],
        ['met1'],
        ['metPhi'],
        ['ptak4jet1'],
        ['etaak4jet1'],
        ['ptak4jet2'],
        ['etaak4jet2'],
        ['phi_jet1MET'],

       # ['mass_zmumu'],
       # ['mass_Zmumu'],
        ['dr_mumu'],
        ['pt_zmumu'],
        ['pt_mu1'],
        ['eta_mu1'],
        ['pt_mu2'],
        ['eta_mu2'],
      

       # ['mass_zelel'],
       # ['mass_Zelel'],
        ['dr_elel'],
        ['pt_zelel'],
        ['pt_el1'],
        ['eta_el1'],
        ['pt_el2'],
        ['eta_el2'],


        ['nbjets'],
        ['ptbjetleading'],
        ['etabjetleading'],
        ['nak8'],
        ['nwjet'],
        ['nhjet'],
        ['ntjet'],
        ['ptak8leading'],
        ['etaak8leading'],
        ['mak8leading'],
        ['prunedmak8leading'],
        ['trimmedmak8leading'],
        ['softdropmak8leading'],
        ['ptak82nd'],
        ['etamak82nd'],
        ['mak82nd'],
        ['purnedmak82nd'],
        ['trimmedmak82nd'],
        ['softdropmak82nd'],
       # ['ptTprime'],
       # ['yTprime'],
       # ['mTprime'],
       # ['ptBprime'],
       # ['yBprime'],
       # ['mBprime'],
       # ['ZJetMasslep'],
       # ['chi2_chi'],
       # ['sqrtChi2'],
       # ['chi_mass'],

        ]
elif dir == 'cat':
    options = [
        #['ST_sig'],
        #['ST_sigT1Z1H1b1'],
        # ['ST_sigT1Z1H1b2'],
        # ['catC'],


       # ['cutflow1'],
       # ['cutflow2'],
       # ['cutflow3'],
        ['cutflow4'],

        # ['ST_sig'],       
        # ['ST_sig1b'],
        # ['ST_sig2b'],
        # ['ST_sigT1Z1'],
        # ['ST_sigT0Z1'],
        # ['ST_sigT1Z0'],
        # ['ST_sigT0Z0'],
        

        #  ['ST_sigT1Z1H1'],
        #  ['ST_sigT1Z1H0'],
        #  ['ST_sigT0Z1H1'],
        #  ['ST_sigT0Z1H0'],
        #  ['ST_sigT1Z0H1'],
        #  ['ST_sigT1Z0H0'],
        #  ['ST_sigT0Z0H1'],
        #  ['ST_sigT0Z0H0'],
        
        # ['ST_sigT1Z1H1b1'],
        # ['ST_sigT1Z1H1b2'],
        # ['ST_sigT1Z1H0b1'],
        # ['ST_sigT1Z1H0b2'],
         ['ST_sigT0Z1H1b1'],
         ['ST_sigT0Z1H1b2'],
         ['ST_sigT0Z1H0b1'],
         ['ST_sigT0Z1H0b2'],
         ['ST_sigT1Z0H1b1'],
         ['ST_sigT1Z0H1b2'],
         ['ST_sigT1Z0H0b1'],
         ['ST_sigT1Z0H0b2'],
         ['ST_sigT0Z0H1b1'],
         ['ST_sigT0Z0H1b2'],
         ['ST_sigT0Z0H0b1'],
         ['ST_sigT0Z0H0b2'],
        
        
        # ['ST_sigT1Z1b1'],
        # ['ST_sigT1Z1b2'],
        # ['ST_sigT0Z1b1'],
        # ['ST_sigT0Z1b2'],
        # ['ST_sigT1Z0b1'],
        # ['ST_sigT1Z0b2'],
        # ['ST_sigT0Z0b1'],
        # ['ST_sigT0Z0b2'],
        #  ['ZHmass-boosted'],                                                                                                                                           
        # ['ZHPt-boosted'],                                                                                                                                             
        # ['nZHcandidate-boosted'],                                                                                                                                     
        # ['ZHmassnb'],                                                                                                                                                 
        # ['ZHPtnb'],                                                                                                                                                   
        # ['nZHcandidatesnb'],                                                                                                                                          
        #['nZHcandidates-tot'],  
        # ['nZHcandidates1-tot'],

        
        ['Hmass-boosted-cnt'],
        ['HPt-boosted-cnt'],
        ['nHcandidate-boosted-cnt'],
        ['Hmassnb-cnt'],
        ['HPtnb-cnt'],
        ['nHcandidatesnb-cnt'],
        ['nHcandidates-tot-cnt'],                                                                                                                                             
        ['nHcandidates1-tot-cnt'],
        
        
        
        
        ['Zmass-boosted-cnt'],
        ['ZPt-boosted-cnt'],
        ['nzcandidate-boosted-cnt'],
        
        ['Zmass-cnt'],
        ['ZPt-cnt'],
        ['nzcandidates-cnt'],
        ['nzcandidates-tot-cnt'],
        ['nzcandidates1-tot-cnt'], 
   
        ['topmass-D-cnt'],
        ['topPt-D-cnt'],
        ['ntopcandidate-D-cnt'],
        
        
        ['Wmass-BC-cnt'],
        ['nWcandidate-BC-cnt'],
        ['lightjetmass-BC-cnt'],
        ['nlightjetcandidate-cnt'],
        
        ['topmas-A-cnt'],
        ['topPt-A-cnt'],
        ['ntopcandidate-A-cnt'],
        
        ['topmass-Bc-cnt'],
        ['topPt-BC-cnt'],
        ['ntopcandidate-BC-cnt'],
        
        ['ntopcandidate-tot-cnt'],
        ['ntopcandidate1-tot-cnt'],
        
        
        ['Hmass-boosted-sig'],
        ['HPt-boosted-sig'],
        ['nHcandidate-boosted-sig'],
        ['Hmassnb-sig'],
        ['HPtnb-sig'],
        ['nHcandidatesnb-sig'],
        ['nHcandidates-tot-sig'],                                                                                                                                       

        ['nHcandidates1-tot-sig'],




        ['Zmass-boosted-sig'],
        ['ZPt-boosted-sig'],
        ['nzcandidate-boosted-sig'],

        ['Zmass-sig'],
        ['ZPt-sig'],
        ['nzcandidates-sig'],
        ['nzcandidates-tot-sig'],
        ['nzcandidates1-tot-sig'],

        ['topmass-D-sig'],
        ['topPt-D-sig'],
        ['ntopcandidate-D-sig'],


        ['Wmass-BC-sig'],
        ['nWcandidate-BC-sig'],
        ['lightjetmass-BC-sig'],
        ['nlightjetcandidate-sig'],

        ['topmas-A-sig'],
        ['topPt-A-sig'],
        ['ntopcandidate-A-sig'],

        ['topmass-Bc-sig'],
        ['topPt-BC-sig'],
        ['ntopcandidate-BC-sig'],

        ['ntopcandidate-tot-sig'],
        ['ntopcandidate1-tot-sig'],

             

        ]

elif dir == 'gen':
    options = [
        ['mass_zelel_pre'],
        ['mass_Zelel_pre'],
        ['dr_elel_pre'],
        ['pt_zelel_pre'],
        ['pt_el1_pre'],
        ['eta_el1_pre'],
        ['pt_el2_pre'],
        ['eta_el2_pre'],

        ['mass_zelel'],
        ['mass_Zelel'],
        ['dr_elel'],
        ['pt_zelel'],
        ['pt_el1'],
        ['eta_el1'],
        ['pt_el2'],
        ['eta_el2'],


        ['htdr_pre'],
        ['stdr_pre'],
        ['drdr_elel_pre'],
        ['htdr_cnt'],
        ['stdr_cnt'],
        ['drdr_elel_cnt'],
        ['htdr'],
        ['stdr'],
        ['drdr_elel'],
        ['st'],
        ]
command = 'python plot.py --var={0:s}'

for option in options :
    s = command.format(
        option[0]
        )

    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo %s"%s,""]                                                                      , shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( [s, ""]                                                                               , shell=True)
