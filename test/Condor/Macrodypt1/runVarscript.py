#!/bin/python

import subprocess

dir = 'cat'



if dir == 'cat':
    options = [

         # ['cutflow4'],
         # ['cutflow1'],
          #['cutflow2'],
          #['cutflow3'],
          ['cutflow'],

       # ['ST_sig'],       
        # ['ST_sig1b'],
        # ['ST_sig2b'],
        # ['ST_sigT1Z1'],
        # ['ST_sigT0Z1'],
        # ['ST_sigT1Z0'],
        # ['ST_sigT0Z0'],
        

         # ['ST_sigT1Z1H1'],
         # ['ST_sigT1Z1H0'],
         # ['ST_sigT0Z1H1'],
         # ['ST_sigT0Z1H0'],
         # ['ST_sigT1Z0H1'],
         # ['ST_sigT1Z0H0'],
         # ['ST_sigT0Z0H1'],
         # ['ST_sigT0Z0H0'],
        ['st'],
         ['ST_sigT1Z1H1b1'],
         ['ST_sigT1Z1H1b2'],
         ['ST_sigT1Z1H0b1'],
         ['ST_sigT1Z1H0b2'],
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

        
       # ['Hmass-boosted-cnt'],
       # ['HPt-boosted-cnt'],
       # ['nHcandidate-boosted-cnt'],
       # ['Hmassnb-cnt'],
       # ['HPtnb-cnt'],
       # ['nHcandidatesnb-cnt'],
       # ['nHcandidates-tot-cnt'],                                                                                                                                             
       # ['nHcandidates1-tot-cnt'],
        
        
        
        
       # ['Zmass-boosted-cnt'],
       # ['ZPt-boosted-cnt'],
       # ['nzcandidate-boosted-cnt'],
        
       # ['Zmass-cnt'],
       # ['ZPt-cnt'],
       # ['nzcandidates-cnt'],
       # ['nzcandidates-tot-cnt'],
       # ['nzcandidates1-tot-cnt'], 
   
       # ['topmass-D-cnt'],
       # ['topPt-D-cnt'],
       # ['ntopcandidate-D-cnt'],
        
        
       # ['Wmass-BC-cnt'],
       # ['nWcandidate-BC-cnt'],
       # ['lightjetmass-BC-cnt'],
       # ['nlightjetcandidate-cnt'],
        
       # ['topmas-A-cnt'],
       # ['topPt-A-cnt'],
       # ['ntopcandidate-A-cnt'],
        
       # ['topmass-Bc-cnt'],
       # ['topPt-BC-cnt'],
       # ['ntopcandidate-BC-cnt'],
        
       # ['ntopcandidate-tot-cnt'],
       # ['ntopcandidate1-tot-cnt'],
        
        
       # ['Hmass-boosted-sig'],
       # ['HPt-boosted-sig'],
       # ['nHcandidate-boosted-sig'],
       # ['Hmassnb-sig'],
       # ['HPtnb-sig'],
       # ['nHcandidatesnb-sig'],
       # ['nHcandidates-tot-sig'],                                                                                                                                       

        #['nHcandidates1-tot-sig'],




        #['Zmass-boosted-sig'],
        #['ZPt-boosted-sig'],
        #['nzcandidate-boosted-sig'],

        #['Zmass-sig'],
        #['ZPt-sig'],
        #['nzcandidates-sig'],
        #['nzcandidates-tot-sig'],
        #['nzcandidates1-tot-sig'],

        #['topmass-D-sig'],
        #['topPt-D-sig'],
        #['ntopcandidate-D-sig'],


        #['Wmass-BC-sig'],
        #['nWcandidate-BC-sig'],
        #['lightjetmass-BC-sig'],
        #['nlightjetcandidate-sig'],

        #['topmas-A-sig'],
        #['topPt-A-sig'],
        #['ntopcandidate-A-sig'],

        #['topmass-Bc-sig'],
        #['topPt-BC-sig'],
        #['ntopcandidate-BC-sig'],

        #['ntopcandidate-tot-sig'],
        #['ntopcandidate1-tot-sig'],

             

        ]


command = 'python plotscript.py --var={0:s}'

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
