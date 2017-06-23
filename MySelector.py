from ROOT import *; import glob, numpy as n; from array import array
from variables import *
isMC = 0

_fileOUT = 'B_JPsi_Kstar_2171.root'

MyFileNamesMC = glob.glob( MCpath(1) + "*.root")
MyFileNamesDA = glob.glob("/afs/cern.ch/work/o/ofilatov/CMSSW_5_3_24/src/XbFrame/Xb_frame/crab_projects_Bfinder_B_JPsi_Kstar_v1/crab_Bfinder_*/results/*.root")

##__aa = 0;    __bb = 10
__aa = 0;  __bb =  len(MyFileNamesDA);
MyFileNames = (MyFileNamesMC if isMC else MyFileNamesDA[__aa: __bb]); ch = TChain('mkcands/ntuple');

for fName in  MyFileNames:
    ii = ch.Add(fName);

print 'get ', len(MyFileNames), 'files from', __aa,'to',__bb,';  chain created'

fileOUT  = TFile (_fileOUT, "recreate");    mytree = TTree("mytree","mytree");

nEvt = ch.GetEntries(); print "entries: from", 0, 'to', nEvt-1; 
H_cuts = TH1F("H_cuts", "H_cuts", 40, 0, 20)

###  declaration and connecting to the branches of my new variables {{{1
NOUT, NOUT_evt, BBB, ibs = [int(0) for i in range(4)];
PV, PVE, JPV, JPVE, JPP3, BsV_Cjp, BsP3_Cjp, BsVE_Cjp, _TV3zero = [TVector3() for i in range(9)]
JP_Kaon_P4_cjp, B_P4_c0, B_P4_Cjp, kaonP_P4_cjp, kaonP_P4_0c, pion_P4, kaonStar_P4, MU1P4_cjp, MU2P4_cjp = [TLorentzVector() for i in range(9)];
_TV3zero  = TVector3(0,0,0)

_MY_VARS_ = [

'kaonP_pt_0c', 'kaonP_pt_cjp',
'kaonP_track_normchi2', 'kaonP_Hits',  'kaonP_PHits',
'kaonP_dxy_Bsdecay', 'kaonP_dz_Bsdecay', 'kaonP_NTrackerLayers',  'kaonP_NPixelLayers',

#-----~-----
'pfMass', 'pion_pt',
'pi_maxDeltaR', 'pi_numberOfGammas',

'deltaR_K_pi',
#-----~-----

'kStar_mass', 'kStar_pt',
'deltaR_Jpsi_kStar',

#-----~-----
'areSoft', 'areTight_def', 'areTight_HM', 'areMyGlobal',

'mum_relIso', 'mum_NMuonStations', 'mum_dxy_Bsdecay', 'mum_dz_Bsdecay', 
'mum_isGlobalMuon', 'mum_isTrackerMuon', 'mum_isGoodLS_OptimT',

'mup_relIso', 'mup_NMuonStations', 'mup_dxy_Bsdecay', 'mup_dz_Bsdecay', 
'mup_isGlobalMuon', 'mup_isTrackerMuon', 'mup_isGoodLS_OptimT',
'deltaR_mupmum_cjp',

#-----~-----
##"JP_Eta_cjp", "JP_Phi_cjp", "JP_pvcos2_Cmumu",
'JP_mass_c0',
'Jpsi_VtxProb', "JP_DS_2D_Cmumu", 

#-----~-----

'JP_Kaon_mass_cjp',
"B_mass_c0", 'B_mass_cjp',
"B_pt_Cjp", "B_pvdistsignif2_Cjp", 
"B_pvcos2_Cjp", "B_vtxprob_Cjp",
"B_Eta_cjp", "B_Phi_cjp", 

'PV_refit_prob',

"SAMEEVENT"]
_MC_VARS = ["MC_mu", "MC_k1"];
if isMC: _MY_VARS_ += _MC_VARS

for _var_ in _MY_VARS_:
    exec(_var_ + ' = n.zeros(1, dtype=float)')

for _var_ in _MY_VARS_:
    #print 'executing ' + 'mytree.Branch("' + _var_ + '"' + ' '*(25-len(_var_)) + ',' + _var_ + ' '*(25-len(_var_)) + ', "'+ _var_ + '/D")'
    exec('mytree.Branch("' + _var_ + '"' + ' '*(25-len(_var_)) + ',' + _var_ + ' '*(25-len(_var_)) + ', "'+ _var_ + '/D")')

###  declaration and connecting to the branches of my new variables }}}1

match_i = 0
for evt in range(0, nEvt):
    ##
    if (ch.GetEntry(evt) <= 0) : break;
    BInfo_size  = ch.nB
    if len(ch.mum_dxy_Bsdecay) != BInfo_size:
        print 'Sizes do not match!', 'array len = ', len(ch.mum_dxy_Bsdecay), ' nB = ', BInfo_size
	match_i += 1
        continue
    
    for Bj in range(BInfo_size):
        ibs = Bj

#####~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~Muons~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~#####

	# Soft J/psi muons #
	areSoft[0] = 0   if (  ch.mum_NTrackerLayers[ibs] <= 5 or ch.mup_NTrackerLayers[ibs] <= 5 or 
	      			      ch.mum_NPixelLayers[ibs] <= 0 or ch.mup_NPixelLayers[ibs] <= 0 or
	      			      abs(ch.mum_dxy_Bsdecay[ibs]) >= 0.3 or abs(ch.mup_dxy_Bsdecay[ibs]) >= 0.3 or 
				      abs(ch.mum_dz_Bsdecay[ibs]) >= 20. or abs(ch.mup_dz_Bsdecay[ibs]) >= 20. or
	      			      ch.mumAngT[ibs] == 0 or ch.mupAngT[ibs] == 0  ) else 1
	
	# Default Tight J/psi muons #	
	areTight_def[0] = 0   if ( ch.mum_isTight[ibs] <= 0 or ch.mup_isTight[ibs] <= 0) else 1


	# Handmade Tight J/psi muons #
	areTight_HM[0] = 0	  if ( ch.mum_isGlobalMuon[ibs] == 0 or ch.mup_isGlobalMuon[ibs] == 0 or 
	    			       ch.mum_normChi2[ibs] >= 10. or ch.mup_normChi2[ibs] >= 10. or 
	    			       ch.mum_normChi2[ibs] < 0. or ch.mup_normChi2[ibs] < 0. or
	    			       ch.mum_NMuonHits[ibs] <= 0 or ch.mup_NMuonHits[ibs] <= 0 or 
	    			       ch.mum_NMuonStations[ibs] <= 1 or ch.mup_NMuonStations[ibs] <= 1 or
	    			       abs(ch.mum_dxy_Bsdecay[ibs]) >= 0.2 or abs(ch.mup_dxy_Bsdecay[ibs]) >= 0.2 or abs(ch.mum_dz_Bsdecay[ibs]) >= 0.5 or abs(ch.mup_dz_Bsdecay[ibs]) >= 0.5 or
	  			       ch.mumNPHits[ibs] <= 0 or ch.mupNPHits[ibs] <= 0 or
	  			       ch.mum_NTrackerLayers[ibs] <= 5 or ch.mup_NTrackerLayers[ibs] <= 5 or 
				       ch.mum_normChi2[ibs] < 0 or ch.mum_NMuonHits < 0 or
				       ch.mup_normChi2[ibs] < 0 or ch.mup_NMuonHits < 0 ) else 1


	mum_relIso[0] = ch.mum_relIso[ibs]; mup_relIso[0] = ch.mup_relIso[ibs];
	mum_isGlobalMuon[0] = ch.mum_isGlobalMuon[ibs]; mup_isGlobalMuon[0] = ch.mup_isGlobalMuon[ibs];
	mum_NMuonStations[0] = ch.mum_NMuonStations[ibs]; mup_NMuonStations[0] = ch.mup_NMuonStations[ibs]; 


	#   Global muon requirements from CMS AN-2008/098   # 
	#                 (without d0 cut)                  #
	areMyGlobal[0] = 0    if (   ch.mum_isGlobalMuon[ibs] == 0 or ch.mup_isGlobalMuon[ibs] == 0 or 
	     				 ch.mum_normChi2[ibs] >= 10. or ch.mup_normChi2[ibs] >= 10. or 
	     				 ch.mum_normChi2[ibs] < 0. or ch.mup_normChi2[ibs] < 0. or
	     				 ch.mumNHits[ibs] <= 10 or ch.mupNHits[ibs] <= 10)   else 1		



	mum_dxy_Bsdecay[0]   = ch.mum_dxy_Bsdecay[ibs];    mum_dz_Bsdecay[0]      = ch.mum_dz_Bsdecay[ibs]; 
	mum_isTrackerMuon[0] = ch.mum_isTrackerMuon[ibs];  mum_isGoodLS_OptimT[0] = ch.mum_isGoodLS_OptimT[ibs];  

	mup_dxy_Bsdecay[0]   = ch.mup_dxy_Bsdecay[ibs];    mup_dz_Bsdecay[0]      = ch.mup_dz_Bsdecay[ibs];  
	mup_isTrackerMuon[0] = ch.mup_isTrackerMuon[ibs];  mup_isGoodLS_OptimT[0] = ch.mup_isGoodLS_OptimT[ibs];  

	#
##        if (not 'HLT_DoubleMu4_Jpsi_Displaced' in ch.triggersMuPL[ibs]) or (not 'HLT_DoubleMu4_Jpsi_Displaced' in ch.triggersMuML[ibs])  :continue
        #



#####~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~J/psi~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~#####

        MU1P4_cjp   .SetXYZM(ch.B_mu_px1_cjp[ibs], ch.B_mu_py1_cjp[ibs], ch.B_mu_pz1_cjp[ibs], PDG_MUON_MASS)
        MU2P4_cjp   .SetXYZM(ch.B_mu_px2_cjp[ibs], ch.B_mu_py2_cjp[ibs], ch.B_mu_pz2_cjp[ibs], PDG_MUON_MASS)
        MUMUP4_cjp = MU1P4_cjp + MU2P4_cjp
        
        JPV     = TVector3( ch.B_J_DecayVtxX[ibs],  ch.B_J_DecayVtxY[ibs],  ch.B_J_DecayVtxZ[ibs]   )
        JPVE    = TVector3( sqrt(ch.B_J_DecayVtxXE[ibs]), sqrt(ch.B_J_DecayVtxYE[ibs]), sqrt(ch.B_J_DecayVtxZE[ibs])  )
        JPP3    = TVector3( ch.B_J_px[ibs],         ch.B_J_py[ibs],         ch.B_J_pz[ibs])
        PV          = TVector3( ch.PV_bestBang_RF_X[ibs],   ch.PV_bestBang_RF_Y[ibs],   ch.PV_bestBang_RF_Z[ibs]    )
        PVE         = TVector3( sqrt(ch.PV_bestBang_RF_XE[ibs]),  sqrt(ch.PV_bestBang_RF_YE[ibs]),  sqrt(ch.PV_bestBang_RF_ZE[ibs])  )


        PDG_JPSI_MASS       =   3.096916
        if abs(ch.B_J_mass[ibs] - PDG_JPSI_MASS)    > 0.15  :continue

        if MU1P4_cjp.Pt() < 4.0 or MU2P4_cjp.Pt() < 4.0:
            H_cuts.Fill(11)
            continue

        if sqrt (ch.B_J_px[ibs]**2 + ch.B_J_py[ibs]**2) < 6.9:
            H_cuts.Fill(12)
	    continue
	
        if ch.B_J_Prob[ibs] < 0.1:
            H_cuts.Fill(13)
            continue

        if DirectionCos2 ( JPV - PV, JPP3 ) < 0.9:
            H_cuts.Fill(14)
##            continue

        if DetachSignificance2( JPV - PV, PVE, JPVE) < 3.0:
            H_cuts.Fill(15)
            continue
        
        if abs(MUMUP4_cjp.Eta()) > 2.2  :continue

##        JP_Eta_cjp[0] = MUMUP4_cjp.Eta()
##        JP_Phi_cjp[0] = MUMUP4_cjp.Phi()
##        JPSI_masc_Cmumu[0]          = (MU2P4 + MU1P4).M()
##        JPSI_pvcos2_Cmumu_va        = DirectionCos2 ( JPV - PV, JPP3 )

	deltaR_mupmum_cjp[0] = MU1P4_cjp.DeltaR(MU2P4_cjp)
        Jpsi_VtxProb[0]      = ch.B_J_Prob[ibs]
        JP_DS_2D_Cmumu[0]    = DetachSignificance2( JPV - PV, PVE, JPVE)
        JP_mass_c0[0]           = ch.B_J_mass[ibs]

#####~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~Kaon~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~#####

        kaonP_P4_0c  .SetXYZM(ch.kaonP_px_0c[ibs], ch.kaonP_py_0c[ibs], ch.kaonP_pz_0c[ibs], PDG_KAON_MASS)
        kaonP_P4_cjp .SetXYZM(ch.kaonP_px_cjp[ibs], ch.kaonP_py_cjp[ibs], ch.kaonP_pz_cjp[ibs], PDG_KAON_MASS)
        pion_P4      .SetXYZM(ch.pion_px[ibs], ch.pion_py[ibs], ch.pion_pz[ibs], ch.pfMass[ibs])

	kaonP_pt_cjp[0] = kaonP_P4_cjp.Pt()
	kaonP_pt_0c[0] = kaonP_P4_0c.Pt()

	if kaonP_pt_0c[0] < 0.7 or abs(kaonP_P4_0c.Eta()) > 2.5 :continue


        kaonP_track_normchi2[0] = ch.kaonP_track_normchi2[ibs]
        kaonP_Hits[0] = ch.kaonP_Hits[ibs]
        kaonP_PHits[0] = ch.kaonP_PHits[ibs]
        kaonP_NTrackerLayers[0] = ch.kaonP_NTrackerLayers[ibs]
        kaonP_NPixelLayers[0] = ch.kaonP_NPixelLayers[ibs]
        kaonP_dxy_Bsdecay[0] = ch.kaonP_dxy_Bsdecay[ibs]
        kaonP_dz_Bsdecay[0] = ch.kaonP_dz_Bsdecay[ibs]
        

#####~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~Pi~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~#####
	
	pfMass[0] = ch.pfMass[ibs]

	pion_pt[0] = ch.pion_pt[ibs]
	pi_maxDeltaR[0] = ch.pi_maxDeltaR[ibs]
	pi_numberOfGammas[0] = ch.pi_numberOfGammas[ibs]
	
	deltaR_K_pi[0] = kaonP_P4_cjp.DeltaR(pion_P4)


#####~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~Kstar~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~#####
	
        kaonStar_P4 = pion_P4 + kaonP_P4_cjp
        kStar_mass[0] = kaonStar_P4.M()
        kStar_pt = kaonStar_P4.Pt()

        deltaR_Jpsi_kStar[0] = kaonStar_P4.DeltaR(MUMUP4_cjp)
        	

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~B and misc.~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

        JP_Kaon_P4_cjp = MUMUP4_cjp + kaonP_P4_cjp
        B_P4_Cjp     = pion_P4 + JP_Kaon_P4_cjp        
        
        BsV_Cjp     = TVector3(ch.B_DecayVtxX[ibs],  ch.B_DecayVtxY[ibs],  ch.B_DecayVtxZ[ibs]   )
        BsVE_Cjp    = TVector3( sqrt(ch.B_DecayVtxXE[ibs]), sqrt(ch.B_DecayVtxYE[ibs]), sqrt(ch.B_DecayVtxZE[ibs])  )
        BsP3_Cjp    = B_P4_Cjp.Vect()

        B_mass_c0[0]          = ch.Bs_mass_c0[ibs]
        B_mass_cjp[0]          = B_P4_Cjp.M()
        JP_Kaon_mass_cjp[0] = ch.JPsi_K_mass_cjp[ibs]
##        PDG_BU_MASS         = 5.27932
##        if abs(PDG_BU_MASS - B_mass_c0[0])  > 0.05   :continue # were 5.32 - 5.41


        B_pt_Cjp[0]            = B_P4_Cjp.Pt()

        B_pvdistsignif2_Cjp[0] = DetachSignificance2( BsV_Cjp - PV, PVE, BsVE_Cjp)
        if B_pvdistsignif2_Cjp[0] < 3. :continue

        B_pvcos2_Cjp[0]        = DirectionCos2 ( BsV_Cjp - PV, BsP3_Cjp )  ## note that the PV was selected for 'bare' B momentum vector, w/o cjp components
                                                                           ## here you are trying to calculate cosine for cjp B meson vector
        if B_pvcos2_Cjp[0] < 0.9 :
            H_cuts.Fill(9)
##            continue

        B_vtxprob_Cjp[0]       = ch.B_Prob[ibs]
        if B_vtxprob_Cjp[0] < 0.1 :  ## 0.05 in Bfinder
            H_cuts.Fill(10)
            continue

        if abs(B_P4_Cjp.Eta())  > 2.5   :continue

        B_Phi_cjp[0]            = B_P4_Cjp.Phi()
        B_Eta_cjp[0]            = B_P4_Cjp.Eta()

	PV_refit_prob[0] = ch.PV_bestBang_RF_CL[ibs]
#---------------------------------------------------

        _mctr = 0
        if isMC:
            _mctr   =   1 if abs(ch.MCID_k1[ibs])==321      else 0;
            _mctr   +=  2 if abs(ch.MCID_pk1[ibs])==521     else 0;
            MC_k1[0] =  _mctr
            _mctr   =   1 if (abs(ch.MCID_mu1[ibs])==13     and abs(ch.MCID_mu2[ibs])==13)      else 0;
            _mctr   +=  2 if (abs(ch.MCID_pmu1[ibs])==443   and abs(ch.MCID_pmu2[ibs])==443)    else 0;
            _mctr   +=  4 if (abs(ch.MCID_ppmu1[ibs])==521  and abs(ch.MCID_ppmu2[ibs])==521)   else 0;
            MC_mu[0] = _mctr
        #
        SAMEEVENT[0] = 0;
        if (BBB > -1) : 
            SAMEEVENT[0] = 1
            NOUT_evt -= 1;
        #
        mytree.Fill(); NOUT += 1; NOUT_evt +=1; BBB = Bj; 
   
    BBB = -1
    if (evt % 2000 == 0) :    ## printout progress
        _perc = str(TMath.Nint(100*(evt+1)/(nEvt+0.0)));
        print "["+_perc+(' ' * (3 - len(_perc)))+"%];evt",evt,' '*(6-len(str(evt))),";saved ["+str(__aa)+":"+str(__bb)+"]", NOUT, ' ', NOUT_evt

fileOUT.Write();
print NOUT, ' ', NOUT_evt
print 'nEvt = ', nEvt, '; match_i = ', match_i
