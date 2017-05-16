from ROOT import *; import glob, numpy as n; from array import array
from math import sqrt
from variables import *
isMC = 0

_fileOUT = 'B0_parked_D_notall.root'

MyFileNamesMC = glob.glob( MCpath(1) + "*.root")
MyFileNamesDA = (glob.glob("/afs/cern.ch/work/o/ofilatov/CMSSW_5_3_24/src/XbFrame/Xb_frame/crab_projects_Bfinder_B0_parked/crab_Bfinder_B0_parked_parked_D/results/*.root"))

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
PV, PVE, JPV, JPVE, JPP3, B0V, B0VE, B0P3, _TV3zero = [TVector3(0,0,0) for i in range(9)]
pionP4_0c, pionP4_wo_mass, MU1P4_cjp, MU2P4_cjp, MUMUP4_C0, MUMUP4_cjp, B0P4 = [TLorentzVector(0,0,0,0) for i in range(7)];
_TV3zero  = TVector3(0,0,0)

_MY_VARS_ = [

#-----~-----
'areSoft', 'areTight_def', 'areTight_HM', 'areMyGlobal',
'mu1_pt_cjp','mu2_pt_cjp',
'mu1_pt_0c','mu2_pt_0c',

'mum_relIso', 'mum_NMuonStations', 'mum_dxy_Bsdecay', 'mum_dz_Bsdecay',
'mum_isGlobalMuon', 'mum_isTrackerMuon',
'mum_LastStationOptimizedLowPtT', 'mum_LastStationT', 'mum_OneStationT',
'mum_LastStationAngT', 'mum_OneStationAngT', 'mum_2DCompatibilityT',

'mup_relIso', 'mup_NMuonStations', 'mup_dxy_Bsdecay', 'mup_dz_Bsdecay',
'mup_isGlobalMuon', 'mup_isTrackerMuon',
'mup_LastStationOptimizedLowPtT', 'mup_LastStationT', 'mup_OneStationT',
'mup_LastStationAngT', 'mup_OneStationAngT', 'mup_2DCompatibilityT',

'deltaR_mupmum_cjp',

#-----~-----
"JP_Eta_C0", "JP_Phi_C0",
'JPsi_mass_C0', 'JPsi_pt_C0', 
"JPsi_pvcos2_C0", "JPsi_pv_DS_2D_C0", 'JPsi_VtxProb', 

###-----~-----
##"Bs_mass_Cjp",
##"Bs_pt_Cjp", "Bs_bcvtxDS2d_Cjp", "Bs_bcvtxDS2d_vtxfit", 'Bs_pvDS2d_Cjp', 'Bs_pv_detach_2D',
##"Bs_bcvtx_cos2_Cjp", 'Bs_bcvtx_cos2_vtxfit', 'Bs_pv_cos2_Cjp',
##"Bs_Eta_cjp", "Bs_Phi_cjp",
##
##'Bs_Bcdecay_weight', "Bs_vtxprob_Cjp", 'BsVtx_Chi2',
#-----~-----
'pion_pt_0c', 'pion_mass',
'pion_track_normchi2', 'pion_Hits',  'pion_PHits',
'pion_dxy_Bcdecay', 'pion_dz_Bcdecay', 'pion_NTrackerLayers',  'pion_NPixelLayers',
'pion_Bcdecay_weight',
'deltaR_pi_JP',
'mva_gamma_nh', 'mva_nothing_nh',
#-----~-----
"B0_mass",
"B0_pt", "B0_pvDS2d", "B0_pvDS2d_vtxfit", 'B0_pv_detach_2D', 'B0_pv_detach_2D_vtxfit',
"B0_pvcos2", "B0_pvcos2_vtxfit",
"B0_Eta", "B0_Phi",

'B0Vertex_isValid', "B0_vtxprob", 'B0Vtx_Chi2_kinfit', 'B0Vtx_Chi2_vtxfit', 'B0Vtx_normChi2_vtxfit',
'PV_refit_prob', 'PV_bestBang_RF_X',

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
    if len(ch.B_mu_px1_cjp) != BInfo_size:
##        print 'Sizes do not match!', 'array len = ', len(ch.mum_dxy_Bsdecay), ' nB = ', BInfo_size
        match_i += 1
        continue

    for Bj in range(BInfo_size):
        ##
        ibs = Bj
        ##


    #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
    ###~~~~~~~~~~Vectors and Vertices initialization~~~~~~~~~~###
    #####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

        MU1P4_cjp   .SetXYZM(ch.B_mu_px1_cjp[ibs], ch.B_mu_py1_cjp[ibs], ch.B_mu_pz1_cjp[ibs], PDG_MUON_MASS)
        MU2P4_cjp   .SetXYZM(ch.B_mu_px2_cjp[ibs], ch.B_mu_py2_cjp[ibs], ch.B_mu_pz2_cjp[ibs], PDG_MUON_MASS)
        MU1P4_0c = MU1P4_cjp
        MU2P4_0c = MU2P4_cjp
        MUMUP4_C0    .SetXYZM(ch.B_J_px[ibs], ch.B_J_py[ibs], ch.B_J_pz[ibs], ch.B_J_mass[ibs])
        MUMUP4_cjp = MU1P4_cjp + MU2P4_cjp

        pionP4_wo_mass.SetPtEtaPhiE(ch.pionPF_Pt[ibs], ch.pionPF_Eta[ibs], ch.pionPF_Phi[ibs], ch.pionPF_E[ibs])
        pionP4_0c  .SetPtEtaPhiM(ch.pionPF_Pt[ibs], ch.pionPF_Eta[ibs], ch.pionPF_Phi[ibs], PDG_PI0_MASS)

##        BsP4_Cjp    .SetXYZM  ( ch.Bs_px_cjp[ibs], ch.Bs_py_cjp[ibs], ch.Bs_pz_cjp[ibs], ch.Bs_mass_cjp[ibs])
##        B0P4        .SetXYZM  ( ch.B0_px[ibs], ch.B0_py[ibs], ch.B0_pz[ibs], ch.B0_mass[ibs])
        B0P4 = MUMUP4_C0 + pionP4_0c
        
        PV          = TVector3( ch.PV_bestBang_RF_X[ibs],   ch.PV_bestBang_RF_Y[ibs],   ch.PV_bestBang_RF_Z[ibs]    )
        PVE         = TVector3( sqrt(ch.PV_bestBang_RF_XE[ibs]),  sqrt(ch.PV_bestBang_RF_YE[ibs]),  sqrt(ch.PV_bestBang_RF_ZE[ibs])  )

        JPV     = TVector3( ch.B_J_DecayVtxX[ibs],  ch.B_J_DecayVtxY[ibs],  ch.B_J_DecayVtxZ[ibs]   )
        JPVE    = TVector3( sqrt(ch.B_J_DecayVtxXE[ibs]), sqrt(ch.B_J_DecayVtxYE[ibs]), sqrt(ch.B_J_DecayVtxZE[ibs])  )
        JPP3    = TVector3( ch.B_J_px[ibs],         ch.B_J_py[ibs],         ch.B_J_pz[ibs])

##        BsV_Cjp     = TVector3(ch.Bs_DecayVtxX[ibs],  ch.Bs_DecayVtxY[ibs],  ch.Bs_DecayVtxZ[ibs]   )
##        BsVE_Cjp    = TVector3( sqrt(ch.Bs_DecayVtxXE[ibs]), sqrt(ch.Bs_DecayVtxYE[ibs]), sqrt(ch.Bs_DecayVtxZE[ibs])  )
##        BsP3_Cjp    = BsP4_Cjp.Vect()
##
##        B0V    = TVector3(ch.B0_DecayVtxX[ibs],  ch.B0_DecayVtxY[ibs],  ch.B0_DecayVtxZ[ibs]   )
        B0V = JPV
        B0VE = JPVE
##        B0VE   = TVector3( sqrt(ch.B0_DecayVtxXE[ibs]), sqrt(ch.B0_DecayVtxYE[ibs]), sqrt(ch.B0_DecayVtxZE[ibs])  )
##        B0V_vtxfit     = TVector3(ch.B0_DecayVtx_vtxfit_X[ibs],  ch.B0_DecayVtx_vtxfit_Y[ibs],  ch.B0_DecayVtx_vtxfit_Z[ibs]   )
##        B0VE_vtxfit    = TVector3( sqrt(ch.B0_DecayVtx_vtxfit_XE[ibs]), sqrt(ch.B0_DecayVtx_vtxfit_YE[ibs]), sqrt(ch.B0_DecayVtx_vtxfit_ZE[ibs])  )
        B0P3   = B0P4.Vect()



    #####~~~~~~~~~~~~~~~~~~~~~#####
    ###~~~~~~~~~~Muons~~~~~~~~~~###
    #####~~~~~~~~~~~~~~~~~~~~~#####

##	# Soft J/psi muons #
##    	areSoft[0] = 0   if (  ch.mum_NTrackerLayers[ibs] <= 5 or ch.mup_NTrackerLayers[ibs] <= 5 or
##    	      			      ch.mum_NPixelLayers[ibs] <= 0 or ch.mup_NPixelLayers[ibs] <= 0 or
##    	      			      abs(ch.mum_dxy_Bsdecay[ibs]) >= 0.3 or abs(ch.mup_dxy_Bsdecay[ibs]) >= 0.3 or
##    				      abs(ch.mum_dz_Bsdecay[ibs]) >= 20. or abs(ch.mup_dz_Bsdecay[ibs]) >= 20. or
##    	      			      ch.mumAngT[ibs] == 0 or ch.mupAngT[ibs] == 0  ) else 1
##
##    	# Default Tight J/psi muons #
##    	areTight_def[0] = 0   if ( ch.mum_isTight[ibs] <= 0 or ch.mup_isTight[ibs] <= 0) else 1
##
##
##    	# Handmade Tight J/psi muons #
##    	areTight_HM[0] = 0	  if ( ch.mum_isGlobalMuon[ibs] == 0 or ch.mup_isGlobalMuon[ibs] == 0 or
##    	    			       ch.mum_normChi2[ibs] >= 10. or ch.mup_normChi2[ibs] >= 10. or
##    	    			       ch.mum_normChi2[ibs] < 0. or ch.mup_normChi2[ibs] < 0. or
##    	    			       ch.mum_NMuonHits[ibs] <= 0 or ch.mup_NMuonHits[ibs] <= 0 or
##    	    			       ch.mum_NMuonStations[ibs] <= 1 or ch.mup_NMuonStations[ibs] <= 1 or
##    	    			       abs(ch.mum_dxy_Bsdecay[ibs]) >= 0.2 or abs(ch.mup_dxy_Bsdecay[ibs]) >= 0.2 or abs(ch.mum_dz_Bsdecay[ibs]) >= 0.5 or abs(ch.mup_dz_Bsdecay[ibs]) >= 0.5 or
##    	  			       ch.mumNPHits[ibs] <= 0 or ch.mupNPHits[ibs] <= 0 or
##    	  			       ch.mum_NTrackerLayers[ibs] <= 5 or ch.mup_NTrackerLayers[ibs] <= 5 or
##    				       ch.mum_normChi2[ibs] < 0 or ch.mum_NMuonHits[ibs] < 0 or
##    				       ch.mup_normChi2[ibs] < 0 or ch.mup_NMuonHits[ibs] < 0 ) else 1
##

    	mum_relIso[0] = ch.mum_relIso[ibs]; mup_relIso[0] = ch.mup_relIso[ibs];
    	mum_isGlobalMuon[0] = ch.mum_isGlobalMuon[ibs]; mup_isGlobalMuon[0] = ch.mup_isGlobalMuon[ibs];
    	mum_NMuonStations[0] = ch.mum_NMuonStations[ibs]; mup_NMuonStations[0] = ch.mup_NMuonStations[ibs];


    	#   Global muon requirements from CMS AN-2008/098   #
    	#                 (without d0 cut)                  #
    	areMyGlobal[0] = 0    if (   ch.mum_isGlobalMuon[ibs] == 0 or ch.mup_isGlobalMuon[ibs] == 0 or
    	     				 ch.mum_normChi2[ibs] >= 10. or ch.mup_normChi2[ibs] >= 10. or
    	     				 ch.mum_normChi2[ibs] < 0. or ch.mup_normChi2[ibs] < 0. or
    	     				 ch.mumNHits[ibs] <= 10 or ch.mupNHits[ibs] <= 10)   else 1



##    	mum_dxy_Bsdecay[0]   = ch.mum_dxy_Bsdecay[ibs];    mum_dz_Bsdecay[0]      = ch.mum_dz_Bsdecay[ibs];
    	mum_isTrackerMuon[0] = ch.mum_isTrackerMuon[ibs]

        mum_LastStationOptimizedLowPtT[0] = ch.mum_LastStationOptimizedLowPtT[ibs]
        mum_LastStationT[0] = ch.mum_LastStationT[ibs]
        mum_OneStationT[0] = ch.mum_OneStationT[ibs]
        mum_LastStationAngT[0] = ch.mum_LastStationAngT[ibs]
        mum_OneStationAngT[0] = ch.mum_OneStationAngT[ibs]
        mum_2DCompatibilityT[0] = ch.mum_2DCompatibilityT[ibs]


##    	mup_dxy_Bsdecay[0]   = ch.mup_dxy_Bsdecay[ibs];    mup_dz_Bsdecay[0]      = ch.mup_dz_Bsdecay[ibs];
    	mup_isTrackerMuon[0] = ch.mup_isTrackerMuon[ibs];

        mup_LastStationOptimizedLowPtT[0] = ch.mup_LastStationOptimizedLowPtT[ibs]
        mup_LastStationT[0] = ch.mup_LastStationT[ibs]
        mup_OneStationT[0] = ch.mup_OneStationT[ibs]
        mup_LastStationAngT[0] = ch.mup_LastStationAngT[ibs]
        mup_OneStationAngT[0] = ch.mup_OneStationAngT[ibs]
        mup_2DCompatibilityT[0] = ch.mup_2DCompatibilityT[ibs]


        mu1_pt_cjp[0] = MU1P4_cjp.Pt()
        mu2_pt_cjp[0] = MU2P4_cjp.Pt()
        mu1_pt_0c[0] = MU1P4_0c.Pt()
        mu2_pt_0c[0] = MU2P4_0c.Pt()
    	#
##        if (not 'HLT_DoubleMu4_Jpsi_Displaced' in ch.triggersMuPL[ibs]) or (not 'HLT_DoubleMu4_Jpsi_Displaced' in ch.triggersMuML[ibs])  :continue
            #


    #####~~~~~~~~~~~~~~~~~~~~~#####
    ###~~~~~~~~~~J/psi~~~~~~~~~~###
    #####~~~~~~~~~~~~~~~~~~~~~#####


        if MU1P4_cjp.Pt() < 4. or MU2P4_cjp.Pt() < 4.:
            H_cuts.Fill(11)
            continue

        if  MUMUP4_C0.Pt() < 6.9:
            H_cuts.Fill(12)
            continue
        
        if ch.B_J_Prob[ibs] < 0.1:
            H_cuts.Fill(13)
            continue
        
        PDG_JPSI_MASS       =   3.096916
        if abs(ch.B_J_mass[ibs] - PDG_JPSI_MASS)    > 0.15  :continue  # were 3.04 and 3.15

        if DirectionCos2 ( JPV - PV, JPP3 ) < 0.9:
            H_cuts.Fill(14)
            continue

        if DetachSignificance2( JPV - PV, PVE, JPVE) < 3.0:
            H_cuts.Fill(15)
            continue
        if abs(MUMUP4_C0.Eta()) > 2.5  :continue

        JP_Eta_C0[0] = MUMUP4_C0.Eta()
        JP_Phi_C0[0] = MUMUP4_C0.Phi()
        
        JPsi_pt_C0[0]     = MUMUP4_C0.Pt()
        JPsi_mass_C0[0]   = MUMUP4_C0.M()
        JPsi_pvcos2_C0[0] = DirectionCos2 ( JPV - PV, JPP3 )
        JPsi_pv_DS_2D_C0[0]  = DetachSignificance2( JPV - PV, PVE, JPVE)
        
    	deltaR_mupmum_cjp[0] = MU1P4_cjp.DeltaR(MU2P4_cjp)
        JPsi_VtxProb[0]       = ch.B_J_Prob[ibs]

    #####~~~~~~~~~~~~~~~~~~~~~#####
    ###~~~~~~~~~~Kaons~~~~~~~~~~###
    #####~~~~~~~~~~~~~~~~~~~~~#####

##    	kaonP_pt_cjp[0] = kaonP_P4_cjp.Pt()
##    	kaonM_pt_cjp[0] = kaonM_P4_cjp.Pt()
##    	kaonP_pt_0c[0] = kaonP_P4_0c.Pt()
##    	kaonM_pt_0c[0] = kaonM_P4_0c.Pt()
##
##    	if kaonP_pt_0c[0] < 0.7  or kaonM_pt_0c[0] < 0.7 :continue
##        if abs(kaonP_P4_0c.Eta()) > 2.5 or abs(kaonM_P4_0c.Eta()) > 2.5 :continue
##
##
##        kaonP_track_normchi2[0] = ch.kaonP_track_normchi2[ibs]
##        kaonP_Hits[0] = ch.kaonP_Hits[ibs]
##        kaonP_PHits[0] = ch.kaonP_PHits[ibs]
##        kaonP_NTrackerLayers[0] = ch.kaonP_NTrackerLayers[ibs]
##        kaonP_NPixelLayers[0] = ch.kaonP_NPixelLayers[ibs]
##        kaonP_dxy_Bsdecay[0] = ch.kaonP_dxy_Bsdecay[ibs]
##        kaonP_dz_Bsdecay[0] = ch.kaonP_dz_Bsdecay[ibs]
##
##        kaonM_track_normchi2[0] = ch.kaonM_track_normchi2[ibs]
##        kaonM_Hits[0] = ch.kaonM_Hits[ibs]
##        kaonM_PHits[0] = ch.kaonM_PHits[ibs]
##        kaonM_NTrackerLayers[0] = ch.kaonM_NTrackerLayers[ibs]
##        kaonM_NPixelLayers[0] = ch.kaonM_NPixelLayers[ibs]
##        kaonM_dxy_Bsdecay[0] = ch.kaonM_dxy_Bsdecay[ibs]
##        kaonM_dz_Bsdecay[0] = ch.kaonM_dz_Bsdecay[ibs]
##
##    	deltaR_KpKm[0] = kaonP_P4_0c.DeltaR(kaonM_P4_0c)

    #####~~~~~~~~~~~~~~~~~~~~~#####
    ###~~~~~~~~~~Phi~~~~~~~~~~###
    #####~~~~~~~~~~~~~~~~~~~~~#####

##    	Phi_VtxProb[0] = ch.B_Phi_Prob[ibs]
##
##    	phi_mass_0c[0] = (kaonP_P4_0c + kaonM_P4_0c).M()
##    	phi_mass_cjp[0] = (kaonP_P4_cjp + kaonM_P4_cjp).M()
##
##        PDG_PHI_MASS        =   1.019455
##    	if abs(PDG_PHI_MASS - phi_mass_0c[0] ) > 0.015   :continue  #were 1.01-1.03
##
##    	phi_pt_0c[0] = (kaonP_P4_0c + kaonM_P4_0c).Pt()
##    	phi_pt_cjp[0] = (kaonP_P4_cjp + kaonM_P4_cjp).Pt()

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~Bs and misc.~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

##        Bs_mass_Cjp[0]          = ch.Bs_mass_cjp[ibs]
##        PDG_BS_MASS         =   5.36679
##        if abs(PDG_BS_MASS - Bs_mass_Cjp[0])  > 0.05   :continue # were 5.32 - 5.41
##
##        Bs_pt_Cjp[0]            = BsP4_Cjp.Pt()
####        if Bs_pt_Cjp[0] <   10.0    :continue
####        Bs_p_Cjp[0]             = BsP4_Cjp.Vect().Mag()
##
##        Bs_bcvtxDS2d_Cjp[0] = DetachSignificance2( BsV_Cjp - BcV, BcVE, BsVE_Cjp)
##        Bs_bcvtxDS2d_vtxfit[0] = DetachSignificance2( BsV_Cjp - BcV_vtxfit, BcVE_vtxfit, BsVE_Cjp)
##
##        Bs_pvDS2d_Cjp[0] = DetachSignificance2( BsV_Cjp - PV, PVE, BsVE_Cjp)
##        Bs_pv_detach_2D[0] =  sqrt( (BsV_Cjp - PV).X()**2 + (BsV_Cjp - PV).Y()**2)
####        if Bs_bcvtxDS2d_Cjp[0] < 3. :continue
##
##        Bs_bcvtx_cos2_Cjp[0]        = DirectionCos2 ( BsV_Cjp - BcV, BsP3_Cjp )
##        Bs_bcvtx_cos2_vtxfit[0]        = DirectionCos2 ( BsV_Cjp - BcV_vtxfit, BsP3_Cjp )
##        if Bs_bcvtx_cos2_Cjp[0] < 0.9 :
##            H_cuts.Fill(9)
##            continue
##
##        Bs_pv_cos2_Cjp[0] = DirectionCos2(BsV_Cjp - PV, BsP3_Cjp)
##
##        Bs_vtxprob_Cjp[0]       = ch.Bs_Prob[ibs]
##        if Bs_vtxprob_Cjp[0] < 0.01 :
##            H_cuts.Fill(10)
##            continue
##
##        if abs(BsP4_Cjp.Eta())  > 2.5   :continue
##
##        Bs_Phi_cjp[0]            = BsP4_Cjp.Phi()
##        Bs_Eta_cjp[0]            = BsP4_Cjp.Eta()
##
##        BsVtx_Chi2[0] = ch.Bs_Chi2[ibs]; Bs_Bcdecay_weight[0] = ch.Bs_Bcdecay_weight[ibs];


#####~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~Pion~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~#####

        pion_pt_0c[0] = ch.pionPF_Pt[ibs]
        if pion_pt_0c[0] < 2. :continue
##        pion_track_normchi2[0] = ch.pion_track_normchi2[ibs]
##        pion_Hits[0] = ch.pion_Hits[ibs]
##        pion_PHits[0] = ch.pion_PHits[ibs]
##        pion_NTrackerLayers[0] = ch.pion_NTrackerLayers[ibs]
##        pion_NPixelLayers[0] = ch.pion_NPixelLayers[ibs]
##        pion_dxy_Bcdecay[0] = ch.pion_dxy_Bcdecay[ibs]
##        pion_dz_Bcdecay[0] = ch.pion_dz_Bcdecay[ibs]
##        pion_Bcdecay_weight[0] = ch.pion_Bcdecay_weight[ibs]
	deltaR_pi_JP[0] = pionP4_0c.DeltaR(MUMUP4_C0)
        mva_nothing_nh[0] = ch.mva_nothing_nh[ibs]
        mva_gamma_nh[0] = ch.mva_gamma_nh[ibs]
        pion_mass[0] = pionP4_wo_mass.M()
        
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
###~~~~~~~~~~B0 and misc.~~~~~~~~~~###
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

        B0_mass[0] = ch.B0_mass[ibs]
##        Bc_mass[0]          = ch.Bc_mass[ibs]
##        Bc_mass_delta[0]    = ch.Bc_mass[ibs] - ch.Bs_mass_cjp[ibs] + PDG_BS_MASS
##        PDG_B0_MASS         =   5.2796
##        if abs(PDG_B0_MASS - B0_mass[0])   >   0.3    :continue
##
        B0_pt[0]            = B0P4.Pt()
        B0_Phi[0]            = B0P4.Phi()
        B0_Eta[0]            = B0P4.Eta()

        B0_pvDS2d[0] = DetachSignificance2( B0V - PV, PVE, B0VE)
##        B0_pvDS2d_vtxfit[0] = DetachSignificance2( B0V_vtxfit - PV, PVE, B0VE_vtxfit)
##        B0_pv_detach_2D[0] =  sqrt( (B0V - PV).X()**2 + (B0V - PV).Y()**2)
##        B0_pv_detach_2D_vtxfit[0] =  sqrt( (B0V_vtxfit - PV).X()**2 + (B0V_vtxfit - PV).Y()**2)

        B0_pvcos2[0]        = DirectionCos2 ( B0V - PV, B0P3 )
##        B0_pvcos2_vtxfit[0]        = DirectionCos2 ( B0V_vtxfit - PV, B0P3 )

##        B0_vtxprob[0]       = ch.B0_Prob[ibs]
##        if B0_vtxprob[0] < 0.01 :
##            H_cuts.Fill(10)
##            continue
##
        if B0_pvDS2d[0] < 3. :continue
        if B0_pvcos2[0] < 0.9 :
            H_cuts.Fill(9)
            continue
        if abs(B0P4.Eta())  > 2.5   :continue
##
##        B0Vertex_isValid[0] = ch.B0Vertex_isValid[ibs]
##        B0Vtx_Chi2_kinfit[0] = ch.B0_Chi2[ibs];   B0Vtx_Chi2_vtxfit[0] = ch.B0Vertex_Chi[ibs]; B0Vtx_normChi2_vtxfit[0] = ch.B0Vertex_normChi[ibs];
##
##
        if ch.PV_bestBang_RF_CL[ibs] < 0.8 :continue
        PV_refit_prob[0] = ch.PV_bestBang_RF_CL[ibs]
##        PV_bestBang_RF_X[0] = ch.PV_bestBang_RF_X[ibs]

           
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
