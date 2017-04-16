from ROOT import *; import glob, numpy as n; from array import array
from variables import *
isMC = 0


#__aa = 0;       __bb =  10 ;
# __aa = 0;       __bb =  14000 ;

# __aa = 0;       __bb =  4000 ; ## 
# __aa = 6000;    __bb =  14000 ; ## 

_fileOUT = 'Bs_star.root'


MyFileNamesMC = glob.glob( MCpath(1) + "*.root")
MyFileNamesDA = glob.glob('/afs/cern.ch/work/o/ofilatov/CMSSW_5_3_24/src/XbFrame/Xb_frame/crab_projects_Bfinder_Bs_star_v1/crab_Bfinder_Bs_star_v1_*/results/*.root')

#MyFileNamesDA = glob.glob("/afs/cern.ch/work/o/ofilatov/CMSSW_5_3_24/src/XbFrame/Xb_frame/crab_projects_Bfinder_Bc_star_v_pi_wo_point/crab_Bfinder_Bc_star_v_pi_wo_point_2012D*/results/*.root") ## 1
#MyFileNamesDA1 = glob.glob('/afs/cern.ch/work/o/ofilatov/CMSSW_5_3_24/src/XbFrame/Xb_frame/crab_projects_Bfinder_Bc_star_v_pi_wo_point/crab_Bfinder_Bc_star_v_pi_wo_point_2012C_201604-202324/results/*.root')
#MyFileNamesDA2 = glob.glob('/afs/cern.ch/work/o/ofilatov/CMSSW_5_3_24/src/XbFrame/Xb_frame/crab_projects_Bfinder_Bc_star_v_pi_wo_point/crab_Bfinder_Bc_star_v_pi_wo_point_2012C_202325-203045/results/*.root')
#MyFileNamesDA3 = glob.glob('/afs/cern.ch/work/o/ofilatov/CMSSW_5_3_24/src/XbFrame/Xb_frame/crab_projects_Bfinder_Bc_star_v_pi_wo_point/crab_Bfinder_Bc_star_v_pi_wo_point_2012C_203046-203742/results/*.root')

__aa = 0;  __bb =  len(MyFileNamesDA);
#__aa = 0;       __bb =  100 ;

MyFileNames = (MyFileNamesMC if isMC else MyFileNamesDA[__aa: __bb]); ch = TChain('mkcands/ntuple');
for fName in  MyFileNames:
    ii = ch.Add(fName);

print 'get ', len(MyFileNames), 'files from', __aa,'to',__bb,';  chain created'

fileOUT  = TFile (_fileOUT, "recreate");    mytree = TTree("mytree","mytree");

nEvt = ch.GetEntries(); print "entries: from", 0, 'to', nEvt-1; 
H_cuts = TH1F("H_cuts", "H_cuts", 40, 0, 20)

###  declaration and connecting to the branches of my new variables {{{1
NOUT, NOUT_evt, BBB, ibs = [int(0) for i in range(4)];
##pt_pi1_eq0, pt_pi2_eq0, pt_pi3_eq0, pt_JP_eq0, pt_Bc_eq0, pt_Bcstar_eq0 = [int(0) for i in range(6)];
PV, PVE, JPV, JPVE, JPP3, BcV_Cjp, BcP3_Cjp, BcVE_Cjp, _TV3zero = [TVector3() for i in range(9)]
BcP4_Cjp, BsP4_Cjp, Kplus_P4, Kminus_P4, MU1P4_cjp, MU2P4_cjp = [TLorentzVector() for i in range(6)];
_TV3zero  = TVector3(0,0,0)

_MY_VARS_ = [

"Pi1_pt", "Pi1_Eta",  'pion1_dxy_Bcdecay', 'pion1_dz_Bcdecay',                                                   'pi1_Bcdecay_weight', 'pi1_CommonVtx_weight',  'pion1_Hits', 'pion1_PHits', 'pion1_NTrackerLayers', 'pion1_NPixelLayers', 'pion1_trackchi2',
"Pi2_pt", "Pi2_Eta", 'pion2_dxy_Bcdecay', 'pion2_dz_Bcdecay', 'pion2_dxy_PV', 'pion2_dz_PV', 'pi2_PV_weight', 'pi2_CommonVtx_weight',  'pion2_Hits', 'pion2_PHits', 'pion2_NTrackerLayers', 'pion2_NPixelLayers', 'pion2_trackchi2',

"Bc_mass_Cjp", "Bc_pt_Cjp", "Bc_pvdistsignif2_Cjp", 
"Bc_pvcos2_Cjp", "Bc_vtxprob_Cjp",
"B_Eta_cjp", 

'DeltaM', "Bc2S_mass_Cjp", 'JP_pipi_mass',
"Bc2S_pt_Cjp", "Bc2S_Eta_cjp",
##"Bc2S_Phi_cjp",
##'PV_prob', 'PV_refit_prob',
'BcVertex_isValid', 'BcVertex_Chi', 'BcVertex_normChi', 'CommonVtx_isValid', 'CommonVtx_Chi', 'CommonVtx_normChi', 
'comp_BcVtx_CommonVtx', 'dist_BcVtx_CommonVtx', 'DS_BcVtx_CommonVtx',
'PV_refit_prob',

"DeltaR_JP_Pi_cjp",
##'DeltaR_Bc_Pi2_cjp', 'DeltaR_Bc_Pi3_cjp', 
## "Bc_p_Cjp", "Bc_pvdist_Cjp", 
## "Bc_properLT_Cjp", "Bc_properLTxy_Cjp",
## "Bc_chi2ndf_Cjp", "Bc_pvcossignif3_Cjp",

'areSoft', 'areTight_def', 'areTight_HM',
##'mum_relIso', 'mup_relIso',
'mum_isGlobalMuon', 'mup_isGlobalMuon', 
##'mum_normChi2', 'mup_normChi2',
##'mum_NMuonHits', 'mup_NMuonHits',
'mum_NMuonStations', 'mup_NMuonStations',

'JP_Eta_cjp',
##"JP_mass_Cmumu", "JP_Eta_cjp", "JP_Phi_cjp", 
"SAMEEVENT"]
_MC_VARS = ["MC_mu", "MC_k1"];
if isMC: _MY_VARS_ += _MC_VARS

for _var_ in _MY_VARS_:
    exec(_var_ + ' = n.zeros(1, dtype=float)')

for _var_ in _MY_VARS_:
    #print 'executing ' + 'mytree.Branch("' + _var_ + '"' + ' '*(25-len(_var_)) + ',' + _var_ + ' '*(25-len(_var_)) + ', "'+ _var_ + '/D")'
    exec('mytree.Branch("' + _var_ + '"' + ' '*(25-len(_var_)) + ',' + _var_ + ' '*(25-len(_var_)) + ', "'+ _var_ + '/D")')

###  declaration and connecting to the branches of my new variables }}}1

for evt in range(0, nEvt):
    ##
    if (ch.GetEntry(evt) <= 0) : break;
    BInfo_size  = ch.nB
    if len(ch.B_pion_px1) != BInfo_size :continue
    
    for Bj in range(BInfo_size):
        ##
        ibs = Bj
        ##
        MU1P4_cjp   .SetXYZM(ch.B_mu_px1_cjp[ibs], ch.B_mu_py1_cjp[ibs], ch.B_mu_pz1_cjp[ibs], PDG_MUON_MASS)
        MU2P4_cjp   .SetXYZM(ch.B_mu_px2_cjp[ibs], ch.B_mu_py2_cjp[ibs], ch.B_mu_pz2_cjp[ibs], PDG_MUON_MASS)
        Kplus_P4  .SetXYZM(ch.B_pion_px1[ibs], ch.B_pion_py1[ibs], ch.B_pion_pz1[ibs], PDG_KAON_MASS)
        Kminus_P4  .SetXYZM(ch.B_pion_px2[ibs], ch.B_pion_py2[ibs], ch.B_pion_pz2[ibs], PDG_KAON_MASS)

        BcP4_Cjp    .SetXYZM  ( ch.B_px[ibs], ch.B_py[ibs], ch.B_pz[ibs], ch.B_mass[ibs]) 
      	BsP4_Cjp = BcP4_Cjp + Kminus_P4
        #
        Bc2S_mass_Cjp[0] = BsP4_Cjp.M()

        JP_pipi_mass[0] = (MU1P4_cjp + MU2P4_cjp + Kminus_P4).M()
        #
	DeltaM[0] = BsP4_Cjp.M() - ch.B_mass[ibs] - PDG_KAON_MASS
	if DeltaM[0] > 0.1 :continue

        Pi1_pt[0]        = Kplus_P4.Pt()
        Pi2_pt[0]        = Kminus_P4.Pt()

        if Pi1_pt[0] < 1.  or Pi2_pt[0] < 0.4  :continue
        
	
        Pi1_Eta[0]       = Kplus_P4.Eta()
        Pi2_Eta[0]       = Kminus_P4.Eta()
        
        #
        if abs(Pi1_Eta[0]) > 2.4 or abs(Pi2_Eta[0]) > 2.4 :continue
        #
#	if ch.pion1_Hits[ibs] <= 4 or ch.pion1_PHits[ibs] <= 0 or ch.pion1_NTrackerLayers[ibs] <= 4 or ch.pion1_NPixelLayers[ibs] <= 0: continue
#	if ch.pion2_Hits[ibs] <= 4 or ch.pion2_PHits[ibs] <= 0 or ch.pion2_NTrackerLayers[ibs] <= 4 or ch.pion2_NPixelLayers[ibs] <= 0: continue

	pion1_dxy_Bcdecay[0] = ch.pion1_dxy_Bcdecay[ibs];     pion1_dz_Bcdecay[0] = ch.pion1_dz_Bcdecay[ibs];
	pi1_Bcdecay_weight[0] = ch.pi1_Bcdecay_weight[ibs];     pi1_CommonVtx_weight[0] = ch.pi1_CommonVtx_weight[ibs];
	pion1_trackchi2[0] = ch.pion1_trackchi2[ibs];     pion1_Hits[0] = ch.pion1_Hits[ibs];     pion1_PHits[0] = ch.pion1_PHits[ibs];     pion1_NTrackerLayers[0] = ch.pion1_NTrackerLayers[ibs];     pion1_NPixelLayers[0] = ch.pion1_NPixelLayers[ibs];

	pion2_dxy_Bcdecay[0] = ch.pion2_dxy_Bcdecay[ibs];     pion2_dz_Bcdecay[0] = ch.pion2_dz_Bcdecay[ibs];     pion2_dxy_PV[0] = ch.pion2_dxy_PV[ibs];     pion2_dz_PV[0] = ch.pion2_dz_PV[ibs];
	pi2_PV_weight[0] = ch.pi2_PV_weight[ibs];     pi2_CommonVtx_weight[0] = ch.pi2_CommonVtx_weight[ibs];
	pion2_trackchi2[0] = ch.pion2_trackchi2[ibs];     pion2_Hits[0] = ch.pion2_Hits[ibs];     pion2_PHits[0] = ch.pion2_PHits[ibs];     pion2_NTrackerLayers[0] = ch.pion2_NTrackerLayers[ibs];     pion2_NPixelLayers[0] = ch.pion2_NPixelLayers[ibs];


	# Soft J/psi muons #
	if (  ch.mum_NTrackerLayers[ibs] <= 5 or ch.mup_NTrackerLayers[ibs] <= 5 or 
	      ch.mum_NPixelLayers[ibs] <= 0 or ch.mup_NPixelLayers[ibs] <= 0 or
	      abs(ch.mumdxy[ibs]) >= 0.3 or abs(ch.mupdxy[ibs]) >= 0.3 or abs(ch.mumdz[ibs]) >= 20. or abs(ch.mupdz[ibs]) >= 20. or
	      ch.mumAngT[ibs] == 0 or ch.mupAngT[ibs] == 0  ):
	    areSoft[0] = 0
	else:
	    areSoft[0] = 1

####-------
	if areSoft[0] != 1:  continue
####-------

	# Default Tight J/psi muons #
	if ( ch.mum_isTight[ibs] == 0 or ch.mup_isTight[ibs] == 0):
	    areTight_def[0] = 0
	else:
	    areTight_def[0] = 1


	# Handmade Tight J/psi muons #
	if ( ch.mum_isGlobalMuon[ibs] == 0 or ch.mup_isGlobalMuon[ibs] == 0 or 
	     ch.mum_normChi2[ibs] >= 10. or ch.mup_normChi2[ibs] >= 10. or 
	     ch.mum_NMuonHits[ibs] <= 0 or ch.mup_NMuonHits[ibs] <= 0 or 
	     ch.mum_NMuonStations[ibs] <= 1 or ch.mup_NMuonStations[ibs] <= 1 or
	     abs(ch.mumdxy[ibs]) >= 0.2 or abs(ch.mupdxy[ibs]) >= 0.2 or abs(ch.mumdz[ibs]) >= 0.5 or abs(ch.mupdz[ibs]) >= 0.5 or
	     ch.mumNPHits[ibs] <= 0 or ch.mupNPHits[ibs] <= 0 or
	     ch.mum_NTrackerLayers[ibs] <= 5 or ch.mup_NTrackerLayers[ibs] <= 5	):
	    areTight_HM[0] = 0
	else:
	    areTight_HM[0] = 1

	mum_isGlobalMuon[0] = ch.mum_isGlobalMuon[ibs]; mup_isGlobalMuon[0] = ch.mup_isGlobalMuon[ibs];
	mum_NMuonStations[0] = ch.mum_NMuonStations[ibs]; mup_NMuonStations[0] = ch.mup_NMuonStations[ibs];

        #
        if (not 'HLT_DoubleMu4_Jpsi_Displaced' in ch.triggersMuPL[ibs]) or (not 'HLT_DoubleMu4_Jpsi_Displaced' in ch.triggersMuML[ibs]) :continue
	#
        MUMUP4_cjp = MU1P4_cjp + MU2P4_cjp
        JPV     = TVector3( ch.B_J_DecayVtxX[ibs],  ch.B_J_DecayVtxY[ibs],  ch.B_J_DecayVtxZ[ibs]   )
        JPVE    = TVector3( sqrt(ch.B_J_DecayVtxXE[ibs]), sqrt(ch.B_J_DecayVtxYE[ibs]), sqrt(ch.B_J_DecayVtxZE[ibs])  )
        JPP3    = TVector3( ch.B_J_px[ibs],         ch.B_J_py[ibs],         ch.B_J_pz[ibs])
        PV          = TVector3( ch.PV_bestBang_RF_X[ibs],   ch.PV_bestBang_RF_Y[ibs],   ch.PV_bestBang_RF_Z[ibs]    )
        PVE         = TVector3( sqrt(ch.PV_bestBang_RF_XE[ibs]),  sqrt(ch.PV_bestBang_RF_YE[ibs]),  sqrt(ch.PV_bestBang_RF_ZE[ibs])  )
        #
        if MU1P4_cjp.Pt() < 4.0 or MU2P4_cjp.Pt() < 4.0:
            H_cuts.Fill(11)
            continue
        #
        if MUMUP4_cjp.Pt() < 7.0:
            H_cuts.Fill(12)
	    continue
        if ch.B_J_Prob[ibs] < 0.1:
            H_cuts.Fill(13)
            continue
        #
        if ch.B_J_mass[ibs]   <   3.04    :continue
        if ch.B_J_mass[ibs]   >   3.15    :continue
        JPSI_pvcos2_Cmumu_va        = DirectionCos2 ( JPV - PV, JPP3 )
        JPSI_pvdistsignif2_Cmumu_va = DetachSignificance2( JPV - PV, PVE, JPVE)
        if JPSI_pvcos2_Cmumu_va < 0.9:
            H_cuts.Fill(14)
            continue
        #
        if JPSI_pvdistsignif2_Cmumu_va < 3.0:
            H_cuts.Fill(15)
            continue
        #
        JP_Eta_cjp[0] = MUMUP4_cjp.Eta()
        if abs(JP_Eta_cjp[0]) > 2.2  :continue
        #

        ###
        BcV_Cjp     = TVector3(ch.B_DecayVtxX[ibs],  ch.B_DecayVtxY[ibs],  ch.B_DecayVtxZ[ibs]   )
        BcVE_Cjp    = TVector3( sqrt(ch.B_DecayVtxXE[ibs]), sqrt(ch.B_DecayVtxYE[ibs]), sqrt(ch.B_DecayVtxZE[ibs])  )
        BcP3_Cjp    = BcP4_Cjp.Vect()
        #
        Bc_mass_Cjp[0]          = ch.B_mass[ibs]
	if ch.B_mass[ibs] < 5.1 or ch.B_mass[ibs] > 5.34 :continue
	#
        Bc_pt_Cjp[0]            = BcP4_Cjp.Pt()
        Bc_pvdistsignif2_Cjp[0] = DetachSignificance2( BcV_Cjp - PV, PVE, BcVE_Cjp)
	if Bc_pvdistsignif2_Cjp[0] < 3.   :continue
        Bc_pvcos2_Cjp[0]        = DirectionCos2 ( BcV_Cjp - PV, BcP3_Cjp )
	if Bc_pvcos2_Cjp[0] < 0.99   :continue
        #
        Bc_vtxprob_Cjp[0]       = ch.B_Prob[ibs]
        if Bc_vtxprob_Cjp[0] < 0.05 :
            H_cuts.Fill(10)
            continue
        #
        DeltaR_JP_Pi_cjp[0] = MUMUP4_cjp.DeltaR(Kplus_P4)
        #
        B_Eta_cjp[0]            = BcP4_Cjp.Eta()
        if abs(B_Eta_cjp[0])  > 2.4   :continue

        ###
        Bc2S_pt_Cjp[0] = BsP4_Cjp.Pt()
        Bc2S_Eta_cjp[0] = BsP4_Cjp.Eta()


        BcVertex_isValid[0] = ch.BcVertex_isValid[ibs];     BcVertex_Chi[0] = ch.BcVertex_Chi[ibs];     BcVertex_normChi[0] = ch.BcVertex_normChi[ibs];
        CommonVtx_isValid[0] = ch.CommonVtx_isValid[ibs];     CommonVtx_Chi[0] = ch.CommonVtx_Chi[ibs];     CommonVtx_normChi[0] = ch.CommonVtx_normChi[ibs];
        comp_BcVtx_CommonVtx[0] = ch.comp_BcVtx_CommonVtx[ibs];    dist_BcVtx_CommonVtx[0] = ch.dist_BcVtx_CommonVtx[ibs];     DS_BcVtx_CommonVtx[0] = ch.DS_BcVtx_CommonVtx[ibs];
        PV_refit_prob[0] = ch.PV_bestBang_RF_CL[ibs]
        ###
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
