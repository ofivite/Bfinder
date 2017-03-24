from ROOT import *
import glob 
from math import sqrt
import numpy as np

REMUC = True
##REMUC = False

ch = TChain("mytree");
ch.Add('~/BsJpsiPhi_v3_notall.root');

print "Adding chain done", ch.GetNtrees(), 'files '

mass_min = 5.2; mass_max = 5.5   
Bs_mass_Cjp = RooRealVar('Bs_mass_Cjp', 'm(B_{s}^{0}), GeV', mass_min, mass_max)
##Bs_varset = RooArgSet( Bc_pvcos2_Cjp, Bc_mass_Cjp, Bc_pt_Cjp, Bc_pvdistsignif2_Cjp, Bc_vtxprob_Cjp, DeltaR_JP_Pi_cjp, JP_Eta_cjp, Pi_Eta_cjp, Pi_pt_cjp)
Bs_varset = RooArgSet(Bs_mass_Cjp)
dataset = RooDataSet("ds","Dataset", Bs_varset)

a = []
nEvt = ch.GetEntries();
print "entries:", nEvt ;
par = 0.; par_0 = 9999.; flag_empty = 1;

cut = 1
if cut:
    for evt in range(nEvt):
        if ch.GetEntry(evt) <= 0 :  break
    ##    if (evt % 10000 == 0) :    ## printout progress
    ##        _perc = str(TMath.Nint(100*(evt-0)/(nEvt-0+0.0)));
    ##        print "["+_perc+(' ' * (3 - len(_perc)))+"%];evt",evt,";saved",dataset.numEntries()
    ##    
        if ((ch.SAMEEVENT < 1 or (not REMUC)) and flag_empty != 1): ## if new & !0, Write
            dataset.add(Bs_varset); a.append(b); flag_empty = 1;              ## now empty and reset par
            par_0 = 0;
    #----------------------------------------------------------

    #--------------------------------#
    ###------ Muons and J/Psi------###
    #--------------------------------#
             
    ##    if ch.mum_isGlobalMuon      != 1     :continue      
    ##    if ch.mup_isGlobalMuon      != 1     :continue      
##        if ch.areMyGlobal           != 1     :continue
    ##
        if ch.areSoft               != 1     :continue  
    ##    if ch.areTight_HM           != 1     :continue      
    ##    if ch.areTight_def          != 1     :continue    
##        if ch.mum_isTrackerMuon == 0 or ch.mup_isTrackerMuon == 0 or ch.mum_isGoodLS_OptimT == 0 or ch.mup_isGoodLS_OptimT == 0   :continue

    ##    if ch.mum_NMuonStations     <= 2.    :continue
    ##    if ch.mup_NMuonStations     <= 2.    :continue
        
   ##    if abs(ch.mum_dxy_Bcdecay) > 0.01 or abs(ch.mup_dxy_Bcdecay) > 0.01   :continue
    ##    if abs(ch.mum_dz_Bcdecay) > 0.05 or abs(ch.mup_dz_Bcdecay) > 0.05   :continue
    ##    
    ##    if ch.mum_relIso < -0.15 or ch.mup_relIso < -0.15    :continue
    ##    if ch.deltaR_mupmum_cjp > 0.5   :continue
##        if ch.JP_Bsdecay_weight < 0.5   :continue  ## -1 is included here
##        if ch.Jpsi_VtxProb < 0.1   :continue



    #----------------------#
    #####--- 1st kaon ---###
    #----------------------#
        
    ##    if abs(ch.pion2_dxy_Bsdecay)     > 0.02      :continue    
    ##    if abs(ch.pion2_dz_Bsdecay)      > 0.03      :continue

##        if ch.kaonP_pt_0c            <= 0.7   :continue    # 0.7 in MySel
##            
##        if ch.kaonP_Hits             <= 7.     :continue    
##        if ch.kaonP_PHits            <= 1.     :continue    
##        if ch.kaonP_NTrackerLayers   <= 4.     :continue    
##        if ch.kaonP_NPixelLayers     <= 1.     :continue    
##        if ch.kaonP_track_normchi2   >= 2.     :continue



    #----------------------#
    #####--- 2nd kaon ---###
    #----------------------#
        
    ##    if abs(ch.pion2_dxy_Bsdecay)     > 0.02      :continue    
    ##    if abs(ch.pion2_dz_Bsdecay)      > 0.03      :continue

##        if ch.kaonM_pt_0c            <= 0.7   :continue    # 0.7 in MySel
##            
##        if ch.kaonM_Hits             <= 7.     :continue    
##        if ch.kaonM_PHits            <= 1.     :continue    
##        if ch.kaonM_NTrackerLayers   <= 4.     :continue    
##        if ch.kaonM_NPixelLayers     <= 1.     :continue    
##        if ch.kaonM_track_normchi2   >= 2.     :continue


    #-----------------#
    #####--- Phi ---###
    #-----------------#

        if ch.phi_mass_0c < 1.01 or ch.phi_mass_0c > 1.03   :continue   # 1.01 and 1.03 in MySel
##        if ch.Phi_VtxProb < 0.1          :continue
##        if ch.phi_Bsdecay_weight < 0.5   :continue  ## -1 is included here
    ##    if ch.deltaR_KpKm < 0.7   :continue 
    ##    if ch.phi_pt_cjp < 1.4           :continue



    #--------------------------#
    ###------ Bs meson ------###
    #--------------------------#
        
        if ch.Bs_mass_Cjp < mass_min or ch.Bs_mass_Cjp > mass_max   :continue
        if ch.Bs_vtxprob_Cjp          < 0.05           :continue       # 0.05 in MySel  
        if ch.Bs_pvdistsignif2_Cjp    < 3.                :continue       # 3 in MySel              
        if ch.Bs_pvcos2_Cjp           < 0.99        :continue        # 0.9 in MySel                             
        
    ##    if ch.Bs_pt_Cjp               < 10.             :continue        # none in MySel
        if abs(ch.Bs_Eta_cjp)          > 2.5             :continue        # 2.5 in MySel

    ##    if ch.BsVertex_normChi   > 1.   :continue  ## -1 is included here
    ##    if ch.BsVertex_isValid   != 1     :continue  
    ##    if ch.PV_refit_prob      <= 0.5    :continue 


    #--------------------------------------------------------------
    #----------------------------------
    ### Choosing the best candidate ###
        par = ch.Bs_vtxprob_Cjp
        if  (not REMUC) or ((ch.SAMEEVENT == 1 and par > par_0) or ch.SAMEEVENT == 0) :
            par_0 = par;                            ### if better than was
            
            Bs_mass_Cjp   .setVal( ch.Bs_mass_Cjp ) ### and flag non empty
            b = ch.Bs_mass_Cjp
    ##        Bc_vtxprob_Cjp.            setVal(ch.Bc_vtxprob_Cjp)       
    ##        Bc_pvdistsignif2_Cjp.   setVal(ch.Bc_pvdistsignif2_Cjp)     
    ##        Bc_pt_Cjp.                     setVal(ch.Bc_pt_Cjp)   
    ##        Bc_pvcos2_Cjp.             setVal(ch.Bc_pvcos2_Cjp)
    ##        Pi_pt_cjp.                       setVal(ch.Pi_pt_cjp)
    ##        DeltaR_JP_Pi_cjp.         setVal(ch.DeltaR_JP_Pi_cjp)
    ##        JP_Eta_cjp.                     setVal(ch.JP_Eta_cjp)
    ##        Pi_Eta_cjp.                     setVal(ch.Pi_Eta_cjp)
            
            flag_empty = 0; 


    if flag_empty != 1 :
        dataset.add(Bs_varset);  ## write last event if needed
        a.append(b)
    dataset.Print();
    print 'dataset entries = ', dataset.sumEntries()

    fileOUT = TFile('~/cernbox/Backup/Bs to JPsi Phi/DataSet_Bs_v3_notall.root', 'recreate')
    dataset.Write()
    fileOUT.Close()


#---------------------------#
#### --- Let's plot! --- ####
#---------------------------#

f = TFile('~/cernbox/Backup/Bs to JPsi Phi/DataSet_Bs_v3_notall.root')
dataset = f.Get('ds')

frame = Bs_mass_Cjp.frame(RooFit.Title(""), RooFit.Bins(60))

mean_gauss = RooRealVar ("mean_1", "mean_1", 5.366, 5.35, 5.39)
sigma_1= RooRealVar ("sigma_1", "sigma_1", 0.01, 0.001, 0.03)
sigma_2= RooRealVar ("sigma_2", "sigma_2", 0.01, 0.001, 0.03)
exp_par = RooRealVar('exp_par', 'exp_par', -3, -10, -0.000000001)

N_gauss = RooRealVar('N_gauss', 'N_gauss', 30000, 0, 70000)
fr = RooRealVar('fr', 'fr', 0.5 , 0, 1)
N_gauss_1 = RooFormulaVar('N_gauss_1', 'N_gauss * fr', RooArgList(N_gauss, fr))
N_gauss_2 = RooFormulaVar('N_gauss_2', 'N_gauss * (1-fr)', RooArgList(N_gauss, fr))
N_backgr = RooRealVar('N_backgr', 'N_backgr', 3000, 0, 50000)

gauss_1 = RooGaussian('gauss_1', 'gauss_1', Bs_mass_Cjp, mean_gauss, sigma_1)
gauss_2 = RooGaussian('gauss_2', 'gauss_2', Bs_mass_Cjp, mean_gauss, sigma_2)
backgr = RooExponential('backgr', 'backgr', Bs_mass_Cjp, exp_par)

N_gauss.setPlotLabel('N_{sig}')
mean_gauss.setPlotLabel('m(B_{s}^{0})')
sigma_1.setPlotLabel('#sigma_{1}')
sigma_2.setPlotLabel('#sigma_{2}')
N_backgr.setPlotLabel('N_{bkg}')
exp_par.setPlotLabel('#lambda')
fr.setPlotLabel('fr')

model = RooAddPdf('model', 'model', RooArgList(gauss_1, gauss_2, backgr), RooArgList(N_gauss_1, N_gauss_2, N_backgr))



mean_gauss.setConstant(1)
model.fitTo(dataset)
model.fitTo(dataset)

mean_gauss.setConstant(0)
model.fitTo(dataset)
model.fitTo(dataset)


dataset.plotOn(frame)
model.paramOn(frame, RooFit.Layout(0.65))
model.plotOn(frame, RooFit.Name('model'), RooFit.LineColor(kRed-6), RooFit.LineWidth(5) )
model.plotOn(frame, RooFit.Name('backgr'), RooFit.Components('backgr'), RooFit.LineStyle(kDashed), RooFit.LineColor(kBlue-8), RooFit.LineWidth(4) )
model.plotOn(frame, RooFit.Name('gauss_1'), RooFit.Components('gauss_1'),RooFit.LineStyle(kDashed), RooFit.LineColor(47), RooFit.LineWidth(4), RooFit.Range(mean_gauss.getValV() - 4 * sigma_1.getValV(), mean_gauss.getValV() + 4 * sigma_1.getValV()))
model.plotOn(frame, RooFit.Name('gauss_2'), RooFit.Components('gauss_2'),RooFit.LineStyle(kDashed), RooFit.LineColor(48), RooFit.LineWidth(4), RooFit.Range(mean_gauss.getValV() - 4 * sigma_2.getValV(), mean_gauss.getValV() + 4 * sigma_2.getValV()))

frame.Draw()

##my_tkde = TKDE(len(a), np.array(a))
##my_tkde.Draw()
