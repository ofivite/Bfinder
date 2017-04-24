from ROOT import *
from RooVarSpace import *
import glob 
from math import sqrt
import math


f = TFile('~/Bc_v2_parked_notall_44886df.root') 
##f = TFile('~/BcBspi_v1_a067.root')   
tree = f.Get('mytree')

PDG_BC_MASS   = 6.2751
PDG_BS_MASS   = 5.36679
PDG_PHI_MASS  = 1.019455
PDG_JPSI_MASS = 3.096916

## Change it also in RooVarSpace !
## ------ mass_min = PDG_BC_MASS - 0.3; mass_max = PDG_BC_MASS + 0.3; nbins = 50

cuts = ('1 > 0'

    #--------------------#
    ###------ Bc ------###
    #--------------------#
     
    + '&&' + 'Bc_pvcos2_vtxfit       > 0.999'   #b# 0.999 #m# 0.9
    + '&&' + 'Bc_vtxprob      > 0.05'           #b# 0.05 #m# 0.01
    + '&&' + 'Bc_pvDS2d_vtxfit       > 3'       #b# 3 #m# none
##    + '&&' + 'Bc_pv_detach_2D_vtxfit > 0.03'
##    + '&&' + 'Bc_pt           > 10'           #b#  #m# none

    #--------------------#
    ###------ Bs ------###
    #--------------------#
 
    + '&&' + 'Bs_bcvtx_cos2_vtxfit        > 0.999'   #b#  #m# 0.9
##    + '&&' + 'Bs_pv_cos2_Cjp        < 0.9999999'
    + '&&' + 'Bs_vtxprob_Cjp           > 0.05'    #b#  #m# 0.01
##    + '&&' + 'Bs_bcvtxDS2d_vtxfit         > 5'  
##    + '&&' + 'Bs_pvDS2d_Cjp         > 10'    
    + '&&' + 'Bs_pv_detach_2D         > Bc_pv_detach_2D_vtxfit' #b# +0.02
##    + '&&' + 'Bs_Bcdecay_weight   > 0.9'            #
    + '&&' + 'Bs_pt_Cjp                > 10'  #b#  #m# none        

    #---------------------#
    ###------ Phi ------###
    #---------------------#

    + '&&' + 'phi_pt_0c                > 2.'    #b# 2 #m# none     
##    + '&&' + 'deltaR_KpKm              < 0.5'  # -
##    + '&&' + 'kaonP_pt_0c > .7 && kaonM_pt_0c > .7' #b# ? #m# 0.7

    #------------------------#
    ###------ Masses ------###
    #------------------------#

    + '&&' + 'JPsi_mass_C0   > (' + str(PDG_JPSI_MASS) + ' - 0.15) && JPsi_mass_C0 < (' + str(PDG_JPSI_MASS) + ' + 0.15)' #b# ? #m# 0.15                                 
    + '&&' + 'phi_mass_0c   > (' + str(PDG_PHI_MASS) + ' - 0.015) && phi_mass_0c < (' + str(PDG_PHI_MASS) + ' + 0.015)' #b# ? #m# 0.015                                 
    + '&&' + 'Bs_mass_Cjp   > (' + str(PDG_BS_MASS) + ' - 0.03) && Bs_mass_Cjp < (' + str(PDG_BS_MASS) + ' + 0.03)'   #b# ? #m# 0.050                           #  5.32 and 5.41        
    + '&&' + 'Bc_mass       >' + str(mass_min) + ' && Bc_mass < ' + str(mass_max)                             #b# ? #m# 0.5        
##    + '&&' + 'Bc_mass_delta       >' + str(mass_min) + ' && Bc_mass_delta < ' + str(mass_max)             

    #--------------------------------#
    ###------ Muons and J/psi------###
    #--------------------------------#

##    + '&&' + 'mup_isGlobalMuon == 1 && mum_isGlobalMuon == 1' # +
##    + '&&' + 'mup_isTrackerMuon == 1 && mum_isTrackerMuon == 1' # -
##    + '&&' + 'mup_isTrackerMuon == 1 && mum_isTrackerMuon == 1 && mup_LastStationOptimizedLowPtT == 1 && mum_LastStationOptimizedLowPtT == 1'
##    + '&&' + 'mup_isTrackerMuon == 1 && mum_isTrackerMuon == 1 && mup_2DCompatibilityT == 1 && mum_2DCompatibilityT == 1'
        
##    + '&&' + 'areMyGlobal == 1' 
##    + '&&' + 'areTight_HM == 1' 
##    + '&&' + 'areSoft == 1' #+

    + '&&' + 'mu1_pt_cjp > 4. && mu2_pt_cjp > 4.' #b# 4 #m# 3.5
        
    + '&&' + 'JPsi_pt_C0 > 10' #b# 10 #m# 5.5 #t# 6.9 ?
    + '&&' + 'JPsi_VtxProb > 0.05' #b# 0.05 #m# 0.01
        
##    + '&&' + 'JPsi_pvcos2_C0 > .99'
##    + '&&' + 'JPsi_pv_DS_2D_C0 > 10.'
        
        
    #----------------------#
    ###------ Pion ------###
    #----------------------#

    + '&&' + 'pion_pt_0c            > 1.'   #b# 1 #m# none  
##    + '&&' + 'pion_Bcdecay_weight   > 0.9'   #+  (95)        
##    + '&&' + 'deltaR_piBs   < 1.'            #        
             
##    + '&&' + 'pion_track_normchi2   < 2'            # -  
##    + '&&' + 'pion_Hits   >= 7'            #   
##    + '&&' + 'pion_PHits   >= 3'            #   
##    + '&&' + 'pion_NTrackerLayers   >= 10'            #   
##    + '&&' + 'pion_NPixelLayers   >= 1'            #
##    + '&&' + 'TMath::Abs(pion_dxy_Bcdecay)   <= 0.005'   #-

)


##tree.Draw('Bc_mass - Bs_mass_Cjp + 5.3668 >> hist', cuts)
##tree.Draw('Bc_mass >> hist', cuts)

dataset = RooDataSet('ds', 'ds', tree, RooArgSet(BcSet), cuts)
print 'dataset: ', dataset.sumEntries()


#-------------------------------------------------------------
frame = var.frame(RooFit.Title(""), RooFit.Bins(nbins))

model = RooAddPdf('model', 'model', RooArgList(gauss_1, gauss_2, backgr), RooArgList(N_gauss_1, N_gauss_2, N_backgr))
##model = RooAddPdf('model', 'model', RooArgList(gauss_1, backgr), RooArgList(N_gauss, N_backgr))

#-------------#
### sig+bkg fit
#-------------#

mean_gauss.setConstant(1)
sigma_1.setConstant(1)
sigma_2.setConstant(1)
fr.setConstant(1)

model.fitTo(dataset)
model.fitTo(dataset)


##mean_gauss.setConstant(0)
##sigma_1.setConstant(0)
##sigma_2.setConstant(0)
##fr.setConstant(0)

model.fitTo(dataset)
rrr_sig = model.fitTo(dataset, RooFit.Save())


dataset.plotOn(frame)
model.paramOn(frame, RooFit.Layout(0.65))
model.plotOn(frame, RooFit.Name('model'), RooFit.LineColor(kRed-6), RooFit.LineWidth(3) )
model.plotOn(frame, RooFit.Name('backgr'), RooFit.Components('backgr'), RooFit.LineStyle(kDashed), RooFit.LineColor(kBlue-8), RooFit.LineWidth(2) )
model.plotOn(frame, RooFit.Name('gauss_1'), RooFit.Components('gauss_1'),RooFit.LineStyle(kDashed), RooFit.LineColor(47), RooFit.LineWidth(2), RooFit.Range(mean_gauss.getValV() - 4 * sigma_1.getValV(), mean_gauss.getValV() + 4 * sigma_1.getValV()))
model.plotOn(frame, RooFit.Name('gauss_2'), RooFit.Components('gauss_2'),RooFit.LineStyle(kDashed), RooFit.LineColor(48), RooFit.LineWidth(2), RooFit.Range(mean_gauss.getValV() - 4 * sigma_2.getValV(), mean_gauss.getValV() + 4 * sigma_2.getValV()))

#--------------#
### bkg-only fit
#--------------#

N_gauss.setVal(0)
N_gauss.setConstant(1)

mean_gauss.setConstant(1)
sigma_1.setConstant(1)
sigma_2.setConstant(1)
fr.setConstant(1)

model.fitTo(dataset)
model.fitTo(dataset)

rrr_null = model.fitTo(dataset, RooFit.Save())

##model.paramOn(frame, RooFit.Layout(0.65))
##model.plotOn(frame, RooFit.Name('model'), RooFit.LineColor(kRed-1) )
##model.plotOn(frame, RooFit.Name('backgr'), RooFit.Components('backgr'), RooFit.LineStyle(kDashed), RooFit.LineColor(kRed-1))
##model.plotOn(frame, RooFit.Name('gauss'), RooFit.Components('gauss'), RooFit.LineStyle(kDashed), RooFit.LineColor(kRed-1), RooFit.Range(mean_gauss.getValV() - 4 * sigma.getValV(), mean_gauss.getValV() + 4 * sigma.getValV()))
frame.Draw()


nll_sig  = rrr_sig.minNll()
nll_null = rrr_null.minNll()
P = TMath.Prob(nll_null - nll_sig, 1)## !!! Change delta of ndf appropriately
S = TMath.ErfcInverse(P) * math.sqrt(2)

print 'P=', P, ' nll_sig=', nll_sig, ' nll_null=', nll_null, '\n', 'S=', S


