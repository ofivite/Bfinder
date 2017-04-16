from ROOT import *
from RooVarSpace import *
import glob 
from math import sqrt
import math


f = TFile('~/BcBspi_v1_a067.root') 
##f = TFile('~/BcBspi_v1_refTTfit_notall.root')   
tree = f.Get('mytree')

## Change it also in RooVarSpace !
mass_min = 6.275 - 0.2; mass_max = 6.275 + 0.2; nbins = 40


cuts = ('1 > 0'

    #--------------------#
    ###------ Bc ------###
    #--------------------#
     
    + '&&' + 'Bc_pvcos2       > 0.999'           #  0.9 in MySel
    + '&&' + 'Bc_vtxprob      > 0.05'           #  0.01 in MySel
    + '&&' + 'Bc_pvDS2d       > 3'              #  no
##    + '&&' + 'Bc_pv_detach_2D > 0.03'              #  no
##    + '&&' + 'Bc_pt           > 10'             #  no

    #--------------------#
    ###------ Bs ------###
    #--------------------#
 
    + '&&' + 'Bs_bcvtx_cos2_Cjp        > 0.999'           #  0.9 in MySel
##    + '&&' + 'Bs_pv_cos2_Cjp        < 0.9999999'           #  0.9 in MySel
    + '&&' + 'Bs_vtxprob_Cjp           > 0.05'         #  0.01 in MySel
##    + '&&' + 'Bs_bcvtxDS2d_Cjp         > 5'              #  no
##    + '&&' + 'Bs_pvDS2d_Cjp         > 3'              #  no
    + '&&' + 'Bs_pv_detach_2D         > Bc_pv_detach_2D '          #!!!  'Bs_pv_detach_2D         > Bc_pv_detach_2D + 0.02' with 'Bc_pvDS2d       > 3'
##    + '&&' + 'Bs_Bcdecay_weight   > 0.9'            #
##    + '&&' + 'Bs_pt_Cjp                > 10'             #  no        

    #---------------------#
    ###------ Phi ------###
    #---------------------#

##    + '&&' + 'phi_pt_0c                > 1'    # -         #  no        
##    + '&&' + 'deltaR_KpKm              < 0.5'  # -           #  no,
##    + '&&' + 'kaonP_pt_0c > 1. && kaonM_pt_0c > 1.'

    #------------------------#
    ###------ Masses ------###
    #------------------------#

##    + '&&' + 'phi_mass_0c  > 1.01 && phi_mass_0c < 1.03'  #  1.01 and 1.03 for 0c !!
    + '&&' + 'Bs_mass_Cjp   > (5.366 - 0.03) && Bs_mass_Cjp < (5.366 + 0.03)'                             #  5.32 and 5.41        
    + '&&' + 'Bc_mass       >' + str(mass_min) + ' && Bc_mass < ' + str(mass_max)                             #  5.32 and 5.41        
    + '&&' + 'Bc_mass_delta       >' + str(mass_min) + ' && Bc_mass_delta < ' + str(mass_max)                             #  5.32 and 5.41        

    #-----------------------#
    ###------ Muons ------###
    #-----------------------#

##    + '&&' + 'mup_isGlobalMuon == 1 && mum_isGlobalMuon == 1' # +
##    + '&&' + 'mup_isTrackerMuon == 1 && mum_isTrackerMuon == 1' # -
##    + '&&' + 'mup_isTrackerMuon == 1 && mum_isTrackerMuon == 1 && mup_LastStationOptimizedLowPtT == 1 && mum_LastStationOptimizedLowPtT == 1'
##    + '&&' + 'mup_isTrackerMuon == 1 && mum_isTrackerMuon == 1 && mup_2DCompatibilityT == 1 && mum_2DCompatibilityT == 1'
        
##    + '&&' + 'areMyGlobal == 1' # -
##    + '&&' + 'areSoft == 1'

    #----------------------#
    ###------ Pion ------###
    #----------------------#

    + '&&' + 'pion_pt_0c            > 1.'            #   
##    + '&&' + 'pion_Bcdecay_weight   > 0.95'   #+  (95)        
##    + '&&' + 'deltaR_piBs   < .8'            #        deltaR_piBs
             
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


