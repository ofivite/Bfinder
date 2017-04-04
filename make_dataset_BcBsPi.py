from ROOT import *
import glob 
from math import sqrt


f = TFile('~/BcBspi_v1.root')   
tree = f.Get('mytree')

mass_min = 6.; mass_max = 6.6; nbins = 30

cuts = (
             'Bc_pvcos2       > 0.99'           #  0.9 in MySel
##    + '&&' + 'Bc_vtxprob      > 0.05'           #  0.01 in MySel
    + '&&' + 'Bc_pvDS2d       > 2'              #  no
##    + '&&' + 'Bc_pt           > 10'             #  no        

    + '&&' + 'Bs_bcvtx_cos2_Cjp        > 0.99'           #  0.9 in MySel
##    + '&&' + 'Bs_pv_cos2_Cjp        < 0.99999'           #  0.9 in MySel
##    + '&&' + 'Bs_vtxprob_Cjp           > 0.1'           #  0.01 in MySel
    + '&&' + 'Bs_bcvtxDS2d_Cjp         > 1'              #  no
##    + '&&' + 'Bs_pt_Cjp                > 10'             #  no        


##    + '&&' + 'phi_mass_cjp  > 1.015 && phi_mass_cjp < 1.025'  #  1.01 and 1.03 for 0c !!
    + '&&' + 'Bs_mass_Cjp   > 5.33 && Bs_mass_Cjp < 5.41'                             #  5.32 and 5.41        
##    + '&&' + 'Bc_mass       >' + str(mass_min) + ' && Bc_mass < ' + str(mass_max)                             #  5.32 and 5.41        

    + '&&' + 'pion_pt_0c            > 1.'            #   
    + '&&' + 'pion_Bcdecay_weight   > 0.95'            #
             
##    + '&&' + 'pion_track_normchi2   < 2'            #   
##    + '&&' + 'pion_Hits   >= 5'            #   
##    + '&&' + 'pion_PHits   >= 3'            #   
##    + '&&' + 'pion_NTrackerLayers   >= 5'            #   
##    + '&&' + 'pion_NPixelLayers   >= 1'            #   

)


##tree.Draw('Bc_mass - Bs_mass_Cjp + 5.3668 >> hist', cuts)
tree.Draw('Bc_mass >> hist', cuts)

##
##Bc_mass = RooRealVar('Bc_mass', 'm(B_{c}), GeV', mass_min, mass_max)
##Bc_pvcos2 = RooRealVar('Bc_pvcos2', 'cos', 0., 1.)
##Bc_pvDS2d = RooRealVar('Bc_pvDS2d', 'DS1', 0., 100.)
##Bc_vtxprob = RooRealVar('Bc_vtxprob', 'prob1', 0., 1.)
##Bc_Set = RooArgSet(Bc_mass, Bc_pvcos2, Bc_pvDS2d, Bc_pvcos2, Bc_vtxprob)
##
##Bs_bcvtx_cos2_Cjp = RooRealVar('Bs_bcvtx_cos2_Cjp', 'cos2', 0., 1.)
##Bs_bcvtxDS2d_Cjp = RooRealVar('Bs_bcvtxDS2d_Cjp', 'DS2', 0., 100.)
##Bs_mass_Cjp = RooRealVar('Bs_mass_Cjp', 'm(B_{s}), GeV', 5., 6.)
##Bs_vtxprob_Cjp = RooRealVar('Bs_vtxprob_Cjp', 'prob2', 0., 1.)
##BsSet = RooArgSet(Bs_mass_Cjp, Bs_bcvtxDS2d_Cjp, Bs_bcvtx_cos2_Cjp, Bs_vtxprob_Cjp)
##
##pion_pt_0c = RooRealVar('pion_pt_0c', 'pi pt', 0., 10.)
##pion_Bcdecay_weight = RooRealVar('pion_Bcdecay_weight', 'weight', 0., 1.)
##piSet = RooArgSet(pion_pt_0c, pion_Bcdecay_weight)
##Bs_pi_Set = RooArgSet(piSet, BsSet)
##
##
##BcSet = RooArgSet(Bs_pi_Set, Bc_Set)
##dataset = RooDataSet('ds', 'ds', tree, RooArgSet(BcSet), cuts)
##
##frame = Bc_mass.frame(RooFit.Title(""), RooFit.Bins(nbins))
##
##### --- Model
##mean_gauss = RooRealVar ("mean_1", "mean_1", 6.27, 6.15, 6.55)
##sigma= RooRealVar ("sigma_1", "sigma_1", 0.02, 0.001, 0.05)
##exp_par = RooRealVar('exp_par', 'exp_par', -3, -10, -0.000000001)
##
##N_gauss = RooRealVar('N_gauss', 'N_gauss', 15, 0, 50)
##N_backgr = RooRealVar('N_backgr', 'N_backgr', 100, 0, 700)
##
##gauss = RooGaussian('gauss', 'gauss', Bc_mass, mean_gauss, sigma)
##backgr = RooExponential('backgr', 'backgr', Bc_mass, exp_par)
##
##N_gauss.setPlotLabel('N_{sig}')
##mean_gauss.setPlotLabel('m(B_{c})')
##sigma.setPlotLabel('#sigma')
##N_backgr.setPlotLabel('N_{bkg}')
##exp_par.setPlotLabel('#lambda')
##
##model = RooAddPdf('model', 'model', RooArgList(gauss, backgr), RooArgList(N_gauss, N_backgr))
##
##
##mean_gauss.setConstant(1)
##model.fitTo(dataset)
##model.fitTo(dataset)
##
##mean_gauss.setConstant(0)
##model.fitTo(dataset)
##model.fitTo(dataset)
##
##
##dataset.plotOn(frame)
##model.paramOn(frame, RooFit.Layout(0.65))
##model.plotOn(frame, RooFit.Name('model'), RooFit.LineColor(kRed-6), RooFit.LineWidth(3) )
##model.plotOn(frame, RooFit.Name('backgr'), RooFit.Components('backgr'), RooFit.LineStyle(kDashed), RooFit.LineColor(kBlue-8), RooFit.LineWidth(2) )
##model.plotOn(frame, RooFit.Name('gauss'), RooFit.Components('gauss'),RooFit.LineStyle(kDashed), RooFit.LineColor(47), RooFit.LineWidth(2), RooFit.Range(mean_gauss.getValV() - 4 * sigma.getValV(), mean_gauss.getValV() + 4 * sigma.getValV()))
##
##frame.Draw()
