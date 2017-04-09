from ROOT import *
import glob 
from math import sqrt
import math

f = TFile('~/BcBspi_v1.root')   
tree = f.Get('mytree')

mass_min = 6.; mass_max = 6.6; nbins = 40

cuts = ('1 > 0'
    + '&&' + 'Bc_pvcos2       > 0.999'           #  0.9 in MySel
##    + '&&' + 'Bc_vtxprob      > 0.4'           #  0.01 in MySel
    + '&&' + 'Bc_pvDS2d       > 2'              #  no
##    + '&&' + 'Bc_pt           > 10'             #  no

    + '&&' + 'Bs_bcvtx_cos2_Cjp        > 0.999'           #  0.9 in MySel
##    + '&&' + 'Bs_pv_cos2_Cjp        < 0.99999'           #  0.9 in MySel
##    + '&&' + 'Bs_vtxprob_Cjp           > 0.05'           #  0.01 in MySel
    + '&&' + 'Bs_bcvtxDS2d_Cjp         > 2'              #  no
##    + '&&' + 'Bs_pt_Cjp                > 10'             #  no        

##    + '&&' + 'phi_mass_cjp  > 1.015 && phi_mass_cjp < 1.025'  #  1.01 and 1.03 for 0c !!
    + '&&' + 'Bs_mass_Cjp   > 5.34 && Bs_mass_Cjp < 5.39'                             #  5.32 and 5.41        
##    + '&&' + 'Bc_mass       >' + str(mass_min) + ' && Bc_mass < ' + str(mass_max)                             #  5.32 and 5.41        

    + '&&' + 'pion_pt_0c            > 1.'            #   
    + '&&' + 'pion_Bcdecay_weight   > 0.98'            #
             
##    + '&&' + 'pion_track_normchi2   < 3'            #   
##    + '&&' + 'pion_Hits   >= 5'            #   
##    + '&&' + 'pion_PHits   >= 1'            #   
##    + '&&' + 'pion_NTrackerLayers   >= 5'            #   
##    + '&&' + 'pion_NPixelLayers   >= 1'            #   

)


##tree.Draw('Bc_mass - Bs_mass_Cjp + 5.3668 >> hist', cuts)
tree.Draw('Bc_mass >> hist', cuts)


Bc_mass = RooRealVar('Bc_mass', 'm(B_{c}), GeV', mass_min, mass_max)
Bc_pvcos2 = RooRealVar('Bc_pvcos2', 'cos', 0., 1.)
Bc_pvDS2d = RooRealVar('Bc_pvDS2d', 'DS1', 0., 100.)
Bc_vtxprob = RooRealVar('Bc_vtxprob', 'prob1', 0., 1.)
Bc_Set = RooArgSet(Bc_mass, Bc_pvcos2, Bc_pvDS2d, Bc_pvcos2, Bc_vtxprob)

Bs_pt_Cjp = RooRealVar('Bs_pt_Cjp', 'ptBs', 0., 100.)
Bs_bcvtx_cos2_Cjp = RooRealVar('Bs_bcvtx_cos2_Cjp', 'cos2', 0., 1.)
Bs_pv_cos2_Cjp = RooRealVar('Bs_pv_cos2_Cjp', 'cos3', 0., 1.)
Bs_bcvtxDS2d_Cjp = RooRealVar('Bs_bcvtxDS2d_Cjp', 'DS2', 0., 100.)
Bs_mass_Cjp = RooRealVar('Bs_mass_Cjp', 'm(B_{s}), GeV', 5., 6.)
Bs_vtxprob_Cjp = RooRealVar('Bs_vtxprob_Cjp', 'prob2', 0., 1.)
Bs_Bcdecay_weight = RooRealVar('Bs_Bcdecay_weight', 'weight2', 0., 1.)
BsSet = RooArgSet(Bs_mass_Cjp, Bs_bcvtxDS2d_Cjp, Bs_bcvtx_cos2_Cjp, Bs_vtxprob_Cjp, Bs_pv_cos2_Cjp, Bs_pt_Cjp)

pion_pt_0c = RooRealVar('pion_pt_0c', 'pi pt', 0., 10.)
pion_Bcdecay_weight = RooRealVar('pion_Bcdecay_weight', 'weight', 0., 1.)
pion_Hits = RooRealVar('pion_Hits', 'hits', 0, 50)
piSet = RooArgSet(pion_pt_0c, pion_Bcdecay_weight, pion_Hits)
Bs_pi_Set = RooArgSet(piSet, BsSet)


BcSet = RooArgSet(Bs_pi_Set, Bc_Set)
dataset = RooDataSet('ds', 'ds', tree, RooArgSet(BcSet), cuts)

frame = Bc_mass.frame(RooFit.Title(""), RooFit.Bins(nbins))

### --- Model
mean_gauss = RooRealVar ("mean_1", "mean_1", 6.275, 6.1, 6.4)
sigma_1= RooRealVar ("sigma_1", "sigma_1", 0.0125, 0.001, 0.03)
sigma_2= RooRealVar ("sigma_2", "sigma_2", 0.0262, 0.001, 0.03)
exp_par = RooRealVar('exp_par', 'exp_par', -3, -10, -0.000000001)

N_gauss = RooRealVar('N_gauss', 'N_gauss', 10000, 0, 70000)
fr = RooRealVar('fr', 'fr', 0.603 , 0, 1)
N_gauss_1 = RooFormulaVar('N_gauss_1', 'N_gauss * fr', RooArgList(N_gauss, fr))
N_gauss_2 = RooFormulaVar('N_gauss_2', 'N_gauss * (1-fr)', RooArgList(N_gauss, fr))
N_backgr = RooRealVar('N_backgr', 'N_backgr', 60000, 0, 200000)

gauss_1 = RooGaussian('gauss_1', 'gauss_1', Bc_mass, mean_gauss, sigma_1)
gauss_2 = RooGaussian('gauss_2', 'gauss_2', Bc_mass, mean_gauss, sigma_2)
##backgr = RooExponential('backgr', 'backgr', Bs_mass_Cjp, exp_par)
backgr = RooPolynomial('backgr', 'backgr', Bc_mass, RooArgList())

N_gauss.setPlotLabel('N_{sig}')
mean_gauss.setPlotLabel('m(B_{c})')
sigma_1.setPlotLabel('#sigma_{1}')
sigma_2.setPlotLabel('#sigma_{2}')
N_backgr.setPlotLabel('N_{bkg}')
exp_par.setPlotLabel('#lambda')
fr.setPlotLabel('fr')

model = RooAddPdf('model', 'model', RooArgList(gauss_1, gauss_2, backgr), RooArgList(N_gauss_1, N_gauss_2, N_backgr))

mean_gauss.setConstant(1)
sigma_1.setConstant(1)
sigma_2.setConstant(1)
fr.setConstant(1)
model.fitTo(dataset)
model.fitTo(dataset)

##mean_gauss.setConstant(0)
model.fitTo(dataset)
rrr_sig = model.fitTo(dataset, RooFit.Save())


dataset.plotOn(frame)
model.paramOn(frame, RooFit.Layout(0.65))
model.plotOn(frame, RooFit.Name('model'), RooFit.LineColor(kRed-6), RooFit.LineWidth(3) )
model.plotOn(frame, RooFit.Name('backgr'), RooFit.Components('backgr'), RooFit.LineStyle(kDashed), RooFit.LineColor(kBlue-8), RooFit.LineWidth(2) )
model.plotOn(frame, RooFit.Name('gauss'), RooFit.Components('gauss'),RooFit.LineStyle(kDashed), RooFit.LineColor(47), RooFit.LineWidth(2), RooFit.Range(mean_gauss.getValV() - 4 * sigma_2.getValV(), mean_gauss.getValV() + 4 * sigma_2.getValV()))


### bkg-only fit

N_gauss.setVal(0)
N_gauss.setConstant(1)
##mean_gauss.setConstant(1)
##sigma_1.setConstant(1)
##sigma_2.setConstant(1)
##fr.setConstant(1)
model.fitTo(dataset)
model.fitTo(dataset)

rrr_null = model.fitTo(dataset, RooFit.Save())

##model.paramOn(frame, RooFit.Layout(0.65))
##model.plotOn(frame, RooFit.Name('model'), RooFit.LineColor(kRed-1) )
##model.plotOn(frame, RooFit.Name('backgr'), RooFit.Components('backgr'), RooFit.LineStyle(kDashed), RooFit.LineColor(kRed-1))
##model.plotOn(frame, RooFit.Name('gauss'), RooFit.Components('gauss'), RooFit.LineStyle(kDashed), RooFit.LineColor(kRed-1), RooFit.Range(mean_gauss.getValV() - 4 * sigma.getValV(), mean_gauss.getValV() + 4 * sigma.getValV()))
frame.Draw()

##, RooFit.NumEvents(dataset_CS.numEntries()), RooFit.DrawOption('F'), RooFit.FillStyle(3006), RooFit.LineColor(kWhite), RooFit.FillColor(kRed), 

nll_sig  = rrr_sig.minNll()
nll_null = rrr_null.minNll()
P = TMath.Prob(nll_null - nll_sig, 1)## !!! Change delta of ndf appropriately
S = TMath.ErfcInverse(P) * math.sqrt(2)

print 'P=', P, ' nll_sig=', nll_sig, ' nll_null=', nll_null, '\n', 'S=', S


