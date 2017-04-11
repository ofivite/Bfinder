from ROOT import *


mass_min = 6.275 - 0.2; mass_max = 6.275 + 0.2; nbins = 40
Bc_mass = RooRealVar('Bc_mass', 'm(B_{c}), GeV', mass_min, mass_max)
f1 = TFile('~/BcBspi_MC_signal_modif.root')
tree = f1.Get('mytree')
##f2 = TFile('~/Study/Machine Learning/PP scrypts/Bs_forest.root')
##hist_2 = f2.Get('mass_linreg')
data_1 = RooDataSet('data_1', 'data_1', tree, RooArgSet(Bc_mass))
##data_2 = RooDataHist('data_2', 'data_2', RooArgList(Bs_mass_Cjp), hist_2)

frame_1 = Bc_mass.frame(RooFit.Title(""), RooFit.Bins(nbins))
##frame_2 = Bs_mass_Cjp.frame(RooFit.Title("Random Forest"), RooFit.Bins(50))

mean_gauss = RooRealVar ("mean_1", "mean_1", 6.275, 6.1, 6.4)
sigma_1= RooRealVar ("sigma_1", "sigma_1", 0.01, 0.001, 0.03)
sigma_2= RooRealVar ("sigma_2", "sigma_2", 0.01, 0.001, 0.03)
exp_par = RooRealVar('exp_par', 'exp_par', -3, -10, -0.000000001)

N_gauss = RooRealVar('N_gauss', 'N_gauss', 10000, 0, 70000)
fr = RooRealVar('fr', 'fr', 0.5 , 0, 1)
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
##model = RooAddPdf('model', 'model', RooArgList(gauss_1, backgr), RooArgList(N_gauss, N_backgr))

### ML dataset
mean_gauss.setConstant(1)
model.fitTo(data_1)

mean_gauss.setConstant(0)
model.fitTo(data_1)
model.fitTo(data_1)


data_1.plotOn(frame_1)
model.paramOn(frame_1, RooFit.Layout(0.65))
model.plotOn(frame_1, RooFit.Name('model'), RooFit.LineColor(kRed-6), RooFit.LineWidth(3) )
model.plotOn(frame_1, RooFit.Name('backgr'), RooFit.Components('backgr'), RooFit.LineStyle(kDashed), RooFit.LineColor(kBlue-8), RooFit.LineWidth(2) )
model.plotOn(frame_1, RooFit.Name('gauss_1'), RooFit.Components('gauss_1'),RooFit.LineStyle(kDashed), RooFit.LineColor(47), RooFit.LineWidth(2), RooFit.Range(mean_gauss.getValV() - 4 * sigma_1.getValV(), mean_gauss.getValV() + 4 * sigma_1.getValV()))
model.plotOn(frame_1, RooFit.Name('gauss_2'), RooFit.Components('gauss_2'),RooFit.LineStyle(kDashed), RooFit.LineColor(48), RooFit.LineWidth(2), RooFit.Range(mean_gauss.getValV() - 4 * sigma_2.getValV(), mean_gauss.getValV() + 4 * sigma_2.getValV()))

frame_1.Draw()

##### original data
##mean_gauss.setConstant(1)
##model.fitTo(data_2)
##model.fitTo(data_2)
##
##mean_gauss.setConstant(0)
##model.fitTo(data_2)
##model.fitTo(data_2)
##
##frame_1.GetXaxis().SetLabelSize(0.03)
##frame_1.GetYaxis().SetLabelSize(0.03)
##frame_1.GetYaxis().SetTitleOffset(1.5)
##frame_1.GetXaxis().SetTitleOffset(1.1)
##
##frame_2.GetXaxis().SetLabelSize(0.03)
##frame_2.GetYaxis().SetLabelSize(0.03)
##frame_2.GetYaxis().SetTitleOffset(1.5)
##frame_2.GetXaxis().SetTitleOffset(1.1)
##
##data_2.plotOn(frame_2)
##model.paramOn(frame_2, RooFit.Layout(0.65))
##model.plotOn(frame_2, RooFit.Name('model'), RooFit.LineColor(kRed-6), RooFit.LineWidth(3) )
##model.plotOn(frame_2, RooFit.Name('backgr'), RooFit.Components('backgr'), RooFit.LineStyle(kDashed), RooFit.LineColor(kBlue-8), RooFit.LineWidth(2) )
##model.plotOn(frame_2, RooFit.Name('gauss_1'), RooFit.Components('gauss_1'),RooFit.LineStyle(kDashed), RooFit.LineColor(47), RooFit.LineWidth(2), RooFit.Range(mean_gauss.getValV() - 4 * sigma_1.getValV(), mean_gauss.getValV() + 4 * sigma_1.getValV()))
##model.plotOn(frame_2, RooFit.Name('gauss_2'), RooFit.Components('gauss_2'),RooFit.LineStyle(kDashed), RooFit.LineColor(48), RooFit.LineWidth(2), RooFit.Range(mean_gauss.getValV() - 4 * sigma_2.getValV(), mean_gauss.getValV() + 4 * sigma_2.getValV()))
##
##c = TCanvas("c","ML vs cut-based",2700,800)
##c.Divide(2,1)
##c.cd(1); frame_1.Draw()
##c.cd(2); frame_2.Draw()

##my_tkde = TKDE(len(a), np.array(a))
##my_tkde.Draw()
