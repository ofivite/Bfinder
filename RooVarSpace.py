from ROOT import *

PDG_BC_MASS = 6.2751
mass_min = PDG_BC_MASS - 0.25; mass_max = PDG_BC_MASS + 0.25; nbins = 50

Bc_mass = RooRealVar('Bc_mass', 'm(B_{c}), GeV', mass_min, mass_max)
Bc_mass_delta = RooRealVar('Bc_mass_delta', 'm(B_{c}) - m(B_{s}^{0}) + m_{PDG}(B_{s}^{0}), GeV', mass_min, mass_max)

Bc_pvcos2 = RooRealVar('Bc_pvcos2', 'cos', 0., 1.)
Bc_pvcos2_vtxfit = RooRealVar('Bc_pvcos2_vtxfit', 'Bc_pvcos2_vtxfit', 0., 1.)
Bc_pvDS2d = RooRealVar('Bc_pvDS2d', 'DS1', 0., 1000.)
Bc_pvDS2d_vtxfit = RooRealVar('Bc_pvDS2d_vtxfit', 'Bc_pvDS2d_vtxfit', 0., 1000.)
Bc_pv_detach_2D = RooRealVar('Bc_pv_detach_2D', 'Bc_pv_detach_2D', 0., 100.)
Bc_pv_detach_2D_vtxfit = RooRealVar('Bc_pv_detach_2D_vtxfit', 'Bc_pv_detach_2D_vtxfit', 0., 100.)
Bc_vtxprob = RooRealVar('Bc_vtxprob', 'prob1', 0., 1.)
Bc_pt = RooRealVar('Bc_pt', 'Bc_pt', 0., 1000.)
##Bc_Set = RooArgSet(Bc_mass, Bc_pvcos2, Bc_pvDS2d, Bc_pvcos2, Bc_vtxprob, Bc_pt, Bc_mass_delta, Bc_pv_detach_2D)
Bc_Set = RooArgSet(Bc_mass, Bc_pvcos2_vtxfit, Bc_pvDS2d_vtxfit, Bc_vtxprob, Bc_pt, Bc_mass_delta, Bc_pv_detach_2D_vtxfit)

Bs_pt_Cjp = RooRealVar('Bs_pt_Cjp', 'ptBs', 0., 1000.)
Bs_bcvtx_cos2_Cjp = RooRealVar('Bs_bcvtx_cos2_Cjp', 'cos2', 0., 1.)
Bs_bcvtx_cos2_vtxfit = RooRealVar('Bs_bcvtx_cos2_vtxfit', 'Bs_bcvtx_cos2_vtxfit', 0., 1.)
Bs_pv_cos2_Cjp = RooRealVar('Bs_pv_cos2_Cjp', 'cos3', 0., 1.)
Bs_bcvtxDS2d_Cjp = RooRealVar('Bs_bcvtxDS2d_Cjp', 'DS2', 0., 1000.)
Bs_bcvtxDS2d_vtxfit = RooRealVar('Bs_bcvtxDS2d_vtxfit', 'Bs_bcvtxDS2d_vtxfit', 0., 1000.)
Bs_pv_detach_2D = RooRealVar('Bs_pv_detach_2D', 'Bs_pv_detach_2D', 0., 100.)
Bs_mass_Cjp = RooRealVar('Bs_mass_Cjp', 'm(B_{s}), GeV', 5., 6.)
Bs_vtxprob_Cjp = RooRealVar('Bs_vtxprob_Cjp', 'prob2', 0., 1.)
Bs_Bcdecay_weight = RooRealVar('Bs_Bcdecay_weight', 'weight2', 0., 1.)
Bs_pvDS2d_Cjp = RooRealVar('Bs_pvDS2d_Cjp', 'Bs_pvDS2d_Cjp', 0., 1000.)
##BsSet = RooArgSet(Bs_mass_Cjp, Bs_bcvtxDS2d_Cjp, Bs_bcvtx_cos2_Cjp, Bs_vtxprob_Cjp, Bs_pv_cos2_Cjp, Bs_pt_Cjp, Bs_pvDS2d_Cjp, Bs_pv_detach_2D)
BsSet = RooArgSet(Bs_mass_Cjp, Bs_bcvtxDS2d_vtxfit, Bs_bcvtx_cos2_vtxfit, Bs_vtxprob_Cjp, Bs_pv_cos2_Cjp, Bs_pt_Cjp, Bs_pvDS2d_Cjp, Bs_pv_detach_2D)

mu1_pt_cjp = RooRealVar('mu1_pt_cjp', 'mu1_pt_cjp', 0., 1000.)
mu2_pt_cjp = RooRealVar('mu2_pt_cjp', 'mu2_pt_cjp', 0., 1000.)
mup_isGlobalMuon = RooRealVar('mup_isGlobalMuon', 'mup_isGlobalMuon', 0., 1.)
mum_isGlobalMuon = RooRealVar('mum_isGlobalMuon', 'mum_isGlobalMuon', 0., 1.)
areSoft = RooRealVar('areSoft', 'areSoft', 0., 1.)
areTight_HM = RooRealVar('areTight_HM', 'areTight_HM', 0., 1.)
areMyGlobal = RooRealVar('areMyGlobal', 'areMyGlobal', 0., 1.) 
mup_LastStationOptimizedLowPtT = RooRealVar('mup_LastStationOptimizedLowPtT', 'mup_LastStationOptimizedLowPtT', 0., 1.)
mum_LastStationOptimizedLowPtT = RooRealVar('mum_LastStationOptimizedLowPtT', 'mum_LastStationOptimizedLowPtT', 0., 1.)
mup_2DCompatibilityT = RooRealVar('mup_2DCompatibilityT', 'mup_2DCompatibilityT', 0., 1.)
mum_2DCompatibilityT = RooRealVar('mum_2DCompatibilityT', 'mum_2DCompatibilityT', 0., 1.)
mup_isTrackerMuon = RooRealVar('mup_isTrackerMuon', 'mup_isTrackerMuon', 0., 1.)
mum_isTrackerMuon = RooRealVar('mum_isTrackerMuon', 'mum_isTrackerMuon', 0., 1.)

JPsi_mass_C0 = RooRealVar('JPsi_mass_C0', 'JPsi_mass_C0', 2.95, 3.25)
JPsi_pt_C0 = RooRealVar('JPsi_pt_C0', 'JPsi_pt_C0', 0, 1000.)
JPsi_pvcos2_C0 = RooRealVar('JPsi_pvcos2_C0', 'JPsi_pvcos2_C0', 0., 1.)
JPsi_pv_DS_2D_C0 = RooRealVar('JPsi_pv_DS_2D_C0', 'JPsi_pv_DS_2D_C0', 0., 1000.)
JPsi_VtxProb = RooRealVar('JPsi_VtxProb', 'JPsi_VtxProb', 0., 1.)
JP_Set = RooArgSet(JPsi_mass_C0, JPsi_pt_C0, JPsi_pvcos2_C0, JPsi_pv_DS_2D_C0)

muonSet = RooArgSet(mu1_pt_cjp, mu2_pt_cjp, JPsi_mass_C0, JPsi_pt_C0, JPsi_pvcos2_C0, JPsi_pv_DS_2D_C0, JPsi_VtxProb,
                    areTight_HM, areSoft)

pion_pt_0c = RooRealVar('pion_pt_0c', 'pi pt', 0., 100.)
pion_Bcdecay_weight = RooRealVar('pion_Bcdecay_weight', 'weight', 0., 1.)
pion_Hits = RooRealVar('pion_Hits', 'hits', 0, 100)
pion_PHits = RooRealVar('pion_PHits', 'pion_PHits', 0, 20)
pion_NTrackerLayers = RooRealVar('pion_NTrackerLayers', 'pion_NTrackerLayers', 0, 100)
pion_NPixelLayers = RooRealVar('pion_NPixelLayers', 'pion_NPixelLayers', 0, 30)
deltaR_piBs = RooRealVar('deltaR_piBs', 'deltaR_piBs', 0., 5.)
pion_track_normchi2 = RooRealVar('pion_track_normchi2', 'pion_track_normchi2', 0., 100.)
pion_dxy_Bcdecay = RooRealVar('pion_dxy_Bcdecay', 'pion_dxy_Bcdecay', -10., 10.)
piSet = RooArgSet(pion_pt_0c, pion_Bcdecay_weight, pion_Hits, pion_track_normchi2, pion_PHits, pion_NTrackerLayers, pion_NPixelLayers, pion_dxy_Bcdecay)

phi_mass_0c = RooRealVar('phi_mass_0c', 'pi phi_mass_0c', 1., 1.05)
phi_pt_0c = RooRealVar('phi_pt_0c', 'pi phi_pt_0c', 0., 100.)
Phi_VtxProb = RooRealVar('Phi_VtxProb', 'pi Phi_VtxProb', 0., 1.) 
deltaR_KpKm = RooRealVar('deltaR_KpKm', 'deltaR_KpKm', 0., 10.)
kaonP_pt_0c = RooRealVar('kaonP_pt_0c', 'kaonP_pt_0c', 0., 100.)
kaonM_pt_0c = RooRealVar('kaonM_pt_0c', 'kaonM_pt_0c', 0., 100.)
kaonP_dxy_Bsdecay = RooRealVar('kaonP_dxy_Bsdecay', 'kaonP_dxy_Bsdecay', -10., 10.)
kaonM_dxy_Bsdecay = RooRealVar('kaonM_dxy_Bsdecay', 'kaonM_dxy_Bsdecay', -10., 10.)
phiSet = RooArgSet (phi_pt_0c, deltaR_KpKm, deltaR_piBs, kaonP_pt_0c, kaonM_pt_0c, phi_mass_0c, Phi_VtxProb, kaonP_dxy_Bsdecay, kaonM_dxy_Bsdecay)

pi_phiSet = RooArgSet(piSet, phiSet)
Bs_pi_Set = RooArgSet(pi_phiSet, BsSet)
Bs_pi_muSet = RooArgSet(muonSet, Bs_pi_Set)


BcSet = RooArgSet(Bs_pi_muSet, Bc_Set)
var = Bc_mass

### --- Model
mean_gauss = RooRealVar ("mean_1", "mean_1", 6.27606, 6.27606 - 0.00020, 6.27606 + 0.00020) # 6.27588 +- 0.00021 for one
sigma_1= RooRealVar ("sigma_1", "sigma_1", 0.013, 0.013 - 0.00053, 0.013 + 0.00053) #0.01578 +- 0.00017 - for one, 0.013 +- 0.00053 - for two gaussians 
sigma_2= RooRealVar ("sigma_2", "sigma_2", 0.026, 0.0259 - 0.0024, 0.0259 + 0.0024) #this is for two gaussians
exp_par = RooRealVar('exp_par', 'exp_par', -3, -10, -0.000000001)

N_gauss = RooRealVar('N_gauss', 'N_gauss', 100, 0, 700)
fr = RooRealVar('fr', 'fr', 0.761 , 0.761 - 0.063, 0.761 + 0.063)
N_gauss_1 = RooFormulaVar('N_gauss_1', 'N_gauss * fr', RooArgList(N_gauss, fr))
N_gauss_2 = RooFormulaVar('N_gauss_2', 'N_gauss * (1-fr)', RooArgList(N_gauss, fr))
N_backgr = RooRealVar('N_backgr', 'N_backgr', 600, 0, 2000)

##mean_gauss = RooRealVar ("mean_1", "mean_1", 6.27606, 6.27606 - 0.10020, 6.27606 + 0.10020) # 6.27588 +- 0.00021 for one
##sigma_1= RooRealVar ("sigma_1", "sigma_1", 0.013, 0.013 - 0.10053, 0.013 + 0.10053) #0.01578 +- 0.00017 - for one, 0.013 +- 0.00053 - for two gaussians 
##sigma_2= RooRealVar ("sigma_2", "sigma_2", 0.026, 0.0259 - 0.1024, 0.0259 + 0.1024) #this is for two gaussians
##exp_par = RooRealVar('exp_par', 'exp_par', -3, -10, -0.000000001)
##
##N_gauss = RooRealVar('N_gauss', 'N_gauss', 100, 0, 700)
##fr = RooRealVar('fr', 'fr', 0.761 , 0., 1.)
##N_gauss_1 = RooFormulaVar('N_gauss_1', 'N_gauss * fr', RooArgList(N_gauss, fr))
##N_gauss_2 = RooFormulaVar('N_gauss_2', 'N_gauss * (1-fr)', RooArgList(N_gauss, fr))
##N_backgr = RooRealVar('N_backgr', 'N_backgr', 600, 0, 2000)

gauss_1 = RooGaussian('gauss_1', 'gauss_1', var, mean_gauss, sigma_1)
gauss_2 = RooGaussian('gauss_2', 'gauss_2', var, mean_gauss, sigma_2)
backgr = RooExponential('backgr', 'backgr', var, exp_par)
##backgr = RooPolynomial('backgr', 'backgr', var, RooArgList())

N_gauss.setPlotLabel('N_{sig}')
mean_gauss.setPlotLabel('m(B_{c})')
sigma_1.setPlotLabel('#sigma_{1}')
sigma_2.setPlotLabel('#sigma_{2}')
N_backgr.setPlotLabel('N_{bkg}')
exp_par.setPlotLabel('#lambda')
fr.setPlotLabel('fr')

