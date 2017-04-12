from ROOT import *

mass_min = 6.275 - 0.2; mass_max = 6.275 + 0.2; nbins = 40

Bc_mass = RooRealVar('Bc_mass', 'm(B_{c}), GeV', mass_min, mass_max)
Bc_mass_delta = RooRealVar('Bc_mass_delta', 'm(B_{c}) - m(B_{s}^{0}) + m_{PDG}(B_{s}^{0}), GeV', mass_min, mass_max)

Bc_pvcos2 = RooRealVar('Bc_pvcos2', 'cos', 0., 1.)
Bc_pvDS2d = RooRealVar('Bc_pvDS2d', 'DS1', 0., 100.)
Bc_pv_detach_2D = RooRealVar('Bc_pv_detach_2D', 'Bc_pv_detach_2D', 0., 10.)
Bc_vtxprob = RooRealVar('Bc_vtxprob', 'prob1', 0., 1.)
Bc_pt = RooRealVar('Bc_pt', 'Bc_pt', 0., 100.)
Bc_Set = RooArgSet(Bc_mass, Bc_pvcos2, Bc_pvDS2d, Bc_pvcos2, Bc_vtxprob, Bc_pt, Bc_mass_delta, Bc_pv_detach_2D)

Bs_pt_Cjp = RooRealVar('Bs_pt_Cjp', 'ptBs', 0., 100.)
Bs_bcvtx_cos2_Cjp = RooRealVar('Bs_bcvtx_cos2_Cjp', 'cos2', 0., 1.)
Bs_pv_cos2_Cjp = RooRealVar('Bs_pv_cos2_Cjp', 'cos3', 0., 1.)
Bs_bcvtxDS2d_Cjp = RooRealVar('Bs_bcvtxDS2d_Cjp', 'DS2', 0., 100.)
Bs_pv_detach_2D = RooRealVar('Bs_pv_detach_2D', 'Bs_pv_detach_2D', 0., 10.)
Bs_mass_Cjp = RooRealVar('Bs_mass_Cjp', 'm(B_{s}), GeV', 5., 6.)
Bs_vtxprob_Cjp = RooRealVar('Bs_vtxprob_Cjp', 'prob2', 0., 1.)
Bs_Bcdecay_weight = RooRealVar('Bs_Bcdecay_weight', 'weight2', 0., 1.)
Bs_pvDS2d_Cjp = RooRealVar('Bs_pvDS2d_Cjp', 'Bs_pvDS2d_Cjp', 0., 100.)
BsSet = RooArgSet(Bs_mass_Cjp, Bs_bcvtxDS2d_Cjp, Bs_bcvtx_cos2_Cjp, Bs_vtxprob_Cjp, Bs_pv_cos2_Cjp, Bs_pt_Cjp, Bs_pvDS2d_Cjp, Bs_pv_detach_2D)

pion_pt_0c = RooRealVar('pion_pt_0c', 'pi pt', 0., 10.)
pion_Bcdecay_weight = RooRealVar('pion_Bcdecay_weight', 'weight', 0., 1.)
pion_Hits = RooRealVar('pion_Hits', 'hits', 0, 50)
pion_PHits = RooRealVar('pion_PHits', 'pion_PHits', 0, 20)
pion_NTrackerLayers = RooRealVar('pion_NTrackerLayers', 'pion_NTrackerLayers', 0, 30)
pion_NPixelLayers = RooRealVar('pion_NPixelLayers', 'pion_NPixelLayers', 0, 10)
deltaR_piBs = RooRealVar('deltaR_piBs', 'deltaR_piBs', 0., 5.)
pion_track_normchi2 = RooRealVar('pion_track_normchi2', 'pion_track_normchi2', 0., 10.)
piSet = RooArgSet(pion_pt_0c, pion_Bcdecay_weight, pion_Hits, pion_track_normchi2, pion_PHits, pion_NTrackerLayers, pion_NPixelLayers)

phi_mass_0c = RooRealVar('phi_mass_0c', 'pi phi_mass_0c', 1., 1.05)
phi_pt_0c = RooRealVar('phi_pt_0c', 'pi phi_pt_0c', 0., 50.)
deltaR_KpKm = RooRealVar('deltaR_KpKm', 'deltaR_KpKm', 0., 10.)
kaonP_pt_0c = RooRealVar('kaonP_pt_0c', 'kaonP_pt_0c', 0., 50.)
kaonM_pt_0c = RooRealVar('kaonM_pt_0c', 'kaonM_pt_0c', 0., 50.)
phiSet = RooArgSet (phi_pt_0c, deltaR_KpKm, deltaR_piBs, kaonP_pt_0c, kaonM_pt_0c, phi_mass_0c)

pi_phiSet = RooArgSet(piSet, phiSet)
Bs_pi_Set = RooArgSet(pi_phiSet, BsSet)


BcSet = RooArgSet(Bs_pi_Set, Bc_Set)
