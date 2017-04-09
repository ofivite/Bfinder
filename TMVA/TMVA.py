from ROOT import *

f_out = TFile('/home/jogurt/Study/Particles/Bfinder/TMVA/LearningOutput.root', 'recreate')
##f_sig = TFile('/home/jogurt/Study/Particles/B_c+/mva_sig.root')
f_sig = TFile('/home/jogurt/BcBspi_MC_signal.root')
f_bkg = TFile('/home/jogurt/BcBspi_v1.root')
sig_tree = f_sig.Get('mytree')
bkg_tree = f_bkg.Get('mytree')

var_list = ['kaonP_pt_0c',
'kaonP_track_normchi2', 'kaonP_Hits','kaonP_NTrackerLayers', 'kaonP_dxy_Bsdecay',
#'kaonP_PHits', 'kaonP_NPixelLayers',
        
'kaonM_pt_0c',
'kaonM_track_normchi2', 'kaonM_Hits', 'kaonM_NTrackerLayers', 'kaonM_dxy_Bsdecay',
#'kaonM_PHits', 'kaonM_NPixelLayers',
'deltaR_KpKm',
               
#'phi_Bsdecay_weight', 'Phi_VtxProb', 
'phi_pt_0c',

#'JP_Bsdecay_weight', 'Jpsi_VtxProb',

"Bs_pt_Cjp", "Bs_bcvtxDS2d_Cjp", "BsVtx_Chi2", 
'Bs_bcvtx_cos2_Cjp', 'Bs_Bcdecay_weight',
"Bs_Eta_cjp", "Bs_Phi_cjp", 

'pion_pt_0c',
'pion_track_normchi2', 'pion_Hits', 'pion_NTrackerLayers', 'pion_dxy_Bcdecay',
'pion_Bcdecay_weight',
#'kaonM_PHits', 'kaonM_NPixelLayers',
#'deltaR_piBs',

"Bc_pt", "Bc_pvDS2d", 'Bc_pvcos2', 
#"BcVtx_Chi2_kinfit", 
#'Bc_vtxprob',
"Bc_Eta", "Bc_Phi"]

factory = TMVA.Factory('Factory', f_out, 'AnalysisType=Classification')
for v in var_list:
    factory.AddVariable(v, 'F')
    
factory.AddSignalTree(sig_tree)
factory.AddBackgroundTree(bkg_tree)

sigCut   =  TCut ('1 > 0')#&& Pi_pt_cjp > 2 && Bc_vtxprob_Cjp > 0.15 && Bc_pvcos2_Cjp > 0.99')
bkgCut  =  TCut ('Bc_mass < 6.15 || Bc_mass > 6.4')
numb_str = 'nTrain_Signal=25000:nTest_Signal=25000:nTrain_Background=5000:nTest_Background=5000:'

factory.PrepareTrainingAndTestTree(sigCut, bkgCut, 'SplitMode=Random:NormMode=NumEvents:!V')
##factory.BookMethod(TMVA.Types.kSVM, "SVM", "C=1.0:Gamma=0.005:Tol=0.001:VarTransform=None")
##factory.BookMethod( TMVA.Types.kCuts, "CutsGA", "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" )
factory.BookMethod( TMVA.Types.kCuts, "Cuts", "H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" )
##factory.BookMethod( TMVA.Types.kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )

factory.TrainAllMethods() 
factory.TestAllMethods() 
factory.EvaluateAllMethods()

##Histo_training_S = TH1D('Histo_training_S','S (Train)',40,0.0,1.0) 
##Histo_training_B = TH1D('Histo_training_B','B (Train)',40,0.0,1.0) 
##Histo_testing_S = TH1D('Histo_testing_S','S (Test)',40,0.0,1.0) 
##Histo_testing_B = TH1D('Histo_testing_B','B (Test)',40,0.0,1.0)
##
##TrainTree = f_out.Get("TrainTree") 
##TestTree = f_out.Get("TestTree")
##
##SCut_Tree = 'classID>0.5'
##BCut_Tree = 'classID<0.5'
##TrainTree.Project("Histo_training_S","SVM",SCut_Tree)
##TrainTree.Project("Histo_training_B","SVM",BCut_Tree)
##TestTree.Project("Histo_testing_S","SVM",SCut_Tree)
##TestTree.Project("Histo_testing_B","SVM",BCut_Tree)
##
### Create the color styles
##Histo_training_S.SetLineColor(2)
##Histo_training_S.SetMarkerColor(2)
##Histo_training_S.SetFillColor(2)
##Histo_testing_S.SetLineColor(2)
##Histo_testing_S.SetMarkerColor(2)
##Histo_testing_S.SetFillColor(2)
## 
##Histo_training_B.SetLineColor(4)
##Histo_training_B.SetMarkerColor(4)
##Histo_training_B.SetFillColor(4)
##Histo_testing_B.SetLineColor(4)
##Histo_testing_B.SetMarkerColor(4)
##Histo_testing_B.SetFillColor(4)
## 
### Histogram fill styles
##Histo_training_S.SetFillStyle(3001)
##Histo_training_B.SetFillStyle(3001)
##Histo_testing_S.SetFillStyle(0)
##Histo_testing_B.SetFillStyle(0)
## 
### Histogram marker styles
##Histo_testing_S.SetMarkerStyle(20)
##Histo_testing_B.SetMarkerStyle(20)
##
### Set titles
##Histo_training_S.GetXaxis().SetTitle("Classifier, SVM [rbf kernel, C=1, gamma=0.005]")
##Histo_training_S.GetYaxis().SetTitle("Counts/Bin")
##
## # Draw the objects
##c1 = TCanvas("c1","",800,600)
##gStyle.SetOptStat(0)
##gStyle.SetOptTitle(0)
##Histo_training_S.Draw("HIST")
##Histo_training_B.Draw("HISTSAME")
##Histo_testing_S.Draw("EPSAME")
##Histo_testing_B.Draw("EPSAME")

f_sig.Close(); f_bkg.Close()
f_out.Close()

if gROOT.IsBatch() == 0:
     TMVA.TMVAGui('/home/jogurt/Study/Particles/Bfinder/TMVA/LearningOutput.root')

