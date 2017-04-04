from ROOT import *
import glob 
from math import sqrt


f = TFile('~/BcBspi_v1.root')   
tree = f.Get('mytree')


cuts = (
             'Bc_pvcos2       > 0.999'           #  0.9 in MySel
    + '&&' + 'Bc_vtxprob      > 0.05'           #  0.01 in MySel
    + '&&' + 'Bc_pvDS2d       > 3'              #  no
##    + '&&' + 'Bc_pt           > 10'             #  no        

    + '&&' + 'Bs_bcvtx_cos2_Cjp        > 0.999'           #  0.9 in MySel
    + '&&' + 'Bs_vtxprob_Cjp           > 0.05'           #  0.01 in MySel
    + '&&' + 'Bs_bcvtxDS2d_Cjp         > 5'              #  no
##    + '&&' + 'Bs_pt_Cjp                > 10'             #  no        


##    + '&&' + 'phi_mass_cjp  > 1.015 && phi_mass_cjp < 1.025'  #  1.01 and 1.03 for 0c !!
    + '&&' + 'Bs_mass_Cjp   > 5.34 && Bs_mass_Cjp < 5.40'                             #  5.32 and 5.41        
##    + '&&' + 'Bc_mass       > 6. && Bc_mass < 6.5'                             #  5.32 and 5.41        

##    + '&&' + 'pion_pt_0c      > .4'            #   0.4 in MySel

)


tree.Draw('Bc_mass >> hist', cuts)
