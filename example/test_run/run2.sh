#16. montage analysis 
#all interactions
cLoops2 montage -f ./gm_samp/chr21-chr21.ixy -bed ../data/runx1.bed -o gm_runx1_all -ext 2 -simple -ppmw 0.05 -vmax 500 
cLoops2 montage -f ./k562_samp/chr21-chr21.ixy -bed ../data/runx1.bed -o k562_runx1_all -ext 2 -simple -ppmw 0.05 -vmax 500 
#viewpoints mode for E19 and promoter
cLoops2 montage -f ./gm_samp/chr21-chr21.ixy -bed ../data/runx1.bed -o gm_runx1_vp -ext 2 -simple -ppmw 0.05 -vmax 500 -vp Promoter,E19 
cLoops2 montage -f ./k562_samp/chr21-chr21.ixy -bed ../data/runx1.bed -o k562_runx1_vp -ext 2 -simple -ppmw 0.05 -vmax 500 -vp Promoter,E19 

