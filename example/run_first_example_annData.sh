wget "https://usegalaxy.eu/datasets/11ac94870d0bb33ac1a474def9f936f9/display?to_ext=h5ad" -O /tmp/test.h5ad

# 1 gauss
# 0:00:43 of MCMC
# 0:00:05 of pdf
baredSC_1d \
      --inputAnnData /tmp/test.h5ad \
      --geneColName PIN3 \
      --output example/annData_1d_1gauss \
      --nnorm 1 \
      --xmax 4 \
      --metadata1ColName leiden \
      --metadata1Values columella+QC+NC \
      --figure example/annData_1d_1gauss.png \
      --title "annData_1d_1gauss" \
      --logevidence example/annData_1d_1gauss_logevid.txt

# 2 gauss
# 0:00:43 of MCMC
# 0:00:07 of pdf
baredSC_1d \
      --inputAnnData /tmp/test.h5ad \
      --geneColName PIN3 \
      --output example/annData_1d_2gauss \
      --nnorm 2 \
      --xmax 4 \
      --metadata1ColName leiden \
      --metadata1Values columella+QC+NC \
      --figure example/annData_1d_2gauss.png \
      --title "annData_1d_2gauss" \
      --logevidence example/annData_1d_2gauss_logevid.txt

baredSC_2d \
      --inputAnnData /tmp/test.h5ad \
      --geneXColName PIN3 \
      --geneYColName PIN7 \
      --output example/annData_2d_1gauss \
      --nnorm 1 \
      --xmax 4 \
      --ymax 4 \
      --metadata1ColName leiden \
      --metadata1Values columella+QC+NC \
      --figure example/annData_2d_1gauss.png \
      --title "annData_2d_1gauss" \
      --logevidence example/annData_2d_1gauss_logevid.txt
