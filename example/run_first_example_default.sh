########
## 1d ##
########

# 1 gauss
# 0:00:43 of MCMC
# 0:00:05 of pdf
baredSC_1d \
      --input example/nih3t3_generated_2d_2.txt \
      --geneColName 0.5_0_0_0.5_x \
      --output example/first_example_1d_1gauss \
      --nnorm 1 \
      --figure example/first_example_1d_1gauss.png \
      --title "first gene 1 gauss" \
      --logevidence example/first_example_1d_1gauss_logevid.txt

# 2 gauss
# 0:00:43 of MCMC
# 0:00:07 of pdf
baredSC_1d \
      --input example/nih3t3_generated_2d_2.txt \
      --geneColName 0.5_0_0_0.5_x \
      --output example/first_example_1d_2gauss \
      --nnorm 2 \
      --figure example/first_example_1d_2gauss.png \
      --title "first gene 2 gauss" \
      --logevidence example/first_example_1d_2gauss_logevid.txt

# 3 gauss
# 0:00:46 of MCMC
# 0:00:08 of pdf
baredSC_1d \
      --input example/nih3t3_generated_2d_2.txt \
      --geneColName 0.5_0_0_0.5_x \
      --output example/first_example_1d_3gauss \
      --nnorm 3 \
      --figure example/first_example_1d_3gauss.png \
      --title "first gene 3 gauss" \
      --logevidence example/first_example_1d_3gauss_logevid.txt

# 3 gauss 1M
# 0:08:43 of MCMC
# 0:00:12 of pdf
baredSC_1d \
      --input example/nih3t3_generated_2d_2.txt \
      --geneColName 0.5_0_0_0.5_x \
      --output example/first_example_1d_3gauss_1M \
      --nnorm 3 --nsampMCMC 1000000 \
      --figure example/first_example_1d_3gauss_1M.png \
      --title "first gene 3 gauss 1M" \
      --logevidence example/first_example_1d_3gauss_1M_logevid.txt

# Combine the default
combineMultipleModels_1d \
      --input example/nih3t3_generated_2d_2.txt \
      --geneColName 0.5_0_0_0.5_x \
      --outputs example/first_example_1d_1gauss \
      example/first_example_1d_2gauss \
      example/first_example_1d_3gauss_1M \
      --figure example/first_example_1d_1-3gauss.png \
      --title "first gene 1, 2, and 3 gauss"

########
## 2d ##
########
for nnorm in 1 2 3; do
        baredSC_2d \
        --input example/nih3t3_generated_2d_2.txt \
        --geneXColName 0.5_0_0_0.5_x \
        --geneYColName 0.5_0_0_0.5_y \
        --metadata1ColName group \
        --metadata1Values group1 \
        --output example/first_example_2d_cellgroup1_${nnorm}gauss \
        --nnorm ${nnorm} \
        --figure example/first_example_2d_cellgroup1_${nnorm}gauss.png \
        --title "first example 2d cell group 1 ${nnorm} gauss" \
        --logevidence example/first_example_2d_cellgroup1_${nnorm}gauss_logevid.txt
done

# Redo the plot for splity
# pdf 0:02:23
baredSC_2d \
        --input example/nih3t3_generated_2d_2.txt \
        --geneXColName 0.5_0_0_0.5_x \
        --geneYColName 0.5_0_0_0.5_y \
        --metadata1ColName group \
        --metadata1Values group1 \
        --output example/first_example_2d_cellgroup1_2gauss \
        --nnorm 2 \
        --figure example/first_example_2d_cellgroup1_2gauss.png \
        --title "first example 2d cell group 1 2 gauss" \
        --splity 0.7 1.5
# Increase the nb of sample for 3 gauss:
nnorm=3
baredSC_2d \
        --input example/nih3t3_generated_2d_2.txt \
        --geneXColName 0.5_0_0_0.5_x \
        --geneYColName 0.5_0_0_0.5_y \
        --metadata1ColName group \
        --metadata1Values group1 \
        --nsampMCMC 1000000 \
        --nsampInPlot 75000 \
        --output example/first_example_2d_cellgroup1_1M_${nnorm}gauss \
        --nnorm ${nnorm} \
        --figure example/first_example_2d_cellgroup1_1M_${nnorm}gauss.png \
        --title "first example 2d cell group 1 1M ${nnorm} gauss" \
        --logevidence example/first_example_2d_cellgroup1_1M_${nnorm}gauss_logevid.txt

combineMultipleModels_2d \
          --input example/nih3t3_generated_2d_2.txt \
          --geneXColName 0.5_0_0_0.5_x \
          --geneYColName 0.5_0_0_0.5_y \
          --metadata1ColName group \
          --metadata1Values group1 \
          --outputs example/first_example_2d_cellgroup1_1gauss \
          example/first_example_2d_cellgroup1_2gauss \
          example/first_example_2d_cellgroup1_1M_3gauss \
          --figure example/first_example_2d_cellgroup1_1-3gauss.png \
          --getPVal \
          --title "first example cell group 1 1,2,3 gauss"

# Test:
# nnorm=1
# baredSC_2d \
#         --input example/nih3t3_generated_2d_2.txt \
#         --geneXColName 0.5_0_0_0.5_x \
#         --geneYColName 0.5_0_0_0.5_y \
#         --metadata1ColName group \
#         --metadata1Values group1 \
#         --output /tmp/test \
#         --xmax 2 --ymax 3 --nnorm $nnorm --nx 12 --ny 20 \
#         --prettyBinsx 60 --prettyBinsy 100 \
#         --figure /tmp/test.png \
#         --logevidence /tmp/test.txt
