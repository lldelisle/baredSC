##############
##    nx    ##
##############
nnorm=1
baredSC_2d \
    --input example/nih3t3_generated_second.txt \
    --geneXColName 1_0.5_0.5_0.5_0.5_0.5_x \
    --geneYColName 1_0.5_0.5_0.5_0.5_0.5_y \
    --metadata1ColName group \
    --metadata1Values group1 \
    --output example/second_example_2d_cellgroup1_${nnorm}gauss \
    --nnorm ${nnorm} \
    --figure example/second_example_2d_cellgroup1_${nnorm}gauss.png \
    --title "second example 2d cell group 1 ${nnorm} gauss"

baredSC_2d \
    --input example/nih3t3_generated_second.txt \
    --geneXColName 1_0.5_0.5_0.5_0.5_0.5_x \
    --geneYColName 1_0.5_0.5_0.5_0.5_0.5_y \
    --metadata1ColName group \
    --metadata1Values group1 \
    --output example/second_example_2d_cellgroup1_${nnorm}gauss_nx20 \
    --nnorm ${nnorm} \
    --nx 20 --ny 20 \
    --figure example/second_example_2d_cellgroup1_${nnorm}gauss_nx20.png \
    --title "second example 2d cell group 1 ${nnorm} gauss 20 bins"

baredSC_2d \
    --input example/nih3t3_generated_second.txt \
    --geneXColName 1_0.5_0.5_0.5_0.5_0.5_x \
    --geneYColName 1_0.5_0.5_0.5_0.5_0.5_y \
    --metadata1ColName group \
    --metadata1Values group1 \
    --output example/second_example_2d_cellgroup1_${nnorm}gauss_nx20 \
    --nnorm ${nnorm} \
    --nx 20 --ny 20 --prettyBinsx 50 --prettyBinsy 50 \
    --figure example/second_example_2d_cellgroup1_${nnorm}gauss_nx20_pretty.png \
    --title "second example 2d cell group 1 ${nnorm} gauss 20 bins pretty"
