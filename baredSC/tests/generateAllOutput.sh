
# baredSC_1d
baredSC_1d \
    --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt \
    --geneColName 0.5_0_0_0.5_x \
    --output example/first_example_1d_1gauss \
    --nnorm 1 \
    --figure example/first_example_1d_1gauss.png \
    --title "first gene 1 gauss" \
    --logevidence example/first_example_1d_1gauss_logevid.txt
baredSC_1d \
    --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt \
    --metadata1ColName 0_0.5_0.5_0_group \
    --metadata1Values 1.0 \
    --xmin -15 \
    --xmax -7 \
    --nx 25 \
    --xscale log \
    --minNeff 400 \
    --geneColName 0.5_0_0_0.5_x \
    --output example/first_example_1d_2gauss_log \
    --nnorm 2 \
    --figure example/first_example_1d_2gauss_log.pdf \
    --title "first gene 2 gauss log scale" \
    --logevidence example/first_example_1d_2gauss_log_logevid.txt

# Combined_1d
for nnorm in 1 2; do
    baredSC_1d --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt --geneColName 0.5_0_0_0.5_x --nnorm ${nnorm} --output baredSC/tests/test_data/small_${nnorm}gauss --nx 10 --nsampMCMC 20000 --force
done
combineMultipleModels_1d --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt --geneColName 0.5_0_0_0.5_x \
    --outputs baredSC/tests/test_data/small_1gauss baredSC/tests/test_data/small_2gauss \
    --nx 10 --title 'first gene combine 1 and 2 gauss' --prettyBins 100 --figure baredSC/tests/test_data/combine_test1.png

# baredSC 2d prep
for nnorm in 1 2; do
    baredSC_2d  --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt --geneXColName '0.5_0_0_0.5_x' --geneYColName '0.5_0_0_0.5_y' --nnorm ${nnorm} --output baredSC/tests/test_data/2d_small_${nnorm}gauss --nx 10 --ny 12 --nsampMCMC 20000 --force
done

# baredSC 2d plot
nnorm=2
baredSC_2d  --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt --geneXColName '0.5_0_0_0.5_x' --geneYColName '0.5_0_0_0.5_y' --nnorm ${nnorm} --output baredSC/tests/test_data/2d_small_${nnorm}gauss --nx 10 --ny 12 --nsampMCMC 20000 --figure baredSC/tests/test_data/2d_small_${nnorm}gauss.png --prettyBinsx 50 --prettyBinsy 50

# combined_2d
combineMultipleModels_2d --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt --geneXColName '0.5_0_0_0.5_x' --geneYColName '0.5_0_0_0.5_y' \
    --outputs baredSC/tests/test_data/2d_small_1gauss baredSC/tests/test_data/2d_small_2gauss \
    --nx 10 --ny 12 --prettyBinsx 50 --prettyBinsy 50 \
    --figure baredSC/tests/test_data/2d_small_combined.pdf
