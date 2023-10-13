##############
## minScale ##
##############
nnorm=1
baredSC_1d \
            --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt \
            --metadata1ColName 0.5_0_0_0.5_group \
            --metadata1Values 3.0 \
            --geneColName 0.5_0_0_0.5_x \
            --output example/first_example_1d_group3_${nnorm}gauss_ms0 \
            --nnorm ${nnorm} --minScale 0 \
            --minNeff 200 \
            --figure example/first_example_1d_group3_${nnorm}gauss_ms0.png \
            --title "first gene ${nnorm} gauss group3 minScale 0"

baredSC_1d \
            --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt \
            --metadata1ColName 0.5_0_0_0.5_group \
            --metadata1Values 3.0 \
            --geneColName 0.5_0_0_0.5_x \
            --output example/first_example_1d_group3_${nnorm}gauss_ms0.05 \
            --nnorm ${nnorm} --minScale 0.05 \
            --minNeff 200 \
            --figure example/first_example_1d_group3_${nnorm}gauss_ms0.05.png \
            --title "first gene ${nnorm} gauss group3 minScale 0.05"
