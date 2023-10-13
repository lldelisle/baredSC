################
## 1d compare ##
################

for group in 0 3; do
  for nnorm in 1 2; do
    baredSC_1d \
          --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt \
          --metadata1ColName 0.5_0_0_0.5_group \
          --metadata1Values ${group}.0 \
          --geneColName 0.5_0_0_0.5_x \
          --output example/first_example_1d_group${group}_${nnorm}gauss_25_neff200 \
          --nnorm ${nnorm} --nx 25 --minNeff 200 \
          --figure example/first_example_1d_group${group}_${nnorm}gauss_25_neff200.png \
          --title "first gene ${nnorm} gauss group${group} 25 bins neff200" \
          --logevidence example/first_example_1d_group${group}_${nnorm}gauss_25_neff200_logevid.txt
  done
  combineMultipleModels_1d \
          --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt \
          --metadata1ColName 0.5_0_0_0.5_group \
          --metadata1Values ${group}.0 \
          --geneColName 0.5_0_0_0.5_x --nx 25 \
          --outputs example/first_example_1d_group${group}_1gauss_25_neff200 \
           example/first_example_1d_group${group}_2gauss_25_neff200 \
          --figure example/first_example_1d_group${group}_1-2gauss_25.png \
          --title "first gene group$group 25 bins 1 and 2 gauss"
done

for group in 1 2; do
        for nnorm in 1 2 3; do
          baredSC_1d \
                --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt \
                --metadata1ColName group \
                --metadata1Values group${group} \
                --geneColName 0.5_0_0_0.5_x \
                --output example/first_example_1d_cellgroup${group}_${nnorm}gauss_25_neff200 \
                --nnorm ${nnorm} --nx 25 --minNeff 200 \
            --figure example/first_example_1d_cellgroup${group}_${nnorm}gauss_25_neff200.png \
            --title "first gene ${nnorm} gauss cellgroup${group}" \
            --logevidence example/first_example_1d_cellgroup${group}_${nnorm}gauss_25_neff200_logevid.txt
        done
        combineMultipleModels_1d \
                --input baredSC/tests/test_data/nih3t3_generated_2d_2.txt \
                --metadata1ColName group \
                --metadata1Values group${group} \
                --geneColName 0.5_0_0_0.5_x --nx 25 \
                --outputs example/first_example_1d_cellgroup${group}_1gauss_25_neff200 \
                example/first_example_1d_cellgroup${group}_2gauss_25_neff200 \
                example/first_example_1d_cellgroup${group}_3gauss_25_neff200 \
                --figure example/first_example_1d_cellgroup${group}_1-3gauss_25.png \
                --title "first gene cellgroup$group 25 bins 1, 2 and 3 gauss"
      done
