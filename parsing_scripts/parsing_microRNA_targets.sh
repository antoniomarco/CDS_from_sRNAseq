#!/usr/bin/bash

# Mirtarbase 9.0 (2021) Human
wget https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/9.0/hsa_MTI.xlsx
ssconvert hsa_MTI.xlsx hsa_MTI.csv # Requires Gnumeric
cat hsa_MTI.csv | awk -F',' '{print $2 "\t" $4}' | sort | uniq > mirtarbase_hsa_9.tab
rm hsa_MTI.xlsx hsa_MTI.csv
cd ..

# Targetscan 8.0 Human
wget https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Targets_Info.default_predictions.txt.zip
wget https://www.targetscan.org/vert_80/vert_80_data_download/miR_Family_Info.txt.zip
zcat Predicted_Targets_Info.default_predictions.txt.zip | awk -F'\t' '$5 == "9606" {print $1 "\t" $3}' | sort | uniq | sort -k 1,1 > temp_targets
zcat miR_Family_Info.txt.zip | awk -F'\t' '$3 == 9606 {print $1 "\t" $4}' | sort | uniq | sort -k 1,1 > temp_families
join -1 1 -2 1 temp_targets temp_families | awk '{print $3 "\t" $2}' > targetscan_hsa_8.tab
rm temp_targets temp_families Predicted_Targets_Info.default_predictions.txt.zip  miR_Family_Info.txt.zip
cd ..

exit 0
