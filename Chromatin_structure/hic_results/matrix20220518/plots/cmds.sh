awk -F'\t' 'BEGIN{pChr="";pline=""}{if($1!=pChr){if(pChr!="") {print pline} pChr=$1;print} pline=$0}END{print}' /mnt/ds28b/gongyh/Nanno/HiC/30/ordered/corrected/hicpro/imet1_hicpro_out3/hic_results/matrix/NannoH0/raw/10000/NannoH0_10000_abs.bed | cut -f1,4 > pos_start_end.txt


