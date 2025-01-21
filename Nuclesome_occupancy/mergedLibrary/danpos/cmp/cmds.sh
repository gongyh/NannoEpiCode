cut -f1-4,13,18,23 H0-H24.positions.integrative.xls | awk -F'\t' 'NR>1{if(($5<0.01)||($6<0.01)||($7<0.01)) print}' | cut -f1-3 > DPNs.bed

cat genes.tss.bed | awk -F'\t' '{if($6=="+"){start=$2-1000;if(start<0){start=0;} print $1"\t"start"\t"$3"\t"$4}else{print $1"\t"$2"\t"$3+1000"\t"$4}}' | grep chr | grep "\.1$" > genes.tss1k.bed

bedtools intersect -wa -a genes.tss1k.bed -b DPNs.bed | uniq | cut -f4 | sed 's/\.1//g' > DPNs_genes.txt

$ cat DEGs_Up.txt DPNs_genes.txt | sort | uniq -d | wc -l
332
$ cat DEGs_Down.txt DPNs_genes.txt | sort | uniq -d | wc -l      
332

cat DEGs_Up.txt DPNs_genes.txt | sort | uniq -d > DPNs_DEGs_Up.txt
cat DEGs_Down.txt DPNs_genes.txt | sort | uniq -d > DPNs_DEGs_Down.txt
cat DPNs_DEGs_Up.txt DPNs_DEGs_Down.txt > DPNs_DEGs_all.txt
cat DPNs_DEGs_all.txt DPNs_genes.txt | sort | uniq -u > DPNs_nDEGs.txt

$ wc -l *.txt
  1512 DEGs_Down.txt
  1873 DEGs_Up.txt
   664 DPNs_DEGs_all.txt
   332 DPNs_DEGs_Down.txt
   332 DPNs_DEGs_Up.txt
  1892 DPNs_genes.txt
  1228 DPNs_nDEGs.txt
 10334 H0_H24_DEG.txt

$ cat DPNs_DEGs_Up.txt dna_binding.gids | sort | uniq -d
NO03G00600
NO18G02680

$ cat DPNs_DEGs_Down.txt dna_binding.gids | sort | uniq -d
NO02G05290
NO06G02300
NO12G03440
NO14G00340
NO14G01300
NO22G01340
NO29G01790


$ cut -f1,4,13,18,23 H0-H24.positions.integrative.xls | head
chr     center  smt_diff_FDR    point_diff_FDR  fuzziness_diff_FDR
chr13   2191    0.0326753330124 0.0326753330124 0.0465414861178
chr13   2970    0.358112662494  0.29160648371   0.500208634248
chr13   4013    0.654918953619  0.658481784625  0.209789760873
chr13   4202    0.668399935805  0.669427058257  0.290386775798
chr13   4556    0.417108008345  0.417108008345  0.308682394479
chr13   4917    0.478863745787  0.478863745787  0.171433156797
chr13   5103    0.495361900177  0.486759749639  0.782506820735
chr13   5275    0.648306852833  0.646990852191  0.290771946718
chr13   5463    0.653602952977  0.653827636013  0.647504413417

cut -f1,4,13,18,23 H0-H24.positions.integrative.xls | awk -F'\t' 'NR>1{if($3<0.01) print}' | wc -l #smt_diff_FDR 2031
cut -f1,4,13,18,23 H0-H24.positions.integrative.xls | awk -F'\t' 'NR>1{if($4<0.01) print}' | wc -l #point_diff_FDR 2030
cut -f1,4,13,18,23 H0-H24.positions.integrative.xls | awk -F'\t' 'NR>1{if($5<0.01) print}' | wc -l #fuzziness_diff_FDR 3201

cut -f1,4,13,18,23 H0-H24.positions.integrative.xls | awk -F'\t' 'NR>1{if(($3<0.01)&&($4<0.01)) print}' | wc -l #smt&point 2005
cut -f1,4,13,18,23 H0-H24.positions.integrative.xls | awk -F'\t' 'NR>1{if(($3<0.01)&&($5<0.01)) print}' | wc -l #smt&fuzziness 753
cut -f1,4,13,18,23 H0-H24.positions.integrative.xls | awk -F'\t' 'NR>1{if(($4<0.01)&&($5<0.01)) print}' | wc -l #point&fuzziness 755
cut -f1,4,13,18,23 H0-H24.positions.integrative.xls | awk -F'\t' 'NR>1{if(($3<0.01)&&($4<0.01)&&($5<0.01)) print}' | wc -l #all 747

bedtools intersect -F 0.2 -a ../consensus_peaks.bed  -b DPNs.bed  -wa  > DPNs_consensus.bed
cut -f1-4,13,18,23 H0-H24.positions.integrative.xls | awk -F'\t' 'NR>1{if($5<0.01) print}' | cut -f1-3 > DPNs_smt.bed 
bedtools intersect -F 0.2 -a ../consensus_peaks.bed  -b DPNs_smt.bed -wa  > DPNs_smt_consensus.bed

