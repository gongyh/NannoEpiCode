#!/bin/bash

## GO enrichment of genes in each cluster for each histone mark under HC and LC conditons

# HC-H3K9ac-Cluster1
Rscript enrich_cluster.R  -g H0 -m H3k9ac -c cluster_1 -o H0_H3k9ac_c1_GO.pdf
# HC-H3K9ac-Cluster2
Rscript enrich_cluster.R  -g H0 -m H3k9ac -c cluster_2 -o H0_H3k9ac_c2_GO.pdf
# HC-H3K9ac-Cluster3
Rscript enrich_cluster.R  -g H0 -m H3k9ac -c cluster_3 -o H0_H3k9ac_c3_GO.pdf
# HC-H3K9ac-Cluster4
Rscript enrich_cluster.R  -g H0 -m H3k9ac -c cluster_4 -o H0_H3k9ac_c4_GO.pdf
# HC-H3K9ac-Cluster5
Rscript enrich_cluster.R  -g H0 -m H3k9ac -c cluster_5 -o H0_H3k9ac_c5_GO.pdf
# HC-H3K27ac-Cluster1

# HC-H3K27ac-Cluster2

# HC-H3K27ac-Cluster3

# HC-H3K27ac-Cluster4

# HC-H3K27ac-Cluster5

# HC-Kcr-Cluster1

# HC-Kcr-Cluster2

# HC-Kcr-Cluster3

# HC-Kcr-Cluster4

# HC-Kcr-Cluster5

# HC-H3K4me2-Cluster1

# HC-H3K4me2-Cluster2

# HC-H3K4me2-Cluster3

# HC-H3K4me2-Cluster4

# HC-H3K4me2-Cluster5

# LC-H3K9ac-Cluster1
Rscript enrich_cluster.R  -g H24 -m H3k9ac -c cluster_1 -o H24_H3k9ac_c1_GO.pdf
# LC-H3K9ac-Cluster2
Rscript enrich_cluster.R  -g H24 -m H3k9ac -c cluster_2 -o H24_H3k9ac_c2_GO.pdf
# LC-H3K9ac-Cluster3
Rscript enrich_cluster.R  -g H24 -m H3k9ac -c cluster_3 -o H24_H3k9ac_c3_GO.pdf
# LC-H3K9ac-Cluster4
Rscript enrich_cluster.R  -g H24 -m H3k9ac -c cluster_4 -o H24_H3k9ac_c4_GO.pdf
# LC-H3K9ac-Cluster5
Rscript enrich_cluster.R  -g H24 -m H3k9ac -c cluster_5 -o H24_H3k9ac_c5_GO.pdf
# LC-H3K27ac-Cluster1

# LC-H3K27ac-Cluster2

# LC-H3K27ac-Cluster3

# LC-H3K27ac-Cluster4

# LC-H3K27ac-Cluster5

# LC-Kcr-Cluster1

# LC-Kcr-Cluster2

# LC-Kcr-Cluster3

# LC-Kcr-Cluster4

# LC-Kcr-Cluster5

# LC-H3K4me2-Cluster1

# LC-H3K4me2-Cluster2

# LC-H3K4me2-Cluster3

# LC-H3K4me2-Cluster4

# LC-H3K4me2-Cluster5


