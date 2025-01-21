#!/bin/bash

## call TAD for NannoH0

# 1. produce matrix
python sparseToDense3.py -c -d -b matrix/NannoH0/raw/5000/NannoH0_5000_abs.bed -o NannoH0_5000.matrix matrix/NannoH0/iced/5000/NannoH0_5000_iced.matrix

# 2. calc IS
for i in {1..30}; do
  perl ~/progs/crane-nature-2015/scripts/matrix2insulation.pl -i chr${i}_NannoH0_5000.matrix -is 30000 -ids 20000 -im mean -bmoe 3 -nt 0.1 -v
done

## call TAD for NannoH24

# 1. produce matrix
python sparseToDense3.py -c -d -b matrix/NannoH24/raw/5000/NannoH24_5000_abs.bed -o NannoH24_5000.matrix matrix/NannoH24/iced/5000/NannoH24_5000_iced.matrix

# 2. calc IS
for i in {1..30}; do
  perl ~/progs/crane-nature-2015/scripts/matrix2insulation.pl -i chr${i}_NannoH24_5000.matrix -is 30000 -ids 20000 -im mean -bmoe 3 -nt 0.1 -v
done


## call TAD for H0rep

# 1. produce matrix
python sparseToDense3.py -c -d -b matrix/H0rep/raw/5000/H0rep_5000_abs.bed -o H0rep_5000.matrix matrix/H0rep/iced/5000/H0rep_5000_iced.matrix

# 2. calc IS
for i in {1..30}; do
  perl ~/progs/crane-nature-2015/scripts/matrix2insulation.pl -i chr${i}_H0rep_5000.matrix -is 30000 -ids 20000 -im mean -bmoe 3 -nt 0.1 -v
done

## call TAD for H24rep

# 1. produce matrix
python sparseToDense3.py -c -d -b matrix/H24rep/raw/5000/H24rep_5000_abs.bed -o H24rep_5000.matrix matrix/H24rep/iced/5000/H24rep_5000_iced.matrix

# 2. calc IS
for i in {1..30}; do
  perl ~/progs/crane-nature-2015/scripts/matrix2insulation.pl -i chr${i}_H24rep_5000.matrix -is 30000 -ids 20000 -im mean -bmoe 3 -nt 0.1 -v
done

