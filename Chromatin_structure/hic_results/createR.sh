#!/bin/bash

for i in {1..30}; do
  echo "mapC(hiCH0\$chr${i}chr${i},hiCH24\$chr${i}chr${i},title=\"H0 vs H24 for chr${i}\")"
  echo "mapC(hiCH0rep\$chr${i}chr${i},hiCH24rep\$chr${i}chr${i},title=\"H0rep vs H24rep for chr${i}\")"
done

echo ""

for i in {1..30}; do
  echo "pc <- pca.hic(hiCH0\$chr${i}chr${i}, normPerExpected=TRUE, npc=1)"
  echo "plot(start(pc\$PC1), score(pc\$PC1), type=\"h\", xlab=\"H0 chr${i}\", ylab=\"PC1vec\", frame=FALSE)"
  echo "pc <- pca.hic(hiCH0rep\$chr${i}chr${i}, normPerExpected=TRUE, npc=1)"
  echo "plot(start(pc\$PC1), score(pc\$PC1), type=\"h\", xlab=\"H0rep chr${i}\", ylab=\"PC1vec\", frame=FALSE)"
  echo "pc <- pca.hic(hiCH24\$chr${i}chr${i}, normPerExpected=TRUE, npc=1)"
  echo "plot(start(pc\$PC1), score(pc\$PC1), type=\"h\", xlab=\"H24 chr${i}\", ylab=\"PC1vec\", frame=FALSE)"
  echo "pc <- pca.hic(hiCH24rep\$chr${i}chr${i}, normPerExpected=TRUE, npc=1)"
  echo "plot(start(pc\$PC1), score(pc\$PC1), type=\"h\", xlab=\"H24rep chr${i}\", ylab=\"PC1vec\", frame=FALSE)"
  echo ""
done

for i in {1..30}; do
  echo "di<-directionalityIndex(hiCH0\$chr${i}chr${i},winup=2e+4,windown=2e+4)"
  echo "barplot(di, col=ifelse(di>0,\"darkred\",\"darkgreen\"),main=\"H0 chr${i} DI\")"
  echo "di<-directionalityIndex(hiCH0rep\$chr${i}chr${i},winup=2e+4,windown=2e+4)"
  echo "barplot(di, col=ifelse(di>0,\"darkred\",\"darkgreen\"),main=\"H0rep chr${i} DI\")"
  echo "di<-directionalityIndex(hiCH24\$chr${i}chr${i},winup=2e+4,windown=2e+4)"
  echo "barplot(di, col=ifelse(di>0,\"darkred\",\"darkgreen\"),main=\"H24 chr${i} DI\")"
  echo "di<-directionalityIndex(hiCH24rep\$chr${i}chr${i},winup=2e+4,windown=2e+4)"
  echo "barplot(di, col=ifelse(di>0,\"darkred\",\"darkgreen\"),main=\"H24rep chr${i} DI\")"
  echo ""
done
