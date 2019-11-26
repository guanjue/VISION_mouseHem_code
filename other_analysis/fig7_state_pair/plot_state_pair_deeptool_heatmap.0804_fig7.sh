cd /storage/home/g/gzx103/group/projects/vision/figure_7_state_transition_ct

head -100 /storage/home/g/gzx103/group/projects/vision/snapshot20_reproduce_2_16lim/ideas_list.all.txt

cat /storage/home/g/gzx103/group/projects/vision/snapshot20_reproduce_2_16lim/atac_20cell.function.matrix.no0.txt | awk -F '\t' -v OFS='\t' '{if ($7==9 && $12==12) print $1,int(($2+$3)/2)-100,int(($2+$3)/2)+100}' | sort -k1,1 -k2,2n > cmp9_eryad12.bed

cat /storage/home/g/gzx103/group/projects/vision/snapshot20_reproduce_2_16lim/atac_20cell.function.matrix.no0.txt | awk -F '\t' -v OFS='\t' '{if ($7==9 && $12==3) print $1,int(($2+$3)/2)-100,int(($2+$3)/2)+100}' | sort -k1,1 -k2,2n > cmp9_eryad3.bed


time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormLskAtac.bw ../bw/ideasVisionV20p8NormLskK4me1.bw ../bw/ideasVisionV20p8NormLskK27ac.bw ../bw/ideasVisionV20p8NormLskK27me3.bw -R cmp9_eryad12.bed cmp9_eryad3.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o cmp9_eryad12_eryad3.LSK.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #b3b300' 'white, #fa9600' 'white, #0000e1' -m cmp9_eryad12_eryad3.LSK.gz -out cmp9_eryad12_eryad3.LSK.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 14 --heatmapWidth 3 

time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormCmpAtac.bw ../bw/ideasVisionV20p8NormCmpK4me1.bw ../bw/ideasVisionV20p8NormCmpK27ac.bw ../bw/ideasVisionV20p8NormCmpK27ac.bw -R cmp9_eryad12.bed cmp9_eryad3.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o cmp9_eryad12_eryad3.CMP.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #b3b300' 'white, #fa9600' 'white, #0000e1' -m cmp9_eryad12_eryad3.CMP.gz -out cmp9_eryad12_eryad3.CMP.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 14 --heatmapWidth 3

time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormMepAtac.bw ../bw/ideasVisionV20p8NormMepK4me1.bw ../bw/ideasVisionV20p8NormMepK27ac.bw ../bw/ideasVisionV20p8NormMepK27ac.bw -R cmp9_eryad12.bed cmp9_eryad3.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o cmp9_eryad12_eryad3.MEP.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #b3b300' 'white, #fa9600' 'white, #0000e1' -m cmp9_eryad12_eryad3.MEP.gz -out cmp9_eryad12_eryad3.MEP.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 14 --heatmapWidth 3

time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormCfueAtac.bw ../bw/ideasVisionV20p8NormCfueK4me1.bw ../bw/ideasVisionV20p8NormCfueK27ac.bw ../bw/ideasVisionV20p8NormCfueK27ac.bw -R cmp9_eryad12.bed cmp9_eryad3.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o cmp9_eryad12_eryad3.CFUE.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #b3b300' 'white, #fa9600' 'white, #0000e1' -m cmp9_eryad12_eryad3.CFUE.gz -out cmp9_eryad12_eryad3.CFUE.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 14 --heatmapWidth 3

time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormEryadAtac.bw ../bw/ideasVisionV20p8NormEryadK4me1.bw ../bw/ideasVisionV20p8NormEryadK27ac.bw ../bw/ideasVisionV20p8NormEryadK27me3.bw -R cmp9_eryad12.bed cmp9_eryad3.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o cmp9_eryad12_eryad3.ERY_ad.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #b3b300' 'white, #fa9600' 'white, #0000e1' -m cmp9_eryad12_eryad3.ERY_ad.gz -out cmp9_eryad12_eryad3.ERY_ad.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 14 --heatmapWidth 3


cat /storage/home/g/gzx103/group/projects/vision/snapshot20_reproduce_2_16lim/atac_20cell.function.matrix.no0.txt | awk -F '\t' -v OFS='\t' '{if ($5==13 && $12==13) print $1,int(($2+$3)/2)-100,int(($2+$3)/2)+100}' | sort -k1,1 -k2,2n > lsk13_eryad13.bed
cat /storage/home/g/gzx103/group/projects/vision/snapshot20_reproduce_2_16lim/atac_20cell.function.matrix.no0.txt | awk -F '\t' -v OFS='\t' '{if ($5==13 && $12==7) print $1,int(($2+$3)/2)-100,int(($2+$3)/2)+100}' | sort -k1,1 -k2,2n > lsk13_eryad7.bed


time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormLskAtac.bw ../bw/ideasVisionV20p8NormLskCtcf.bw -R lsk13_eryad13.bed lsk13_eryad7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_eryad13_eryad7.LSK.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #C800FA' -m lsk13_eryad13_eryad7.LSK.gz -out lsk13_eryad13_eryad7.LSK.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 

time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormCmpAtac.bw ../bw/ideasVisionV20p8NormCmpAtac.bw -R lsk13_eryad13.bed lsk13_eryad7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_eryad13_eryad7.CMP.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #C800FA' -m lsk13_eryad13_eryad7.CMP.gz -out lsk13_eryad13_eryad7.CMP.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 

time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormMepAtac.bw ../bw/ideasVisionV20p8NormMepAtac.bw -R lsk13_eryad13.bed lsk13_eryad7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_eryad13_eryad7.MEP.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #C800FA' -m lsk13_eryad13_eryad7.MEP.gz -out lsk13_eryad13_eryad7.MEP.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 

time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormCfueAtac.bw ../bw/ideasVisionV20p8NormCfueAtac.bw -R lsk13_eryad13.bed lsk13_eryad7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_eryad13_eryad7.CFUE.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #C800FA' -m lsk13_eryad13_eryad7.CFUE.gz -out lsk13_eryad13_eryad7.CFUE.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 

time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormEryadAtac.bw ../bw/ideasVisionV20p8NormEryadCtcf.bw -R lsk13_eryad13.bed lsk13_eryad7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_eryad13_eryad7.ERY_ad.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #C800FA' -m lsk13_eryad13_eryad7.ERY_ad.gz -out lsk13_eryad13_eryad7.ERY_ad.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 

time computeMatrix scale-regions -S ../bw/ideasVisionV20p8NormEryflAtac.bw ../bw/ideasVisionV20p8NormEryflCtcf.bw -R lsk13_eryad13.bed lsk13_eryad7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_eryad13_eryad7.ERY_fl.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #c83296' 'white, #C800FA' -m lsk13_eryad13_eryad7.ERY_fl.gz -out lsk13_eryad13_eryad7.ERY_fl.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 






