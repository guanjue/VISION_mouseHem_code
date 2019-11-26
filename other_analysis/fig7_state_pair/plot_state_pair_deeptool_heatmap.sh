cat pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{if ($17==13 && $6==7) print $2,$3,$4}' > mep13_cfue7.bed
cat pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{if ($16==13 && $11==7) print $2,$3,$4}' > lsk13_ery7.bed
cat pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{if ($16==13 && $11==13) print $2,$3,$4}' > lsk13_ery13.bed


cat pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{if ($16==13 && $13==7) print $2,$3,$4}' > lsk13_g1e7.bed
cat pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{if ($16==13 && $13==13) print $2,$3,$4}' > lsk13_g1e13.bed

time computeMatrix scale-regions -S LSK_BM.atac.0_16.bw G1E.atac.0_16.bw LSK_BM.ctcf.0_16.bw G1E.ctcf.0_16.bw -R lsk13_g1e13.bed lsk13_g1e7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_g1e13_7.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #C83296' 'white, #C83296' 'white, #C800FA' 'white, #C800FA' -m lsk13_g1e13_7.gz -out lsk13_g1e13_7.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 

time computeMatrix scale-regions -S LSK_BM.atac.0_16.bw CMP.atac.0_16.bw MEP.atac.0_16.bw CFU_E_ad.atac.0_16.bw ERY_ad.atac.0_16.bw ERY_fl.atac.0_16.bw -R lsk13_g1e13.bed lsk13_g1e7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_g1e13_7_atac.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' -m lsk13_g1e13_7_atac.gz -out lsk13_g1e13_7_atac.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 

time computeMatrix scale-regions -S LSK_BM.h3k4me1.0_16.bw CMP.h3k4me1.0_16.bw MEP.h3k4me1.0_16.bw CFU_E_ad.h3k4me1.0_16.bw ERY_ad.h3k4me1.0_16.bw ERY_fl.h3k4me1.0_16.bw -R lsk13_g1e13.bed lsk13_g1e7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_g1e13_7_h3k4me1.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' -m lsk13_g1e13_7_h3k4me1.gz -out lsk13_g1e13_7_h3k4me1.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 

time computeMatrix scale-regions -S LSK_BM.h3k27ac.0_16.bw CMP.h3k27ac.0_16.bw MEP.h3k27ac.0_16.bw CFU_E_ad.h3k27ac.0_16.bw ERY_ad.h3k27ac.0_16.bw ERY_fl.h3k27ac.0_16.bw -R lsk13_g1e13.bed lsk13_g1e7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_g1e13_7_h3k27ac.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, black' 'white, black' 'white, black' 'white, black' 'white, black' 'white, black' -m lsk13_g1e13_7_h3k27ac.gz -out lsk13_g1e13_7_h3k27ac.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 




time computeMatrix scale-regions -S LSK_BM.atac.0_16.bw ERY_ad.atac.0_16.bw LSK_BM.ctcf.0_16.bw ERY_ad.ctcf.0_16.bw -R lsk13_ery13_7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_ery13_7.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #C83296' 'white, #C83296' 'white, #C800FA' 'white, #C800FA' -m lsk13_ery13_7.gz -out lsk13_ery13_7.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 



time computeMatrix scale-regions -S LSK_BM.atac.0_16.bw ERY_ad.atac.0_16.bw LSK_BM.ctcf.0_16.bw ERY_ad.ctcf.0_16.bw -R lsk13_ery13.bed lsk13_ery7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_ery13_7.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #C83296' 'white, #C83296' 'white, #C800FA' 'white, #C800FA' -m lsk13_ery13_7.gz -out lsk13_ery13_7.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 


time computeMatrix scale-regions -S LSK_BM.atac.0_16.bw LSK_BM.ctcf.0_16.bw ERY_ad.atac.0_16.bw ERY_ad.ctcf.0_16.bw ERY_fl.atac.0_16.bw ERY_fl.ctcf.0_16.bw -R lsk13_ery13.bed lsk13_ery7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_ery13_7_eryfl.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #C83296' 'white, #C800FA' 'white, #C83296' 'white, #C800FA' 'white, #C83296' 'white, #C800FA' -m lsk13_ery13_7_eryfl.gz -out lsk13_ery13_7_eryfl.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 


time computeMatrix scale-regions -S LSK_BM.atac.0_16.bw CMP.atac.0_16.bw MEP.atac.0_16.bw CFU_E_ad.atac.0_16.bw ERY_ad.atac.0_16.bw ERY_fl.atac.0_16.bw LSK_BM.ctcf.0_16.bw ERY_ad.ctcf.0_16.bw ERY_fl.ctcf.0_16.bw -R lsk13_ery13.bed lsk13_ery7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_ery13_7_eryfl_atac.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C800FA' 'white, #C800FA' 'white, #C800FA' -m lsk13_ery13_7_eryfl_atac.gz -out lsk13_ery13_7_eryfl_atac.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 


time computeMatrix scale-regions -S LSK_BM.atac.0_16.bw CMP.atac.0_16.bw MEP.atac.0_16.bw CFU_E_ad.atac.0_16.bw ERY_ad.atac.0_16.bw ERY_fl.atac.0_16.bw LSK_BM.ctcf.0_16.bw ERY_ad.ctcf.0_16.bw ERY_fl.ctcf.0_16.bw -R lsk13_ery13.bed lsk13_ery7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_ery13_7_eryfl_atac.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C800FA' 'white, #C800FA' 'white, #C800FA' -m lsk13_ery13_7_eryfl_atac.gz -out lsk13_ery13_7_eryfl_atac.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 

time computeMatrix scale-regions -S LSK_BM.h3k27ac.0_16.bw CMP.h3k27ac.0_16.bw MEP.h3k27ac.0_16.bw CFU_E_ad.h3k27ac.0_16.bw ERY_ad.h3k27ac.0_16.bw ERY_fl.h3k27ac.0_16.bw -R lsk13_ery13.bed lsk13_ery7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_ery13_7_eryfl_h3k27ac.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #FA9600' 'white, #FA9600' 'white, #FA9600' 'white, #FA9600' 'white, #FA9600' 'white, #FA9600' -m lsk13_ery13_7_eryfl_h3k27ac.gz -out lsk13_ery13_7_eryfl_h3k27ac.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 

time computeMatrix scale-regions -S LSK_BM.h3k4me1.0_16.bw CMP.h3k4me1.0_16.bw MEP.h3k4me1.0_16.bw CFU_E_ad.h3k4me1.0_16.bw ERY_ad.h3k4me1.0_16.bw ERY_fl.h3k4me1.0_16.bw -R lsk13_ery13.bed lsk13_ery7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_ery13_7_eryfl_h3k4me1.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, black' 'white, black' 'white, black' 'white, black' 'white, black' 'white, black' -m lsk13_ery13_7_eryfl_h3k4me1.gz -out lsk13_ery13_7_eryfl_h3k4me1.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 


time computeMatrix scale-regions -S LSK_BM.atac.0_16.bw CMP.atac.0_16.bw MEP.atac.0_16.bw CFU_E_ad.atac.0_16.bw ERY_ad.atac.0_16.bw ERY_fl.atac.0_16.bw LSK_BM.ctcf.0_16.bw ERY_ad.ctcf.0_16.bw ERY_fl.ctcf.0_16.bw -R lsk13_ery13.bed lsk13_ery7.bed --beforeRegionStartLength 5000 --regionBodyLength 200 --afterRegionStartLength 5000 -o lsk13_ery13_7_eryfl_atac.gz --binSize 200 --numberOfProcessors 8 --sortRegions keep --missingDataAsZero --averageTypeBins mean
time plotHeatmap --colorList 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C83296' 'white, #C800FA' 'white, #C800FA' 'white, #C800FA' -m lsk13_ery13_7_eryfl_atac.gz -out lsk13_ery13_7_eryfl_atac.gz.pdf --sortRegions no --zMax 16 --zMin 0 --yMin 0 --yMax 16 --startLabel s --endLabel e --heatmapHeight 28 --heatmapWidth 6 


~/group/software/ucsc/bedGraphToBigWig CMP.atac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes CMP.atac.0_16.bw

~/group/software/ucsc/bedGraphToBigWig CFU_E_ad.atac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes CFU_E_ad.atac.0_16.bw

~/group/software/ucsc/bedGraphToBigWig MEP.atac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes MEP.atac.0_16.bw

~/group/software/ucsc/bedGraphToBigWig ERY_ad.atac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes ERY_ad.atac.0_16.bw
~/group/software/ucsc/bedGraphToBigWig ERY_fl.atac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes ERY_fl.atac.0_16.bw
~/group/software/ucsc/bedGraphToBigWig LSK_BM.atac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes LSK_BM.atac.0_16.bw
~/group/software/ucsc/bedGraphToBigWig G1E.atac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes G1E.atac.0_16.bw


~/group/software/ucsc/bedGraphToBigWig ERY_ad.ctcf.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes ERY_ad.ctcf.0_16.bw
~/group/software/ucsc/bedGraphToBigWig ERY_fl.ctcf.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes ERY_fl.ctcf.0_16.bw
~/group/software/ucsc/bedGraphToBigWig LSK_BM.ctcf.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes LSK_BM.ctcf.0_16.bw
~/group/software/ucsc/bedGraphToBigWig G1E.ctcf.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes G1E.ctcf.0_16.bw
~/group/software/ucsc/bedGraphToBigWig LSK_BM.ctcf.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes LSK_BM.ctcf.0_16.bw

~/group/software/ucsc/bedGraphToBigWig MEP.ctcf.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes MEP.ctcf.0_16.bw


~/group/software/ucsc/bedGraphToBigWig ERY_ad.h3k27ac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes ERY_ad.h3k27ac.0_16.bw
~/group/software/ucsc/bedGraphToBigWig ERY_fl.h3k27ac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes ERY_fl.h3k27ac.0_16.bw
~/group/software/ucsc/bedGraphToBigWig LSK_BM.h3k27ac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes LSK_BM.h3k27ac.0_16.bw
~/group/software/ucsc/bedGraphToBigWig CMP.h3k27ac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes CMP.h3k27ac.0_16.bw
~/group/software/ucsc/bedGraphToBigWig MEP.h3k27ac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes MEP.h3k27ac.0_16.bw
~/group/software/ucsc/bedGraphToBigWig CFU_E_ad.h3k27ac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes CFU_E_ad.h3k27ac.0_16.bw
~/group/software/ucsc/bedGraphToBigWig G1E.h3k27ac.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes G1E.h3k27ac.0_16.bw


~/group/software/ucsc/bedGraphToBigWig ERY_ad.h3k4me1.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes ERY_ad.h3k4me1.0_16.bw
~/group/software/ucsc/bedGraphToBigWig ERY_fl.h3k4me1.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes ERY_fl.h3k4me1.0_16.bw
~/group/software/ucsc/bedGraphToBigWig LSK_BM.h3k4me1.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes LSK_BM.h3k4me1.0_16.bw
~/group/software/ucsc/bedGraphToBigWig CMP.h3k4me1.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes CMP.h3k4me1.0_16.bw
~/group/software/ucsc/bedGraphToBigWig MEP.h3k4me1.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes MEP.h3k4me1.0_16.bw
~/group/software/ucsc/bedGraphToBigWig CFU_E_ad.h3k4me1.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes CFU_E_ad.h3k4me1.0_16.bw
~/group/software/ucsc/bedGraphToBigWig G1E.h3k4me1.0_16.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes G1E.h3k4me1.0_16.bw







