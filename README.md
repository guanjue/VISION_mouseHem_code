# An integrative view of the regulatory and transcriptional landscapes in mouse hematopoiesis

# VISION_mouseHem_code
Code used in analysis of the Genome Research paper (2020) on the regulatory landscape and transcriptomes of mouse blood cells

In support of the manuscript "An integrative view of the regulatory and transcriptional landscapes in mouse hematopoiesis" by Guanjue Xiang, Cheryl A. Keller et al., we provide the code for the analyses in this GitHub repository.

The manuscript is under review after re-submission to Genome Research.

#### The scripts for the analyses in each figures is saved in each of the folder based on the folder name. It was compiled by Guanjue Xiang. Due to the size of the input file, some of the input files are saved in the following link: https://psu.box.com/s/3bmc2zndpm5kiwkinspr7a0chso9u1bm 

### Other packages that are used in the analysis:
### VISION project dataset portal: usevision.org
### IDEAS GitHub repository: https://github.com/seqcode/IDEAS
### S3norm GitHub repository: https://github.com/guanjue/S3norm

### This study is done in mm10

#### Abstract: Thousands of epigenomic datasets have been generated in the past decade, but it is difficult for researchers to effectively utilize all the data relevant to their projects. Systematic integrative analysis can help meet this need, and the VISION project was established for ValIdated Systematic IntegratiON of epigenomic data in hematopoiesis. Here, we systematically integrated extensive data recording epigenetic features and transcriptomes from many sources, including individual laboratories and consortia, to produce a comprehensive view of the regulatory landscape of differentiating hematopoietic cell types in mouse. By employing IDEAS as our Integrative and Discriminative Epigenome Annotation System, we identified and assigned epigenetic states simultaneously along chromosomes and across cell types, precisely and comprehensively. Combining nuclease accessibility and epigenetic states produced a set of over 200,000 candidate cis-regulatory elements (cCREs) that efficiently capture enhancers and promoters. The transitions in epigenetic states of these cCREs across cell types provided insights into mechanisms of regulation, including decreases in numbers of active cCREs during differentiation of most lineages, transitions from poised to active or inactive states, and shifts in nuclease accessibility of CTCF-bound elements. Regression modeling of epigenetic states at cCREs and gene expression produced a versatile resource to improve selection of cCREs potentially regulating target genes. These resources are available from our VISION website (usevision.org) to aid research in genomics and hematopoiesis. 



<img src="https://github.com/hardisonlab/vision_pipeline_mouse_mm10/blob/master/figures_for_github/vision_mouse.png" width="800"/>

##### Hematopoietic cell types and datasets used for integrative analysis. Schematic representation of the main lineage commitment steps in hematopoiesis, along with three immortalized cell lines (HPC7, G1E, G1E-ER4) and their approximate position relative to the primary cell populations shown. Abbreviations for cell populations are LSK = Lin-Sca1+Kit+ (which includes hematopoietic stem cells and multipotent progenitor cells), CMP = common myeloid progenitor cells, GMP = granulocyte monocyte progenitor cells, MEP = megakaryocyte erythrocyte progenitor cells, CLP = common lymphoid progenitor cells, CFUE = colony forming unit erythroid, ERY = erythroblasts, RBC = red blood cells, CFUMK = colony forming unit megakaryocyte, iMK = immature megakaryocytes, MK_fl = maturing megakaryocytes from fetal liver, PLTS = platelets, EOS = eosinophils, MAS = mast cells, NEU = neutrophils, MON = monocytes, T_CD8 = CD8+ T-cells, T_CD4 = CD4+ T-cells, B = B-cells, NK = natural killer cells.

##### Figure: Segmentation of the epigenomes of hematopoietic cells after integrative modeling with IDEAS. A. Heatmap of the emission frequencies of each of the 27 states discovered by IDEAS, with state number and function-associated labels. Each letter in the label indicates a function associated with the combination of features in each state, defined in the box. The indicator for transcribed is H3K36me3, active is H3K27ac, enhancer-like is H3K4me1>H3K4me3, promoter-like is H3K4me3>H3K4me1, heterochromatin is H3K9me3, and polycomb is H3K27me3. B. Example of normalized epigenetic data from ERY in fetal liver around the Gfi1b locus, covering 70kb from chr2:28,565,001-28,635,000 in GRCm38/mm10, used as input to IDEAS for segmentation. The signal tracks are colored distinctively for each feature, with the inferred epigenetic states shown on the last track. The upper limit for signal in each normalized track is given at the right. C. Bar graphs of the average coverage of genomes by each state. The top graph emphasizes the high abundance of state Q, and the second graph shows the abundances of the 26 non-quiescent states. The key for annotated colors is the same order as the states in the bar graph. D. Segmentation pattern across cell types around the Gfi1b exemplar locus. Signal tracks for EP300 (ENCSR982LJQ, ENCODE consortium) and CTCF from mouse fetal liver were included for validation and confirmation, along with the locations of enhancers shown to be active (Enh_vald; Moignard et al. 2013).



## References
Xiang, G., Keller, C. A., Heuston, E., Giardine, B. M., An, L., Wixom, A. Q., ... & Lichtenberg, J. (2020). An integrative view of the regulatory and transcriptional landscapes in mouse hematopoiesis. bioRxiv, 731729.

Zhang, Y., An, L., Yue, F., & Hardison, R. C. (2016). Jointly characterizing epigenetic dynamics across multiple human cell types. Nucleic acids research, 44(14), 6721-6731.

Xiang, G., Keller, C. A., Giardine, B., An, L., Li, Q., Zhang, Y., & Hardison, R. C. (2019). S3norm: simultaneous normalization of sequencing depth and signal-to-noise ratio in epigenomic data. bioRxiv, 506634.

Zhang, Y., & Mahony, S. (2019). Direct prediction of regulatory elements from partial data without imputation. PLoS computational biology, 15(11).


