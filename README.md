# Unlocking the Epigenetic Switch of Regeneration by Targeted Suppression of DNMT3a in Retinal Ganglion Cells Restores Vision in Adult Mice with Optic Nerve Injury

<br />
Collaboration project with Dong Feng Chen lab on post optic nerve crush RGC regeneration. <br />
Limited regenerative potential of the optic nerve in adult mammals hinders vision restoration after optic nerve trauma or disease. We identified DNA methyltransferase 3a (DNMT3a) as a potent inhibitor of optic nerve regeneration. Selective suppression of DNMT3a in the retinal ganglion cells (RGCs) promoted robust and functional optic nerve regeneration into the brain, restoring vision in adult mice with traumatic nerve injury. Analyses of DNMT3a deficient RGCs with single-nuclei RNA sequencing, genome-wide bisulfite and bulk RNA profiling revealed demethylation and reactivation of axon growth-related gene networks, while simultaneously deactivating the neurodegenerative signals in various subsets of RGCs.

The analysis includes: <br />
1. Datasets processing - standard Seurat processing, [GSEA pathway analysis](https://github.com/mcrewcow/RGC_Dnmt3a_pONC2d_paper/blob/main/ssGSEA_heatmap.R), [Pathways into patterns normalization and 3D](https://github.com/mcrewcow/RGC_Dnmt3a_pONC2d_paper/blob/main/3d_normalized_pathways.R), [Volcano plot](https://github.com/mcrewcow/RGC_Dnmt3a_pONC2d_paper/blob/main/volcano_sc.R)  <br />
2. RGC subtypes annotation and visualization - [Subtypes heatmap](https://github.com/mcrewcow/RGC_Dnmt3a_pONC2d_paper/blob/main/heatmap_subtypes.R) <br />
3. CellChat analysis for predicting cell communications - [Cellchat analysis](https://github.com/mcrewcow/RGC_Dnmt3a_pONC2d_paper/blob/main/cellchat.R) <br />
4. Bisulfite data visualization - [Plots](https://github.com/mcrewcow/RGC_Dnmt3a_pONC2d_paper/blob/main/bisulfite.R) <br />
5. Bulk data analysis - [GSEA](https://github.com/mcrewcow/RGC_Dnmt3a_pONC2d_paper/blob/main/GSEA_bulk.R) <br />

The sequencing datasets are available under the following GEO numbers: <br />
-[Bulk](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229033) <br />
-[Bisulfite](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229034) <br />
-[Single-nuclei RNA seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228627) <br />

If you would like to get the access to the datasets before they get public, or\
If you have any questions regarding the data, analysis or further improvements, here are the contacts: <br />
ekriukov@meei.harvard.edu - Bioinformatics + data availability\
wailydia_tai@meei.harvard.edu - Contributing author\
dongfeng_chen@meei.harvard.edu - PI
