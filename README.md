## Python
scRNA_practice.ipynb

This is personal practice of scRNA analysis following steps and code in https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb. The dataset can be downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92332&format=file. Some steps might be modified compared to the original ones. Starting from batch correction, sub-sampled dataset (1/4 of the original size) is used so that it will not take forever to run the code on an old personal computer... More detailed explanations for each step can be referred to the original notebook above.

## R
As a comparison, scRNA_practice.Rmd analyzes the same dataset using Seurat in R. It includes steps from QC to cluster annotation. The dataset is not sub-sampled here. scRNA_practice.nb.html is HTML document from R Markdown. To view plots directly from GitHub, scRNA_practice.md and folder scRNA_practice_files were uploaded.

## References
1. [Luecken, M.D. and Theis, F.J., 2019. Current best practices in single‐cell RNA‐seq analysis: a tutorial. Molecular systems biology, 15(6), p.e8746.](https://www.embopress.org/doi/full/10.15252/msb.20188746)

2. [Haber, A.L., Biton, M., Rogel, N., Herbst, R.H., Shekhar, K., Smillie, C., Burgin, G., Delorey, T.M., Howitt, M.R., Katz, Y. and Tirosh, I., 2017. A single-cell survey of the small intestinal epithelium. Nature, 551(7680), pp.333-339.](https://www.nature.com/articles/nature24489)

3. [Seurat](https://satijalab.org/seurat/vignettes.html)