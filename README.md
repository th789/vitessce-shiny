# vitessce-shiny

[R Shiny app](https://gehlenborglab.shinyapps.io/vitessce-shiny/) demonstrating integrative and interactive visualization of single-cell data with [Vitessce](http://beta.vitessce.io/) (Visual integration tool for exploration of spatial single cell experiments).



## Usage

### Basic demo

The *Basic demo* page shows examples of Vitessce visualizations for sample datasets (for dataset information, see Datasets section below).

### Run analysis

The *Run analysis* page allows for tailored analysis of Vitessce visualizations for sample single-cell datasets and user-uploaded, single-cell datasets.

**Select example dataset:** Upon choosing an example dataset, the analysis can be tailored by filtering the dataset (using min.cells, min.features, and percent.mt) and by customizing the Vitessce visualization (selecting the analyses, summaries, descriptions, and view options to display for the dataset). If dataset filtering conditions result in a dataset that is too large, the Shiny app will disconnect from the <shinyapps.io> server. Occasionally, when toggling between tabs, visualizations may appear distorted. This can be fixed by zooming in and then zooming out of the page.

**Upload dataset: ** Users can also upload a single-cell dataset to analyze in the form of an .rds file containing the single-cell dataset as a SeuratObject. For examples of the file format, see **_Run analysis_ datasets:** (ending in "_full.rds") under Datasets. Upon uploading the single-cell dataset, the analysis can be performed the same way as for example datasets described above. To further customize single-cell data analysis and visualization, see [Vitessce] (http://beta.vitessce.io/).



## Datasets

Sample datasets were obtained from the following sources:
- Peripheral blood mononuclear cells (PBMC) -- 10X Genomics <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k>
- CD4 T cells -- Zheng, G., Terry, J., Belgrader, P. et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8, 14049 (2017).
- CD8 T cells -- Zheng, G., Terry, J., Belgrader, P. et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8, 14049 (2017).
- Non-small cell lung cancer -- 10X Genomics <https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_nsclc_5gex>

Code used to process the raw data can be found in the data-processing.R file. 

The processed data sets can be found at the following links. **_Basic demo_ datasets:** <https://vitessce-export-examples.s3.amazonaws.com/shiny-app/data_pbmc_results.rds>, <https://vitessce-export-examples.s3.amazonaws.com/shiny-app/data_tcellcd4_results.rds>, <https://vitessce-export-examples.s3.amazonaws.com/shiny-app/data_tcellcd8_results.rds>, and <https://vitessce-export-examples.s3.amazonaws.com/shiny-app/data_nsclc_results.rds>. **_Run analysis_ datasets:** <https://vitessce-export-examples.s3.amazonaws.com/shiny-app/data_pbmc_full.rds>, <https://vitessce-export-examples.s3.amazonaws.com/shiny-app/data_tcellcd4_full.rds>, <https://vitessce-export-examples.s3.amazonaws.com/shiny-app/data_tcellcd8_full.rds>, and <https://vitessce-export-examples.s3.amazonaws.com/shiny-app/data_nsclc_full.rds>.



## Using Vitessce in a Shiny app

The Vitessce R Shiny app embeds Vitessce in an R Shiny app (code in app.R). For a simpler example of using Vitessce in an R Shiny app, see [here](https://vitessce.github.io/vitessce-r/articles/shiny.html#shiny-apps-on-remote-servers).



## App deployment

The Vitessce R Shiny app is hosted by the <shinyapps.io> server at <https://gehlenborglab.shinyapps.io/vitessce-shiny/>.


