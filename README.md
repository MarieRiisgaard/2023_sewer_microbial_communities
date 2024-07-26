# The Sewer Microbiome: A Missing Link to Understanding Community Composition and Dynamics in Wastewater Treatment Plants

R scripts and shell scripts used in "The Sewer Microbiome: A Missing Link to Understanding Community Composition and Dynamics in Wastewater Treatment Plants"

## Code Structure

### Data

Raw data files for the scripts are available in the 'data' folder.

Raw reads processing was done using the AmpProc 5.0 workflow (https://github.com/eyashiro/AmpProc), which generated ASV sequences. The ASVs were then mapped to full-length ASVs from the MiDAS 5.2 database (Dueholm et al., 2023) using `usearch` with a `-sintax_cutoff` of 0.6. Raw sequencing data can be found in the SRA archive bioprojects: PRJNA946374 and PRJNA1139651.

The study used the MiDAS database for taxonomic classification, which is available at https://www.midasfieldguide.org/guide/downloads.

### Scripts

The main script is the R Markdown file found at *scripts/markdowns/sewer_microbiome_article.Rmd*. This script imports functions from the *scripts/functions* directory to process the data and generate plots.

### Plots

All plots used in the data analysis are supplied in the following directories:

- Plots in the main article: *output/plots/1_results*
- All supplementary plots: *output/plots/2_supplementary_figures*

## License

MIT
