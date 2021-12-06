# NanoDam_analysis

This repository comprises R markdowns & python scripts accompanying the preview article:

**Tang JLY, Hakes AE, Krautz R, Suzuki T, Contreras EG, Fox PM, Brand AH**
*NanoDam identifies novel temporal transcription factors conserved between the Drosophila central brain and visual system*
doi: [https://doi.org/10.1101/2021.06.07.447332](https://doi.org/10.1101/2021.06.07.447332)

## Running damMer

The suite of damMer python3 scripts automates the application of the damidseq_pipeline across a multitude of TaDa- or NanoDam-samples by performing all pairwise comparisons between (1.) one 'Dam-fusion'- (TaDa) or one 'NanoDam-tagged protein'-sample (NanoDam) and (2.) one 'Dam-only'- (TaDa) or one 'NanoDam-only'-sample (NanoDam). All damMer parts expect the 'slurm workload manager' on the local system as scripts will be submitted as individual slurm jobs for efficient parallelisation. ![See schematic overview.](https://github.com/robertkrautz/NanoDam_analysis/blob/master/20180123_workflow_damMer_v2.pdf)

### [1.] 'damMer.py'

'damMer' requires [samtools](http://www.htslib.org/), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), the damidseq_pipeline (vR.1), a GATC-fragment file (i.e., 'gff'-format, generated with '[gatc.track.maker.pl](https://github.com/AHBrand-Lab/DamID_scripts)'), and a bowtie2 index for the genome corresponding to the analysed samples (see arguments). 'damMer.py' attempts to identify the location of samtools, bowtie2, the damidseq_pipeline (vR.1) on the system itself and will suggest alternatives, if present.

Names of all '\*.fastq.gz'-files need to be provided as filename lists corresponding to 'Dam-fusion'- ('--experiment') and 'Dam-only'-samples ('--control'). They should simply be added one after the other or provided as a shell array.

#### [1.1.] 'damMer.py' usage
```
dam=($(find . -type f -iname "*.fastq.gz" -and -iname "dam_*"))
exp=($(find . -type f -iname "*.fastq.gz" -and -iname "experiment_*"))
python3 damMer.py -e ${exp[@]} -c ${dam[@]} -i /path/to/index -g /path/to/file.GATC.gff -b /path/to/bowtie2 -s /path/to/samtools -q /path/to/damidseq_pipeline_vR.1.pl
```
#### [1.2.] 'damMer.py' arguments
```
-e / --experiment  List of experimental '*.fastq.gz'-files for the Dam-fusion samples.
-c / --control     List of control '*.fastq.gz'-files for the Dam-only samples.
-i / --index       'bowtie2_build'-derived genome index.
-g / --gatcfrag     '*.GATC.gff'-file listing coordinates of GATC-fragments.
-b / --bow2dir       Path to bowtie2 executables.
-s / --samdir        Path to samtools executables.
-q / --damidseq      Path to damidseq_pipeline executable.
-f / --feedback    Complete mail address to receive slurm feedback.
-d / --defaults    Load defaults for species of interest.
```

#### [1.3.] 'damMer.py' output

For every pairwise comparison, one subdirectory with the results from the damidseq_pipeline_vR.1.pl will be created, including individual '\*.bedgraph'-files enlisting the normalized, genomewide binding intensities of the DNA/chromatin-binding protein of interest. As filenames will be changed while running 'damMer_tracks.py', it is recommended to not change them manually. All arguments & parameters will be logged in a '*.log'-file. All shell scripts submitted by 'damMer.py' are kept for full transparency, too.

### [2.] 'damMer_tracks.py'

The second part of the workflow - 'damMer_tracks.py' - ensures successful and complete execution of 'damMer.py' before streamlining the filenames in the subdirectories. All '\*.bedgraph'-files will be copied into a separate subdirectory, quantile normalized to each other, averaged and converted into '\*.bigwig'-files. In addition, MACS2-dependent peak calling for all pairwise comparisons will also be initialised. 'damMer_tracks.py' requires [bedGraphToBigWig](https://www.encodeproject.org/software/bedgraphtobigwig/), [MACS2](https://pypi.org/project/MACS2/), as well as '[quantile_norm_bedgraph.pl](https://github.com/AHBrand-Lab/DamID_scripts)', ['average_tracks.pl'(https://github.com/AHBrand-Lab/DamID_scripts)] and a file enlisting the chromosome sizes (e.g., [dm6.chrom.sizes](https://www.encodeproject.org/files/dm6.chrom.sizes/)) for downstream processing.

Names of all subdirectories to be included in downstream analysis need to be provided ('--repos') either one after the other or as a shell array. A prefix for the track output folders needs to be specified ('--out') as well as common strings in the names of the 'Dam-fusion'- ('--exppre', e.g., gene symbol of the DNA/chromatin binding protein) and 'Dam-only'-samples ('--ctrlpre', e.g., 'Dam', 'ctrl').

#### [2.1.] 'damMer.py' usage
```
dirs=($(find . -type f -iname "*_vs_*"))
python3 damMer_tracks.py -r ${dirs[@]} -o *output_folder_name* -p *Dam_fusion_protein* -c *Dam* -m /path/to/MACS2 -n /path/to/quantile_norm_bedgraph.pl -a /path/to/average_tracks.pl -b /path/to/bedGraphToBigWig -l /path/to/*genome*.chrom.sizes
```
#### [2.2.] 'damMer_tracks.py' arguments
```
-r / --repos    List of repositories (i.e., directories).
-o / --out      Directory for output.
-p / --exppre   Common string in experimental samples.
-c / --ctrlpre  Common string in control Dam-samples.
-m / --macs2    Path to 'MACS2'.
-n / --quantile Path to 'quantile_norm_bedgraph.pl'.
-a / --average  Path to 'average_tracks.pl'.
-b / --bgToBw   Path to 'bedGraphToBigWig'.
-l / --chrSize  List of chromsome sizes.
-d / --defaults Load defaults for species of interest.
-f / --feedback Complete mail address to receive slurm feedback.
```

#### [2.3.] 'damMer_tracks.py' output

Two output folders will be generated with names based on the indicated prefix ('--out') preceded by either '\*\_DamOnly\_tracks' or '\*\_tracks'. The latter includes '\*.bedgraph' files copied from all included subdirectories ('--repos', i.e., '\*\_vs\_\*') before and after quantile normalization, an average of all quantile normalized files in '\*.bedgraph'-format as well as '\*.bigwig'-files after conversion of all '\*.bedgraph's. In parallel, all chosen subdirectories will include the results from MACS2, e.g., '\*_peaks.broadPeak'.

#### [3.] 'damMer_peaks.py'

In the third part of the workflow - 'damMer_peaks.py' - the presence of all '\*.broadPeak'-files in the chosen subdirectories is ensured before copying them in separate '\*\_DamOnly\_peaks' and '\*\_peaks'-folders. Peaks will be thresholded according to their __*false discovery rate*__, sorted and merged. At last, reproducible peaks are identified based on their appearance in ≥50% across all initial '\*.broadPeak'-files based on all individual, pairwise comparisons.

The list of directories used in 'damMer_tracks.py' ('--repos') either in consecutive order or as a shell array need to be specified together with a prefix for the output folder ('--out') when invoking 'damMer_peaks.py'.

#### [2.1.] 'damMer__peaks.py' usage
```
dirs=($(find . -type f -iname "*_vs_*"))
python3 damMer_peaks.py -r ${dirs[@]} -o *output_folder_name*
```

#### [3.2.] 'damMer_peaks.py' arguments
```
-r / --repos  List of repositories (i.e., directories).
-o / --out    Directory for output.
```

#### [3.3.] 'damMer_peaks.py' output

Similar to 'damMer_tracks.py', two subdirectories will be generated, named according to the indicated prefix ('--out') preceded by '\*\_DamOnly\_peaks' or '\*\_peaks'. They include copies of the '\*.broadPeak'-files from all chosen subdirectories ('--repos') and their '\*.mergePeak'- and  '\*.reproPeak'-derivatives. For each of the 41 predefined FDR-thresholds (i.e., 0 - 2000; -log10-normalized), all merged peaks are enlisted in the corresponding '\*FDR\*.mergePeak'-files (e.g., '75.mergePeak'; bed-format) and the reproducible peaks, present in ≥50% of all pairwise comparisons, are enlisted in the '\*FDR\*.reproPeak' (e.g., '75.reproPeak'; bedgraph-format). '\*.mergePeak'-files allow discrimination of non-/reproducible peaks by color (i.e., red/blue) when loaded into a Genome Browser.

## R markdowns

All custom R markdowns are based on tidyverse to ensure transparency in the analytic workflows and exceptions are only made when absolutely necessary. Code is written explicitly with e.g., functions preceding R library names (e.g., 'dplyr::pull()') and all markdowns are as self-contained as possible. Crucial objects are saved as '\*.rds' to avoid rerunning time-consuming calculations or to load the objects in subsequent R markdowns.

The set of markdowns should be contained in a separate R project directory with the corresponding '\*.Rproj'-file in the root as this is required for the relative paths specified with the here package. Additional subdirectories, i.e., 'resources', 'results', 'data', should be created before running through the analyses. Processed files from running the suite of damMer scripts or from the NCBI Gene Expression Omnibus (accession code: GSE190210) should be placed in the 'data'-subdirectory.

## [4.] 'genomewide_correlation.Rmd'

Pearson correlation coefficients for pairwise comparisons of the genomewide, binned Chronophage-signal from all individual TaDa and NanoDam libraries (i.e., '\*.bgr'-files) was calculated and visualized as correlation matrix.

## [6.] 'create_annotations.Rmd'

Employs biomaRt to acquire a complete set of annotations for transcription start sites (TSSs) in the _Drosophila melanogaster_ genome belonging to protein-coding genes ("BDGP6.bm.TssBiomart.ProteinCoding").

## [7.] 'cluster_peaks.Rmd'

Binding intensities are aggregated over peaks derived from DichaeteGFP, GrainyheadGFP and EyelessGFP via the custom 'aggregator()'-function. After identifying the optimal number of clusters and an optimisation of the clustering procedure based on silhouette widths as read-out, peaks were clustered via stats::kmeans() into 6 clusters based on the quantified binding intensities of Dichaete, Grainyhead and Eyeless. Peaks were sorted according to their respective silhouettes in the visualizations. For convenient tracking of the peaks and their cluster association in a Genome Browser, a bed-file with cluster-based color coding was created.

## [8.] 'annotate_peaks.Rmd'

Peaks were annotated to the nearest TSS of protein-coding genes on the linear genome and subsetted for peaks associated with genes encoding for transcription factors  (TFs). The latter was based on the curated list of TFs from the [_Drosophila_ Transcription Factor Database](https://www.mrc-lmb.cam.ac.uk/genomes/FlyTF/old_index.html).

## [9.] 'scRNAseq_analysis.Rmd'

The two single cell replicates were first processed individually with [Seurat v2.3.4.](https://cran.r-project.org/src/contrib/Archive/Seurat/) before integrating them via _canonical correlation analysis_ (CCA). Dimensionality reduction (i.e., tSNE, UMAP) and cluster definitions were based on 10 dimensions of the aligned CCA. Identification of cluster-specific transcription factors especially for the Type II Neural Stem Cell-derived Intermediate Progenitors (INPs) preceded an overlap of the functional, NanoDam-derived TFs with the TFs detected to be expressed in INPs.
