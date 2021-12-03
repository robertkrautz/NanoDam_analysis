# NanoDam_analysis

This repository comprises R markdowns & python scripts accompanying the preview article:

**Tang JLY, Hakes AE, Krautz R, Suzuki T, Contreras EG, Fox PM, Brand AH**
*NanoDam identifies novel temporal transcription factors conserved between the Drosophila central brain and visual system*
doi: [https://doi.org/10.1101/2021.06.07.447332](https://doi.org/10.1101/2021.06.07.447332)

## Running damMer

The suite of damMer python3 scripts automates the application of the damidseq_pipeline across a multitude of TaDa- or NanoDam-samples by performing all pairwise comparisons between (1.) one 'Dam-fusion'- (TaDa) or one 'NanoDam-tagged protein'-sample (NanoDam) and (2.) one 'Dam-only'- (TaDa) or one 'NanoDam-only'-sample (NanoDam). All damMer parts expect the 'slurm workload manager' on the local system as scripts will be submitted as individual slurm jobs for efficient parallelisation. See schematic overview further down.

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

![Overview of damMer scripts.](https://github.com/robertkrautz/NanoDam_analysis/blob/master/20180123_workflow_damMer_v2.pdf)