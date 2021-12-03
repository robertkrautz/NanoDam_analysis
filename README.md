# NanoDam_analysis

Rmarkdowns & python scripts accompanying the preview article:

**Jocelyn L.Y. Tang, Anna E. Hakes, Robert Krautz, Takumi Suzuki, Esteban G. Contreras, Paul M. Fox, Andrea H. Brand**
*NanoDam identifies novel temporal transcription factors conserved between the Drosophila central brain and visual system*
doi: [https://doi.org/10.1101/2021.06.07.447332](https://doi.org/10.1101/2021.06.07.447332)

## Running damMer

The suite of damMer python3 scripts automates the application of the damidseq_pipeline across a multitude of TaDa- or NanoDam-samples by performing all pairwise comparisons between (1.) one 'Dam-fusion'- (TaDa) or one 'NanoDam-tagged protein'-sample (NanoDam) and (2.) one 'Dam-only'- (TaDa) or one 'NanoDam-only'-sample (NanoDam). All damMer parts expect the 'slurm workload manager' on the local system as scripts will be submitted as individual slurm jobs for efficient parallelisation. See schematic overview further down.

'damMer' requires samtools, bowtie2, the damidseq_pipeline (vR.1), a GATC-fragment file (i.e., 'gff'-format, generated with 'gatc.track.maker.pl'), and a bowtie2 index for the genome corresponding to the analysed samples (see arguments).

#### 'damMer.py' arguments
```
-e , --experiment  List of experimental '*.fastq.gz'-files for the Dam-fusion samples.
-c , --control     List of control '*.fastq.gz'-files for the Dam-only samples.
-d , --defaults    Load defaults for species of interest.
-f , --feedback    Complete mail address to receive slurm feedback.
-i , --index       'bowtie2_build'-derived genome index.
-g , --gatcfrag     '*.GATC.gff'-file listing coordinates of GATC-fragments.
-b, --bow2dir       Path to bowtie2 executables.
-s, --samdir        Path to samtools executables.
-q, --damidseq      Path to damidseq_pipeline executable.
```

```
dam=($(find . -type f -iname "*.fastq.gz" -and -iname "dam_*"))
exp=($(find . -type f -iname "*.fastq.gz" -and -iname "experiment_*"))
python3 ~/Desktop/damMer.py -e "${exp[@]}" -c "${dam[@]}"
```
