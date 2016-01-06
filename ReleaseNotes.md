
# Versions Notes

## Release 0.3, the 2th of December, 2015

### Modules

- Adding the picardTools AddOrReplaceReadGroup CleanSam ValidateSamFile SamFormatConverter
- Adding the GATK UnifiedGenotyper
- Adding the samTools sortSam flagstats depth idxstats
- Adding the SGE and Slurm schedulers

### Functions

- *On the fly* creation of pipeline (see MANUAL.md for a detailled explanation)
- Use of **Graphviz** to generate a visual output of the pipeline structure
- Adding the scheduler aware launching system (automatically)
- Adding a cleaning system to remove chosen intermediate data in order to save hard drive space (see MANUAL.md)
- Use of compressed gzipped files
- Possibility of starting from FASTQ, SAM/BAM or VCF (gzipped or not)
- Possibility of working in relative path
- Providing an output folder in order to not modify the input data (working with symbolic links)
- Single creation of index/dictionary for reference, once per pipeline and no more once per sample, only if they do not exist.
