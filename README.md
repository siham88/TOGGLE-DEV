TOGGLE : Toolbox for generic NGS analyses  - TEST VERSION FOR ONTHEFLY creation of PIPELINES
===========

TOGGLE (TOolbox for Generic nGs anaLysEs) is a suite of 10 packages and more than 110 modules able to manage a large set of NGS softwares
and utilities to easily design pipelines able to handle hundreds of samples. Moreover, TOGGLE offers an easy way to manipulate the various
options of the different softwares through the pipelines in using a single basic configuration file, that can be changed for each assay without
having to change the code itself.

We present also the implementation of TOGGLE in a complete analysis pipeline designed for SNP discovery for large sets of NGS data, ready to use
in different environments (single machine to HPC clusters).


##  Contributing

* Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
* Intellectual property belongs to IRD, CIRAD, ADNid and SouthGreen development platform
* Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Enrique Ortega-Abboud, Sébastien Ravel, Julie Orjuela-Bouniol, Souhila Amanzougarene, Gauthier Sarah, Marilyne Summo, and Francois Sabot
* Copyright 2014-2015

##  Citation
**TOGGLE: Toolbox for generic NGS analyses**. Cécile Monat, Christine Tranchant-Dubreuil, Ayité Kougbeadjo, Cédric Farcy, Enrique
Ortega-Abboud, Souhila Amanzougarene, Sébastien Ravel, Mawussé Agbessi, Julie Orjuela-Bouniol, Maryline Summo and François Sabot. In press, *BMC Bioinformatics*.

##  INSTALLATION

see https://github.com/SouthGreenPlatform/TOGGLE-DEV/blob/onTheFlyv2/INSTALL.md

## MANUAL

see https://github.com/SouthGreenPlatform/TOGGLE-DEV/blob/onTheFlyv2/MANUAL.md

## REQUIREMENTS

#### Perl


* [Data::Translate](http://search.cpan.org/~davieira/Data_Translate-0.3/Translate.pm)
* [Data::Dumper](http://search.cpan.org/~smueller/Data-Dumper-2.154/Dumper.pm)
* [Test::More](http://search.cpan.org/~exodist/Test-Simple-1.001014/lib/Test/More.pm)
* [Test::Deep](http://search.cpan.org/~rjbs/Test-Deep-0.119/lib/Test/Deep.pm)
* [Capture::Tiny](http://search.cpan.org/~dagolden/Capture-Tiny-0.30/lib/Capture/Tiny.pm)
* [List::Compare](http://search.cpan.org/~jkeenan/List-Compare-0.53/lib/List/Compare.pm)
* Switch


#### Bioinformatics software (minimal version)

* [java 1.7](https://www.java.com/fr/)
* [BWA 0.7.2](http://bio-bwa.sourceforge.net/)
* [SAMtools 0.1.18](http://samtools.sourceforge.net/)
* [picardTools 1.63](http://broadinstitute.github.io/picard/)
* [gatk 3.3](https://www.broadinstitute.org/gatk/)
* [fastQC v0.10.1](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [cutadapt 1.2.1](https://pypi.python.org/pypi/cutadapt)
* [FastxToolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [Tophat](https://ccb.jhu.edu/software/tophat/index.shtml)
* [Snpeff](http://snpeff.sourceforge.net/)

#### Bioinformatics tools included

#####BWA (http://bio-bwa.sourceforge.net/)

- bwaAln
- bwaSampe
- bwaSamse
- bwaIndex
- bwaMem
      
#####SamTools (http://samtools.sourceforge.net/)

- samToolsFaidx
- samToolsIndex
- samToolsView
- samToolsSort
- mergeHeader
- samToolsMerge
- samToolsIdxstats
- samToolsDepth
- samToolsFlagstat
       
#####PicardTools (http://broadinstitute.github.io/picard/)

- picardToolsMarkDuplicates
- picardToolsCreateSequenceDictionary
- picardToolsSortSam
        
#####Gatk (https://www.broadinstitute.org/gatk/)

- gatkBaseRecalibrator
- gatkRealignerTargetCreator
- gatkIndelRealigner
- gatkHaplotypeCaller
- gatkSelectVariants
- gatkVariantFiltration
- gatkReadBackedPhasing
    
#####Fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- fastqc

#####FastxToolkit (http://hannonlab.cshl.edu/fastx_toolkit/)

- fastxTrimmer

#####Tophat (https://ccb.jhu.edu/software/tophat/index.shtml)

- bowtiebuild
- bowtie2build
- tophat2

#####Snpeff (http://snpeff.sourceforge.net/)

- snpeffAnnotation
    
#####Cutadapt (https://pypi.python.org/pypi/cutadapt)

- cutadapt

#### OPTIONAL
* Graphviz v2.xx (www.graphviz.org/)


##  Versions Notes - onTheFlyv2 branch

Release 0.3, xxxxxxx

Third beta version release
