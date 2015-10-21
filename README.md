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
Ortega-Abboud, Souhila Amanzougarene, Sébastien Ravel, Mawussé Agbessi, Julie Orjuela-Bouniol, Maryline Summo and François Sabot. In Revision for *BMC Bioinformatics*.

##  INSTALLATION

see https://github.com/SouthGreenPlatform/TOGGLE-DEV/blob/onTheFlyv2/INSTALL.md

## MANUAL

see https://github.com/SouthGreenPlatform/TOGGLE-DEV/blob/onTheFlyv2/MANUAL.md

## REQUIREMENTS

#### Perl


* Data::Translate
* Data::Dumper
* Test::More
* Test::Deep
* Capture::Tiny
* List::Compare
* Switch


#### Bioinformatics software (minimal version)

* java 1.7
* fastQC v0.10.1
* cutadapt 1.2.1
* BWA 0.7.2
* gatk 3.3
* picardTools 1.63
* SAMtools 0.1.18

#### Bioinformatics tools included

#####BWA

>bwaAln

>bwaSampe

>bwaSamse

>bwaIndex

>bwaMem
      
#####SamTools

>samToolsFaidx

>samToolsIndex

>samToolsView

>samToolsSort

>mergeHeader

>samToolsMerge

>samToolsIdxstats

>samToolsDepth

>samToolsFlagstat

        
#####PicardTools

>picardToolsMarkDuplicates

>picardToolsCreateSequenceDictionary

>picardToolsSortSam
        
#####Gatk

>gatkBaseRecalibrator

>gatkRealignerTargetCreator

>gatkIndelRealigner

>gatkHaplotypeCaller

>gatkSelectVariants

>gatkVariantFiltration

>gatkReadBackedPhasing
        
#####Fastqc

>fastqc


#####FastxToolkit

>fastxTrimmer

#####Tophat

>bowtiebuild

>bowtie2build

>tophat2

#####Snpeff

>snpeffAnnotation
    
#####Cutadapt

>cutadapt

#### OPTIONAL
* Graphviz v2.xx


##  Versions Notes - onTheFlyv2 branch

Release 0.3, xxxxxxx

Third beta version release
