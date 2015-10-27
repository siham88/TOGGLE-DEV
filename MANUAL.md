MANUAL for TOGGLEv3 using ONTHEFLY creation of PIPELINES
===========

The *onTheFly* version allows users to create their own customized pipelines.
You can modify not only the different options for each software but also define specific organization for your analysis.

You can therefore remove some steps compared to the previous version, starting from SAM/BAM or VCF files instead of FASTQ only, asking for individual treatments only, individual then common (such as mapping followed by common calling), or even only common treatment.

#Launching an analysis
The current version is based on the **toggleGenerator.pl** script

````
$toggleGenerator.pl -d initialDirectory -r referenceFile -c configurationFile -o outputFolder
````

with
> -d initialDirectory: a folder with raw data to be treated (FASTQ, FASTQ.GZ, SAM, BAM, VCF) 

>>!!! THE FILE NAMES MUST BE UNDER THE FORM individual1_1.fastq OR mapping1.sam !! THE ONLY DOT ('.') ALLOWED IN NAMES IS FOR THE EXTENSION !!

> -r referenceFile: the reference FASTA file to be used. If no index exists it will be created accordingly to the pipeline requested index. If the index exist, they will not be re-created UNLESS the pipeline order (see below) expressively requests it (updating the index e.g.)

> -c configurationFile: generally it is the *software.config.txt* file but it can be any text file structured as shown below.

> -o outputFolder: the current version of TOGGLE will not modify the initial data folder but will create an output directory with all analyses in.

All the locations (reference, config file, folder in and out) can be provided as absolute (/home/mylogin/data/myRef.fasta) or in relative (../data/myRef.fasta).

The *software.config.txt* file is an example of how to customize your pipeline.

Any software configuration will start as follows:
 ````
 $mySoftware
 option1
 option2

 $mySecond soft
 option1
 option2
 ````

# Sending options
As for the previous version, you can address any option to any given software (as soon as the given option exists for the given software ^^) as follows:
````
$bwaAln
-e=1
-n=5
````

You can also write as follows
````
$BWA ALN
-e=1
-n=5
````

The software name is not case sensitive and the subprogram can be "glued" to the name (bwaALN is recognized, as well as bwa aln).

If you plan to use the same software with different configuration in different steps of the pipeline (see below), please specify as follows
````
$samtoolsView 1
OPTIONS

$samtools view 2
OTHEROPTIONS
````


#Creating a pipeline

The **toggleGenerator.pl** script will use the configuration file informations (generally named as *software.config.txt* file, but not mandatory) to create specific pipeline scripts in the output folder called **toggleBzz.pl** for individual treatments (individual mapping e.g.) and **toggleMultiple.pl** for global treatment (global SNP calling e.g.).

The order will be verified before creating the scripts, by checking if the output files from the step n-1 will be accepted as input for the step n.

###Providing an order
The order of a specific pipeline is provided in the *software.config.txt* as the software *order*

Thus, if you provide option such as:
````
$order
1=bwa aln
2=BWASAMPE
3=samtools view
````
You will obtain a pipeline taking fastq files (plain text of gzipped), aligning them upon a given reference using bwa aln/sampe, and finally cleaning the resulting sam and providing a BAM file (if samtools view option are as such).

Note that the way you write the name of the step is not really important, e.g. it is not case-sensitive, and can have space between the different words. A dedicated module will adjust everything.

### Same software repeated multiple times

If you want to adress multiple times the same step BUT with different configurations (e.g. a first samtools view to obtain a BAM then a second samtools view to clean this BAM using the -f option), you have to indicate as follows
````
$order
1=bwa aln
2=BWASAMPE
3=samtools view 1
4= samtools view 2
````

In the same time you have to provide the same informations in your configuration:
````
$samtoolsView1
-Sb

$samtools View 2
-f 0x02
````

###Giving a common step to all individuals

If you want for instance to map all individuals separately then perform a common SNP calling, please add a step number higher than 999.

````
$order
1= bwaAln
2=bwaSampe
3=picardTools SortSam
1000=gatkUnifiedGenotyper
````
This pipeline will map FASTQ then sort per coordinates the SAM for each individuals, and then will launch a common step for all as a global SNP calling. You can add after this calling other steps (1001=gatkVariantFiltrator for example).

###Starting only as a common treatment
If you want only a global treatment (you may have all your formatted BAM and you want to perform the calling and subsequent steps), you can create the following pipeline:
````
$order
1000=gatk UnifiedGenotyper
1001=gatk VariantFiltrator
1002=gatkSelectVariants
````

The pipeline will treat the data as a global pipeline in this case, without separating the individual files.

#Cleaning steps
In order to gain hard drive space, you may want to remove some steps from your complete analysis.

In this case, you can specify in the configuration file which step must be REMOVED along the pipeline (as soon as the step following them has been correctly completed), using the key *cleaner*.

As an example
````
$order
1=bwaAln
2=bwaSampe
3=samtools sort

$cleaner
1
2
````

There only the last step result (samtools sort) will be conserved. The folders and logs of the cleaned steps are conserved, but not the resulting files.

#Drawing the pipeline
If *Graphviz* is installed on the running server, the **toggleGenerator.pl** script will also generate a graphical view of the pipeline in the output folder.
If *Graphviz* is not installed, a .dot file is nevertheless created, allowing the user to create a graphical view on another machine having *Graphviz*
