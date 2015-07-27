#!/usr/bin/env perl


###################################################################################################################################
#
# Copyright 2014-2015 IRD-CIRAD-INRA-ADNid
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform for all versions also for ADNid for v2 and v3 and INRA for v3
# Version 1 written by Cecile Monat, Ayite Kougbeadjo, Christine Tranchant, Cedric Farcy, Mawusse Agbessi, Maryline Summo, and Francois Sabot
# Version 2 written by Cecile Monat, Christine Tranchant, Cedric Farcy, Enrique Ortega-Abboud, Julie Orjuela-Bouniol, Sebastien Ravel, Souhila Amanzougarene, and Francois Sabot
# Version 3 written by Cecile Monat, Christine Tranchant, Cedric Farcy, Maryline Summo, Julie Orjuela-Bouniol, Sebastien Ravel, Gautier Sarah, and Francois Sabot
#
###################################################################################################################################




use strict;
use warnings;
use lib qw(./Modules);
use localConfig;
use Data::Dumper;

use pairing;
use toolbox;
use onTheFly;





##########################################
# recovery of parameters/arguments given when the program is executed
##########################################
my $cmd_line=$0." @ARGV";
my ($nomprog)=$0=~/([^\/]+)$/;
unless ($#ARGV>=0)                                                                                          # if no argument given
{

  print <<"Mesg";

  perldoc $nomprog display the help

Mesg

  exit;
}

my %param = @ARGV;                                                                                          # get the parameters 
if (not defined($param{'-d'}) or not defined($param{'-c'}) or not defined($param{'-r'}) or not defined ($param{'-o'}))
{
  print <<"Mesg";

  ERROR: Parameters -d and -c and -r and -o are required.
  perldoc $nomprog display the help

Mesg
  exit;
}


##########################################
# recovery of initial informations/files
##########################################
my $initialDir = $param{'-d'};        # recovery of the name of the directory to analyse
my $fileConf = $param{'-c'};          # recovery of the name of the software.configuration.txt file
my $refFastaFile = $param{'-r'};      # recovery of the reference file
my $outputDir = $param{'-o'};      # recovery of the output folder

##########################################
# Creation of the output folder
##########################################

if (not -d $outputDir) #The output folder is not existing yet
{
    #creating the folder
    my $createOutputDirCommand = "mkdir -p $outputDir";
    system ("$createOutputDirCommand") and die ("\n$0 cannot create the output folder $outputDir: $!\nExiting...\n");
}

chdir $outputDir;

my $infosFile = "individuSoft.txt";

my ($sec, $min, $h, $mois_jour, $mois, $an, $sem_jour, $cal_jour, $heure_ete) = localtime(time);
$mois+=1;
$mois = $mois < 10 ? $mois = "0".$mois : $mois;
$mois_jour = $mois_jour < 10 ? $mois_jour = "0".$mois_jour : $mois_jour;
$h = $h < 10 ? $h = "0".$h : $h;
$min = $min < 10 ? $min = "0".$min : $min;
$an+=1900;
my $date="$mois_jour-$mois-$an-$h"."_"."$min";


#my $infosFile = "$pathIndividu[1]/individuSoft.txt";
open (F1, ">",$infosFile) or die ("$0 : open error of $infosFile .... $!\n");
print F1 "GLOBAL\n";
print F1 "ANALYSIS_$date\n";

toolbox::exportLog("#########################################\nINFOS: Global analysis \n#########################################\n",1);
toolbox::exportLog("----------------------------------------",1);
toolbox::exportLog("INFOS: $0 : Command line : $cmd_line\n",1);
toolbox::exportLog("INFOS: Your output folder is $outputDir\n",1);
toolbox::exportLog("----------------------------------------",1);
toolbox::checkFile($fileConf);                              # check if this file exists
toolbox::existsDir($initialDir);                            # check if this directory exists
toolbox::checkFile($refFastaFile);                          #Retriving the configuration

my $configInfo=toolbox::readFileConf($fileConf);


#Verifying the correct ordering for the experiment, based on input output files and recovering the last step value
my $lastStep = onTheFly::checkOrder($configInfo);

#Generation of Index required for the analysis to work (on the reference only)
toolbox::exportLog("#########################################\nINFOS: Generating reference index requested \n#########################################\n",1);
toolbox::exportLog("----------------------------------------",1);
onTheFly::indexCreator($configInfo,$refFastaFile);

#Generate script

my $script = "toggleBzz.pl";

onTheFly::generateScript($configInfo,$script);


##########################################
# Transferring data in the output folder and organizing
#########################################

#Checking if inly regular files are present in the initial directory
my $initialDirContent=toolbox::readDir($initialDir);

my $initialDirFolder=toolbox::checkInitialDirContent($initialDir);

if ($initialDirFolder != 0)#The initial dir contains subdirectories, so dying
{
    toolbox::exportLog("ERROR : $0 : The initial data directory $initialDir contains subdirectories and not only regular files.\n",0);
}

#Checking input data homogeneity

my $previousExtension=0;
foreach my $file (@{$initialDirContent})
{
    $file =~ m/\.(\w+)$/;
    my $extension = $1;
    if ($extension eq "fq" or $extension eq "fastq") #homogeneisation for fastq extension
    {
        $extension = "fastq";
    }
    if ($previousExtension == 0) #first round
    {
        $previousExtension = $extension;
        next;
    }
    if ($previousExtension ne $extension) #not the same extension
    {
        toolbox::exportLog("ERROR : $0 : The file type in the initial directory are not homogeneous : $previousExtension and $extension are not compatible in the same analysis.\n",0);
    }
    next;
}


#Linking the original data to the output dir

my $workingDir = $outputDir."/Results";
toolbox::makeDir($workingDir);

foreach my $file (@{$initialDirContent})
{    
    my ($shortName)=toolbox::extractPath($file);
    my $lnCommand = "ln -s $file $workingDir/$shortName";
    ##DEBUG print $lnCommand,"\n";
    
    if(toolbox::run($lnCommand)==1)       #Execute command
    {
        toolbox::exportLog("INFOS: $0 : Transferring $file to $workingDir\n",1);
    }
}
toolbox::exportLog("----------------------------------------",1);


if ($previousExtension eq "fastq")               # the data are all in FASTQ format
{
    #########################################
    # recognition of pairs of files and create a folder for each pair
    #########################################
    my $pairsInfos = pairing::pairRecognition($workingDir);            # from files fasta recognition of paired files
    pairing::createDirPerCouple($pairsInfos,$workingDir);              # from infos of pairs, construction of the pair folder
    
    my $listOfFiles = toolbox::readDir($workingDir);                     # read it to recover files in it
    toolbox::exportLog("INFOS: $0 : toolbox::readDir : $workingDir after create dir per couple: @$listOfFiles\n",1);
    
}

#Other Data are not always treated singlely, but can work together => check if order hash steps higher than 1000 using the $lastStep value





#exit;

#########################################
# Launching the generated script on all subfolders
#########################################

my $listSamples=toolbox::readDir($workingDir);

foreach my $currentDir(@{$listSamples})
{
    next unless $currentDir =~ m/:$/; # Will work only on folders
    $currentDir =~ s/:$//;
    my $launcherCommand="$script -d $currentDir -c $fileConf -r $refFastaFile";

    if(toolbox::run($launcherCommand)==1)       #Execute command
    {
        toolbox::exportLog("INFOS: $0 : Correctly launched $script\n",1);
    }
}
exit;