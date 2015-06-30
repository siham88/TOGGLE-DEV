#!/usr/bin/env perl

###################################################################################################################################
#
# Copyright 2014 IRD-CIRAD
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
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform
# Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Marilyne Summo, Julie Orjuela, Sebastien Ravel and Francois Sabot
#
###################################################################################################################################


use strict;
use warnings;
use lib qw(./Modules);
use localConfig;
use Data::Dumper;

use toolbox;
use radseq;

##########################################
# recovery of parameters/arguments given when the program is executed
##########################################
my $cmd_line=$0." @ARGV";
# print "LIGNE COMMANDE $cmd_line\n";
my ($nomprog)=$0=~/([^\/]+)$/;
unless ($#ARGV>=0)                                                                      # if no argument given
{
  print <<"Mesg";

  perldoc $nomprog display the help

Mesg

  exit;
}

my %param = @ARGV;                                                                      # get the parameters 
if (not defined($param{'-d'}) or not defined($param{'-c'}) or not defined($param{'-r'}) or not defined($param{'-k'}))
{
  print <<"Mesg";

  ERROR: Parameters -d or -c or -r or -k are required.
  perldoc $nomprog display the help

Mesg
  exit;
}


##########################################
# recovery of initial informations/files
##########################################
my $initialDir = $param{'-d'};                                                          # recovery of the name of the directory to analyse
my $fileConf = $param{'-c'};                                                            # recovery of the name of the software.configuration.txt file
my $refFastaFile = $param{'-r'};                                                        # recovery of the reference file
my $keyFile = $param {'-k'}; 

my $fileAdaptator = defined($param{'-a'})? $param{'-a'} : "$toggle/adaptator.txt";      # recovery of the adaptator file

my $optionref = toolbox::readFileConf($fileConf);
#print Dumper ($optionref); 
my $softParameters = toolbox::extractHashSoft($optionref, "radseq");
#print Dumper ($softParameters); 
my $options="";
$options=toolbox::extractOptions($softParameters);                                       #Get radseq options given in sofwareConfig file
#print Dumper ($options);

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
print F1 "GBS\n";
print F1 "ANALYSIS_$date\n";

my $outDir=$initialDir;
$outDir=~s/(^.+\/).+\/$/$1/;
$outDir="$outDir/outputRadseq";
mkdir $outDir;
#print "dirOut: $outDir\n";

toolbox::exportLog("#########################################\nINFOS: GBS analysis \n#########################################\n",1);
toolbox::exportLog("----------------------------------------",1);
toolbox::exportLog("INFOS: $0 : Command line : $cmd_line\n",1);
toolbox::exportLog("----------------------------------------",1);
toolbox::checkFile($fileConf);                                                              # check if this file exists
toolbox::existsDir($initialDir);                                                            # check if this directory exists
toolbox::checkFile($refFastaFile);                                                          # check if the reference file exists
toolbox::checkFile($keyFile);                                                               # check if the key file exists
toolbox::checkFile($fileAdaptator);                                                         # check if adaptator file exists

#radseq::splitKeyFile($keyFile, $outDir);
my @laneDirectoriesOK=radseq::processRadtags($keyFile,$initialDir,$options, $outDir);       # run radseq::processRadtags
print Dumper @laneDirectoriesOK;

#for my $i (0 .. $#laneDirectoriesOK)                                                       # foreach directory run globalAnalysisSGE.pl 
#{
#  print "$laneDirectoriesOK[$i]\n";
#}


