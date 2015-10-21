#!/usr/bin/perl

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
use Data::Dumper;

use Test::More tests => 7;
use lib qw(../Modules/);


########################################
#use of radseq module ok
########################################
use_ok('toolbox') or exit;                                                                      # Check if toolbox is usable
use_ok('radseq') or exit;                                                                       # Check if radseq is usable
can_ok('radseq','splitKeyFile');                                                                # Check if radseq::splitKeyFile is find
can_ok('radseq','parseDirectory');                                                              # Check if radseq::parseDirectory is find
can_ok('radseq','processRadtags');                                                              # Check if radseq::processRadtags is find

use toolbox;
use radseq;


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"radseq\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0 : Cannot create the individuSoft.txt file with the command $creatingCommand\n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -f radseq_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0 : Cannot remove the previous log files with the command $cleaningCommand \n$!\n");

#########################################
#Remove the files and directory created by the previous test
#########################################
$cleaningCommand="rm -Rf ../DATA-TEST/radseqTestDir";
system($cleaningCommand) and die ("ERROR: $0 : Cannot remove the previous test dir with the command $cleaningCommand \n$!\n");

########################################
#Creation of test directory
########################################
my $testingDir="../DATA-TEST/radseqTestDir";
my $makeDirCom = "mkdir $testingDir";
system ($makeDirCom) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCom\n$!\n");

########################################
#Input keyfile and directory
########################################

#Copying keyFile into testingDir
my $originalKeyFile = "../DATA/expectedData/radseq/keyFileTest.txt";        # Keyfile
my $keyFile = "$testingDir/keyFileTest.txt";                                # Keyfile for test
my $FileCopyCom = "cp $originalKeyFile $keyFile";                           # command to copy the original keyFile into the test directory
system ($FileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalKeyFile in the test directory with the command $FileCopyCom\n$!\n");    # RUN the copy command

#Copying fastq directory into testingDir
my $originalInitialDir = "../DATA/expectedData/radseq/initialDir/";         # initialDir with two fastq files
my $initialDir = "$testingDir/";                                            # initialDir for test
$FileCopyCom = "cp -r $originalInitialDir $initialDir";                     # command to copy the original initialDir into the test directory
system ($FileCopyCom) and die ("ERROR: $0 : Cannot copy the file $originalInitialDir in the test directory with the command $FileCopyCom\n$!\n");    # RUN the copy command
$initialDir = "$testingDir/initialDir/";
my $options="-e apeKI ";                                                             #options de radseq::processRadtags are empty by default
######################


#######################################

### Test of radseq::processRadtags ###
is ((radseq::processRadtags($keyFile, $initialDir, $options)),1, 'radseq::processRadtags');         # TEST IF FONCTION WORKS
my $expectedOutput = `ls ../DATA/expectedData/radseq/outputRadseq/` or die ("ERROR: $0 : Cannot list the directory ../DATA/expectedData/radseq/outputRadseq with the command ls \n$!\n");
my @expectedOutput = split(/\n/,$expectedOutput);

my $observedOutput = `ls $testingDir/outputRadseq/` or die ("ERROR: $0 : Cannot list the directory .$testingDir/outputRadseq/ with the command ls \n$!\n");
my @observedOutput = split(/\n/,$observedOutput);

is_deeply (\@expectedOutput, \@observedOutput, "radseq output checkout");                             # TEST IF THE STRUCTURE OF THE FILE OUT IS GOOD
##############################

exit;