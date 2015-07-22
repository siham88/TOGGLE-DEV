#/usr/bin/perl

###################################################################################################################################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
#
# Intellectual property belongs to IRD, CIRAD and South Green developpement plateform 
# Written by Cecile Monat, Ayite Kougbeadjo, Mawusse Agbessi, Christine Tranchant, Marilyne Summo, Cedric Farcy, Francois Sabot
#
###################################################################################################################################

#Will test if onTheFly works correctly

use strict;
use warnings;
use Test::More 'no_plan'; #tests => 19; #Number of tests, to modify if new tests implemented. Can be changed as 'no_plan' instead of tests=>11 .
use Data::Dumper;
use Test::Deep;
use lib qw(../Modules/);

use localConfig;
my $configFile='software.config.txt';

use_ok('onTheFly');
can_ok('onTheFly','checkOrder');
can_ok('onTheFly','generateScript');

use onTheFly;


#######################################
#Creating the IndividuSoft.txt file
#######################################
my $creatingCommand="echo \"onTheFly\nTEST\" > individuSoft.txt";
system($creatingCommand) and die ("ERROR: $0: Cannot create the individuSoft.txt file with the command $creatingCommand \n$!\n");

#######################################
#Cleaning the logs for the test
#######################################
my $cleaningCommand="rm -Rf onTheFly_TEST_log.*";
system($cleaningCommand) and die ("ERROR: $0: Cannot clean the previous log files for this test with the command $cleaningCommand \n$!\n");



#########################################
#Remove the files and directory created by the previous test
#########################################
my $testingDir="../DATA-TEST/onTheFlyTestDir";
$cleaningCommand="rm -Rf ../DATA-TEST/$testingDir"; 
system ("rm -Rf $testingDir") and die ("ERROR: $0 : Cannot remove the previous test directory with the command $cleaningCommand \n$!\n");


########################################
#Creation of test directory
########################################
my $makeDirCom = "mkdir $testingDir";
system ($makeDirCom) and die ("ERROR: $0 : Cannot create the new directory with the command $makeDirCom\n$!\n");

########################################
#checkOrder test
#######################################

#testing the correct rendering
#Adding a configHash
my $hashOrderOk =   {
                        "order"=>   {
                                        "1" => "bwaSampe",
                                        "2" => "samtools view",
                                        "3" => "samtools sort"
                                    }
                    };
is (onTheFly::checkOrder($hashOrderOk),'1','Test for correct pipeline onTheFly::checkOrder');

#testing the uncorrect rendering TEST Ok in dev, but die so cannot be tested...
#Adding a configHash
#my $hashOrderNotOk =   {
#                        "order"=>   {
#                                        "2" => "bwaSampe",
#                                        "1" => "samtools view",
#                                        "3" => "samtools sort"
#                                    }
#                        };
#is (onTheFly::checkOrder($hashOrderNotOk),'0','Test for uncorrect pipeline onTheFly::checkOrder');

#testing for dead-end program beginning
my $hashOrderNAOk =   {
                        "order"=>   {
                                        "1" => "fastqc",
                                        "2" => "bwaSampe",
                                        "3" => "samtools view",
                                        "4" => "samtools sort"
                                    }
                    };
is (onTheFly::checkOrder($hashOrderNAOk),'1','Test for correct pipeline with \'dead-end\' software beginning onTheFly::checkOrder');

#Testing for dead-end program in the middle
my $hashOrderNAOkBis =   {
                        "order"=>   {
                                        "3" => "fastqc",
                                        "1" => "bwaSampe",
                                        "2" => "samtools view",
                                        "4" => "samtools sort"
                                    }
                    };
is (onTheFly::checkOrder($hashOrderNAOkBis),'1','Test for correct pipeline with \'dead-end\' software in the middle onTheFly::checkOrder');

#testing the single program
my $hashOrderSingle =   {
                        "order"=>   {
                                        "3" => "bwaSampe"
                                    }
                    };
is (onTheFly::checkOrder($hashOrderSingle),'1','Test for correct pipeline with a single software onTheFly::checkOrder');

#testing multiple call of the same program
my $hashOrderMultiple =   {
                        "order"=>   {
                                        "1" => "bwaSampe",
                                        "2" => "samtools view",
                                        "3" => "samtools sort",
                                        "4" => "samtools view"
                                    }
                    };
is (onTheFly::checkOrder($hashOrderMultiple),'1','Test for correct pipeline with multiple call of the same software onTheFly::checkOrder');



########################################
#generateScript test
#######################################

#testing the correct rendering
#Adding a configHash
my $hashConf =   {
                        "order"=>   {
                                        "1" => "bwaSampe",
                                        "2" => "samtools view",
                                        "3" => "samtools sort"
                                    }
                    };
is (onTheFly::generateScript($hashConf,"$testingDir/ToggleBzzz.pl"),'1','Test for correct pipeline onTheFly::generateScript');

