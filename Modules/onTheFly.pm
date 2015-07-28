package onTheFly;



###################################################################################################################################
#
# Copyright 2014-2015 IRD-CIRAD
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
# Written by Cecile Monat, Christine Tranchant, Ayite Kougbeadjo, Cedric Farcy, Mawusse Agbessi, Marilyne Summo, and Francois Sabot
#
###################################################################################################################################




use strict;
use warnings;
use Data::Dumper;
use List::Compare;

use lib qw(.);
use localConfig;
use toolbox;
use bwa;
use samTools;
use picardTools;

################################################################################################
# sub checkOrder =>  Will verify the order of the softwares in the pipeline 
#                                           are Ok (ie outfile 1 in Ok as infile 2) 
################################################################################################
# arguments :
# 	- hash of complete configuration
################################################################################################
sub checkOrder
{
    my ($hashConf)=@_;
    my $hashOrder=toolbox::extractHashSoft($hashConf,"order"); #Picking up the options for the order of the pipeline
    
    #Picking up input output for each software
    my $hashInOut=toolbox::readFileConf("$toggle/softwareFormats.txt");
    
    #Verifying the coherence of input/output
    my ($previousSoft,$previousFormat,$currentFormat,$initialStep,$lastStep);
    foreach my $step (sort {$a<=> $b} keys %{$hashOrder})
    {
	my $currentSoft=$$hashOrder{$step};
        $currentSoft =~ s/\d+$//; # Removing number if a software is used more than once with different options
	##DEBUG print $previousSoft,"->",$currentSoft,"\n";
	#if first round
	if (!defined $previousFormat && $$hashInOut{$currentSoft}{"OUT"} ne "NA")
	{ 
	    $previousFormat=$$hashInOut{$currentSoft}{"OUT"};
	    $previousSoft=$currentSoft;
            $initialStep=$step unless $initialStep;
            $lastStep = $step;
	    #print "prout\n";
	    next;
	}
	elsif (!defined $previousFormat && $$hashInOut{$currentSoft}{"OUT"} eq "NA")
	{
            $initialStep = $step;
	    ##DEBUG print "pas prout\n";
	    next;
	}
	
	#For other rounds
	$currentFormat=$$hashInOut{$currentSoft}{"IN"};
	
	#Comparison IN/OUT
	my @listCurrentFormat = split (",", $currentFormat);
	my @listPreviousFormat = split (",", $previousFormat);
	
	## DEBUG print "++",@listCurrentFormat,"++",@listPreviousFormat,"\n";
	
	my $compareList = List::Compare->new(\@listCurrentFormat,\@listPreviousFormat);
	my @intersection = $compareList->get_intersection;
	##DEBUG print "**",@intersection,"\n";
	unless (scalar (@intersection)) #No element are Ok...
	{
	    #Print Error
	    toolbox::exportLog("ERROR: onTheFly::checkOrder : The $previousSoft outputs ($previousFormat) are not compatible with the $currentSoft inputs ($currentFormat).\nPlease check your pipeline order.\n",0);
	}

	#Preparing for the next round

	next if ($$hashInOut{$currentSoft}{"OUT"} eq "NA"); #for a brick such as FastQC which is a 'dead-end'
	
	$previousSoft=$currentSoft;
	$previousFormat=$$hashInOut{$currentSoft}{"OUT"};
        $lastStep = $step;
        ##DEBUG print $lastStep,"\n";

    }
    ##DEBUG print $initialStep,"--",$lastStep,"\n";
    return ($initialStep,$lastStep); #Will return the last step number
}

################################################################################################
# sub generateScript =>  will generate scripts on the fly 
################################################################################################
# arguments :
# 	- hash of complete configuration
#       - script name
################################################################################################
sub generateScript
{
    my ($hashOrder,$script,$limit)=@_;
    
    #Picking up input output for each software
    my $hashSoftware=toolbox::readFileConf("$toggle/softwareFormats.txt");
    
    my $catCommand = "cat $toggle/onTheFly/startBlock.txt"; #Will create the preambule for the pipeline code, including paths, use modules, etc...
    
    foreach my $step (sort {$a<=> $b} keys %{$hashOrder})
    {
        my $currentSoft=$$hashOrder{$step}; #Picking up the step name
        $currentSoft =~ s/ /_/g; #Removing extraspace
        $currentSoft =~ s/\d+$//; #Removing numbers if a soft is called more than once
        $catCommand .= " ".$toggle."/onTheFly/".$currentSoft."Block.txt"; #Adding the code block for the corresponding step in the cat command, as all txt files with these blocks are name as softBlock.txt
	
    }
    
    $catCommand .= " $toggle/onTheFly/endBlock.txt > $script && chmod 775 $script"; #Adding the end of the script and rending it executable
    
    ##DEBUG print $catCommand,"\n";
    
    if(toolbox::run($catCommand)==1)       #Execute command
    {
        toolbox::exportLog("INFOS: onTheFly::generateScript : Correctly done\n",1);
        return 1;
    }
    
    else       
    {
        toolbox::exportLog("ERROR: onTheFly::generateScript : The script $script cannot be created using the following command:\n $catCommand\n",0);       # ... and return an error message
        return 0;
    }
}

################################################################################################
# sub indexCreator =>  will create the different index on the reference needed for the analysis
################################################################################################
# arguments :
# 	- hash of complete configuration
#       - reference file
################################################################################################

sub indexCreator
{
    my ($hashConf,$reference)=@_;
    my $hashOrder=toolbox::extractHashSoft($hashConf,"order"); #Picking up the options for the order of the pipeline
    
    my @listConfig = keys %{$hashConf}; #Picking up all softwares with any option declared
    
    foreach my $step (sort {$a <=> $b}  keys %{$hashOrder})
    {
        my $currentSoft = $$hashOrder{$step};
        
        #INDEXING for BWA
        if ($currentSoft =~ m/bwa/i) #Any step involving BWA
        {
            if ($currentSoft eq "bwaIndex") # If the index is expressely asked
            {
                my $softParameters = toolbox::extractHashSoft($hashConf,"bwaIndex");                                  # recovery of specific parameters of bwa index
                bwa::bwaIndex($reference,$softParameters);
            }
            else #We check if the index is present or not
            {
                my $refIndexedFile = $reference.".ann";
                if (-e $refIndexedFile)
                {# The index is already created
                    toolbox::exportLog("INFOS: onTheFly::indexCreator : The reference index for bwa index already exists, skipped...\n",1);
                    next;
                }
                my $softParameters = toolbox::extractHashSoft($hashConf,"bwaIndex");                                  # recovery of specific parameters of bwa index
                bwa::bwaIndex($reference,$softParameters);
            }
        }
        #INDEXING for PICARDTOOLS
        if ($currentSoft eq "picardToolsCreateSequenceDictionary" or $currentSoft =~ m/gatk/i) #Any step involving GATK
        {
            my $dictFileOut=$reference;   # name of dictionary file    
            $dictFileOut =~ s/\.[^\.]*$/.dict/;
            if ($currentSoft eq "picardToolsCreateSequenceDictionary")
            {
                my $softParameters = toolbox::extractHashSoft($hashConf,"picardToolsCreateSequenceDictionary"); # recovery of specific parameters of picardToolsCreateSequenceDictionary
                my $rmCommand = "rm -f $dictFileOut";
                toolbox::run($rmCommand); #Removing of any existing previous dictionary
                picardTools::picardToolsCreateSequenceDictionary($reference,$dictFileOut,$softParameters);
            }
            else #We check if the dict is present or not
            {
                if (-e $dictFileOut)
                {# The index is already created
                    toolbox::exportLog("INFOS: onTheFly::indexCreator : The reference index for picardTools CreateSequenceDictionary already exists, skipped...\n",1);
            
                }
                else
                {
                    my $softParameters = toolbox::extractHashSoft($hashConf,"picardToolsCreateSequenceDictionary");# recovery of specific parameters of picardToolsCreateSequenceDictionary
                    picardTools::picardToolsCreateSequenceDictionary($reference,$dictFileOut,$softParameters);
                }
            }
        }
        
        #INDEXING for SAMTOOLS
        if ($currentSoft eq "samToolsFaidx" or $currentSoft =~ m/gatk/i) #Any step involving GATK
        {
            if ($currentSoft eq "samToolsFaidx")
            {
                my $softParameters = toolbox::extractHashSoft($hashConf,"samToolsFaidx"); # recovery of specific parameters of samToolsFaidx
                samTools::samToolsFaidx($reference);
            }
            else #We check if the dict is present or not
            {
                my $indexFileOut=$reference.".fai";
                ##DEBUG print $reference,"--",$indexFileOut,"\n";
                if (-e $indexFileOut)
                {# The index is already created
                    toolbox::exportLog("INFOS: onTheFly::indexCreator : The reference index for samtools faidx already exists, skipped...\n",1);
                    ##DEBUG print "skipped faidx\n";
                }
                else
                {                
                    ##DEBUG print "samtools faidx\n";
                    my $softParameters = toolbox::extractHashSoft($hashConf,"samToolsFaidx");# recovery of specific parameters of samToolsFaidx
                    samTools::samToolsFaidx($reference);
                }
            }
        }
    }
    
    
    return 1;
}


1;

=head1 NAME

package I<onTheFly> 

=head1 SYNOPSIS

	use onTheFly;


	
=head1 DESCRIPTION



=head2 Functions
