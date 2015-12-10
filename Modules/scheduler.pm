package scheduler;


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

use localConfig;

use toolbox;
use Switch;

##################################
#
#LAUNCHING
#
##################################

our ($commandLine, $requirement, $sample, $configInfo);

sub launcher { #Global function for launching, will recover the command to be launch and will re-dispatch to normal or other scheduler
    
    ($commandLine,$requirement, $sample, $configInfo) = @_;
    
    my $hashCapability = &checkingCapability;
    my $runOutput;
    switch (1)
    {
        case ($hashCapability->{"sge"} && $configInfo->{"sge"}){$runOutput = &sgeRun} #For SGE running
        case ($hashCapability->{"slurm"} && $configInfo->{"slurm"}){$runOutput = &slurmRun} #For SLURM running
        
        #If no scheduler available or configurated in the config info file, let's run it in a normal way
        else {$runOutput = &normalRun};
    }
    
    return $runOutput;
}

sub checkingCapability { #Will test the capacity of launching using various schedulers on the current system
    
    my %capabilityValue;
    
    #SGE test
    $capabilityValue{"sge"} = `qsub -help 2>/dev/null | grep usage`; #Will provide a not-empty output if SGE is installed
    chomp $capabilityValue{"sge"};
    
    #SLURM test
    $capabilityValue{"slurm"}=""; #Will provide a not-empty output if SLURM is installed
    chomp $capabilityValue{"slurm"};
    
    
    #Returning infos as a reference
    return \%capabilityValue;
}


sub normalRun { #For running in normal mode, ie no scheduler
    
    #BASED ON TOOLBOX::RUN, but will not stop the whole pipeline for an error
    use Capture::Tiny qw(capture);
        
    exportLog("INFOS: scheduler::normalRun : $commandLine\n",1);
    
    ##Execute the command
    my ($result,$stderr)=capture {` $commandLine `};
    
    ##Log export according to the error
    if ($?==0)
    {
	##DEBUG
	exportLog("INFOS: scheduler::normalRun on $sample: ".$result."\n--".$stderr."\n",1);
	return 1;
    }
    else
    {
	##DEBUG
	exportLog("WARNING: scheduler::normalRun on $sample: ".$result."\n--".$stderr."\nThe $sample data set has provoked and error, and was not analyzed anymore.\n",2);
	return 2;
    }
    
}

sub sgeRun { #For SGE cluster, running using qsub
    
    my $sgeOptionsHash=toolbox::extractHashSoft($configInfo,"sge");
    my $sgeOptions=toolbox::extractOptions($sgeOptionsHash);
}

sub slurmRun{ #for SLURM cluster, running using sbatch
    
    my $slurmOptionsHash=toolbox::extractHashSoft($configInfo,"slurm");
    my $slurmOptions=toolbox::extractOptions($slurmOptionsHash);
    
}

##################################
#
#WAITING for schedulers ONLY!
#
##################################

sub sgeWait {
    
}

sub slurmWait{
    
}

1;
