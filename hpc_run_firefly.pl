#!/usr/bin/perl -w
use strict;

#Usage - first argument is number of nodes this instance will run, e.g. 16
#We assume it's consecutive - second argument is the first file number we will 
#start at. So './make_and_submit.pl 16 48' runs files 48 through 64

my $firstarg=$ARGV[0];
my $lastarg=$ARGV[1];
my $modname="py";
my $workdir="/users/goddardd/firefly_original";
chdir $workdir or die "Couldn't go to $workdir";
my $procpernode=3; #cores/node requested

open PYOUT, "> ../jobs/j$modname.$firstarg.$lastarg.py" or die "Couldn't open j$modname.$firstarg.$lastarg.py for writing"; 
print PYOUT "from firefly_job import *\nfirefly_job($firstarg,$lastarg)\nprint 'complete'\nexit()";
close PYOUT;

#Now develop the qsub script

my $qsubscript=<<ENDSCRIPT; 
#!/bin/sh
#PBS -N FFLY-mo-$firstarg-$lastarg
#PBS -l nodes=1:ppn=$procpernode,walltime=300:00:00,pmem=1500mb
#PBS -q sciama1.q

cd $workdir
export OMP_NUM_THREADS=1
/opt/apps/python/2.7.8/gcc-4.4.7/bin/python ../jobs/j$modname.$firstarg.$lastarg.py
ENDSCRIPT
open QSUB,"> ../jobs/$modname-job$firstarg-$lastarg.sh" or die "Couldn't open $modname-job$firstarg-$lastarg.sh for writing"; 
print QSUB $qsubscript;
print "Wrote $modname-job$firstarg-$lastarg.sh\n";
#submit job
print "Submitting...";
system("qsub ../jobs/$modname-job$firstarg-$lastarg.sh");
print "Done\n";

# #PBS -M david\@okgb.net
# #PBS -m abe - place before cd workdir
# 
