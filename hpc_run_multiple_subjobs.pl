#!/usr/bin/perl
use strict;

my $modname="py";
my $numproc=$ARGV[0]; #how many to fork
my $command="/opt/gridware/apps/gcc/python/2.7.8/bin/python j$modname.!!num!!.pro"; #what to run - !!num!! is argument 
my $firstfile=$ARGV[1]; #we need the first file to figure out where to start
my $lastfile=$ARGV[2]; #we take note of the last of the batches so we don't run anything that doesn't exist

my @pids;
my $amchild=0;

for (my $i=0; $i<$numproc; $i++) {
	print "and i = ";
	print $i;
	if (my $pid = fork) { 
	#I am the parent 
	print "parent initiated...\n";
	push @pids, $pid;
} else {
	#I am a child
	$amchild=1;
	print "child initiated...\n";
	my $thiscom=$command;
	my $n=$i+$firstfile;
	print "my n is ";
	print $n;
	if ($n > $lastfile) {exit;}
	$thiscom =~ s/!!num!!/$n/g;
	print "Executing $thiscom\n";
	system($thiscom); #sends output to the right place
	$i=$numproc+1; #forces loop exit
	} 
}

unless ($amchild) {
	foreach (@pids) {waitpid($_,0);}
}



#PBS -M david\@okgb.net
#PBS -m abe

export PY_STARTUP="/users/wilkinda/.pystartup"