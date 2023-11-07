#!/usr/bin/perl 
####################################################
# parse PSIBLAST result file;
# INPUT: The psiblast output result.and the outputlist file. 
# OUPUT: A list of id and evalue, and score.
###################################################
use strict;

if (@ARGV != 3)
{
        die "INPUT: The psiblast output result.and the outputlist file..\n";
}

my $id_format = 30; # like AT1G18790.1  the formate of the ids of the sequence; suggested to be longer than the real id, just not too long to reach the score and evalue part
my $input_file = $ARGV[0];
my $output_file_2 = $ARGV[1];
my $iteration = $ARGV[2];

my @evalues;
my @identities;
my @names;
my $align = "false";
open(READ, "<$input_file");   #this time find the hits;
my $flag = "false";
while(<READ>){
	my $line = $_;
	if($line =~ /Results from round $iteration/){
		$flag = "true";
		
	}
	if ($flag eq "true"){
                my @items = split(/\s+/, $line);
		if($items[6] eq "Expect" && $align eq "true"){			
				$items[8] =~ s/,//;
				push(@evalues, $items[8]);
		}
                if($items[1] eq "Identities" && $align eq "true"){
                                my @nums = split(/\//, $items[3]);
				my $div = $nums[0] / $nums[1];
				push(@identities, $div);
				$align = "false";
                }
                if(substr($line, 0, 1) eq ">"){
				$line =~ s/\n//;
				$line =~ s/>//;
                                push(@names, $line);
				$align = "true";
                }
	}
	if($flag eq "true" && $line =~ /Searching..................................................done/){
		$flag = "false";
	}
}
close(READ);
my $length = @evalues;
if($length == 0){	#"CONVERGED!" before reaching the iteration
	open(READ, "<$input_file");   #this time find the hits;
	open(WRITE, ">$output_file_2");
	my $flag = "false";
	while(<READ>){
        	my $line = $_;
        	if($line =~ /CONVERGED!/){
                	$flag = "true";

        	}
        	if ($flag eq "true"){
                       	        my @items = split(/\s+/, $line);
                               	if($items[6] eq "Expect" && $align eq "true"){
                                       	$items[8] =~ s/,//;
                                       	push(@evalues, $items[8]);
                               	}
                               	if($items[1] eq "Identities" && $align eq "true"){
                                       	my @nums = split(/\//, $items[3]);
                                       	my $div = $nums[0] / $nums[1];
                                       	push(@identities, $div);
                                       	$align = "false";
                               	}
                               	if(substr($line, 0, 1) eq ">"){
                                       	$line =~ s/\n//;
                                       	$line =~ s/>//;
                                       	push(@names, $line);
                                       	$align = "true";
                               	}
               	}
       	}
	close(READ);
}
open(WRITE, ">$output_file_2");
$length = @evalues;
for(my $i = 0; $i < $length; $i++){
	print WRITE "$names[$i] $evalues[$i] $identities[$i]\n";
}
close(WRITE);
