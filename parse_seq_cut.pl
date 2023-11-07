#######################################
#parse the sequence file a singel file and save each sequence in one file
#specificly for Sorbi1, as it's ID is too long and waste space.
#Author: Zheng Wang
#######################################

use strict;

my $seq_file = $ARGV[0];
my $output_dir = $ARGV[1];

open(READ, "<$seq_file");

my $counter = 1;
my $sequence;
my $title;
my $first = "true";
while(<READ>){
  my $line = $_;
	$line =~ s/\s+$//;
	$line =~ s/\n//;
	if(length($line) == 0){
		next;
	}
  if(substr($line, 0, 1) eq ">"){
    if($first eq "false"){
      open FASTA, ">$output_dir/$title.fasta" or die $!;
      print FASTA ">".$title."\n";
      print FASTA $sequence;
      close FASTA;
			`chmod 774 $output_dir/$title.fasta`;
      $counter++;
      $sequence = "";
    }
		my @items = split(/\s+/, $line);
		$title = $items[0];	
    $title =~ s/>//;
  }
  if(substr($line, 0, 1) ne ">" && length($line) > 1){
    $sequence = $sequence.$line;
    $sequence =~ s/\n//g;
    $sequence =~ s/ //g;
    $sequence =~ s/^A-Z//g;
    $sequence =~ s/\*//g;
    $first = "false";
  }
}
open(FASTA, ">$output_dir/$title.fasta");
print FASTA ">".$title."\n";
print FASTA $sequence;
close FASTA;
`chmod 774 $output_dir/$title.fasta`;
print "$counter";
