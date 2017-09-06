#!/usr/bin/perl  

# notes
# 20131101, bug fixed: when there is no available description, 
# just output the name without any tracing space character.

my $Arg_list = join " ", @ARGV;
if (@ARGV < 0 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit; }

	 This script is to rename the scaffold ID of the assembly
	 produced by SSPACE, i.e., to remove the |sizeNNN part.
	 
	 Example: gunzip -c test.fa.gz | rename4sspace.pl >test2.fa
	 
USAGES

while(my $line=<STDIN>){
	if($line=~m/^>scaffold(\d+)/){
		my $tt= '0' x (7-length($1));
		$line='>sc'.$tt.$1."\n";
	}
	print $line;
}