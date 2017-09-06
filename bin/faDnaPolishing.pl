#!/usr/bin/perl  

# notes
# 

my $Arg_list = join " ", @ARGV;
if (@ARGV < 0 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit; }

	 This a robust script is to do these four tasks stepwise (can be skipped):
	 1. change illegal characters to n;
	 2. if a sequence portion between two n/N 
	    is short than a given length, replace it with n+. 
	 3. remove leading and ending n/N+;
	 4. remove a sequence if short than a given length;
   
   --legalizing - change illegal characters to n; def=no
   --maskShortPortion=N - task 2; def=0
   --noLeadingN - task 3; def=no
   --removeShortSeq=N - task 4; def=0

   Example 1: gunzip -c test.fa.gz | faDnaPolishing.pl --legalizing >test2.fa
	 
USAGES

print STDERR "\n\n========== Start at " . qx(date);
print STDERR "$0 $Arg_list\n\n";
my $timer = time();


my %para;

print STDERR "\n";

$para{'--legalizing'}=0;
if($Arg_list =~ m/--legalizing/){
	$para{'--legalizing'}=1;
}
print STDERR "--legalizing=$para{'--legalizing'} .\n";

$para{'--maskShortPortion'}=0;
if($Arg_list =~ m/--maskShortPortion=(\d+)/){
	$para{'--maskShortPortion'}=$1;
}
print STDERR "--maskShortPortion=$para{'--maskShortPortion'} .\n";

$para{'--noLeadingN'}=0;
if($Arg_list =~ m/--noLeadingN/){
	$para{'--noLeadingN'}=1;
}
print STDERR "--noLeadingN=$para{'--noLeadingN'} .\n";

$para{'--removeShortSeq'}=0;
if($Arg_list =~ m/--removeShortSeq=(\d+)/){
	$para{'--removeShortSeq'}=$1;
}
print STDERR "--removeShortSeq=$para{'--removeShortSeq'} .\n";


#### read and output sequences
my ($line,$name,$desc,$seq,$is_the_end)=("","","","",0);

while($line=<STDIN>){
	last if $line =~ m/^>/;
}	
die "Problems with the fasta file; no symbol > found !\n" unless defined($line);

while($is_the_end==0){
	
	($name,$desc,$seq)=('','','');
	
	if($line =~ m/>\s*(\S+)/){
		$name=$1;
	}
	if($line =~ m/>\s*(\S+)\s+(.+)/){
		$desc=$2;
		chomp $desc;
	}	

	while($line=<STDIN>){
		last if $line =~ m/^>/;
		chomp $line;
		$seq .= $line;
		
	}
	$is_the_end=1 unless defined($line);
	
	print STDERR ">$name: "; #### output processing info
	
	if($para{'--legalizing'}>0){
  	if($seq=~s/[^acgtACGTnN]/n/g){
    	print STDERR "illegal characters found; ";  #### output processing info
  	}
  }
  
	if($para{'--maskShortPortion'}>0){
		my $flag=0;
		while($seq=~m/[nN]([^nN]+)[nN]/g){
			my $pp=pos($seq)-1;			
			if(length($1)<=$para{'--maskShortPortion'}){
				if($flag<1){ $flag=1; print STDERR "masking short portion; "; }  #### output processing info
				my $tt="n" x length($1);
				substr($seq,$-[1],length($1),$tt);
			}
			pos($seq)=$pp; #### restore the pos altered by substr
		}
		if($seq=~m/^([^nN]+)[nN]/){
			if(length($1)<=$para{'--maskShortPortion'}){
				if($flag<1){ $flag=1; print STDERR "masking short portion; "; }  #### output processing info
				my $tt="n" x length($1);
				substr($seq,$-[1],length($1),$tt);
			}
		}
		if($seq=~m/[nN]([^nN]+)$/){
			if(length($1)<=$para{'--maskShortPortion'}){
				if($flag<1){ $flag=1; print STDERR "masking short portion; "; }  #### output processing info
				my $tt="n" x length($1);
				substr($seq,$-[1],length($1),$tt);
			}
		}				
	}

	if($para{'--noLeadingN'}>0){
  	if($seq=~s/^[nN]+//){
    	print STDERR "leading nN found; ";  #### output processing info
  	}
  	if($seq=~s/[nN]+$//){
    	print STDERR "ending nN found; ";  #### output processing info
  	}  	
  }	
	

	if( length($seq) > $para{'--removeShortSeq'} ){
		
		if(length($desc) == 0){		
		  print STDOUT ">$name\n";
		}else{
		  print STDOUT ">$name $desc\n";		
		}

		while(length($seq)>0){
			print STDOUT substr($seq,0,60,'')."\n";
		}
	}else{
		print STDERR "short seq removed; ";  #### output processing info
	}
	
	print STDERR "finished.\n";  #### output processing info
}


print STDERR "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";