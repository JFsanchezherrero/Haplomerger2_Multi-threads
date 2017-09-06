#!/usr/bin/perl  


## notes
## The agp file's coordinates are 0-base; the cabog6.1 uses agp specification v1.1.
## NOTE: the original scf should be adjusted - no leading Ns for scf and ctg sequences.
## NOTE: the original scf file and the gap filling files should be match with each other!!
## NOTE: the ca6 agp (and our g4f) have no scf prefix for the scf_id!!!!!
## NOTE: the original g4f has only 3 types of record, "W", "WN", "N"
## NOTE: gap-finishing sequences from different configuration/dataset are not compatible 
##       and cannot be pulled together. Namely, multiple fctgs should be merged.

use strict; 
use warnings;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 3 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

   This script is to incorporate the gap-filling cotigs into the 
   present genome assemblies.
   At the point, the gapcloser result is used (in *.gc_fctgs format).
   A new scaffold fasta files and a corresponding agp files were created.
   
   Gzip files are supported.
   
   --fctgs=<gc_fctgs_file> - from GF_extract_filling_ctgs.pl 
   --g4f=<gap4fill file> - as indicated (based on agp v1.1)
   --scf=<scf fasta> - should matech the above two files
   --symbol=N - symbol for the gap-filling contig; def=GF_
   
   Outputs include (placed in the same dir of the source file):
   *_gf.fa.gz, 
   *_gf.agp.
    	 
USAGES

print STDERR "\n\n========== Start at " . qx(date);
print STDERR "$0 $Arg_list\n\n";
my $timer = time();

#### ################################
#### reading arguments

my ($in_scf_name,$out_scf_name);
my ($out_agp_name); 
if ($Arg_list =~ m/--scf=(\S+)/) {
	$in_scf_name=$1;
	$in_scf_name=~m/(.+)\.(fasta|fa)(\.gz)*$/;
	$out_scf_name=$1.'_gf.fa.gz';
	$out_agp_name=$1.'_gf.agp';
}else{
	die "no --scf arg.\n";
}

my $in_g4f_name; 
if ($Arg_list =~ m/--g4f=(\S+)/) {
	$in_g4f_name=$1;
}else{
	die "no --g4f arg.\n";
}

my $in_fctgs_name; 
if ($Arg_list =~ m/--fctgs=(\S+)/) {
	$in_fctgs_name=$1;
}else{
	die "no --fctgs arg.\n";
}

my $symbol="GF_"; 
if ($Arg_list =~ m/--symbol=(\S+)/) {
	$symbol=$1;
}

#### ################################
#### reading fasta  
my %scfs;

{
	my $inFH;
	if($in_scf_name=~m/\.gz$/){
		open($inFH,"gunzip -c $in_scf_name |") or die "Can not open gunzip -c $in_scf_name |.\n";
	}else{
		open($inFH,"<$in_scf_name") or die "Can not open $in_scf_name.\n"; 
	}
	
	my ($line,$name,$seq,$is_the_end)=("","","",0);
	
	while($line=<$inFH>){
		last if $line =~ m/^>/;
	}	
	die "Problems with the fasta file; no symbol > found !\n" unless defined($line);
	
	while($is_the_end==0){
		
		($name,$seq)=('','');
		
		if($line =~ m/>\s*(\S+)/){
			$name=$1;
		}
	
		while($line=<$inFH>){
			last if $line =~ m/^>/;
			chomp $line;
			$seq .= $line;
		}
		$is_the_end=1 unless defined($line);
		
		#$name=~s/^scf//; # the ca6 agp have no scf prefix but the scf fasta has, fix it
		$scfs{$name}=$seq;
	}
	close $inFH;	
}

print STDERR "Finished reading scf.\n";

#### ################################
#### reading fasta  
my %fctgs;

{
	my $inFH;
	if($in_fctgs_name=~m/\.gz$/){
		open($inFH,"gunzip -c $in_fctgs_name |") or die "Can not open gunzip -c $in_fctgs_name |.\n";
	}else{
		open($inFH,"<$in_fctgs_name") or die "Can not open $in_fctgs_name.\n"; 
	}
	
	my ($count1,$count2, $count3, $count4, $count5)=(0,0,0, 0, 0);
	while(<$inFH>){
		next if m/^#|^\s+/;
		chomp;
		my @tt=split /\t/,$_;
		#$tt[0]=~s/^scf//; # the ca6 agp have no scf prefix but the scf fasta has, fix it
		$tt[2]='' unless defined($tt[2]);
		$fctgs{$tt[1].'@'.$tt[0]}=$tt[2];
		$count1++;
	}
	close $inFH;
	
	foreach my $tt (keys %fctgs){
		$count2++;
		$count3++ if $fctgs{$tt}=~m/[ATCGatcg]/;
		$count4++ if length($fctgs{$tt})<1;
		$count5++ if $fctgs{$tt}=~m/\d/;
	}	
	print STDERR "A total of $count2 gaps can be finished,\n"; 
	print STDERR "  with $count3 gaps filled with sequences, $count4 with zero length and $count5 with negative length.\n";
}

print STDERR "Finished reading fctgs.\n";

#### ################################
#### reading and filling the g4f data

#the fields in the g4f
#"#agpv1.1; 0-based. ctg_start/end refers to the strand used in the scaffold, not the original strand.\n";
#"#for gaps between cotig, 4N, 5len, 6fragment, 7yes, 8null\n";
#"#0scf_id\t1scf_start\t2scf_end\t3count\t4type\t5ctg_id\t6ctg_start\t7ctg_end\t8ctg_strand\t9gap_size\t10gap_count\t11filling_seq\n";

my %g4f;

{
	my $inFH;
	if($in_g4f_name=~m/\.gz$/){
		open($inFH,"gunzip -c $in_g4f_name |") or die "Can not open gunzip -c $in_g4f_name |.\n";
	}else{
		open($inFH,"<$in_g4f_name") or die "Can not open $in_g4f_name.\n"; 
	}

  my $ptt;
	while(<$inFH>){
		next if m/^#|^\s+/;
		chomp;
		my @tt=split /\t/,$_;
		if($tt[4] eq 'W'){
			die "can not find the corresponding scf.\n" unless exists $scfs{$tt[0]};
			my $seq;
			if($tt[2]-$tt[1]>0){
				$seq = substr($scfs{$tt[0]},$tt[1],$tt[2]-$tt[1]);
				die "can not get the filling seq.\n" unless length($seq)==$tt[2]-$tt[1];
			}else{
				$seq='';
			}
			$tt[11]=$seq;			
		}elsif($tt[4] eq 'N' or $tt[4] eq 'WN'){
			if(exists $fctgs{$tt[10].'@'.$tt[0]}){
			  $tt[11]=$fctgs{$tt[10].'@'.$tt[0]};
				$tt[4]='FN' if $tt[4] eq 'N'; # filled Ngaps
			}else{
				$tt[11]='N' x $tt[9];
			}
		}
		push @{$g4f{$tt[0]}},\@tt;
		$ptt=\@tt;
	}
	close $inFH;
}

print STDERR "Finished reading g4f.\n";

#### ################################
#### adjusting the g4f data

#the fields in the g4f
#"#agpv1.1; 0-based. ctg_start/end refers to the strand used in the scaffold, not the original strand.\n";
#"#for gaps between cotig, 4N, 5len, 6fragment, 7yes, 8null\n";
#"#0scf_id\t1scf_start\t2scf_end\t3count\t4type\t5ctg_id\t6ctg_start\t7ctg_end\t8ctg_strand\t9gap_size\t10gap_count\t11filling_seq\n";

{
  #processing negative gaps
	foreach my $scf_id (keys %g4f){
		my $tt=$g4f{$scf_id};
		my $count=scalar(@$tt);
		for(my $i=0;$i<$count-1;$i++){
			if( $tt->[$i][11]=~m/^-(\d+)$/ ){
			  $tt->[$i][11]='';
				my $eat=$1;
				if(exists($tt->[$i-1][11]) and length($tt->[$i-1][11])>=$eat){
					$tt->[$i-1][11]=substr($tt->[$i-1][11],0,length($tt->[$i-1][11])-$eat);
				}elsif(exists($tt->[$i+1][11]) and length($tt->[$i+1][11])>=$eat){
					$tt->[$i+1][11]=substr($tt->[$i+1][11],0,length($tt->[$i+1][11])-$eat);
				}else{
					die "Negative gap is large than both the adjacent contigs.\n";
				}
				#while($eat>0){
				#  if($k<0){ die "Negative gap is large than the whole preceding seq.\n"; }
				#  if(length($tt->[$k][11])<1){ $k--; next; }
				#  if(length($tt->[$k][11])>=$eat){ $tt->[$k][11]=substr($tt->[$k][11],0,length($tt->[$k][11])-$eat); $eat=0; next; }
				#  $eat-=length($tt->[$k][11]);$tt->[$k][11]=''; $k--;  
				#}				
			}
		}
	}
	
  #remove records of zero length
	foreach my $scf_id (keys %g4f){
		my $tt=$g4f{$scf_id};
		my $count=scalar(@$tt);
		for(my $i=$count-1;$i>=0;$i--){
			if( length($tt->[$i][11])<1 ){
				splice(@$tt,$i,1); # $count--; $i--;	
			}
		}
	}
	
	# merge contigs
	foreach my $scf_id (keys %g4f){
		my $tt=$g4f{$scf_id};
		my $count=scalar(@$tt);
		for(my $i=0;$i<$count-1;$i++){
			if( ($tt->[$i][4] eq 'W' or $tt->[$i][4] eq'WN') 
						and ($tt->[$i+1][4] eq 'W' or $tt->[$i+1][4] eq'WN')
						and ($tt->[$i][5] eq $tt->[$i+1][5]) ){
				$tt->[$i][11] .= $tt->[$i+1][11];
				splice(@$tt,$i+1,1); $count--; $i--;				
			}elsif(($tt->[$i][4] eq 'N' or $tt->[$i][4] eq 'FN') and length($tt->[$i][11]) == 0){
				splice(@$tt,$i,1); $count--; $i--;	
			}
		}
	}
	
	# adjust the coordinates
	foreach my $scf_id (keys %g4f){
		my $tt=$g4f{$scf_id};
		my ($scf_start,$scf_end, $ctg_end, $count)=(0,0,0,0);
		for(my $i=0;$i<scalar(@$tt);$i++){
			$count++;
			$scf_start=$scf_end;
			$scf_end+=length($tt->[$i][11]);
			($tt->[$i][1],$tt->[$i][2],$tt->[$i][3])=($scf_start,$scf_end,$count);
			if($tt->[$i][4] eq "W" or $tt->[$i][4] eq "WN"){
				$tt->[$i][4]="W";
				($tt->[$i][6],$tt->[$i][7])=(0, length($tt->[$i][11]));
			}elsif($tt->[$i][4] eq "FN"){
				$tt->[$i][4]="W";
				$tt->[$i][5]=$symbol.$tt->[$i][0]."_".$tt->[$i][10];
				($tt->[$i][6],$tt->[$i][7],$tt->[$i][8])=(0, length($tt->[$i][11]),'+');
			}
		}
	}	
}

print STDERR "Finished adjusting g4f.\n";

#### ################################
#### output agp
#### 

{
	$"="\t"; #"
	open(my $outFH," >$out_agp_name") or die "can not create  $out_agp_name.\n";
	foreach my $scf_id (keys %g4f){
		my $tt=$g4f{$scf_id};
		for(my $i=0;$i<scalar(@$tt);$i++){
			$tt->[$i][1]++; 
			if($tt->[$i][4] eq 'W'){
				$tt->[$i][6]++;
				print $outFH "@{$tt->[$i]}[0..8]\n";
				$tt->[$i][6]--;
			}else{
				print $outFH "@{$tt->[$i]}[0..7]\n";
			}
			$tt->[$i][1]--;
		}
	}
	close $outFH;
	$"=" "; #"
}

print STDERR "finished agp output.\n";

#### ################################
#### output scaffolds 
#### 

{
	open(my $outFH,"| gzip -c >$out_scf_name") or die "can not create | gzip -c >$out_scf_name.\n";
	foreach my $scf_id (keys %g4f){
		my $tt=$g4f{$scf_id};
		my $seq;
		for(my $i=0;$i<scalar(@$tt);$i++){
			$seq .= $tt->[$i][11];
		}
		print $outFH ">$scf_id\n";
		while(length($seq)>0){
			print $outFH substr($seq,0,60,'')."\n";
		}
	}
	close $outFH;
}

print STDERR "Finished outputing scaffolds.\n";




print STDERR "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";
