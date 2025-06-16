#! perl -w

my $new_name = shift;
my $pep_infile = "protein.faa";
my $pep_outfile = "$new_name.pep.fa";
my $cds_infile = "rna.fna";
my $cds_outfile = "$new_name.cds.fa";
my $gff_file ="genomic.gff";

open IN, "<", $gff_file;
while(<IN>){
  next if /^#/;
  my @in = split /\t/;
  if($in[2] eq 'mRNA'){
	if($in[8] =~ m/ID=(\S+?);Parent=(\S+?);[\d\D]+?product=([\d\D]+?);/){
	  $gid{$1} = $2;
	  $anno{$1} = $3;
          print "$1\t$2\t$3\n";
	}
  }
  elsif($in[2] eq 'exon'){
	if($in[8] =~ m/Parent=(\S+?);[\d\D]+?transcript_id=(\S+)[;\s]/){
	  $mid{$2} = $1;
	  $rid{$1} = $2;
          #print "$1\t$2\n";
	}
  }
  elsif($in[2] eq 'CDS'){
	if($in[8] =~ m/Parent=(\S+?);[\d\D]+?protein_id=(\S+)[;\s]/){
	  my $m = (split /;/,$2)[0];
          $mid{$m} = $1;
	  $pid{$1} = $m;
          #print "$1\t$m\n";
	}
  }
}
close IN;

my $flag = 0;
open IN, "<", $cds_infile;
open BUG, ">", "not_used.list";
while(<IN>){
  if(/^>(\S+)/){ 
	if(defined($mid{$1})){ 
	  $id = $mid{$1};
          print "$1\t$id\n";
	  $flag = 1;
	}
	else{ 
	  $flag = 0;
	  print BUG "$1\n"; 
	}
  }
  elsif($flag){ 
    s/\s+//g;
    $cds{$id} .= $_;
  }
}
close IN;

$flag = 0;
open IN, "<", $pep_infile;
while(<IN>){
  if(/^>(\S+)/){ 
	if(defined($mid{$1})){
  	  $id = $mid{$1}; 
	  $flag = 1;
	}
	else{
	  $flag = 0;
	  print BUG "$1\n";
	}
  }
  elsif($flag){
	s/\s+//g;
	$pep{$id} .= $_;
	$len{$gid{$id}}{$id} += length($_) if defined $gid{$id};
  }	 
}
close IN;
close BUG;

open OUT1, ">", $pep_outfile;
open OUT2, ">", $cds_outfile;
open DUP, ">", "duplicate.list";
foreach my $gid (keys %len){
  my $flag = 0;
  foreach my $id (sort {$len{$gid}{$b} <=> $len{$gid}{$a}} keys %{$len{$gid}}){
	$flag ++;
	if($flag == 1){
  	  print OUT1 ">$gid{$id} $pid{$id} $anno{$id}\n$pep{$id}\n" if defined $pep{$id};
	  print OUT2 ">$gid{$id} $pid{$id} $anno{$id}\n$cds{$id}\n" if defined $cds{$id};
	}
	else{ print DUP "$id\t$gid\t$pid{$id}\t$rid{$id}\n"; }
  }
}
close OUT1;
close OUT2;
close DUP;

