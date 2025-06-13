#!perl -w

my $name = shift;
my $pep_infile = "protein.faa";
my $pep_outfile = "$name.pep.fa";
my $cds_infile = "cds_from_genomic.fna";
my $cds_outfile = "$name.cds.fa";
my $gff_file = "genomic.gff";
my $bed_file = "$name.bed";

open IN, "<", $gff_file;
while(<IN>){
  chomp;
  next if /^#/;
  my @in = split /\t/;
  warn "Invalid line: $_" unless @in >= 9;
  if($in[2] =~ 'mRNA'){
	if($in[8] =~ m/locus_tag=(\S+?);[\d\D]+?/){
	  $pos = "$in[0]\t$in[3]\t$in[4]\t";
	}
  }
  elsif($in[2] =~ 'CDS'){
    my $attr = $in[8];
    my ($name) = $attr =~ /Name=([^;]+)/;
    my ($locus) = $attr =~ /locus_tag=([^;]+)/;
    my ($product) = $attr =~ /product=([^;]+)/;
    
    if($name && $locus && $product){
        $anno{$name} = $product;
        $gid{$name} = $locus;
        $bed{$name} = "$pos\t$name";
	#print "$name\t$locus\t$product\n";
    }
  }
}
close IN;

my $flag = 0;
open IN, "<", $cds_infile;
open BUG, ">", "not_used.list";
while(<IN>){
  if(/^>(\S+)/){
	my $tmp = (split /_/,$1)[-2];
	if(defined($gid{$tmp})){ 
	  $id = $tmp;
	  $flag = 1;
	}
	else{ 
	  $flag = 0;
	  print BUG "$tmp\n"; 
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
	my $tmp = $1;
	if(defined($gid{$tmp})){
  	  $id = $tmp; 
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
open OUT3, ">", $bed_file;
open DUP, ">", "duplicate.list";
foreach my $gid (keys %len){
  my $flag = 0;
  foreach my $id (sort {$len{$gid}{$b} <=> $len{$gid}{$a}} keys %{$len{$gid}}){
	$flag ++;
	if($flag == 1){
  	  print OUT1 ">$id $anno{$id}\n$pep{$id}\n" if defined $pep{$id};
	  print OUT2 ">$id $anno{$id}\n$cds{$id}\n" if defined $cds{$id};
	  print OUT3 "$bed{$id}\n" if defined $bed{$id};
	}
	else{ print DUP "$id\t$gid\n"; }
  }
}
close OUT1;
close OUT2;
close OUT3;
close DUP;
