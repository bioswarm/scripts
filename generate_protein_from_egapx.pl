#!perl -w

# complete.proteins.faa complete.cds.fna complete.genomic.gff
# remove prefix "gnl|WGS:ZZZZ|"
my $name = shift;
my $pep_infile = "complete.proteins.faa";
my $pep_outfile = "$name.pep.fa";
my $cds_infile = "complete.cds.fna";
my $cds_outfile = "$name.cds.fa";
my $gff_file = "complete.genomic.gff";
my $bed_file = "$name.bed";
my $gff_outfile ="$name.gff";

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
    $name =~ s/gnl\|WGS:ZZZZ\|//;
    my ($locus) = $attr =~ /locus_tag=([^;]+)/;
    my ($rna) = $attr =~ /orig_transcript_id=([^;]+)/;
    my ($product) = $attr =~ /product=([^;]+)/;
    
    if($name && $locus && $product && $rna){
	$rna =~ s/gnl\|WGS:ZZZZ\|//;
	$rna = "rna-$rna";
	$mid{$name} = $rna;
        $anno{$name} = $product;
        $gid{$name} = $locus;
        $bed{$name} = "$pos\t$name";
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
	$tmp = "$name\_$tmp";
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
# gnl|WGS:ZZZZ|Hp_009919-P1 uncharacterized protein Hp_009919
while(<IN>){
  if(/^>(\S+)/){
	my $tmp = $1;
        $tmp =~ s/gnl\|WGS:ZZZZ\|//;
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
	  $keep{$mid{$id}} = 1;
	  $keep{$id} = 1;
  	  print OUT1 ">$id $mid{$id} $anno{$id}\n$pep{$id}\n" if defined $pep{$id};
	  print OUT2 ">$id $mid{$id} $anno{$id}\n$cds{$id}\n" if defined $cds{$id};
	  print OUT3 "$bed{$id}\n" if defined $bed{$id};
	}
	else{ print DUP "$id\t$gid\n"; }
  }
}
close OUT1;
close OUT2;
close OUT3;
close DUP;


open GFF,">",$gff_outfile;
open IN, "<", $gff_file;
while(<IN>){
  next if /^#/;
  my @in = split /\t/;
  if($in[2] eq 'mRNA'){
        if($in[8] =~ m/ID=(\S+?);Parent=(\S+?);[\d\D]+?product=([^;]+)/){
		#print "$1\t$2\t$3\n";
                if(defined $keep{$1}){
                        $n_exon = 0;
                        $n_cds = 0;
                        print GFF "$in[0]\t$in[1]\tgene\t$in[3]\t$in[4]\t$in[5]\t$in[6]\t$in[7]\tID=$2\n";
                        print GFF "$in[0]\t$in[1]\tmRNA\t$in[3]\t$in[4]\t$in[5]\t$in[6]\t$in[7]\tID=$1;Parent=$2\n";
                }
        }
  }
  elsif($in[2] eq 'exon'){
        if($in[8] =~ m/Parent=(\S+?);[\d\D]+?;orig_transcript_id=(\S+)[;\s]/){
                if(defined $keep{$1}){
                        $n_exon ++;
                        print GFF "$in[0]\t$in[1]\texon\t$in[3]\t$in[4]\t$in[5]\t$in[6]\t$in[7]\tID=exon-$1-$n_exon;Parent=$1\n";
                }
        }
  }
  elsif($in[2] eq 'CDS'){
        if($in[8] =~ m/Parent=(\S+?);[\d\D]+?;protein_id=([^;]+)/){
                if(defined $keep{$1}){
                        $n_cds ++;
                        print GFF "$in[0]\t$in[1]\tCDS\t$in[3]\t$in[4]\t$in[5]\t$in[6]\t$in[7]\tID=cds-$1-$n_cds;Parent=$1\n";
                }
        }
  }
}
close IN;
close GFF;
