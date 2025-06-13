#! perl -w

my $pep_infile = 'dmel-all-translation-r6.63.fasta';
my $pep_outfile = 'dmel.pep.fa';
my $cds_infile = 'dmel-all-transcript-r6.63.fasta';
my $cds_outfile = 'dmel.cds.fa';

open IN, "<", $cds_infile;
while(<IN>){
  if(/^>(\S+)/){ $id = $1; }
  else{ 
    s/\s+//g;
    $cds{$id} .= $_;
  }
}
close IN;

open IN, "<", $pep_infile;
while(<IN>){
  if(/^>(\S+)/){
    my $id = $1;
    next unless /type\=polypeptide/;
    print $_ unless /length\=(\d+)/;
    my $len = $1;
    print $_ unless /parent\=(\w+)\,(\w+);/;
    my $gid = $1;
    my $tid = $2;
    print $_ unless /name\=(\S+);/;
    my $name = $1;
    $len{$gid}{$id} = $len;
    $name{$id} = $name;
    $gid{$id} = $gid;
    $tid{$id} = $tid;
  }
}
close IN;

foreach my $gid (keys %len){
  my $cnt = 0;
  foreach my $id (sort {$len{$gid}{$b} <=> $len{$gid}{$a}} keys %{$len{$gid}}){
    $cnt ++;
    if($cnt == 1){
      $flag{$id} = 1;
      last;
    }
  }
}

open IN, "<", $pep_infile;
open OUT1, ">", $pep_outfile;
open OUT2, ">", $cds_outfile;
open OUT3, ">", "fly.annotations";
while(<IN>){
  if(/^>(\S+)/){
    $flag = 0;
    my $id = $1;
    if(defined($flag{$id})){
      $flag = 1;
      print OUT1 ">$gid{$id} $name{$id}\n";
      print OUT2 ">$gid{$id} $name{$id}\n$cds{$tid{$id}}\n";
      print OUT3 "$gid{$id}\t$name{$id}\n";
    }
  }
  elsif($flag){ print OUT1 $_; }
}
close IN;
close OUT1;
close OUT2;
close OUT3;
