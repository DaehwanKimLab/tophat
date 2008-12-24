#!/usr/bin/perl


# Author: lh3
# Version: 0.1.7

use strict;
use warnings;
use Getopt::Std;

my $usage = qq(
Usage:   fq_all2std.pl <command> <in.txt>

Command: scarf2std      Convert SCARF format to the standard/Sanger FASTQ
         fqint2std      Convert FASTQ-int format to the standard/Sanger FASTQ
         sol2std        Convert Solexa/Illumina FASTQ to the standard FASTQ
         std2sol        Convert standard FASTQ to Solexa/Illumina FASTQ (simplified)
         fa2std         Convert FASTA to the standard FASTQ
         seqprb2std     Convert .seq and .prb files to the standard FASTQ
         fq2fa          Convert various FASTQ-like format to FASTA
         export2sol     Convert Solexa export format to Solexa FASTQ
         export2std     Convert Solexa export format to Sanger FASTQ
         csfa2std       Convert AB SOLiD read format to Sanger FASTQ
         std2qual       Convert standard FASTQ to .seq+.qual
         instruction    Explanation to different format
         example        Show examples of various formats

Note:    Read/quality sequences MUST be presented in one line.
\n);

die($usage) if (@ARGV < 1);

# Solexa->Sanger quality conversion table
my @conv_table;
for (-64..64) {
  $conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
}

# parsing command line
my $cmd = shift;
my %cmd_hash = (scarf2std=>\&scarf2std, fqint2std=>\&fqint2std, sol2std=>\&sol2std, fa2std=>\&fa2std,
				fq2fa=>\&fq2fa, example=>\&example, instruction=>\&instruction, export2sol=>\&export2sol,
				export2std=>\&export2std, csfa2std=>\&csfa2std, seqprb2std=>\&seqprb2std, std2sol=>\&std2sol,
				std2qual=>\&std2qual);
if (defined($cmd_hash{$cmd})) {
  &{$cmd_hash{$cmd}};
} else {
  die("** Unrecognized command $cmd");
}

sub fa2std {
  my %opts = (q=>25);
  getopts('q:', \%opts);
  die("Usage: fq_all2std.pl fa2std [-q $opts{q}] <in.fa>\n") if (-t STDIN && @ARGV == 0);
  my $q = chr($opts{q} + 33);
  while (<>) {
	if (/^>(\S+)/) {
	  print "\@$1\n";
	  $_ = <>;
	  print "$_+\n", $q x (length($_)-1), "\n";
	}
  }
}

sub csfa2std {
  my %opts = (q=>25, Q=>'', l=>0);
  getopts('q:Q:l:', \%opts);
  die("
Usage:   fq_all2std.pl csfa2std [options] <in.csfa>\n
Options: -q INT     default base quality [$opts{q}]
         -Q FILE    quality file [null]
         -l INT     output read length, 0 for auto [$opts{l}]

Note:    For paired-end alignment, Maq requires two sequence files as the
         input. The n-th read in the first file should forms a read pair with
         the n-th read in the second file. However, SOLiD reads may be
         singletons and therefore further prepocessing is needed.
\n") if (-t STDIN && @ARGV == 0);
  my ($fh, $name, $seq);
  my $len = $opts{l};
  my $q = chr($opts{q} + 33);
  if ($opts{Q}) {
	open($fh, $opts{Q}) || die("** fail to open quality file '$opts{Q}'");
  }
  while (1) {
	while (<>) { last if (/^>/); }
	last unless ($_);
	/^>(\S+)/; $name = $1;
	$_ = substr(<>, 2); chomp;
	tr/0123./ACGTN/;
	$seq = $_;
	if ($fh) { # .qual file is available
	  while (<$fh>) { last if (/^>(\S+)/); }
	  /^>(\S+)/;
	  die("** unmatched seq-qual name: '$name' ne '$1'") unless ($1 eq $name);
	  $_ = <$fh>;
	  s/(\s*\d+\s*)/chr(int($1) + 33)/eg;
	  $_ = substr($_, 1);
	} else {
	  $_ = $q x length($seq);
	}
	if ($name =~ /^(\S+)_F\d$/) { # change read name for maq
	  $name = "$1/1";
	} elsif ($name =~ /^(\S+)_R\d$/) {
	  $name = "$1/2";
	}
	if ($len) { # chop the sequence if required
	  $seq = substr($seq, 1, $len);
	  $_ = substr($_, 1, $len);
	}
	print "\@$name\n$seq\n+\n$_\n";
  }
  close($fh) if ($fh);
}

sub fq2fa {
  while (<>) {
	if (/^@(\S+)/) {
	  print ">$1\n";
	  $_ = <>; print;
	  <>; <>;
	}
  }
}

sub scarf2std {
  while (<>) {
	my @t = split(':', $_);
	my $name = join('_', @t[0..4]);
	print "\@$name\n$t[5]\n+\n";
	my $qual = '';
	@t = split(/\s/, $t[6]);
	$qual .= $conv_table[$_+64] for (@t);
	print "$qual\n";
  }
}

sub seqprb2std {
  die("Usage: fq_all2std.pl seqprb2std <in.seq.txt> <in.prb.txt>\n") if (@ARGV != 2);
  my ($fhs, $fhq);
  open($fhs, $ARGV[0]) || die;
  open($fhq, $ARGV[1]) || die;
  while (<$fhs>) {
	my @t = split;
	my $name = join(":", @t[0..3]);
	$t[4] =~ tr/./N/;
	print "\@$name\n$t[4]\n+\n";
	$_ = <$fhq>;
	@t = split;
	my $q = '';
	my $max = -100;
	for (0 .. $#t) {
	  $max = $t[$_] if ($t[$_] > $max);
	  if (($_&0x3) == 3) {
		$q .= $conv_table[$max+64];
		$max = -100;
	  }
	}
	print "$q\n";
  }
  close($fhs); close($fhq);
}

sub export2sol {
  while (<>) {
	chomp;
	my @t = split("\t", $_);
	if ($t[21] eq 'Y') {
	  my $x = (defined($t[7]) && ($t[7] == 1 || $t[7] == 2))? "/$t[7]" : '';
	  $t[0] =~ s/^(SLXA|HWI)-//; $t[0] =~ s/_Human//i; $t[0] =~ s/_PhiX//i; $t[0] =~ s/_R1//;
	  my $rn_head = ($t[0] =~ /(^[A-Z]+\d+_\d+)/)? $1 : ($t[1]? "$t[0]_$t[1]" : $t[0]);
	  print "\@$rn_head:$t[2]:$t[3]:$t[4]:$t[5]$x\n$t[8]\n+\n$t[9]\n";
	}
  }
}

sub export2std {
  while (<>) {
	chomp;
	my @t = split("\t", $_);
	if ($t[21] eq 'Y') {
	  my $x = (defined($t[7]) && ($t[7] == 1 || $t[7] == 2))? "/$t[7]" : '';
	  $t[0] =~ s/^SLXA-//; $t[0] =~ s/_Human//i; $t[0] =~ s/_PhiX//i; $t[0] =~ s/_R1//;
	  my $rn_head = ($t[0] =~ /(^[A-Z]+\d+_\d+)/)? $1 : ($t[1]? "$t[0]_$t[1]" : $t[0]);
	  print "\@$rn_head:$t[2]:$t[3]:$t[4]:$t[5]$x\n$t[8]\n";
	  my @s = split('', $t[9]);
	  my $qual = '';
	  $qual .= $conv_table[ord($_)] for (@s);
	  print "+\n$qual\n";
	}
  }
}

sub fqint2std {
  while (<>) {
	if (/^@/) {
	  print;
	  $_ = <>; print; $_ = <>; $_ = <>;
	  my @t = split;
	  my $qual = '';
	  $qual .= $conv_table[$_+64] for (@t);
	  print "+\n$qual\n";
	}
  }
}

sub sol2std {
  my $max = 0;
  while (<>) {
	if (/^@/) {
	  print;
	  $_ = <>; print; $_ = <>; $_ = <>;
	  my @t = split('', $_);
	  my $qual = '';
	  $qual .= $conv_table[ord($_)] for (@t);
	  print "+\n$qual\n";
	}
  }
}

sub std2sol {
  my $max = 0;
  while (<>) {
	if (/^@/) {
	  print;
	  $_ = <>; print; $_ = <>; $_ = <>;
	  chomp; tr/!-]/@-|/;
	  print "+\n$_\n";
	}
  }
}

sub std2qual {
  die("Usage fq_all2std.pl std2qual <out.prefix> <in.fastq>\n") if (@ARGV == 0);
  my $pre = shift(@ARGV);
  my ($fhs, $fhq);
  open($fhs, ">$pre.seq") || die;
  open($fhq, ">$pre.qual") || die;
  while (<>) {
	s/^@/>/;
	print $fhs $_; print $fhq $_;
	$_ = <>; print $fhs $_; <>;
	$_ = <>; s/([!-~])/" ".(ord($1)-33)/eg; $_ = substr($_, 1); print $fhq $_;
  }
  close($fhs); close($fhq);
}

sub instruction {

  print "
FASTQ format is first used in the Sanger Institute, and therefore
we take the Sanger specification as the standard FASTQ. Although
Solexa/Illumina reads file looks pretty much like the standard
FASTQ, they are different in that the qualities are scaled
differently. In the quality string, if you can see a character
with its ASCII code higher than 90, probably your file is in the
Solexa/Illumina format.

Sometimes we also use an integer, instead of a single character,
to explicitly show the qualities. In that case, negative
qualities indicates that Solexa/Illumina qualities are used.

";

}

sub example {
  my $exam_scarf = '
USI-EAS50_1:4:2:710:120:GTCAAAGTAATAATAGGAGATTTGAGCTATTT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 19 23 23 23 18 23 23 23
USI-EAS50_1:4:2:690:87:GTTTTTTTTTTTCTTTCCATTAATTTCCCTTT:23 23 23 23 23 23 23 23 23 23 23 23 12 23 23 23 23 23 16 23 23 9 18 23 23 23 12 23 18 23 23 23
USI-EAS50_1:4:2:709:32:GAGAAGTCAAACCTGTGTTAGAAATTTTATAC:23 23 23 23 23 23 23 23 20 23 23 23 23 23 23 23 23 23 23 23 23 12 23 18 23 23 23 23 23 23 23 23
USI-EAS50_1:4:2:886:890:GCTTATTTAAAAATTTACTTGGGGTTGTCTTT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23
USI-EAS50_1:4:2:682:91:GGGTTTCTAGACTAAAGGGATTTAACAAGTTT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 20 23 23 23 23 23 23 23 23 23 23 23 18 23 23 23 23
USI-EAS50_1:4:2:663:928:GAATTTGTTTGAAGAGTGTCATGGTCAGATCT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23
';

  my $exam_fqint = '
@4_1_912_360
AAGGGGCTAGAGAAACACGTAATGAAGGGAGGACTC
+4_1_912_360
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 21 40 40 40 40 40 40 40 40 40 26 40 40 14 39 40 40
@4_1_54_483
TAATAAATGTGCTTCCTTGATGCATGTGCTATGATT
+4_1_54_483
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 16 40 40 40 28 40 40 40 40 40 40 16 40 40 5 40 40
@4_1_537_334
ATTGATGATGCTGTGCACCTAGCAAGAAGTTGCATA
+4_1_537_334
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 21 29 40 40 33 40 40 33 40 40 33 31 40 40 40 40 18 26 40 -2
@4_1_920_361
AACGGCACAATCCAGGTTGATGCCTACGGCGGGTAC
+4_1_920_361
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 39 40 40 40 40 40 40 40 40 31 40 40 40 40 40 40 15 5 -1 3
@4_1_784_155
AATGCATGCTTCGAATGGCATTCTCTTCAATCACGA
+4_1_784_155
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 31 40 40 40 40 40
@4_1_595_150
AAAGACGTGGCCAGATGGGTGGCCAAGTGCCCGACT
+4_1_595_150
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 30 40 40 40 40 40 40 40 40 40 20 40 40 40 40 40 14 40 40
';

  my $exam_sol = '
@SLXA-B3_649_FC8437_R1_1_1_610_79
GATGTGCAATACCTTTGTAGAGGAA
+SLXA-B3_649_FC8437_R1_1_1_610_79
YYYYYYYYYYYYYYYYYYWYWYYSU
@SLXA-B3_649_FC8437_R1_1_1_397_389
GGTTTGAGAAAGAGAAATGAGATAA
+SLXA-B3_649_FC8437_R1_1_1_397_389
YYYYYYYYYWYYYYWWYYYWYWYWW
@SLXA-B3_649_FC8437_R1_1_1_850_123
GAGGGTGTTGATCATGATGATGGCG
+SLXA-B3_649_FC8437_R1_1_1_850_123
YYYYYYYYYYYYYWYYWYYSYYYSY
@SLXA-B3_649_FC8437_R1_1_1_362_549
GGAAACAAAGTTTTTCTCAACATAG
+SLXA-B3_649_FC8437_R1_1_1_362_549
YYYYYYYYYYYYYYYYYYWWWWYWY
@SLXA-B3_649_FC8437_R1_1_1_183_714
GTATTATTTAATGGCATACACTCAA
+SLXA-B3_649_FC8437_R1_1_1_183_714
YYYYYYYYYYWYYYYWYWWUWWWQQ
';

  print qq(
solexa
======
$exam_sol
scarf
=====
$exam_scarf
fqint
=====
$exam_fqint
);
}
