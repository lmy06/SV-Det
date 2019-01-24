#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;

my $usage = <<'USAGE';
Discription:
Author:leimengyue
contact:leimengyue@genomics.cn
Usage
-i	     input directory [require]
-samp        the name of tumor sample    [require]
-bp	     maxium bp distence to consider a overlap [defult:50]
-n	     at least detected in how many tools to consider a SV [defult:2]
-o           output directory [require]

Example
perl a03multi.SV.overlap.pl  -samp tumor.sample -i input directory -o output directory -bp 50 -n 2
USAGE
my ( $samp , $bp, $tools, $indir, $outdir);
GetOptions(
	'samp=s'      => \$samp,
	'bp=n'	      => \$bp,
	'n=n'	      => \$tools,
	'i=s'         => \$indir,
	'o=s'         => \$outdir,
);
die "$usage\n" if ( !$samp || !$indir || !$outdir);
$bp ||= 50;
$tools ||= 2;

my $key;
my ($chr1,$pos1,$chr2,$pos2,$type);
my (%hash,%chr1,%pos1,%chr2,%pos2,%type);

`less -SN $indir/$samp.lumpy.sv.txt |sort|uniq >$indir/$samp.lumpy.uniq.sv.txt`;
`mv $indir/$samp.lumpy.uniq.sv.txt $indir/$samp.lumpy.sv.txt`;
`less -SN $indir/$samp.manta.sv.txt |sort|uniq >$indir/$samp.manta.uniq.sv.txt`;
`mv $indir/$samp.manta.uniq.sv.txt $indir/$samp.manta.sv.txt`;
`ls $indir/$samp.*.sv.txt >$indir/$samp.sv.list`;

open LST,"$indir/$samp.sv.list" or die $!;
while(<LST>){
        chomp;
        open IN,$_ or die $!;
	while(my $line=<IN>){
        	chomp $line;
	       	$key++;
        	my @line=split "_",$line;
        	($chr1,$pos1,$chr2,$pos2,$type)=@line[0,1,2,3,4];
	        $chr1{$key}=$chr1;
	        $chr2{$key}=$chr2;
	        $pos1{$key}=$pos1;
        	$pos2{$key}=$pos2;
		$type{$key}=$type;
		$hash{$key}=$key;
	}
        close IN;
}
close LST;

open SV,">$outdir/$samp.$tools-tools.overlap.txt" or die $!;
open LST,"$indir/$samp.sv.list" or die $!;
while(<LST>){
        chomp;
        open IN,$_ or die $!;
        while(my $line=<IN>){
                chomp $line;
		my $count=0;	
                my @line=split "_",$line;
 	        ($chr1,$pos1,$chr2,$pos2,$type)=@line[0,1,2,3,4];
	        foreach $key (sort {$a<=>$b} keys %hash){
#			print "$line\t$chr1{$key}_$pos1{$key}_$chr2{$key}_$pos2{$key}_$type{$key}\n";
        		if (($chr1{$key} eq $chr1) and ($chr2{$key} eq $chr2) and ($type{$key} eq $type)){
                		if((abs($pos1-$pos1{$key}) <$bp) and (abs($pos2-$pos2{$key}) <$bp)){
		                        $count++;
                		}
	        	}
        	}

		if($count ge $tools){
			print SV "$line\n";
		}
#		print "$count\t$line\n";
	}
        close IN;
}
close LST;
close SV;

`less -SN $outdir/$samp.$tools-tools.overlap.txt |sort|uniq >$outdir/$samp.$tools-tools.overlap.uniq.txt`;
#`mv $outdir/$samp.$tools-tools.overlap.uniq.txt $outdir/$samp.$tools-tools.overlap.txt`;

my (%uniq,%uchr1,%uchr2,%upos1,%upos2,%utype);
open SV,"$outdir/$samp.$tools-tools.overlap.uniq.txt" or die $!;
open UNIQSV,">$outdir/$samp.$tools-tools.overlap.txt" or die $!;
        while(my $line=<SV>){
                chomp $line;
                my $count=0;
		my $c=0;
		$c++;
                my @line=split "_",$line;
                ($chr1,$pos1,$chr2,$pos2,$type)=@line[0,1,2,3,4];
                $uniq{$c}="1";
		if(not exists $uchr1{$c}){
		print UNIQSV "$line\n";
		}else{
			foreach my $k (sort keys %uniq){
        	                if (($uchr1{$k} eq $chr1) and ($uchr2{$k} eq $chr2) and ($utype{$k} eq $type)){
                	                if((abs($pos1-$upos1{$k}) <$bp) and (abs($pos2-$upos2{$k}) <$bp)){
                        	                $count++;
                                	}
	                        }
        	        }
		if($count==0){print UNIQSV "$line\n"};
		}
		$uchr1{$c}=$chr1;
		$uchr2{$c}=$chr2;
                $upos1{$c}=$pos1;
                $upos2{$c}=$pos2;
		$utype{$c}=$type;
		
}
close SV;
close UNIQSV;
