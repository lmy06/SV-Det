#!/usr/bin/perl -w

my $usage = <<'USAGE';
Discription: to filt the structure variations called by multi-tools
Author:leimengyue
contact:leimengyue@genomics.cn
Usage

-i           <in.input>sample_infomation file in vcf format [require]
-o           output directory [require]
-t           the tool that be used to call mutation "crest|delly|lumpy|manta|novoBreak" [require]
-sam         the name of tumor sample    [require]
Example
perl a02multi.SV.filt.pl  -i sample.vcf  -t novoBreak -sam tumor.sample -o out.dir
USAGE
use strict;
use File::Basename;
use Getopt::Long;
use lib 'lib/perl5/site_perl/5.8.8';
my (    $infile, $out, $tool, $sam );
GetOptions(
        'i=s'         => \$infile,
        'o=s'         => \$out,
        't=s'         => \$tool,
        'sam=s'       => \$sam,
);
die "$usage\n" if ( !$infile || !$out || !$tool || !$sam );
if( not -d $out){`mkdir -p $out`;}
open OUT,">$out/$sam.$tool.sv.txt" or die $!;

if (! $infile) {print "you must input the vcf file!\n";}
if($infile =~ /gz$/){open IN,"gzip -dc $infile | " or die $!;}
else{open IN, $infile or die $!;}
if(! -d $out){`mkdir $out`;}
while(<IN>){
	chomp;
	next if ($_=~/#/);
	my @line=split"\t",$_;
	my ($chr1,$chr2,$pos1,$pos2,$svtype);
	if($tool eq "crest"){ 
		($chr1,$chr2,$pos1,$pos2,$svtype)=@line[0,4,1,5,8]; 
	}elsif($tool eq "delly"){
	        next if ($line[6] ne "PASS");
	        my ($i)=$_=~/;CHR2=([\w\.]+);/;
	        my ($j)=$_=~/;END=([\d\.]+);/;
        	my $t=(split"<",$line[4])[1];
	        $svtype=(split">",$t)[0];
		$chr2=$i;
		next if($chr2=~/GL000/);			
	        $pos1=$line[1];
        	$pos2=$j;
		$chr1=$line[0];
                next if($chr1=~/GL000/);
	}elsif($tool eq "lumpy"){
	        my ($i)=$_=~/SVTYPE=([\w\.]+);/;
        	my ($j)=$_=~/END=([\d\.]+);/;
	        next if ($line[-1] ne "./.:0:0:0");
	        my $alt=(split":",$line[-2])[-3];
		next if ($alt le 5);
		$chr1=$line[0];	
		$pos1=$line[1];
		$svtype=$i;		
		if($svtype ne "BND"){
                	$chr2=$chr1;
	                $pos2=$j;
	        }else{
        	        my $info=$line[4];
	               	my @info;
        	        my @a=split":",$info;
	                my $c1=(split"\\[",$a[0])[-1];
        	        my $c2=(split"\\]",$c1)[-1];
                	my $p1=(split"\\[",$a[1])[0];
	                my $p2=(split"\\]",$p1)[0];
			$chr2=$c2;
			$pos2=$p2;
		}
	}elsif($tool eq "manta"){
                my ($i)=$_=~/SVTYPE=([\w\.]+);/;
                my ($j)=$_=~/END=([\d\.]+);/;
                $chr1=$line[0];
                $pos1=$line[1];
                $svtype=$i;
                if($svtype ne "BND"){
                        $chr2=$chr1;
                        $pos2=$j;
                }else{
                        my $info=$line[4];
                        my @info;
                        my @a=split":",$info;
                        my $c1=(split"\\[",$a[0])[-1];
                        my $c2=(split"\\]",$c1)[-1];
                        my $p1=(split"\\[",$a[1])[0];
                        my $p2=(split"\\]",$p1)[0];
                        $chr2=$c2;
                        $pos2=$p2;
		}		
        }elsif($tool eq "novoBreak"){
               	my ($i)=$_=~/;CHR2=([\w\.]+);/;
              	my ($j)=$_=~/;END=([\d\.]+);/;	
	    	my $t=(split"<",$line[4])[1];
        	$svtype=(split">",$t)[0];		
		$chr2=$i;
                next if($chr2=~/GL000/);
                $pos1=$line[1];
                $pos2=$j;
                $chr1=$line[0];
		next if($chr1=~/GL000/);
	}
	if($svtype eq "CTX"||$svtype eq "ITX"||$svtype eq "TRA"){$svtype = "BND"};	
	if($svtype eq "INS"){$svtype = "DUP"};
	&SITES($chr1,$chr2,$pos1,$pos2,$svtype);
}
close IN;

sub SITES{
        my($chr1,$chr2,$pos1,$pos2,$svtype) = @_;
        if ($chr1 eq $chr2){
        if($pos1 > $pos2){my $t=$pos2;$pos2=$pos1;$pos1=$t};
        	print OUT "$chr1\_$pos1\_$chr2\_$pos2\_$svtype\n";	
	}else{
                if($chr1 le $chr2){
                        print OUT "$chr2\_$pos2\_$chr1\_$pos1\_$svtype\n";
	        }else{
                        print OUT "$chr1\_$pos1\_$chr2\_$pos2\_$svtype\n";                }
        }
}
