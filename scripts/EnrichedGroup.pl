use strict;
use Math::GSL::CDF qw/:all/;
use Math::GSL::Randist qw/:all/;
use Statistics::Basic qw(:all);
use Statistics::Multtest qw(bonferroni holm hommel hochberg BH BY qvalue);
use Statistics::Multtest qw(:all);

my $stFolder=$ARGV[0];					###	SupplementaryTable12.Data/
$stFolder =~ s/\/$//;
my @stData=glob("$stFolder\/*.$ARGV[1].all.txt");	###	1Mb
my $nCut=$ARGV[2];					###	0.05

my %stCheck;
my %nTotal;
my %nList;
my %nGroup;
my $nTotal=0;
my @stRegion=();
for(my $i=0; $i<@stData; $i++)
{
	open(DATA, "$stData[$i]");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo=split(/\t/,$stLine);
		my $stWin="$stInfo[0]	$stInfo[1]	$stInfo[2]";	###	region
		if($stCheck{$stWin} eq "")
		{
			push(@stRegion, "$stWin");
			$stCheck{$stWin}="OK";
		}
		if($stInfo[3] != 0)
		{
			$nTotal+=$stInfo[3];					###	Total number of genes
			$nTotal{$stWin}+=$stInfo[3];				###	Total number of genes in each region
			for(my $j=4; $j<@stInfo; $j++)
			{
				my @nG=split(/=/,$stInfo[$j]);
				$nList{$stWin}{$nG[0]}+=$nG[1];		###	group number in each region
				$nGroup{$nG[0]}+=$nG[1];		###	each group total
			}
		}
	}
	close DATA;
}

open FH, ">", "NLRsPerWindow.Total.$ARGV[1].all.enriched.txt";
open FH2, ">", "NLRsPerWindow.Total.$ARGV[1].enriched.txt";
my @stFinal=();
my @nPval=();

foreach my $stKey (@stRegion)
{
	foreach my $nKey (sort{$nList{$stKey}{$b} <=> $nList{$stKey}{$a}}keys %{$nList{$stKey}})
	{
		my $nA=$nList{$stKey}{$nKey};
		my $nB=$nGroup{$nKey};
		my $nC=$nTotal-$nB;
		my $nD=$nTotal{$stKey};
		my $nPval=gsl_cdf_hypergeometric_Q($nA,$nB,$nC,$nD);
		push(@stFinal, "$stKey	$nD	$nKey=$nA	$nA,$nB,$nC,$nD");
		push(@nPval, "$nPval");
	}
}
my $nRefP=\@nPval;
my $nFDR=BH($nRefP);
my @nValue=values(%nTotal);
my $nRefV=\@nValue;
my $nMedian=median($nRefV);
my @stTemp=();
my $stID="";
for(my $i=0; $i<@stFinal; $i++)
{
	print FH "$stFinal[$i]	$nPval[$i]	$$nFDR[$i]\n";
	if($stID ne "$stFinal[$i]" && $stID ne "")
	{
		@stTemp=sort{($a =~ /	([^\s]+)$/)[0] <=> ($b =~ /	([^\s]+)$/)[0] || ($b =~ /=([0-9]+)/)[0] <=> ($a =~ /=([0-9]+)/)[0]}@stTemp;
		if($stTemp[0] =~ /=([0-9]+)/)
		{
			my $nNum=$1;
			if($nNum > $nMedian)
			{
				my @stInfo=split(/\t/,$stTemp[0]);
				$stInfo[4] =~ s/=.+//;
				if($stInfo[$#stInfo] < $nCut)
				{
					print FH2 "$stInfo[0]	$stInfo[1]	$stInfo[2]	$stInfo[4]\n";
				}
			}
		}
		@stTemp=();
	}
	push(@stTemp, "$stFinal[$i]	$nPval[$i]	$$nFDR[$i]");
	$stID=$stFinal[$i];
}
if($#stTemp != -1)
{
	@stTemp=sort{($a =~ /	([^\s]+)$/)[0] <=> ($b =~ /	([^\s]+)$/)[0] || ($b =~ /=([0-9]+)/)[0] <=> ($a =~ /=([0-9]+)/)[0]}@stTemp;
	if($stTemp[0] =~ /=([0-9]+)/)
	{
		my $nNum=$1;
		if($nNum > $nMedian)
		{
			my @stInfo=split(/\t/,$stTemp[0]);
			$stInfo[4] =~ s/=.+//;
			if($stInfo[$#stInfo] < $nCut)
			{
				print FH2 "$stInfo[0]	$stInfo[1]	$stInfo[2]	$stInfo[4]\n";
			}
		}
	}
}
close FH;
close FH2;
