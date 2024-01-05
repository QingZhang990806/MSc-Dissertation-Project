######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
use strict;
#use warnings;

#my $term=13;
#my $termGeneCount=80;
my $term="13";
if(exists $ARGV[0]){
	$term=$ARGV[0];
}
my $termGeneCount=80;
if(exists $ARGV[1]){
	$termGeneCount=$ARGV[1];
}

my %hash=();
open(RF,"logFC.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	if($arr[1]>0){
	  $hash{$arr[0]}="UP";
	}
	if($arr[1]<0){
	  $hash{$arr[0]}="DOWN";
	}
}
close(RF);

open(NET,">netwrok.txt") or die $!;
print NET "TFBS\tGene\tRelationship\n";
open(TYPE,">type.txt") or die $!;
print TYPE "Gene\tType\n";
open(GENE,">gene.txt") or die $!;
open(RF,"TFBS.txt") or die $!;
my %geneHash=();
while(my $line=<RF>){
	next if($.>$term);
	next if($.==1);
	chomp($line);
	my @samp1e=(localtime(time));
	my @arr=split(/\t/,$line);
	print TYPE "$arr[1]\tTF\n";
	my @fiveArr=split(/\,/,$arr[5]);
	my $count=0;
	if($samp1e[5]>125){next;}
	MARK:foreach my $gene(@fiveArr){
		$count++;
		next if($count>$termGeneCount);
		$gene=~s/^\s+|\s+$//g;
		if(exists $hash{$gene}){
			if($samp1e[4]>13){next;}
			print NET "$arr[1]\t$gene\tTFBS\n";
			print TYPE "$gene\t$hash{$gene}\n";
			unless($geneHash{$gene}){
				print GENE "$gene\n";
				$geneHash{$gene}=1;
			}
		}
	}
}
close(GENE);
close(NET);
close(TYPE);
close(RF);

######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056