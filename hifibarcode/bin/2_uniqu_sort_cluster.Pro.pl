#! usr/bin/perl -w
use strict;
use File::Basename qw(dirname basename);
use FindBin qw($Bin);

die "perl $0 <prefix>" unless (@ARGV==1);

my $pre= shift;


open IN, "<$pre\.endList" || die $!;

open OUT, ">$pre\.2lis" || die $!;
my $cut=0.1;

while(my $nam = <IN>){
	my %hash;
	my %count;
	chomp($nam);
    my $na = (split /[\.\/\_]/, $nam)[-2];
	open FH1, "<$nam" || die $!;
	while (my $t=<FH1>){
		chomp($t);
		chomp(my $seq=<FH1>);
        my ($tt, $file) = (split /\s+/, $t)[0,-1];
        $tt = "$tt\_$file";
		#seq number
		if (exists $hash{$seq}){ 
			$hash{$seq}.=$tt;
			$count{$seq}++;	
		}else{
			$hash{$seq}=$tt;
			$count{$seq}=1;
		}
	}
	close FH1;
	
    ##usearch delet identity>=98% 
	my $tout="temp.fa";
	my $teuc="temp.uc";
	open (IFA,">",$tout) || die $!;
	for my $keyi  (sort {$count{$b} <=> $count{$a}} keys %hash){
		print IFA ">$keyi\n$keyi\n";
	}
    close IFA;
	system "$Bin/vsearch --cluster_smallmem temp.fa --uc temp.uc --id 0.98";
	open (UC,"<",$teuc) || die $!;
	while (my $line = <UC>){
		chomp($line);
		next unless ($line=~/^H/);
		my ($q,$s)=(split (/\t/,$line))[-2,-1];
		die "seqs meet problem\n$q\n$s\n" unless (exists $count{$q} and exists $count{$s});
        $count{$s}+=$count{$q};
        $hash{$s}.=$hash{$q};
	    delete $hash{$q};
        delete $count{$q};
    }
	close UC;
	unlink $tout;
	unlink $teuc;

	my $topseq;
	my @title;
	my $acum=0;
	my $max;
	TTT:for my $key (sort {$count{$b} <=> $count{$a}} keys %hash){
		$acum++;
        last TTT if ($acum>5);
		if ($acum==1){
			$max=$count{$key};
		}
		my @title=split /\>/,$hash{$key};
		my $rate=$count{$key}/$max;
		last TTT if ($rate<$cut);
		last TTT if ($acum>=3 && $count{$key}<10);
		print OUT ">$na\t$count{$key}\t$key\n";
        foreach my $ele (@title){
            next if ($ele eq "");
            print OUT "$ele\n";
        }
	}
}
close IN;
close OUT;
