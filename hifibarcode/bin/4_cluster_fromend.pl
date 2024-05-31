#! usr/bin/perl -w
use strict;

die "perl <prefix> <minimum cov for last position|5>" unless (@ARGV==2);

my $pre = shift;
my $posCut = shift;

open FA, "<$pre\_endConnected.fasta" || die $!;
open OU, ">$pre\_cluster.fasta" or die $!;

open IP, "<$pre\_idxPrimerLength" || die $!;
chomp(my $ipline = <IP>);
my ($pfl, $prl) = (split /\s+/, $ipline)[0,1];
die "primer length not well formated\n" unless (defined $pfl and $pfl> 0 and defined $prl and $prl > 0);
close IP;

my %seqs;
$/="\>";
while(<FA>){
	chomp;
    next if ($_ eq "");
    my ($id, $seq) = (split /\n/, $_, 2)[0,1];
	my $idx = (split /[\s\_]+/, $id)[1];
    $seq=~s/\n//g;
    if ($idx=~/Forward/){
        $seq = substr $seq, $pfl;
    }elsif($idx=~/Reverse/){
        $seq = substr $seq, $prl;
    }else{
        die "did not find correct ids for $idx";
    }
    push @{$seqs{$idx}}, $seq;
}
close FA;

for my $key (keys %seqs){
    my ($cov, $seq) = &cluster($key);
    my $len = length ($seq);
    print OU ">$key\t$cov\t$len\n$seq\n";
}
close OU;

#####sub####


sub cluster{
	my $file = $_[0];
	my %hash;
	my $maxlen=0;
	foreach my $rr (@{$seqs{$file}}){
		chomp($rr);
		next if ($rr eq "");
		my @s=split /\s*/, $rr;
		$maxlen = ($#s+1) if ($maxlen < ($#s+1));
		for my $i (0..$#s){
			if (exists $hash{$i} and exists $hash{$i}{$s[$i]}){
				$hash{$i}{$s[$i]}++;     #numbers of ATCG locate in $i
			}else{
				$hash{$i}{$s[$i]}=1;
			}
		}
	}
	my $eve;
	my $fin=0;
	my $seq;
	TTT:for (my $n = $maxlen-1; $n >0; $n--){
		my ($cov,$max) = &total (\%{$hash{$n}});
		if ($max >= $posCut){   #ends trimmed of < 5 read coverag
			$fin=$n+1;
			last TTT;
		}
	}

	for my $n (0..($fin-1)){
		my ($cov,$max,$base) = &total (\%{$hash{$n}});
		$eve+=$max;
        $seq.=$base;
	}
	my $ecov;
	if ($fin!=0){
		$ecov=int($eve/$fin);
	}else{
		$ecov=0;
		$seq="O";
	}
	return ($ecov,$seq);
}

sub total {
	my $su = $_[0];
	my $tn=0;
	my $mn=0;
    my $base;
	for my $sk (keys %{$su}){
        $tn+= $su -> {$sk};
		if ($mn < $su -> {$sk}){
            $mn = $su -> {$sk};
            $base = $sk;
        }
	}
	return ($tn,$mn, $base)
}
