#! usr/bin/perl -w
use strict;
use File::Basename qw(dirname basename);
use FindBin qw($Bin);

die "perl <pre> <cpuNum>" unless (@ARGV==2);

my $pre=shift;
my $cpu=shift;

my %seqs;
my %deps;
open FA, "$pre\_cluster.fasta" || die $!;
$/="\>";
while (<FA>){
	chomp;
    next if ($_ eq "");
    my ($id, $seq) = split /\n/, $_, 2;
    $seq=~s/\n//g;
    my ($ori, $idx, $rank, $dep, $len) = (split /[\s\-]+/, $id)[0,1,2,3,4];
    die "not well formated of $id by direction\n" unless ($ori eq "Reverse" or $ori eq "Forward");
    my $seqL = length $seq;
    die "not well formated of $id by seq length of $seqL\n" unless ($len == $seqL);
    next unless ($len >= 130);
    next unless ($rank <= 5);
    my $i = $rank -1;
    my $j = 0;
    if ($ori eq "Reverse"){
        $j = 1;
        $seq = reverse $seq;
        $seq =~ tr/ATCGN/TAGCN/;
    }
    $seqs{$idx} -> [$i][$j] = $seq;
    $deps{$idx} -> [$i][$j] = $dep;
}
close FA;
$/="\n";

open OUT, ">$pre\_ends.fasta" or die $!;
open DEP, ">$pre\_ends.depth" or die $!;

for my $key (keys %seqs){
    my (@dfs, @drs);
    TNT:for my $i (0..4){
        next TNT unless (defined $seqs{$key} -> [$i][0] and defined $seqs{$key} -> [$i][1]);
        next TNT unless (defined $deps{$key} -> [$i][0] and defined $deps{$key} -> [$i][1]);
        my $sf = $seqs{$key} -> [$i][0];
        my $sr = $seqs{$key} -> [$i][1];
        my $df = $deps{$key} -> [$i][0];
        my $dr = $deps{$key} -> [$i][1];
        if ($df/$dr >=10 or $dr/$df >=10){
            print DEP "$key\t$i\-$df\t$i\-$dr\tfailed\n";
            next TNT;
        }
        if ($i <= 1){
            print OUT ">$key\_$i\-$i\n$sf\n$sr\n";
            print DEP "$key\t$i\-$df\t$i\-$dr\tanalyzed\n";
        }
        push @dfs, $deps{$key} -> [$i][0];
        push @drs, $deps{$key} -> [$i][1];
    }
    my @matches = &match (\@dfs, \@drs);
    foreach my $ele (@matches){
        next if ($ele eq "0-0" or $ele eq "1-1");
        my ($fi, $ri) = (split /\-/, $ele)[0,1];
        next unless (defined $seqs{$key} -> [$fi][0] and defined $seqs{$key} -> [$ri][1]);
        my $sfa = $seqs{$key} -> [$fi][0];
        my $sra = $seqs{$key} -> [$ri][1];
        my $dfa = $deps{$key} -> [$fi][0];
        my $dra = $deps{$key} -> [$ri][1];
        print OUT ">$key\_$fi\-$ri\-possible\n$sfa\n$sra\n";
        print DEP "$key\t$fi\-$dfa\t$ri\-$dra\tanalyzed\n";
    }
}

close OUT;
close DEP;

open (MIDLIS,">","$pre\_mid.lis") || die $!;
print MIDLIS ">\nf=$pre\_midMerged.fa\n";
close MIDLIS;

unless ((-e "$pre\_midMerged.fa") or (-z "$pre\_midMerged.fa")){
    `$Bin/cmr -a $pre\_mid_1.fasta  -b $pre\_mid_2.fasta  -o $pre\_midMerged.fa -2 $pre\_midFail.1 -3 $pre\.midFail.2 -l 15 -u 120 -c 0.95 -m 0`;
}

`$Bin/barcode -e $pre\_ends.fasta -r $pre\_mid.lis -o $pre\_barcode -x 350 -n 200 -l 119 -v 8 -k 127 -t $cpu`;


sub match {
    my ($f,$r) = @_;
    my @mats;
    FIR:for my $i (0..2){
        next FIR unless (defined $f -> [$i]);
        SEC:for my $j (0..2){
            next SEC unless (defined $r -> [$j]);
            my ($x, $y) = $f -> [$i] > $r -> [$j] ? ( $f -> [$i] , $r -> [$j]) : ( $r -> [$j] , $f -> [$i]);
            if (($y*2) >= $x){
                push @mats, "$i\-$j";
            }
        }
    }
    return @mats;
}
