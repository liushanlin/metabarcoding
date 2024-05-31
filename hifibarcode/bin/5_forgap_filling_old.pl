#! usr/bin/perl -w
use strict;
use File::Basename qw(dirname basename);
use FindBin qw($Bin);

die "perl <pre> <cpuNum>  " unless (@ARGV==2);

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

for my $key (keys %seqs){
    TNT:for my $i (0..4){
        next TNT unless (defined $seqs{$key} -> [$i][0] and defined $seqs{$key} -> [$i][1]);
        my $sf = $seqs{$key} -> [$i][0];
        my $sr = $seqs{$key} -> [$i][1];
        print OUT ">$key\_$i\n$sf\n$sr\n";
    }
}

close OUT;


open (MIDLIS,">","$pre\_mid.lis") || die $!;
print MIDLIS ">\nf=$pre\_midMerged.fa\n";
close MIDLIS;

unless ((-e "$pre\_midMerged.fa") or (-z "$pre\_midMerged.fa")){
    `$Bin/cmr -a $pre\_mid_1.fasta  -b $pre\_mid_2.fasta  -o $pre\_midMerged.fa -2 $pre\_midFail.1 -3 $pre\.midFail.2 -l 15 -u 120 -c 0.95 -m 0`;
}

`$Bin/barcode -e $pre\_ends.fasta -r $pre\_mid.lis -o $pre\_barcode -x 350 -n 200 -l 119 -v 8 -k 127 -t $cpu`;

=begin comment

open (BA,">>","$dir/../final.sh") || die $!;
#obtain mid.fa
print BA "$Bin/cmr -a $dir/../step1/$ARGV[1]\_mid_1.fasta  -b $dir/../step1/$ARGV[1]\_mid_2.fasta  -o $dir/mid/mid.fa -2 $dir/mid/$ARGV[1].midfail1 -3 $dir/mid/$ARGV[1].midfail2 -l 30 -u 120 -c 0.95 -m 0\n";
#gap close
for (my $i=1;$i<=4;$i++){
	print BA "$Bin/barcode  -e $dir/ends/$ARGV[1].$i.fa -r $dir/mid/mid.lis  -o $dir/asm/bar$i  -x 470 -n 220 -l 100 -v 8 -k 127 -t $cpu \n" ;
	print BA "perl $Bin/6_rename_kmer.pl $dir/asm/bar$i.contig $dir/ends/$ARGV[1].$i.fa $dir/asm\n";
}

print BA "perl $Bin/7_final.pl $dir/asm  bar $dir/../Result\n" ;
close BA;

sub val{
	my ($sor1,$sor2)=@_;
	my $val;
	if ($sor1==1 && $sor2==2){$val=3;}
	if ($sor1==2 && $sor2==1){$val=4;}
	if ($sor1==2 && $sor2==2){$val=2;}
	if ($sor1==1 && $sor2==1){$val=1;}
	return $val;
}
=end comment
=cut
