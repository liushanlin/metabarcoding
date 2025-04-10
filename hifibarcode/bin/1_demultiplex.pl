#! usr/bin/perl -w
use strict;
=head1 Description

        Usage: perl extract.pl <parameter>

        --fq1           "xx.1.fq"
        --fq2           "xx.2.fq"
        --pri           "primer set, format in "forward \t Seq \n reverse \t Seq \n";
        --ind           "index file in format of "id\t sequence\t forward|reverse";
        --max            ""maximum bp trimmed from 5' end | 1 by default";
        --out           "prefix of the output file"
        --help          "print out this information"

=cut

use Getopt::Long;
use FindBin qw($Bin $Script);
use strict;


my ($Fq1,$Fq2,$Pri,$Out,$Ind,$MaxT,$Help);

my%primer=(
        "V" => "(A|C|G)",
        "D" => "(A|T|G)",
        "B" => "(T|G|C)",
        "H" => "(A|T|C)",
        "W" => "(A|T)",
        "S" => "(C|G)",
        "K" => "(T|G)",
        "M" => "(A|C)",
        "Y" => "(C|T)",
        "R" => "(A|G)",
        "N" => "(A|T|C|G)",
);

GetOptions(
        "fq1:s"=>\$Fq1,
        "fq2:s"=>\$Fq2,
        "out:s"=>\$Out,
        "pri:s"=>\$Pri,
        "ind:s"=>\$Ind,
        "max=i"=>\$MaxT,
        "help"=>\$Help
);

die `pod2text $0` if ($Help || !defined ($Fq1) || !defined ($Fq2)|| !defined ($Pri) || !defined ($Ind));

unless (defined $MaxT){
        $MaxT = 1;
}

if ($Fq1=~/gz$/){
    open (FFQ, "gzip -dc $Fq1 |") || die $!;
}else{
    open (FFQ, $Fq1) || die $!;
}
if ($Fq2=~/gz$/){
    open (RFQ,"gzip -dc $Fq2 |") || die $!;
}else{
    open (RFQ, $Fq2) || die $!;
}

open my $fhEL, ">$Out\.endList" or die $!;

my %indf;
my %indr;
open IND, "<$Ind" || die $!;
my $index_len=&checkIdx;
die "unless not correct index in length\n" unless (defined $index_len);
seek IND,0,0;
while(<IND>){
    chomp;
    s/\r//g;
    my @a=split /\s+/;
    die "index file should be formated as:sampleID index_sequence forward/reverse" unless (@a==3);
    my @mis=&misRev($a[1]);
    if ($a[2] eq "forward"){
        foreach my $ele (@mis){
            if (exists $indf{$ele}){
                if ("$a[0]" ne "$indf{$ele}"){
                    warn "forward barcode $ele of $_ conflict to $indf{$ele}";
                    delete $indf{$ele};
                }
            }else{
                $indf{$ele}=$a[0];
            }
        }
    }elsif($a[2] eq "reverse"){
        foreach my $ele (@mis){
                if (exists $indr{$ele}){
                        if ("$a[0]" ne "$indr{$ele}"){
                            warn "reverse barcode $ele of $_ conflict to $indr{$ele}";
                            delete $indr{$ele};
                        }
                }else{
                        $indr{$ele}=$a[0];
                }
        }
    }else{
        die "the thrid column should be labeled as either forward or reverse\n";
    }
}
close IND;

open PRI, "$Pri" || die $!;
my ($primerF, $primerR, %kh, %khr);
my $k =4; #kmer length used for primer matching
my $pril=0;
while (<PRI>){
    chomp;
    my @a=split /\s+/;
    my $len = length($a[1]);
    $pril = $len >= $pril ? $len : $pril;
    die "primer set is not in correct format\n" unless (@a==2);
    if ($a[0] eq "forward"){
        $primerF = $a[1];
    }elsif($a[0] eq "reverse"){
        $primerR = $a[1];
    }else{
        die "primer set not correct $a[0]\n";
    }
}
close PRI;
die "primer set not defined\n" unless (defined $primerF and defined $primerR);
my $prilf = length ($primerF);
my $prilr = length ($primerR);

my $checkL = $index_len + $pril;
die "unless not correct index in length\n" unless (defined $checkL and $checkL > $index_len and $checkL > $pril and $checkL < 100);

# write out the index + primer length for the following analyses
open my $fhidx, ">$Out\_idxPrimerLength" or die $!;
my $outprif = $index_len + $prilf;
my $outprir =  $index_len + $prilr;
print $fhidx "$outprif\t$outprir\n";
close $fhidx;
 
my $dup =0;
for (my $i = 0; $i <= length($primerF) - $k; $i++) {
    my $kmer = substr($primerF, $i, $k);
    if (exists $kh{$kmer}){
        $dup++;
    }else{
        $kh{$kmer} = $i;
    }
}
my $fullMatchF = length($primerF) - $k - $dup;

$dup =0;
for (my $i = 0; $i <= length($primerR) - $k; $i++) {
    my $kmer = substr($primerR, $i, $k);
    if (exists $khr{$kmer}){
        $dup++;
    }else{
        $khr{$kmer} = $i;
    }
}
my $fullMatchR = length($primerR) - $k - $dup;


## generate output files
my %st;
my %fh;
for my $key (keys %indf){
    my $id=$indf{$key};
    my $fhf=$id."F";
    if (exists $st{$fhf}){
        next;
    }
    $st{$fhf}=0;
    open ($fh{"$fhf"}, ">$Out\_Forward-$id\.fasta") || die $!;
    print $fhEL "$Out\_Forward-$id\.fasta\n";
}

for my $key (keys %indr){
    my $id=$indr{$key};
    my $fhr=$id."R";
    if (exists $st{$fhr}){
        next;
    }
    $st{$fhr}=0;
    open ($fh{"$fhr"}, ">$Out\_Reverse-$id\.fasta") || die $!;
    print $fhEL "$Out\_Reverse-$id\.fasta\n";
}

close $fhEL;

## read, match and split;

open ERR, ">$Out\_err.fasta" || die $!;
open MIF, ">$Out\_mid_1.fasta" || die $!;
open MIR, ">$Out\_mid_2.fasta" || die $!;

my %hash;
my $seqnum =0;
TNT:while(my $t1=<FFQ>){
    $seqnum++;
    chomp($t1);
    chomp(my $seqf=<FFQ>);
    <FFQ>;
    <FFQ>;
    chomp (my $t2=<RFQ>);
    chomp (my $seqr=<RFQ>);
    <RFQ>;
    <RFQ>;
    my $seqLenf = length $seqf;
    my $seqLenr = length $seqr;
    next TNT unless ($seqLenf == $seqLenr and $seqLenr == 150);
    my ($pmf, $dirf, $posf) = &kmerMatch($seqf);
    my ($pmr, $dirr, $posr) = &kmerMatch($seqr);
    if ($pmf == 0 and $pmr == 0){
        print MIF ">$seqnum\_f\n$seqf\n";
        print MIR ">$seqnum\_r\n$seqr\n";
        next TNT;
    }elsif (($pmf == 1 and $pmr == 1) or ($dirf == 3 or $dirr == 3)){
        print ERR ">$seqnum\_f\n$seqf\n>$seqnum\_r\n$seqr\n";
        next TNT;
    }
    unless ($posf == $index_len or $posr == $index_len){
        print ERR ">errIdx_$posf\_$posr\_$seqnum\_f\n$seqf\n>errIdx_$posf\_$posr\_$seqnum\_r\n$seqr\n";
        next TNT;
    }
    my ($trim, $priIdx, $hmd);
    if ($pmf == 1){
        ($trim, $priIdx,$hmd) = &findSeq($seqf, $dirf, $posf);
        if ($hmd > 1){
            print MIF ">ham_$hmd\_$posf\_$seqnum\_f\n$seqf\n";
            print MIR ">ham_$hmd\_$posf\_$seqnum\_r\n$seqr\n";
            next TNT;
        }
        next TNT unless (exists $fh{$priIdx});
        print {$fh{$priIdx}} ">$t1\tfq1\n$trim\n";
    }elsif($pmr == 1){
        ($trim, $priIdx, $hmd) = &findSeq($seqr, $dirr, $posr);
        if ($hmd > 1){
            print MIF ">ham_$hmd\_$posr\_$seqnum\_f\n$seqf\n";
            print MIR ">ham_$hmd\_$posr\_$seqnum\_r\n$seqr\n";
            next TNT;
        }
        next TNT unless (exists $fh{$priIdx});
        print {$fh{$priIdx}} ">$t2\tfq2\n$trim\n";
    }
}	
close FFQ;
close RFQ;

for my $key (sort{$hash{$b} <=> $hash{$a}} keys %hash){
    print "$key\t$hash{$key}\n"
}

sub findSeq{
    my ($seq, $dir, $pos) = @_;
    my $idx = substr $seq, 0, $pos;
    my $dis;
    my $start = $pos + $prilf;
    my $trim = substr ($seq, $start);
    my $id;
    my $prm;
    if ($dir == 1){
        if (exists $indf{$idx}){
            $id = $indf{$idx}."F";
        }else{
            $id = "noFind";
        }
        $dis = &hamming_distance (substr ($seq, $pos, $prilf), $primerF);
    }else{
        if (exists $indr{$idx}){
            $id = $indr{$idx}."R";
        }else{
            $id = "noFind";
        }
        $dis = &hamming_distance (substr ($seq, $pos, $prilr), $primerR);
    }    
    return ($trim, $id, $dis);
}


sub checkIdx {
    my $len=0;
    my %dupCheck;
    while(my $line=<IND>){
        chomp($line);
        my @a = split /\s+/, $line;
        my ($id,$key)=@a[0,1];
        my $leni=length $key;
        if (exists $dupCheck{$id}){
            warn "index ID $id appears more than once\n" unless (@a > 2);
        }else{
            $dupCheck{$id}=1;
        }
        if ($len!=0){
            die "lenth of barcode is various\n" unless ($len==$leni);
        }else{
            $len=$leni;
            }
    }
    return $len;
}

sub misRev {
        my $s=shift @_;
        my @set;
        my @a=split /\s*/,$s;
        for my $i (0..$#a){
                my $f=substr $s,0,$i;
                my $r=substr $s,$i+1,$#a;
                foreach my $b ("A","T","C","G","N"){
                        my $idx=sprintf "$f$b$r";
                        push @set,$idx;
                }
        }
        for my $i (1..$MaxT){
                my $trims = substr $s, $i;
                push @set, $trims;
        }
        return @set;
}


sub kmerMatch{
    my $seq2 = shift @_;
    my %alnf;
    my %alnr;
    my $matchf = 0;
    my $posf = 0;
    my $matchr = 0;
    my $posr = 0;
    for (my $i = 0; $i <= $checkL - $k; $i++) {
        my $kmer = substr($seq2, $i, $k);
        if (exists $kh{$kmer}){
            my $score = $i - $kh{$kmer};
            if (exists $alnf{$score}){
                $alnf{$score}++;
                if ($alnf{$score} > $matchf) {
                    $matchf = $alnf{$score};
                    $posf = $score;
                }
            }else{
                $alnf{$score}=1;
            }
        }
        if (exists $khr{$kmer}){
            my $score = $i - $khr{$kmer};
            if (exists $alnr{$score}){
                $alnr{$score}++;
                if ($alnr{$score} > $matchr) {
                    $matchr = $alnr{$score};
                    $posr = $score;
                }
            }else{
                $alnr{$score}=1;
            }
        }
    }
    my $fr = ($matchf+$k)/$fullMatchF;
    my $rr = ($matchr+$k)/$fullMatchR;
    if($fr >= 0.9 and $rr >= 0.9){
        return (1, 3, 0);
    }elsif ($fr >= 0.90){
        return (1, 1, $posf);
    }elsif($rr >= 0.90){
        return (1, 2, $posr);
    }else{
        return (0, 0, 0);
    }
}

sub hamming_distance {
    my ($str1, $str2) = @_;
    my $distance = 0;
    for my $i (0 .. length($str1) - 1) {
        $distance++ if substr($str1, $i, 1) ne substr($str2, $i, 1);
    }
    return $distance;
}
