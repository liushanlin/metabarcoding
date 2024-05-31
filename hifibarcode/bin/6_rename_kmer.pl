#! usr/bin/perl -w
use strict;
use File::Basename qw(dirname basename);

die "perl <prefix> <target length> <clean intermedia files or not 1/0>" unless (@ARGV==3);
my $pre = shift;
my $target = shift;
my $clean = shift;

#clean the intermedia files 

if ($clean == 1 ){
    open my $fh, "$pre\.endList" or die $!;
    while(<$fh>){
        chomp;
        unlink("$_");
    }
    close $fh;
}

open ASS, "$pre\_barcode.contig" or die $!;
open OLD, "$pre\_ends.fasta" or die $!;

open ASF, ">$pre\_barcode.contig\.F" || die $!;
open POI, ">$pre\_barcode.contig.F.add" || die $!;

my %name;
my $count=0;
while(<OLD>){
	chomp;
	next unless (s/^>//);
	$count++;
	$name{$count}=$_;	
}
close OLD;
my %aha;
$/="\>";
while(<ASS>){
    chomp;
    next if ($_ eq "");
    my ($id, $seq) = (split /\n/, $_, 2)[0,1];
    my ($k,$n,$m) = (split /[\_\s]+/, $id)[0,1,2];
    die "not well splited for id $id\n" unless (exists $name{$n});
    $seq=~s/\n//g;
    my $len=length $seq;
    $aha{$n}=0 unless(exists $aha{$n});
    if ($m==1){
            $aha{$n}++;
            print ASF ">$name{$n};k=$k\n$seq\n";
    }elsif($aha{$n}==0){
            $aha{$n}++;
            print ASF ">$name{$n};k=$k\n$seq\n";
    }else{	
            my $lend = abs($len - $target);
			next unless ($lend < 20);
			my $outnum=$m-1;
            print POI ">$name{$n};k=$k;AddSeq$outnum\n$seq\n";
    }
}
close ASS;
close ASF;
close POI;
