#! usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);
die "perl $0 <prefix> <fq1> <fq2>" unless (@ARGV==3);
my $pre = $ARGV[0];
open IN, "<$pre\.2lis" || die $!;

if ($ARGV[1]=~/gz$/){
    open (FQ1, "gzip -dc $ARGV[1] |") || die $!;
}else{
    open (FQ1, $ARGV[1]) || die $!; 
}

if ($ARGV[2]=~/gz$/){
    open (FQ2,"gzip -dc $ARGV[2] |") || die $!;
}else{
    open (FQ2, $ARGV[2]) || die $!;
}


open (FFQ, ">$pre\_overlap_1.fasta") || die $!;
open (FRQ, ">$pre\_overlap_2.fasta") || die $!;

my %hash;
my %ha;
$/ = "\>";
while (<IN>){
	chomp;
    next if ($_ eq "");
    my ($id, $titles) = split /\n/, $_, 2;
	my $name=(split /\s+/, $id)[0];
	if (exists $ha{$name}){
		$ha{$name}++;
    }else{
        $ha{$name}=1;
    }
    my $num = sprintf "%.2i", $ha{$name};
    my $uid = "$name\-$num";
	my @a=split /\n/,$titles;
	foreach my $ele (@a){
        my $file;
		if ($ele=~s/\_fq2$//){
            $file = "fq2";
        }elsif($ele=~s/\_fq1$//){
            $file = "fq1";
        }
        if (exists $hash{$ele}){
			warn "$ele exists more than once\n"
		}else{
			$hash{$ele}[0]=$uid;
            $hash{$ele}[1]=$file;
		}
	}	
}
close IN;
$/="\n";

while (my $t1=<FQ1>){
	chomp($t1);
    $t1 = (split /\s+/, $t1)[0];
	chomp (my $f=<FQ1>);
	<FQ1>;
	<FQ1>;
	chomp (my $t2=<FQ2>);
	chomp (my $r=<FQ2>);
    $t2 = (split /\s+/, $t2)[0];
	<FQ2>;
	<FQ2>;
    unless ($t2 eq $t1){
        warn "fastq not well ordered\n";
        next;
    }
	next unless (exists $hash{$t1});
    my $file = $hash{$t1}[1];
    if ($file eq "fq1"){
        print FFQ ">$hash{$t1}[0]\n$f\n";
        print FRQ ">$hash{$t1}[0]\n$r\n"
    }elsif($file eq "fq2"){
        print FFQ ">$hash{$t1}[0]\n$r\n";
        print FRQ ">$hash{$t1}[0]\n$f\n"
    }
}
close FQ1;
close FQ2;
close FFQ;
close FRQ;
`$Bin/cmr -a $pre\_overlap_1.fasta -b $pre\_overlap_2.fasta -o $pre\_endConnected.fasta -2 $pre\_endFail.1 -3 $pre\_endFail.2 -l 15 -u 120 -c 0.95 -m 0`; 
print "done\n"
