#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($SN,$SP);
my $prefix;
my ($Base,$MaskFa);
my $Sub_Base;
GetOptions(

        "SN:i"  =>      \$SN,
        "SP:i"  =>      \$SP,

        "p:s"   =>      \$prefix,

        "Base:s"=>      \$Base,
        "MaskFa"=>      \$MaskFa,

        "SB:s"  =>      \$Sub_Base,
);
$SN ||= 100;
$SP ||= 10;
$prefix ||= "Xxxx";

$Base ||= "ATCG";

$Sub_Base ||= "N";

@ARGV || die "Usage:perl $0 <fa>\n";

my $fa = shift;

#my @BASE = ("A","T","C","G");
my @BASE = $Base ? (split //,$Base) : ("A","T","C","G");

foreach my $b(@BASE){
        $b eq $Sub_Base && die "It will keep going on forever\n";
}

#print join("\n",@BASE),"\n";exit;

if($fa=~/gz$/){
        open IN,"<::gzip",$fa;
}else{
        open IN,$fa;
}
open O,">$prefix.simple.mask.list";
$/=">";<IN>;$/="\n";
while(<IN>){
        my $id = (split /\s+/,$_)[0];
        $/=">";
        my $seq = <IN>;
        $seq =~ s/\s+//g;
        $seq =~ s/>//;
        $/="\n";

        $seq = uc $seq;
        my $len = length $seq;

        my $seq_tmp = $seq;
        foreach my $b(@BASE){
#print "$b\n";exit;
                while($seq_tmp =~ /($b+)/g){
                        my $N = $1;
                        my $Nlen = length $N;

                        $Nlen >= $SN || next;

                        my $e = pos $seq_tmp;
                        my $s = $e - $Nlen + 1;

#                       print "$id\t$s\t$e\t$Nlen\t$N\n";

                        my $ms = $s + $SP - 1;
                        my $me = $e - $SP - 1;
                        $ms > $me && next;
                        my $mlen = $me - $ms + 1;

#                       my $mstr = "N" x $mlen;
                        my $mstr = $Sub_Base x $mlen;

#                       https://www.jianshu.com/p/16f82d799e2c
                        substr($seq_tmp,$ms,$mlen) = $mstr;

                        my $Nbase = (split //,$N)[0];
#                       my $Nbase = $1 if $N =~ /([ATCG])+/;
                        print O "$id\t$len\t$s\t$e\t$Nlen\t$ms\t$me\t$Nbase\n";

                }
        }

if($MaskFa){
#       print ">$id\n$seq\n$seq_tmp\n";
        print ">$id\n$seq_tmp\n";
}

}
close IN;
close O;
