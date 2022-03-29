#!/usr/bin/perl  -w

=head1 Name

fasta_split.pl

=head1 Description

to split the fasta file into seval

 perl fasta_split.pl  <infile> [fadir] [-Option]
 infile            input fasta file
 1 common options:
 --fadir <str>      directory to store sub fasta file, default ARGV[1] or {basename}_Fa
 --cuts <num>       sequence number per subfile, default 1
 --cutf <num>       specified subfile number to split infile averagely, default notset
 --cutd <num>       maximum file number allow in 2-rank directory, default no 2-rank-dir
 --one <num>        the number of sequence to output one time, default=20
 --selst <file>     only fasta sequence in list can output

 2 option about nib(only used when -nib):
 --nib              change fasta sequence into nib form, then -cuts==1
 --mkgff <file>     masked repeat wiht N, use repeat maksed gff
 --nibfa <file>     set the outfile combine all the outseq in fasta form, default no output
 --bitfa <file>     set the outfile combine all the outseq then change into twobit form
 --nibdir <dir>     output unmasked nib to specified directory when -mgff, defaut no output
 --flst <file>      nib file list pathway, default no output
 --dlst <file>      nib directory pathway list, default no output
 --umflst <file>    unmasked nib file list pathway, default no output
 --umdlst <file>    unmasked nib directory way list, default no ouput
 --help             output help information to screen

=head1 Notice

 1 When -cutf defined, -cuts will be failure.
 2 When -cutf bigger then sequence number, it will equal to sequence number.
 3 When budgetary -cutf bigger than 500 without -cutd set, -cutd will set to be 200

=head1 Example

 perl faspa_split.pl  beer.fa
 perl fasta_split.pl  beer.fa  Beer_Fa/

=cut

use strict;
use lib "";
use COMM qw(abs_path get_path);
use File::Temp;
use Getopt::Long;

use PerlIO::gzip;

my ($fafile,$fadir,$snpfile,$splitn,$fnpdir,$one,$umflst,$bitfa,
$selst,$nib,$dlst,$flst,$umdlst,$nibdir,$nibfa,$mkgff,$help);
GetOptions(
        "fadir:s"=>\$fadir,
    "cuts:i"=>\$snpfile,
    "cutf:i"=>\$splitn,
    "cutd:i"=>\$fnpdir,
    "one:i"=>\$one,
    "nib"=>\$nib,
    "dlst:s"=>\$dlst,
    "flst:s"=>\$flst,
    "umflst:s"=>\$umflst,
    "umdlst:s"=>\$umdlst,
    "nibdir:s"=>\$nibdir,
    "nibfa:s"=>\$nibfa,
    "mkgff:s"=>\$mkgff,
    "selst:s"=>\$selst,
    "bitfa:s"=>\$bitfa,
    "help"=>\$help
);
($help || @ARGV==0) && (die `pod2text $0`);
$fafile = $ARGV[0];
$fadir ||= ($ARGV[1] || 0);
my ($fatonib,$msort,$commlib,$fa_twobit) = get_path("fatonib","msort","commlib","fa_twobit");
if(!$fadir){
        my $bname = (split/\./,(split/\//,$fafile)[-1])[0];
        $fadir = $bname . ($nib ? '_Nib' : '_Fa');
}
$fadir = abs_path($fadir);
$nib && ($snpfile = 1, $splitn = 0, $nibdir &&= abs_path($nibdir), `$commlib`);
fasplit($fafile,$fadir,$snpfile,$splitn,$fnpdir,$one,$nib,$fatonib,$dlst,$selst,
$flst,$umflst,$umdlst,$nibdir,$nibfa,$mkgff,$msort,$fa_twobit,$bitfa,$commlib);#sub1

#sub 1
###########
sub fasplit
###########
{
        my ($fa,$fadir,$snpf,$splitn,$fnpd,$one,$nib,$fatonib,$dlst,$selst,
        $flst,$umflst,$umdlst,$nibdir,$nibfa,$mgff,$msort,$fa_twobit,$bitfa,$commlib) =@_;
  $one ||= 20;
        (-d $fadir) || mkdir"$fadir";
        my $dir = $fadir;
        my $ndir = ($nibdir || "");
        $nibdir && !(-d $nibdir) && mkdir"$nibdir";
        $snpf ||= 1;
        my $seqn = (split/\s+/,($selst ? `wc -l $selst` : `grep -c \">\" $fa`))[0];
        if(!$splitn){
                $splitn = int($seqn/$snpf - 1e-7) + 1;
        }elsif($splitn > $seqn){
                $splitn = $seqn;
        }
        ($splitn > 500) && ($fnpd ||= 200);
  if($fnpd){
        foreach(0..int($splitn / $fnpd - 1e-7)){
                (-d "$fadir/$_") || mkdir"$fadir/$_";
                $nibdir && !(-d "$nibdir/$_") && mkdir"$nibdir/$_";
        }
  }
        my ($j,$n) = (0,0);
        my @outprint;
        my $oneseq_onef = ($splitn>=$seqn) ? 1 : 0;
        chomp(my $temdir = `pwd`);
        my %selfa = ($selst && (-s $selst)) ? (split/\s+/,`awk '{print \$1,1}' $selst`) : ();
        $dlst && (open DLS,">$dlst" || die"$!");
        $flst && (open FLS,">$flst" || die"$!");
        $umflst && (open UFLS,">$umflst" || die"$!");
        $umdlst && (open UDLS,">$umdlst" || die"$!");
        $nibfa && (open NFA,">$nibfa" || die"$!");
        my %repeat = repeat_region($mgff,$msort);#sub1.1
        my $combit;
        $nib && !$nibfa && $bitfa && ($combit = new File::Temp( DIR => $temdir));

if($fa=~/\.gz$/){
        open IN,"<::gzip",$fa;
}else{
        open IN,$fa;
}

        $/=">";<IN>;$/="\n";
        while(<IN>){
                /^(\S+)/ || next;
                chomp(my $id = $_);
                $selst && (-s $selst) && !$selfa{$1} && ($/=">",<IN>,$/="\n",next);
                $/=">";
                chomp(my $seq = <IN>);
                $/="\n";
                my $k = ($j % $splitn);
                $j++;
                if($oneseq_onef){
                        if($fnpd && !($k % $fnpd)){
        my $i=int($k/$fnpd);
        $dir = "$fadir/$i";
        $ndir &&= "$nibdir/$i";
                        }
                        if($nib){
                                if($mgff && $ndir){
                                        my $tmp2 = new File::Temp( DIR => $temdir);
                                        print $tmp2 ">$id\n$seq";
                                        $tmp2->close();
                                        system"$commlib;$fatonib $tmp2 $ndir/$1.nib >& /dev/null";
                                        $tmp2 = undef;
                                        $umdlst && (print UDLS "$1\t$ndir\n");
                                        $umflst && (print UFLS "$1\t$ndir/$1.nib\n");
                                }
                                my $norp = %repeat ? rp_masked(\$seq,$repeat{$1}) : 1;#sub1.2
                                my $tmp = new File::Temp( DIR => $temdir);
                                print $tmp ">$id\n$seq";
                                $tmp->close();
                                ($mgff && $ndir && $norp) ? system"ln -s $ndir/$1.nib $dir/$1.nib" : system"$commlib;$fatonib $tmp $dir/$1.nib >& /dev/null";
                                $tmp = undef;
                                $dlst && (print DLS "$1\t$dir\n");
                                $flst && (print FLS "$1\t$dir/$1.nib\n");
                                $nibfa ? (print NFA ">$id\n$seq") : ($bitfa && (print $combit ">$id\n$seq"));
                        }else{
                                open FAF,">$dir/$1.fa";
                                print FAF ">$id\n$seq";
                                close FAF;
                        }
                }else{
                        $outprint[$k] ? ($outprint[$k] .= ">$id\n$seq") : ($outprint[$k] = ">$id\n$seq");
                        ($j % ($one*$splitn)) ||  outprint($splitn,\$n,$fadir,$dir,$fnpd,\@outprint);#sub1.3
                }
        }
        close IN;
        $nib && !$nibfa && $bitfa && ($combit->close());
        $nib && $bitfa && ($nibfa ? system"$commlib;$fa_twobit $nibfa $bitfa" : system"$commlib;$fa_twobit $combit $bitfa");
        $nib && !$nibfa && $bitfa && ($combit=undef);
        $oneseq_onef || (@outprint && outprint($splitn,\$n,$fadir,$dir,$fnpd,\@outprint));#sub1.3
        $dlst && (close DLS);
        $flst && (close FLS);
        $umflst && (close UFLS);
        $umdlst && (close UDLS);
        $nibfa && (close NFA);
}
#sub1.1
##################
sub repeat_region
##################
{
        my ($gff,$msort) = @_;
        my $del_gff_overlap = 'awk \'(x){if(n==$1){if(e>$4+2){e=(e>$5)?e:$5}else{y=y" "s" "e;n=$1;s=$4;e=$5}}';
        $del_gff_overlap .= 'else{print y;n=$1;s=$4;e=$5;y=n"\t"s-1" "e-1}}(!x && !/^#/)';
        $del_gff_overlap .= '{n=$1;s=$4;e=$5;y=n"\t"s-1" "e-1;x=1}END{print y,s-1,e-1}\'';
        $gff ? split/\t|\n/,`$msort -k 'm1,n4,n5' $gff | $del_gff_overlap` : ();
}
#sub1.2
############
sub rp_masked
############
{
        my ($seq,$region) = @_;
        $region || return(1);
        my @regions = split/\s+/,$region;
        foreach(0..@regions/2-1){
                my @l = @regions[2*$_,2*$_+1];
                my $leng = $l[1] - $l[0] + 1;
                substr($$seq,$l[0],$leng) = 'N' x $leng;
        }
        return(0);
}

#sub1.3
#############
sub outprint
#############
{
        my ($splitn,$noo,$fadir,$dir,$fnpd,$outprint) = @_;
        foreach(0..$splitn-1){
                ${$outprint}[$_] || next;
                if($fnpd && !($_ % $fnpd)){my $i=int($_/$fnpd);$dir = "$fadir/$i";}
                $$noo ? (open FAF,">>$dir/fa$_.fa") : (open FAF,">$dir/fa$_.fa");
                print FAF ${$outprint}[$_];
                close FAF;
                ${$outprint}[$_] = "";
                ($_ == $splitn-1) && ($$noo = 1);
        }
        @{$outprint} = ();
}
