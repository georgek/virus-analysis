#!/usr/bin/perl

# print "NAME\tPOS\tCOV\tREF\tA\tC\tG\tT\tD\n";

# samtools mpileup <bam> | ./pileup2sqlite.pl | awk '{print "animal\t1\t" $0;}' | sqlite3 <db> -separator "\t" -init <(echo -e "CREATE TABLE if not exists pileup(animal, day, chromosome, position, A, C, G, T, D);\n.mode tabs\n.import /dev/stdin pileup\n")

$bam = $ARGV[0];
$animal = $ARGV[1];
$day = $ARGV[2];
$db = $ARGV[3];

use Fcntl;

# $^F = 255;

# print "samtools mpileup " . $bam . " |";

open PILEUP, "samtools mpileup -BQ0 -d10000000 " . $bam . " |";

pipe DBPIPEREAD, DBPIPEWRITE;

$flags = fcntl(DBPIPEREAD, F_GETFL, 0) or die "Can't get flags for the fd: $!\n";
fcntl(DBPIPEREAD, F_SETFD, $flags & (~FD_CLOEXEC)) or die "Can't set flags for the fd: $!\n";

open DBCMD, "| sqlite3 -echo " . $db;

binmode( DBCMD, ":unix" );

print DBCMD "CREATE TABLE if not exists
pileupdir(animal, day INT, chromosome INT, position INT,
          Af INT, Cf INT, Gf INT, Tf INT,
          Ar INT, Cr INT, Gr INT, Tr INT, D INT);\n";
print DBCMD "CREATE temp TABLE if not exists
pileup_temp(chromosome, position INT,
            Af INT, Cf INT, Gf INT, Tf INT,
            Ar INT, Cr INT, Gr INT, Tr INT, D INT);\n";
print DBCMD ".mode tabs\n";
print DBCMD ".import /dev/fd/" . fileno(DBPIPEREAD) . " pileup_temp\n";
close DBPIPEREAD;

my %nuc;

while(<PILEUP>) {
    @list = split /\t/, $_;
    my $pileupnucleotides = $list[4];
    $nuc{"A"} = 0;
    $nuc{"C"} = 0;
    $nuc{"G"} = 0;
    $nuc{"T"} = 0;
    $nuc{"a"} = 0;
    $nuc{"c"} = 0;
    $nuc{"g"} = 0;
    $nuc{"t"} = 0;
    $nuc{"*"} = 0;
    while (length $pileupnucleotides) {
        my @matches = ($pileupnucleotides  =~ m/\.|,|\*|\+[0-9]+[ACGTNacgtn]+|-[0-9]+[ACGTNacgtn]+|[ACGTNacgtn]/g);
        my %ins = ();
        $pileupnucleotides = "";
        foreach $m (@matches) {
            if ($m eq ".") {
                $nuc{uc($list[2])}++;
            } elsif ($m eq ",") {
                $nuc{lc($list[2])}++;
            } elsif ($m =~ m/\+/) {
                #print "insertion here: " . $m . "\n";
                my @number = ($m =~ m/[0-9]+/g);
                my @sequence = ($m =~ m/[ACGTNacgtn]+/g);
                $pileupnucleotides .= substr $sequence[0], int($number[0]);
            } elsif ($m =~ m/\-/) {
                #print "deletion here: " . $m . "\n";
                my @number = ($m =~ m/[0-9]+/g);
                my @sequence = ($m =~ m/[ACGTNacgtn]+/g);
                $pileupnucleotides .= substr $sequence[0], int($number[0]);
            } else {
                $nuc{$m}++;
            }
        }
    }
    print DBPIPEWRITE $list[0] . "\t" . $list[1] . "\t" . $nuc{"A"} . "\t" . $nuc{"C"} . "\t" . $nuc{"G"} . "\t" . $nuc{"T"} . "\t" . $nuc{"a"} . "\t" . $nuc{"c"} . "\t" . $nuc{"g"} . "\t" . $nuc{"t"} . "\t" . $nuc{"*"} . "\n" ;
}

close DBPIPEWRITE;
close PILEUP;

print DBCMD "begin exclusive transaction;\n";
print DBCMD "CREATE TABLE if not exists chromosomes(name, length INT);\n";

print DBCMD 'insert into chromosomes
select distinct chromosome, 0
from pileup_temp
left join chromosomes on (chromosome = name)
where name is null;
';

printf DBCMD 'insert into pileupdir
select "%s", %d, chromosomes.rowid, position, Af, Cf, Gf, Tf, Ar, Cr, Gr, Tr, D
from pileup_temp
left join chromosomes on (chromosome = name);
', $animal, $day;

print DBCMD "update chromosomes set length = (select max(position) from pileup where chromosome = chromosomes.rowid);\n";
print DBCMD "end transaction;\n";

print DBCMD ".q\n";

close DBCMD;

exit(0);

# Local Variables:
# mode: perl
# End:
