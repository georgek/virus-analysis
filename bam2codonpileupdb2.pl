#!/usr/bin/perl

# for i in alignments/F*.bam; do echo cd `pwd` "&& source samtools-0.1.18 && ./bam2codonpileupdb.pl $i `echo $i | cut -d / -f 2 | cut -d _ -f 1` `echo $i | cut -d / -f 2 | cut -d _ -f 2 | cut -d d -f 2` ferretCCC.db"; done

$BAM = $ARGV[0];
$animal = $ARGV[1];
$day = $ARGV[2];
$db = $ARGV[3];

use Switch;
use Fcntl;

open SAM, "samtools view $BAM |";

my $maxpos = 0;

pipe DBPIPEREAD1, DBPIPEWRITE1;
# pipe DBPIPEREAD2, DBPIPEWRITE2;

$flags = fcntl(DBPIPEREAD1, F_GETFL, 0) or die "Can't get flags for the fd: $!\n";
fcntl(DBPIPEREAD1, F_SETFD, $flags & (~FD_CLOEXEC)) or die "Can't set flags for the fd: $!\n";
# $flags = fcntl(DBPIPEREAD2, F_GETFL, 0) or die "Can't get flags for the fd: $!\n";
# fcntl(DBPIPEREAD2, F_SETFD, $flags & (~FD_CLOEXEC)) or die "Can't set flags for the fd: $!\n";

# $ENV{HOME}

open DBCMD, "| sqlite3 -echo " . $db;
binmode( DBCMD, ":unix" );

print DBCMD "CREATE temp TABLE temp_codonpileup ( chromosome TEXT , position INT, codon TEXT, count INT );\n";
print DBCMD ".mode csv\n";
print DBCMD ".import /dev/fd/" . fileno(DBPIPEREAD1) . " temp_codonpileup\n";
close DBPIPEREAD1;

my @codonpileup;
my $chromosome = undef;
my $numreads = 0;

my $pi=0;

# open DEBUG, ">> /dev/ttys003";

while(<SAM>) {
	# print DEBUG $_;
	my $samline = $_;
	@sam = split /\t/, $_;

	# $start = $sam[3];
	
	if ($sam[5] eq "*") { next; }
	
	my $name = $sam[0];
	my $read = $sam[9];
	my $pos = int($sam[3]);
	
	my $seqoffset = $codon_start - $sam[3];
	my @cigar = ($sam[5]  =~ m/[0-9]+[MIDNSHPX=]/g);
	
	my $codon = "";
	
	my $debugmsg = "";
	
	my $part = 1;
	
	
	$numreads++;
	if ($chromosome ne undef) {
		for(; $pi < scalar(@codonpileup) && (($pi < $pos) || ($chromosome ne $sam[2])); $pi++) {
			foreach my $codon (keys %{$codonpileup[$pi]}) {
				printf DBPIPEWRITE1 "%s,%d,%s,%d\n", $chromosome, $pi, $codon, $codonpileup[$pi]{$codon};
			}
			$codonpileup[$pi] = undef;
		}
	} 
	if ($chromosome ne $sam[2]) {
		$chromosome = $sam[2];
		$pi = 0;
	}
	
	
	FOR: foreach my $cstr (@cigar) {
		
		# print $cstr;
		
		my $op = substr $cstr, length($cstr) - 1;
		my $olen = int(substr $cstr, 0, length($cstr) - 1);
		
		$debugmsg .= $op . "\t" . $olen . "\n";
		
		switch ($op) {
			case "M"	{
				# printf DBPIPEWRITE1 "%d,%d,%s,%s\n", $pos, $olen, $sam[2], substr($read, 0, $olen);
				
				for ($i = 0; $i < $olen - 2; $i++) {
					$codonpileup[$pos + $i]{substr($read, $i, 3)}++;
				}
				
				$read = substr $read, $olen;
				$pos += $olen;
				$part++;
				
			}
			case "I"	{$read = substr $read, $olen;}
			case "S"	{$read = substr $read, $olen;}
			case "D"	{$pos += $olen;}
			# case "H"	{print STDERR "Unknown CIGAR option " . $op . ", skipped read " . $sam[0] . "\n"; last;}
			# case "P"	{print STDERR "Unknown CIGAR option " . $op . ", skipped read " . $sam[0] . "\n"; last;}
			# case "X"	{print STDERR "Unknown CIGAR option " . $op . ", skipped read " . $sam[0] . "\n"; last;}
			# case "="	{print STDERR "Unknown CIGAR option " . $op . ", skipped read " . $sam[0] . "\n"; last;}
			# case "N"	{print STDERR "Unknown CIGAR option " . $op . ", skipped read " . $sam[0] . "\n"; last;}
			else	{print STDERR "Unknown CIGAR option " . $op . " in CIGAR " . $sam[5] . ", skipped read " . $sam[0] . "\n"; last;}
		}
	}

}

if ($chromosome ne undef) {
	for(; $pi < scalar(@codonpileup); $pi++) {
		foreach my $codon (keys %{$codonpileup[$pi]}) {
			printf DBPIPEWRITE1 "%s,%d,%s,%d\n", $chromosome, $pi, $codon, $codonpileup[$pi]{$codon};
		}
		$codonpileup[$pi] = undef;
	}
}

close DBPIPEWRITE1;

printf DBCMD '
CREATE TABLE IF NOT EXISTS codonpileup (animal, day INT, position INT, codon TEXT, count INT, protein_id int);

INSERT INTO codonpileup
SELECT "%s" , %d , temp_codonpileup.position, codon, count , proteins.rowid
FROM temp_codonpileup
join chromosomes on (temp_codonpileup.chromosome = chromosomes.name)
join proteins on (proteins.chromosome = chromosomes.rowid
and temp_codonpileup.position >= proteins.start
and temp_codonpileup.position <= proteins.end
and (temp_codonpileup.position - proteins.start) % 3 = 0
);
.q
', $animal, $day;

close DBCMD;

exit(0);
