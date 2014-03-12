#!/bin/sh

DB=$1
PROTEIN=$2

cat << EOF
<?xml version="1.0"?>
<?xml-stylesheet type="text/xsl" href="mutation_sites_trans.xsl"?>
<animalset>
<alldays>
EOF


sqlite3 -separator "\t" $DB "select distinct day from codonpileup join proteins on (protein_id = proteins.rowid and proteins.name = '"${PROTEIN}"')" | awk '{printf "<day id =\"%d\" />\n", $1}'

cat << EOF
</alldays>
<allsites>
EOF

sqlite3 -separator ',' $DB "select name, codonpileup.position, count(distinct codonpileup.animal) from codonpileup join proteins on (protein_id = proteins.rowid and proteins.name = '"${PROTEIN}"') join seqmutprobs ON (codonpileup.animal = seqmutprobs.animal and seqmutprobs.position in (codonpileup.position, codonpileup.position + 1, codonpileup.position + 2) AND proteins.chromosome = seqmutprobs.chromosome) group by protein_id, codonpileup.position" |\
awk -F , '{printf "<site site_id=\"%s_%d\" protein=\"%s\" codon_start=\"%d\" count=\"%d\" />\n", $1, $2, $1, $2, $3}'

cat << EOF
</allsites>
EOF

sqlite3 $DB << EOF | awk '{
	if (animal != $1) {if (day != "") print "</day>"; if (site_id != "") print "</site>"; if (animal != "") print "</animal>"; print "<animal id=\""$1"\">"; animal = $1; site_id = ""; day = ""}
	if (site_id != $5) {if (day != "") print "</day>"; if (site_id != "") print "</site>"; print "<site site_id=\""$5"\" protein=\""$3"\" codon_start=\""$4"\" own=\""$9"\">"; site_id = $5; day = ""}
	if (day != $2) {if (day != "") print "</day>"; print "<day id=\""$2"\">"; day = $2}
	print "<codon seq=\""$6"\" aminoacid=\""$8"\" count=\""$7"\" />"
}'


.mode tabs

CREATE temp TABLE codon_translation (codon, aminoacid);
insert into codon_translation (codon, aminoacid) values
('TTT','F'),
('TCT','S'),
('TAT','Y'),
('TGT','C'),
('TTC','F'),
('TCC','S'),
('TAC','Y'),
('TGC','C'),
('TTA','L'),
('TCA','S'),
('TAA','.'),
('TGA','.'),
('TTG','L'),
('TCG','S'),
('TAG','.'),
('TGG','W'),
('CTT','L'),
('CCT','P'),
('CAT','H'),
('CGT','R'),
('CTC','L'),
('CCC','P'),
('CAC','H'),
('CGC','R'),
('CTA','L'),
('CCA','P'),
('CAA','Q'),
('CGA','R'),
('CTG','L'),
('CCG','P'),
('CAG','Q'),
('CGG','R'),
('ATT','I'),
('ACT','T'),
('AAT','N'),
('AGT','S'),
('ATC','I'),
('ACC','T'),
('AAC','N'),
('AGC','S'),
('ATA','I'),
('ACA','T'),
('AAA','K'),
('AGA','R'),
('ATG','M'),
('ACG','T'),
('AAG','K'),
('AGG','R'),
('GTT','V'),
('GCT','A'),
('GAT','D'),
('GGT','G'),
('GTC','V'),
('GCC','A'),
('GAC','D'),
('GGC','G'),
('GTA','V'),
('GCA','A'),
('GAA','E'),
('GGA','G'),
('GTG','V'),
('GCG','A'),
('GAG','E'),
('GGG','G');

select c.animal, day, proteins.name protein, c.position, proteins.name || "_" || c.position, codon, count, aminoacid, sum(c.animal = s.animal)
from codonpileup c
join proteins on (protein_id = proteins.rowid)
join seqmutprobs s ON (s.position in (c.position, c.position + 1, c.position + 2) AND proteins.chromosome = s.chromosome)
join codon_translation using (codon)
where proteins.name = '${PROTEIN}'
group by c.animal, day, proteins.name, c.position, codon
order by c.animal, proteins.name, c.position, day, count desc;

EOF

cat << EOF
</day>
</site>
</animal>
</animalset>
EOF

