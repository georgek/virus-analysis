CREATE TABLE chromosomes(
       id INTEGER PRIMARY KEY AUTOINCREMENT,
       name TEXT UNIQUE,
       length INTEGER);
CREATE TABLE genes(
       id INTEGER PRIMARY KEY AUTOINCREMENT,
       name TEXT UNIQUE,
       product TEXT);
CREATE TABLE cds(
       id INTEGER PRIMARY KEY AUTOINCREMENT,
       chromosome INTEGER,
       gene INTEGER,
       FOREIGN KEY(chromosome) REFERENCES chromosomes(id),
       FOREIGN KEY(gene) REFERENCES genes(id));
CREATE TABLE cds_regions(
       cds INTEGER,
       number INTEGER,
       strand INTEGER,
       start INTEGER,
       end INTEGER,
       FOREIGN KEY(cds) REFERENCES cds(id));
CREATE TABLE animals(
       id INTEGER PRIMARY KEY AUTOINCREMENT,
       name TEXT UNIQUE,
       ndays INTEGER);
CREATE TABLE read_data(
       animal INTEGER,
       day INTEGER,
       nreads INTEGER,
       paired INTEGER,
       FOREIGN KEY(animal) REFERENCES animals(id));
CREATE TABLE nucleotides(
       animal INTEGER,
       day INTEGER,
       chromosome INTEGER,
       position INTEGER,
       Af INTEGER, Cf INTEGER, Gf INTEGER, Tf INTEGER,
       Ar INTEGER, Cr INTEGER, Gr INTEGER, Tr INTEGER,
       D INTEGER,
       FOREIGN KEY(animal) REFERENCES animals(id),
       FOREIGN KEY(chromosome) REFERENCES chromosomes(id));
CREATE TABLE codons(
       animal INTEGER,
       day INTEGER,
       chromosome INTEGER,
       position INTEGER,
       AAA INTEGER, AAC INTEGER, AAG INTEGER, AAT INTEGER, ACA INTEGER,
       ACC INTEGER, ACG INTEGER, ACT INTEGER, AGA INTEGER, AGC INTEGER,
       AGG INTEGER, AGT INTEGER, ATA INTEGER, ATC INTEGER, ATG INTEGER,
       ATT INTEGER, CAA INTEGER, CAC INTEGER, CAG INTEGER, CAT INTEGER,
       CCA INTEGER, CCC INTEGER, CCG INTEGER, CCT INTEGER, CGA INTEGER,
       CGC INTEGER, CGG INTEGER, CGT INTEGER, CTA INTEGER, CTC INTEGER,
       CTG INTEGER, CTT INTEGER, GAA INTEGER, GAC INTEGER, GAG INTEGER,
       GAT INTEGER, GCA INTEGER, GCC INTEGER, GCG INTEGER, GCT INTEGER,
       GGA INTEGER, GGC INTEGER, GGG INTEGER, GGT INTEGER, GTA INTEGER,
       GTC INTEGER, GTG INTEGER, GTT INTEGER, TAA INTEGER, TAC INTEGER,
       TAG INTEGER, TAT INTEGER, TCA INTEGER, TCC INTEGER, TCG INTEGER,
       TCT INTEGER, TGA INTEGER, TGC INTEGER, TGG INTEGER, TGT INTEGER,
       TTA INTEGER, TTC INTEGER, TTG INTEGER, TTT INTEGER,
       FOREIGN KEY(animal) REFERENCES animals(id),
       FOREIGN KEY(chromosome) REFERENCES chromosomes(id));

CREATE INDEX nuc_animal_chromosome_position on nucleotides (animal,chromosome,position);
CREATE INDEX nuc_animal_chromosome on nucleotides (animal,chromosome);
CREATE INDEX nuc_animal on nucleotides (animal);

CREATE INDEX cod_animal_chromosome_position on codons (animal,chromosome,position);
CREATE INDEX cod_animal_chromosome on codons (animal,chromosome);
CREATE INDEX cod_animal on codons (animal);

CREATE VIEW nucleotides_nd AS
       SELECT animal, day, chromosome, position, 
       Af+Ar AS A, Cf+Cr AS C, Gf+Gr AS G, Tf+Tr AS T, D FROM nucleotides;

CREATE VIEW proteins AS
       SELECT animal, day, chromosome, position,
       GCA+GCC+GCG+GCT AS Ala,
       AGA+AGG+CGA+CGC+CGG+CGT AS Arg,
       AAC+AAT AS Asn,
       GAC+GAT AS Asp,
       TGC+TGT AS Cys,
       CAA+CAG AS Gln,
       GAA+GAG AS Glu,
       GGA+GGC+GGG+GGT AS Gly,
       CAC+CAT AS His,
       ATA+ATC+ATT AS Ile,
       CTA+CTC+CTG+CTT+TTA+TTG AS Leu,
       AAA+AAG AS Lys,
       ATG AS Met,
       TTC+TTT AS Phe,
       CCA+CCC+CCG+CCT AS Pro,
       AGC+AGT+TCA+TCC+TCG+TCT AS Ser,
       ACA+ACC+ACG+ACT AS Thr,
       TGG AS Trp,
       TAC+TAT AS Tyr,
       GTA+GTC+GTG+GTT AS Val,
       TAA+TAG+TGA AS STOP
       FROM codons;
