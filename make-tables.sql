-- this turns the pileup_temp table into pileup, animals and chromosome
-- indexed tables

create table animals(id integer primary key autoincrement, name);

insert into animals(name) select distinct animal from pileup_temp;

create table chromosomes(id integer primary key autoincrement, name, length);

insert into chromosomes(name) select distinct chromosome from pileup_temp;

create table pileup(animal int, day int, chromosome int, position int,
Af int, Cf int, Gf int, Tf int, Ar int, Cr int, Gr int, Tr int, D int);

insert into pileup select a.id, day, c.id, position, Af, Cf, Gf, Tf, Ar, Cr, Gr, Tr, D 
from pileup_temp as p 
join animals as a on (p.animal = a.name) 
join chromosomes as c on (p.chromosome = c.name);

drop table pileup_temp;

create index pileup_animal_chromosome_position on pileup (animal,chromosome,position);
create index pileup_animal_chromosome on pileup (animal,chromosome);
create index pileup_animal on pileup (animal);

