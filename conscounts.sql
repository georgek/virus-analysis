-- gives counts of consensus base and 
create temp table conscounts as
select p.position, p.chromosome,
case c.nuc
     when 'A' then Af
     when 'C' then Cf
     when 'G' then Gf
     when 'T' then Tf
end af,
case c.nuc
     when 'A' then Ar
     when 'C' then Cr
     when 'G' then Gr
     when 'T' then Tr
end ar,
case a.nuc
     when 'A' then Af
     when 'C' then Cf
     when 'G' then Gf
     when 'T' then Tf
end df,
case a.nuc
     when 'A' then Ar
     when 'C' then Cr
     when 'G' then Gr
     when 'T' then Tr
end dr
from pileup as p
join consensus as c using (position, chromosome)
join alternative as a using (position, chromosome)
where animal like $animal and day like $day;
