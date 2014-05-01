create temp table consensus as
select position, case max(sum(A), sum(C), sum(T), sum(G))
        WHEN sum(A) then 'A'
        WHEN sum(C) then 'C'
        WHEN sum(G) then 'G'
        WHEN sum(T) then 'T'
        END nuc, chromosome
from pileupnd
where animal like $animal and day like $day
group by chromosome, position;
