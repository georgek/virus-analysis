create temp table alternative as
select position,
case nuc
    when 'A' then case max(sum(C),sum(G),sum(T))
                      when sum(C) then 'C'
                      when sum(G) then 'G'
                      when sum(T) then 'T'
                  end
    when 'C' then case max(sum(A),sum(G),sum(T))
                      when sum(A) then 'A'
                      when sum(G) then 'G'
                      when sum(T) then 'T'
                  end
    when 'G' then case max(sum(C),sum(A),sum(T))
                      when sum(C) then 'C'
                      when sum(A) then 'A'
                      when sum(T) then 'T'
                  end
    when 'T' then case max(sum(C),sum(G),sum(A))
                      when sum(C) then 'C'
                      when sum(G) then 'G'
                      when sum(A) then 'A'
                  end
end nuc, chromosome
from pileupnd
join consensus using (chromosome, position)
where animal like $animal and day like $day
group by chromosome, position;
