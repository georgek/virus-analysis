sbscore <- function(af,ar,df,dr) {
    (abs(df/(af+df)-dr/(ar+dr))/((df+dr)/(af+ar+df+dr)))*(df+dr)
}

args <- commandArgs(TRUE)
args

databaseFile <- args[1]
outputFileFormat <- args[2]
titleAddition <- args[3]

animals <- FALSE
days <- FALSE
if (length(args) > 3) {
    animals <- args[4]
}
if (length(args) > 4) {
    days <- as.integer(args[5])
}

library("DBI")
library("RSQLite")
library("reshape2")
library("ggplot2")

db <- dbConnect(dbDriver("SQLite"), dbname=databaseFile)

if (animals == FALSE) {
    animals <- dbGetQuery(db,"select id,name from animals;")
}

for (i in 1:length(animals$id)) {

    if (days == FALSE) {
        days <- dbGetQuery(db,sprintf('select distinct day from pileup where animal = %d;', animals$id[i]))$day
    }

    for (j in 1:length(days)) {
        chromosomes <- dbGetQuery(db,sprintf('select id,name from chromosomes;', animals$id[i], days[j]))
        aliases <- chromosomes$name
        chromids <- chromosomes$id
        dat <- list()
        lengths <- list()

        dbSendQuery(db, sprintf("create temp table consensus as
select position, case max(sum(A), sum(C), sum(T), sum(G))
	WHEN sum(A) then 'A'
	WHEN sum(C) then 'C'
	WHEN sum(G) then 'G'
	WHEN sum(T) then 'T'
	END nuc, chromosome
from pileupnd
where animal = %d and day = %d
group by chromosome, position;", animals$id[i], days[j]));

        dbSendQuery(db, sprintf("create temp table alternative as
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
where animal = %d and day = %d
group by chromosome, position;", animals$id[i], days[j]));

        ## set up data frames
        for (k in 1:length(chromids)) {
            print(sprintf("Animal: %s, day: %d, segment: %s", animals$name[i], days[j], aliases[k]))
##             chr <- dbGetQuery(db, sprintf('select position,Af,Ar,Cf,Cr,Gf,Gr,Tf,Tr from div
## where animal = %d and day = %d and chromosome = %d;',
##                                           animals$id[i], days[j], chromids[k]))

            chr <- dbGetQuery(db, sprintf("select position, c.nuc as consensus, a.nuc as alternative,
case c.nuc when 'A' then Af
           when 'C' then Cf
           when 'G' then Gf
           when 'T' then Tf
end af,
case c.nuc when 'A' then Ar
           when 'C' then Cr
           when 'G' then Gr
           when 'T' then Tr
end ar,
case a.nuc when 'A' then Af
           when 'C' then Cf
           when 'G' then Gf
           when 'T' then Tf
end df,
case a.nuc when 'A' then Ar
           when 'C' then Cr
           when 'G' then Gr
           when 'T' then Tr
end dr
from pileup
join consensus as c using (chromosome, position)
join alternative as a using (chromosome, position)
where animal = %d and day = %d and chromosome = %d;", animals$id[i], days[j], chromids[k]))

            sbs <- sbscore(chr$af, chr$ar, chr$df, chr$dr)
            dat[[k]] <- data.frame(pos=c(chr$position,chr$positon),
                                   base=c(chr$consensus, chr$alternative),
                                   sb=c(sbs,sbs))
            lengths[[k]] <- max(dat[[k]]$pos) - min(dat[[k]]$pos)
        }

        ## set up PDF
        maxlength <- max(unlist(lengths))
        filename <- sprintf(outputFileFormat, animals$name[i], days[j])
        print(filename)
        pdf(file=filename, height=8.3, width=(maxlength / 10), onefile=TRUE)

        ## print plots
        for (k in 1:length(chromids)) {
            c <- ggplot()
            c <- c + geom_bar(data=dat[[k]], aes(x=pos, y=sb, fill=base), stat="identity")
            ## c <- c + geom_bar(data=dat[[k]], aes(x=pos, y=sb, fill=alt), stat="identity")
            c <- c + scale_x_discrete(expand=c(0, (maxlength-lengths[[k]])/2 + 1))
            c <- c + labs(title = sprintf("Strand bias in %s day %d, %s segment", animals$name[i], days[j], aliases[k]),
                          x = "Position", y = "SB score", fill = "Base")
            c <- c + theme(axis.text.x = element_text(size=8,angle=90, hjust=1))
            suppressWarnings(print(c))  # ggplot warns about negative bars, but it's ok
        }
        dev.off()
        dbSendQuery(db, "drop table consensus;");
        dbSendQuery(db, "drop table alternative;");
    }
    days <- FALSE
}

dbDisconnect(db)
