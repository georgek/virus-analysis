sbscore <- function(af,ar,df,dr) {
    (abs(df/(af+df)-dr/(ar+dr))/((df+dr)/(af+ar+df+dr)))
}

minmaxmean <- function(x) {
    print(sprintf("Min: %f", min(x)))
    print(sprintf("Max: %f", max(x)))
    print(sprintf("Mean: %f", mean(x)))
}

args <- commandArgs(TRUE)
args

databaseFile <- args[1]
outputFileFormat <- args[2]
titleAddition <- args[3]

animal <- FALSE
day <- FALSE
if (length(args) > 3) {
    animal <- args[4]
}
if (length(args) > 4) {
    day <- as.integer(args[5])
}

library("DBI")
library("RSQLite")
library("reshape2")
library("ggplot2")

db <- dbConnect(dbDriver("SQLite"), dbname=databaseFile)

if (animal == FALSE) {
    animals <- dbGetQuery(db,"select id,name from animals;")
} else {
    animals <- dbGetQuery(db,sprintf("select id, name from animals where id = %s;", animal))
}

distance <- 15

nbias <- 0
nnorm <- 0
mbias <- c()
malts <- c()
matchbias <- c()
matchnorm <- c()
for (dif in 1:distance) {
    matchbias[dif] <- 0
    matchnorm[dif] <- 0
}

for (i in 1:length(animals$id)) {

    if (day == FALSE) {
        days <- dbGetQuery(db,sprintf('select distinct day from pileup where animal = %d;', animals$id[i]))$day
    } else {
        days <- day
    }

    for (j in 1:length(days)) {
        chromosomes <- dbGetQuery(db,sprintf('select id,name from chromosomes;', animals$id[i], days[j]))
        aliases <- chromosomes$name
        chromids <- chromosomes$id
        sbs <- list()
        norms <- list()
        con <- list()
        alt <- list()
        dir <- list()
        lengths <- list()
        maxsbs <- list()

        res <- dbGetQuery(db, sprintf("create temp table consensus as
select position, case max(sum(A), sum(C), sum(T), sum(G))
	WHEN sum(A) then 'A'
	WHEN sum(C) then 'C'
	WHEN sum(G) then 'G'
	WHEN sum(T) then 'T'
	END nuc, chromosome
from pileupnd
where animal = %d and day = %d
group by chromosome, position;", animals$id[i], days[j]));

        res <- dbGetQuery(db, sprintf("create temp table alternative as
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
where animal = %d and day = %d and chromosome = %d
and Af+Ar+Cf+Cr+Gf+Gr+Tf+Tr > 2000
and df > 0 and dr > 0 and af > 0 and ar > 0;",
                                          animals$id[i], days[j], chromids[k]))

            if (nrow(chr) > 0) {
                sbs[[k]] <- sbscore(chr$af, chr$ar, chr$df, chr$dr)*1000*((chr$df+chr$dr)/(chr$af+chr$ar))
                norms[[k]] <- ((chr$df+chr$dr)/(chr$af+chr$ar))/sbscore(chr$af, chr$ar, chr$df, chr$dr)
                mbias[[k]] <- abs(chr$df - chr$dr)
                malts[[k]] <- chr$df + chr$dr
                con[[k]] <- chr$consensus
                alt[[k]] <- chr$alternative
                dir[[k]] <- chr$df > chr$dr
                maxsbs[[k]] <- max(sbs[[k]])
                ## minmaxmean(norms[[k]])
            } else {
                sbs[[k]] <- 0
                norms[[k]] <- 0
                con[[k]] <- 0
                alt[[k]] <- 0
                dir[[k]] <- 0
                maxsbs[[k]] <- 0
            }
        }

        for (k in 1:length(chromids)) {
            positions <- which(sbs[[k]] > 10)
            if (length(positions) > 0) {
                for (p in 1:length(positions)) {
                    pos <- positions[p]
                    b <- mbias[[k]][pos]
                    if (pos > distance && pos < length(con[[k]])-distance) {
                        nbias <- nbias+b
                        for (dif in 1:distance) {
                            if (dir[[k]][pos]) {
                                if (alt[[k]][pos] == con[[k]][pos-dif]) {
                                    matchbias[dif] <- matchbias[dif]+b
                                }
                            } else {
                                if (alt[[k]][pos] == con[[k]][pos+dif]) {
                                    matchbias[dif] <- matchbias[dif]+b
                                }
                            }
                        }
                    }
                }
            }
            positions <- which(norms[[k]] > 0.1)
            if (length(positions) > 0) {
                for (p in 1:length(positions)) {
                    pos <- positions[p]
                    b <- malts[[k]][pos]
                    if (pos > distance && pos < length(con[[k]])-distance) {
                        nnorm <- nnorm+b
                        for (dif in 1:distance) {
                            if (dir[[k]][pos]) {
                                if (alt[[k]][pos] == con[[k]][pos-dif]) {
                                    matchnorm[dif] <- matchnorm[dif]+b
                                }
                            } else {
                                if (alt[[k]][pos] == con[[k]][pos+dif]) {
                                    matchnorm[dif] <- matchnorm[dif]+b
                                }
                            }
                        }
                    }
                }
            }
        }

        ## ## set up PDF
        ## maxlength <- max(unlist(lengths))
        ## maxmaxsbs <- max(unlist(maxsbs))
        ## print(maxmaxsbs)
        ## filename <- sprintf(outputFileFormat, animals$name[i], days[j])
        ## print(filename)
        ## pdf(file=filename, height=8.3, width=(maxlength / 10), onefile=TRUE)

        ## ## print plots
        ## for (k in 1:length(chromids)) {
        ##     c <- ggplot()
        ##     c <- c + geom_bar(data=dat[[k]], aes(x=pos, y=sb, fill=base), stat="identity")
        ##     ## c <- c + geom_bar(data=dat[[k]], aes(x=pos, y=sb, fill=alt), stat="identity")
        ##     c <- c + scale_x_discrete(expand=c(0, (maxlength-lengths[[k]])/2 + 1))
        ##     c <- c + scale_y_continuous(limits=c(0,maxmaxsbs*2))
        ##     c <- c + labs(title = sprintf("Strand bias in %s day %d, %s segment", animals$name[i], days[j], aliases[k]),
        ##                   x = "Position", y = "SB score", fill = "Base")
        ##     c <- c + theme(axis.text.x = element_text(size=8,angle=90, hjust=1))
        ##     suppressWarnings(print(c))  # ggplot warns about negative bars, but it's ok
        ## }

        ## dev.off()
        dbGetQuery(db, "drop table consensus;")
        dbGetQuery(db, "drop table alternative;")
    }
    days <- FALSE
}

print("With bias:")
for (dif in 1:distance) {
    print(sprintf("Dis %d: %f", dif, matchbias[dif]/nbias))
}
print("Without bias:")
for (dif in 1:distance) {
    print(sprintf("Dis %d: %f", dif, matchnorm[dif]/nnorm))
}

dbDisconnect(db)
