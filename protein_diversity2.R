## this is pretty inefficient because it rebuilds the protein, consensus and
## diversity tables for each animal/day

rbinde <- function(f1,f2) {
    if (length(f1)==0) {
        return(f2)
    } else {
        return(rbind(f1,f2))
    }
}

args <- commandArgs(TRUE)
args

databaseFile <- args[1]
outputFileFormat <- args[2]
titleAddition <- args[3]

library("DBI")
library("RSQLite")
library("reshape2")
library("ggplot2")

db <- dbConnect(dbDriver("SQLite"), dbname=databaseFile)

genes <- dbGetQuery(db, "select id, name, product from genes;")

animals <- dbGetQuery(db, "select id, name from animals;")

for (i in 1:length(animals$id)) {
    days <- dbGetQuery(db, sprintf("select distinct day from nucleotides where animal = %d;", animals$id[i]))$day
    if (length(days) == 0) {
        next
    }

    for (j in 1:length(days)) {
        dat.cov <- list()
        dat.var <- list()
        lengths <- list()

        for (k in 1:length(genes$id)) {
            print(sprintf("Gene %s (%s)", genes$name[k], genes$product[k]))
            dat.cov[[k]] <- list()
            dat.var[[k]] <- list()
            lengths[[k]] <- 0
            next_start <- 0
            cds_regions <- dbGetQuery(db, sprintf("select chromosome, start, end, strand from cds_regions
join cds on (cds.id = cds) where gene = %s;", genes$id[k]))
            for (l in 1:length(cds_regions$start)) {
                print(sprintf("cds: %s: %s - %s (%s)", cds_regions$chromosome[l], cds_regions$start[l],
                              cds_regions$end[l], cds_regions$strand[l]))
                dbGetQuery(db, (sprintf("create temp view protein as select animal, day, (position-%s)/3+1 as pos,
Ala,Arg,Asn,Asp,Cys,Gln,Glu,Gly,His,Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val,STOP
from amino_acids where chromosome = %s and position >= %s and position <= %s and (position-%s)%%3 = 0;",
                                        cds_regions$start[l], cds_regions$chromosome[l],
                                        cds_regions$start[l], cds_regions$end[l],
                                        cds_regions$start[l])))
                ## make consensus protein
                dbGetQuery(db, "create temp view consensus as
select pos,case max(sum(Ala),sum(Arg),sum(Asn),sum(Asp),sum(Cys),
                    sum(Gln),sum(Glu),sum(Gly),sum(His),sum(Ile),
                    sum(Leu),sum(Lys),sum(Met),sum(Phe),sum(Pro),
                    sum(Ser),sum(Thr),sum(Trp),sum(Tyr),sum(Val),
                    sum(STOP))
   when sum(Ala) then 'Ala'
   when sum(Arg) then 'Arg'
   when sum(Asn) then 'Asn'
   when sum(Asp) then 'Asp'
   when sum(Cys) then 'Cys'
   when sum(Gln) then 'Gln'
   when sum(Glu) then 'Glu'
   when sum(Gly) then 'Gly'
   when sum(His) then 'His'
   when sum(Ile) then 'Ile'
   when sum(Leu) then 'Leu'
   when sum(Lys) then 'Lys'
   when sum(Met) then 'Met'
   when sum(Phe) then 'Phe'
   when sum(Pro) then 'Pro'
   when sum(Ser) then 'Ser'
   when sum(Thr) then 'Thr'
   when sum(Trp) then 'Trp'
   when sum(Tyr) then 'Tyr'
   when sum(Val) then 'Val'
   when sum(STOP) then 'STOP'
   end acid
from protein
group by pos;")

                dbGetQuery(db, "create temp view div as
select animal,day,pos,
case acid when 'Ala' then Ala else -Ala end Ala,
case acid when 'Arg' then Arg else -Arg end Arg,
case acid when 'Asn' then Asn else -Asn end Asn,
case acid when 'Asp' then Asp else -Asp end Asp,
case acid when 'Cys' then Cys else -Cys end Cys,
case acid when 'Gln' then Gln else -Gln end Gln,
case acid when 'Glu' then Glu else -Glu end Glu,
case acid when 'Gly' then Gly else -Gly end Gly,
case acid when 'His' then His else -His end His,
case acid when 'Ile' then Ile else -Ile end Ile,
case acid when 'Leu' then Leu else -Leu end Leu,
case acid when 'Lys' then Lys else -Lys end Lys,
case acid when 'Met' then Met else -Met end Met,
case acid when 'Phe' then Phe else -Phe end Phe,
case acid when 'Pro' then Pro else -Pro end Pro,
case acid when 'Ser' then Ser else -Ser end Ser,
case acid when 'Thr' then Thr else -Thr end Thr,
case acid when 'Trp' then Trp else -Trp end Trp,
case acid when 'Tyr' then Tyr else -Tyr end Tyr,
case acid when 'Val' then Val else -Val end Val,
case acid when 'STOP' then STOP else -STOP end STOP
from protein
join consensus using (pos);")

                prot <- dbGetQuery(db, sprintf("select pos,Ala,Arg,Asn,Asp,Cys,Gln,Glu,Gly,His,Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val,STOP
from div where animal = %d and day = %d;", animals$id[i], days[j]))
                prot$pos <- prot$pos + next_start
                next_start <- max(prot$pos)+1
                prot <- melt(prot, id.vars=c("pos"), variable.name="acid", value.name="count")
                dat.cov[[k]] <- rbinde(dat.cov[[k]], subset(prot, count > 0))
                dat.var[[k]] <- rbinde(dat.var[[k]], subset(prot, count < 0))
                lengths[[k]] <- lengths[[k]] + (max(prot$pos) - min(prot$pos))

                dbGetQuery(db, "drop view div;")
                dbGetQuery(db, "drop view consensus;")
                dbGetQuery(db, "drop view protein;")
            }
        }

        ## set up PDF
        maxlength <- max(unlist(lengths))
        filename <- sprintf(outputFileFormat, animals$name[i], days[j])
        print(filename)
        pdf(file=filename, height=8, width=(maxlength / 10), onefile=TRUE)

        ## print plots
        for (k in 1:length(genes$id)) {
            c <- ggplot()
            c <- c + geom_bar(data=dat.cov[[k]], aes(x=pos, y=count, fill=acid), stat="identity")
            c <- c + geom_bar(data=dat.var[[k]], aes(x=pos, y=count, fill=acid), stat="identity")
            c <- c + scale_x_discrete(expand=c(0, (maxlength-lengths[[k]])/2 + 1))
            c <- c + labs(title = sprintf("Protein diversity in %s day %s, %s", animals$name[i], days[j], genes$product[k]),
                          x = "Position", y = "Coverage", fill = "Acid")
            c <- c + theme(axis.text.x = element_text(size=8,angle=90, hjust=1))
            suppressWarnings(print(c))  # ggplot warns about negative bars, but it's ok
        }
        dev.off()
    }
}
