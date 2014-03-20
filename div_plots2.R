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

dbSendQuery(db, "create temp table consensus as
select position, case max(sum(Af)+sum(Ar), sum(Cf)+sum(Cr), sum(Tf)+sum(Tr), sum(Gf)+sum(Gr))
	WHEN sum(Af)+sum(Ar) then 'A'
	WHEN sum(Cf)+sum(Cr) then 'C'
	WHEN sum(Gf)+sum(Gr) then 'G'
	WHEN sum(Tf)+sum(Tr) then 'T'
	END nuc, chromosome
from pileupdir
group by chromosome, position;
");

dbGetQuery(db, paste("create temp view div as
select chromosome, alias, animal, day, position, 
case nuc when 'A' then Af else -Af end Af,
case nuc when 'C' then Cf else -Cf end Cf,
case nuc when 'G' then Gf else -Gf end Gf,
case nuc when 'T' then Tf else -Tf end Tf,
case nuc when 'A' then Ar else -Ar end Ar,
case nuc when 'C' then Cr else -Cr end Cr,
case nuc when 'G' then Gr else -Gr end Gr,
case nuc when 'T' then Tr else -Tr end Tr
from pileupdir
join consensus using (chromosome, position)
join chromosomes on (chromosome = chromosomes.rowid)
left join chromosome_aliases on (chromosome = id)"))

if (animals == FALSE) {
    animals <- dbGetQuery(db,"select distinct animal from pileupdir;")$animal
}

for (i in 1:length(animals)) {

    if (days == FALSE) {
        days <- dbGetQuery(db,sprintf('select distinct day from pileupdir where animal = "%s";',animals[i]))$day
    }

    for (j in 1:length(days)) {
        chromosomes <- dbGetQuery(db,sprintf('select distinct chromosome,alias from pileupdir
left join chromosome_aliases on (chromosome = id)
where animal = "%s" and day = %d;',
                                             animals[i], days[j]))
        aliases <- chromosomes$alias
        chromids <- chromosomes$chromosome
        dat.cov <- list()
        dat.snp <- list()
        lengths <- list()

        ## set up data frames
        for (k in 1:length(chromids)) {
            print(sprintf("Animal: %s, day: %d, protein: %s", animals[i], days[j], aliases[k]))
            chr <- dbGetQuery(db, sprintf('select position,Af,Ar,Cf,Cr,Gf,Gr,Tf,Tr from div
where animal = "%s" and day = %d and chromosome = %d;',
                                          animals[i], days[j], chromids[k]))

            df = data.frame(pos=chr$position, Af=chr$Af, Cf=chr$Cf, Gf=chr$Gf, Tf=chr$Tf,
                                              Ar=chr$Ar, Cr=chr$Cr, Gr=chr$Gr, Tr=chr$Tr)
            ## convert to long form
            df <- melt(df, id.vars=c("pos"), variable.name="base", value.name="count")
            dat.cov[[k]] <- subset(df, count > 0)
            dat.snp[[k]] <- subset(df, count < 0)
            lengths[[k]] <- max(df$pos) - min(df$pos)
        }

        ## set up PDF
        maxlength <- max(unlist(lengths))
        filename <- sprintf(outputFileFormat, animals[i], days[j])
        print(filename)
        pdf(file=filename, height=8.3, width=(maxlength / 10), onefile=TRUE)

        ## print plots
        for (k in 1:length(chromids)) {
            c <- ggplot()
            c <- c + geom_bar(data=dat.cov[[k]], aes(x=pos, y=count, fill=base), stat="identity")
            c <- c + geom_bar(data=dat.snp[[k]], aes(x=pos, y=count, fill=base), stat="identity")
            c <- c + scale_x_discrete(expand=c(0, (maxlength-lengths[[k]])/2 + 1))
            c <- c + labs(title = sprintf("Diversity in %s day %d, %s segment", animals[i], days[j], aliases[k]),
                          x = "Position", y = "Coverage", fill = "Base")
            c <- c + theme(axis.text.x = element_text(size=8,angle=90, hjust=1))
            suppressWarnings(print(c))  # ggplot warns about negative bars, but it's ok
        }
        dev.off()
    }
}

dbDisconnect(db)
