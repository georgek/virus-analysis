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
select position, case max(sum(A), sum(C), sum(T), sum(G))
	WHEN sum(A) then 'A'
	WHEN sum(C) then 'C'
	WHEN sum(G) then 'G'
	WHEN sum(T) then 'T'
	END nuc, chromosome
from pileup
group by chromosome, position;
");

dbGetQuery(db, paste("create temp view div as
select chromosome, alias, animal, day, position, 
case nuc when 'A' then A else -A end A,
case nuc when 'C' then C else -C end C,
case nuc when 'G' then G else -G end G,
case nuc when 'T' then T else -T end T
from pileup
join consensus using (chromosome, position)
join chromosomes on (chromosome = chromosomes.rowid)
left join chromosome_aliases on (chromosome = id)"))

if (animals == FALSE) {
    animals <- dbGetQuery(db,"select distinct animal from pileup;")$animal
}

for (i in 1:length(animals)) {

    if (days == FALSE) {
        days <- dbGetQuery(db,sprintf('select distinct day from pileup where animal = "%s";',animals[i]))$day
    }

    for (j in 1:length(days)) {
        chromosomes <- dbGetQuery(db,sprintf('select distinct chromosome,alias from pileup
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
            chr <- dbGetQuery(db, sprintf('select position,A,C,G,T from div
where animal = "%s" and day = %d and chromosome = %d;',
                                          animals[i], days[j], chromids[k]))

            df = data.frame(pos=chr$position, A=chr$A, C=chr$C, G=chr$G, T=chr$T)
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
