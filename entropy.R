## entropy.R db outputname titleaddition [threshold]

## output name is sprintfed, positions with covereage under threshold are
## ignored

args <- commandArgs(TRUE)
args

databaseFile <- args[1]
outputFileFormat <- args[2]
titleAddition <- args[3]

if (length(args) > 3) {
    threshold <- as.integer(args[4])
}

entpart <- function(x) {
    if (is.nan(x) || x == 0) {
        0
    } else {
        -x*log(x)
    }
}

library("DBI")
library("RSQLite")
library("ggplot2")

db <- dbConnect(dbDriver("SQLite"), dbname=databaseFile)

animals <- dbGetQuery(db,"select distinct animal from pileup;")$animal

for (i in 1:length(animals)) {

    days <- dbGetQuery(db,sprintf('select distinct day from pileup where animal = "%s";',animals[i]))$day

    for (j in 1:length(days)) {
        chromosomes <- dbGetQuery(db,sprintf('select distinct chromosome,alias from pileup
left join chromosome_aliases on (chromosome = id)
where animal = "%s" and day = %d;',
                                             animals[i], days[j]))
        aliases <- chromosomes$alias
        chromids <- chromosomes$chromosome

        filename <- sprintf(outputFileFormat, animals[i], days[j])
        print(paste(i, ": ", filename))
        pdf(file=filename, width=10, height=8, paper="a4r", onefile=TRUE)

        dat <- list()
        dat.lab <- list()
        ymax <- 0

        for (k in 1:length(chromids)) {
##             limits <- dbGetQuery(db,sprintf('select min(position) as min, max(position) as max from pileup
## where animal = "%s" and day = %d and chromosome = %d;',
##                                             animals[i], days[j], chromids[k]))
            ## start <- limits$min + trimFront
            ## end <- limits$max - trimEnd
##             chr <- dbGetQuery(db,sprintf('select position,A,C,G,T,A+C+G+T as cov from pileup
## where animal = "%s" and day = %d and chromosome = %d and position >= %d and position <= %d;',
##                                          animals[i], days[j], chromids[k], start, end))
            chr <- dbGetQuery(db,sprintf('select position,A,C,G,T,A+C+G+T as cov from pileup
where animal = "%s" and day = %d and chromosome = %d;',
                                         animals[i], days[j], chromids[k]))
            probA = chr$A/chr$cov
            probC = chr$C/chr$cov
            probG = chr$G/chr$cov
            probT = chr$T/chr$cov
            entropies <- sapply(probA, entpart) + sapply(probC, entpart)
                         + sapply(probG, entpart) + sapply(probT, entpart)

            df <- data.frame(pos=chr$position,ent=entropies,cov=chr$cov)
            df <- df[df$cov > min(threshold,max(df$cov)-1),]
            dat[[k]] <- df
            dat.lab[[k]] <- df[df$ent>0.1,]
            ymax <- max(ymax,max(df$ent))
        }

        for (k in 1:length(chromids)) {
            p <- qplot(pos,ent,data=dat[[k]],colour=cov)
            if (nrow(dat.lab[[k]]) > 0) {
                p <- p + geom_text(data = dat.lab[[k]], aes(pos,ent, label = pos), hjust = 2,size=3)
            }
            p <- p + labs(title = sprintf("Nucleotide entropy %s %s Day %d", titleAddition, animals[i], days[j]),
                          x = aliases[k], y = "Entropy", colour = "Coverage")
            p <- p + ylim(0,ymax)
            print(p)
        }

        dev.off()
    }
}
