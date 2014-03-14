## entropy.R db outputname titleaddition threshold wholestudy

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
if (length(args) > 4) {
    wholestudy = TRUE
} else {
    wholestudy = FALSE
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
if (wholestudy) {
    animals <- animals[1]
}

for (i in 1:length(animals)) {

    days <- dbGetQuery(db,sprintf('select distinct day from pileup where animal = "%s";',animals[i]))$day
    if (wholestudy) {
        days <- days[1]
    }

    for (j in 1:length(days)) {

        chromosomes <- dbGetQuery(db,sprintf('select distinct chromosome,alias from pileup
left join chromosome_aliases on (chromosome = id)
where animal = "%s" and day = %d;',
                                             animals[i], days[j]))

        aliases <- chromosomes$alias
        chromids <- chromosomes$chromosome

        if (wholestudy) {
            filename <- sprintf(outputFileFormat)
        } else {
            filename <- sprintf(outputFileFormat, animals[i], days[j])
        }
        print(paste(i, ": ", filename))
        pdf(file=filename, width=10, height=8, paper="a4r", onefile=TRUE)

        dat <- list()
        dat.lab <- list()
        ymax <- 0

        for (k in 1:length(chromids)) {
            if (wholestudy) {
                chr <- dbGetQuery(db,sprintf('select position,sum(A) as A,sum(C) as C,sum(G) as G,sum(T) as T,
sum(A)+sum(C)+sum(G)+sum(T) as cov from pileup where chromosome = %d group by position;',
                                             chromids[k]))
            } else {
                chr <- dbGetQuery(db,sprintf('select position,A,C,G,T,A+C+G+T as cov from pileup
where animal = "%s" and day = %d and chromosome = %d;',
                                             animals[i], days[j], chromids[k]))
            }
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
            if (wholestudy) {
                p <- p + labs(title = sprintf("Nucleotide entropy %s whole study", titleAddition, animals[i], days[j]),
                              x = aliases[k], y = "Entropy", colour = "Coverage")
            } else {
                p <- p + labs(title = sprintf("Nucleotide entropy %s %s Day %d", titleAddition, animals[i], days[j]),
                              x = aliases[k], y = "Entropy", colour = "Coverage")
            }
            p <- p + ylim(0,ymax)
            print(p)
        }

        dev.off()
    }
}

dbDisconnect(db)
