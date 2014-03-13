args <- commandArgs(TRUE)
args

databaseFile <- args[1]
outputFileFormat <- args[2]
titleAddition <- args[3]

entpart <- function(x) {
    if (x == 0) {
        0
    } else {
        -x*log(x)
    }
}

library("DBI")
library("RSQLite")
library("ggplot2")

db <- dbConnect(dbDriver("SQLite"), dbname=databaseFile)

refs <- dbGetQuery(db,"select distinct chromosome, animal, day from pileup;")
dim <- length(refs$animal)

chr <- dbGetQuery(db,sprintf('select alias,position,A,C,G,T,A+C+G+T as cov from pileup
left join chromosome_aliases on (chromosome = id)
where chromosome = %d and animal = "%s" and day = %d;',
                             refs$chromosome[1], refs$animal[1], refs$day[1]))
alias <- chr$alias[1]

probA = chr$A/chr$cov
probC = chr$C/chr$cov
probG = chr$G/chr$cov
probT = chr$T/chr$cov
entropies <- sapply(probA, entpart) + sapply(probC, entpart) + sapply(probG, entpart) + sapply(probT, entpart)

filename <- sprintf(outputFileFormat, refs$animal[1], refs$day[1], alias)
pdf(file=filename, width=10, height=8, paper="a4r", onefile=FALSE)

df <- data.frame(pos=chr$position,ent=entropies,cov=chr$cov)
p <- qplot(pos,ent,data=df,colour=cov)
p + labs(title = sprintf("Nucleotide entropy %s %s Day %d", titleAddition, refs$animal[1], refs$day[1]),
         x = alias, y = "Entropy", colour = "Coverage")

dev.off()
