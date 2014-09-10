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
} else {
    threshold <- 100
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

strand_bias <- function(Af,Cf,Gf,Tf,Ar,Cr,Gr,Tr) {
    forward <- c(Af,Cf,Gf,Tf)
    reverse <- c(Ar,Cr,Gr,Tr)
    both <- forward + reverse
    pri <- which.max(both)
    both[pri] <- 0
    var <- which.max(both)
    return(1 - fisher.test(matrix(c(forward[pri],reverse[pri],
                                    forward[var],reverse[var]),
                                  nrow=2))$p.value)
}
sbcutoff <- 0.5

library("DBI")
library("RSQLite")
library("ggplot2")

db <- dbConnect(dbDriver("SQLite"), dbname=databaseFile)

animals <- dbGetQuery(db,"select id, name from animals;")
print(animals$name)
if (wholestudy) {
    animals <- data.frame(id=c(1),name=c(1))
}

for (i in 1:length(animals$id)) {

    days <- dbGetQuery(db,sprintf('select distinct day from nucleotides where animal = "%s";',
                                  animals$id[i]))$day
    if (wholestudy) {
        days <- days[1]
    }

    for (j in 1:length(days)) {

        chromosomes <- dbGetQuery(db,'select distinct id, name from chromosomes;')

        if (wholestudy) {
            filename <- sprintf(outputFileFormat)
        } else {
            filename <- sprintf(outputFileFormat, animals$name[i], days[j])
        }
        print(paste(i, ": ", filename))
        pdf(file=filename, width=10, height=8, paper="a4r", onefile=TRUE)

        dat <- list()
        dat.lab <- list()
        ymax <- 0

        for (k in 1:length(chromosomes$id)) {
            if (wholestudy) {
                chr <- dbGetQuery(db,sprintf('select position,
sum(Af) as Af,sum(Cf) as Cf,sum(Gf) as Gf,sum(Tf) as Tf,
sum(Ar) as Ar,sum(Cr) as Cr,sum(Gr) as Gr,sum(Tr) as Tr,
sum(Af)+sum(Cf)+sum(Gf)+sum(Tf)+sum(Ar)+sum(Cr)+sum(Gr)+sum(Tr) as cov
from nucleotides where chromosome = %d group by position;',
                                             chromosomes$id[k]))
            } else {
                chr <- dbGetQuery(db,sprintf('select position,
Af, Cf, Gf, Tf, Ar, Cr, Gr, Tr, Af+Cf+Gf+Tf+Ar+Cr+Gr+Tr as cov from nucleotides
where animal = %d and day = %d and chromosome = %d;',
                                             animals$id[i], days[j], chromosomes$id[k]))
            }
            probA = (chr$Af+chr$Ar)/chr$cov
            probC = (chr$Cf+chr$Cr)/chr$cov
            probG = (chr$Gf+chr$Gr)/chr$cov
            probT = (chr$Tf+chr$Tr)/chr$cov
            entropies <- sapply(probA, entpart) + sapply(probC, entpart)
                       + sapply(probG, entpart) + sapply(probT, entpart)
            sb <- apply(matrix(c(chr$Af,chr$Cf,chr$Gf,chr$Tf,chr$Ar,chr$Cr,chr$Gr,chr$Tr),ncol=8),
                        1, function (x) do.call(strand_bias, as.list(x)))
            sb[sb < sbcutoff] <- sbcutoff

            df <- data.frame(pos=chr$position,ent=entropies,cov=chr$cov,sb=sb)
            df <- df[df$cov > min(threshold,max(df$cov)-1),]
            dat[[k]] <- df
            dat.lab[[k]] <- df[df$ent>0.1,]
            ymax <- max(ymax,max(df$ent))
        }

        for (k in 1:length(chromosomes$id)) {
            p <- qplot(pos,ent,data=dat[[k]],colour=sb)
            p <- p + scale_colour_gradient(limits=c(sbcutoff,1), low="green", high="red")
            if (nrow(dat.lab[[k]]) > 0) {
                p <- p + geom_text(data = dat.lab[[k]], aes(pos,ent, label = pos), hjust = 2,size=3)
            }
            if (wholestudy) {
                p <- p + labs(title = sprintf("Nucleotide entropy %s whole study",
                                  titleAddition),
                              x = chromosomes$name[k], y = "Entropy", colour = "Strand bias")
            } else {
                p <- p + labs(title = sprintf("Nucleotide entropy %s %s Day %d",
                                  titleAddition, animals$name[i], days[j]),
                              x = chromosomes$name[k], y = "Entropy", colour = "Strand bias")
            }
            p <- p + ylim(0,ymax)
            print(p)
        }

        dev.off()
    }
}

dbDisconnect(db)
