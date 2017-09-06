args <- commandArgs(TRUE)
args

databaseFile <- args[1]
outputFileFormat <- args[2]
titleAddition <- args[3]
splitLength <- Inf

if (length(args) > 3) {
    splitLength <- as.integer(args[4])
}

showgc <- FALSE
windowlength <- 10

gccontent <- function(a) {
    (length(a[a=='G']) + length(a[a=='g']) +
     length(a[a=='C']) + length(a[a=='c'])) /
         length(a)
}

## returns windows of length k from a as a matrix with windows in rows
windows <- function(a, k) {
    l <- length(a)
    t(sapply(1:(l-k+1), function(x){a[x:(x+k-1)]}))
}

## returns list of n colours equally spaced in HSL, like ggplot default
gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)[1:n]
    c(hcl(h=hues, l=60, c=110), hcl(h=hues, l=70, c=90))
}

library("DBI")
library("RSQLite")
library("reshape2")
library("ggplot2")

db <- dbConnect(dbDriver("SQLite"), dbname=databaseFile)

dbGetQuery(db, "create temp table consensus as
select position, case max(sum(Af)+sum(Ar), sum(Cf)+sum(Cr), sum(Tf)+sum(Tr), sum(Gf)+sum(Gr))
	WHEN sum(Af)+sum(Ar) then 'A'
	WHEN sum(Cf)+sum(Cr) then 'C'
	WHEN sum(Gf)+sum(Gr) then 'G'
	WHEN sum(Tf)+sum(Tr) then 'T'
	END nuc, chromosome
from nucleotides
group by chromosome, position;
");

dbGetQuery(db, paste("create temp view div as
select animal, day, chromosome, position,
case nuc when 'A' then Af else -Af end Af,
case nuc when 'C' then Cf else -Cf end Cf,
case nuc when 'G' then Gf else -Gf end Gf,
case nuc when 'T' then Tf else -Tf end Tf,
case nuc when 'A' then Ar else -Ar end Ar,
case nuc when 'C' then Cr else -Cr end Cr,
case nuc when 'G' then Gr else -Gr end Gr,
case nuc when 'T' then Tr else -Tr end Tr,
D
from nucleotides
join consensus using (chromosome, position);
"))

animals <- dbGetQuery(db,"select id,name from animals;")

for (i in 1:length(animals$id)) {
    days <- dbGetQuery(db,sprintf('select distinct day from nucleotides where animal = %d;', animals$id[i]))$day
    if (length(days) == 0) {
        print(sprintf("Skipping animal %s.", animals$name[i]))
        next
    }

    for (j in 1:length(days)) {
        chromosomes <- dbGetQuery(db,sprintf('select id,name from chromosomes;',
                                             animals$id[i], days[j]))
        aliases <- chromosomes$name
        chromids <- chromosomes$id
        dat.cov <- list()
        dat.snp <- list()
        dat.del <- list()
        dat.gc <- list()
        lengths <- list()
        coverages <- list()

        ## set up data frames
        for (k in 1:length(chromids)) {
            print(sprintf("Animal: %s, day: %d, segment: %s", animals$name[i], days[j], aliases[k]))
            chr <- dbGetQuery(db, sprintf('select position,Af,Ar,Cf,Cr,Gf,Gr,Tf,Tr,D from div
where animal = %d and day = %d and chromosome = %d;',
                                          animals$id[i], days[j], chromids[k]))
            if (nrow(chr) == 0) {
                lengths[[k]] <- 0
                next
            }

            df <- data.frame(pos=chr$position, Af=chr$Af, Cf=chr$Cf, Gf=chr$Gf, Tf=chr$Tf,
                                              Ar=chr$Ar, Cr=chr$Cr, Gr=chr$Gr, Tr=chr$Tr)
            df.del <- data.frame(pos=chr$position, D=chr$D)
            ## convert to long form
            df <-     melt(df,     id.vars=c("pos"), variable.name="base", value.name="count")
            df.del <- melt(df.del, id.vars=c("pos"), variable.name="base", value.name="count")
            dat.cov[[k]] <- rbind(subset(df, count >= 0),df.del)
            dat.snp[[k]] <- subset(df, count <= 0)
            dat.del[[k]] <- df.del
            lengths[[k]] <- max(df$pos) - min(df$pos)
            coverages[[k]] <- max(chr$Af + chr$Cf + chr$Gf + chr$Tf +
                                  chr$Ar + chr$Cr + chr$Gr + chr$Tr + chr$D)

            consensus <- dbGetQuery(db, sprintf('select position,nuc from consensus
where chromosome = %d;',
                                                chromids[k]))
            conswindows <- windows(consensus$nuc, windowlength)
            gcconts <- apply(conswindows, 1, gccontent)
            dat.gc[[k]] <- data.frame(pos=consensus$position[1:(length(consensus$position)-windowlength+1)],
                                      gc=gcconts)
        }

        ## set up PDF
        maxlength <- max(unlist(lengths))
        maxcoverage <- max(unlist(coverages))
        filename <- sprintf(outputFileFormat, animals$name[i], days[j])
        print(filename)
        pdf(file=filename, height=8, width=(min(splitLength, maxlength) / 10), onefile=TRUE)

        ## print plots
        for (k in 1:length(chromids)) {
            if (lengths[[k]] == 0) {
                next
            }

            if (splitLength > lengths[[k]]) {
                nsegs <- 1
                seglength <- lengths[[k]]
            } else {
                nsegs <- ceiling(lengths[[k]]/splitLength)
                seglength <- floor(lengths[[k]]/nsegs)
            }
            print(sprintf("split length: %.0f, nsegs: %d", splitLength, nsegs))
            start = 1
            end = seglength
            for (l in 1:nsegs) {
                if (l == nsegs) {
                    seglength <- seglength + 1
                }
                print(sprintf("Start: %d, end: %d", start, end))
                c <- ggplot() + scale_fill_manual(values=c(gg_color_hue(4),"#444444"))
                seg.cov <- subset(dat.cov[[k]], pos>=start & pos <=end)
                seg.snp <- subset(dat.snp[[k]], pos>=start & pos <=end)
                c <- c + geom_bar(data=seg.cov, aes(x=pos, y=count, fill=base), stat="identity")
                c <- c + geom_bar(data=seg.snp, aes(x=pos, y=count, fill=base), stat="identity")
                ## c <- c + geom_bar(data=dat.del[[k]], aes(x=pos, y=count), fill="black", stat="identity")
                ## c <- c + scale_x_discrete(expand=c(0, (maxlength-lengths[[k]])/2 + 1))
                ## print(sprintf("splitLength: %d, seglength: %d", splitLength, seglength))
                if (nsegs > 1) {
                    c <- c + scale_x_continuous(expand=c(0, (splitLength - seglength)/2 + 1),
                                                breaks=seq(start, end, by = 1))
                } else {
                    c <- c + scale_x_continuous(expand=c(0, (maxlength - lengths[[k]])/2 + 1),
                                                breaks=seq(start, end, by = 1))
                }
                ## c <- c + expand_limits(y=c(-maxcoverage,maxcoverage))
                if (nsegs > 1) {
                    plottitle <- sprintf("Diversity in %s day %d, %s segment part %d", animals$name[i], days[j], aliases[k], l)
                } else {
                    plottitle <- sprintf("Diversity in %s day %d, %s segment", animals$name[i], days[j], aliases[k])
                }
                c <- c + labs(title = plottitle,
                              x = "Position", y = "Coverage", fill = "Base")
                c <- c + theme(axis.text.x = element_text(size=8,angle=90, hjust=1))
                suppressWarnings(print(c))  # ggplot warns about negative bars, but it's ok
                if (showgc) {
                    gc <- ggplot()
                    seg.gc <- subset(dat.gc[[k]], pos>=start & pos<=end)
                    gc <- gc + geom_line(data=seg.gc, aes(x=pos, y=gc))
                    if (nsegs > 1) {
                        c <- c + scale_x_continuous(expand=c(0, (splitLength - seglength)/2 + 1),
                                                    breaks=seq(start, end, by = 1))
                    } else {
                        c <- c + scale_x_continuous(expand=c(0, (maxlength - lengths[[k]])/2 + 1),
                                                    breaks=seq(start, end, by = 1))
                    }
                    gc <- gc + labs(title = sprintf("GC content in %s segment", aliases[k]),
                                    x = "Position", y = "GC content")
                    gc <- gc + theme(axis.text.x = element_text(size=8,angle=90, hjust=1))
                    print(gc)
                }
                start <- start + seglength
                end <- end + seglength
            }
        }
        dev.off()
    }
    days <- FALSE
}

dbDisconnect(db)
