args <- commandArgs(TRUE)
args

databaseFile <- args[1]
outputFileFormat <- args[2]   # sprintf(args[2],samples$animal[sid], samples$day[sid], alias)
								# "%s_d%d_%s.pdf"
titleAddition <- args[3]

# CREATE TABLE chromosomes(name, length);
# CREATE TABLE chromosome_aliases (name, alias);
# CREATE TABLE pileup(animal, day, chromosome, position INT, A INT, C INT, G INT, T INT, D INT);

where <- "1"
if (length(args) > 3) {
    where <- paste("1", args[-c(1,2,3)], sep=" AND ", collapse=" AND ")
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
left join chromosome_aliases on (chromosome = id)
where ", where))

refs <- dbGetQuery(db,'select distinct chromosome, animal, day from pileup;')

dim<-length(refs$animal)

for (i in 1:dim) {
    chr <- dbGetQuery(db,sprintf('select * from div where chromosome = %d and animal = "%s" and day = %d;',
                                 refs$chromosome[i], refs$animal[i], refs$day[i]))
    animal <- chr$animal[1]
    day <- chr$day[1]
    alias <- chr$alias[1]
    length <- max(chr$position)

    filename <- sprintf(outputFileFormat, animal, day, alias)
    print(filename)

    pdf(file=filename,height=8.3, width=(length / 10),onefile = FALSE)
    print(c(min(c(chr$A,chr$C,chr$G,chr$T)),max(c(chr$A,chr$C,chr$G,chr$T))))

    df = data.frame(pos=chr$position, A=chr$A, C=chr$C, G=chr$G, T=chr$T)
    ## convert to long form
    df <- melt(df, id.vars=c("pos"), variable.name="base", value.name="count")
    df.cov <- subset(df, count >= 0)
    df.snp <- subset(df, count < 0)

    c <- ggplot()
    c <- c + geom_bar(data=df.cov, aes(x=pos, y=count, fill=base), stat="identity")
    c <- c + geom_bar(data=df.snp, aes(x=pos, y=count, fill=base), stat="identity")
    c <- c + scale_x_discrete(expand=c(0, 1))
    c <- c + labs(title = sprintf("Diversity in %s day %d, %s segment", animal, day, alias),
                  x = "Position", y = "Coverage", fill = "Base")
    c <- c + theme(axis.text.x = element_text(size=8,angle=90, hjust=1))
    print(c)

    dev.off()
}

dbDisconnect(db)
