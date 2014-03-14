args <- commandArgs(TRUE)
args

database_filename <- args[1]
output_file_format <- args[2]   # sprintf(args[2],samples$animal[sid], samples$day[sid], alias)
								# "%s_d%d_%s.pdf"
title_addition <- args[3]

nuc_colours <- c("springgreen", "plum", "slateblue", "seagreen")

# CREATE TABLE chromosomes(name, length);
# CREATE TABLE chromosome_aliases (name, alias);
# CREATE TABLE pileup(animal, day, chromosome, position INT, A INT, C INT, G INT, T INT, D INT);

where <- "1"
if (length(args) > 3) {
    where <- paste("1", args[-c(1,2,3)], sep=" AND ", collapse=" AND ")
}


library("DBI")
library("RSQLite")

db <- dbConnect(dbDriver("SQLite"), dbname=args[1])

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

# data

# data$COV = data$A + data$C + data$G + data$T

# max(c(data$A,data$C,data$G,data$T))

refs <- dbGetQuery(db,'select distinct chromosome, animal, day from pileup;')

# refs<-unique(data.frame(name=data$chromosome, alias=data$alias))
refs
dim<-length(refs$animal)

for (i in 1:1) {
    ## print(name)

    ## chr<-subset(sdata,chromosome==name)
    chr <- dbGetQuery(db,sprintf('select * from div where chromosome = %d and animal = "%s" and day = %d;',
                                 refs$chromosome[i], refs$animal[i], refs$day[i]))
    alias <- chr$alias[1]

    tempfilename <- sprintf(args[2],refs$animal[i], refs$day[i], alias)

    pdf(file=tempfilename,height=8.3, width=((length(chr$position) / 10) + 2.5),onefile = FALSE)
    print(c(min(c(chr$A,chr$C,chr$G,chr$T)),max(c(chr$A,chr$C,chr$G,chr$T))))
    par(lty=0)
    barplot(t(matrix(pmax(c(chr$A,chr$C,chr$G,chr$T),0),ncol=4)),
            ylim=c(min(c(chr$A,chr$C,chr$G,chr$T)),max(c(chr$A,chr$C,chr$G,chr$T))),
            xlab="Position",
            ylab="Coverage",
            main=paste("Diversity in",refs$animal[i],"day", refs$day[i],args[3], sep=" "),
            sub=alias,
            names.arg=chr$position,
            col=nuc_colours,
            space=0)

    barplot(t(matrix(pmin(c(chr$A,chr$C,chr$G,chr$T),0),ncol=4)),col=nuc_colours,add=TRUE)

    par(xpd=TRUE)
    legend(x="left",c("A","C","G","T"),fill=nuc_colours)

    ## print(paste("STATS",args[1],name, sum(chr$COV),sum(chr$COV)/length(chr$COV), sep="      ") )

    dev.off()
}
