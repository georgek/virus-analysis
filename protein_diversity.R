args <- commandArgs(TRUE)
args

databaseFile <- args[1]
outputFileFormat <- args[2]   # sprintf(args[2],samples$animal[sid], samples$day[sid], alias)
								# "%s_d%d_%s.pdf"
titleAddition <- args[3]

library("DBI")
library("RSQLite")
library("reshape2")
library("ggplot2")

db <- dbConnect(dbDriver("SQLite"), dbname=databaseFile)

print("Translating codons...")
result <- dbGetQuery(db, "create temp table codon_translation (codon, aminoacid);")
result <- dbGetQuery(db, "insert into codon_translation (codon, aminoacid) values
('TTT','F'),('TCT','S'),('TAT','Y'),('TGT','C'),('TTC','F'),('TCC','S'),('TAC','Y'),
('TGC','C'),('TTA','L'),('TCA','S'),('TAA','.'),('TGA','.'),('TTG','L'),('TCG','S'),
('TAG','.'),('TGG','W'),('CTT','L'),('CCT','P'),('CAT','H'),('CGT','R'),('CTC','L'),
('CCC','P'),('CAC','H'),('CGC','R'),('CTA','L'),('CCA','P'),('CAA','Q'),('CGA','R'),
('CTG','L'),('CCG','P'),('CAG','Q'),('CGG','R'),('ATT','I'),('ACT','T'),('AAT','N'),
('AGT','S'),('ATC','I'),('ACC','T'),('AAC','N'),('AGC','S'),('ATA','I'),('ACA','T'),
('AAA','K'),('AGA','R'),('ATG','M'),('ACG','T'),('AAG','K'),('AGG','R'),('GTT','V'),
('GCT','A'),('GAT','D'),('GGT','G'),('GTC','V'),('GCC','A'),('GAC','D'),('GGC','G'),
('GTA','V'),('GCA','A'),('GAA','E'),('GGA','G'),('GTG','V'),('GCG','A'),('GAG','E'),
('GGG','G');")

## translate codons
result <- dbGetQuery(db, "create temp view proteinpileup as
select animal, day, protein_id, position, aminoacid, sum(count) as count from codonpileup
left join codon_translation on codonpileup.codon = codon_translation.codon
group by animal,day,protein_id,position,aminoacid;")

## find overall consensus proteins
print("Finding consensus proteins...")
result <- dbGetQuery(db, "create temp table consensus as
select protein_id,position,aminoacid,max(count) as count from
(select protein_id,position,aminoacid,sum(count) as count
 from proteinpileup group by protein_id,position,aminoacid)
group by protein_id,position;")

animals <- dbGetQuery(db, "select distinct animal from codonpileup;")$animal

for (i in 1:length(animals)) {

    days <- dbGetQuery(db, sprintf('select distinct day from codonpileup where animal = "%s";',animals[i]))$day

    for (j in 1:length(days)) {
        proteins <- dbGetQuery(db, sprintf('select distinct protein_id,name from codonpileup
left join proteins on (chromosome = protein_id)
where animal = "%s" and day = %d;',
                                            animals[i], days[j]))
        names <- proteins$name
        proteids <- proteins$protein_id
        dat.cov <- list()
        dat.var <- list()
        lengths <- list()

        ## set up data frames
        for (k in 1:length(proteids)) {
            print(sprintf("Animal: %s, day: %d, protein: %s", animals[i], days[j], names[k]))
            result <- dbGetQuery(db, sprintf('create temp table agree as
select p.position, p.aminoacid, p.count from proteinpileup as p
join consensus as c on p.protein_id = c.protein_id and p.position = c.position and p.aminoacid = c.aminoacid
where animal = "%s" and day = %d and p.protein_id = %d;', animals[i], days[j], proteids[k]))

            df <- dbGetQuery(db, "select * from agree;")
            df$position <- (df$position - min(df$position))/3
            dat.cov[[k]] <- df

            df <- dbGetQuery(db, sprintf("select position, aminoacid, count from proteinpileup
where animal = '%s' and day = %d and protein_id = %d except select * from agree;",
                                                  animals[i], days[j], proteids[k]))
            df$position <- (df$position - min(df$position))/3
            df$count <- -df$count
            dat.var[[k]] <- df

            lengths[[k]] <- max(df$position)

            result <- dbGetQuery(db, "drop table agree;")
        }

        ## set up PDF
        maxlength <- max(unlist(lengths))
        filename <- sprintf(outputFileFormat, animals[i], days[j])
        print(paste("Writing ", filename, "..."))
        pdf(file=filename, height=8.3, width=(maxlength / 10), onefile=TRUE)

        ## print plots
        for (k in 1:length(proteids)) {
            c <- ggplot()
            c <- c + geom_bar(data=dat.cov[[k]], aes(x=position, y=count, fill=aminoacid), stat="identity")
            c <- c + geom_bar(data=dat.var[[k]], aes(x=position, y=count, fill=aminoacid), stat="identity")
            c <- c + scale_x_discrete(expand=c(0, (maxlength-lengths[[k]])/2 + 1))
            c <- c + labs(title = sprintf("Diversity in %s day %d, %s protein", animals[i], days[j], names[k]),
                          x = "Position", y = "Coverage", fill = "Base")
            c <- c + theme(axis.text.x = element_text(size=8,angle=90, hjust=1))
            suppressWarnings(print(c))  # ggplot warns about negative bars, but it's ok
        }
        dev.off()
    }
}

dbDisconnect(db)
