library(ggplot2)
library(grid)
library(gridExtra)

args <- commandArgs(TRUE)

tabs <- list()
plots <- list()
maxy <- 0

outputname <- args[1]
cds.begin <- strtoi(args[2])
cds.end <- strtoi(args[3])
args <- args[4:length(args)]

for (i in 1:length(args)) {
    print(args[i])
    tab <- read.table(args[i],header=T,sep=",")
    tab$amount <- tab$amount * 100
    tab$amount[tab$amount < 1] <- 1
    maxy <- max(maxy, max(tab$amount))
    tabs[[i]] <- tab
}

coding.region <- data.frame(xstart = cds.begin, xend = cds.end, col = "CDS")

for (i in 1:length(tabs)) {
    p <- ggplot()
    p <- p + guides(fill=FALSE)
    p <- p + theme(axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   plot.margin=unit(c(0,1,0,1),"cm"),
                   text = element_text(size=6))
    p <- p + scale_fill_hue(c=100,l=50)
    p <- p + labs(title=args[[i]])
    p <- p + scale_y_log10()
    p <- p + geom_rect(data = coding.region, aes(xmin = xstart, xmax = xend, ymin = 1, ymax = Inf, fill = col), alpha = 0.2)
    p <- p + geom_bar(data=tabs[[i]], aes(x=position,y=amount,fill=nucleotide),stat="identity")
    plots[[i]] <- p
}

pdf(file=outputname, height=length(args), width=20, onefile=TRUE)
plots <- c(plots, list(ncol=1))
do.call(grid.arrange, plots)
dev.off()
