library(ggplot2)
library(grid)
library(gridExtra)

args <- commandArgs(TRUE)

tabs <- list()
plots <- list()
maxy <- 0

outputname <- args[1]
args <- args[2:length(args)]

for (i in 1:length(args)) {
    print(args[i])
    tab <- read.table(args[i],header=T,sep=",")
    tab$amount <- tab$amount * 100
    tab$amount[tab$amount < 1] <- 1
    maxy <- max(maxy, max(tab$amount))
    tabs[[i]] <- tab
}

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
    p <- p + geom_bar(data=tabs[[i]], aes(x=position,y=amount,fill=nucleotide),stat="identity")
    ## p <- p + scale_y_continuous(limits=c(0,maxy))
    p <- p + scale_y_log10()
    plots[[i]] <- p
}

pdf(file=outputname, height=length(args), width=20, onefile=TRUE)
plots <- c(plots, list(ncol=1))
do.call(grid.arrange, plots)
dev.off()
