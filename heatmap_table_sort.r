#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
infile = paste("", args[1], sep="")
outfile = paste("", args[2], sep="")
if(nchar(infile) >= 1 && !infile == "NA" && nchar(outfile) >= 1 && !outfile == "NA") {
    x = read.table(infile, sep="\t", header=T)
    xm1 = as.matrix(x[,-1])
    rownames(xm1) = x[,1]
    xm = xm1
    xm[is.na(xm)] = 0
    hcr = hclust(dist(xm), "ave")
    hcc = hclust(dist(t(xm)), "ave")
    xmo = xm1[hcr$order,hcc$order]
    write.table(xmo, outfile, quote=F, sep="\t", row.names=T, col.names=T, na="NA")
} else {
    warning("Usage: heatmap_table_sort_row_col.r [input table] [output sorted table]")
}
