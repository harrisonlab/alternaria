#!/usr/bin/Rscript

# Plot a 3-way Venn diagram from a tab delimited file containing a matrix showing
 # presence /absence of orthogroups between 3 sets of genomes.

 # This is intended to be used on the output of the orthoMCL pipeline following
 # building of the matrix using:
 # ~/git_repos/emr_repos/tools/Atenogen/orthology/orthoMCL/orthoMCLgroups2tab.py

 # The script also requires the colorspace package. This can be downloaded by
 # opening R and running the following command:
 # options(download.file.method = "wget")
 # install.packages("colorspace")

 #get config options
library(optparse)
library(colorspace)
library(VennDiagram)
opt_list = list(
    make_option("--inp", type="character", help="tab seperated file containing matrix of presence of orthogroups"),
    make_option("--out", type="character", help="output venn diagram in pdf format")
)
opt = parse_args(OptionParser(option_list=opt_list))
f = opt$inp
o = opt$out

orthotabs <-data.frame()
orthotabs <- read.table(f)
df1 <- t(orthotabs)
summary(df1)

df2 <- data.frame()
df2 <- subset(df1, df1[,"Aa_1"] == 1 & df1[,"Aa_2"] == 1 & df1[,"Aa_3"] == 1 | df1[,"At_1"] == 1 & df1[,"At_2"] == 1 & df1[,"At_3"] == 1 & df1[,"At_4"] == 1 | df1[,"At_5"] == 1 & df1[,"At_6"] == 1 & df1[,"At_7"] == 1 & df1[,"At_8"] == 1 | df1[,"Ag_1"] == 1)


Aarb=subset(df1, df1[,"Aa_1"] == 1 & df1[,"Aa_2"] == 1 & df1[,"Aa_3"] == 1 & df1[,"At_1"] == 0 & df1[,"At_2"] == 0 & df1[,"At_3"] == 0 & df1[,"At_4"] == 0 & df1[,"At_5"] == 0 & df1[,"At_6"] == 0 & df1[,"At_7"] == 0 & df1[,"At_8"] == 0 & df1[,"Ag_1"] == 0)
Aten_nonpath=subset(df1, df1[,"Aa_1"] == 0 & df1[,"Aa_2"] == 0 & df1[,"Aa_3"] == 0 & df1[,"At_1"] == 1 & df1[,"At_2"] == 0 & df1[,"At_3"] == 1 & df1[,"At_4"] == 1 & df1[,"At_5"] == 0 & df1[,"At_6"] == 0 & df1[,"At_7"] == 0 & df1[,"At_8"] == 0 & df1[,"Ag_1"] == 0)
Aten_path=subset(df1, df1[,"Aa_1"] == 0 & df1[,"Aa_2"] == 0 & df1[,"Aa_3"] == 0 & df1[,"At_1"] == 0 & df1[,"At_2"] == 0 & df1[,"At_3"] == 0 & df1[,"At_4"] == 0 & df1[,"At_5"] == 1 & df1[,"At_6"] == 1 & df1[,"At_7"] == 1 & df1[,"At_8"] == 1 & df1[,"Ag_1"] == 0)
Agai_path=subset(df1, df1[,"Aa_1"] == 0 & df1[,"Aa_2"] == 0 & df1[,"Aa_3"] == 0 & df1[,"At_1"] == 0 & df1[,"At_2"] == 0 & df1[,"At_3"] == 0 & df1[,"At_4"] == 0 & df1[,"At_5"] == 0 & df1[,"At_6"] == 0 & df1[,"At_7"] == 0 & df1[,"At_8"] == 0 & df1[,"Ag_1"] == 1)

# orthologs=subset(df1, df1[,"Aa_1"] == 1 & df1[,"Aa_2"] == 1 & df1[,"Aa_3"] == 1 & df1[,"At_1"] == 1 & df1[,"At_2"] == 1 & df1[,"At_3"] == 1 & df1[,"At_4"] == 1 & df1[,"At_5"] == 1 & df1[,"At_6"] == 1 & df1[,"At_7"] == 0& df1[,"At_8"] == 0)

# area1=(nrow(Aarb) + nrow(orthologs))
# area2=(nrow(Aten_nonpath) + nrow(orthologs))
# area3=(nrow(Aten_path) + nrow(orthologs))
# area3
# area2
# area1



# Print labels
label1 <- paste('Aarb', sep="" )
label2 <- paste('Aten_nonpath', sep="" )
label3 <- paste('Aten_path', sep="" )
label4 <- paste('Agai_path', sep="" )

#label1 <- paste('Aarb', ' (', area2, ')', sep="" )
#label2 <- paste('Aarb', ' (', area1, ')', sep="" )
#label3 <- paste('FoL', ' (', area3, ')', sep="" )
#label4 <- paste('FoL', ' (', area3, ')', sep="" )

# No labels
#label1 <- paste("", sep="" )
#label2 <- paste("", sep="" )
#label3 <- paste("", sep="" )
#label4 <- paste("", sep="" )

#n123=nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1 & df2[,"At_1"] == 1 & df2[,"At_2"] == 1 & df2[,"At_3"] == 1 & df2[,"At_4"] == 1 & df2[,"At_5"] == 1 & df2[,"At_6"] == 1))
#n12=n123 + nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1 & df2[,"At_1"] == 1 & df2[,"At_2"] == 1 & df2[,"At_3"] == 1 & df2[,"At_4"] == 1 & df2[,"At_5"] == 1 & df2[,"At_6"] == 0))
#n13=n123 + nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1 & df2[,"At_1"] == 1 & df2[,"At_2"] == 1 & df2[,"At_3"] == 0 & df2[,"At_4"] == 0 & df2[,"At_5"] == 0 & df2[,"At_6"] == 1))
#n23=n123 + nrow(subset(df2, df2[,"Aa_1"] == 0 & df2[,"Aa_2"] == 0 & df2[,"Aa_3"] == 0 & df2[,"At_1"] == 0 & df2[,"At_2"] == 0 & df2[,"At_3"] == 1 & df2[,"At_4"] == 1 & df2[,"At_5"] == 1 & df2[,"At_6"] == 1))

n12=nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1 & df2[,"At_1"] == 0 & df2[,"At_2"] == 0 & df2[,"At_3"] == 0 & df2[,"At_4"] == 0 & df2[,"At_5"] == 0 & df2[,"At_6"] == 0 & df2[,"At_7"] == 0 & df2[,"At_8"] == 0 & df2[,"Ag_1"] == 0))
n13=nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1 & df2[,"At_1"] == 0 & df2[,"At_2"] == 0 & df2[,"At_3"] == 0 & df2[,"At_4"] == 0 & df2[,"At_5"] == 1 & df2[,"At_6"] == 1 & df2[,"At_7"] == 1 & df2[,"At_8"] == 1 & df2[,"Ag_1"] == 0))
n14=nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1 & df2[,"At_1"] == 0 & df2[,"At_2"] == 0 & df2[,"At_3"] == 0 & df2[,"At_4"] == 0 & df2[,"At_5"] == 0 & df2[,"At_6"] == 0 & df2[,"At_7"] == 0 & df2[,"At_8"] == 0 & df2[,"Ag_1"] == 1))
n23=nrow(subset(df2, df2[,"Aa_1"] == 0 & df2[,"Aa_2"] == 0 & df2[,"Aa_3"] == 0 & df2[,"At_1"] == 1 & df2[,"At_2"] == 1 & df2[,"At_3"] == 1 & df2[,"At_4"] == 1 & df2[,"At_5"] == 1 & df2[,"At_6"] == 1 & df2[,"At_7"] == 1 & df2[,"At_8"] == 1 & df2[,"Ag_1"] == 0))
n24=nrow(subset(df2, df2[,"Aa_1"] == 0 & df2[,"Aa_2"] == 0 & df2[,"Aa_3"] == 0 & df2[,"At_1"] == 1 & df2[,"At_2"] == 1 & df2[,"At_3"] == 1 & df2[,"At_4"] == 1 & df2[,"At_5"] == 0 & df2[,"At_6"] == 0 & df2[,"At_7"] == 0 & df2[,"At_8"] == 0 & df2[,"Ag_1"] == 1))
n34=nrow(subset(df2, df2[,"Aa_1"] == 0 & df2[,"Aa_2"] == 0 & df2[,"Aa_3"] == 0 & df2[,"At_1"] == 0 & df2[,"At_2"] == 0 & df2[,"At_3"] == 0 & df2[,"At_4"] == 0 & df2[,"At_5"] == 1 & df2[,"At_6"] == 1 & df2[,"At_7"] == 1 & df2[,"At_8"] == 1 & df2[,"Ag_1"] == 1))
n123=nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1 & df2[,"At_1"] == 1 & df2[,"At_2"] == 1 & df2[,"At_3"] == 1 & df2[,"At_4"] == 1 & df2[,"At_5"] == 1 & df2[,"At_6"] == 1 & df2[,"At_7"] == 1 & df2[,"At_8"] == 1 & df2[,"Ag_1"] == 0))
n124=nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1 & df2[,"At_1"] == 1 & df2[,"At_2"] == 1 & df2[,"At_3"] == 1 & df2[,"At_4"] == 1 & df2[,"At_5"] == 0 & df2[,"At_6"] == 0 & df2[,"At_7"] == 0 & df2[,"At_8"] == 0 & df2[,"Ag_1"] == 1))
n134=nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1 & df2[,"At_1"] == 0 & df2[,"At_2"] == 0 & df2[,"At_3"] == 0 & df2[,"At_4"] == 0 & df2[,"At_5"] == 1 & df2[,"At_6"] == 1 & df2[,"At_7"] == 1 & df2[,"At_8"] == 1 & df2[,"Ag_1"] == 1))
n234=nrow(subset(df2, df2[,"Aa_1"] == 0 & df2[,"Aa_2"] == 0 & df2[,"Aa_3"] == 0 & df2[,"At_1"] == 1 & df2[,"At_2"] == 1 & df2[,"At_3"] == 1 & df2[,"At_4"] == 1 & df2[,"At_5"] == 1 & df2[,"At_6"] == 1 & df2[,"At_7"] == 1 & df2[,"At_8"] == 1 & df2[,"Ag_1"] == 1))
n1234=nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1 & df2[,"At_1"] == 1 & df2[,"At_2"] == 1 & df2[,"At_3"] == 1 & df2[,"At_4"] == 1 & df2[,"At_5"] == 1 & df2[,"At_6"] == 1 & df2[,"At_7"] == 1 & df2[,"At_8"] == 1 & df2[,"Ag_1"] == 1))


summary(n12)
summary(n13)
summary(n23)
summary(n123)

#
#area1=(nrow(subset(df2, df2[,"Aa_1"] == 1 & df2[,"Aa_2"] == 1 & df2[,"Aa_3"] == 1)))
#area2=(nrow(subset(df2, df2[,"At_1"] == 1 & df2[,"At_2"] == 1 & df2[,"At_3"] == 1 & df2[,"At_4"] == 1)))
#area3=(nrow(subset(df2, df2[,"At_5"] == 1 & df2[,"At_6"] == 1 & df2[,"At_7"] == 1 & df2[,"At_8"] == 1)))
#area4=(nrow(subset(df2, df2[,"Ag_1"] == 1)))
#+n12 +n13 +n14 +n23 +n24 +n34 +n123 +n124 +n134 +n234 +n1234

area1=(nrow(Aarb) +n12 +n13 +n14 +n123 +n124 +n134 +n1234)
area2=(nrow(Aten_nonpath) +n12 +n23 +n24 +n123 +n124 +n234 +n1234)
area3=(nrow(Aten_path) +n13 +n23 +n34 +n123 +n134 +n234 +n1234)
area4=(nrow(Agai_path) +n14 +n24 +n34 +n124 +n134 +n234 +n1234)


#nrow(Aarb)
nrow(Aarb)
nrow(Aten_nonpath)
nrow(Aten_path)
n12
n13
n23
n123
area1
#area1 - n12 - n13 + n123
area2
area3

 +n123 +n124 +n134 +n234 +n1234

pdf(o)
draw.quad.venn(area1, area2, area3, area4,
    (n12 +n123 +n124 +n1234), (n13 +n123 +n134 +n1234), (n14 +n124 +n134 +n1234), (n23 +n123 +n234 +n1234), (n24 +n124 +n234 +n1234), (n34 +n134 +n234 +n1234),
    (n123 +n1234), (n124 +n1234), (n134 +n1234), (n234 +n1234),
    n1234,
    category = c(label1, label2, label3, label4),
#    rep("", 4),
    rotation = 1,
    reverse = FALSE,
    lwd = rep(2, 4),
    lty = rep("solid", 4),
    col = rep("black", 4),
    fill = c(rainbow_hcl(4)),
    alpha = rep(0.5, 4),
    label.col = rep("black", 15),
    cex = rep(1, 15),
    fontface = rep("plain", 15),
    fontfamily = rep("serif", 15),
    cat.pos = c(-15, 15, 0, 0),
    cat.dist = c(0.22, 0.22, 0.11, 0.11),
    cat.col = rep("black", 4),
    cat.cex = rep(1, 4),
    cat.fontface = rep("plain", 4),
    cat.fontfamily = rep("serif", 4),
    cat.just = rep(list(c(0.5, 0.5)), 4),
    cat.default.pos = "outer",
    cat.prompts = FALSE,
    rotation.degree = 0,
    rotation.centre = c(0.5, 0.5),
    ind = TRUE, sep.dist = 0.05, offset = 0,
    )


dev.off()
singles = df1[grepl("single*", rownames(df1)), ]
print("Aa_1")
total_1 = nrow(subset (df1, df1[,"Aa_1"] == 1))
missing_1 = (total_1 - area2)
uniq_1=sum(singles[, "Aa_1"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_1)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_1)
paste('The total number of singleton genes not in the venn diagram: ', uniq_1)
print("Aa_2")
total_2 = nrow(subset (df1, df1[,"Aa_2"] == 1))
missing_2 = (total_2 - area2)
uniq_2=sum(singles[, "Aa_2"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_2)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_2)
paste('The total number of singleton genes not in the venn diagram: ', uniq_2)
print("Aa_3")
total_3 = nrow(subset (df1, df1[,"Aa_3"] == 1))
missing_3 = (total_3 - area2)
uniq_3=sum(singles[, "Aa_3"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_3)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_3)
paste('The total number of singleton genes not in the venn diagram: ', uniq_3)
print("At_3")
total_4 = nrow(subset (df1, df1[,"At_3"] == 1))
missing_4 = (total_4 - area1)
uniq_4=sum(singles[, "At_3"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_4)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_4)
paste('The total number of singleton genes not in the venn diagram: ', uniq_4)
print("At_4")
total_5 = nrow(subset (df1, df1[,"At_4"] == 1))
missing_5 = (total_5 - area1)
uniq_5=sum(singles[, "At_4"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_5)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_5)
paste('The total number of singleton genes not in the venn diagram: ', uniq_5)
print("At_5")
total_6 = nrow(subset (df1, df1[,"At_5"] == 1))
missing_6 = (total_6 - area1)
uniq_6=sum(singles[, "At_5"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_6)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_6)
paste('The total number of singleton genes not in the venn diagram: ', uniq_6)
print("At_6")
total_7 = nrow(subset (df1, df1[,"At_6"] == 1))
missing_7 = (total_6 - area3)
uniq_7=sum(singles[, "At_6"])
paste('The total number of orthogroups and singleton genes in this isolate: ', total_6)
paste('The total number of orthogroups and singleton genes not in the venn diagram: ', missing_6)
paste('The total number of singleton genes not in the venn diagram: ', uniq_6)

#inpara_2 = sum(orthogroups[,"Aa_1"] == 0 & orthogroups[,"Aa_2"] == 1)
#label1
#uniq_1
#inpara_1
#label2
#uniq_2
#inpara_2

warnings()
q()
