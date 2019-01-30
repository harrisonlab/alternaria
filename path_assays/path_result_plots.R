

df1 <- read.csv("~/Downloads/14d_Alt_data.csv")
View(df1)

library(ggplot2)
library(readr)
#setwd('~/Downloads/Aalt')


# --- Summary SE function ---
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
# --- /Summary SE function ends ---


df1$name <- gsub("-ve.*","control",df1$name)
df1$Leaf.type <- gsub("Bramley","Bramley's seedling",df1$Leaf.type)
df1$Leaf.type <- factor(df1$Leaf.type, levels=c("Spartan", "Bramley's seedling"))


df2 <- summarySE(df1, measurevar="no..leisions", groupvars=c("name", 'Leaf.type'))
df2$name <- factor(df2$name, levels=c("635", "743", "1166", "648", "1082", "1164", "control"))


# Create plot
p<-ggplot(data=subset(df2), aes(x=name, y=no..leisions))
p <- p + geom_bar(stat="identity")
p <- p + theme(axis.text.x=element_text(angle = -45, hjust = 0))
p <- p + ylab('Number of lesions') + xlab('')
p <- p + geom_errorbar(aes(ymin=no..leisions-se, ymax=no..leisions+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
p <- p + theme(plot.margin=unit(c(1,3,0.5,0.5),"cm"))
p <- p + facet_grid(Leaf.type ~ .)
p


df3 <- summarySE(df1, measurevar="no..leisions", groupvars=c("Treatment", 'Leaf.type'))
df3$name[df3$Treatment == 1] <- "apple pathotype"
df3$name[df3$Treatment == 2] <- "non-pathotype"
df3$name[df3$Treatment == 3] <- "control"
df3$name <- factor(df3$name, levels=c("apple pathotype", "non-pathotype", "control"))
p<-ggplot(data=subset(df3), aes(x=name, y=no..leisions))
p <- p + geom_bar(stat="identity")
p <- p + theme(axis.text.x=element_text(angle = -45, hjust = 0))
p <- p + ylab('Number of lesions') + xlab('')
p <- p + geom_errorbar(aes(ymin=no..leisions-se, ymax=no..leisions+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
p <- p + theme(plot.margin=unit(c(1,3,0.5,0.5),"cm"))
p <- p + facet_grid(Leaf.type ~ .)
p


df2$level <- "by isolate"
df3$level <- "by pathotype"
df2$sig <- 'c'
df2$sig[df2$name == '635'] <- 'a'
df2$sig[df2$name == '743'] <- 'b'
df2$sig[df2$name == '1166'] <- 'a'

df3$sig[df3$Treatment == 1] <- 'a'
df3$sig[df3$Treatment == 2] <- 'b'
df3$sig[df3$Treatment == 3] <- 'b'

df3 <- df3[,-which(names(df3) == "Treatment")]

df4 <- rbind(df2,df3)
df4$name <- factor(df4$name, levels=c("apple pathotype", "non-pathotype", "635", "743", "1166", "648", "1082", "1164", "control"))
df4$name <- revalue(df4$name, c("635"="FERA 635", "743"="FERA 743", "1166"="FERA 1166", "648"="FERA 648", "1082"="FERA 1082", "1164"="FERA 1164"))
df4$name <- factor(df4$name, levels=c("apple pathotype", "non-pathotype", "FERA 635", "FERA 743", "FERA 1166", "FERA 648", "FERA 1082", "FERA 1164", "control"))
df4$sig[df4$Leaf.type == 'Spartan']
df4$sig[df4$Leaf.type == "Bramley's seedling"]

# df4 <- df4[order(df4$name, df4$Leaf.type, df4$level),]
p<-ggplot(data=subset(df4), aes(x=name, y=no..leisions))
p <- p + geom_bar(stat="identity")
p <- p + theme(axis.text.x=element_text(angle = -45, hjust = 0))
p <- p + scale_y_continuous(limits = c(0, 6))
p <- p + ylab('Number of lesions') + xlab('')
p <- p + geom_errorbar(aes(ymin=no..leisions-se, ymax=no..leisions+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
# p <- p + geom_text(aes(x=name, y=1.2 * max(no..leisions), label=df4$sig),
#                                    data=)
 p <- p + geom_text(size = 3,aes(x=name, y=5.9,
                  label=c("c", "c", "a", "a", "c", "b", "c",
                    "a", "b", "b",
                    "c", "c", "a", "a", "c", "b", "c",
                      "a", "b", "b")),
                                    data=df4)
p <- p + facet_grid(Leaf.type ~ level, scales="free_x")
# p <- p + facet_grid(Leaf.type ~ level)
p

# filename <- '~/Downloads/Aalt/Fig5_path_results.pdf'
filename <- '~/Downloads/Aalt/Fig5_path_results.tiff'
ggsave(filename, plot = p, width =10, height = 10, units = "cm", limitsize = FALSE)
