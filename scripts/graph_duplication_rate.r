# Read car and truck values from tab-delimited autos.dat
#args <- commandArgs(trailingOnly = TRUE)
#print(args)

filename <- "P:\\data\\preprocessed\\TEST\\TEST_duplication_rate_reads.txt"

duplication_rate = read.table(filename, header=T, sep="\t") 

# Compute the largest y value used in the data (or we could
# just use range again)

# Define colors to be used for cars, trucks, suvs
plot_colors = c("black","black","darkgray")
line_types = c("dotted","solid","solid");

# Start PNG device driver to save output to figure.png
pdf("duplication_rate.pdf",)

# Graph autos using y axis that ranges from 0 to max_y.
# Turn off axes and annotations (axis labels) so we can 
# specify them ourself
par(mar=c(5, 4, 4, 6) + 0.1)
max_y = max(duplication_rate)
plot(duplication_rate$rate[c(1:30)], duplication_rate$non.unique[c(1:30)]+duplication_rate$unique[c(1:30)], log="y",type="l", lty=line_types[1], col=plot_colors[1], axes=FALSE, ann=FALSE)

axis(1, at=2*(1:15), labels=2*(1:15))

#ylims <- c(0.2, max_y)
#get_axp <- function(x) 10^c(ceiling(x[1]), floor(x[2]))
#usr.i <- log10(ylims)
#(aty <- axTicks(side = 2, usr = usr.i,
                 #axp = c(get_axp(usr.i), n = 3), log = TRUE, nintLog = Inf))
aty <- axTicks(2,log=TRUE)
labels <- sapply(log(aty,10),function(i)
            as.expression(bquote(10^.(i))))
#axis(2, aty, labels = labels)
axis(2, at=axTicks(2,log=TRUE), labels=format(aty, scientific=TRUE))
#axis(2, at=axTicks(2,log=TRUE), labels=labels)

# Graph trucks with red dashed line and square points
lines(duplication_rate$rate[c(1:30)],duplication_rate$unique[c(1:30)],type="l", lty=line_types[2], col=plot_colors[2])
title(xlab= "Duplication Rate")
title(ylab= "Number of Reads")
box()


duplication_rate_difference <-  duplication_rate$non.unique / (duplication_rate$unique + duplication_rate$non.unique)
max_y_difference = max(duplication_rate_difference)
par(new=TRUE)
plot(duplication_rate$rate[c(1:30)],duplication_rate_difference[c(1:30)], type="l", lty=line_types[3], col=plot_colors[3],ann=F, axes=F)
mtext("Percentage PCR artifacts",side=4,col="black",line=4) 

#aty_diff <- seq(0,1,by=0.2)
#labels_diff <- sapply(aty_diff,function(i)
#  as.expression(bquote(10^.(i))))
axis(4, at=axTicks(4), labels=format(axTicks(4), scientific=FALSE),las=1)

#axis(4, ylim=c(min(duplication_rate_difference),max_y_difference))
title(main="Reads per Duplication Rate", col.main="black", font.main=4)

# Label the x and y axes with dark green text


legend("topright", c("non unique","unique", "PCR artifacts"), cex=0.8, col=plot_colors, lty=line_types)
   
# Turn off device driver (to flush output to png)
#dev.off()
