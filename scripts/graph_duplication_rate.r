# Read car and truck values from tab-delimited autos.dat
args <- commandArgs(trailingOnly = TRUE)
print(args)
autos_data = read.table(args[1], header=T, sep="\t") 

# Compute the largest y value used in the data (or we could
# just use range again)
max_y = max(autos_data)

# Define colors to be used for cars, trucks, suvs
plot_colors = c("blue","red")

# Start PNG device driver to save output to figure.png
png(filename="CRfigure.png", height=295, width=300, 
 bg="white")

# Graph autos using y axis that ranges from 0 to max_y.
# Turn off axes and annotations (axis labels) so we can 
# specify them ourself
plot(autos_data$rate[c(1:30)], log(autos_data$unique[c(1:30)]), type="o", col=plot_colors[1], 
   ylim=log(range(autos_data$unique[c(1:30)])), axes=FALSE, ann=FALSE,)

# Make x axis using Mon-Fri labels
axis(1, at=2*c(1:15), lab=2*c(1:15))

# Make y axis with horizontal labels that display ticks at 
# every 4 marks. 40max_y is equivalent to c(0,4,8,12).
#axis(2, las=1, at=c(1,10,100,1000,100000,100000))
aty <- axTicks(2)
labels <- sapply(aty,function(i)
            as.expression(bquote(10^.(i)))
			)
axis(2,at=aty)
# Create box around plot
box()

# Graph trucks with red dashed line and square points
lines(log(autos_data$non_unique), type="o", pch=22, lty=1, 
   col=plot_colors[2])

# Graph suvs with green dotted line and diamond points
#lines(autos_data$suvs, type="o", pch=23, lty=3, 
#   col=plot_colors[3])

# Create a title with a red, bolditalic font
title(main="Mapping Distribution", col.main="red", font.main=4)

# Label the x and y axes with dark green text
title(xlab= "Rate", col.lab=rgb(0,0.5,0))
title(ylab= "Number of Positions", col.lab=rgb(0,0.5,0))

# Create a legend at (1, max_y) that is slightly smaller 
# (cex) and uses the same line colors and points used by 
# the actual plots
legend(15, max_y, c("unique","non unique"), cex=0.8, col=plot_colors, 
   pch=21:23, lty=13);
   
# Turn off device driver (to flush output to png)
dev.off()