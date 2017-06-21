#####################################
##	REPARARTION: Ribosome Profiling Assisted (Re-) Annotation of Bacterial genomes
##
##	Copyright (C) 2017 Elvis Ndah
##
##	This program is free software: you can redistribute it and/or modify
##	it under the terms of the GNU General Public License as published by
##	the Free Software Foundation, either version 3 of the License, or
##	(at your option) any later version.
##	
##	This program is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##	GNU General Public License for more details.
##
##	You should have received a copy of the GNU General Public License
##	along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##	contact: elvis.ndah@gmail.com
#####################################

# R script to plot Profiles


args <- commandArgs(TRUE)
start_file <- as.character(args[1])
stop_file <- as.character(args[2])
work_dir <- as.character(args[3])
count <- as.numeric(args[4])

start <- read.table(file=start_file,sep='\t',h=T)
stop <- read.table(file=stop_file,sep='\t',h=T)

ymax <- max(start$Average_read,stop$Average_read)

cex <- 2
pdf(file=paste(work_dir,"metagene_profile.pdf",sep=""), width=12, height=9)
mar.default = c(5, 4, 4, 2) + 0.1
par(mar=mar.default + c(0,4,0,0)) 
plot(start$position,start$Average_read,col="blue", lwd=2.5, type="l", main = paste("Start Profile ",count), xlab="Distance from start codon", ylab="Average RPF Reads",cex.lab=cex, cex.axis = cex, ylim=c(0,ymax))
abline(v=0,col="red")

plot(stop$position,stop$Average_read, col="blue", lwd=2.5, type="l", main = paste("Stop Profile ",count), xlab="Distance from stop codon", ylab="Average RPF Reads",cex.lab=cex, cex.axis = cex, ylim=c(0,ymax))
abline(v=0,col="red")
invisible(dev.off()) 
<<<<<<< HEAD

=======
>>>>>>> origin/master
