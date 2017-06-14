#####################################
##	REPARARTION: Ribosome Profiling Assisted (Re-) Annotation of Bacterial proteome
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

# R script for random forest ORF delineation

# functions
run_randomforest <- function(train,pos,neg){
	set.seed(1)	# set seed
	rf.output <- randomForest(class~., data=train, mtry=4, ntree=3001, nodesize=4, classwt = c(P=pos, N=neg))
	return(rf.output)
}

# Library
suppressMessages(library(ggplot2))
suppressMessages(library(randomForest))
suppressMessages(library(ROCR))
suppressMessages(library(PRROC))
suppressMessages(library(SiZer))

# Input varibles
args <- commandArgs(TRUE)
ORFs_file <- as.character(args[1])
Positive_set <- as.character(args[2])
work_dir <- as.character(args[3])
codons <- as.character(args[4])
ncodon <- as.character(args[5])
sd <- as.character(args[6])
cutoff <- as.numeric(args[7])

# split string to get start codons
valid_start = strsplit(codons, ",")[[1]]
ncodons = strsplit(ncodon, ",")[[1]]

# ORFs
ORFs <- read.table(file=ORFs_file, sep="\t", h=T)
ORFs <- ORFs[which(ORFs$start_rpkm > 0),]
ORFs$log_rpkm <- log(ORFs$rpkm)

# Positive Set
prodigal <- read.table(file=Positive_set, sep="\t", h=T)
positive <- ORFs[which(as.vector(ORFs$orf_id) %in% as.vector(prodigal$orf_id)),]
positive <- positive[which(as.vector(positive$start_codon) %in% valid_start),]
positive$class <- "P"

# S-CURVE Threshold Estimation
model <- nls(coverage ~ SSfpl(log_rpkm, D, A, C, B), data=positive)
x <- fitted(model)
y <- predict(model)

# SiZer bent point prediction
bent <- bent.cable(x, y, grid.size=100)
MINRPKM <- round(bent$alpha, digits = 2)
MINCOV <- round(predict(model, newdata = data.frame(log_rpkm = MINRPKM))[1], digits = 2)

scurve <- paste(work_dir,"S_Curve.pdf",sep="")
pdf(file=scurve, width=10, height=7)
par(cex.main = 0.75)
head <- paste(paste(paste("Min Log RPKM ", MINRPKM), " Min Coverage "), MINCOV)
plot(positive$coverage ~ positive$log_rpkm, main=head, xlab="Log Read density", ylab="RPF Coverage", cex.lab=1.5, cex.axis=1.25)
points(positive$log_rpkm,predict(model),lty=0.5, cex=0.5, pch=20,col='red')
abline(h=MINCOV, col="blue")
abline(v=MINRPKM, col="blue")
dev.off()

cat("Minimum Log Read density ", MINRPKM,"\n", sep=" ")
cat("Minimum ORF coverage ", MINCOV,"\n", sep=" ")

# keep all ORFs above threshold
ORFs <- ORFs[which(ORFs$log_rpkm >= MINRPKM),]
ORFs <- ORFs[which(ORFs$coverage >= MINCOV),]

# positive set thresholding
positive <- positive[which(positive$log_rpkm >= MINRPKM),]
positive <- positive[which(positive$coverage >= MINCOV),]
cat("Total number of ORFs in positive set ",dim(positive)[1],"\n", sep=" ")

min.pos <- min(positive$length)

# Negative Set
negative.all <- ORFs[which(ORFs$start_codon == ncodon),]
negative.all <- negative.all[which(negative.all$length >= min.pos),]
negative.all$class <- "N"

# keep longest ORF
family.neg <- unique(as.vector(negative.all$gene))
family.neg.vector <- character(length = length(family.neg))
for (i in 1:length(family.neg)) {
	g <- family.neg[i]
	family <- negative.all[which(negative.all$gene == g),]
	longest.orf <- family[which(family$length == max(family$length)),]
	family.neg.vector[i] <- as.vector(longest.orf$orf_id)[1]
}
negative <- negative.all[which(as.vector(negative.all$orf_id) %in% family.neg.vector),]
cat("Number of ORFs in negative Set ",nrow(negative),"\n", sep=" ")

# Number of Valid ORFs for prediction
ORFs <- ORFs[which(ORFs$start_codon %in% valid_start),]
cat("Total number of valid ORFs for prediction ",nrow(ORFs),"\n", sep=" ")
cat("Total number of valid ORF families ",length(unique(as.vector(ORFs$gene))),"\n",sep=" ")

# Training Set
trainset.rf <- rbind(positive,negative)
trainset.rf$class <- as.factor(trainset.rf$class)

if (sd == 'N') {
	feat <- c("class","start_coverage","start_rpkm","coverage","stop_rpkm","accumulation_proportion")
} else {
	feat <- c("class","start_coverage","start_rpkm","coverage","stop_rpkm","accumulation_proportion","SD_score") 
}

trainset <- trainset.rf[,feat]

#	RandomForest
pos.prob <- nrow(positive)/nrow(negative)
neg.prob <- 1 - pos.prob

if (pos.prob < 0.25) {
	pos.prob <- 0.25
	neg.prob <- 1 - pos.prob
} 

if (neg.prob < 0.25) {
	neg.prob <- 0.25
	pos.prob <- 1 - neg.prob
}

# run random forest model
rf_output <- run_randomforest(trainset,pos.prob,neg.prob)

preds.train = predict(rf_output, type="response")
conf.matrix <- table(trainClass=trainset$class,predClass=preds.train)

TP <- conf.matrix[2,2]
FP <- conf.matrix[1,2]
FN <- conf.matrix[2,1]

precision <- TP/(TP + FP)
recall <- TP/(TP + FN)
cat("Precision ",precision,"\n",sep=" ")
cat("Recall ",recall,"\n",sep=" ")

# Variable Importance
vimp <- paste(work_dir,"variable_importance.pdf",sep="")
pdf(file=vimp)
varImpPlot(rf_output, type=2, n.var=length(feat), scale=FALSE, main="Variable Importance (Gini) predictors")
dev.off()

# Area Under the Curve
ROC_curve <- paste(work_dir,"PR_and_ROC_curve.pdf",sep="")
pdf(file=ROC_curve)
prediction.rf <- as.vector(rf_output$votes[,2])
pred.rf=prediction(prediction.rf,trainset[,1])
perf.rf <- performance(pred.rf,"tpr","fpr")
auc.rf <- performance(pred.rf, measure = "auc")
auc.rf <- auc.rf@y.values[[1]]
roc.data.rf <- data.frame(fpr=unlist(perf.rf@x.values),tpr=unlist(perf.rf@y.values),model="GLM")
ggplot(roc.data.rf, aes(x=fpr, ymin=0, ymax=tpr)) +
    geom_ribbon(alpha=0.2) +
    geom_line(aes(y=tpr)) +
    ggtitle(paste0("ROC Curve w/ AUC=", auc.rf)) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold")
	)

# precision
pred.rf=prediction(prediction.rf,trainset$class)
prec.rf <- performance(pred.rf, "prec", "rec")
plot(prec.rf, colorize=T, cex.lab=1.75, cex.axis=2, cex=2)

dev.off()


# PREDICTION
ORFs$pred <- predict(rf_output, ORFs, type="response")
ORFs$prob <- predict(rf_output, ORFs, type="prob")[,2]

head <- c("orf_id","gene","strand","length","start_codon","count","rpkm","coverage","SD_score","SD_pos","prob")
result <- ORFs[which(ORFs$prob >= cutoff),head]

cat("Total number of ORF families predicted ",length(unique(as.vector(result$gene))),"\n",sep=" ")

result_file <- paste(work_dir,"tmp/RF_predicted_ORFs.txt",sep="")
write.table(result, file=result_file, sep = "\t",col.names = TRUE,row.names = F, quote = FALSE)

result_file_all <- paste(work_dir,"tmp/RF_predicted_all.txt",sep="")
write.table(ORFs, file=result_file_all, sep = "\t",col.names = TRUE,row.names = F, quote = FALSE)


# export thresholds
cutoff <- c(MINCOV,exp(MINRPKM))
threshold <- paste(work_dir,"tmp/threshold.txt",sep="")
write(cutoff,threshold,sep="\n")



