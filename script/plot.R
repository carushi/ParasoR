args <- commandArgs(trailingOnly = T)
NUM <- args[1]
for (i in 1:NUM) {
	a <- read.table(paste("stem_", i, ".txt", sep=""), header=F)
	postscript(paste("stem_", i, ".ps", sep=""), horizontal=FALSE)
	plot(1:dim(a)[2], a, type="l", lwd=3, xlab="position", ylab="stem probability")
	dev.off()
}

for (i in 1:NUM) {
	a <- read.table(paste("acc_", i, ".txt", sep=""), header=F)
	postscript(paste("acc_", i, ".ps", sep=""), horizontal=FALSE)
	plot(1:dim(a)[2]*10, a, type="l", lwd=3, xlab="position", ylab="accessibility (kcal/mol)")
	dev.off()
}