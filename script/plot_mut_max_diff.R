library(dplyr)
library(ggpubr)

check.trans <- function(data) {
    vec <- apply(data, c(1), function(x){
        mut = c(x[1], x[3])
        if (mut[1] > mut[2])
            mut = c(mut[2], mut[1])
        if (mut[1] == 1) {
            if (mut[2] == 3) return('transition')
            return('transversion')
        } else if (any(mut == 3)) {
            return('transversion')
        } else {
            if (mut[1] == 2 && mut[2] == 4)
                return('transition')
            else
                return('transversion')
        }
    })
    return(vec)
}
compute.pvalue <- function(data) {
    index=order(data[,4], decreasing=TRUE)
    data <- cbind(data[index,], 1:dim(data)[1]/dim(data)[1])
    print(min(data[,4]))
    print(max(data[,4]))
    data <- data[order(data[,2], data[,3]),]
    return(as.vector(data[,dim(data)[2]]))
}
compute.pvalue.base <- function(data) {
    vec <- NULL
    for (base_type in 1:4) {
        part <- data[data[,1] == base_type,]
        index=order(part[,4], decreasing=TRUE)
        print(c(base_type, min(part[,4]), max(part[,4])))
        temp <- cbind(part[index,], p.value.base=1:dim(part)[1]/dim(part)[1])
        if (is.null(vec)) {
            vec <- temp
        } else {
            vec <- rbind(vec, temp)
        }
    }
    vec <- vec[order(vec[,2], vec[,3]),]
    return(as.vector(vec[,'p.value.base']))
}
input=commandArgs(trailingOnly=TRUE)
a <- read.table(input, sep="\t", header=T)
before <- unlist(sapply(1:max(a[,1]), function(x){
    temp <- a[a[,1] == x,2]
    return(rep(c(1:4)[unlist(sapply(1:4, function(x){return(!any(x == temp))}))][1], 3))
}))

a <- cbind(before=as.vector(before), a)
trans <- check.trans(a)
a <- cbind(a, mutation_type=trans)
bases=c('A', 'C', 'G', 'U')
a <- cbind(a, substitution=unlist(apply(a, c(1), function(x){return(paste0(bases[as.integer(x[1])], '->', bases[as.integer(x[3])]))})))

pdf('mutation_max_diff.pdf')
g <- ggboxplot(a, x='substitution', y='max_diff', fill='mutation_type', palette=c('grey', 'white'), ylab="Maximum difference of stem probabilities")
plot(g)
dev.off()

colnames(a)[dim(a)[2]-2]="Pearson correlation coefficient"
pdf('mutation_cor_diff.pdf')
g <- ggboxplot(a, x='substitution', y='Pearson correlation coefficient', fill='mutation_type', palette=c('grey', 'white'))
plot(g)
dev.off()


a <- cbind(a, p.value=compute.pvalue(a))
a <- cbind(a, p.value.base=compute.pvalue.base(a))

a[,1] <- unlist(sapply(a[,1], function(x){return(bases[x])}))
pdf('compare_pvalue_for_all_and_each_base.pdf')
g <- ggscatter(a, x='p.value', y='p.value.base', color='before', palette=rainbow(4))
plot(g)
dev.off()
pdf('compare_pvalue_for_all_and_each_base.pdf')
g <- ggscatter(a, x='p.value', y='p.value.base', color='mutation_type', palette=rainbow(16))
plot(g)
dev.off()
write.table(a, file='pvalue.tsv', sep="\t")

