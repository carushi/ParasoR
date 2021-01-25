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
    print(vec)
    return(vec)
}

a <- read.table("test", sep="\t", header=T)
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

