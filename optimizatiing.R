# this procedure is for generating the drug profile of relationships with lncRNAs through unconstrained minimization

# C is coefficient vector for drug d with gene, c has length of number of genes
# P is predicted score matrix decription the relationship between genes and lncRNAs with row as genes, column as lncRNAs
# the optimized variable d is the profile of relationships with lncRNAs, so x has the length of number of lncRNAs

f<-function(d,C,P) sum(abs(C)*unlist(lapply(1:length(C),function(x) sum(unlist(lapply(1:ncol(P),function(y) (sign(C[x])*P[x,y]-d[y])^2))))))

opt<-nlm(f,rep(0,ncol(P)),C,P)

optim_d<-opt$estimate
