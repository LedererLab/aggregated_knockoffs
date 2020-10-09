# TreeFDR
require(ape)
require(nlme)
require(cluster)
require(StructFDR)

# X
X.tree <- t(scale(X.raw))

# Pairwise distance matrix
rownames(X.tree)[1] <- "V1"
D <- as.matrix(dist(X.tree))
colnames(D) <- rownames(X.tree)
rownames(D) <- rownames(X.tree)

# Define testing and permutation function 
test.func <- function (X, Y) {
  obj <- apply(X, 1, function(x) {
    ttest.obj <- t.test(x ~ Y)
    c(ttest.obj$p.value, sign(ttest.obj$statistic)) })
  return(list(p.value=obj[1, ], e.sign=obj[2, ])) 
}
perm.func <- function (X, Y) { 
  return(list(X=X, Y=sample(Y)))
}
# TreeFDR
tree.fdr.obj <- StructFDR(X.tree, Y, D, test.func, perm.func)
rownames(X.tree)[which(tree.fdr.obj$p.adj <= fdr)]

# BH
BH.p.adj <- p.adjust(tree.fdr.obj$p.unadj, 'fdr')
rownames(X.tree)[which(BH.p.adj <= fdr)]


