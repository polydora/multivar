foram <- read.table("data/Golikova_etal_2020_cluster_data.csv", sep = "\t", header = TRUE)
row.names(foram) <- foram$rowname
library(vegan)
foram_rel <- decostand(foram[, c(12:21)], method = 'total')
foram_ord <- metaMDS(foram_rel)
foram_ord$stress
library(ggplot2)
theme_set((theme_bw(base_size = 10)))
update_geom_defaults('point', list(size = 5))
ggplot(as.data.frame(foram_ord$points), aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = foram$Level, shape = foram$Site))
d <- dist(x = foram_rel, method = "euclidean")

library(ape)
library(dendextend)
hc_single <- hclust(d, method = "single")
ph_single <- as.phylo(hc_single)
c_single <- cophenetic(ph_single)
cor(d, as.dist(c_single))

hc_compl <- hclust(d, method = "complete")
ph_compl <- as.phylo(hc_compl)
c_compl <- cophenetic(ph_compl)
cor(d, as.dist(c_compl))

hc_avg <- hclust(d, method = "average")
ph_avg <- as.phylo(hc_avg)
c_avg <- cophenetic(ph_avg)
cor(d, as.dist(c_avg))

hc_w2 <- hclust(d, method = "ward.D2")
ph_w2 <- as.phylo(hc_w2)
c_w2 <- cophenetic(ph_w2)
cor(d, as.dist(c_w2))

plot(ph_avg, type = "phylogram", cex = 0.7)
axisPhylo()
den_avg <- as.dendrogram(hc_avg)
plot(den_avg)

library(cluster)
complete3 <- cutree(tree = hc_avg, k = 10)
plot(silhouette(x = complete3, dist = d))

library(pvclust)
set.seed(123)
cl_boot <- pvclust(t(foram_rel), method.hclust = "average",
                   nboot = 10000, method.dist = "euclidean",
                   parallel = TRUE, iseed = 42)
plot(cl_boot)

# методы неиерархической кластеризации работают не очень хорошо
set.seed(333)
f_kmeans <- kmeans(foram_rel, centers = 5, nstart = 100)
f_cascade <- cascadeKM(foram_rel, inf.gr = 2, sup.gr = 10,
                       iter = 100, criterion = 'calinski')
plot(f_cascade, sortg = TRUE)

library(factoextra)
fviz_cluster(f_kmeans, data = foram_rel,
             ggtheme = theme_bw())

f_10k <- kmeans(foram_rel, centers = 10,
               nstart = 100)

fviz_cluster(f_10k, data = foram_rel,
             ggtheme = theme_bw())

plot(silhouette(f_kmeans$cluster, d))
plot(silhouette(f_10k$cluster, d))

library(dbscan)
kNNdistplot(foram_rel, k = 3)
abline(h = 0.28, lty = 2)

f_dbscan <- dbscan(foram_rel, 0.28, 3)
fviz_cluster(f_dbscan, data = foram_rel)

f_cmeans <- fanny(d, k = 2, memb.exp = 2)
summary(f_cmeans)

plot(silhouette(f_cmeans), cex.names = 0.6)

library(scatterpie)
f_clust <- cbind(foram_ord$points, as.data.frame(f_cmeans$membership))

colnames(f_clust) <- c("MDS1", "MDS2", "1", "2")
ggplot() + geom_scatterpie(data = f_clust, aes(x = MDS1, y = MDS2),
                           cols = c("1", "2"))
