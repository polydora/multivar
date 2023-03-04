##### Загружаем данные ####
library(candisc)
data("Wolves")


##### Анализ данных ####
dim(Wolves)
colnames(Wolves)
head(rownames(Wolves))
any(is.na(Wolves))
table(Wolves$group)


##### Построение MDS-ординации ####
library(vegan)
library(ggplot2)
theme_set(theme_bw(base_size = 16))
update_geom_defaults("point", list(size = 4))
update_geom_defaults("text", list(size = 5))
st_w <- scale(Wolves[, ]) ## стандартизируем
ord_w <- metaMDS()
dfr_w <- data.frame(, Group = Wolves$group)
gg_w <- ggplot(dfr_w, aes()) +
  geom_point()


##### Пакеты для визуализации кластеризации ####
library(ape)
library(dendextend)

# Матрица расстояний
d <- dist(x = st_w, method = "euclidean")


##### Метод ближайшего соседа ####

# Базовая графика
hc_single <- hclust(d, method = "single")
plot(hc_single)

# С помощью пакета ape
ph_single <- as.phylo(hc_single)
plot(ph_single, type = "phylogram")
axisPhylo()


# С помощью dendextend
den_single <- as.dendrogram(hc_single)
plot(den_single)


##### Метод отдалённого соседа ####

# Ape
hc_compl <- hclust(d, method = "complete")
ph_compl <- as.phylo(hc_compl)
plot(ph_compl, type = "phylogram")
axisPhylo()


# Dendextend
den_compl <- as.dendrogram(hc_compl)
plot(den_compl)


##
hc_avg <- hclust(d, method = "average")
ph_avg <- as.phylo(hc_avg)
plot(ph_avg, type = "phylogram")
axisPhylo()


##
den_avg <- as.dendrogram(hc_avg)
plot(den_avg)


##
hc_w2 <-hclust(d, method = "ward.D2")
ph_w2 <- as.phylo(hc_w2)
plot(ph_w2, type = "phylogram")
axisPhylo()


##
den_w2 <- as.dendrogram(hc_w2)
plot(den_w2)


##
# Матрица кофенетических расстояний
c_single <- cophenetic(ph_single)

# Кофенетическая корреляция =
# = корреляция матриц кофенетич. и реальн. расстояний
cor(d, as.dist(c_single))


##
cor(d, as.dist(c_single))

c_compl <- cophenetic(ph_compl)
cor(d, as.dist(c_compl))

c_avg <- cophenetic(ph_avg)
cor(d, as.dist(c_avg))

c_w2 <- cophenetic(ph_w2)
cor(d, as.dist(c_w2))


##
library(cluster)
complete3 <- cutree(tree = hc_avg, k = 3)
plot(silhouette(x = complete3, dist = d), cex.names = 0.6)


##
library(pvclust)
# итераций должно быть 1000 и больше, здесь мало для скорости
cl_boot <- pvclust(t(st_w), method.hclust = "average", nboot = 100,
                   method.dist = "euclidean", parallel = TRUE, iseed = 42)


##
plot(cl_boot)
pvrect(cl_boot) # достоверные ветвления


##
seplot(cl_boot)
print(cl_boot) # все значения


##
set.seed(395)
untang_w <- untangle_step_rotate_2side(den_compl, den_w2, print_times = F)

# танглграмма
tanglegram(untang_w[[1]], untang_w[[2]],
           highlight_distinct_edges = FALSE,
           common_subtrees_color_lines = F,
           main = "Tanglegram",
           main_left = "Left tree",
           main_right = "Right tree",
           columns_width = c(8, 1, 8),
           margin_top = 3.2, margin_bottom = 2.5,
           margin_inner = 4, margin_outer = 0.5,
           lwd = 1.2, edge.lwd = 1.2,
           lab.cex = 1.5, cex_main = 2)


# Раскраска дендрограмм ====================================

# Раскраска происходит в порядке следования веток на дендрограмме

# а) Вручную
# Здесь в примере просто произвольные цвета
cols <- rainbow(30)
den_single_manual <- color_labels(dend = den_single, col = cols)
plot(den_single_manual, horiz = TRUE)

# б) При помощи функции
# Функция для превращения лейблов в цвета
# (группы определяются по `n_chars` первых букв)
get_colours <- function(dend, n_chars, palette = "Dark2"){
  labs <- get_leaves_attr(dend, "label")
  group <- substr(labs, start = 0, stop = n_chars)
  group <- factor(group)
  cols <- brewer.pal(length(levels(group)), name = palette)[group]
  return(cols)
}

# Применяем функцию
library(RColorBrewer)
cols <- get_colours(dend = den_single, n_chars = 3)
den_single_c <- color_labels(dend = den_single, col = cols)
plot(den_single_c, horiz = TRUE)


##
webpage <- "http://evolution.genetics.washington.edu/book/primates.dna"
primates.dna <- read.dna(webpage)
d_pri <- dist.dna(primates.dna)
hc_pri <- hclust(d_pri, method = "average")
ph_pri <- as.phylo(hc_pri)
plot.phylo(ph_pri)
axisPhylo()


##
plot.phylo(ph_pri, type = "radial")

##### Самостоятельная работа ####
