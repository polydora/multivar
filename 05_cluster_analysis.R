# install.packages(c("candisc", "ape", "dendextend", "fpc", "pvclust"))

## Пример: Волки ===========================================
# Морфометрия черепов у волков в Скалистых горах и в Арктике (Jolicoeur, 1959)

# Загружаем данные #
library(candisc)
data("Wolves")


##### Анализ данных ####
dim(Wolves)
colnames(Wolves)
head(rownames(Wolves))
any(is.na(Wolves))
table(Wolves$group)

## Задание 1 -------------------------------------------------
# - Постройте ординацию nMDS данных о морфометрии волков
# - Оцените качество ординации
# - Обоснуйте выбор коэффициента
# - Раскрасьте точки на ординации волков в зависимости от
# географического происхождения (`group`)

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


# Кластерный анализ в R =====================================

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

plot(ph_single, type = "phylogram")
plot(ph_single, type = "cladogram")
plot(ph_single, type = "fan")
plot(ph_single, type = "unrooted")
plot(ph_single, type = "radial")

axisPhylo()


# С помощью dendextend
den_single <- as.dendrogram(hc_single)
plot(den_single, horiz = TRUE)


##### Метод отдалённого соседа ####

# Визаулизация средства ape
hc_compl <- hclust(d, method = "complete")
ph_compl <- as.phylo(hc_compl)
plot(ph_compl, type = "phylogram")
axisPhylo()


# Визуализация средствами dendextend
den_compl <- as.dendrogram(hc_compl)
plot(den_compl, horiz = TRUE)


##### Метод невзвешенного попарного среднего (UPGMA) ####

# Средствами ape
hc_avg <- hclust(d, method = "average")
ph_avg <- as.phylo(hc_avg)
plot(ph_avg, type = "phylogram")
axisPhylo()


# Средствами dendextend
den_avg <- as.dendrogram(hc_avg)
plot(den_avg, horiz = TRUE)


##### Метод Варда ####

# Средствами ape
hc_w2 <- hclust(d, method = "ward.D2")
ph_w2 <- as.phylo(hc_w2)
plot(ph_w2, type = "phylogram")
axisPhylo()

# Средствами dendextend
den_w2 <- as.dendrogram(hc_w2)
plot(den_w2, horiz = TRUE)


## Кофенетическая корреляция ===============================
# Матрица кофенетических расстояний
c_single <- cophenetic(ph_single)

# Кофенетическая корреляция =
# = корреляция матриц кофенетич. и реальн. расстояний
cor(d, as.dist(c_single))

## Задание 2 ------------------------------------------------
# Оцените при помощи кофенетической
# корреляции качество кластеризаций, полученных разными
# методами.
# Какой метод дает лучший результат?



# Качество и количество кластеров ==========================

## Ширина силуэта ==========================================

# Оценка ширины силуэта для 3 кластеров
library(cluster)
avg3 <- cutree(tree = hc_avg, k = 3) # делим дерево на нужное количество кластеров
plot(silhouette(x = avg3, dist = d), cex.names = x.names = 0.6)


## Бутстреп поддержка ветвей ===============================
library(pvclust)
set.seed(389)
# итераций должно быть 1000 и больше, здесь мало для скорости
cl_boot <- pvclust(t(st_w), method.hclust = "average", nboot = 100,
                   method.dist = "euclidean", parallel = TRUE, iseed = 42)


plot(cl_boot)
pvrect(cl_boot) # достоверные ветвления


## Для диагностики качества оценок AU
seplot(cl_boot)
print(cl_boot) # все значения


##### Танглграмма ####
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

## Задание 3 ===============================================
# Постройте танглграмму из дендрограмм,
# полученных методом ближайшего соседа и методом Варда.



# Раскраска дендрограмм ====================================

# Раскраска происходит в порядке следования веток на дендрограмме

# а) Вручную
# Здесь в примере просто произвольные цвета
cols <- rainbow(30)
den_avg_manual <- color_labels(dend = den_avg, col = cols)
plot(den_avg_manual, horiz = TRUE)

# б) С помощью функции
# Функция для превращения лейблов в цвета
# (группы определяются по `n_chars` первых букв в лейбле)
get_colours <- function(dend, n_chars, palette = "Dark2"){
  labs <- get_leaves_attr(dend, "label")
  group <- substr(labs, start = 0, stop = n_chars)
  group <- factor(group)
  cols <- brewer.pal(length(levels(group)), name = palette)[group]
  return(cols)
}

cols <- get_colours(dend = den_avg, n_chars = 3)
den_avg_c <- color_labels(dend = den_avg, col = cols)
plot(den_avg_c, horiz = TRUE)

##### Метод k-средних ####

set.seed(333)
w_kmeans <- kmeans(st_w, centers = 3, nstart = 100)

# Сравниваем полученные кластеры с результатом UPGMA
table(avg3, w_kmeans$cluster)

# Выбор нужного количества кластеров

w_cascade <- cascadeKM(st_w, inf.gr = 2, sup.gr = 10,
                       iter = 100, criterion = 'calinski')
plot(w_cascade, sortg = TRUE) # чтобы объекты, относящиеся к одному кластеру, рисовались вместе

# Визуализация
library(factoextra)
fviz_cluster(w_kmeans, data = st_w,
             ggtheme = theme_bw()) # 3 кластера

w_2k <- kmeans(st_w, centers = 2,
               nstart = 100) # 2 кластера

fviz_cluster(w_2k, data = st_w,
             ggtheme = theme_bw())

## Оценка качества кластеризации по ширине силуэта

plot(silhouette(w_kmeans$cluster, d)) # 3 кластера
plot(silhouette(w_2k$cluster, d)) # 2 кластера

##### DBSCAN: Density-based spatial clustering of applications with noise ####
library(dplyr)
data('multishapes')

multi_circle <- multishapes %>% filter(shape == c(1, 2))
multi <- multi_circle[, 1:2]

ggplot(multi, aes(x, y)) + geom_point()

# когда метод K-средних не справляется...
set.seed(123)
circle_kmeans <- kmeans(multi, centers = 2, nstart = 20)
my_col_circle <- c("#2E9FDF", "#E7B800")
fviz_cluster(circle_kmeans, data = multi,
             palette = my_col_circle)

## DBSCAN
library(dbscan)
circle_dbscan <- dbscan(multi, 0.23, 5)
fviz_cluster(circle_dbscan, data = multi,
             palette = my_col_circle)

## Выбор нужного радиуса
kNNdistplot(multi, k = 5)
abline(h = 0.23, lty = 2)

circle_dbscan <- dbscan(multi, 0.23, 5)
fviz_cluster(circle_dbscan, data = multi,
             palette = my_col_circle)

## Задание 4 ===============================================
# Проанализируйте данные по волкам, используя DBSCAM алгоритм



##### Метод нечёткой кластеризации C-средних (C-means, fuzzy clustering) ####

w_cmeans <- fanny(d, k = 4, memb.exp = 2)
summary(w_cmeans)

## Оценка качества кластеризации
plot(silhouette(w_cmeans), cex.names = 0.6)

## Задание 5 ===============================================
# Попробуйте оценить ширину силуэта для C-means кластеризации с более подходящим числом кластеров.

## Визуализация кластеризации методом C-средних ####
library(scatterpie)
w_clust <- cbind(dfr_w, $membership)
ggplot() + geom_scatterpie(data = w_clust, aes(x = MDS1, y = MDS2),
                           cols = c())

##### Самостоятельная работа ####
# Проанализируйте данные об относительных обилиях фораминифер в пробах на Белом
# море (данные из Golikova et al. 2020 из файла Golikova_etal_2020_cluster_data.csv).
#
# - Выберите и обоснуйте трансформацию данных и расстояние.
#
# - Постройте ординацию nMDS по относительным обилиям фораминифер: - цвет
# значков --- растение-доминант, - форма значков --- точка сбора.
#
# - Постройте дендрограмму проб по сходству отн. обилий фораминифер. - оцените
# при помощи кофенетической корреляции, какой метод аггрегации лучше,
# - Постройте визуализацию для методов неиерархической кластеризации
#
# - Опишите получившиеся кластеры при помощи различных параметров: - ширина
# силуэта - бутстреп-поддержка ветвлений

foram <- read.table("data/Golikova_etal_2020_cluster_data.csv",
                    sep = "\t", header = TRUE)
