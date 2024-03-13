# ---
# title: "Анализ главных компонент"
# subtitle: "Анализ и визуализация многомерных данных с использованием R"
# author: Марина Варфоломеева, Вадим Хайтов, Анастасия Лянгузова
# ---

# ## Тренировочные данные: медузы из реки Хоксбери ####
# Данные из Lunn & McNeil, 1991

# Исходные данные
jelly <- read.delim("data/jellyfish.csv")
X_raw <- jelly[, 2:3]
X <- scale() # Центрируем
X_mat <- as.matrix()
A <- t() %*%/(nrow()) # Матрица ковариаций

E <-            # Спектральное разложение
U <-           # Собственные векторы
Lambda <-     # Собственные числа
# Координаты точек в новом пространстве
Y <- X %*% U

# график исходных данных

gg_jelly_raw <- ggplot(as.data.frame(X_raw),
                       aes(x = width, y = length)) +
  geom_point(size = 2)

# график в главных осях
(gg_jelly_rotated <- ggplot(as.data.frame(Y),
                            aes(x = V1, y = V2)) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    labs(x = "PC1", y = "PC2"))

# графики вместе


# график
gg_rotated +
  aes(colour = jelly$location) +
  scale_color_brewer("Location", palette = "Set1") +
  labs(title = "Результаты PCA",
       x = paste0("PC1, ", round(Explained[1] * 100, 1), "%"),
       y = paste0("PC2, ", round(Explained[2] * 100, 1), "%"))

# ## Пример: Потребление белков в странах Европы с разными видами продуктов питания ####
# Данные из Weber, 1973

# ## Открываем данные
protein <- read.table(file="data/protein.csv", sep="\t", dec=".", header=TRUE)
protein$region <- factor(protein$region)
rownames(protein) <- protein$country
head(protein)

# ## Делаем PCA  ####
library(vegan)
prot_pca <- rda(protein[, -c(1, 2)], scale = TRUE)
op <- par(mar = c(5, 4, 0, 1) + 0.1)
biplot(prot_pca)
par(op)


# ## Разбираемся с результатами PCA  ####
tmp <- summary(prot_pca)

# # 1. Сколько компонент нужно оставить? ----

eigenvals(prot_pca) # собственные числа

bstick(prot_pca) # ожидаемое по Broken Stick Model

screeplot(prot_pca, type = "lines", bstick = TRUE) # График собственных чисел

prop_expl <- eigenvals(prot_pca)/sum(eigenvals(prot_pca))
cumsum(prop_expl*100)

# # 2. Графики факторных н агрузок и ординации ----
# ## Параметр `scaling`
# Внимание! Координаты объектов или переменных можно получить в нескольких вариантах, отличающихся масштабом. От этого масштаба будет зависеть интерпретация.

# Графики в vegan
op <- par(mfrow = c(1, 2))
# График факторных координат
biplot(prot_pca, display = "sites", scaling = "sites")
# График факторных нагрузок
biplot(prot_pca, display = "species", scaling = "species")
par(op)

# ## Те же самые графики можно построить в ggplot2
library(ggplot2)
# Данные для графиков
df_scores <- data.frame(protein[, 1:2],
                        scores(prot_pca, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
## График ординации в ggplot

p_scores <- ggplot(df_scores, aes(x = PC1, y = PC2, colour = region)) +
  geom_text(aes(label = country)) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))

## Данные для графика нагрузок
df_load <- as.data.frame(scores(prot_pca, display = "species", choices = c(1, 2, 3), scaling = "species"))
# поправки для размещения подписей
df_load$hjust <- ifelse(df_load$PC1 >= 0, -0.1, 1)
df_load$vjust <- ifelse(df_load$PC2 >= 0, -0.1, 1)
library(grid) # для стрелочек
ar <- arrow(length = unit(0.25, "cm"))
## График нагрузок в ggplot
p_load <- ggplot(df_load) +
  geom_text(aes(x = PC1, y = PC2, label = rownames(df_load)),
            size = 3, vjust = df_load$vjust, hjust = df_load$hjust) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               colour = "grey40", arrow = ar) +
  coord_equal(xlim = c(-2, 2), ylim = c(-2, 2))



# Все вместе
library(cowplot)
plot_grid(p_load, p_scores, align = "hv", rel_widths = c(0.4, 0.64))


# # 3. Интерпретация компонент -----

# Факторные нагрузки оценивают вклады переменных в изменчивость по главной компоненте
scores(prot_pca, display = "species", choices = c(1, 2, 3), scaling = 0)




# # Создание составных переменных при помощи PCA  ####


# ## При помощи дисперсионного анализа можно проверить, различается ли значение первой главной компоненты ("Мясо -- злаки и орехи") между разными регионами Европы

# Значения факторов (= факторные координаты)
df <- data.frame(region = protein$region,
  scores(prot_pca, display = "sites", choices = c(1, 2, 3), scaling = "sites"))
mod <- lm(PC1 ~ region, data = df)
anova(mod)

# ## Проверка условий применимости дисперсионного анализа
mod_diag <- fortify(mod)
res_p <- ggplot(data = mod_diag, aes(x = .fitted, y = .stdresid)) + geom_point(aes(size = .cooksd)) + geom_hline(yintercept = 0) + geom_smooth(method="loess", se=FALSE)
mean_val <- mean(mod_diag$.stdresid)
sd_val <- sd(mod_diag$.stdresid)
norm_p <- ggplot(mod_diag, aes(sample = .stdresid)) + geom_point(stat = "qq") + geom_abline(intercept = mean_val, slope = sd_val)
plot_grid(res_p, norm_p, ncol = 2, rel_widths = c(0.55, 0.45))


# ## График значений первой компоненты по регионам
df$region <- reorder(df$region, df$PC1, FUN=mean)
ggplot(df, aes(x = region, y = PC1, colour = region)) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 1) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

# ## Пост-хок тест
TukeyHSD(aov(mod))

# # Классический подход к морфометрии ####

# ## Пример: морфометрия черепах ####

# Черепахи - единственные живые представители анапсид (череп не имеет височных окон). Морфология черепа важна для их систематики (Claude et al., 2004).
# Данные - 24 разных измерения черепов черепах 122 ныне живущих пресноводных, морских и наземных видов и одного ископаемого. (Из Zuur et al. 2007)

# ## Читаем данные
turt <- read.table("data/turtles.txt", header = TRUE)
turt$Environment3 <- factor(turt$Environment3, levels = c(0, 1, 2, 9), labels = c("Freshwater", "Terrestrial", "Marine", "Fossil"))
colnames(turt)


# Нужно ли стандартизовать исходные данные?
boxplot(x = turt[8:31])

# ## Делаем анализ главных компонент
library(vegan)
turt_pca <- rda(turt[, 8:31], scale = TRUE)

# ## Сколько компонент достаточно для описания данных?
eig <- eigenvals(turt_pca)[1:5]
eig*100/sum(eig) # доля объясненной изменчивости
screeplot(turt_pca, bstick = TRUE)

# ## Что странного в этой картинке?
biplot(turt_pca, display = "species", scaling = 2)

# - Как вы думаете, почему у всех переменных большие нагрузки по первой компоненте?
# - Что отражает первая компонента?


# ## Двойное центрирование - один из классических способов избавиться от влияния размера

# Функция, которая может центрировать вектор
center <- function(x){
  x - mean(x, na.rm = TRUE)
}
# Почему для двойного центрирования перед PCA достаточно применить эту функцию к строкам?
dbcent <- t(apply(turt[, 8:31], 1, center))

# PCA
turt_db_pca <- rda(dbcent, scale = TRUE)
eig_db <- eigenvals(turt_db_pca)[1:5]
eig_db*100/sum(eig_db)
screeplot(turt_db_pca, bstick = TRUE)
biplot(turt_db_pca, display = "species", scaling = 2)


# ## Код для графика ординации черепах по морфометрии черепов
op <- par(mfrow = c(1, 2), mar = c(4, 4, 0.5, 0.5), cex = 1)
biplot(turt_db_pca, display = "species", scaling = 2)
# цвета для графика факторных координат
colvec <- c("orange2", "limegreen", "steelblue", "red3")
# пустой график
plot(turt_db_pca, type = "n", scaling = 1)
# точки, раскрашенные по уровням фактора turt$Environment3
points(turt_db_pca, display = "sites", scaling = 1, pch = 21,
       col = colvec[turt$Environment3], bg = colvec[turt$Environment3])
# легенда
legend("topright", legend = levels(turt$Environment3), bty = "n", pch = 21,
       col = colvec, pt.bg = colvec)
par(op)


# # Геометрическая морфометрия  ####

# ## Пример: Форма головы Апалачских саламандр рода _Plethodon_ ####
#
# _Plethodon jordani_ и _P.teyahalee_ встречаются вместе и раздельно.
# В совместно обитающих популяциях меняется форма головы обоих видов. В разных группах популяций этот процесс параллельно приводит к одинаковым результатам. По-видимому, одной из причин параллельной эволюции может быть межвидовая конкуренция (Adams, 2004, 2010).

# install.packages("geomorph", dependencies = TRUE)
library(geomorph)
data(plethodon)
str(plethodon, vec.len = 2, give.attr = F)

# ## Сырые морфометрические данные еще не выравнены
# Два образца для примера
plotRefToTarget(plethodon$land[, , 1], plethodon$land[, ,10],
                method = "points", mag = 1,
                links = plethodon$links)
# Слева - три образца, справа - все. Жирные точки - центроиды соответствующих меток
op <- par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plotAllSpecimens(plethodon$land[, , 1:3], links=plethodon$links)
plotAllSpecimens(plethodon$land,links=plethodon$links)
par(op)


## Шаг 1. Выравниваем данные при помощи обобщенного прокрустова анализа ----

gpa <- gpagen(plethodon$land, print.progress = FALSE)
plotAllSpecimens(gpa$coords,links=plethodon$links)

# ## Усредненная форма
ref <- mshape(gpa$coords)
plotRefToTarget(ref, ref, method = "TPS", links = plethodon$links)

# ## Можем посмотреть, как отличается любой из образцов от усредненной формы

# Изменение формы можно представить графически несколькими способами


# матрица, в которой хранится разметка общего графика
m <- matrix(data = c(1, 2,
                     3, 3),
            nrow = 2, ncol = 2, byrow = TRUE)
l <- layout(m, heights = c(1, 1), widths = c(1, 1))
# layout.show(l) # можно просмотреть разметку

# Графики
op <- par( mar = c(0, 0, 0, 0))
# 1) изменение конфигурации обозначено векторами
plotRefToTarget(ref, gpa$coords[, , 11],
                method = "vector", mag = 1,
                links = plethodon$links)
# 2) формы обозначены точками
plotRefToTarget(ref, gpa$coords[, , 11],
                method = "points", mag = 1,
                links = plethodon$links)
# 3) сплайн
plotRefToTarget(ref, gpa$coords[, , 11],
                method = "TPS", mag = 1,
                links = plethodon$links)
par(op)



# ## Шаг 2. Создаем морфопространство ----
op <- par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
ord <- gm.prcomp(gpa$coords)
plot(ord, main = "PCA")


# ## Можно раскрасить по группам
gp <- as.factor(paste(plethodon$species, plethodon$site)) # группа должна быть фактором
# задаем соответствие цветов уровням фактора
colvec <- c("Jord Allo" = "yellow2",
            "Jord Symp" = "orange",
            "Teyah Allo" = "green4",
            "Teyah Symp" = "green1")
# вектор цветов в порядке заданном фактором gp
colvec <- colvec[match(gp, names(colvec))]
# график
plot(ord, bg = colvec, pch = 21, col = "grey20")
# легенда
legend("topright", legend = levels(gp),
       bty = "n", pch = 21,
       col = "grey20",
       pt.bg = levels(as.factor(colvec)))
par(op)


# ## Доля объясненной изменчивости и факторные координаты
expl <- round(ord$d[1:5]/sum(ord$d) * 100, 1) # Доля изменчивости объясненной 1-5 компонентами
head(ord$x[, 1:5]) # Факторные координаты по 1-5 компонентам

# ## Чтобы легко рисовать изменения формы вдоль главной компоненты нам понадобится функция
plot_shape_change <- function(ord, ref_shape, PC,
                              horiz = TRUE,
                              gridPars = NULL, ...){
  if(horiz){
    op <- par(mfrow = c(1, 2), mar = c(0, 0 , 0, 0))
    plotRefToTarget(M1 = ref_shape, M2 = ord$shapes[[PC]]$min,
                    gridPars = gridPars,  ...)
    plotRefToTarget(M1 = ref_shape, M2 = ord$shapes[[PC]]$max,
                    gridPars = gridPars, ...)
    par(op)
  } else {
    op <- par(mfrow = c(2, 1), mar = c(0, 0 , 0, 0))
    plotRefToTarget(M1 = ref_shape, M2 = ord$shapes[[PC]]$max,
                    gridPars = gridPars,  ...)
    plotRefToTarget(M1 = ref_shape, M2 = ord$shapes[[PC]]$min,
                    gridPars = gridPars, ...)
    par(op)
  }
}

# ## Изменение формы вдоль главных компонент относительно средней формы
plot_shape_change(ord, ref_shape = gpa$consensus, PC = 1, links = plethodon$links, method = "TPS")

plot_shape_change(ord, ref_shape = gpa$consensus, PC = 2, links = plethodon$links, method = "TPS", horiz = FALSE)

# ## Можно нарисовать одновременно изменение формы вдоль обеих компонент и ординацию
library(cowplot)
library(gridGraphics)
my_gridPar <- gridPar(tar.pt.size = 0.6, grid.lwd = 0.7)

gg_pca <- plot_grid(
  ~ plot_shape_change(ord, ref_shape = gpa$consensus, PC = 2,
                      horiz = FALSE, links = plethodon$links,
                      method = "TPS", gridPars = my_gridPar),
  ~ {plot(ord, bg = colvec, pch = 21, col = "grey20")
    legend("topright", legend = levels(gp),  bty = "n",
           pch = 21, col = "grey20",
           pt.bg = levels(as.factor(colvec)))},
  NULL,
  ~ plot_shape_change(ord, ref_shape = gpa$consensus, PC = 1,
                      links = plethodon$links,
                      method = "TPS", gridPars = my_gridPar),
  ncol = 2, rel_heights = c(5, 1), rel_widths = c(1, 4)
)

gg_pca


# # Эволюционные изменения формы #####

# ## Фило-морфопространство саламандр рода Plethodon
data(plethspecies)
str(plethspecies, vec.len = 2, give.attr = F)

# ## Выравниваем средние формы для видов
species_gpa <- gpagen(plethspecies$land) #GPA-alignment

# ## Наложение филогенетического дерева и анцестральных форм на график PCA ординации

# Филоморфопространство
pca_with_phylo <- gm.prcomp(species_gpa$coords, phy = plethspecies$phy)
plot(pca_with_phylo, phylo = TRUE)

