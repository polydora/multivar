# Анализ датасета

# - Зообентос рек Тасмании --- Grazing_Magierowski_et_al_2015.xls

##### Загружаем данные ####
library(readxl)
library(car)
library(vegan)
library(tidyverse)
library(ggplot2)
library(DAAG)

env <- as.data.frame(read_xls("data/Grazing_Magierowski_et_al_2015.xls", sheet = 'env'))
rownames(env) <- env$SITE
env <- env[, 2:18]
head(env)

coord <- read_xls("data/Grazing_Magierowski_et_al_2015.xls", sheet = 'coord')
head(coord)
coord <- coord[, 2:3]

fauna <- as.data.frame(read_xls("data/Grazing_Magierowski_et_al_2015.xls", sheet = 'fauna'))
rownames(fauna) <- fauna$SITE
fauna <- fauna[, 2:200]
head(fauna)

zoobentos <- na.omit(data.frame(env, coord, fauna))

fauna <- zoobentos[, 20:218]
env <- zoobentos[, 1:17]
env_geo <- zoobentos[, 1:19]
env_geo$GrazingRank <- as.factor(env_geo$GrazingRank)
head(env_geo)



########## Подготовка данных ######

##### Трансформация данных ####
View(fauna)
fauna <- fauna[, which(colSums(fauna) != 0)]

fauna_log <- log(fauna + 1)

##### Линейность связей и мультиколлинеарность ####
pairs(env_geo)

vif(lm(fauna_log$Nematoda ~ ., data = env_geo))
vif(lm(fauna_log$Nematoda ~ . -Regulation -Alkalinity.Total..mg.CaCO3.L.
       -P.Total..mg.P.L. -northing_GDA94 -easting_GDA94 -GrazingRank -Nitrate.Nitrite..mg.N.L.,
       data = env_geo))


##### RDA анализ ####

zoo_rda <- rda(fauna_log ~ Abstraction + Grazing...proportion.of.total.catchment.area.
               + fines..proportion.substrata. + Temperature..oC. + Conductivity..uS.cm.
               + average.turbidity..NTU. + pH + N.total..mg.N.L. + DRP..mg.P.L.
               + Average...shading + Average.algae.cover.... + Chl.a..mg.m2., data = env_geo)

summary(zoo_rda)

##### Оценка значимости модели ####
anova(zoo_rda, permutations = 9999)

##### Type I эффекты ####
anova(zoo_rda, by = "term", permutations = 9999)

##### Type III эффекты ####
anova(zoo_rda, by = "mar", permutations = 9999)

##### Значимость осей ####
anova(zoo_rda, by = "axis", permutations = 9999)


##### Подбор оптимальной модели ####

full_mod <- rda(fauna_log ~ Abstraction + Grazing...proportion.of.total.catchment.area.
                + fines..proportion.substrata. + Temperature..oC. + Conductivity..uS.cm.
                + average.turbidity..NTU. + pH + N.total..mg.N.L. + DRP..mg.P.L.
                + Average...shading + Average.algae.cover.... + Chl.a..mg.m2., data = env_geo)



null_model <- rda(fauna_log ~ 1, data = env_geo)

model <- ordistep(null_model, scope = formula(full_mod), permutations = 9999)
model$anova

m_radj <- ordiR2step(null_model, scope = formula(full_mod), permutations = 9999)
m_radj$anova

##### Проверка значимости оптимальной модели ####

optim_rda <- rda(fauna_log ~ Abstraction + Temperature..oC. + Chl.a..mg.m2.,
                 data = env_geo)

##### Оценка значимости модели ####
anova(optim_rda, permutations = 9999)

##### Type I эффекты ####

anova(optim_rda, by = "term", permutations = 9999)

##### Type III эффекты ####
anova(optim_rda, by = "mar", permutations = 9999)

##### Значимость осей ####
anova(optim_rda, by = "axis", permutations = 9999)


##### Изучаем summary от модели ####
sum_opt <- summary(optim_rda)

partit <- function(smr){
  if(!is.null(smr$partial.chi)){
    Inertia <- c(smr$tot.chi, smr$partial.chi, smr$constr.chi, smr$unconst.chi)
    nms <- c("Total", "Conditioned", "Constrained", "Unconstrained")
  } else {
    Inertia <- c(smr$tot.chi, smr$constr.chi, smr$unconst.chi)
    nms <- c("Total", "Constrained", "Unconstrained")
  }
  part <- data.frame(Inertia = round(Inertia, 2), Proportion = round(Inertia/smr$tot.chi, 3))
  rownames(part) <- nms
  cat("Partitioning of variance:\n")
  part
}

##### Структура общей изменчивости ####
partit(sum_opt)

# Влияние компонент
sum_opt$cont

# Распределение изменчивости, потенциально объяснимой факторами
sum_opt$concont

# Собственные векторы, нагрузки переменных = “species scores”
scores_opt_rda <- as.data.frame(scores(optim_rda, display = "species", choices = 1:2))
View(scores_opt_rda)
scores_opt_rda[which.max(scores_opt_rda$RDA1), ]
scores_opt_rda %>% slice_max(RDA2)

##### Визуализация результатов ####
op <- par(mar = c(4, 4, 0, 1))

## Триплот корреляций
plot(optim_rda, scaling = 2, display = c("sites", "bp"))


## Триплот расстояний
plot(optim_rda, scaling = 1)

plot_rda <- plot(optim_rda, scaling = 1)
splen <- abs(rowSums(plot_rda$species))
text(plot_rda, "sites", pch=21, col="black", cex=0.8)
text(plot_rda, "species", select = splen > 1.8, length = 0.05,
     col = "red")
par(op)

##### Частный RDA ####

## Влияние среды без учёта географии
zoo_prda_env <- rda(fauna_log ~ Abstraction + Temperature..oC. + Chl.a..mg.m2.
                    + Condition(northing_GDA94 + easting_GDA94), data = env_geo)

anova(zoo_prda_env, permutations = 9999)
sum_env <- summary(zoo_prda_env)

##### Триплоты ####
op <- par(mfrow = c(1, 2))

cor_prda <- plot(zoo_prda_env, main = "Correlation triplot", scaling = 2)
cor_splen <- abs(rowSums(cor_prda$species))
text(cor_prda, "sites", pch=21, col="black", cex=0.8)
text(cor_prda, "species", select = cor_splen > 0.8, length = 0.05,
     col = "red")

dist_prda <- plot(zoo_prda_env, main = "Distance triplot", scaling = 1)
dist_splen <- abs(rowSums(dist_prda$species))
text(dist_prda, "sites", pch=21, col="black", cex=0.8)
text(dist_prda, "species", select = dist_splen > 2, length = 0.05,
     col = "red")

par(op)


## Влияние географии без учёта среды
zoo_prda_geo <- rda(fauna_log ~ northing_GDA94 + easting_GDA94
                     + Condition(Abstraction + Temperature..oC.
                                 + Chl.a..mg.m2.), data = env_geo)
sum_geo <- summary(zoo_prda_geo)

## Совместное влияние среды и географии
zoo_prda_full <- rda(fauna_log ~ Abstraction + Temperature..oC. + Chl.a..mg.m2.
                   + northing_GDA94 + easting_GDA94, data = env_geo)
sum_full <- summary(zoo_prda_full)

##### Инерция, потенциально объяснённая...####

## потенциально объяснённая инерция
(I_total <- sum_full$constr.chi)

## средой, но не географией
(I_env <- sum_env$constr.chi)

## географией, но не средой
(I_geo <- sum_geo$constr.chi)

## средой и географией
(I_env_geo <- I_total - I_env - I_geo)

comp <- data.frame(Inertia = c(I_env, I_geo, I_env_geo, I_total))
rownames(comp) <- c('Только среда',
                    'Только география',
                    'Среда и география вместе',
                    'Общая объяснимая инерция')
comp$Proportion <- comp$Inertia/sum(comp$Inertia[1:3]) * 100
colnames(comp) <- c('Инерция', '%')
comp

# varpart(fauna_log, ~ northing_GDA94 + easting_GDA94,
#         ~ Abstraction + Temperature..oC. + Chl.a..mg.m2., data = env_geo)

##### CCA: канонический корреспондентный анализ на тех же данных ####

##### Подбор оптимальной модели для CCA

##### Подбор оптимальной модели ####

full_mod_cca <- cca(fauna_log ~ Abstraction + Grazing...proportion.of.total.catchment.area.
                    + fines..proportion.substrata. + Temperature..oC. + Conductivity..uS.cm.
                    + average.turbidity..NTU. + pH + N.total..mg.N.L. + DRP..mg.P.L.
                    + Average...shading + Average.algae.cover.... + Chl.a..mg.m2., data = env_geo)



null_model_cca <- cca(fauna_log ~ 1, data = env_geo)

m_radj_cca <- ordiR2step(null_model, scope = formula(full_mod),
                         direction = "forward", permutations = 9999)
m_radj_cca$anova

zoo_cca <- cca(fauna_log ~ Abstraction + Temperature..oC. + Chl.a..mg.m2., data = env_geo)

vif.cca(zoo_cca)

sum_cca <- summary(zoo_cca)
partit(sum_cca)


##### CCA вручную ####

## Вероятности
f_ij <- fauna_log #Частота встречи данного вида в данной пробе, то есть это первичные даные!

Ft <- sum(fauna_log) #Общее количество найденных животных

f_i <- apply(fauna_log, 1, FUN = sum) #Общее количество особей в каждой пробе

p_i <- f_i/Ft #Вектор вероятностей встретить какую-либо особь в данной пробе

f_j <- apply(fauna_log, 2, FUN = sum) #Общее количество особей в каждом виде

p_j <- f_j/Ft #Вектор вероятностей встретить особь данного вида

## Матрица вкладов
Q <- (f_ij*Ft - f_i %*% t(f_j))/(Ft*sqrt(f_i %*% t(f_j)))
Q <- as.matrix(Q)

## Сингулярное разложение матрицы вкладов (Q)
U <- svd(Q)$u
D <- diag(svd(Q)$d)
V <- svd(Q)$v

## Расчёт инерции
Inertia_total <- sum(D^2) #Общая инерция в системе
CA_number <- sum(round(D, 2) != 0) # Главных осей на 1 меньше, чем в матрице

## Модельная матрица (связь между предикторами - параметрами среды - и зависимыми переменными)
# видами)

X <- model.matrix(~ Abstraction + Grazing...proportion.of.total.catchment.area.
                  + fines..proportion.substrata. + Temperature..oC. + Conductivity..uS.cm.
                  + average.turbidity..NTU. + pH + N.total..mg.N.L. + DRP..mg.P.L.
                  + Average...shading + Average.algae.cover.... + Chl.a..mg.m2., data = env)

## Матрица коэффициентов
betas <- solve(t(X) %*% diag(p_i) %*% X) %*% (t(X) %*% diag(p_i)^(1/2) %*% Q)

## Матрица предсказанных значенй
Q_pred <- diag(p_i)^(1/2) %*% X %*% betas

## Сингулярное разложение матрицы предсказанных значений
U_pred <- svd(Q_pred)$u
D_pred <- diag(svd(Q_pred)$d) # канонические оси
V_pred <- svd(Q_pred)$v

## Инерция в пространстве канонических осей
Inertia_constrained <- sum(D_pred^2) #Инерция в ограниченной ординации

CCA_number <- sum(round(D_pred, 2) != 0) #Количество канонических осей в СCA

## Матрица остатков
Q_resid <- Q - Q_pred

## Сингулярное разложение матрицы остатков
U_res <- svd(Q_resid)$u
D_res <- diag(svd(Q_resid)$d)
V_res <- svd(Q_resid)$v

## Инерция в пространстве неканонических осей - остатков
Inertia_unconstrained <- sum(D_res^2)

##### Результат вручную сделанного CCA ####
cca_inertia <- data.frame(Inertia = c(Inertia_total, Inertia_constrained,
                                      Inertia_unconstrained))
rownames(cca_inertia) <- c("Total",
                    "Constrained",
                    "Unconstrained")
cca_inertia$Proportion <- cca_inertia$Inertia/Inertia_total
cca_inertia

partit(sum_cca)

##### Ординация проб в канонических осях CCA ####

## Данные, полученные вручную
constr_CA_samples <- diag(p_i^(-1/2))%*% U_pred

## Данные, полученные в результате cca()
cca_ord <- zoo_cca$CCA$u

qplot(constr_CA_samples[ ,1], cca_ord[ ,1]) + geom_abline(slope = 1)

ggplot(as.data.frame(constr_CA_samples), aes(x=V1, y=V2)) +
  geom_text(label = rownames(fauna_log)) + labs(x = "CCA1", y = "CCA2") +
  theme_bw() + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + ggtitle("Результаты, полученные вручную")

plot(zoo_cca, display = "lc", main = "Результаты сca()", scaling = "sites")

##### Ординация видов в канонических осях CCA ####

# Координаты видов канонических осях полученные вручную
constr_CA_species <- diag(p_j^(-1/2))%*% V_pred

# Координаты видов в канонических осях по версии cca()
cca_sp_constr <- zoo_cca$CCA$v

qplot(constr_CA_species[ ,1], cca_sp_constr[ ,1]) + geom_abline(slope = 1)

ggplot(as.data.frame(constr_CA_species), aes(x=V1, y=V2)) +
  geom_text(label = colnames(fauna_log)) + labs(x = "CCA1", y = "CCA2") +
  theme_bw() + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + ggtitle("Результаты, полученные вручную")

dev.off()
plot(zoo_cca, display = "sp", main = "Результаты сca()", scaling = "species")


##### Значимость ординации ####
anova(zoo_cca)

## Type I эффекты
anova(zoo_cca, by="term")

## Type III эффекты
anova(zoo_cca, by="mar")

## Значимость осей
anova(zoo_cca, by="axis")

##### Интерпретация результатов ####

## Влияние компонент
sum_cca$cont

## Распределение изменчивости, потенциально объяснимой факторами
sum_cca$concont

## Нагрузки species scores
# Влияние компонент
cca_zoo_scores <- scores(zoo_cca, display = "species", choices = 1:5)
View(cca_zoo_scores)

spenvcor(zoo_cca)

##### Визуализация ####
plot(zoo_cca, scaling = "sites",
     main = "scaling 1")

plot(zoo_cca, scaling = "species",
     main = "scaling 2")

## Триплот корреляций
op <- par(mar = c(4, 4, 0, 1))
plot(zoo_cca, scaling = 2, display = c("sites", "bp"))


## Триплот расстояний
plot_cca <- plot(zoo_cca, scaling = 1)
cca_splen <- abs(rowSums(plot_cca$species))
text(plot_cca, "sites", pch=21, col="black", cex=0.8)
text(plot_cca, "species", select = cca_splen > 3, length = 0.05,
     col = "red")
par(op)

##### Данные для самостоятельной работы ####

# 1. Фауна Долгой губы (обилие видов)

abund <- read.table("data/dolg_abundance.txt", skip = 1,
                    header = TRUE, sep = ";", row.names = 1)
env <- read.table("data/dolg_hydrology.txt", skip = 1,
                  header = TRUE, sep = ";", row.names = 1)

# 2. Деревья на острове Barro Colorado
bci_species <- read.table('data/BCI_species.csv',
                          header=TRUE, sep=',', row.names = 1)

bci_env <- read.table('data/BCI_env.csv',
                      header=TRUE, sep=',', row.names = 1)

# 3. Морфометрия поссумов в Австралии
data("possum")
data("possumsites")
possum$site <- as.factor(possum$site)
levels(possum$site) <- rownames(possumsites)
