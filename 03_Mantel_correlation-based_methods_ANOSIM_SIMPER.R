# Анализ связи между наборами данных
# Вадим Хайтов
# Многомерные методы на R

# Пакеты
library(vegan)
library(ggplot2)


# Загружаем данные
data(varespec)
data(varechem)

head(varespec)
str(varespec)


head(varechem)
str(varechem)


sum(is.na(varechem))
sum(is.na(varespec))


# Код для построения ординаци в осях nMDS

set.seed(12345)
ord <- metaMDS(comm = varespec, distance = "bray", k = 2, trymax = 40)

plot(ord)

stressplot(ord)

stres_ord <- ord$stress

ord$points

scores(ord)

mds_points <- as.data.frame(ord$points)

# код для построения графика ординации









# Применяем функцию envfit()
env_fit <- envfit(  ~ ., data = )
env_fit

# Визуализация результатов
ordiplot( , display = )
plot( )


# Анализ связи с переменными c помощью функции `ordisurf()`

ordiplot( , display = "sites")

ordisurf( , varechem$Al,
         add = TRUE, col="blue", method = "REML")
ordisurf(  ,varechem$Mn,
         add = TRUE, col="green")



# Задание: Отразите связь ординации растительности со значениями концентрации гумуса.




# Вычисление мантеловской корреляции
str(varespec)
str(varechem)


dist_com <- vegdist(varespec, method = "bray")

dist_chem <- vegdist(varechem, method = "euclidean")


x <- as.vector(dist_com)
y <- as.vector(dist_chem)

qplot(x, y) + geom_smooth(se = F, method = "lm")


xy <- data.frame(x, y)

R <- round(cor(x, y, method = "spearman"), 3)

mant <- ggplot(xy, aes(x = x, y = y))

mant + geom_point(size=3) + xlab("Biological dissimilarity") +
  ylab("Chemical dissimilarity") +
  annotate("text", x = 0.25, y = 0.35, label=paste("Rspearmen =", R, sep=" ")) +
  theme_bw() + geom_smooth(method = "lm", se = FALSE)




cor.test(x, y, method = "pearson") # Это неправильное действие! Так делать нельзя!








## Пермутационный метод

set.seed(12345)

male <- rnorm(100, 130, 5)
female <- rnorm(100, 129, 5)

t.test(male, female)



SE_m <- sd(male) / sqrt(length(male))
SE_f <- sd(female) / sqrt(length(female))

t_initial <- (mean(male) - mean(female))/sqrt(SE_m^2 + SE_f^2)
t_initial

f <- female
m <- male

num_perm <- sample(1:100,1)
order_m <- sample(1:100, num_perm)
order_f <- sample(1:100, num_perm)

f[order_f] <- male[order_f]
m[order_m] <- female[order_f]

SE_m <- sd(m) / sqrt(length(m))
SE_f <- sd(f) / sqrt(length(f))

t_p <- (mean(m) - mean(f))/sqrt(SE_m^2 + SE_f^2)
t_p


Nperm <- 10000

tperm <- rep(NA, Nperm)

set.seed(12345)

for (i in 1:(Nperm - 1))
{
  BOX <- c(male,female)
  ord <-sample(1:200, 200)
  f <- BOX[ord[1:100]]
  m <- BOX [ord[101:200]]
  SE_m <- sd(m) / sqrt(length(m))
  SE_f <- sd(f) / sqrt(length(f))
  tperm[i]=(mean(m) - mean(f))/sqrt(SE_m^2 + SE_f^2)
}


head(tperm)
tail(tperm)

tperm[Nperm] <- t_initial

tdf <- data.frame(t = tperm)

ggplot(tdf, aes(x =t)) + geom_histogram(fill="blue", color = "black") +
  geom_vline(xintercept = c(t_initial, -t_initial))


p_value <- mean(tperm <= -t_initial |tperm >= t_initial )
p_value


# Пермутационная оценка значимости корреляции
library(coin)
library(MASS)

set.seed(1234567)

mu <- c(10, 20) #Вектор средних значений

Sigma <- matrix(.7, nrow=2, ncol=2) #Ковариационная матрица

diag(Sigma) <- c(1, 3)

dat <- as.data.frame(mvrnorm(n=100, mu=mu, Sigma=Sigma))

qplot(dat$V1, dat$V2)

cor.test(dat$V1, dat$V2, method = "spearman")

spearman_test(V1 ~ V2, data = dat,
              distribution = approximate(nresample = 99999))



# Пермутационная оценка значимости мантеловской корреляции

mant <- mantel(dist_com, dist_chem, method="pearson", permutations = 9999)
mant


# Частная мантеловская корреляция

# Матрица координат описаний
geo <- read.table("data/Coordinates.txt",header = TRUE, sep = "\t")

# Матрица расстояний между точками
dist_geo <- vegdist(geo[, -1], method = "euclidean")


mantel_partial <- mantel.partial(dist_com, dist_chem,
                                 dist_geo, method = "pearson",
                                 permutations = 9999)
mantel_partial

# Функция `bioenv()`из пакета `vegan`


# Количество комбинаций

2^ncol(varechem) - 1


BioEnv <- bioenv(varespec, varechem, method = "spearman", index = "bray")

BioEnv


# Код для вмзуализации связи ординации растительности с парамтерами среды, которые вошли в финальную модель BioEnv






# Оценка статистической значимости результатов BIO-ENV
# ЗАПУСКАТЬ КОД МЕЖДУ ДВУМЯ ЛИНИЯМИ ТОЛЬКО ЕСЛИ НЕ ЖАЛКО ВРЕМЕНИ!!!
#------------------------------------
perm_binv <- c(1:100)
for (i in 1:99)
{
  perm_num <- sample(1:nrow(varespec))
  perm_com <- varespec[perm_num,]
  perm_i <- bioenv(perm_com, varechem)

  manteltop <- length(perm_i$models)
  est <- c(1:manteltop)
  for (j in 1:ncol(varechem)) est[j] <- perm_i$models[[j]]$est
  perm_binv[i] <- max(est)
  cat("interation", i, "\n")
}

manteltop <- length(BioEnv$models)


estim <- c(1:manteltop)
for (j in 1:ncol(varechem)) estim[j] <- BioEnv$models[[j]]$est
perm_binv[100] <- max(estim)

perm_binv <- data.frame(perm_i = perm_binv)

p = length(perm_binv[perm_binv[,1] >= perm_binv[100,1],1]) /
  nrow(perm_binv)

# png("BIOENV.png")

hist <- ggplot(perm_binv, aes(x=perm_i))
hist + geom_histogram (bin=0.1, fill="blue", colour="black") +
  geom_vline(xintercept=perm_binv[100,1], linetype=2) + theme_bw() +
  xlab("Мантеловские корреляции, \nполученные при пермутациях \nпроцедуры BIO-ENV ") +
  annotate("text", x=0.8, y=20, label=(paste("P=", p, sep=" ")))

# dev.off()
#------------------------------------



# Модельные матрицы

com <- read.csv("data/mussel_beds.csv",
                sep=';', header = T)

ascam <- read.csv("data/ASCAM.csv",
                  sep=';', header = T)

library(dplyr)

log_com <- com %>%
  filter(Bank == "Vor2") %>%
  select(-c(1:3)) %>% decostand(., method = "log")


log_ascam <- ascam %>%
  filter(Bank == "Vor2") %>%
  select(-c(1:2)) %>% decostand(.,method = "log")


mds_vor2_com <- as.data.frame(metaMDS(vor2_log_com)$points)


mds_vor2_ascam <- as.data.frame(metaMDS(log_ascam,
                                        distance = "euclid")$points)

ggplot(mds_vor2_ascam, aes(x = MDS1, y = MDS2)) +
  geom_path() +
  geom_text(label = 1997:2011)



dist_com <- vegdist(log_com, method = "bray")

dist_ascam <- vegdist(log_ascam, method = "euclidean")


mantel(dist_com, dist_ascam)



## Задание
# 1. Выясните есть ли многолетний градиент в динамике размерной струтуры и структуры
#сообщества на банке Vor4.
# 2. Оцените связь между размерной структурой мидий и структурой сообщества.



## Градиентная модельная матрица
gradient_model <- vegdist(com$Year[com$Bank == "Vor2"], method="euclidian")
gradient_model

## Тестируем гипотезу о наличии градиента с помощью теста Мантела

# Получаем матарицы расстояний для двух наборов данных

dist_vor2_com <- vegdist(log_com, method = "bray")
dist_vor2_ascam <- vegdist(log_ascam, method = "euclidean")



### 1) Наличие градиента в структуре сообщества
mantel(dist_com, gradient_model)



### 2) Наличие градиента в размерной структуре мидий
mantel(dist_ascam, gradient_model)

## Прослеживается ли связь между размерной структурой мидий и структурой сообщества?

### Не самое правильное решение
mantel(dist_vor2_com, dist_vor2_ascam)

### Более корректное решение
mantel.partial(dist_ascam , dist_com , gradient_model )


## Циклическая модельная матрица

cycmod <- function(x){
  points <- data.frame(X=c(1:x), Y=c(1:x))
  for (i in 1:x) {
    points$X[i] <- cos(2*pi*(i-1)/x)
    points$Y[i] <- sin(2*pi*(i-1)/x)
  }
  return(points)
}

qplot(cycmod(nrow(mds_com))$X, cycmod(nrow(mds_com))$Y, xlab="X", ylab="Y", geom = "point", size = 4)

cycl_model <- round(vegdist(cycmod(nrow(mds_vor2_ascam)), method = "euclidean"))
cycl_model

## Выявляется ли циклическая составляющая в динамике размерной структуры?
mantel(dist_ascam, cycl_model)

## Более корректная оценка
mantel.partial(dist_ascam, cycl_model, gradient_model)

str(varechem)


# ANOSIM, SIMPER

com <- read.csv("data/mussel_beds.csv", sep=';', header = T)
ascam <- read.csv("data/ASCAM.csv", sep=';', header = T)
com$Mussel_size <- factor(com$Mussel_size)


str(com)



# ANOSIM: Analysis Of Similarity

## Задание
# - Постройте ординацию всех описаний датасета `com` (логарифмированные данные) в осях nMDS на основе матрицы Брея-Куртиса
# - Раскрасьте разными цветами точки, относящиеся к двум разным группам: "Large-dominated" и "Small-dominated"


library(vegan)
library(ggplot2)

log_com <- log(com[,-c(1:3)] + 1)

ord_log_com <- metaMDS(log_com, distance = "bray", k=2)

# Обратите внимание, что при разных способах извлечения данных из объеута 0ординации будут разные заголовки столбцов

MDS <- data.frame(scores(ord_log_com)[[1]])


ggplot(MDS, aes(x = NMDS1, y = NMDS2, fill = com$Mussel_size)) +
  geom_point(shape = 21, size = 4) +
  scale_fill_manual(values = c("red", "blue")) +
  ggtitle(paste("Stress = ", round(ord_log_com$stress, 2))) +
  theme_bw()





# **Задание:**
#
# # 1. Вычислите матрицу коэффициентов Брея-Куртиса на основе матрицы `log_com`
# # 2. Разверните полученную матрицу в вектор.
# # 3. На основе полученного вектора создайте вектор, содержащий ранги расстояний.
#
dist_com <- vegdist(log_com, method = "bray")
#
# write.table(as.data.frame(dist_com), "clipboard", sep = "\t", row.names = F)

unfold_dist_com <- as.vector(dist_com)


rank_dist_com <-




# **Задание 4. Создайте треугольную матрицу `dummy_dist`, той же размерности, что и матрица `dist_com`, в которой `0` будет с стоять в тех ячейках, которые соответствуют межгрупповым расстояниями, а   `1` -- внутригрупповым.

dummy_dist <- dist()

dummy_dist <- ifelse(dummy_dist == 0, 0, 1)




# **Задание:**
#
# 5. Вычислите средние значения рангов внутри- и межгрупповых расстояний
# 6. Вычислите R-статистику

dists <- data.frame(rank = rank_dist_com, dummy = as.vector(dummy_dist))

library(dplyr)


mean_dists <- dists %>%
  group_by() %>% summarize(rank_mean = )

n <- nrow(log_com)

R_glob <-





# **Задание:**

# 7. Напишите пользовательскую функцию (пусть она будет называться `R_perm`), которая пермутировала бы принадлежность каждого объекта к той или иной группе и вычисляла бы значение R-статистики для новой комбинации.
# 8. Используя функцию `for()` вычислите 10000 значений  R-статистики и запишите их в вектор.


R_perm <- function(comm, group){
  require(vegan)
  dist_com <- vegdist(comm)
  rank_dist_com <- rank(dist_com)
  dummy_dist <- dist() #Придумайте как перемешать группы
  dummy_dist <- ifelse(dummy_dist == 0, 0, 1)
  dists <- data.frame(rank = rank_dist_com, dummy = as.vector(dummy_dist))
  require(dplyr)
  mean_dists <- dists %>% group_by(dummy) %>% summarize(rank_type = mean(rank))
  n <- nrow(log_com)
  R_p <- (mean_dists$rank_type[2] - mean_dists$rank_type[1])/(n * (n - 1)/4)
  R_p
}


R_perm(comm = log_com,
       group = com$Mussel_size)

R_perms <- rep(NA, 10000)

for(i in 1:9999) R_perms[i] <- R_perm(comm = log_com, group = com$Mussel_size)

R_perms[10000] <-



# **Задание:**
#
# 9. Постройте частотную гистограмму, характеризующую распределение пермутационных оценок.
# 10. Нанесите на нее полученное значение $R_{global}$.
# 11. Вычислите уровень значимости.

R_perms <- data.frame()

Pl_our <- ggplot(R_perms, aes(x = R_perms)) +
  geom_histogram() + geom_vline(xintercept = ) +
  xlim(-0.2, 0.2)

Pl_our

#P-value




## Процедура ANOSIM в пакете `vegan`
com_anosim <- anosim(log_com,
                     grouping = com$Mussel_size,
                     permutations = 9999,
                     distance = "bray")

## Задание
# Изучите структуру объекта `com_anosim` и постройте частотное распределение
#значений $R_{global}$, полученных при каждом акте пермутации

R_perms_vegan <- data.frame(vegan_R = )

Pl_vegan <- ggplot(R_perms_vegan, aes(x = vegan_R)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "blue") +
  geom_vline(xintercept = , linetype = 2) +
  xlim(-0.2, 0.2)

Pl_vegan

library(gridExtra)

grid.arrange(Pl_our, Pl_vegan)



## Ограничения (Assumptions) для применения ANOSIM

# Внутригрупповые расстояния (ранги)
plot(com_anosim, main = "Dissimilarity ranks \n between and within classes")


## Задание

# + Постройте ординацию в осях nMDS, раскрасив точки в разные цвета в зависимости от номера мидиевой банки
# + Проверьте гипотезу о различиях в структуре сообществ на разных банках
# + Проверьте условия применимости ANOSIM
# + Проведите попарное сравнение всех банок

ggplot( , aes(x = MDS1, y = MDS2, fill = )) +
  geom_point(shape = 21, size = 4) +
  scale_fill_manual(values = c("red", "blue", "green")) +
  labs(fill = "Mussel beds") +
  ggtitle(paste("Stress = ", round(ord_log_com$stress, 3), sep = " "))


# Проверка гипотезы о различиях в структуре сообществ на разных банках
bank_anosim <- anosim(log_com, grouping = )

# Применимость
plot(bank_anosim)

# Попарное сравнение банок
# Vor2 vs Vor4
anosim()

# Vor2 vs Vor5
anosim()

#Vor4 vs Vor5
anosim()


# SIMPER: Similarity Percentages

## Какие признаки зависимой матрицы вносят наибольший вклад в формирование различий между группами?
log_com_simper <- simper(log_com, group = com$Mussel_size, permutations = 999)
summary(log_com_simper)


## Задание
# Выявитие виды, отвечающие за различия в сообществах разых банок

log_com_simper2 <- simper(log_com, group = com$Bank, permutations = 9999)

summary(log_com_simper2)




### Самостоятельная работа

#### Тест Мантела
#### Растительные сообщества во Франции La Mafragh (Annaba, Algérie)

#### ANOSIM, SIMPER
#### Данные по мидиям [mussel_patches.csv](data/mussel_patches.csv)
