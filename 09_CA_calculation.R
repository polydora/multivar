protein <- read.table(file="data/protein.csv", sep="\t", dec=".", header=TRUE)

protein$region <- factor(protein$region)
rownames(protein) <- protein$country
head(protein)

library(vegan)

prot_pca <- rda(protein[, -c(1,2)])

plot(prot_pca)

eigenvals(prot_pca)
summary(prot_pca)

scores(prot_pca, choices = 1:ncol(protein[, -c(1,2)]))

prot_cent <- scale(protein[, -c(1,2)], center = TRUE, scale = FALSE)

D <- svd(prot_cent)$d

diag(D)

U <- svd(prot_cent)$u

V <- svd(prot_cent)$v

library(ggplot2)

qplot(y = D^2, x = eigenvals(prot_pca)) + geom_line()

dim(U)


plot(x = U, y= prot_pca$CA$u)

plot(x = V, y= prot_pca$CA$v)




########vegan#######################################
# МЕХАНИКА КОРРЕСПОНДЕНТНОГО  АНАЛИЗА #
#######################################

#Вадим Хайтов, Марина Варфоломеева

# Проблемы PCA

library(readxl)
birds <- read_excel(path = "data/macnally.xlsx")

# имена переводим в нижний регистр
colnames(birds) <- tolower(colnames(birds))

# Проведите анализ главных компонент и визуализируйте его результаты для данных приведенных в датасете macnally.xlsx

#Код для анализа главных компонент

library(vegan)

bird_pca <- rda(birds[,-c(1:2)], scale = TRUE)
screeplot(bird_pca, bstick = T)

eigenvals(bird_pca)


# Код для вывода информации об информативности PC

summary(bird_pca)


plot(bird_pca, display = "sites")

# Код для построения биплота

biplot(bird_pca, scaling = "sites")

biplot(bird_pca, scaling = "species")

biplot(bird_pca, scaling = "species", display = "species")







# Механика Корреспондентного анализа

library(vegan)
data(mite)
data(mite.env)
data(mite.xy)
head(mite[ , 1:6], 2)
str(mite)


str(mite.xy)
str(mite.env)


mite_pca <- rda(mite, scaling = TRUE)

screeplot(mite_pca, bstick = T)

biplot(mite_pca,  scaling = "sites", type = "t")

mite$LCIL


plot(mite_pca,  display = "sites", type = "t")


mite_mds <- metaMDS(mite)

plot(mite_mds, display = "site")

summary(mite_pca)

# CA в темную ######


mite_ca <- cca(mite)

summary(mite_ca)


biplot(mite_ca)

biplot(mite_pca)


screeplot(mite_ca, bstick = T)

plot(mite_ca, display = "sites")

plot(mite_ca, display = "species")




#### Матрицы сопряженности #####

peas <- matrix(c(99, 42, 29, 13), byrow = T, ncol = 2)



Ft <- sum(peas)



f_i <- apply(peas, MARGIN = 1, FUN = sum)

f_j <- apply(peas, 2, FUN = sum)



p_i <- f_i / Ft #Вектор вероятностей для формы
p_j <- f_j / Ft #Вектор вероятностей для цвета




q <- p_i %*% t(p_j)




E <- q * Ft
E <- round(E, 1)


O <- peas


sum((O-E)^2/E)

chisq.test(x = O, p = q, correct = F)



### CA Вручную ##############


Ft <- sum(mite)

f_ij <- mite #Частота встречи данного вида в данной пробе, то есть это первичные даные!

p_ij <- mite/Ft #вероятность встречи данного вида в данной пробе


Ft <- sum(mite) #Общее количество найденных животных

f_i <- apply(mite, MARGIN = 1, FUN = sum) #Общее количество особей в каждой пробе

p_i <- f_i/Ft #Вектор вероятностей встретить какую-либо особь в данной пробе

f_j <- apply(mite, MARGIN = 2, FUN = sum) #Общее количество особей в каждом виде

p_j <- f_j/Ft #Вектор вероятностей встретить особь данного вида


q <- p_i %*% t(p_j) #вероятность встретить особь в данной пробе.


E <- (p_i %*% t(p_j) * Ft)

O <- mite

Chi2 <- sum((O-E)^2/E)

summary(mite_ca)

Inertia <- Chi2/Ft




mite_cca <- cca(mite)
summary(mite_cca)


Q1 <- (p_ij - p_i %*% t(p_j))/sqrt(p_i %*% t(p_j))

#Та же матрица, вычисленная через частоты
Q <- (f_ij*Ft - f_i %*% t(f_j))/(Ft*sqrt(f_i %*% t(f_j)))

round(sum(Q1 - Q))



sum(Q^2)




### Сингулярное разложеине матрицы сопряженности ######

U <- svd(Q)$u
D <- svd(Q)$d
V <- svd(Q)$v

dim(U)

dim(V)

round(D, 2)

str(D)


Qsvd <- U %*% diag(D) %*% t(V) #матрица "восстановленная" из "вспомогательных" матриц

round(sum(Q - Qsvd)) #разность между исходной и "восстановленной" матрицами


# Связь SVD и собственных значений

D <- diag(D)


dim(D)

round(t(Q) %*% Q -   V %*% t(D) %*% D %*% t(V))

Q <- as.matrix(Q)

A <- t(Q) %*% Q


eig_values <- eigen(A)$values #Собственные числа матрицы A
eig_vectors <- eigen(A)$vectors #Матрица собственных векторов для матрицы A


plot(eig_values, diag(D))

round(eig_values, 4)

eigen(t(Q) %*% Q)$values #Собственные значения для матрицы ковариации Q'Q

svd((t(Q) %*% Q))$d #Это те же собственные значения.

diag(D)^2 #Квадраты сингулярных чисел


plot(eigen(t(Q) %*% Q)$values, diag(D))

sum(eig_values)


Information <- data.frame(
  CA = 1:length(eig_values),
  Eigenval =round(eig_values, 5),
  Prop_Explained = round(eig_values/sum(eig_values), 5),
  Cumul_Prop=round(cumsum(eig_values/sum(eig_values)),5)
)


summary(mite_ca)


CA_samples <- diag(p_i^(-1/2))%*% U[,1:2]


library(ggplot2)
Pl_CA_st <-
  ggplot(as.data.frame(CA_samples), aes(x=V1, y=V2) ) +
  geom_text(label = rownames(mite)) +
  geom_hline(yintercept=0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(x= "CA1", y = "CA2")


Pl_CA_st


CA <- as.data.frame(CA_samples)

mite.env

plot(x = mite.env$SubsDens, y = CA$V2)




CA_species <- diag(p_j^(-1/2))%*% V[,1:2]




Pl_CA_sp <-
  ggplot(as.data.frame( CA_species), aes(x = V1, y = V2) )  +
  geom_hline(yintercept=0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(x= "CA1", y = "CA2") +
  geom_text(label = names(mite))

Pl_CA_sp






# Задание: Проведите корреспондентный анализ данных по птицам Австралии, используя "ручной" метод обработки
# После этого сравните полученные результаты с теми результатами, которые получаются с использованием функции из пакета vegan.

library(ggvegan)

birds_ca <- cca(birds[, -c(1,2)])

autoplot(birds_ca, scaling = "sites") +
  theme_bw()


birds_2 <- birds[, -c(1,2)] #только численность

Ft <- sum(birds_2) #Общее количество
f_ij <- birds_2 #Частота встречи данного вида в данной пробе, то есть это первичные даные!
p_ij <- birds_2/Ft #вероятность встречи данного вида в данной пробе


f_i <- apply(birds_2, MARGIN = 1, FUN = sum) #Общее количество особей в каждой пробе
p_i <- f_i/Ft #Вектор вероятностей встретить какую-либо особь в данной пробе

f_j <- apply(birds_2, MARGIN = 2, FUN = sum) #Общее количество особей в каждом виде
p_j <- f_j/Ft #Вектор вероятностей встретить особь данного вида

dim(birds_2)



Q <- (p_ij - p_i %*% t(p_j))/sqrt(p_i %*% t(p_j))

sum(Q^2)


U <- svd(Q)$u
D <- svd(Q)$d
V <- svd(Q)$v

# Q <- as.matrix(Q)
# A <- t(Q) %*% Q #матрица ковариации
#
# eig_values <- eigen(A)$values #Собственные числа матрицы A
# eig_vectors <- eigen(A)$vectors #Матрица собственных векторов для матрицы A


CA_samples <- diag(p_i^(-1/2))%*% U[,1:2]

CA_samples <- as.data.frame(CA_samples)

ggplot(CA_samples, aes(V1, V2)) +
  geom_point()

CA_species <- diag(p_j^(-1/2))%*% V[,1:2]

CA_species <- as.data.frame(CA_species)




