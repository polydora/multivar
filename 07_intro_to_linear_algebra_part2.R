#' title: "Краткое введение в мир линейной алгебры. Часть 3"
#' subtitle: "Анализ и визуализация многомерных данных с использованием R"
#' author: Вадим Хайтов, Марина Варфоломеева

####################################################################
# Level 8: Собственные значения, собственные векторы и главные оси #
####################################################################
library(ggplot2)


#' Во многих многомерных методах требуется найти оси максимального варьирования

set.seed(123456789)

x1 <- rnorm(500, 30, 4)
y1 <- rnorm(500, 700, 50)
x2 <- rnorm(500, 40, 5)
y2 <- 10 * x2 + 200 + rnorm(500, 0, 100)

XY <-data.frame(x = c(x1, x2), y = c(y1, y2) )

ggplot(XY, aes(x, y)) +
  geom_point() +
  labs(x = "Переменная 1", y = "Переменная 2") +
  geom_point(aes(x = mean(x), y = mean(y)), size = 4, shape = 21, fill = "yellow")


XY_cent <- as.data.frame(scale(XY,
                               center = T, scale = F))

head(XY_cent)


# Угол между векорами в многомерном пространстве

Cos_alpha <-
  with(XY_cent,
       (x %*% y) /
         (norm(t(x), type = "F") *
            norm(t(y), type = "F"))  )

Cos_alpha

cor(XY$x, XY$y)



#' ## Нормализуем векторы

x_norm <- XY$x/sqrt(sum(XY$x)^2)
y_norm <- XY$y/sqrt(sum(XY$y)^2)


XY_norm <- data.frame(x = x_norm, y = y_norm)


ggplot(XY_norm , aes(x = x, y = y)) + geom_point() +
  geom_point(aes(x = mean(x), y = mean(y)), size = 4, color = "yellow")


#' ## Центрируем нормализованные векторы

XY_norm_cent <- as.data.frame(scale(XY_norm,  center = TRUE, scale = FALSE))

ggplot(XY_norm_cent , aes(x = x, y = y)) + geom_point() +
  geom_point(aes(x = mean(x), y = mean(y)), size = 4, color = "yellow")


#' ## Находим ковариационную матрицу {.smaller}

mXY_norm_cent <- as.matrix(XY_norm_cent)

Sxy_norm_cent <- t(mXY_norm_cent) %*% mXY_norm_cent/(nrow(mXY_norm_cent) - 1)

Sxy_norm_cent

#' ## Находим собственные числа и собственные векторы {.smaller}

eig <- eigen(Sxy_norm_cent) # Стандартная функция R для извлечения собственных чисел и собственных векторов

Lambda <- eig$values # Собственные числа


# Превращаем вектор собственный чисел в матрицу
diag(Lambda)


U <- eig$vectors # Собственные векторы

# Должна получиться ковариационная матрица
U %*% diag(Lambda)  %*% solve(U)

# Проверим
Sxy_norm_cent

U[,1] %*% U[,2]





#' ## Стандартизованные собственные векторы {.smaller}

U_scaled <- U %*% sqrt(diag(Lambda)) #


(U %*% sqrt(diag(Lambda))) %*% t(U %*% sqrt(diag(Lambda)))




#' ## Рисуем собственные векторы {.smaller}

PC1 <- data.frame(x = c(mean(XY_norm_cent$x), U_scaled[1, 1]),
                  y = c(mean(XY_norm_cent$y),  U_scaled[2,1]))

PC2 <- data.frame(x = c(mean(XY_norm_cent$x),  U_scaled[1, 2]),
                  y = c(mean(XY_norm_cent$y),  U_scaled[2,2]))

ggplot(XY_norm_cent, aes(x = x, y = y)) + geom_point() +
  geom_point(aes(x = mean(x), y = mean(y)), size = 4, color = "yellow") +
  geom_line(data = PC1, aes(x = x, y = y), color = "yellow", size = 1.5)  +
  geom_line(data = PC2, aes(x = x, y = y), color = "yellow", size = 1.5) +
  coord_equal()



#' ## Рисуем главные оси {.smaller .columns-2}

ggplot(XY_norm_cent, aes(x = x, y = y)) + geom_point() +
  geom_point(aes(x = mean(x), y = mean(y)), size = 4, color = "yellow") +
  geom_line(data = PC1, aes(x = x, y = y), color = "yellow", size = 1.5)  +
  geom_line(data = PC2, aes(x = x, y = y), color = "yellow", size = 1.5) +
  coord_equal() + geom_abline(slope = tan(acos(U[1,1])), color = "blue") +
  geom_abline(slope = (U[2,2])/(U[1,2]), color = "blue")

#' ## Вращение осей {.smaller .columns-2}

#' Вращающая матрица

angle <- -1 * acos(U[1,1]) #Отрицательный угол, так как поворачиваем оси по часовой стрелке

Rot <- matrix(c(cos(angle), sin(angle),
                -sin(angle), cos(angle)), nrow = 2)
Rot


XY_norm_cent_rot <- as.data.frame(t(Rot %*% t(mXY_norm_cent)))

ggplot(XY_norm_cent, aes(x = x, y = y)) +
  geom_point(color = "gray") +
  geom_point(data = XY_norm_cent_rot, aes(x = V1, y = V2)) +
  labs(x = "Первая главная ось", y = "Вторая главная ось") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


## Задание
# Исследуйте структуру матрицы *X*.


set.seed(123456789)

x1 <- c(rnorm(250, 30, 4), rnorm(250, 60, 4))
x2 <- rnorm(500, 70, 50)
x3 <- rnorm(500, 40, 5)
x4 <- c(rnorm(100, 10, 5), rnorm(100, 40, 5), rnorm(100, 70, 5), rnorm(200, 100, 5))
x5 <- c(rnorm(250, 50, 5), rnorm(250, 100, 5))


X <-data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)



# Постройте методами матричной алгебры ковариационную матрицу
# Проведите ее спектральное разложение (вычислите ее собственные числа и собственные векторы).
# Оцените информативность главных осей.
# Изобразите точки в пространстве первой и второй главной оси
X_cent <- scale(X, center = T, scale = F)

Cov_X <- (t(X_cent) %*% X_cent)/(nrow(X_cent) - 1)

eig <- eigen(Cov_X)

Lambda <- eig$values

U <- eig$vectors

sum(as.vector(Lambda/(sum(Lambda)))[1:2])

X_projected <- as.matrix(X) %*% U[,1:2]

qplot(X_projected[,1], X_projected[,2] )


#########################################################################
# Level 9: Сингулярное разложение матриц (Singular value decomposition) #
#########################################################################



# Работа с изображениями, как с матричными объектами #
######################################################



################## Поворот изображения ##########################3


load("data/face.rda")

faceData
dim(faceData)

library(reshape2)

faceData_XY <- melt(faceData) ## Переводим матрицу в два вектора координат и вектор значений интенсивности заливки

names(faceData_XY) <- c("X1", "X2", "value")


ggplot(faceData_XY, aes(X1, X2)) + geom_tile(aes(fill = value)) + scale_fill_gradient(low = "darkblue",   high =  "white" ) + coord_equal()



# Повернем изображение на угол 30 градусов


angle <-  -30*pi/180 #Задаем угол поворота в радианах

# Вращающая матрица
Rot <- matrix(c(cos(angle), sin(angle),
                -sin(angle), cos(angle)), nrow = 2)

Image_rot <-   data.frame(t((Rot) %*% t(faceData_XY[, 1:2] )), value = faceData_XY[3]) #Надо заполнить пропуски

ggplot(Image_rot, aes(X1, X2)) + geom_point(aes(color = value), size = 5) + scale_fill_gradient(low = "darkblue",   high =  "white" )


# Проведем масштабирование полученного изображения

Scale <- matrix(c(3, 0, 0, 1), nrow = 2)

Image_trans <-   data.frame(t((Scale) %*% t(Image_rot[,1:2])), value = faceData_XY$value)

ggplot(Image_trans, aes(X1, X2)) + geom_point(aes(color = value), size = 5) + scale_fill_gradient(low = "darkblue",   high =  "white" ) + coord_equal()



# Задание: посмотрите, что за изображение спрятано в следующей матрице


img <- read.csv("data/a_matrix.csv")

img <-as.matrix(img)

img_XY <- melt(img) ## Переводим матрицу в два вектора координат и вектор значений интенсивности заливки

names(img_XY) <- c("X1", "X2", "value")


# Изобразите эту матрицу






# Сингулярное разложение матриц #
#################################

set.seed(123456789)
B <- matrix(round(runif(50, 1, 5))  , byrow = T, ncol=5) #Некоторая матрица
SVD <- svd(B) #Сингулярное Разложение матрицы B с помощью функции svd()
V <- SVD$v #"Вспомогательная" матрица - правые сингулярные векторы
D <- SVD$d #Вектор сингулярных чисел
U <- SVD$u #"Вспомогательная" матрица - левые сингулярные векторы

# Восстановим исходную матрицу

B_reconstructed <- U %*% diag(D) %*% t(V)



# Задание
# Вычислите матрицу, которая получится при использовании только 1 и 2 сингулярного числа для матрицы B, использованной на предыдущем слайде.

B_reconstructed <- U[,1:5] %*% diag(D[1:5]) %*% t(V[,1:5])

qplot(B, B_reconstructed) +
  geom_abline()









#Применение сингулярного разложения матриц  в сжатии изображений #
##################################################################


load("data/face.rda")

gg_face <- function(x) {
  library(reshape2)
  library(ggplot2)
  rotate <- function(x) t(apply(x, 2, rev))
  dd <- rotate(x)
  ddd <- melt(dd)
  ggplot(ddd, aes(Var1, Var2)) + geom_tile(aes(fill = value)) + scale_fill_gradient(low = "darkblue",   high =  "white" ) + coord_equal()
}

gg_face(faceData)


SVD_face <- svd(faceData)

U_face <- SVD_face$u
D_face <- SVD_face$d
V_face <- SVD_face$v




#' ##Рекоструируем изображение, используя только часть информации

reduction <- function(x, U, D, V) U[,1:x] %*% diag(D[1:x]) %*% t(V[, 1:x])

gg_face(reduction(20, U_face, D_face, V_face))


library(jpeg)

ander1 <- readJPEG("images/Anderson.jpg")[,,1]

str(ander)

gg_face(ander[, , 3])


### Самостоятельная работа ############################


# Задание про волков


# Восстановление изображения


# Сжатие изображения




