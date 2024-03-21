# ---
# title: "Анализ избыточности (Redundancy analysis, RDA)"
# subtitle: "Анализ и визуализация многомерных данных с использованием R"
# author: Марина Варфоломеева, Вадим Хайтов, Анастасия Лянгузова


# ### Пример: генетика бабочек Euphydryas editha ####
#
# Частоты разных аллелей фосфоглюкоизомеразы и данные о факторах среды для 16 колоний бабочек Euphydryas editha в Калифорнии и Орегоне (данные McKechnie et al., 1975)

# # Подготовка данных к RDA ####

library(ade4)
data(butterfly)
str(butterfly)

# ## Создадим переменные с более короткими названиями для удобства
gen <- butterfly$genet
head(gen, 1)
env_geo <- cbind(butterfly$envir, butterfly$xy)
head(env_geo, 1)

# ## Нужно ли стандартизовать зависимые переменные?
summary(gen)

# ## Линейны ли связи между переменными?
pairs(data.frame(gen, env_geo))


# ## Удаляем колинеарные предикторы

library(car)
vif(lm(gen$`0.4` ~ ., data = env_geo))
vif(lm(gen$`0.4` ~ . -Temp_Min, data = env_geo))

# ## RDA в vegan ####

library(vegan)
bf_rda <- rda(gen ~ Altitude + Precipitation + Temp_Max, data = env_geo)
summary(bf_rda)

# Обратите внимание на
# Структуру общей изменчивости
# Важность различных компонент
# Как распределяется изменчивость между осями
# Распределение изменчивости, потенциально объяснимой факторами

# ## Корреляции между откликами и предикторами ####
spenvcor(bf_rda)

# # Визуализация ординации #######################

# ## Триплот корреляций (scaling = 2): Какие переменные среды сильнее всего определяют сходство объектов? {.columns-2}
plot(bf_rda, scaling = 2)


# ## Триплот расстояний (scaling = 1): Насколько похожи друг на друга объекты? {.columns-2}
plot(bf_rda, scaling = 1)


# # Проверка значимости ординации ################

# ## Общий тест: Влияют ли факторы на зависимые переменные?
anova(bf_rda)


# ## Тест факторов, type I эффекты: Какие факторы влияют на зависимые переменные?
anova(bf_rda, by = "term")

# ## Тест факторов, type III эффекты: Какие факторы влияют на зависимые переменные?
anova(bf_rda, by = "mar")

# ## Тест значимости осей, ограниченных факторами: Вдоль какой из осей значимо меняется генетическая структура?
anova(bf_rda, by = "axis")


# # Выбор оптимальной модели ####################

m1 <- rda(gen ~ Altitude + Precipitation + Temp_Max, data = env_geo)
m0 <- rda(gen ~ 1, data = env_geo)

m <- ordistep(m0, scope = formula(m1), permutations = 99999)

m$anova

m_radj <- ordiR2step(m0, scope = formula(m1), permutations = 9999)

m_radj$anova


# # Частный анализ избыточности ##################


###### pRDA в виде функции
pRDA <- function(Y, X = NULL, W = NULL, scale.Y = FALSE)
{
  Y <- scale(as.matrix(Y), center = TRUE, scale = scale.Y)
  if (!is.null(W)) {
    # При наличии ковариат W
    W <- scale(as.matrix(W), center = TRUE, scale = FALSE)
    Y <- qr.resid(qr(W), Y)
  }
  if (!is.null(X)) {
    # При начилии матрицы X
    X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
    X <- cbind(X, W)
    Q <- qr(X)
    RDA <- svd(qr.fitted(Q, Y))
    RDA$w <- Y %*% RDA$v %*% diag(1/RDA$d)
    Y <- qr.resid(Q, Y)
  } else {
    # Если нет ни X, ни W
    RDA <- NULL
  }
  RES <- svd(Y)
  # Анализ главныз компонент по остаткам
  list(RDA = RDA, RES = RES)
}

######################################


# ## Делаем частный RDA: зависимость генетической структуры от среды с учетом географического положения
bf_prda_1 <- rda(gen ~ Altitude + Condition(x + y), data = env_geo)
anova(bf_prda_1, permutations = 99999) ## Пермутационный тест

plot(bf_prda_1, main = "Partial RDA", scaling = )


# # Компоненты объясненной инерции ###############
showvarparts(2)

# ## Подбираем модели RDA, нужные для поиска компонентов инерции
# bf_prda_1 уже есть
bf_prda_2 <- rda(gen ~ x + y + Condition(Altitude), data = env_geo)
bf_rda_full <- rda(gen ~ x + y + Altitude, data = env_geo)

# ## Задание: Найдите компоненты инерции -----------------------
# По результатам трех RDA найдите всю потенциально объяснимую инерцию, а так же долю инерции, объясненную:
# - средой и географией в сумме (a + b + c)
# - средой, но не с географией (a)
# - географией, но не со средой (b)
# - совместно средой и географией (c) ???
# - а если выразить все (a, b, c) в процентах?
# В процентах от чего???


# ## Данные для самостоятельной работы -----------------------
# Влияние выпаса скота на зообентос рек о.Тасмания. Данные из работы Magierowski, RH, Davies, PE, Read, SM (2015). Biological and physical data collected from river sites selected across a gradient of catchment area under grazing in northern Tasmania., Version 1. http://doi.org/10.4227/05/54BEEE248D22B. ÆKOS Data Portal, rights owned by University of Tasmania. Accessed 02 Mar 2017.
# В таблице содержится несколько листов:
# fauna - обилия разных групп животных на разных реках;
# env - параметры среды в каждой из рек:
##### Abstraction - накопленый индекс водозабора (сумма всех водозаборов и водоотводов вверх по течению, делённая на средний многолетний (25-30 лет) годовой сток),
##### Regulation - накопленный индекс регулирования (сумма водосборных водохранилищ вверх по течению, делённая на средний многолетний годовой сток),
##### Grazing  (доля от общей площади) - выпас скота,
##### fines (proportion substrata) - мелкие отложения,
##### Temperature (oC) - температура,
##### Conductivity (uS/cm) - проводимость,
##### average turbidity (NTU) - средняя мутность,
##### pH,
##### Alkalinity Total (mg CaCO3/L) - общая щелочность,
##### Nitrate+Nitrite (mg-N/L),
##### DRP (mg-P/L) - растворённый реактивный фосфор,
##### N total (mg-N/L),
##### P Total (mg-P/L),
##### Average % shading - средняя затенённость,
##### Average algae cover (%) - средний покров водорослей,
##### Chl a (mg/m2),
##### GrazingRank - уровень выпаса скота;
# coord - географические координаты каждой из рек.


library(readxl)
fauna <- read_xls('data/Grazing_Magierowski_et_al_2015.xls', sheet = 'fauna')
env <- read_xls('data/Grazing_Magierowski_et_al_2015.xls', sheet = 'env')
coord <- read_xls('data/Grazing_Magierowski_et_al_2015.xls', sheet = 'coord')



# Морфометрия поссумов в Австралии. Данные из работы Lindenmayer, D. B., Viggers, K. L., Cunningham, R. B., and Donnelly, C. F. 1995. Morphological variation among columns of the mountain brushtail possum, Trichosurus caninus Ogilby (Phalangeridae: Marsupiala). Australian Journal of Zoology 43: 449-458.
#
library(DAAG)
data(possum)
data(possumsites)
