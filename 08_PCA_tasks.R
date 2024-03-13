# Данные из Machine Learning Repository
# https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)

# Результаты тонкоигольной аспирационной пункционной биопсии. Описание ядер клеток

# Переменные:
# 1) id - идентификационный номер пациента
# 2) diagnosis - диагноз (M = malignant, B = benigh)
# 3-12) среднее, 13-22) - стандартное отклонение и 23-32) - худшее значение следующих признаков:
# a) radius (mean of distances from center to points on the perimeter; среднее расстояние от центра до точек на периметре)
# b) texture (standard deviation of gray-scale values; стандартное отклонение значений серого)
# c) perimeter (периметр)
# d) area (площадь)
# e) smoothness (local variation in radius lengths; локальное изменение длины радиуса)
# f) compactness (perimeter^2 / area - 1.0)
# g) concavity (severity of concave portions of the contour; выраженность вогнутых участков по контору клетки)
# h) concave points (number of concave portions of the contour; количество вогнутых участков по контору клетки)
# i) symmetry
# j) fractal dimension ("coastline approximation" - 1)
#

# # Задание:
# - Проведите анализ главных компонент. Какую долю общей изменчивости объясняют первые две главные компоненты?
# - Постройте график ординации объектов в пространстве первых двух компонент и раскрасьте точки в зависимости от диагноза.
# - При помощи таблицы или графика факторных нагрузок определите, какие признаки вносят вклад в изменчивость вдоль первых двух главных компонент.
# - Вдоль какой из компонент более выражено разделение облаков точек?

# Открываем данные и создаем названия признаков
brc <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data", header = F)
features <- c("radius", "texture", "perimeter", "area", "smoothness", "compactness", "concavity", "concave_points", "symmetry", "fractal_dimension")
names(brc) <- c("id", "diagnosis", paste0(features,"_mean"), paste0(features,"_se"), paste0(features,"_worst"))

# Геометрическая морфометрия тела рыб Cyprindon pecosensis.
# - coords - координаты лендмарок, выравненные при помощи gpa
# - CS - размер центроида
# - Sex - пол ("F","M")
# - Pop - популяция ("Marsh","Sinkhole")
# Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described by high-dimensional data. Heredity. 115: 357-365.

# Задание:
# Сделайте PCA по выравненным координатам лендмарок.
# Сколько изменчивости объясняют первые две главных компоненты?
# Нарисуйте график главных компонент:
# - раскрасьте точки в зависимости от пола и популяции рыб,
# - приведите графики изменения формы вдоль главных компонент.

library(geomorph)
data("pupfish")
