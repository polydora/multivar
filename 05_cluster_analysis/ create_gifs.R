## Как работает метод ближайшего соседа

# Данные
dat <- data.frame(x = c(9, 19, 30, 32, 38, 50, 20, 50), y = c(42, 40, 10, 30, 60, 35, 42, 31))
rownames(dat) <- LETTERS[1:8]
plot(dat, type = "n"); text(dat, rownames(dat))
# Кластеризация
hc_s <- hclust(dist(dat), method = "single")
hc_c <- hclust(dist(dat), method = "complete")
hc_a <- hclust(dist(dat), method = "average")

cluster_ani <- function(dat, gg_dat, dist_fun = "vegdist", dist_method = "euclidean", hclust_method = "average", k = nrow(dat)){
  library(vegan)
  library(ggplot2)
  library(ggforce)
  library(dendextend)
  library(cowplot)
  if (dist_fun == "vegdist") {
    d <- vegdist(dat, method = dist_method)
  } else if (dist_fun == "dist") {
    d <- vegdist(dat, method = dist_method)
  } else {
    stop("dist_fun should be either `vegdist` or `dist`")
  }

  hc <- hclust(d, method = hclust_method)
  den <- as.dendrogram(hc)
  # ordination plot
  gg_ord <- ggplot(data = gg_dat, aes(x = MDS1, y = MDS2, label = rownames(gg_dat))) +
    coord_fixed() +
    geom_point(size = 2) +
    geom_text(hjust = 1.2, vjust = 1.2) +
    geom_mark_ellipse(aes(group = cutree(hc, k = k)), colour = "red", expand = 0.001) +
    scale_y_continuous(expand=c(0.1,0.1))

  # dendrogram plot
  gg_tree <- function(){
    par(mar = c(2, 2, 0, 0))
    if (k == 1) {
      plot(den)
    } else {
      plot(den)
      rect.dendrogram(den, k = k, lty = 1, lwd = 1, border = "red")
    }
  }
  # together
  plot_grid(gg_ord,
            gg_tree,
            nrow = 1, rel_widths = c(0.65, 0.35),
            hjust = 0, vjust = 1, scale = c(0.8, 0.9))
}

suppressWarnings(ord <- metaMDS(dat, distance = "euclidean", autotransform = FALSE))
gg_dat <- data.frame(ord$points)

gg_list_s <- lapply(8:1, function(x) cluster_ani(dat, gg_dat, hclust_method = "single", k = x))
gg_list_c <- lapply(8:1, function(x) cluster_ani(dat, gg_dat, hclust_method = "complete", k = x))
gg_list_a <- lapply(8:1, function(x) cluster_ani(dat, gg_dat, hclust_method = "average", k = x))
gg_list_w <- lapply(8:1, function(x) cluster_ani(dat, gg_dat, hclust_method = "ward.D2", k = x))


ggplot(data = gg_dat, aes(x = MDS1, y = MDS2, label = rownames(gg_dat))) +
  geom_point(size = 2) +
  geom_text(hjust = 1.2, vjust = 1.2) +
  scale_y_continuous(expand=c(0.1,0.1))

library(gifski)

gifski::save_gif(
  expr = sapply(gg_list_s, plot),
  gif_file = "cluster_analysis/animation_single.gif",
  width = 800,
  height = 600,
  delay = 2  # seconds between frames
)

gifski::save_gif(
  expr = sapply(gg_list_c, plot),
  gif_file = "cluster_analysis/animation_complete.gif",
  width = 800,
  height = 600,
  delay = 2  # seconds between frames
)

gifski::save_gif(
  expr = sapply(gg_list_a, plot),
  gif_file = "cluster_analysis/animation_average.gif",
  width = 800,
  height = 600,
  delay = 2  # seconds between frames
)

gifski::save_gif(
  expr = sapply(gg_list_w, plot),
  gif_file = "cluster_analysis/animation_ward.gif",
  width = 800,
  height = 600,
  delay = 2  # seconds between frames
)


library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

diagram <- DiagrammeR(
  "graph TD;
           A(Набор признаков)-->B(Матрица расстояний или сходств);
           B-->C(Группировка объектов);
           C-->D(Отношения между кластерами); D-->E(Поправки); E-->A;
           style A fill:#C6FFE7,stroke:#C6FFE7;
           style B fill:#C6FFE7,stroke:#C6FFE7;
           style C fill:#C6FFE7,stroke:#C6FFE7;
           style D fill:#C6FFE7,stroke:#C6FFE7;
           style E fill:#C6FFE7,stroke:#C6FFE7;
  ")


webshot::install_phantomjs()
htmlwidgets::saveWidget(diagram, "temp_diagram.html")
webshot::webshot("temp_diagram.html", "images/diagrammer_cluster.png")
