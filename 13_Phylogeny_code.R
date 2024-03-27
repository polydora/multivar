## Филогения с использованием R

## Анастасия Лянгузова

library(ggtree)
library(ggplot2)
library(treeio)

library(ape)
ids <- readLines("data/genbank_ids.txt")
seqs <- read.GenBank(ids) # позволяет извлечь fasta-файлы из базы данныx
str(seqs)

# чтобы было удобнее читать
seqs_ids <- paste(attr(seqs, "species"), names(seqs), sep = "_")
write.dna(seqs, file = "data/rhiz_ids.fasta", format = "fasta", nbcol = 1) # сохранение файла в формате fasta посредство пакета ape
#

library(seqinr)
rhiz_fasta <- read.fasta(file = "data/rhiz_ids.fasta", seqtype = "DNA",
           as.string = TRUE, forceDNAtolower = FALSE) # чтение файла
write.fasta(sequences = rhiz_fasta, names = seqs_ids, file.out = "data/rhiz_not_align.fasta")

## расчет матрицы коэфициентов сходств/различий
library(adegenet)
# выравнивание
rhiz_dnabin <- fasta2DNAbin('data/rhiz_not_align.fasta')
library(ape)
clustal_align <- clustal(rhiz_dnabin)
image(clustal_align)


mafft_fasta <- fasta2DNAbin('data/18S_alignment.fasta') # читаем в рабочем формате
image(mafft_fasta)

library(phangorn)
sub_models <- modelTest(as.phyDat(mafft_fasta), multicore = TRUE,
                        mc.cores = 6)
(best <- sub_models[which.min(sub_models$AIC), ])
(worst <- sub_models[which.max(sub_models$AIC), ])


# чтение деревьев в формате Nexus (баейсовская филогения)
library(treeio)
library(ggtree)
nex_rhiz <- read.beast("data/18S_ver2.nex.con.tre") # GTR+I+G, MAFFT
# просто график по филогении (но с поддержками)
gtr_tree <- ggtree(nex_rhiz) +
  geom_tiplab(size = 4) +
  geom_text(aes(label = prob_percent),
            hjust = 1.2, vjust = -0.3)

# структура данных
str(nex_rhiz)
str(nex_rhiz@phylo)
nex_rhiz@phylo$tip.label



# укоренение
# library(phytools)
# nex_rhiz@phylo <- midpoint.root(nex_rhiz@phylo)
library(ape)
nex_rhiz@phylo <- root(nex_rhiz@phylo, node = 43)

# создаем группы
group_rhiz <- groupOTU(nex_rhiz@phylo, c("Spug_18S_14_5",
                                       "Spug_18S_16_1",
                                       "Spug_18S_19_4",
                                       "Spug_18S_24_2"))

(gg_spug <- ggtree(nex_rhiz) +
  geom_tiplab(size = 4) +
  geom_treescale() +
  geom_tippoint(data = group_rhiz,
                aes(alpha = group), col = "red") +
  geom_text(aes(label = prob_percent),
            hjust = 1.2, vjust = -0.3) +
  theme(legend.position = 'null'))

sac_para_nodes <- list("Sacculina pugettiae" = c("Spug_18S_14_5",
                                                 "Spug_18S_16_1",
                                                 "Spug_18S_19_4",
                                                 "Spug_18S_24_2"),
                       "Parasacculina pilosella" = nex_rhiz@phylo$tip.label[18:24])

sac_para_otus <- groupOTU(nex_rhiz, sac_para_nodes)
ggtree(sac_para_otus) +
  geom_tiplab(size = 4) +
  geom_nodelab() +
  geom_treescale() +
  geom_tippoint(aes(col = group)) +
  scale_color_manual(values=c("black", "blue", "orange")) +
  geom_text(aes(label = prob_percent),
            hjust = 1.2, vjust = -0.3)
  # theme(legend.position = 'null')

# визуализация номеров узлов
ggtree(sac_para_otus) +
  geom_tiplab(size = 4) +
  geom_treescale() +
  geom_text(aes(label=node), hjust=1.2, vjust = -.3)

# окраска нужных узлов
ggtree(sac_para_otus) +
  geom_tiplab(size = 4) +
  geom_nodelab() +
  geom_treescale() +
  geom_hilight(node = 63, fill = "orange") +
  geom_hilight(node = 58, fill = "blue") +
  geom_text(aes(label = prob_percent),
          hjust = 1.2, vjust = -0.3)

# выделение групп на дереве
ggtree(sac_para_otus) +
  geom_tiplab(size = 4) +
  geom_text(aes(label = prob_percent),
            hjust = 1.2, vjust = -0.3) +
  geom_cladelabel(node = 63, label="Sacculina pugettiae",
                color = 'orange', offset = 0.3,
                hjust = -0.1) +
  geom_cladelabel(node = 58, label="Parasacculina pilosella",
                color = 'blue', offset = 0.28,
                hjust = -0.1)



# сравниваем с другим деревом
jc69 <- read.tree("data/18S_JC69_sac_parasac.treefile")

#### Задание ####
## Визуализируйте дерево, полученное с помощью JC69 и укорените его соответственно той же внешней группе, что и построенное по GTR ####



# Построение кофилогении
library(ggpubr)
ggarrange(gtr_tree, jc69_tree)

library(ape)
cophyloplot(nex_rhiz@phylo, jc69, length.line=4, space=40)

library(phytools)
trees.cophylo <- cophylo(nex_rhiz@phylo, jc69, rotate = TRUE)
png("cophyplot.png", width = 1200, height = 800)
plot(trees.cophylo, link.type = "curved", link.lwd = 4,
     link.lty="solid", link.col = "blue", size = 1)
dev.off()


#### Задание ####
### Трематоды из семейства Zoogonidae (https://www.mdpi.com/1424-2818/15/1/121) Kremnev, G., Gonchar, A., Uryadova, A., Krapivin, V., Skobkina, O., Gubler, A., & Krupenko, D. (2023). No Tail No Fail: Life Cycles of the Zoogonidae (Digenea). Diversity, 15(1), 121.
### 1. Выберите нужную модель нуклеотидных замен, основываясь на выравнивании
### 2. Сравните деревья, полученные методом максимального правдоподобия и байесовским методом. Выделите на дереве сем. Zoogonidae
