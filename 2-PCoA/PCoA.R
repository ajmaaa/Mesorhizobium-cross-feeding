#############################################
# PCoA (Bray–Curtis) with group1 ellipse, group2 shape, group3 color
library(vegan)
library(ggplot2)
library(ggsci)


otu_data <- read.csv("B_asv.csv", row.names = 1, header = TRUE, check.names = FALSE)

group <- read.table("group_data.txt",
                    header = TRUE,
                    sep = "\t",        
                    stringsAsFactors = FALSE)


otu_data <- t(otu_data)

#  Bray–Curtis 
distance <- vegdist(otu_data, method = "bray")

#  PCoA
poca <- cmdscale(distance, k = 2, eig = TRUE)
plot_data <- as.data.frame(poca$points)
colnames(plot_data) <- c("PCoA1", "PCoA2")


poca_plot_data <- cbind(plot_data, group)


eig_percent <- round(poca$eig / sum(poca$eig) * 100, 1)

#  Adonis（ group1）
dune.div <- adonis2(distance ~ group1, data = group, permutations = 999)
dune_adonis <- paste0("Adonis R2 = ",
                      round(dune.div$R2[1], 3),
                      ", P = ",
                      dune.div$`Pr(>F)`[1])


p <- ggplot(poca_plot_data,
            aes(x = PCoA1,
                y = PCoA2,
                color = group3,
                shape = group2)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = group1, fill = group3),
               geom = "polygon",
               alpha = 0.15,
               level = 0.95,
               show.legend = FALSE) +
  scale_color_d3() +
  scale_fill_d3() +
  labs(x = paste0("PCoA1 (", eig_percent[1], "%)"),
       y = paste0("PCoA2 (", eig_percent[2], "%)"),
       caption = dune_adonis) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  theme(legend.title = element_blank(),
        panel.grid = element_blank())

print(p)





