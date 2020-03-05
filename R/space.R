library(ggplot2)
library(stringr)

ind.loc <- read.csv(file = "~/ProjetoFinalJB/data/indloc.csv", sep = ",")

max.space=100

for (i in 1:10) {
  individuals.fig=
    ggplot() +  
    geom_point(data = subset(ind.loc, time == i), aes(x = x, y = y,color= factor(sp)), size=1.5, alpha=0.5)+
    guides(fill=FALSE, shape="none") +
    labs(x = "", y="") +  
    xlim(0, 100) +
    ylim(0, 100) +
    theme_bw()+
    theme(text = element_text(size=12, family="Helvetica"),
          panel.grid.minor = element_blank(),panel.grid.major = element_blank(), 
          legend.title=element_blank(),
          legend.position = c(1.12,0.365), 
          legend.background = element_rect(fill="transparent",size=0.01, linetype="solid",colour ="black"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.text=element_text(size=5),
          legend.box.background = element_rect(),
          legend.margin = margin(-4, 4, -1, -1),
          plot.margin = unit(c(0.1,2,0.1,0.1), "cm"))
  
  individuals.fig
  # Set here to a place to save the images
  ggsave(paste("../figs/position/ind_location_", str_pad(i, 2, pad="0"), ".png", sep = ""), individuals.fig,width=6, height=5)
}