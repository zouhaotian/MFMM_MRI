library(tidyverse)
library(gghighlight)
library(gridExtra)

dat = read.csv('dataset/ADNI_long.csv')
uID = unique(dat$ID)

set.seed(2050)
ID.rand = sort(sample(uID, size = 15))
dat2 = dat[which(dat$ID %in% ID.rand), ]

dat2$color = 'Others'
fid1 = 84
fid2 = 144
dat2[which(dat2$ID %in% fid1), ]$color = 'Subject A'
dat2[which(dat2$ID %in% fid2), ]$color = 'Subject B'
p1 = ggplot(dat = dat2, aes(x = Time, y = y1, group = ID)) + 
  geom_line(aes(color = color)) + 
  theme_bw() + 
  xlab('Visit Year') + 
  ylab('ADAS-Cog 13') + 
  scale_colour_manual(values = c('grey', 'blue', 'red')) + 
  theme(legend.position = 'none')

p2 = ggplot(dat = dat2, aes(x = Time, y = y2, group = ID)) + 
  geom_line(aes(color = color)) + 
  theme_bw() + 
  xlab('Visit Year') + 
  ylab('RAVLT immediate') + 
  scale_colour_manual(values = c('grey', 'blue', 'red')) + 
  theme(legend.position = 'none')

p3 = ggplot(dat = dat2, aes(x = Time, y = y3, group = ID)) + 
  geom_line(aes(color = color)) + 
  theme_bw() + 
  xlab('Visit Year') + 
  ylab('RAVLT learning') + 
  scale_colour_manual(values = c('grey', 'blue', 'red')) + 
  theme(legend.position = 'none')

p4 = ggplot(dat = dat2, aes(x = Time, y = y4, group = ID)) + 
  geom_line(aes(color = color)) + 
  theme_bw() + 
  xlab('Visit Year') + 
  ylab('MMSE') + 
  scale_colour_manual(values = c('grey', 'blue', 'red')) + 
  theme(legend.position = 'none')

p5 = ggplot(dat = dat2, aes(x = Time, y = y5, group = ID)) + 
  geom_line(aes(color = color)) + 
  theme_bw() + 
  xlab('Visit Year') + 
  ylab('CDR-SB') + 
  scale_colour_manual(values = c('grey', 'blue', 'red')) + 
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = c(1.5, 1),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

cairo_ps(filename = 'plot/spaghetti.eps', width = 12, height = 8)
grid.arrange(p1, p2, p3, p4, p5, nrow = 2, ncol = 3)
dev.off()