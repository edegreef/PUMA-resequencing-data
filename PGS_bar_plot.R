ggplot(test_individuals, aes(x = decile, y = spring_migration, group= decile)) + 
  geom_jitter(shape = 16, position = position_jitter(0.1), alpha = 0.6, size=0.5, col="gray80") +
  geom_boxplot(col=c("#49417B", "#6055A1", "#786EB2", "#9363AD", "#A47CBA", "#B16BA8", "#BB7FB4", "#AE6583", "#C28CA2", "#CD9FB2"), fill=NA, lwd=1) +
 # geom_boxplot(col=c("#786EB2", "#9363AD","#B16BA8", "#AE6583"), fill=NA, lwd=1) +
  #scale_colour_brewer(palette="PuOr")+
  ylab("Spring migration timing") +
  xlab("Decile") +
  scale_x_continuous(breaks = 1:10, labels = c(as.character(1:10))) +
  theme_classic()
 # xlim(c(3.5,7.5))

ggsave("spring_cov7_PGS_BAR_1-10decile_plot_updated.png", width=4, height=3, dpi=1000)

# point thing
brewer.pal(11, 'PuOr')
ggplot(dat2, aes(x=SCORE, y=spring_migration))+#, colour=Colony))+
  geom_smooth(method="lm",formula = y ~ x, linetype=0, fill="gray60")+
  geom_point(shape=16,size=1.8, alpha=1, aes(color=factor(Colony, levels=c("AB", "MN", "SD", "ON", "PA", "NJ", "VA", "TAM", "FLD", "TX", "FLN"))))+
  theme_classic()+
  scale_color_manual(values=brewer.pal(11, 'PuOr'))+
  labs(x="Polygenic score", y="Spring migration timing")+
  labs(color="Colony")

ggsave("spring_cov7_PGS_plot_updated_LEGENDSPEC.png", width=4.5, height=4, dpi=1000)
