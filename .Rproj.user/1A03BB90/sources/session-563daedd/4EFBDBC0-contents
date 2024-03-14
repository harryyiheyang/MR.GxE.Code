library(ggplot2)
library(data.table)
library(dplyr)
w=fread("~/NewInteraction/newsample_pleiotropy_h2.csv")
w1=w[which(w$population=="Cross"),]
w2=w[which(w$population=="EUR"),]
w1$fill="1"
w1$fill[which(w1$environment=="marginal.effect")]="2"
w2$fill="1"
w2$fill[which(w2$environment=="marginal.effect")]="2"
w1$environment=ordered(w1$environment,levels=c("marginal.effect","current.drinking","regular.drinking","current.smoking","ever.smoking"))
w2$environment=ordered(w2$environment,levels=c("marginal.effect","current.drinking","regular.drinking","current.smoking","ever.smoking"))
w1$annotation <- sprintf("%.3f%%", 100 * w1$h2)
w2$annotation <- sprintf("%.3f%%", 100 * w2$h2)
plotcross <- ggplot(data=w1) +
  geom_bar(aes(x=h2, y=environment, fill=fill), color="grey20", stat="identity") +
  scale_fill_manual(values=c("grey20","#fe4365")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_errorbar(aes(y=environment, xmin=h2-2*h2_se, xmax=h2+2*h2_se), width=0.2, colour="grey20", alpha=0.9, size=0.5) +
  geom_text(aes(y=environment, x=h2+2*h2_se+0.01, label=annotation), hjust=0, vjust=0.5, size=3) +  # 添加注释
  facet_grid(.~trait) +
  theme(plot.title = element_text(size=18)) +
  theme(strip.text.x = element_text(size = 18), strip.text.y = element_text(size = 18)) +
  theme(axis.title.y=element_blank()) +
  ggtitle("A. Heritability Estimates of Cross Population") +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=15)) +
  scale_x_continuous(limits=c(-0.01,0.31), breaks=seq(0,0.3,0.05)) +
  guides(fill="none")


ploteur <- ggplot(data=w2) +
  geom_bar(aes(x=h2, y=environment, fill=fill), color="grey20", stat="identity") +
  scale_fill_manual(values=c("grey20","#fe4365")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_errorbar(aes(y=environment, xmin=h2-2*h2_se, xmax=h2+2*h2_se), width=0.2, colour="grey20", alpha=0.9, size=0.5) +
  geom_text(aes(y=environment, x=h2+2*h2_se+0.01, label=annotation), hjust=0, vjust=0.5, size=3) +  # 添加注释
  facet_grid(.~trait) +
  theme(plot.title = element_text(size=18)) +
  theme(strip.text.x = element_text(size = 18), strip.text.y = element_text(size = 18)) +
  theme(axis.title.y=element_blank()) +
  ggtitle("B. Heritability Estimates of European Population") +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=15)) +
  scale_x_continuous(limits=c(-0.01,0.31), breaks=seq(0,0.3,0.05)) +
  guides(fill="none")

