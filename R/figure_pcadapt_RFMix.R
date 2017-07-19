# nIND <- 25
# textfile <- read.table("HYBRID.LOCALANC.0.6")
# nSNP <- nrow(textfile)
# v <- matrix(0, nrow = nSNP, ncol = 2)
# for (i in 1:nIND){
#   textfile <- read.table(paste0("HYBRID.LOCALANC.", (i - 1), ".6"))
#   v[, 1] <- v[, 1] + (2 * textfile[, 1] + textfile[, 2]) / (2 * nIND)
#   v[, 2] <- v[, 1] + (2 * textfile[, 3] + textfile[, 2]) / (2 * nIND)
# }
# 
# subseq <- seq(1, nSNP, by = 20)
# xx <-  c(pos_map[subseq], pos_map[subseq]) * 1e-6
# yy <- c(stat.rfmix[subseq], stat.pcadapt[subseq, 1])
# populus.df <- data.frame(x = xx, y = yy, Method = c(rep("RFMix", length(stat.rfmix[subseq])),
#                                                     rep("pcadapt", length(stat.pcadapt[subseq, 1])))) 
# 
# 
# vert_lines <- (pos_map[gt.df$beg] + pos_map[gt.df$end]) * 1e-6 / 2
# tp.df <- data.frame(Type = "True Positive", vals = vert_lines)
# 
# ggplot2::ggplot(populus.df, aes(x = xx, y = yy, colour = Method)) + 
#   ggplot2::geom_line(size = 1) + geom_hline(yintercept = 3) +
#   theme_bw() + 
#   scale_x_continuous(breaks = seq(0, max(populus.df$x), 5), limits = c(0, max(populus.df$x))) +
#   xlab("SNP position (Mbp)") +
#   ylab("Normalized Balsamifera ancestry") +
#   theme(axis.text=element_text(size=15),
#         axis.title=element_text(size=15, face="bold"),
#         title=element_text(size = 15, face="bold"),
#         legend.text=element_text(size=15),
#         legend.key.height=unit(1, "line"),
#         legend.key.width=unit(3, "line")
#   )
# 
# populus.df %>% ggplot(aes(x=x, y=y, colour=Method)) + 
#   facet_grid(~Method) +   
#   ggplot2::geom_line(size=1) + 
#   geom_hline(yintercept=3) +
#   theme_bw() +
#   xlab("SNP position (Mbp)") +
#   ylab("Normalized Balsamifera ancestry") +
#   annotate("text", 
#            x = vert_lines,
#            y = max(populus.df$y, na.rm = TRUE),
#            label = "*",
#            size = 10) +
#   theme(axis.text=element_text(size=20, face="bold"),
#         axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)),
#         axis.title=element_text(size=20, face="bold"),
#         title=element_text(size=20, face="bold"),
#         strip.text=element_text(size=20, face="bold"),
#         legend.text=element_text(size=20),
#         legend.key.height=unit(3, "line"),
#         legend.key.width=unit(3, "line"))
# 
