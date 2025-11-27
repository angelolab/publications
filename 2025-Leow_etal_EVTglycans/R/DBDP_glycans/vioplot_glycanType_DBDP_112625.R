# Vioplot of glycans 
# Author: Ke Leow
# Date: 10/01/25

#--------------------------------
# Load packages/functions
#--------------------------------
library(tidyverse)
library(ggpubr)
library(rstatix)
#--------------------------------
# Load data
#--------------------------------
data <- read_csv('data/DBDP_transition/MALDI/library_matched/DBDP_transition_glyRatio_glycanTypes_031825.csv')

#--------------------------------
# Vioplot for each glycan
#--------------------------------
data_int <- data %>% select(R_fucosylated:R_bisecting) 
data_pheno <- data %>% select(mask, Region) 

#organize data for vioplot
data_vioplot <- cbind(data_pheno, data_int) %>% 
  pivot_longer(starts_with(c("R_")), names_to = "glycan_type", values_to = "ratio") %>% 
  mutate(glycan_type = gsub("R_", "", glycan_type))

#plot single gly
input = "highMannose" 
stat.test <- data_vioplot %>%
  filter(glycan_type == input) %>%
  wilcox_test(ratio ~ Region) %>%
  add_xy_position(x = "Region")%>%
  add_significance()

#Compare between groups
pdf(paste0("R_plots/DBDP_transition/dotplot_DBDP_glycanTypes_",input,".pdf"), width = 3, height = 3)
data_vioplot %>% filter(glycan_type == input) %>% 
  ggplot(aes(x = Region, y = ratio))+
  # geom_boxplot(fill = "grey", col = "black", width = 0.5)+
  geom_jitter(position=position_jitter(0.1))+
  theme_classic()+
  ggtitle(input)+
  labs(y = "Relative Intensity"
  )+
  stat_pvalue_manual(stat.test, label = "p.signif")
dev.off()


#########
#plot single gly
input = "bisecting" 
stat.test <- data_vioplot %>%
  filter(glycan_type == input) %>%
  wilcox_test(ratio ~ Region) %>%
  add_xy_position(x = "Region")%>%
  add_significance()

#Compare between groups
pdf(paste0("R_plots/DBDP_transition/dotplot_DBDP_glycanTypes_",input,".pdf"), width = 3, height = 3)
data_vioplot %>% filter(glycan_type == input) %>% 
  ggplot(aes(x = Region, y = ratio))+
  # geom_boxplot(fill = "grey", col = "black", width = 0.5)+
  geom_jitter(position=position_jitter(0.1))+
  theme_classic()+
  ggtitle(input)+
  labs(y = "Relative Intensity"
  )+
  stat_pvalue_manual(stat.test, label = "p.signif")
dev.off()


#########
#plot single gly
input = "biantennary" 
stat.test <- data_vioplot %>%
  filter(glycan_type == input) %>%
  wilcox_test(ratio ~ Region) %>%
  add_xy_position(x = "Region")%>%
  add_significance()

#Compare between groups
pdf(paste0("R_plots/DBDP_transition/dotplot_DBDP_glycanTypes_",input,".pdf"), width = 3, height = 3)
data_vioplot %>% filter(glycan_type == input) %>% 
  ggplot(aes(x = Region, y = ratio))+
  # geom_boxplot(fill = "grey", col = "black", width = 0.5)+
  geom_jitter(position=position_jitter(0.1))+
  theme_classic()+
  ggtitle(input)+
  labs(y = "Relative Intensity"
  )+
  stat_pvalue_manual(stat.test, label = "p.signif")
dev.off()

#########
#plot single gly
input = "tetraantennary" 
stat.test <- data_vioplot %>%
  filter(glycan_type == input) %>%
  wilcox_test(ratio ~ Region) %>%
  add_xy_position(x = "Region")%>%
  add_significance()

#Compare between groups
pdf(paste0("R_plots/DBDP_transition/dotplot_DBDP_glycanTypes_",input,".pdf"), width = 3, height = 3)
data_vioplot %>% filter(glycan_type == input) %>% 
  ggplot(aes(x = Region, y = ratio))+
  # geom_boxplot(fill = "grey", col = "black", width = 0.5)+
  geom_jitter(position=position_jitter(0.1))+
  theme_classic()+
  ggtitle(input)+
  labs(y = "Relative Intensity"
  )+
  stat_pvalue_manual(stat.test, label = "p.signif")
dev.off()

#########
#plot single gly
input = "polylacnac" 
stat.test <- data_vioplot %>%
  filter(glycan_type == input) %>%
  wilcox_test(ratio ~ Region) %>%
  add_xy_position(x = "Region")%>%
  add_significance()

#Compare between groups
pdf(paste0("R_plots/DBDP_transition/dotplot_DBDP_glycanTypes_",input,".pdf"), width = 3, height = 3)
data_vioplot %>% filter(glycan_type == input) %>% 
  ggplot(aes(x = Region, y = ratio))+
  # geom_boxplot(fill = "grey", col = "black", width = 0.5)+
  geom_jitter(position=position_jitter(0.1))+
  theme_classic()+
  ggtitle(input)+
  labs(y = "Relative Intensity"
  )+
  stat_pvalue_manual(stat.test, label = "p.signif")
dev.off()


# #plot single gly
# input = "H7N6F1" 
# # pdf(paste0("R_plots/MALDI_IF_EVT_glycans/vioplot_EVTglycans_mean_int_",input,".pdf"), width = 4, height = 4)
# data_vioplot %>% filter(glycan == input) %>% 
#   ggplot(aes(x = cell_anno_update, y = average_pixel_intensity))+
#   geom_violin(alpha = 0.7) + 
#   geom_jitter(position=position_jitter(0.1))+
#   stat_summary(fun = "mean", geom = "crossbar",  width = 0.5,) +
#   theme_classic()+
#   ggtitle(input)+ theme(legend.position="none")
# dev.off()
# 
# input = "H5N5F1" 
# pdf(paste0("R_plots/MALDI_IF_EVT_glycans/vioplot_EVTglycans_mean_int_",input,".pdf"), width = 4, height = 4)
# data_vioplot %>% filter(glycan == input) %>% 
#   ggplot(aes(x = cell_anno_update, y = average_pixel_intensity))+
#   geom_violin(alpha = 0.7) + 
#   geom_jitter(position=position_jitter(0.1))+
#   stat_summary(fun = "mean", geom = "crossbar",  width = 0.5,) +
#   theme_classic()+
#   ggtitle(input)+ theme(legend.position="none")
# dev.off()
# 
# input = "H6N2" 
# pdf(paste0("R_plots/MALDI_IF_EVT_glycans/vioplot_EVTglycans_mean_int_",input,".pdf"), width = 4, height = 4)
# data_vioplot %>% filter(glycan == input) %>% 
#   ggplot(aes(x = cell_anno_update, y = average_pixel_intensity))+
#   geom_violin(alpha = 0.7) + 
#   geom_jitter(position=position_jitter(0.1))+
#   stat_summary(fun = "mean", geom = "crossbar",  width = 0.5,) +
#   theme_classic()+
#   ggtitle(input)+ theme(legend.position="none")
# dev.off()
# 
# input = "H5N4S1" 
# pdf(paste0("R_plots/MALDI_IF_EVT_glycans/vioplot_EVTglycans_mean_int_",input,".pdf"), width = 4, height = 4)
# data_vioplot %>% filter(glycan == input) %>% 
#   ggplot(aes(x = cell_anno_update, y = average_pixel_intensity))+
#   geom_violin(alpha = 0.7) + 
#   geom_jitter(position=position_jitter(0.1))+
#   stat_summary(fun = "mean", geom = "crossbar",  width = 0.5,) +
#   theme_classic()+
#   ggtitle(input)+ theme(legend.position="none")
# dev.off()
# 
# 
# # # plot all gly into pdf
# # gly <- colnames(data_int)
# # pdf("figures/vioplot_allGlycans_evtVctFv.pdf", width = 8, height = 6)
# # for (i in 1:length(gly)) {
# #   input = gly[i]
# #   
# #   plot = data_vioplot %>% filter(glycan == input) %>% 
# #     ggplot(aes(x = cell_anno_update, y = average_pixel_intensity))+
# #     geom_jitter(position=position_jitter(0.1))+
# #     geom_violin(aes(fill = cell_anno_update), alpha = 0.7) + 
# #     stat_summary(fun = "mean", geom = "crossbar",  width = 0.5,) +
# #     theme_classic()+
# #     ggtitle(input)
# #   
# #   print(plot)
# # }
# # dev.off()