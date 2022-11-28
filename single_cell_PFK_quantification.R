library(tidyverse)
library(ggforce)
source("geom_hpline.R")
library(ggpubr)

rename_nested_columns <- function(.df, .prefix) {
  .df %>% 
    rename_with(tolower) %>% 
    rename(cell = `...1`) %>% 
    rename_with(~paste0(.prefix, .), .cols = 2:6)
}

# Import and join data ----
dpath <- "data/pfk/results_csv/"

dat_cyto <-
  tibble(
    dpaths = list.files(dpath, full.names = T, pattern = "cytoplasma-signal.csv"),
    data = map(dpaths, read_csv)) %>%
  mutate(dpaths = basename(dpaths)) %>%
  separate(col = dpaths,
           into = c("exp_id", "donor", "stiffness", "treatment", 
                    "sample", "image"),
           sep = '_', extra = 'drop') %>% 
  mutate(data = map(data, rename_nested_columns, "cyto_")) %>% 
  unnest(data)

dat_nuc <-
  tibble(
    dpaths = list.files(dpath, full.names = T, pattern = "nuclear-signal.csv"),
    data = map(dpaths, read_csv)) %>%
  mutate(dpaths = basename(dpaths)) %>%
  separate(col = dpaths,
           into = c("exp_id", "donor", "stiffness", "treatment", 
                    "sample", "image"),
           sep = '_', extra = 'drop') %>% 
  mutate(data = map(data, rename_nested_columns, "nuc_")) %>% 
  unnest(data)

dat <- full_join(dat_cyto, dat_nuc)
summary(dat)


## Filter tiny/huge cells and calculate ratio ----
dat_filter <- filter(dat, cyto_area < 500)
ggplot(dat_filter, aes(nuc_area)) +
  geom_freqpoly(bins = 100)

dat_ratio <- dat %>% 
  filter(between(cyto_area, 50, 500)) %>% 
  filter(between(nuc_area, 20, 300)) %>% 
  mutate(nuc_cyto_ratio = nuc_mean/cyto_mean) %>% 
  mutate(stiffness = parse_number(stiffness) %>% as_factor(),
         treatment = as_factor(treatment))

summary(dat_ratio)

## Export results ---
# write_csv(dat_ratio, "results/PFK/MW1-29_results.csv")


# Results ----
## Nuc/Cyto ratio ----
stiff_labs <- c("40kPa", "180kPa")
names(stiff_labs) <- c("40", "180")
y_limits <- c(0, 4000)

### Boxplot ----
ratio_box <- ggplot(dat_ratio, aes(treatment, nuc_cyto_ratio, 
                      group = interaction(stiffness, treatment), 
                      color = treatment)) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~stiffness, 
             labeller = labeller(stiffness = stiff_labs)) +
  theme_classic() +
  scale_color_manual(values = c("gray50", "darkred")) +
  labs(x = "Drug treatment",
       y = "Nuc/Cyto PFK ratio")
ratio_box

### Sinaplot ----
ratio_sina <- ggplot(dat_ratio, aes(treatment, nuc_cyto_ratio, 
                      group = interaction(stiffness, treatment), 
                      color = treatment)) +
  geom_sina(pch = 16, show.legend = FALSE) +
  stat_summary(fun.data = ~mean_sdl(., mult =1), 
               geom = "pointrange", 
               color = "black", 
               show.legend = FALSE, 
               # cex = 1, 
               position = position_dodge(0.9)) +
  facet_wrap(~stiffness, 
             labeller = labeller(stiffness = stiff_labs)) +
  theme_classic() +
  scale_color_manual(values = c("gray50", "darkred")) +
  labs(x = "Drug treatment",
       y = "Nuc/Cyto PFK ratio")
ratio_sina

### Density plot ----
ggplot(dat_ratio, aes(nuc_cyto_ratio, color = treatment)) +
  geom_density(show.legend = TRUE) +
  labs(x = "Nuclear/cytoplasmic PFK ratio",
       y = "Density",
       color = "Drug treatment") +
  facet_wrap(~stiffness, 
             labeller = labeller(stiffness = stiff_labs)) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("gray50", "darkred"))


## Numeric Nuc/Cyto ratio results ----
### Summarized per sample ----
dat_ratio %>% 
  group_by(donor, stiffness, treatment, sample) %>% 
  summarise(n = n(), 
            mean = mean(nuc_cyto_ratio),
            sd = sd(nuc_cyto_ratio))


## Cytoplasmic PFK signal only ----
### Signal in segmented cytoplasma area ----
cyto <- ggplot(dat_ratio, aes(treatment, cyto_mean, 
                      group = interaction(stiffness, treatment),
                      color = treatment)) +
  geom_boxplot(show.legend = FALSE) +
  # geom_sina(pch = 16, show.legend = FALSE) +
  # stat_summary(fun.data = ~mean_sdl(., mult =1),
  #              geom = "pointrange",
  #              color = "black",
  #              show.legend = FALSE,
  #              # cex = 1,
  #              position = position_dodge(0.9)) +
  theme_classic() +
  facet_wrap(~stiffness) +
  scale_color_manual(values = c("gray50", "darkred")) +
  scale_y_continuous(limits = y_limits) +
  labs(x = "Drug treatment",
       y = 'Cytoplasma PFK signal [a.u.]')
cyto

## Nuclear PFK signal only ----
### Nucleus size ----
nuc_size <- ggplot(dat_ratio, aes(treatment, nuc_area, 
                      group = interaction(stiffness, treatment),
                      color = treatment)) +
  geom_boxplot(show.legend = FALSE) +
  # geom_sina(pch = 16, show.legend = FALSE) +
  # stat_summary(fun.data = ~mean_sdl(., mult =1),
  #              geom = "pointrange",
  #              color = "black",
  #              show.legend = FALSE,
  #              # cex = 1,
  #              position = position_dodge(0.9)) +
  theme_classic() +
  facet_wrap(~stiffness, 
             labeller = labeller(stiffness = stiff_labs)) +
  scale_color_manual(values = c("gray50", "darkred")) +
  labs(x = "Drug treatment",
       y = 'Nucleus size [um^2]')
nuc_size

### Nucleus signal ----
nuc_signal <- ggplot(dat_ratio, aes(treatment, nuc_mean, 
                      group = interaction(stiffness, treatment),
                      color = treatment)) +
  geom_boxplot(show.legend = FALSE) +
  # geom_sina(pch = 16, show.legend = FALSE) +
  # stat_summary(fun.data = ~mean_sdl(., mult =1),
  #              geom = "pointrange",
  #              color = "black",
  #              show.legend = FALSE,
  #              # cex = 1,
  #              position = position_dodge(0.9)) +
  theme_classic() +
  facet_wrap(~stiffness, 
             labeller = labeller(stiffness = stiff_labs)) +
  scale_color_manual(values = c("gray50", "darkred")) +
  scale_y_continuous(limits = y_limits) +
  labs(x = "Drug treatment",
       y = 'Nucleus PFK signal [a.u.]')
nuc_signal


## Total signal ----
dat_total <- mutate(dat_ratio, 
                    total_mean = (nuc_intden + cyto_intden)/(nuc_area + cyto_area),
                    total_intden = nuc_intden + cyto_intden)

total_signal <- ggplot(dat_total, aes(treatment, total_mean, 
                      group = interaction(stiffness, treatment),
                      color = treatment)) +
  geom_boxplot(show.legend = FALSE) +
  # geom_sina(pch = 16, show.legend = FALSE) +
  # stat_summary(fun.data = ~mean_sdl(., mult =1),
  #              geom = "pointrange",
  #              color = "black",
  #              show.legend = FALSE,
  #              # cex = 1,
  #              position = position_dodge(0.9)) +
  theme_classic() +
  facet_wrap(~stiffness, 
             labeller = labeller(stiffness = stiff_labs)) +
  scale_color_manual(values = c("gray50", "darkred")) +
  scale_y_continuous(limits = y_limits) +
  labs(x = "Drug treatment",
       y = 'Total PFK signal [a.u.]')
total_signal

## Cell count ----
dat_count <- dat_ratio %>% 
  group_by(across(1:4)) %>% 
  summarise(cell_count = n())

count_plot <- ggplot(dat_count, aes(treatment, cell_count, 
                                      group = interaction(stiffness, treatment),
                                      color = treatment)) +
  geom_jitter(pch = 16, width = 0.1, show.legend = FALSE) +
  stat_summary(fun.data = ~mean_sdl(., mult =1), 
               geom = "hpline", 
               color = "black", 
               show.legend = FALSE, 
               # cex = 1, 
               position = position_dodge(0.9)) +
  theme_classic() +
  facet_wrap(~stiffness, 
             labeller = labeller(stiffness = stiff_labs)) +
  scale_color_manual(values = c("gray50", "darkred")) +
  labs(x = "Drug treatment",
       y = 'Cell count [n]')
count_plot


# ggplot(dat_count, aes(treatment, cell_count, color = treatment)) +
#   stat_summary(geom = "hpline", show.legend = FALSE) +
#   stat_summary(geom = "point", 
#                aes(group = interaction(stiffness, sample)), 
#                position = position_jitter(width = 0.07),
#                pch = 1,
#                show.legend = FALSE) +
#   geom_col(show.legend = FALSE, fill = NA) +
#   theme_classic() +
#   scale_color_manual(values = c("black", "darkred")) +
#   labs(x = "Substrate hardness [kPa]",
#        y = "Cell count [n]") +
#   scale_y_continuous(expand = c(0,0))


# Combined plot -----------------------------------------------------------

fig <- ggpubr::ggarrange(ggpubr::ggarrange(total_signal, nuc_signal, cyto, nrow = 1), 
                         # ggpubr::ggarrange(nuc_signal, cyto), 
                         ggpubr::ggarrange(ratio_box, nuc_size, count_plot, nrow = 1), 
                         nrow = 2)
ggpubr::annotate_figure(fig, 
                top = ggpubr::text_grob("MW1-29 overview", 
                                face = "bold", 
                                size = 16))

# ggsave("results/PFK/MW1-29_overview.svg",
       # width = 550, height = 400, dpi = 72, units = "px")



# By sample ---------------------------------------------------------------

# dat_sum <- dat_ratio %>% 
#   group_by(across(1:4)) %>% 
#   summarise(cell_count = n(),
#             cyto_mean = mean(cyto_mean),
#             nuc_mean = mean(nuc_mean),
#             ratio_mean = mean(nuc_cyto_ratio))
# 
# 
# ggplot(dat_sum, aes(treatment, ratio_mean, 
#                     group = interaction(stiffness, treatment),
#                     color = treatment)) +
#   geom_jitter(pch = 16, width = 0.1, show.legend = FALSE) +
#   stat_summary(fun.data = ~mean_sdl(., mult =1), 
#                geom = "hpline", 
#                color = "black", 
#                show.legend = FALSE, 
#                # cex = 1, 
#                position = position_dodge(0.9)) +
#   theme_classic() +
#   facet_wrap(~stiffness, 
#              labeller = labeller(stiffness = stiff_labs)) +
#   scale_color_manual(values = c("gray50", "darkred")) +
#   labs(x = "Drug treatment",
#        # y = 'Cell count [n]'
#        )
