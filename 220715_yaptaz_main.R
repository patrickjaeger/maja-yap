library(tidyverse)
library(ggforce)
source("geom_hpline.R")

dpath <- "data/results_csv"

dat <-
  tibble(
    dpaths = list.files(dpath, full.names = T),
    data = map(dpaths, read_csv)) %>%
  mutate(dpaths = basename(dpaths)) %>%
  separate(col = dpaths,
           into = c("x1", "x2", "stiffness", "antibody", 
                    "magnification", "sample"),
           sep = '_', extra = 'drop')

dat

datl <- unnest(dat, data) %>% 
  rename(cell_id = "...1", 
         area = "Area", 
         mean_intensity = "Mean",
         median_intensity = "Median",
         normalized_mean = "normalized_mean")

# datl <- mutate(datl, stiffness = as_factor(stiffness) %>% 
#                  fct_relevel("40kPa", "2kPa"))


## YAP ratio ----
datl <- mutate(datl, stiffness = parse_number(stiffness) %>% as_factor())

ggplot(datl, aes(stiffness, normalized_mean, color = stiffness)) +
  geom_boxplot(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("black", "darkred")) +
  labs(x = "PDMS hardness [kPa]",
       y = "Nuc/Cyto YR")

ggplot(datl, aes(stiffness, normalized_mean, color = stiffness)) +
  geom_sina(pch = 16, show.legend = FALSE) +
  stat_summary(fun.data = ~mean_sdl(., mult =1), geom = "pointrange", color = "black", show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("gray50", "darkred")) +
  labs(x = "PDMS hardness [kPa]",
       y = "Nuc/Cyto YR")
  
# ggplot(datl, aes(normalized_mean, color = stiffness)) +
#   geom_density() +
#   theme_classic() +
#   scale_color_manual(values = c("black", "darkred")) +
#   labs(x = "Nuclear/cytoplasmic YAP ratio",
#        y = "Density",
#        color = "Substrate hardness") +
#   theme(legend.position = "top")


### numeric YR results ----
datl %>% 
  group_by(stiffness, sample) %>% 
  summarise(n = n(), 
            mean = mean(normalized_mean),
            sd = sd(normalized_mean))


## Cell count ----
dat_count <- datl %>% group_by(across(1:6)) %>% summarise(cell_count = n())

ggplot(dat_count, aes(stiffness, cell_count, color = stiffness)) +
  stat_summary(geom = "hpline", show.legend = FALSE) +
  stat_summary(geom = "point", 
               aes(group = interaction(stiffness, sample)), 
               position = position_jitter(width = 0.07),
               pch = 1,
               show.legend = FALSE) +
  geom_col(show.legend = FALSE, fill = NA) +
  theme_classic() +
  scale_color_manual(values = c("black", "darkred")) +
  labs(x = "Substrate hardness [kPa]",
       y = "Cell count [n]") +
  scale_y_continuous(expand = c(0,0))
