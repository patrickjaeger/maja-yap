library(tidyverse)
library(ggforce)
source("geom_hpline.R")

rename_nested_columns <- function(.df, .prefix) {
  .df %>% 
    rename_with(tolower) %>% 
    rename(cell = `...1`) %>% 
    rename_with(~paste0(.prefix, .), .cols = 2:6)
}

# Import and join data ----
dpath <- "data/donut_test/results_csv"

dat_cyto <-
  tibble(
    dpaths = list.files(dpath, full.names = T, pattern = "cytoplasmicYAP.csv"),
    data = map(dpaths, read_csv)) %>%
  mutate(dpaths = basename(dpaths)) %>%
  separate(col = dpaths,
           into = c("exp_id", "passage", "stiffness", "antibody", 
                    "magnification", "sample", "type"),
           sep = '_', extra = 'drop') %>% 
  select(-type) %>% 
  mutate(data = map(data, rename_nested_columns, "cyto_")) %>% 
  unnest(data)

dat_nuc <-
  tibble(
    dpaths = list.files(dpath, full.names = T, pattern = "nuclearYAP.csv"),
    data = map(dpaths, read_csv)) %>%
  mutate(dpaths = basename(dpaths)) %>%
  separate(col = dpaths,
           into = c("exp_id", "passage", "stiffness", "antibody", 
                    "magnification", "sample", "type"),
           sep = '_', extra = 'drop') %>% 
  select(-type) %>% 
  mutate(data = map(data, rename_nested_columns, "nuc_")) %>% 
  unnest(data)

dat <- full_join(dat_cyto, dat_nuc)

dat_ratio <- mutate(dat, yap_ratio = nuc_mean/cyto_mean)


# Results ----
## YAP ratio ----
dat_ratio <- mutate(dat_ratio, stiffness = parse_number(stiffness) %>% as_factor())

ggplot(dat_ratio, aes(stiffness, yap_ratio, color = stiffness)) +
  geom_boxplot(show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = c("black", "darkred")) +
  labs(x = "PDMS hardness [kPa]",
       y = "Nuc/Cyto YR")

ggplot(dat_ratio, aes(stiffness, yap_ratio, color = stiffness)) +
  geom_sina(pch = 16, show.legend = FALSE) +
  stat_summary(fun.data = ~mean_sdl(., mult =1), geom = "pointrange", color = "black", show.legend = FALSE, cex = 1) +
  theme_classic() +
  scale_color_manual(values = c("gray50", "darkred")) +
  labs(x = "PDMS hardness [kPa]",
       y = "Nuc/Cyto YR")
  
ggplot(dat_ratio, aes(yap_ratio, color = stiffness)) +
  geom_density() +
  theme_classic() +
  scale_color_manual(values = c("black", "darkred")) +
  labs(x = "Nuclear/cytoplasmic YAP ratio",
       y = "Density",
       color = "Substrate hardness") +
  theme(legend.position = "top")


## numeric YR results ----
dat_ratio %>% 
  group_by(stiffness, sample) %>% 
  summarise(n = n(), 
            mean = mean(yap_ratio),
            sd = sd(yap_ratio))


## Cell count ----
dat_count <- dat_ratio %>% group_by(across(1:6)) %>% summarise(cell_count = n())

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
