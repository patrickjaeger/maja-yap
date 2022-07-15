library(tidyverse)

dpath <- "data/results_csv"

dat <-
  tibble(
    dpaths = list.files(dpath, full.names = T),
    data = map(dpaths, read_csv)) %>%
  mutate(dpaths = basename(dpaths)) %>%
  separate(col = dpaths,
           into = c("x1", "x2", "stiffness", "antibody", 
                    "magnification", "x3"),
           sep = '_', extra = 'drop')

dat

datl <- unnest(dat, data) %>% 
  rename(cell_id = "...1", 
         area = "Area", 
         mean_intensity = "Mean",
         median_intensity = "Median")

ggplot(datl, aes(stiffness, mean_intensity)) +
  geom_jitter(width = 0.2)
