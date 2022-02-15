classify_ATL03 <- function(atl03,
                           atl08) {

  atl03 <- readr::read_csv("E:/R functions/icesat2R_test/processed_ATL03_20200401073328_00880706_003_01_GT1R_STRONG.csv")
  atl08 <- readr::read_csv("E:/R functions/icesat2R_test/processed_ATL08_20200401073328_00880706_003_01_GT1R_STRONG_ATL03_class.csv")

  atl08_seg <- readr::read_csv("E:/R functions/icesat2R_test/processed_ATL08_20200401073328_00880706_003_01_GT1R_STRONG.csv")


  atl03_classified <- dplyr::inner_join(atl03, atl08)


}

color_palette = list(ground = "#7C3626",
                     canopy =  "#90A959",
                     `top of canopy` = "#206500",
                     noise = "#E64DB8")

p <- atl03_classified %>%
  filter(classed_pc_flag != "noise") %>%
  ggplot(mapping = aes(x = cum_segment_length,
                       y = h_ph,
                       color = classed_pc_flag)) +
  geom_point() +
  labs(x = "Along profile distance (m)",
       y = "Height (m)",
       color = "ATL08 Class") +
  scale_color_manual(values = color_palette) +
  theme(legend.position = "bottom")

fig <- ggplotly(p)

fig <- fig %>% toWebGL()


test <- atl03_classified %>%
  slice_head(n = 30000)

fig <- plot_ly(data = test, x = ~cum_segment_length, y = ~h_ph, color = ~classed_pc_flag,
               type = "scatter")
