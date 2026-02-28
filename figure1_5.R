library(ggplot2)
library(dplyr)

df_pie <- data.frame(
  type  = c("Non-cluster", "Cluster"),
  count = c(147304, 429629)
) %>%
  mutate(
    percent = count / sum(count) * 100,
    label = sprintf("%s\n%.1f%%", type, percent)
  )

total_n <- sum(df_pie$count)

ggplot(df_pie, aes(x = "", y = count, fill = type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 5
  ) +
  labs(
    title    = "Cluster vs Non-cluster Cell Distribution",
    subtitle = paste0("Total cells = ", format(total_n, big.mark = ",")),
    fill     = "Category"
  ) +
  theme_void() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.title  = element_text(size = 12),
    legend.text   = element_text(size = 11)
  )
