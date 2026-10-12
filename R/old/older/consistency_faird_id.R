library(tidyverse)
library(ggplot2)
library(ggpubr)


data <- read_csv("etraps_daily_biomass_clean.csv")

# --- Función para análisis de correlación completo ---
analyze_device_consistency <- function(data, device1, device2, metric = "biomass") {
  
  # Filtrar datos para los dos dispositivos
  df <- data %>%
    filter(Device %in% c(device1, device2)) %>%
    select(Date, Device, !!sym(metric)) %>%
    pivot_wider(names_from = Device, values_from = !!sym(metric))
  
  # Remover NAs
  df_clean <- df %>% drop_na()
  
  # Correlación de Pearson
  cor_pearson <- cor.test(df_clean[[device1]], df_clean[[device2]], 
                          method = "pearson")
  
  # Correlación de Spearman
  cor_spearman <- cor.test(df_clean[[device1]], df_clean[[device2]], 
                           method = "spearman")
  
  # Regresión lineal
  model <- lm(df_clean[[device2]] ~ df_clean[[device1]])
  r_squared <- summary(model)$r.squared
  
  # Gráfico de dispersión
  p <- ggplot(df_clean, aes(x = !!sym(device1), y = !!sym(device2))) +
    geom_point(size = 3, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.9) +
    stat_regline_equation(label.x.npc = 0.1, label.y.npc = 0.85) +
    labs(
      title = paste("Consistency Analysis:", device1, "vs", device2),
      x = paste(device1, metric, "(mg)"),
      y = paste(device2, metric, "(mg)")
    ) +
    theme_minimal()
  
  # Resultados
  results <- list(
    pearson_r = cor_pearson$estimate,
    pearson_p = cor_pearson$p.value,
    spearman_rho = cor_spearman$estimate,
    spearman_p = cor_spearman$p.value,
    r_squared = r_squared,
    model = model,
    plot = p,
    data = df_clean
  )
  
  return(results)
}

# --- Ejecutar análisis para cada par ---
# FAIRD en Maize
faird_maize <- analyze_device_consistency(data, "FAIRD1", "FAIRD2", "biomass_mg")
print(faird_maize$plot)

# FAIRD en Meadow
faird_meadow <- analyze_device_consistency(data, "FAIRD3", "FAIRD4", "biomass_mg")
print(faird_meadow$plot)

# ID en Maize
id_maize <- analyze_device_consistency(data, "ID1", "ID2", "biomass_mg")
print(id_maize$plot)