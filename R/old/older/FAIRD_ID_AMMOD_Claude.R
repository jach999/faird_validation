# --- DIAGNÓSTICO DE WARNINGS ---

# 1. Verificar problema de Checkout
cat("=== Problema 1: Checkout que no parseó ===\n")
faird_raw <- read_csv("FAIRD_For_Statistics.csv", show_col_types = FALSE)
id_raw <- read_csv("ID_For_Statistics.csv", show_col_types = FALSE)

# Buscar el registro problemático
faird_problem <- faird_raw %>% 
  mutate(parsed = !is.na(ymd_hms(Checkout))) %>%
  filter(!parsed)

id_problem <- id_raw %>%
  mutate(parsed = !is.na(ymd_hms(Checkout))) %>%
  filter(!parsed)

if(nrow(faird_problem) > 0) {
  cat("FAIRD Checkout problemático:\n")
  print(faird_problem %>% select(Device_ID, Checkin, Checkout))
}
if(nrow(id_problem) > 0) {
  cat("ID Checkout problemático:\n")
  print(id_problem %>% select(Device_ID, Checkin, Checkout))
}

# 2. Verificar problema de parsing general
cat("\n=== Problema 2: Parsing issues ===\n")
ammod_tax_raw <- read_csv("AMMOD_24h_otu_counts_For_Statistics.csv",
                          show_col_types = FALSE)
cat("Primeras fechas de AMMOD taxonomy:\n")
print(head(ammod_tax_raw$Date, 10))
cat("\nClase de Date column:\n")
print(class(ammod_tax_raw$Date))

# Verificar problemas
probs <- problems(ammod_tax_raw)
if(nrow(probs) > 0) {
  cat("\nProblemas encontrados:\n")
  print(probs)
}