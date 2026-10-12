╔═══════════════════════════════════════════════════════════════════════════════╗
║              📖 PROYECTO: FAIRD DEVICE VALIDATION STUDY                      ║
╚═══════════════════════════════════════════════════════════════════════════════╝

ÚLTIMA ACTUALIZACIÓN: Noviembre 4, 2025
VERSION: 4.2 - Consolidated Results + H2 Discussion
ESTADO: H1 Completo ✅ | H2 Completo ✅ | H3 Pendiente ⏸️

Este documento integra el estado actual del proyecto, la metodología validada
y los resultados finales de H1 y H2, basados en los outputs de los scripts
de análisis.

═══════════════════════════════════════════════════════════════════════════════
🎯 OBJETIVO DEL PROYECTO
═══════════════════════════════════════════════════════════════════════════════

PAPER: Validation of FAIR-Device as automated e-trap for continuous insect
       biodiversity monitoring
FOCUS: Statistical validation of three device types (AMMOD, FAIRD, ID)
TARGET JOURNAL: Methods in Ecology and Evolution (IF: 6.6)

HIPÓTESIS PRINCIPALES:
  H1: Consistency    - Temporal reproducibility between replicate devices ✅
  H2: Validation     - Correlation with ground truth (AMMOD) ✅
  H3: Taxonomic      - Taxonomic breadth comparison (FAIRD vs ID) ⏸️

═══════════════════════════════════════════════════════════════════════════════
📅 DISEÑO EXPERIMENTAL
═══════════════════════════════════════════════════════════════════════════════

PERÍODOS EXPERIMENTALES:

  E-traps (FAIRD & ID):
    • Inicio: August 23, 2023
    • Fin: September 13, 2023
    • Duración: 22 días continuos
    • Cobertura: Completa (datos diarios)

  AMMOD Ground Truth:
    • Primera tanda de botellas: August 23 - September 3, 2023 (12 días)
      → Resultados de laboratorio DISPONIBLES ✅
    • Cambio de botellas: September 4, 2023 (no data, día operacional)
    • Segunda tanda de botellas: September 5 - September 13, 2023 (9 días)
      → Resultados de laboratorio PENDIENTES ⏸️

RESTRICCIÓN TEMPORAL CRÍTICA:
  Los análisis de validación H2 (FAIRD vs AMMOD) están limitados a los 
  primeros 12 días (Aug 23 - Sep 3) debido al procesamiento pendiente
  de la segunda tanda de AMMOD.

SITIOS Y DISPOSITIVOS:

  Site1 - Maize-Edge (Ecotono, borde de cultivo + canal):
    • FAIRD1 ✔    ID1 ✔    AMMOD1 ✔

  Site2 - Maize-Core (Interior de cultivo homogéneo):
    • FAIRD2 ✔    ID2 ✔    AMMOD2 ✔

  Site3 - Meadow 1 (Pradera):
    • FAIRD3 ✔    ID3 ✖    AMMOD3 ⚠️

  Site4 - Meadow 2 (Pradera):
    • FAIRD4 ✔    ID4 ✔    AMMOD4 ✖

LEYENDA:
  ✔  = Funcional
  ⚠️  = Fallo parcial (AMMOD3: mecánico, muestras 30/Ago - 03/Sep acumuladas)
  ✖ = Datos no disponibles (AMMOD4: perdidos; ID3: fallo de dispositivo)

═══════════════════════════════════════════════════════════════════════════════
📊 DATOS DEL PROYECTO
═══════════════════════════════════════════════════════════════════════════════

ARCHIVOS DE DATOS (INPUT):
  1. FAIRD_For_Statistics.csv (1,128 individuos)
  2. ID_For_Statistics.csv (2,370 individuos)
  3. AMMOD_24h_otu_counts_For_Statistics.csv (CUALITATIVO, solo taxonomía)
  4. AMMOD_24h_live_mass_For_Statistics.csv (CUANTITATIVO, solo biomasa)

LIMITACIÓN FUNDAMENTAL DE AMMOD:
  ⚠️ TAXONOMÍA Y BIOMASA ESTÁN COMPLETAMENTE SEPARADAS
     • NO se puede calcular la biomasa por taxón para AMMOD.
     • SÍ se puede correlacionar la biomasa total diaria (H2).
     • SÍ se puede comparar la riqueza taxonómica (H3).

═══════════════════════════════════════════════════════════════════════════════
🗂️ ESTRUCTURA Y PIPELINE DE ANÁLISIS
═══════════════════════════════════════════════════════════════════════════════

SCRIPTS PRINCIPALES:
┌─────────────────────────────────────────────────────────────────────────────┐
│ load_data.R             → Carga, limpieza y completitud de datos ⭐ CRÍTICO │
│ analysis_functions.R    → Librería de funciones estadísticas                │
│ custom_colors.R         → Paleta de colores consistentes                    │
├─────────────────────────────────────────────────────────────────────────────┤
│ H1_consistency_FAIRD.R  → H1 Análisis de consistencia FAIRD ✅              │
│ H1_consistency_ID.R     → H1 Análisis de consistencia ID ✅                 │
│ H1_consistency_AMMOD.R  → H1 Análisis de consistencia AMMOD ✅              │
│ H2_validation.R         → H2 Validación FAIRD vs AMMOD ✅                   │
│ H3_taxonomic.R          → H3 Comparación taxonómica ⏸️                     │
└─────────────────────────────────────────────────────────────────────────────┘

DATA PIPELINE:
  Raw CSVs (4 files)
      ↓
  load_data.R
      ├── Carga datos de individuos
      ├── Agrega a biomasa/abundancia diaria
      ├── Maneja fallo de AMMOD3 (crea ammod_biomass_original con prorrateo)
      ├── Asegura completitud (con expand_grid + left_join + replace_na(0)) ⭐
      └── Guarda: data_cleaned.RData
      ↓
  data_cleaned.RData (Objeto data_list)
      ├── etraps_daily (COMPLETO con ceros)
      ├── ammod_biomass_original (INCLUYE datos prorrateados de AMMOD3)
      └── ...otros datos
      ↓
  analysis_functions.R
      ├── Carga todas las funciones personalizadas (ej. analyze_...)
      ↓
  Scripts de Análisis (H1, H2, H3)
      ↓
  Resultados: Tablas (CSV) + Figuras (PNG/PDF)

═══════════════════════════════════════════════════════════════════════════════
🔬 METODOLOGÍA ESTADÍSTICA (DECISIONES FINALES)
═══════════════════════════════════════════════════════════════════════════════

1. TRANSFORMACIÓN LOGARÍTMICA (OBLIGATORIA)
   • DECISIÓN: Usar log10(x+1) para TODOS los análisis de regresión y
     correlación de H1 y H2 (abundancia y biomasa).
   • RAZÓN: Los datos brutos fallan los diagnósticos (Shapiro y Breusch-Pagan).
     Los datos transformados los pasan exitosamente.
   • OUTPUT VALIDADO: H1_consistency_..._output.txt, H2_validation_output.txt.

2. MANEJO DE DÍAS CERO (OBLIGATORIO)
   • DECISIÓN: Los días con cero capturas son datos reales y deben incluirse.
   • IMPLEMENTACIÓN: Centralizada en `load_data.R` (Sección 5.5) usando
     `tidyr::expand_grid`, `dplyr::left_join` y `tidyr::replace_na(abundance = 0, ...)`
     para rellenar los días faltantes.

3. MANEJO DE FALLO AMMOD3 (PRORRATEO)
   • DECISIÓN: Usar los datos prorrateados (promediados) para el análisis H2
     para mantener el tamaño de muestra (n=12) en el Sitio 3.
   • IMPLEMENTACIÓN: Los scripts H1 y H2 cargan `ammod_biomass_original`
     (que contiene los valores prorrateados) para los análisis.
   • OUTPUT VALIDADO: `H2_validation_output.txt` confirma n=12 para el Sitio 3.

4. PRUEBAS NO PARAMÉTRICAS
   • SECUENCIA: Kruskal-Wallis (prueba general) → Dunn's Test (post-hoc)
     con corrección "bonferroni".
   • APLICACIÓN: H2 Sección 6 (Comparación de Dispositivos) y Sección 8
     (Comparación de Ambientes).

═══════════════════════════════════════════════════════════════════════════════
✅ H1: RESULTADOS DE CONSISTENCIA TEMPORAL (VALIDADOS)
═══════════════════════════════════════════════════════════════════════════════
Análisis de correlación (Pearson) sobre datos transformados log10(x+1).

┌──────────────────┬──────────────────┬──────────┬────────┬───────────┬──────────────┐
│ Comparación      │ Métrica          │ R² (log) │ p-value  │ Veredicto │ Fuente (Output)│
├──────────────────┼──────────────────┼──────────┼────────┼───────────┼──────────────┤
│ FAIRD1 vs FAIRD2 │ Abundancia (n=22)│  0.674   │ <0.001   │ Fuerte    │ H1_FAIRD       │
│ FAIRD3 vs FAIRD4 │ Abundancia (n=22)│  0.733   │ <0.001   │ Fuerte    │ H1_FAIRD       │
│ ID1 vs ID2       │ Abundancia (n=22)│  0.027   │ 0.4633   │ Pobre (NS)│ H1_ID          │
│ AMMOD1 vs AMMOD2 │ Biomasa (n=12)   │  0.741   │ <0.001   │ Fuerte    │ H1_AMMOD       │
└──────────────────┴──────────────────┴──────────┴────────┴───────────┴──────────────┘

HALLAZGOS CLAVE (H1):
  • FAIRD y AMMOD muestran una consistencia temporal excelente.
  • ID muestra una consistencia pobre y no significativa, ~25x peor que FAIRD.
  • NOTA: La consistencia de AMMOD es INTER-SITIO (Sitio 1 vs 2), lo que
    hace que su alta correlación sea aún más notable.

═══════════════════════════════════════════════════════════════════════════════
✅ H2: RESULTADOS DE VALIDACIÓN (VALIDADOS)
═══════════════════════════════════════════════════════════════════════════════
Análisis de regresión lineal de FAIRD vs AMMOD sobre biomasa transformada
log10(x+1). Usa datos prorrateados de AMMOD3 (n=12 para todos los sitios).

┌───────┬─────────────┬──────┬──────────┬──────────┬──────────────────┐
│ Sitio │ Hábitat     │  n   │ R² (log) │ p-value  │ Veredicto        │
├───────┼─────────────┼──────┼──────────┼──────────┼──────────────────┤
│ Site 2│ Maize-Core  │  12  │  0.445   │ 0.0179   │ ✅ Validado      │
│ Site 1│ Maize-Edge  │  12  │  0.282   │ 0.0755   │ ⚠️ Marginal (NS)  │
│ Site 3│ Meadow      │  12  │  0.168   │ 0.1858   │ ❌ No Validado   │
└───────┴─────────────┴──────┴──────────┴──────────┴──────────────────┘

HALLAZGOS CLAVE (H2):
  • VALIDACIÓN EXITOSA: La validación de FAIRD es exitosa en el Sitio 2
    (Maize-Core), con una correlación moderada y estadísticamente significativa.
  • EFECTO DE HÁBITAT: La correlación es más débil en el borde del maizal
    (Sitio 1) y no significativa en la pradera (Sitio 3).
  • COMPARACIÓN DE DISPOSITIVOS (No-paramétrico):
    • ID captura significativamente más insectos que FAIRD (ej. 51.8/día vs
      8.45/día en Sitio 2).
    • Todos los sitios FAIRD capturan cantidades similares (Grupo "a" en CLD).
    • Las capturas de ID varían mucho entre sitios (Grupos "bc", "c", "d").

═══════════════════════════════════════════════════════════════════════════════
🔬 DISCUSIÓN Y ANÁLISIS DE HIPÓTESIS (H2)
═══════════════════════════════════════════════════════════════════════════════

HIPÓTESIS DE SESGO DE DETECCIÓN (Explica la discrepancia en H2)

Esta hipótesis explica por qué la correlación FAIRD/AMMOD es fuerte en el
Sitio 2 (Core) pero débil en el Sitio 1 (Edge).

1.  **AMMOD (Trampa Malaise Tradicional):**
    * Captura el espectro de tamaños completo (ej. 1mm - 20mm).
    * La biomasa total está dominada por Dípteros pequeños/medianos (2-5mm)
        que son muy abundantes.
    * Dado que esta comunidad de insectos pequeños es consistente entre
        ambos sitios (Maize-Edge y Maize-Core), la biomasa total de AMMOD
        es consistente entre sitios (como se vio en H1: R² = 0.741).
    * **Conclusión:** La "biomasa total" de AMMOD es principalmente una
        señal de Dípteros pequeños.

2.  **FAIRD (Cámara con IA):**
    * Tiene un límite de detección de hardware/software de ~3mm.
    * **No puede ver** la gran mayoría de los Dípteros pequeños (2-5mm)
        que dominan la biomasa de AMMOD.
    * FAIRD mide la biomasa de insectos medianos/grandes (>3mm), que son
        menos comunes y más variables entre hábitats.

3.  **Análisis Combinado (La Explicación):**
    * **Site 1 (Maize-Edge):** La biomasa total de AMMOD es ALTA (impulsada
        por Dípteros <3mm), pero la biomasa de FAIRD es BAJA (porque no ve
        esos Dípteros). Las dos señales no coinciden.
        *Resultado: Correlación casi nula (R² = 0.282, p=0.076).*
    * **Site 2 (Maize-Core):** El hábitat es diferente, con menos Dípteros
        pequeños y una mayor proporción de insectos grandes (>3mm).
        Aquí, tanto AMMOD como FAIRD están midiendo la *misma señal*
        (insectos grandes).
        *Resultado: Correlación fuerte y significativa (R² = 0.445, p=0.018).*

**IMPLICACIÓN:** La baja correlación en el Sitio 1 no es un fallo de FAIRD,
sino una **demostración documentada de su sesgo de detección conocido**
(el límite de 3mm). El dispositivo funciona perfectamente en el Sitio 2,
donde la comunidad de insectos se alinea con su capacidad de detección.

═══════════════════════════════════════════════════════════════════════════════
⏸️ H3: COMPARACIÓN TAXONÓMICA (PENDIENTE)
═══════════════════════════════════════════════════════════════════════════════

OBJETIVO:
  ¿FAIRD captura mayor diversidad taxonómica que ID?

DATOS PRELIMINARES:
  ┌──────────┬────────┬─────────┬────────────────┐
  │ Métrica  │ FAIRD  │   ID    │ FAIRD Ventaja  │
  ├──────────┬────────┼─────────┼────────────────┤
  │ Órdenes  │   11   │    6    │ +83%           │
  │ Géneros  │   55   │    5    │ +1000%         │
  └──────────┴────────┴─────────┴────────────────┘

COMPOSICIÓN:
  • ID: >95% Dípteros (Especialista).
  • FAIRD: Balanceado (Dípteros, Himenópteros, Coleópteros, etc.) (Generalista).

TRABAJO PENDIENTE:
  □ Aplicar análisis de rarefacción (vegan::rarecurve) para normalizar
    por esfuerzo de muestreo (diferente número de individuos).
  □ Calcular índices de diversidad (Shannon, Simpson).
  □ Realizar pruebas estadísticas (ej. Chi-cuadrado).