# R Namespace Standards — MANDATORY

## Why

This project loads `MASS` after `tidyverse`, which masks `dplyr::select()` and `dplyr::filter()`. Explicit namespaces prevent silent failures.

## Required namespaces

### dplyr (ALWAYS prefix)
```r
dplyr::filter()    dplyr::select()     dplyr::mutate()      dplyr::summarise()
dplyr::arrange()   dplyr::group_by()   dplyr::ungroup()     dplyr::rename()
dplyr::distinct()  dplyr::count()      dplyr::pull()        dplyr::slice()
dplyr::case_when() dplyr::if_else()    dplyr::n()           dplyr::lag()
dplyr::bind_rows() dplyr::bind_cols()  dplyr::across()      dplyr::n_distinct()
dplyr::left_join() dplyr::inner_join() dplyr::right_join()  dplyr::full_join()
dplyr::semi_join() dplyr::anti_join()
```

### tidyr (ALWAYS prefix)
```r
tidyr::pivot_wider()  tidyr::pivot_longer()  tidyr::drop_na()
tidyr::complete()     tidyr::separate()      tidyr::unite()
tidyr::replace_na()   tidyr::nest()          tidyr::unnest()
```

### readr (ALWAYS prefix)
```r
readr::read_csv()  readr::write_csv()
```

### tibble (when used explicitly)
```r
tibble::tibble()  tibble::as_tibble()  tibble::column_to_rownames()  tibble::rownames_to_column()
```

## Package loading order

```r
library(tidyverse)   # 1st — includes dplyr, tidyr, ggplot2, readr, stringr
library(patchwork)   # Safe — no conflicts
library(scales)      # Safe
library(car)         # 2nd — ncvTest(), vif()
library(lmtest)      # 3rd — dwtest()
library(MASS)        # LAST — masks dplyr::select(), dplyr::filter()
library(vegan)       # After MASS — vegdist(), diversity()
```

## ggplot2 — no prefix needed

`ggplot()`, `aes()`, `geom_*()`, `scale_*()`, `theme_*()`, `labs()`, `ggsave()` do not conflict. No namespace required.

## Base R — no prefix needed

`mean()`, `sd()`, `sum()`, `round()`, `paste()`, `sprintf()`, `cat()`, `print()`, `cor.test()`, `lm()`, `shapiro.test()`, `wilcox.test()`, `t.test()`, `setdiff()`, `intersect()`, `unique()`, `nrow()`, `ncol()`, `seq_len()` — all safe without prefix.

## Common error signatures

```
Error in select(...) : unused arguments     → Use dplyr::select()
Error in filter(...) : object not found     → Use dplyr::filter()
could not find function "n"                 → Use dplyr::n() inside dplyr::summarise()
```
