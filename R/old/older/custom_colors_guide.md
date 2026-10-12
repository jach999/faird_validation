# Custom Colors Guide
## Project: FAIRD Field Validation — ggplot2 Color Usage

---

## The Core Rule

**Always pass `custom_colors` directly to `scale_color_manual(values = ...)`.  
Never construct intermediate named vectors by subsetting `custom_colors`.**

ggplot2 automatically selects only the color names that appear in the data.  
No pre-filtering needed.

---

## ✅ CORRECT patterns

```r
# Pattern 1: Device type as aesthetic (FAIRD / ID / AMMOD)
scale_color_manual(values = custom_colors, name = "Device Type")

# Pattern 2: Site as aesthetic (Site1 / Site2 / Site3 / Site4)
scale_color_manual(values = custom_colors, name = "Site")

# Pattern 3: Habitat as aesthetic (Maize / Meadow)
scale_color_manual(values = custom_colors, name = "Habitat")

# Pattern 4: fill variant (e.g., violin, boxplot)
scale_fill_manual(values = custom_colors, name = "Device Type")
```

---

## ❌ WRONG patterns — cause gray/invisible points

```r
# WRONG — subsetting creates compound names (e.g., "FAIRD.FAIRD")
scale_color_manual(
  values = c(FAIRD = custom_colors["FAIRD"],
             AMMOD = custom_colors["AMMOD"])
)

# WRONG — same problem with site colors
site_palette <- c(
  Site1 = custom_colors["Site1"],
  Site2 = custom_colors["Site2"]
)
scale_color_manual(values = site_palette)

# WHY IT FAILS:
# custom_colors["FAIRD"] returns a named vector: c(FAIRD = "springgreen4")
# c(FAIRD = c(FAIRD = "springgreen4")) → name becomes "FAIRD.FAIRD"
# ggplot looks for "FAIRD" but finds only "FAIRD.FAIRD" → no match → gray
```

---

## Custom Colors Reference

```r
# Device types (used when color = Device_type or color = Method)
"FAIRD"  = "springgreen4"
"ID"     = "coral2"
"AMMOD"  = "dodgerblue"

# Sites (used when color = Site)
"Site1"  = "gold"
"Site2"  = "darkorange"
"Site3"  = "greenyellow"
"Site4"  = "darkgreen"

# Habitats (used when color = Habitat or color = Environment)
"Maize"  = "gold"
"Meadow" = "darkgreen"

# Individual devices (used when color = Device)
"FAIRD1" = "gold"       "ID1" = "gold"       "AMMOD1" = "gold"
"FAIRD2" = "darkorange" "ID2" = "darkorange" "AMMOD2" = "darkorange"
"FAIRD3" = "greenyellow""ID3" = "greenyellow""AMMOD3" = "greenyellow"
"FAIRD4" = "darkgreen"  "ID4" = "darkgreen"  "AMMOD4" = "darkgreen"

# Weather variables
"AIRTEMP"       = "red"
"WINDSPEED"     = "cyan"
"PRECIPITATION" = "lightskyblue"
```

---

## Real Examples from Project Scripts

### H2_validation_FAIRD_AMMOD.R (reference implementation)

```r
# Temporal plot — Method levels are "FAIRD" and "AMMOD"
ggplot(site_temporal, aes(x = Date, y = Biomass, color = Method)) +
  geom_line() +
  scale_color_manual(values = custom_colors, name = "Device Type")

# Scatter plot — color by site, full palette passed
geom_smooth(fill = site_color, ...)  # site_color extracted separately for fill only
```

### Extracting a single color (for non-scale uses)

```r
# OK to extract for direct use in geom arguments (not in scale_color_manual)
site_color  <- custom_colors[site]       # e.g., for geom_smooth(fill = site_color)
faird_color <- custom_colors["FAIRD"]   # e.g., for annotate() or geom_hline()

# Just don't put these extracted values back into c() with names
```

---

## Checklist — Before finalizing any ggplot with custom colors

- [ ] `scale_color_manual(values = custom_colors)` — full vector, no subsetting
- [ ] `scale_fill_manual(values = custom_colors)` — same rule for fills
- [ ] The aesthetic variable name in `aes(color = ...)` matches a key in `custom_colors`
      (e.g., if using `color = Method`, ensure Method values are "FAIRD"/"ID"/"AMMOD")
- [ ] `source("custom_colors.R")` called before any ggplot code

---

*Last updated: 2026-03-17*
