# 1. CORAL COMMUNITY ANALYSIS----------------------------------------------------
### Load packages ----------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(vegan)
library(glmmTMB)
library(DHARMa)
library(performance)
library(lme4)
library(lmerTest)
library(ggrepel)
library(patchwork)
library(showtext)
library(tidyverse)
library(sjPlot)
library(brms)
library(bayesplot)
library(tidybayes)
library(ggridges)
library(gridExtra)
library(cowplot)
library(performance)
library(knitr)
library(performance)
library(flextable)
library(broom.mixed)
library(dplyr)
library(officer)
library(readxl)
library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)
library(flextable)
library(officer)
library(ggrepel)
library(patchwork)
### Aesthetics -------------------------------------------------------------------

font_add(
        family = "Montserrat",
        regular = "D:/Descargas/Montserrat/static/Montserrat-Regular.ttf",
        bold = "D:/Descargas/Montserrat/static/Montserrat-Bold.ttf",
        italic = "D:/Descargas/Montserrat/static/Montserrat-Italic.ttf",
        bolditalic = "D:/Descargas/Montserrat/static/Montserrat-BoldItalic.ttf"
)
showtext_auto()

# 1.1 Colony Density Analysis------------------------
### Load and prepare data --------------------------------------------------------
df_com_rel <- read_excel(
        "D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/Corals_data_JR_Forereefs_Terrace_2024_August_May.xlsx",
        sheet = "Commu_sp_rel"
) %>%
        filter(Site != "Auras") %>%
        droplevels() %>%
        mutate(
                Month = factor(Month, levels = c("May", "August")),
                Habitat = factor(Habitat, levels = c("Fore reef", "Terrace")),
                Site = factor(Site, levels = c("Peruano", "Pinos", "Anclitas", 
                                               "Mariflores", "Cruces")),
                Transept = factor(Transept)
        )

### Data exploration
glimpse(df_com_rel)
table(df_com_rel$Site, df_com_rel$Habitat)
View(df_com_rel)
### Data exploration -------------------------------------------------------------
density_stats_habitat <- df_com_rel %>%
        group_by(Habitat) %>%
        summarise(
                Mean = mean(Density, na.rm = TRUE),
                SD = sd(Density, na.rm = TRUE),
                N = n(),
                .groups = 'drop'
        )

print(density_stats_habitat)

# Habitat    Mean    SD     N
# Fore reef  81.4  21.7   100
# Terrace    72    26.8   100

density_stats_month <- df_com_rel %>%
        group_by(Month) %>%
        summarise(
                Mean = mean(Density, na.rm = TRUE),
                SD = sd(Density, na.rm = TRUE),
                N = n(),
                .groups = 'drop'
        )

print(density_stats_month)

# Month   Mean    SD     N
# May     81.9  25.4   100
# August  71.5  23.1   100

density_stats_hm <- df_com_rel %>%
        group_by(Month, Habitat) %>%
        summarise(
                Mean = mean(Density, na.rm = TRUE),
                SD = sd(Density, na.rm = TRUE),
                N = n(),
                .groups = 'drop'
        )

print(density_stats_hm)

# Month     Habitat     Mean  SD    N
# May    x  Fore reef   87.2  21.7  50
# May    x  Terrace     76.5  27.8  50
# August x  Fore reef   75.5  20.3  50
# August x  Terrace     67.5  25.1  50

mean(df_com_rel$Density)  # 76.685
var(df_com_rel$Density)   # 612.96 (Razon varianza media es 8, lo que indica sobredispersion)

# Visualize distributions
ggplot(df_com_rel, aes(x = Density)) +
        geom_histogram(bins = 15) +
        facet_grid(Month ~ Habitat) +
        theme_bw() +
        labs(title = "Density distribution by Habitat and Month")

ggplot(df_com_rel, aes(x = Month, y = Density, fill = Month)) +
        geom_violin(outlier.shape = 21, alpha = 0.3) +
        geom_boxplot(outlier.shape = 21, alpha = 0.5) +
        geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
        facet_wrap(~ Habitat) +
        theme_minimal() +
        labs(
                title = "Colony density by Habitat and Month",
                x = "Month",
                y = expression("Colony density (col / 10 m"^2*")")
        ) +
        theme(legend.position = "none")

ggplot(df_com_rel, aes(sample = Density)) +
        stat_qq() + stat_qq_line() +
        facet_grid(Month ~ Habitat) +
        theme_bw()

# Shapiro test by group. Density is normal 
df_com_rel %>%
        group_by(Month, Habitat) %>%
        summarise(p_value = shapiro.test(Density)$p.value, n = n(), 
                  .groups = "drop") 

### Models Full model with interaction vs additive -----------------------------

# Model 1: Interaction Mixed 
###
model_density_inter_mix <- glmmTMB(
        Density ~ Habitat * Month + (1|Site),
        family = nbinom2,
        data = df_com_rel
)
summary(model_density_inter_mix)
r2(model_density_inter_mix)
# Model 2: Additive Mixed
model_density_add_mix <- glmmTMB(
        Density ~ Habitat + Month + (1|Site),
        family = nbinom2,
        data = df_com_rel
)
summary(model_density_add_mix)
r2(model_density_add_mix)
# Model 3: Additive
model_density_add <- glmmTMB(
        Density ~ Habitat + Month,
        family = nbinom2,
        data = df_com_rel
)
summary(model_density_add)
r2(model_density_add)

model_density_all_mix <- glmmTMB(
        Density ~ Habitat + Month + Habitat * Month + (1|Site),
        family = nbinom2,
        data = df_com_rel
)
summary(model_density_all_mix)
r2(model_density_add_mix)

#Prueba de razón de verosimilitud (Likelihood Ratio Test, LRT)
lrt <-anova(model_density_inter_mix, model_density_add_mix, model_density_add) 

# Crear tabla con los parámetros clave
tabla_lrt <- data.frame(
        Modelo = c("model_density_add", "model_density_add_mix", "model_density_inter_mix"),
        Formula = c(
                deparse(formula(model_density_add)),
                deparse(formula(model_density_add_mix)),
                deparse(formula(model_density_inter_mix))
        ),
        AIC     = lrt$AIC,
        BIC     = lrt$BIC,
        logLik  = as.numeric(lrt$logLik),
        Chi2    = lrt$Chisq,
        gl      = lrt$`Chi Df`,
        p_valor = lrt$`Pr(>Chisq)`,
        Signif  = symnum(
                lrt$`Pr(>Chisq)`,
                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                symbols   = c("***", "**", "*", ".", "ns"),
                na        = ""
        ) |> as.character()
)

# Tabla
knitr::kable(
        tabla_lrt,
        format    = "pipe",
        digits    = 3,
        col.names = c("Modelo", "Fórmula", "AIC", "BIC", "logLik", "χ²", "gl", "p-valor", "Sig."),
        caption   = "Comparación de modelos — Likelihood Ratio Test"
)

# R² and variance partitioning 
r2(model_density_add_mix)
# Conditional R²: 0.291 (29% total variance explained)
# Marginal R²: 0.086 (8.6% by fixed effects)

# Model validation (DHARMa) 
sim_res <- simulateResiduals(model_density_add_mix, n = 1000)
plot(sim_res)
testDispersion(sim_res) # dispersión = 0.860, p = 0.618 → Sin sobre ni subdispersión, la varianza está bien capturada. 
testUniformity(sim_res) # D = 0.061, p = 0.447 → Los residuos se distribuyen uniformemente, el modelo ajusta bien.
testOutliers(sim_res) # 0 outliers en 200 observaciones, p = 1.000 → Sin valores extremos problemáticos.
check_overdispersion(model_density_add_mix) # ratio = 0.847, p = 0.568 → No se detecta sobredispersión.

### Observed vs Predicted 
df_com_rel$pred_final <- predict(model_density_add_mix, type = "response")
cor(df_com_rel$Density, df_com_rel$pred_final) # r = 0.527, ajuste moderado.

### RESULT. TABLE 1. Model parameters table. -----------------------------------------------------

# Extraer resultados del modelo
model_results <- broom.mixed::tidy( model_density_add_mix)
model_results

# Preparar efectos fijos
tabla_fixed <- model_results %>%
        filter(effect == "fixed") %>%
        mutate(
                Parameter = case_when(
                        term == "(Intercept)"     ~ "Intercept",
                        term == "HabitatTerrace"  ~ "Habitat (Terrace)",
                        term == "MonthAugust"     ~ "Month (August)"
                ),
                Estimate  = as.character(round(estimate, 3)),
                SE        = as.character(round(std.error, 3)),
                `z-value` = as.character(round(statistic, 2)),
                `p-value` = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
        ) %>%
        select(Parameter, Estimate, SE, `z-value`, `p-value`)
tabla_fixed
# Preparar efecto aleatorio
tabla_random <- model_results %>%
        filter(effect == "ran_pars") %>%
        mutate(
                Parameter = "Site (SD)",
                Estimate  = as.character(round(estimate, 3)),
                SE        = "—",
                `z-value` = "—",
                `p-value` = "—"
        ) %>%
        select(Parameter, Estimate, SE, `z-value`, `p-value`)

# Unir ambas partes
tabla_df <- bind_rows(tabla_fixed, tabla_random)
tabla_df
# Crear flextable
ft_complete_density <- flextable(tabla_df) %>%
        theme_booktabs() %>%
        autofit() %>%
        #bold(i = 1, part = "body") %>%
        hline(
                i      = nrow(tabla_fixed),
                border = fp_border(color = "gray50", width = 0.75)
        ) %>%
        align(
                j     = c("Estimate", "SE", "z-value", "p-value"),
                align = "center",
                part  = "all"
        ) %>%
        width(width = 6.5 / 6) %>% 
        # Añadir título como header
        add_header_lines(
                values = "Table 1. Parameter estimates from the additive generalized linear mixed-effects model (GLMM; negative binomial error distribution with log link) evaluating the effects of Habitat (Terrace vs. Fore reef) and Month (May vs. August) on colony density. Site was included as a random intercept. R\u00b2 marginal = 0.086, R\u00b2 conditional = 0.291, AIC = 1828.9"
        ) %>%
        # Formato del título
        align(i = 1, align = "left", part = "header")

ft_complete_density
save_as_docx(ft_complete_density, path = "Table_1_Community_Density_Model.docx")

# 1.2 Community Structure Analysis Relative abundance---------------------------------------------

### Load and prepare data --------------------------------------------------------
df_com_rel <- read_excel(
        "D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/Corals_data_JR_Forereefs_Terrace_2024_August_May.xlsx",
        sheet = "Commu_sp_rel"
) %>%
        filter(Site != "Auras") %>%
        droplevels() %>%
        mutate(
                Month = factor(Month, levels = c("May", "August")),
                Habitat = factor(Habitat, levels = c("Fore reef", "Terrace")),
                Site = factor(Site, levels = c("Peruano", "Pinos", "Anclitas", 
                                               "Mariflores", "Cruces")),
                Transept = factor(Transept)
        )
df_com_rel 
### Extract community matrix and metadata
comm_matrix <- df_com_rel %>%
        select(Agaaga:Steint) %>%
        as.matrix()

metadata <- df_com_rel %>%
        select(Month, Habitat, Site, Transept)

# Data summary
dim(comm_matrix)  # 200 samples × 37 species
glimpse(comm_matrix)
View(comm_matrix)
### RESULT. TABLE 2.A. Relative abundance per species for Habitat and Month--------
# Diccionario de nombres científicos
species_names <- c(
        "Sidsid" = "Siderastrea siderea (Ellis & Solander, 1786)",
        "Steint" = "Stephanocoenia intersepta (Esper, 1795)",
        "Dicsto" = "Dichocoenia stokesii Milne Edwards & Haime, 1849",
        "Orbfav" = "Orbicella faveolata (Ellis & Solander, 1786)",
        "Psestr" = "Pseudodiploria strigosa (Dana, 1846)",
        "Moncav" = "Montastraea cavernosa (Linnaeus, 1767)",
        "Porast" = "Porites astreoides Lamarck, 1816",
        "Orbann" = "Orbicella annularis (Ellis & Solander, 1786)",
        "Diplab" = "Diploria labyrinthiformis (Linnaeus, 1758)",
        "Sidrad" = "Siderastrea radians (Pallas, 1766)",
        "Meamea" = "Meandrina meandrites (Linnaeus, 1758)",
        "Porfur" = "Porites furcata Lamarck, 1816",
        "Porpor" = "Porites porites (Pallas, 1766)",
        "Pordiv" = "Porites divaricata Le Sueur, 1820",
        "Maddec" = "Madracis decactis (Lyman, 1859)",
        "Eusfas" = "Eusmilia fastigiata (Pallas, 1766)",
        "Orbfra" = "Orbicella franksi (Gregory, 1895)",
        "Agalam" = "Agaricia lamarcki Milne Edwards & Haime, 1851",
        "Favfra" = "Favia fragum (Esper, 1793)",
        "Scocub" = "Scolymia cubensis (Milne Edwards & Haime, 1848)",
        "Agaaga" = "Agaricia agaricites (Linnaeus, 1758)",
        "Agaund" = "Agaricia undata (Ellis & Solander, 1786)",
        "Madfor" = "Madracis formosa Wells, 1973",
        "Colnat" = "Colpophyllia natans (Houttuyn, 1772)",
        "Dencyl" = "Dendrogyra cylindrus (Ehrenberg, 1834)",
        "Milalc" = "Millepora alcicornis Linnaeus, 1758",
        "Psecli" = "Pseudodiploria clivosa (Ellis & Solander, 1786)",
        "Musang" = "Mussa angulosa (Pallas, 1766)",
        "Agafra" = "Agaricia fragilis Dana, 1846",
        "Meajac" = "Meandrina jacksoni Weil & Pinzón, 2011",
        "Mycali" = "Mycetophyllia aliciae Wells, 1973",
        "Mycfer" = "Mycetophyllia ferox Wells, 1973",
        "Myclam" = "Mycetophyllia lamarckiana Milne Edwards & Haime, 1849",
        "Agahum" = "Agaricia humilis Verrill, 1901",
        "Madmir" = "Madracis mirabilis (Duchassaing & Michelotti, 1860)",
        "Mycdan" = "Mycetophyllia danaana Milne Edwards & Haime, 1849",
        "Scolac" = "Scolymia lacera (Pallas, 1766)",
        "Lepcuc" = "Helioseris cucullata (Ellis & Solander, 1786)"
        )

species_names <- c(
        "Sidsid" = "Siderastrea siderea",
        "Steint" = "Stephanocoenia intersepta",
        "Dicsto" = "Dichocoenia stokesii",
        "Orbfav" = "Orbicella faveolata",
        "Psestr" = "Pseudodiploria strigosa",
        "Moncav" = "Montastraea cavernosa",
        "Porast" = "Porites astreoides",
        "Orbann" = "Orbicella annularis",
        "Diplab" = "Diploria labyrinthiformis",
        "Sidrad" = "Siderastrea radians",
        "Meamea" = "Meandrina meandrites",
        "Porfur" = "Porites furcata",
        "Porpor" = "Porites porites",
        "Pordiv" = "Porites divaricata",
        "Maddec" = "Madracis decactis",
        "Eusfas" = "Eusmilia fastigiata",
        "Orbfra" = "Orbicella franksi",
        "Agalam" = "Agaricia lamarcki",
        "Favfra" = "Favia fragum",
        "Scocub" = "Scolymia cubensis",
        "Agaaga" = "Agaricia agaricites",
        "Agaund" = "Agaricia undata",
        "Madfor" = "Madracis formosa",
        "Dencyl" = "Dendrogyra cylindrus",
        "Colnat" = "Colpophyllia natans",
        "Milalc" = "Millepora alcicornis",
        "Psecli" = "Pseudodiploria clivosa",
        "Musang" = "Mussa angulosa",
        "Agafra" = "Agaricia fragilis",
        "Meajac" = "Meandrina jacksoni",
        "Mycali" = "Mycetophyllia aliciae",
        "Mycfer" = "Mycetophyllia ferox",
        "Myclam" = "Mycetophyllia lamarckiana",
        "Agahum" = "Agaricia humilis",
        "Madmir" = "Madracis mirabilis",
        "Mycdan" = "Mycetophyllia danaana",
        "Scolac" = "Scolymia lacera",
        "Lepcuc" = "Helioseris cucullata"
        )


# 1. Calcular media ± SD GLOBAL por especie (sin agrupar por Habitat/Month)
abundancia_global <- df_com_rel %>%
        pivot_longer(cols = Agaaga:Steint,
                     names_to = "Species_Code",
                     values_to = "Abundance") %>%
        group_by(Species_Code) %>%
        summarise(
                Mean = mean(Abundance, na.rm = TRUE),
                SD = sd(Abundance, na.rm = TRUE),
                .groups = "drop"
        ) %>%
        mutate(
                Label = sprintf("%.2f ± %.2f", Mean, SD),
                Species = species_names[Species_Code]
        ) %>%
        arrange(desc(Mean)) %>%
        select(Species, Label) %>%
        rename(Overall = Label)

abundancia_global

# 2. Calcular media ± SD GLOBAL por especie agrupando por Habitat/Month)
abundancia_por_grupo <- df_com_rel %>%
        pivot_longer(cols = Agaaga:Steint,
                     names_to = "Species_Code",
                     values_to = "Abundance") %>%
        group_by(Species_Code, Habitat, Month) %>%
        summarise(
                Mean = mean(Abundance, na.rm = TRUE),
                SD = sd(Abundance, na.rm = TRUE),
                N = n(),
                .groups = "drop"
        ) %>%
        mutate(
                Label = sprintf("%.2f ± %.2f", Mean, SD),
                Species = species_names[Species_Code]
        ) %>%
        select(Species, Habitat, Month, Label)
abundancia_por_grupo
# 3. Pivotar a formato ancho

abundancia_tabla <- abundancia_global %>%
        left_join(
                abundancia_por_grupo %>%
                        pivot_wider(
                                names_from = c(Habitat, Month),
                                values_from = Label,
                                names_sep = "_"
                        ),
                by = "Species"
        ) %>%
        select(Species, Overall, everything())

abundancia_tabla

# 4. CREAR Tabla 2a de abundancias relativas

# Ancho disponible en página Word con márgenes de 2 cm
# Carta (21.59 cm) - margen izq (2 cm) - margen der (2 cm) = 17.59 cm
page_width_cm <- 21.59 - 2 - 2  # = 17.59 cm
page_width_in <- page_width_cm / 2.54  # convertir a pulgadas (flextable usa pulgadas)
page_width_in

ft_abundancia <- flextable(abundancia_tabla) %>%
        set_header_labels(
                Species = "Species",
                Overall = "Overall\n(%)",
                `Fore reef_May` = "Fore reef\nMay (%)",
                `Fore reef_August` = "Fore reef\nAugust (%)",
                `Terrace_May` = "Terrace\nMay (%)",
                `Terrace_August` = "Terrace\nAugust (%)"
        ) %>%
        theme_booktabs() %>%
        font(fontname = "Times New Roman", part = "all") %>%        
        fontsize(size = 10, part = "all") %>%                       
        compose(
                j = "Species",
                part = "body",
                value = as_paragraph(
                        as_chunk(
                                gsub("^(\\S+\\s+\\S+).*", "\\1", Species),
                                props = fp_text(italic = TRUE, 
                                                font.size = 10,
                                                font.family = "Times New Roman")  # <-- agrega font.family
                        ),
                        as_chunk(
                                ifelse(
                                        str_detect(Species, "^\\S+\\s+\\S+\\s+"),
                                        paste0(" ", gsub("^\\S+\\s+\\S+\\s+(.*)", "\\1", Species)),
                                        ""
                                ),
                                props = fp_text(italic = FALSE, 
                                                font.size = 10,
                                                font.family = "Times New Roman")  # <-- agrega font.family
                        )
                )
        ) %>%
        align(j = 1, align = "left", part = "all") %>%
        align(j = 2:6, align = "center", part = "all") %>%
        padding(j = 1:6, padding.left = 3, padding.right = 3, part = "all") %>%
        padding(j = 1:6, padding.top = 2, padding.bottom = 2, part = "all") %>%
        add_header_lines(
                values = "Table 2.a Mean relative abundance (%) ± standard deviation of coral species in each habitat and month combination, with overall average across all samples."
        ) %>%
        align(i = 1, align = "left", part = "header") %>%
        #autofit() %>%
        width(j = 1:6, width = 2.5) %>%      # <-- ancho col 1 en pulgadas
        width(j = 2:6, width = (page_width_in-2.5)/5 ) %>%
        fit_to_width(max_width = page_width_in)
ft_abundancia

save_as_docx(ft_abundancia, path = "Table_2.a_Community_Species_Abundance_HabitatMonth.docx")

# Calcular la media de cada especie y ordenar de mayor a menor
media_especies <- df_com_rel %>%
        select(-Month, -Habitat, -Site, -Transept, -Density) %>%
        summarise(across(everything(), list(
                Media = ~mean(., na.rm = TRUE),
                SD = ~sd(., na.rm = TRUE)
        ))) %>%
        pivot_longer(cols = everything(), 
                     names_to = c("Especie", ".value"), 
                     names_pattern = "(.+)_(.+)") %>%
        arrange(desc(Media))

# Ver los resultados
media_especies
# Ver los resultados
View (media_especies)


# 1. Media global (ya la tienes)
media_global <- df_com_rel %>%
        select(-Month, -Habitat, -Site, -Transept, -Density) %>%
        summarise(across(everything(), list(
                Media = ~mean(., na.rm = TRUE),
                SD    = ~sd(., na.rm = TRUE)
        ))) %>%
        pivot_longer(cols = everything(),
                     names_to = c("Especie", ".value"),
                     names_pattern = "(.+)_(.+)") %>%
        arrange(desc(Media))
View(media_global)
# 2. Media por Habitat x Month
media_habitat_month <- df_com_rel %>%
        select(-Site, -Transept, -Density) %>%
        group_by(Habitat, Month) %>%
        summarise(across(everything(), list(
                Media = ~mean(., na.rm = TRUE),
                SD    = ~sd(., na.rm = TRUE)
        )), .groups = "drop") %>%
        pivot_longer(cols = -c(Habitat, Month),
                     names_to = c("Especie", ".value"),
                     names_pattern = "(.+)_(.+)") %>%
        mutate(Grupo = paste0(Habitat, "_", Month)) %>%
        select(Especie, Grupo, Media, SD) %>%
        pivot_wider(names_from = Grupo,
                    values_from = c(Media, SD),
                    names_glue = "{Grupo}_{.value}")
View(media_habitat_month)
# 3. Media por Habitat x Site
media_habitat_site <- df_com_rel %>%
        select(-Month, -Transept, -Density) %>%
        group_by(Habitat, Site) %>%
        summarise(across(everything(), list(
                Media = ~mean(., na.rm = TRUE),
                SD    = ~sd(., na.rm = TRUE)
        )), .groups = "drop") %>%
        pivot_longer(cols = -c(Habitat, Site),
                     names_to = c("Especie", ".value"),
                     names_pattern = "(.+)_(.+)") %>%
        mutate(Grupo = paste0(Habitat, "_", Site)) %>%
        select(Especie, Grupo, Media, SD) %>%
        pivot_wider(names_from = Grupo,
                    values_from = c(Media, SD),
                    names_glue = "{Grupo}_{.value}")
View(media_habitat_site)
# 4. Unir todo
media_completa <- media_global %>%
        left_join(media_habitat_month, by = "Especie") %>%
        left_join(media_habitat_site,  by = "Especie")
media_completa
View(media_completa)
media_completa %>%
        slice(1:10) %>%
        select(1:3, 12:ncol(.)) %>%
        print(n = Inf, width = Inf)

media_completa %>%
        select(Especie, Media, SD, contains("Cruces")) %>%
        print(n = Inf, width = Inf)

## 1.3 Community Structure Analysis Density---------------------------------------------
### Load packages 

### Load and prepare data --------------------------------------------------------
df_com <- read_excel(
        "D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/Corals_data_JR_Forereefs_Terrace_2024_August_May.xlsx",
        sheet = "Commu_sp"
) %>%
        filter(Site != "Auras") %>%
        droplevels() %>%
        mutate(
                Month = factor(Month, levels = c("May", "August")),
                Habitat = factor(Habitat, levels = c("Fore reef", "Terrace")),
                Site = factor(Site, levels = c("Peruano", "Pinos", "Anclitas", 
                                               "Mariflores", "Cruces")),
                Transept = factor(Transept)
        )
df_com 

df_com <- df_com %>%
        mutate(across(where(is.numeric), ~replace_na(., 0)))
### Extract community matrix and metadata
comm_matrix <- df_com %>%
        select(Agaaga:Steint) %>%
        as.matrix()

metadata <- df_com_rel %>%
        select(Month, Habitat, Site, Transept)

# Data summary
dim(comm_matrix)  # 200 samples × 37 species
glimpse(comm_matrix)
comm_matrix

### RESULT. TABLE 2.B Colony density per species for Habitat and Month--------

# 1. Calcular media ± SD GLOBAL por especie (sin agrupar por Habitat/Month)
abundancia_global <- df_com %>%
        pivot_longer(cols = Agaaga:Steint,
                     names_to = "Species_Code",
                     values_to = "Abundance") %>%
        group_by(Species_Code) %>%
        summarise(
                Mean = mean(Abundance, na.rm = TRUE),
                SD = sd(Abundance, na.rm = TRUE),
                .groups = "drop"
        ) %>%
        mutate(
                Label = sprintf("%.2f ± %.2f", Mean, SD),
                Species = species_names[Species_Code]
        ) %>%
        arrange(desc(Mean)) %>%
        select(Species, Label) %>%
        rename(Overall = Label)
abundancia_global
# 2. Crear la tabla con datos por grupo (como antes)
abundancia_por_grupo <- df_com %>%
        pivot_longer(cols = Agaaga:Steint,
                     names_to = "Species_Code",
                     values_to = "Abundance") %>%
        group_by(Species_Code, Habitat, Month) %>%
        summarise(
                Mean = mean(Abundance, na.rm = TRUE),
                SD = sd(Abundance, na.rm = TRUE),
                N = n(),
                .groups = "drop"
        ) %>%
        mutate(
                Label = sprintf("%.2f ± %.2f", Mean, SD),
                Species = species_names[Species_Code]
        ) %>%
        select(Species, Habitat, Month, Label)

# 3. Pivotar a formato ancho

abundancia_tabla <- abundancia_global %>%
        left_join(
                abundancia_por_grupo %>%
                        pivot_wider(
                                names_from = c(Habitat, Month),
                                values_from = Label,
                                names_sep = "_"
                        ),
                by = "Species"
        ) %>%
        select(Species, Overall, everything())

abundancia_tabla

# 4. CREAR FLEXTABLE con la nueva columna
ft_abundancia <- flextable(abundancia_tabla) %>%
        set_header_labels(
                Species = "Species",
                Overall = "Overall",
                `Fore reef_May` = "Fore reef\nMay",
                `Fore reef_August` = "Fore reef\nAugust",
                `Terrace_May` = "Terrace\nMay",
                `Terrace_August` = "Terrace\nAugust"
        ) %>%
        theme_booktabs() %>%
        font(fontname = "Times New Roman", part = "all") %>%        
        fontsize(size = 10, part = "all") %>%                       
        compose(
                j = "Species",
                part = "body",
                value = as_paragraph(
                        as_chunk(
                                gsub("^(\\S+\\s+\\S+).*", "\\1", Species),
                                props = fp_text(italic = TRUE, 
                                                font.size = 10,
                                                font.family = "Times New Roman")  # <-- agrega font.family
                        ),
                        as_chunk(
                                ifelse(
                                        str_detect(Species, "^\\S+\\s+\\S+\\s+"),
                                        paste0(" ", gsub("^\\S+\\s+\\S+\\s+(.*)", "\\1", Species)),
                                        ""
                                ),
                                props = fp_text(italic = FALSE, 
                                                font.size = 10,
                                                font.family = "Times New Roman")  # <-- agrega font.family
                        )
                )
        ) %>%
        align(j = 1, align = "left", part = "all") %>%
        align(j = 2:6, align = "center", part = "all") %>%
        padding(j = 1:6, padding.left = 3, padding.right = 3, part = "all") %>%
        padding(j = 1:6, padding.top = 2, padding.bottom = 2, part = "all") %>%
        add_header_lines(
                values = "Table 2.b. Mean density ± standard deviation of coral species in each habitat and month combination, with overall average across all samples."
        ) %>%
        align(i = 1, align = "left", part = "header") %>%
        #autofit() %>%
        width(j = 1:6, width = 2.5) %>%      # <-- ancho col 1 en pulgadas
        width(j = 2:6, width = (page_width_in-2.5)/5 ) %>%
        fit_to_width(max_width = page_width_in)
ft_abundancia

save_as_docx(ft_abundancia, path = "Table_2.b_Community_Species_Density_HabitatMonth.docx")

# Grouped by habitat for analysis
abundancia_habitat <- df_com %>%
        pivot_longer(cols = Agaaga:Steint,
                     names_to = "Species_Code",
                     values_to = "Abundance") %>%
        group_by(Species_Code, Habitat) %>%
        summarise(
                Mean = mean(Abundance, na.rm = TRUE),
                SD = sd(Abundance, na.rm = TRUE),
                .groups = "drop"
        ) %>%
        mutate(
                Label = sprintf("%.2f ± %.2f", Mean, SD),
                Species = species_names[Species_Code]
        ) %>%
        select(Species, Habitat, Label) %>%
        pivot_wider(names_from = Habitat, values_from = Label)

View(abundancia_habitat)


### Permanova --------------------------------------------------------------------

# Calculate Bray-Curtis dissimilarity 
dist_bray <- vegdist(comm_matrix, method = "bray")

# Permanova interaction marginal with Site as blocking factor 
permanova_int_marginal <- adonis2(
        dist_bray ~ Habitat * Month,
        data = metadata,
        strata = metadata$Site,
        permutations = 9999,
        #method = "bray",
        by = "margin"
)

# Permanova additive marginal with Site as blocking factor 
permanova_add_marginal <- adonis2(
        dist_bray ~ Habitat + Month,
        data = metadata,
        strata = metadata$Site,
        permutations = 9999,
        #method = "bray",
        by = "margin"
)

permanova_int_marginal
permanova_add_marginal

# Complementa PERMANOVA
anosim_habitat <- anosim(dist_bray, metadata$Habitat, permutations = 9999)
anosim_month <- anosim(dist_bray, metadata$Month, permutations = 9999)
anosim_habitat
anosim_month

# Test for homogeneity of multivariate dispersions 

# Habitat
disp_habitat <- betadisper(dist_bray, metadata$Habitat)
disp_habitat 
habitat_dispersion <- permutest(disp_habitat, pairwise = TRUE, permutations = 9999)
habitat_dispersion
anova(disp_habitat)
# Month
disp_month <- betadisper(dist_bray, metadata$Month)
month_dispersion <- permutest(disp_month, permutations = 9999)
disp_month
month_dispersion
# Interaction
metadata$Hab_Month <- interaction(metadata$Habitat, metadata$Month)
disp_interaction <- betadisper(dist_bray, metadata$Hab_Month)
habitatxmonth_dispersion <- permutest(disp_interaction, permutations = 9999)

# Visualize dispersion by habitat
par(mfrow = c(1, 2))
plot(disp_habitat, main = "PCoA - Dispersion by Habitat")
boxplot(disp_habitat, main = "Multivariate dispersion by Habitat")
dev.off()

### RESULT. TABLE 3.Permanova table --------------------------------------------------

# Extraer valores de PERMANOVA
perm_habitat_F <- permanova_add_marginal$F[1]
perm_habitat_R2 <- permanova_add_marginal$R2[1]
perm_habitat_p <- permanova_add_marginal$`Pr(>F)`[1]

perm_month_F <- permanova_add_marginal$F[2]
perm_month_R2 <- permanova_add_marginal$R2[2]
perm_month_p <- permanova_add_marginal$`Pr(>F)`[2]

perm_interaction_F <- permanova_int_marginal$F[1]
perm_interaction_R2 <- permanova_int_marginal$R2[1]
perm_interaction_p <- permanova_int_marginal$`Pr(>F)`[1]

perm_residual_R2 <- permanova_add_marginal$R2[3]
perm_total_R2 <- permanova_add_marginal$R2[4]

perm_residual_R2
perm_total_R2

betad_habitat_F <- habitat_dispersion$tab$F[1]
betad_habitat_p <- habitat_dispersion$tab$`Pr(>F)`[1]
betad_habitat_F 
betad_habitat_p 

betad_month_F <- month_dispersion$tab$F[1]
betad_month_p <- month_dispersion$tab$`Pr(>F)`[1]
betad_month_F
betad_month_p

betad_interaction_F <- habitatxmonth_dispersion$tab$F[1]
betad_interaction_p <- habitatxmonth_dispersion$tab$`Pr(>F)`[1]

# Formatear p-values
format_p <- function(p) {
        if (p < 0.001) return("< 0.001")
        else return(sprintf("= %.3f", p))
}

tabla_combined <- data.frame(
        Analysis = c("PERMANOVA", 
                     "  Habitat", 
                     "  Month", 
                     "  Residual", 
                     "  Total",
                     "", 
                     "BETADISPER", 
                     "  Habitat", 
                     "  Month"), 
        Statistic = c("", 
                      sprintf("Pseudo-F = %.2f", perm_habitat_F),
                      sprintf("Pseudo-F = %.2f", perm_month_F),
                      "—",
                      "—",
                      "",
                      "", 
                      sprintf("F = %.2f", betad_habitat_F),
                      sprintf("F = %.2f", betad_month_F)),
        R2_or_Effect = c("", 
                         sprintf("%.3f", perm_habitat_R2),
                         sprintf("%.3f", perm_month_R2),
                         sprintf("%.3f", perm_residual_R2),
                         sprintf("%.3f", perm_total_R2),
                         "",
                         "", 
                         "—", 
                         "—"),
        p_value = c("", 
                    format_p(perm_habitat_p),
                    format_p(perm_month_p),
                    "—", 
                    "—",
                    "",
                    "", 
                    format_p(betad_habitat_p),
                    format_p(betad_month_p)),
        Interpretation = c("", 
                           ifelse(perm_habitat_p < 0.05, "Significant", "Non-significant"),
                           ifelse(perm_month_p < 0.05, "Significant", "Non-significant"),
                           "—", 
                           "—",
                           "",
                           "", 
                           ifelse(betad_habitat_p < 0.05, "Heterogeneous*", "Homogeneous"),
                           ifelse(betad_month_p < 0.05, "Heterogeneous*", "Homogeneous"))
)


tabla_combined

# Calcular el ratio si aún no lo tienes
ratio_value <- perm_habitat_R2 / perm_month_R2
ratio_text <- sprintf("%.1f", ratio_value)

# Crear footer dinámico
footer_note <- sprintf(
        "Note: The interaction term was %s (p %s), thus the %s model was selected. %s effects were approximately %s-fold stronger than %s effects (R²%s/R²%s = %s).",
        ifelse(perm_interaction_p < 0.05, "significant", "non-significant"),
        format_p(perm_interaction_p),
        ifelse(perm_interaction_p < 0.05, "interactive", "additive"),
        ifelse(perm_habitat_R2 > perm_month_R2, "Habitat", "Month"),
        ratio_text,
        ifelse(perm_habitat_R2 > perm_month_R2, "Month", "Habitat"),
        ifelse(perm_habitat_R2 > perm_month_R2, "Habitat", "Month"),
        ifelse(perm_habitat_R2 > perm_month_R2, "Month", "Habitat"),
        ratio_text
)
footer_note

# Crear flextable
ft_combined <- flextable(tabla_combined) %>%
        set_header_labels(
                Analysis = "Analysis/Term",
                Statistic = "Statistic",
                R2_or_Effect = sprintf("R²"),
                p_value = "p-value",
                Interpretation = "Interpretation"
        ) %>%
        theme_booktabs() %>%
        autofit() %>%
        fontsize(size = 10, part = "all") %>%
        align(j = 2:5, align = "center", part = "all") %>%
        align(j = 1, align = "left", part = "all") %>%
        bold(i = c(1, 7), part = "body") %>%  # Solo PERMANOVA y BETADISPER
        italic(j = c(4), i = c(2, 3, 8, 9), part = "body") %>%  # p-values
        hline(i = c(5, 9), border = fp_border(color = "gray50", width = 1)) %>%
        bg(i = 6, bg = "white", part = "body") %>% 
        add_header_lines(
                values = "Table 3. Permutational multivariate analysis of variance (PERMANOVA) results evaluating the effects of Habitat (Terrace vs. Fore reef) and Month (August vs. May) on coral community composition based on Bray-Curtis dissimilarity. Site was included as a blocking factor (strata). Marginal (Type III) tests with 9,999 permutations."
        ) %>%
        align(i = 1, align = "left", part = "header") %>%
        #fontsize(i = 1, size = 10, part = "header") %>%
        add_footer_lines(values = footer_note) %>%
        align(part = "footer", align = "left") %>%
        italic(part = "footer") %>%
        fontsize(part = "footer", size = 9)


ft_combined
# Guardar
save_as_docx(ft_combined, path = "Table_3_Community_PERMANOVA_Betadisper.docx")

### NMDS for community -------------------------------------------------------------------------

# Crear etiqueta formateada de ANOSIM

perm_habitat_F 
perm_habitat_R2 
perm_habitat_p

F_val_text <- round(perm_habitat_F , 2)
F_val_text

R_val_text <- round(perm_habitat_R2, 3)
R_val_text

p_val_text <- if(perm_habitat_p < 0.001) {
        "< 0.001"
} else {
        paste0("= ", round(perm_habitat_p, 3))
}
p_val_text

set.seed(123)
nmds <- metaMDS(
        comm_matrix, 
        #distance = "bray", 
        k = 2, 
        trymax = 100,
        autotransform = F)

# Goodness of fit
nmds_stress <- nmds$stress 
stressplot(nmds)
nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
        bind_cols(metadata)

# NMDS plot

p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, 
                                  color = Habitat, fill = Habitat, shape = Month)) +
        stat_ellipse(data = nmds_scores,
                     aes(x = NMDS1, y = NMDS2,
                         group = Habitat,
                         fill = Habitat),
                     level = 0.95,
                     geom = "polygon",
                     alpha = 0.3,
                     color = NA)+
        geom_point(size = 2.5, alpha = 0.8) +
        scale_color_manual(values = c("Fore reef" = "#1E90FF", "Terrace" = "#ffba08"),
                           breaks = c("Terrace", "Fore reef"),
                           labels = c("Fore reef" = "Fore-reef", "Terrace" = "Terrace")) +
        scale_fill_manual(values = c("Fore reef" = "#1E90FF", "Terrace" = "#ffba08"),
                          breaks = c("Terrace", "Fore reef"),
                          labels = c("Fore reef" = "Fore-reef", "Terrace" = "Terrace")) +
        scale_shape_manual(values = c(16, 17),
                           breaks = c("May", "August")) +  # <-- pon aquí tus niveles reales de Month
        scale_x_continuous(limits = c(-1.1, 1.7)) +
        theme_bw(base_size = 10) +
        theme(
                axis.title = element_text(size = 11, face = "bold"),
                axis.text = element_text(size = 10),
                legend.title = element_text(size = 10, face = "bold"),
                legend.text = element_text(size = 10),
                legend.position =  "bottom",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", size = 0.5)
        ) +
        labs(x = "NMDS Axis 1", y = "NMDS Axis 2") +
        annotate("text", x = - 1.1, y = -0.6,
                 label = paste0("Stress = ", round(nmds$stress, 3),
                                "\nPseudo-F = ", round(F_val_text, 3),
                                ", R² = ", R_val_text,
                                ", p-value ", p_val_text),
                 hjust = 0, vjust = 0, size = 3.3)
p_nmds


# NMDS with species vectors (top 10) 
fit <- envfit(nmds, comm_matrix, permutations = 999)
fit
species_scores <- as.data.frame(scores(nmds, display = "species"))
species_scores$species <- rownames(species_scores)

top_species <- comm_matrix %>%
        colMeans() %>%
        sort(decreasing = TRUE) %>%
        head(10) %>%
        names()
top_species
species_top <- species_scores %>%
        filter(species %in% top_species)

p_nmds_species <- p_nmds + 
        # Vectores de especies
        geom_segment(data = species_top,
                     aes(x = 0, y = 0, 
                         xend = NMDS1, 
                         yend = NMDS2),
                     arrow = arrow(length = unit(0.15, "cm")),
                     color = "gray30", 
                     #alpha = 0.6, 
                     linewidth = 0.6,
                     inherit.aes = FALSE) +
        # Etiquetas con ggrepel (evita solapamiento)
        geom_text_repel(data = species_top,
                        aes(x = NMDS1, y = NMDS2, label = species),
                        size = 3, 
                        color = "black", 
                        fontface = "italic",
                        inherit.aes = FALSE,
                        box.padding = 0.2,        # Espacio alrededor del texto
                        point.padding = 0.3,      # Espacio desde el punto
                        segment.color = "gray30", # Color de líneas conectoras
                        segment.size = 0.3,       # Grosor de líneas conectoras
                        max.overlaps = 20)       # Permitir más etiquetas

p_nmds_species

#ggsave("Fig_NMDS.pdf", p_nmds_species, width = 6, height = 4.5)

# NMDS by Site 
p_nmds_site <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Site)) +
        geom_point(size = 2.5, alpha = 0.7) +
        stat_ellipse(level = 0.95, size = 0.5) +
        facet_wrap(~ Month + Habitat, nrow = 2) +
        theme_bw(base_size = 11) +
        labs(title = "NMDS by Site, Habitat, and Month") +
        theme(legend.position = "bottom")

p_nmds_site

### Simper: Species contribution to dissimilarity --------------------------------

# Species names dictionary
species_names <- c(
        "Sidsid" = "Siderastrea siderea",
        "Steint" = "Stephanocoenia intersepta",
        "Dicsto" = "Dichocoenia stokesi",
        "Orbfav" = "Orbicella faveolata",
        "Psestr" = "Pseudodiploria strigosa",
        "Moncav" = "Montastraea cavernosa",
        "Porast" = "Porites astreoides",
        "Orbann" = "Orbicella annularis",
        "Diplab" = "Diploria labyrinthiformis",
        "Sidrad" = "Siderastrea radians",
        "Meamea" = "Meandrina meandrites",
        "Porfur" = "Porites furcata",
        "Porpor" = "Porites porites",
        "Pordiv" = "Porites divaricata",
        "Maddec" = "Madracis decactis",
        "Eusfas" = "Eusmilia fastigiata",
        "Orbfra" = "Orbicella franksi",
        "Agalam" = "Agaricia lamarcki",
        "Favfra" = "Favia fragum",
        "Scocub" = "Scolymia cubensis",
        "Agaaga" = "Agaricia agaricites",
        "Agaund" = "Agaricia undata",
        "Madfor" = "Madracis formosa",
        "Colnat" = "Colpophyllia natans",
        "Milalc" = "Millepora alcicornis",
        "Psecli" = "Pseudodiploria clivosa",
        "Musang" = "Mussa angulosa",
        "Agafra" = "Agaricia fragilis",
        "Meajac" = "Meandrina jacksoni",
        "Mycali" = "Mycetophyllia aliciae",
        "Mycfer" = "Mycetophyllia ferox",
        "Myclam" = "Mycetophyllia lamarckiana",
        "Agahum" = "Agaricia humilis",
        "Madmir" = "Madracis mirabilis",
        "Mycdan" = "Mycetophyllia danaana",
        "Scolac" = "Scolymia lacera",
        "Lepcuc" = "Helioseris cucullata"
)

# By Month
simper_month <- simper(comm_matrix, metadata$Month, permutations = 999)
simper_month_summary <-summary(simper_month, ordered = TRUE)
simper_month
simper_month_summary 

# By Habitat 
simper_habitat <- simper(comm_matrix, metadata$Habitat, permutations = 999)
simper_habitat_summary <-summary(simper_habitat, ordered = TRUE)
print(simper_habitat)
print(simper_habitat_summary)

# Verificar el nombre exacto del contraste
names(summary(simper_habitat))

library(flextable)
library(officer)
library(dplyr)

simper_hab_df <- summary(simper_habitat)[[1]] %>%
        as.data.frame() %>%
        tibble::rownames_to_column("species") %>%
        arrange(desc(average)) %>%
        head(10) %>%
        mutate(
                species_full = species_names[species],
                habitat_indicator = case_when(
                        p >= 0.05 ~ "Non-significant",
                        ava > avb ~ "Fore reef",
                        avb > ava ~ "Terrace"
                ),
                fold_change = ifelse(
                        ava > avb, 
                        ava / avb, 
                        avb / ava
                        ),
                fold_change_label = ifelse(
                        p < 0.05, paste0(
                                round(fold_change, 1), "×"), 
                        ""
                        ),
                p_label = case_when(
                        p < 0.001 ~ "p < 0.001",
                        p < 0.01 ~ "p < 0.01",
                        p < 0.05 ~ "p < 0.05",
                        TRUE ~ "ns"
                ),
                label_out = case_when(
                        p < 0.05 & !is.na(fold_change) ~ 
                                paste0("×", round(fold_change, 1),",  ", p_label),
                        p < 0.05 ~ p_label,
                        TRUE ~ "")
        )
simper_hab_df

# SIMPER plot
p_simper <- ggplot(simper_hab_df, 
                   aes(x = reorder(species_full, average), y = average,
                       fill = habitat_indicator)) +
        geom_col(alpha = 1, width = 0.6) +
        scale_fill_manual(values = c("Fore reef" = "#1E90FF", "Terrace" = "#ffba08", "Non-significant" = "gray70"),
                          breaks = c("Terrace", "Fore reef", "Non-significant"),
                          labels = c("Terrace" = "Terrace", "Fore reef" = "Fore-reef",  "Non-significant"="Non-significant" )) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
        geom_text(aes(label = label_out),
                  hjust = -0.1,
                  size = 3.3) +
        coord_flip(clip = "off") +
        labs(x = "", 
             y = "Average contribution to dissimilarity",
             fill = "Habitat") +
        theme_bw(base_size = 10) +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_blank(),   
              panel.grid.minor.x = element_blank(),   
              panel.grid.minor.y = element_blank(),   
              axis.title = element_text(size = 11, face = "bold"),    # títulos de eje
              axis.text.y = element_text(face = "italic", size = 10, color = "gray30"),
              axis.text.x = element_text(size = 10.5),
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 10),
              legend.position =  "bottom",
              panel.border = element_rect(color = "black", size = 0.25))
              
p_simper
#ggsave("Fig_SIMPER_Habitat.pdf", p_simper, width = 6, height = 4.5)

### RESULT. FIGURE 3. Community Structure Panel ----------------------------------------------------
library(patchwork)
combined <- p_nmds_species + p_simper +
        plot_layout(ncol = 2, widths = c(1.5, 1.3)) +
        plot_annotation(tag_levels = "A") &
        theme(
                plot.tag = element_text(size = 12, face = "bold"),
                plot.tag.position = c(0.02, 0.98),
                #plot.margin = margin(t = 0, r = 0, b = 0, l = 5, unit = "mm")
                #plot.margin = margin(l = 5, unit = "mm")
        )
combined

ggsave("Figure_3_Community_NMDs_SIMPER.pdf", combined,
       width = 12, height = 5.5, dpi = 400)

### RESULT. TABLE 4. Simper table for community and habitat differences ------------------------------------------------------
# Preparar datos

# Denominador correcto: suma de TODAS las especies del SIMPER, antes de recortar a 10
total_average <- sum(summary(simper_habitat)[[1]]$average)
total_average
simper_table <- simper_hab_df %>%
        mutate(
                Species = species_full,
                `Contribution (%)` = sprintf("%.1f", average * 100 / total_average),  # <-- corregido
                `Fore reef (%)` = sprintf("%.2f", ava),
                `Terrace (%)` = sprintf("%.2f", avb),
                `Fold change` = sprintf("%.1f×", fold_change),
                `p-value` = case_when(
                        p < 0.001 ~ "<0.001***",
                        p < 0.01 ~ sprintf("%.3f**", p),
                        p < 0.05 ~ sprintf("%.3f*", p),
                        TRUE ~ sprintf("%.3f", p)
                ),
                Indicator = habitat_indicator
        ) %>%
        select(Species, `Contribution (%)`, `Fore reef (%)`, 
               `Terrace (%)`, `Fold change`, `p-value`, Indicator)

# Crear flextable
ft_simper <- flextable(simper_table) %>%
        theme_booktabs() %>%
        # Anchos
        width(j = 1, width = 1.8) %>%   # Species
        width(j = 2, width = 1.0) %>%   # Contribution (%)
        width(j = 3, width = 0.7) %>%   # Fore reef (%)
        width(j = 4, width = 0.7) %>%   # Terrace (%)
        width(j = 5, width = 0.65) %>%  # Fold change
        width(j = 6, width = 0.7) %>%   # p-value
        width(j = 7, width = 1.1) %>%   # Indicator
        # Alineación
        align(j = 1, align = "left", part = "all") %>%
        align(j = 2:7, align = "center", part = "all") %>%
        # Formato
        italic(j = 1, part = "body") %>%
        #bold(part = "header") %>% 
        # Fuentes
        fontsize(size = 10, part = "body") %>%
        fontsize(size = 10, part = "header") %>%
        # Padding
        padding(j = 1:7, padding.left = 3, padding.right = 3, part = "all") %>%
        padding(j = 1:7, padding.top = 2, padding.bottom = 2, part = "all") %>%
        # HEADING (dividido en 2 líneas para legibilidad)
        add_header_lines(
                "Table 4. Similarity percentage (SIMPER) analysis identifying coral species contributing most to compositional dissimilarity between fore-reef and terrace habitats. Species are ranked by their contribution to overall Bray-Curtis dissimilarity; the top ten species are shown, of which the top five account for 74.4 % of the dissimilarity and the top six for 78.9 %. Fold change indicates enrichment ratio. Indicator classification based on permutation tests (n = 999). Average dissimilarity = 29.3 %. * p < 0.05, ** p < 0.01, *** p < 0.001."
        ) 

ft_simper
save_as_docx(ft_simper, path = "Table_4_Community_SIMPER_Habitat.docx")

### Rank-Abundance Curves --------------------------------------------------------
library(BiodiversityR)
# Por grupo
rankabundance_fore <- rankabundance(
        comm_matrix[metadata$Habitat == "Fore reef", ]
)
rankabundance_terrace <- rankabundance(
        comm_matrix[metadata$Habitat == "Terrace", ]
)
# Gráfico
par(mfrow = c(1, 2))
rankabunplot(rankabundance_fore, main = "Fore reef")
rankabunplot(rankabundance_terrace, main = "Terrace")

permanova_int_marginal
permanova_add_marginal
habitat_dispersion
month_dispersion
habitatxmonth_dispersion
anosim_habitat
anosim_month
species_top
fit
nmds_stress
simper_habitat_summary
simper_habitat
simper_month_summary
simper_month
### Clean Space -----
rm(list = ls())
graphics.off()
### Summary Results: Community Structure ------------------------------------------------

# CRITICAL METHODOLOGICAL CAVEAT:
# BETADISPER revealed significant heterogeneity in multivariate dispersions between habitats
# (F = 4.11, p = 0.045). Terrace communities are 15.9% more dispersed (0.197 vs. 0.170).
# This violates PERMANOVA's homogeneity of variance assumption and may partially inflate
# the habitat F-statistic. ANOSIM validation (below) corroborates true compositional differences
# exist beyond variance heterogeneity alone.

# PERMANOVA - MODEL SELECTION:
# Interaction model: Habitat × Month not significant (R² = 0.007, F = 1.50, p = 0.093)
# → Additive model selected (Habitat + Month)
# 
# PERMANOVA - ADDITIVE MODEL (Marginal tests, by = "margin"):
# Overall: R² = 0.124, F = 13.97, p < 0.001 (12.4% variance explained)
# - Habitat: R² = 0.107, F = 23.97, p < 0.001 *** (10.7% variance, PRIMARY DRIVER)
#   [NOTE: F-statistic may be partially inflated by dispersion heterogeneity, but effect is genuine]
# - Month:   R² = 0.018, F = 3.97,  p < 0.001 *** (1.8% variance, secondary)
# - Site controlled as blocking factor (strata)
# Habitat effects 5.9× stronger than temporal effects (R² ratio: 0.107/0.018)
#
# ANOSIM (Complementary rank-based validation):
# ANOSIM is more robust to dispersion heterogeneity than PERMANOVA
# - Habitat: R = 0.136, p < 0.001 *** (moderate, statistically significant separation)
# - Month:   R = 0.032, p < 0.001 *** (weak but statistically significant separation)
# Both habitat and month indices 4× larger in PERMANOVA than ANOSIM (as expected given dispersion differences)
# PERMANOVA and ANOSIM AGREE on relative effect sizes: habitat >> month
# Concordance across both methods strengthens evidence for true compositional differences
# 
# BETADISPER (Homogeneity of multivariate dispersions):
# - Habitat: F = 4.11, p = 0.045 * (SIGNIFICANT HETEROGENEITY)
#   → Terrace communities 15.9% more dispersed: 0.197 vs. 0.170
#   → Interpretation: Greater compositional variance within terrace habitat
#      (higher beta-diversity, more heterogeneous microhabitat assembly, or stochastic variability)
# - Month: F = 0.05, p = 0.829 ns (homogeneous dispersion)
#   → Within-month community variance consistent across seasons
# - Habitat×Month:  F = 1.55, p = 0.203 ns (homogeneous dispersion)
#   → No evidence of differential dispersion patterns in habitat × season combinations

# INTERPRETATION SUMMARY:
# Habitat differences are driven by TWO components:
# (1) Compositional centroid shift (different mean species abundances) — CONFIRMED by ANOSIM
# (2) Dispersion increase (higher beta-diversity in terraces) — CONFIRMED by BETADISPER
# Both contribute to multivariate separation; PERMANOVA conflates both effects.
# ANOSIM evidence that centroid shift is statistically robust.

# NMDS ORDINATION:
# - Stress = 0.135 (good fit: below 0.15 threshold)
# - Clear habitat separation along NMDS1 axis
# - Temporal patterns weaker but detectable along secondary axes
# - Dominant species: Sidsid, Steint, Orbfav, Dicsto
#
# SIMPER - HABITAT (Fore reef vs Terrace)
# Average Bray-Curtis dissimilarity = 0.293 (29.3%)
# Top 8 species explain 79.9% of habitat dissimilarity
# 
# GENERALIST (no significant habitat preference):
# - Sidsid (S. siderea): 20.6% contribution to dissimilarity, p = 0.102 ns
#   Abundance: 52.3% (fore reef) vs. 50.4% (terrace), 1.04× ratio
#   → Dominant in both habitats equally; true environmental generalist
# 
# FORE REEF INDICATORS (significantly more abundant, p < 0.05):
# - Orbfav (O. faveolata):    12.6% contribution, p < 0.001 ***, 2.2× enriched
#   Abundance: 10.3% vs. 4.6%
# - Moncav (M. cavernosa):    5.6% contribution, p = 0.003 **, 1.4× enriched
#   Abundance: 4.6% vs. 3.3%
# - Porast (P. astreoides):   5.1% contribution, p < 0.001 ***, 1.6× enriched
#   Abundance: 3.8% vs. 2.4%
# - Porfur (P. furcata):      1.6% contribution, p = 0.005 **, 3.6× enriched
#   Abundance: 0.82% vs. 0.23%
# → ECOLOGICAL TRAIT: Massive, framework-building corals
# → INTERPRETATION: Favor deeper, stable fore reef conditions with reduced wave disturbance
#    and lower chronic physical stress. Space pre-emption and competitive dominance
#    drive assembly under low-disturbance equilibrium conditions.
# 
# TERRACE INDICATORS (significantly more abundant, p < 0.05):
# - Dicsto (D. stokesi):      14.1% contribution, p < 0.001 ***, 5.9× enriched
#   Abundance: 8.8% vs. 1.5%
# - Steint (S. intersepta):   14.6% contribution, p < 0.001 ***, 1.2× enriched
#   Abundance: 20.3% vs. 17.2%
# - Psestr (P. strigosa):     7.3% contribution, p < 0.001 ***, 3.7× enriched
#   Abundance: 5.1% vs. 1.4%
# - Pordiv (P. divaricata):   1.1% contribution, p < 0.001 ***, 4.3× enriched
#   Abundance: 0.56% vs. 0.13%
# → ECOLOGICAL TRAIT: Small-polyped, morphologically flexible species with rapid recovery
# → INTERPRETATION: Dominate shallow, high-energy terrace environments with chronic
#    wave disturbance, variable light, and unstable substrate. Physiological stress-tolerance
#    and rapid growth/reproduction capacity override competitive ability as assembly driver.

# SIMPER - MONTH (May vs August):
# Average Bray-Curtis dissimilarity = 0.281 (28.1%)
# WEAK TEMPORAL VARIATION: Dissimilarity between months (~28%) comparable to habitat
# dissimilarity (~29%), but PERMANOVA variance explained is 5.9× smaller (1.8% vs. 10.7%)
# This indicates temporal changes are modest in magnitude despite statistical significance.
# 
# SIGNIFICANT TEMPORAL CHANGES (p < 0.05):
# Species INCREASING May → August:
# - Orbann (O. annularis):    3.4% contribution, p < 0.001 ***, 2.1× increase
#   Abundance: 0.75% (May) → 1.61% (August)
# - Moncav (M. cavernosa):    5.8% contribution, p = 0.010 **, 1.3× increase
#   Abundance: 3.48% (May) → 4.44% (August)
# - Orbfav (O. faveolata):    11.7% contribution, p = 0.034 *, 1.2× increase
#   Abundance: 6.75% (May) → 8.18% (August)
# - Diplab (D. labyrinthiformis): 2.6% contribution, p = 0.010 **, 1.7× increase
#   Abundance: 0.81% (May) → 1.38% (August)
# → All are competitively dominant framework-builders; may reflect summer recruitment pulses
#    or phenological activity changes from May to August
# 
# Species DECREASING May → August:
# - Sidrad (S. radians):      2.7% contribution, p < 0.001 ***, 4.1× decrease
#   Abundance: 1.30% (May) → 0.32% (August)
# - Porfur (P. furcata):      1.7% contribution, p = 0.008 **, 1.8× decrease
#   Abundance: 0.68% (May) → 0.38% (August)
# → Both are smaller-polyped species; may reflect post-recruitment mortality or
#    summer stress-induced colony loss
# 
# NON-SIGNIFICANT TEMPORAL TRENDS (p > 0.05):
# - Sidsid (S. siderea):      21.4% contribution, p = 0.088 ns
#   Abundance: 52.1% (May) vs. 50.6% (August) — dominant generalist stable over time
# - Steint (S. intersepta):   14.9% contribution, p = 0.051 ns (marginally non-significant)
#   Abundance: 19.8% (May) vs. 17.7% (August) — terrace dominant stable across seasons
# 
# NOTE: Temporal effects are WEAK compared to habitat effects (Month R² = 1.8% vs. Habitat R² = 10.7%)
# Community structure fundamentally stable; seasonal changes affect only subset of species.
# 
# KEY ECOLOGICAL INTERPRETATION
# 1. HABITAT IS THE DOMINANT STRUCTURING FORCE (5.9× stronger than season)
#    - Explains 10.7% variance vs. 1.8% for Month
#    - No Habitat × Month interaction (p = 0.093) indicates habitat signatures stable across seasons
#    - PERMANOVA habitat effect (F = 23.97, p < 0.001) confirmed by ANOSIM (R = 0.136, p < 0.001)
#    - Both statistical approaches agree: habitat >> temporal patterns
# 
# 2. TWO DISTINCT REEF ASSEMBLAGES WITH CONTRASTING FUNCTIONAL STRATEGIES
#    
#    FORE REEF ASSEMBLAGES: Competitive dominants in stable conditions
#    - Dominated by massive framework-builders (Orbicella, Montastraea, Porites spp.)
#    - Slow-growing, space pre-emptive, superior competitors
#    - Reflect deeper, low-disturbance reef conditions
#    - Assembly driven by competitive exclusion and equilibrium conditions
#    
#    TERRACE ASSEMBLAGES: Stress-tolerant opportunists in dynamic conditions
#    - Dominated by small-polyped, morphologically flexible species (Dichocoenia, Stephanocoenia)
#    - Rapid growth/recovery capacity, physiological stress tolerance
#    - Shallow, high-energy, variable light conditions
#    - Assembly driven by disturbance, environmental heterogeneity, and stochasticity
# 
# 3. TERRACE COMMUNITIES ARE MORE COMPOSITIONALLY HETEROGENEOUS
#    - 15.9% higher within-habitat beta-diversity (BETADISPER: p = 0.045)
#    - Suggests mosaic of microhabitats with varying exposure/substrate/hydrodynamics
#    - May reflect stochastic assembly processes more sensitive to local variation
#    - Potential implications for ecological stability and resilience
# 
# 4. SIDERASTREA SIDEREA: ECOLOGICAL GENERALIST SPECIES
#    - 50% relative abundance dominates both habitats equally (1.04× ratio, p = 0.102 ns)
#    - High habitat generality despite contributing 20.6% to dissimilarity
#    - Reflects broad environmental tolerance and consistent ecosystem role
# 
# 5. TEMPORAL PATTERNS ARE SUBTLE, AFFECTING ONLY RARE SPECIES
#    - Weak seasonal shifts in lower-abundance species (Orbicella annularis, Siderastrea radians)
#    - May reflect post-disturbance mortality
#    - Dominant species (Sidsid, Steint) show compositional stability across seasons
#    - Overall assemblage structure fundamentally stable (no interaction effect)
# 
# 6. SPATIAL GRADIENTS OVERRIDE TEMPORAL FLUCTUATIONS
#    - Environmental filtering along reef zone gradients >> seasonal variability
#    - Habitat-driven zonation patterns consistent across both sampling months
#    - Strong implication: Reef zone preservation critical for maintaining community diversity
#    - Different reef zones = distinct functional assemblages with contrasting ecological strategies
# 
# 7. METHODOLOGICAL ROBUSTNESS
#    - BETADISPER identified heterogeneous dispersions between habitats (p = 0.045)
#    - PERMANOVA habitat F-statistic may be partially inflated by variance heterogeneity
#    - ANOSIM validation (rank-based, more robust) confirms significant compositional separation
#    - Concordance between PERMANOVA and ANOSIM provides strong evidence for true differences
#    - Conservative interpretation: Compositional differences are genuine, but dispersion
#      heterogeneity contributes additional variance in terrace communities

# 2. CORAL HEALTH ANALYSIS------------------------------
### 2.1 Health Composition Analysis --------------------------------------------
# Load packages 
library(readxl)
library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)
library(flextable)
library(officer)
library(ggrepel)
library(patchwork)

### Load and prepare data -------------------------------------------------------
df_health_prev <- read_excel(
        "D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/Corals_data_JR_Forereefs_Terrace_2024_August_May.xlsx",
        sheet = "Health_prev"
) %>%
        filter(Site != "Auras") %>%
        droplevels() %>%
        mutate(
                Month = factor(Month, levels = c("May", "August")),
                Habitat = factor(Habitat, levels = c("Fore reef", "Terrace")),
                Site = factor(Site, levels = c("Peruano", "Pinos", "Anclitas", 
                                               "Mariflores", "Cruces")),
                Transept = factor(Transept))
View(df_health_prev)
### Data exploration -------------------------------------------------------------
glimpse(df_health_prev)
table(df_health_prev$Site, df_health_prev$Habitat)

### Extract health matrix and metadata ----------------------------------------

prev_matrix <- df_health_prev %>%
        select(Bl, AM, RM, RTLD, BBD, DSS, BIOERO_SPO, BIOERO_POL, Diseased) %>%
        as.matrix()

colnames(prev_matrix) <- c("Bleaching",
                           "Old Mortality",
                           "Recent Mortality",
                           "Rapid Tissue Loss Disease",
                           "Black Band Disease",
                           "Dark Spot Disease",
                           "Sponge Bioerosion",
                           "Polychaete Bioerosion",
                           "Disease")

metadata <- df_health_prev %>%
        select(Month, Habitat, Site, Transept)

# Data summary
dim(prev_matrix)  # 200 samples × 8 variables
glimpse(prev_matrix)
summary(prev_matrix)

# Opción 1: Crear un dataframe con media y SD para cada variable
summary_stats <- as.data.frame(prev_matrix) %>%
        summarise(
                across(everything(), 
                       list(Mean = mean, SD = sd),
                       .names = "{.col}_{.fn}")
        ) %>%
        pivot_longer(everything(),
                     names_to = c("Variable", "Statistic"),
                     names_sep = "_",
                     values_to = "Value") %>%
        pivot_wider(names_from = Statistic, values_from = Value) %>%
        mutate(Label = sprintf("%.2f ± %.2f", Mean, SD)) %>%
        select(Variable, Mean, SD, Label)

summary_stats


### RESULT. TABLE 5. Prevalences per Habitat and Month--------

# 1. Crear un dataframe con media ± SD para cada condición de salud en cada combinación
prevalencia_por_grupo <- as.data.frame(prev_matrix) %>%
        bind_cols(metadata) %>%
        pivot_longer(cols = Bleaching:Disease,  # Todas las condiciones
                     names_to = "Health_Condition",
                     values_to = "Prevalence") %>%
        group_by(Health_Condition, Habitat, Month) %>%
        summarise(
                Mean = mean(Prevalence, na.rm = TRUE),
                SD = sd(Prevalence, na.rm = TRUE),
                N = n(),
                .groups = "drop"
        ) %>%
        mutate(
                Label = sprintf("%.2f ± %.2f", Mean, SD)
        ) %>%
        select(Health_Condition, Habitat, Month, Label)
prevalencia_por_grupo

# 2. Pivotar a formato ancho (Health_Condition × Habitat-Month combinations)
prevalencia_tabla <- prevalencia_por_grupo %>%
        mutate(
                Health_Condition = factor(Health_Condition, 
                                          levels = c("Bleaching",
                                                     "Recent Mortality",
                                                     "Old Mortality",
                                                     "Disease",
                                                     "Rapid Tissue Loss Disease",
                                                     "Black Band Disease",
                                                     "Dark Spot Disease",
                                                     "Sponge Bioerosion",
                                                     "Polychaete Bioerosion"))
        ) %>%
        pivot_wider(
                names_from = c(Habitat, Month),
                values_from = Label,
                names_sep = "_"
        ) %>%
        arrange(Health_Condition) %>%
        rename(Condition = Health_Condition)

# Ver la tabla
prevalencia_tabla

abundancia_tabla <- abundancia_global %>%
        left_join(
                abundancia_por_grupo %>%
                        pivot_wider(
                                names_from = c(Habitat, Month),
                                values_from = Label,
                                names_sep = "_"
                        ),
                by = "Species"
        ) %>%
        select(Species, Overall, everything())

abundancia_tabla

# 3. CREAR FLEXTABLE con estilo similar a tu Table 2
ft_prevalencia <- flextable(prevalencia_tabla) %>%
        set_header_labels(
                Condition = "Health Condition",
                `Fore reef_May` = "Fore reef\nMay (%)",
                `Fore reef_August` = "Fore reef\nAugust (%)",
                `Terrace_May` = "Terrace\nMay (%)",
                `Terrace_August` = "Terrace\nAugust (%)"
        ) %>%
        theme_booktabs() %>%
        autofit() %>%
        fontsize(size = 9, part = "all") %>%
        align(j = 1, align = "left", part = "all") %>%  # Condition name a la izquierda
        align(j = 2:5, align = "center", part = "all") %>%  # Prevalencias centradas
        padding(j = 1:5, padding.left = 3, padding.right = 3, part = "all") %>%
        padding(j = 1:5, padding.top = 2, padding.bottom = 2, part = "all") %>%
        add_header_lines(
                values = "Table 5. Mean prevalence (%) ± standard deviation of coral health conditions in each habitat and month combination. Values represent the average and variability across all transects within each group. RTLD includes White Plague Disease (WPD), White Band Disease (WBD), and Stony Coral Tissue Loss Disease (SCTLD)."
        ) %>%
        align(i = 1, align = "left", part = "header")
ft_prevalencia

# 4. GUARDAR
save_as_docx(ft_prevalencia, path = "Table_5_Health_Prevalence_HabitatMonth.docx")

prev_matrix <- prev_matrix %>%
        as.data.frame() %>%
        select(-Disease) %>%
        as.matrix()
prev_matrix
### Permanova --------------------------------------------------------------------

#prev_matrix_hell <- decostand(prev_matrix, method = "hellinger")

# Calculate Bray-Curtis dissimilarity 
dist_bray <- vegdist(prev_matrix, method = "bray")
# Permanova with Site as blocking factor 
permanova_prev_int_margin <- adonis2(
        dist_bray ~ Habitat * Month,
        data = metadata,
        strata = metadata$Site,
        permutations = 9999,
        #method = "bray",
        by = "margin"
)

# Permanova with Site as blocking factor 
permanova_prev_add_margin <- adonis2(
        dist_bray ~ Habitat + Month,
        data = metadata,
        strata = metadata$Site,
        permutations = 9999,
        #method = "bray",
        by = "margin"
)

permanova_prev_int_margin
permanova_prev_add_margin

# ANOSIM Complementa PERMANOVA
anosim_habitat <- anosim(dist_bray, metadata$Habitat, permutations = 9999)
anosim_month <- anosim(dist_bray, metadata$Month, permutations = 9999)
metadata$Hab_Month <- interaction(metadata$Habitat, metadata$Month)
anosim_hab_month <- anosim(dist_bray, metadata$Hab_Month, permutations = 9999)
anosim_habitat
anosim_month
anosim_hab_month
# Test for homogeneity of multivariate dispersions 

# Habitat
disp_habitat <- betadisper(dist_bray, metadata$Habitat)
habitat_dispersion <- permutest(disp_habitat, pairwise = TRUE, permutations = 9999)
disp_habitat 
habitat_dispersion
#Fore reef   Terrace 
#0.1964    0.2359 
#F = 4.29, p = 0.0357

# Month
disp_month <- betadisper(dist_bray, metadata$Month)
month_dispersion <- permutest(disp_month, permutations = 9999)
disp_month
month_dispersion
#May August 
#0.2393 0.1880 
#F = 6.9358, p = 0.008

# Interaction
metadata$Hab_Month <- interaction(metadata$Habitat, metadata$Month)
disp_interaction <- betadisper(dist_bray, metadata$Hab_Month)
habitatxmonth_dispersion <- permutest(disp_interaction, permutations = 9999)
disp_interaction
habitatxmonth_dispersion
#Fore reef.May      Terrace.May Fore reef.August   Terrace.August 
#0.1543           0.2966           0.1978           0.1609 
#F = 16.168, p < 0.001


# Visualize dispersion 
par(mfrow = c(1, 2))
plot(disp_habitat, main = "PCoA - Dispersion by Habitat")
boxplot(disp_habitat, main = "Multivariate dispersion by Habitat")
dev.off()
par(mfrow = c(1, 2))
plot(disp_month, main = "PCoA - Dispersion by Month")
boxplot(disp_month, main = "Multivariate dispersion by Month")
dev.off()
par(mfrow = c(1, 2))
plot(disp_interaction, main = "PCoA - Dispersion by Habitat x Month")
boxplot(disp_interaction, main = "Multivariate dispersion by Habitat x Month")
dev.off()

### RESULT. TABLE 5. Permanova table --------------------------------------------------

# Extraer valores de PERMANOVA
permanova_prev_add_margin

perm_habitat_F <- permanova_prev_add_margin$F[1]
perm_habitat_R2 <- permanova_prev_add_margin$R2[1]
perm_habitat_p <- permanova_prev_add_margin$`Pr(>F)`[1]
perm_habitat_F 
perm_habitat_R2
perm_habitat_p

perm_month_F <- permanova_prev_add_margin$F[2]
perm_month_R2 <- permanova_prev_add_margin$R2[2]
perm_month_p <- permanova_prev_add_margin$`Pr(>F)`[2]

perm_residual_R2_add <- permanova_prev_add_margin$R2[3]
perm_total_R2_add <- permanova_prev_add_margin$R2[4]
perm_residual_R2_add
perm_total_R2_add


permanova_prev_int_margin
perm_interaction_F <- permanova_prev_int_margin$F[1]
perm_interaction_R2 <- permanova_prev_int_margin$R2[1]
perm_interaction_p <- permanova_prev_int_margin$`Pr(>F)`[1]

perm_residual_R2 <- permanova_prev_int_margin$R2[2]
perm_total_R2 <- permanova_prev_int_margin$R2[3]
perm_residual_R2
perm_total_R2

betad_habitat_F <- habitat_dispersion$tab$F[1]
betad_habitat_p <- habitat_dispersion$tab$`Pr(>F)`[1]

betad_month_F <- month_dispersion$tab$F[1]
betad_month_p <- month_dispersion$tab$`Pr(>F)`[1]

betad_interaction_F <- habitatxmonth_dispersion$tab$F[1]
betad_interaction_p <- habitatxmonth_dispersion$tab$`Pr(>F)`[1]

# Verificar la estructura primero
str(permanova_prev_int_margin)

# Formatear p-values
format_p <- function(p) {
        if (p < 0.001) return("< 0.001")
        else return(sprintf("= %.3f", p))
}

tabla_combined <- data.frame(
        Model = c("", 
                     "  Additive", 
                     "  Additive", 
                     "  Additive", 
                     "  Additive",
                     "  Interaction",
                     "  Interaction", 
                     "  Interaction",
                     "", 
                     "", 
                     "", 
                     ""),
        Analysis = c("PERMANOVA", 
                     "  Habitat", 
                     "  Month", 
                     "  Residual", 
                     "  Total",
                     "  Habitat × Month",
                     "  Residual", 
                     "  Total",
                     "BETADISPER", 
                     "  Habitat", 
                     "  Month", 
                     "  Habitat × Month"),
        Statistic = c("", 
                      sprintf("Pseudo-F = %.2f", perm_habitat_F),
                      sprintf("Pseudo-F = %.2f", perm_month_F),
                      "—",
                      "—",
                      sprintf("Pseudo-F = %.2f", perm_interaction_F),
                      "—",
                      "—",
                      "", 
                      sprintf("F = %.2f", betad_habitat_F),
                      sprintf("F = %.2f", betad_month_F),
                      sprintf("F = %.2f", betad_interaction_F)),
        R2_or_Effect = c("", 
                         sprintf("%.3f", perm_habitat_R2),
                         sprintf("%.3f", perm_month_R2),
                         sprintf("%.3f", perm_residual_R2_add),
                         sprintf("%.3f", perm_total_R2_add),
                         sprintf("%.3f", perm_interaction_R2),
                         sprintf("%.3f", perm_residual_R2),
                         sprintf("%.3f", perm_total_R2),
                         "", 
                         "—", 
                         "—", 
                         "—"),
        p_value = c("", 
                    format_p(perm_habitat_p),
                    format_p(perm_month_p),
                    "—", 
                    "—",
                    format_p(perm_interaction_p),
                    "—", 
                    "—",
                    "", 
                    format_p(betad_habitat_p),
                    format_p(betad_month_p),
                    format_p(betad_interaction_p)),
        Interpretation = c("", 
                           ifelse(perm_habitat_p < 0.05, "Significant", "Non-significant"),
                           ifelse(perm_month_p < 0.05, "Significant", "Non-significant"),
                           "—", 
                           "—",
                           ifelse(perm_interaction_p < 0.05, "Significant", "Non-significant"),
                           "—", 
                           "—",
                           "", 
                           ifelse(betad_habitat_p < 0.05, "Heterogeneous", "Homogeneous"),
                           ifelse(betad_month_p < 0.05, "Heterogeneous", "Homogeneous"),
                           ifelse(betad_interaction_p < 0.05, "Heterogeneous", "Homogeneous"))
)

tabla_combined

# Calcular el ratio si aún no lo tienes
ratio_value <- perm_month_R2 / perm_habitat_R2
ratio_text <- sprintf("%.1f", ratio_value)
ratio_text

# Crear flextable
ft_combined <- flextable(tabla_combined) %>%
        set_header_labels(
                Model = "Model",
                Analysis = "Analysis/Term",
                Statistic = "Statistic",
                R2_or_Effect = sprintf("R²"),
                p_value = "p-value",
                Interpretation = "Interpretation"
        ) %>%
        theme_booktabs() %>%
        autofit() %>%
        align(j = 2:5, align = "center", part = "all") %>%
        align(j = 1, align = "left", part = "all") %>%
        bold(i = c(1, 9), part = "body") %>%
        fontsize(part = "all", size = 10) %>%
        #italic(i = c(2, 3, 4, 5, 6, 9, 10, 11), part = "body") %>%
        italic(j = c(4), i = c(2,3,4,9,10,11), part = "body") %>%
        hline(i = c(8), border = fp_border(color = "gray40", width = 1)) %>%
        bg(i = 7, bg = "white", part = "body") %>% 
        merge_at(i = 2:5, j = 1, part = "body") %>%  # Filas 1-3, columna 1
        merge_at(i = 6:8, j = 1, part = "body") %>%  # Filas 1-3, columna 1
        hline(i = c(1,5,9), border = fp_border(color = "gray70", width = 0.5)) %>%
        add_header_lines(
                values = "Table 6. Permutational multivariate analysis of variance (PERMANOVA) results evaluating the effects of Habitat (Terrace vs. Fore reef), Month (August vs. May) and their interaction on coral health condition prevalence based on Bray-Curtis dissimilarity prevalence data. Site was included as a blocking factor (strata). Marginal tests with 9,999 permutations."
        ) %>%
        align(i = 1, align = "left", part = "header") %>%
        fontsize(i = 1, size = 10, part = "header")

ft_combined

# Guardar
save_as_docx(ft_combined, path = "Table_6_Health_PERMANOVA_Betadisper.docx")

### NMDS -------------------------------------------------------------------------

# Crear etiqueta formateada de ANOSIM

perm_interaction_F <- permanova_prev_int_margin$F[1]
perm_interaction_R2 <- permanova_prev_int_margin$R2[1]
perm_interaction_p <- permanova_prev_int_margin$`Pr(>F)`[1]

F_val_text <- round(perm_interaction_F , 2)
F_val_text

R_val_text <- round(perm_interaction_R2, 3)
R_val_text

p_val_text <- if(perm_interaction_p < 0.001) {
        "< 0.001"
} else {
        paste0("= ", round(perm_interaction_p, 3))
}
p_val_text


set.seed(123)
nmds <- metaMDS(
        dist_bray, 
        #distance = "bray", 
        k = 2, 
        trymax = 100,
        autotransform = F)

# Goodness of fit
nmds_stress <- nmds$stress  #  (good fit)
stressplot(nmds)
nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
        bind_cols(metadata)

nmds_scores 
unique(nmds_scores$Hab_Month)
nmds_scores <- nmds_scores %>%
        mutate(Hab_Month = gsub("\\.", " ", Hab_Month))

p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
        #xlim(-1.3, 2.1) +
        #ylim(-1.2, 1.3) +
        stat_ellipse(aes(group = Hab_Month, fill = Hab_Month,
                     color = Hab_Month),
                     level = 0.95,
                     geom = "polygon",
                     alpha = 0.05,
                     linewidth = 0.7
                     ) +
        geom_point(aes(color = Hab_Month, shape = Habitat),
                   size = 2, alpha = 1
                   ) +
        scale_color_manual(
                name = "Habitat × Month",
                values = c(
                        "Terrace May" = "#F57F17",
                        "Terrace August" = "#ffba08",
                        "Fore reef May" = "#4169E1",
                        "Fore reef August" = "#4ecdc4"
                )
        ) +
        scale_fill_manual(
                name = "Habitat × Month",
                values = c(
                        "Terrace May" = "#F57F17",
                        "Terrace August" = "#ffba08",
                        "Fore reef May" = "#4169E1",
                        "Fore reef August" = "#4ecdc4"
                )
        ) +
        scale_shape_manual(
                name = "Habitat",
                values = c("Fore reef" = 16, "Terrace" = 17)
        ) +
        guides(color = guide_legend(order = 1, nrow = 2),
               fill = guide_legend(order = 1, nrow = 2),
               shape = guide_legend(order = 2, nrow = 2)) +
        theme_bw(base_size = 10) +
        theme(
                axis.title = element_text(size = 10, face = "bold"),
                axis.text = element_text(size = 10),
                legend.title = element_text(size = 10, face = "bold"),
                legend.text = element_text(size = 9),
                legend.position = "bottom",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", linewidth = 0.5)
        ) +
        annotate("text", x = -1.1, y = -1,
                 label = paste0("Stress = ", round(nmds$stress, 3),
                                "\nPseudo-F = ", round(F_val_text, 3),
                                ", R² = ", R_val_text,
                                ", p-value ", p_val_text),
                 hjust = 0, vjust = 0, size = 3.3)

print(p_nmds)

# NMDS with species vectors (top 10) 
fit <- envfit(nmds, prev_matrix, permutations = 999)
fit
species_scores <- as.data.frame(scores(fit, display = "vectors"))
species_scores$Variable <- rownames(species_scores)
species_scores
p_nmds_species <- p_nmds +
        geom_segment(data = species_scores,
                     aes(x = 0, y = 0,
                         xend = NMDS1,
                         yend = NMDS2),
                     arrow = arrow(length = unit(0.2, "cm")),
                     color = "gray30",
                     linewidth = 0.6,
                     inherit.aes = FALSE) +
        geom_text_repel(data = species_scores,
                        aes(x = NMDS1, y = NMDS2, label = Variable),
                        size = 3,
                        color = "black",
                        segment.color = "darkgray",
                        hjust = 0,
                        inherit.aes = FALSE)

p_nmds_species

#ggsave("Fig_NMDS_Health.pdf", p_nmds_species, width = 6, height = 4.5)

# NMDS by Site 
p_nmds_site <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Site)) +
        geom_point(size = 2.5, alpha = 0.7) +
        stat_ellipse(level = 0.95, size = 0.5) +
        facet_wrap(~ Month + Habitat, nrow = 2) +
        theme_bw(base_size = 11) +
        labs(title = "NMDS by Site, Habitat, and Month") +
        theme(legend.position = "bottom")

p_nmds_site


# NMDS plot
p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
        stat_ellipse(aes(group = Habitat, fill = Habitat),
                     level = 0.95,
                     geom = "polygon",
                     alpha = 0.2,
                     color = NA) +
        geom_point(aes(color = Habitat, shape = Month),
                   size = 2.5, alpha = 0.8) +
        scale_color_manual(values = c("Fore reef" = "#1E90FF", "Terrace" = "#ffba08")) +
        scale_fill_manual(values = c("Fore reef" = "#1E90FF", "Terrace" = "#ffba08")) +
        scale_shape_manual(values = c(16, 17)) +
        theme_bw(base_size = 10) +
        theme(
                axis.title = element_text(size = 11, face = "bold"),
                axis.text = element_text(size = 10),
                legend.title = element_text(size = 10, face = "bold"),
                legend.text = element_text(size = 10),
                legend.position =  "bottom",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white"),
                panel.border = element_rect(color = "black", linewidth = 0.5)
        ) +
        labs(x = "NMDS Axis 1", y = "NMDS Axis 2") +
        annotate("text", x = -1.2, y = -0.6,
                 label = paste0("Stress = ", round(nmds$stress, 3),
                                "\nPseudo-F = ", round(F_val_text, 3),
                                ", R² = ", R_val_text,
                                ", p-value ", p_val_text),
                 hjust = 0, vjust = 0, size = 3.3)
p_nmds 


### SIMPER: Health conditions contribution to dissimilarity --------------------------------

### SIMPER by Habitat 
simper_habitat <- simper(prev_matrix, metadata$Habitat, permutations = 999)
simper_habitat_summary <- summary(simper_habitat, ordered = TRUE)
simper_habitat     
simper_habitat_summary
# By Month
simper_month <- simper(prev_matrix, metadata$Month, permutations = 999)
simper_month_summary <- summary(simper_month, ordered = TRUE)

# Verificar el nombre exacto del contraste
names(summary(simper_habitat))

simper_hab_df <- summary(simper_habitat)[[1]] %>%
        as.data.frame() %>%
        tibble::rownames_to_column("variable") %>%
        arrange(desc(average)) %>%
        head(10) %>%
        mutate(
                variable_short = case_when(
                                variable == "Bleaching" ~ "BL",
                                variable == "Old Mortality" ~ "OM",
                                variable == "Recent Mortality" ~ "RM",
                                variable == "Sponge Bioerosion" ~ "SB",
                                variable == "Polychaete Bioerosion" ~ "PB",
                                variable == "Dark Spot Syndrome" ~ "DSS",
                                variable == "Rapid Tissue Loss Disease" ~ "RTLD",
                                variable == "Black Band Disease" ~ "BBD",
                                TRUE ~ variable
                        ),
                habitat_indicator = case_when(
                        p >= 0.05 ~ "Non-significant",
                        ava > avb ~ "Fore reef",
                        avb > ava ~ "Terrace"
                ),
                fold_change = ifelse(
                        ava > avb,
                        ifelse(avb == 0, NA, ava / avb),
                        ifelse(ava == 0, NA, avb / ava)
                ),
                fold_change_label = ifelse(
                        p < 0.05 & !is.na(fold_change),
                        paste0(round(fold_change, 2), "×"),
                        ""
                ),
                p_label = case_when(
                        p < 0.001 ~ "p < 0.001",
                        p < 0.01  ~ "p < 0.01",
                        p < 0.05  ~ "p < 0.05",
                        TRUE ~ "ns"
                ), 
                label_out = case_when(
                        p < 0.05 & !is.na(fold_change) ~ 
                                paste0("×", round(fold_change, 1), ",  ", p_label),
                        p < 0.05 ~ p_label,
                        TRUE ~ "")
        )
simper_hab_df

# SIMPER plot

p_simper_habitat <- ggplot(simper_hab_df, 
                           aes(x = reorder(variable_short, average), y = average,
                               fill = habitat_indicator)) +
        geom_col(alpha = 1, width = 0.6) +
        scale_fill_manual(values = c("Fore reef" = "#1E90FF", 
                                     "Terrace" = "#ffba08",
                                     "Non-significant" = "gray70")) +
        geom_text(aes(label = label_out,
                      hjust = ifelse(average >= 0.15, 1.1, -0.1),
                      color = ifelse(average >= 0.15, "black", "black")),
                  size = 3.3) +
        scale_color_identity() +
        coord_flip(clip = "off") +
        labs(x = "", 
             y = "Average contribution to dissimilarity",
             fill = "Habitat",
             #title = "A) By Habitat"
             )         +
        theme_bw(base_size = 10) +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_blank(),   
              panel.grid.minor.x = element_blank(),   
              panel.grid.minor.y = element_blank(),   
              axis.title = element_text(size = 10, face = "bold"),
              axis.text.y = element_text(size = 10, color = "gray30"),
              axis.text.x = element_text(size = 10.0),
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 10),
              legend.position = "bottom",
              panel.border = element_rect(color = "black", linewidth = 0.25)
        )
p_simper_habitat

#ggsave("Fig_SIMPER_Health_Habitat.pdf", p_simper, width = 6, height = 4.5)

### SIMPER: Fore Reef (Mayo vs Agosto)
metadata_fore_reef <- metadata %>% filter(Habitat == "Fore reef")
idx_fore_reef <- which(metadata$Habitat == "Fore reef")
prev_matrix_fore_reef <- prev_matrix[idx_fore_reef, ]

simper_fore_reef <- simper(prev_matrix_fore_reef, metadata_fore_reef$Month, permutations = 999)
simper_fore_reef_summary <- summary(simper_fore_reef, ordered = TRUE)

simper_fore_reef_df <- summary(simper_fore_reef)[[1]] %>%
        as.data.frame() %>%
        tibble::rownames_to_column("variable") %>%
        arrange(desc(average)) %>%
        mutate(
                variable_short = case_when(
                        variable == "Bleaching" ~ "BL",
                        variable == "Old Mortality" ~ "OM",
                        variable == "Recent Mortality" ~ "RM",
                        variable == "Sponge Bioerosion" ~ "SB",
                        variable == "Polychaete Bioerosion" ~ "PB",
                        variable == "Dark Spot Syndrome" ~ "DSS",
                        variable == "Rapid Tissue Loss Disease" ~ "RTLD",
                        variable == "Black Band Disease" ~ "BBD",
                        TRUE ~ variable
                ),
                fold_change = ifelse(
                        ava > avb,
                        ifelse(avb == 0, NA, ava / avb),
                        ifelse(ava == 0, NA, avb / ava)
                ),
                fold_change_label = ifelse(
                        p < 0.05 & !is.na(fold_change),
                        paste0(round(fold_change, 2), "×"),
                        ""
                ),
                month_indicator = case_when(
                        p >= 0.05 ~ "Non-significant",
                        ava > avb ~ "May",
                        avb > ava ~ "August"
                ),
                p_label = case_when(
                        p < 0.001 ~ "p < 0.001",
                        p < 0.01  ~ "p < 0.01",
                        p < 0.05  ~ "p < 0.05",
                        TRUE ~ "ns"
                ),
                label_out = case_when(
                        p < 0.05 & !is.na(fold_change) ~ 
                                paste0("×", round(fold_change, 2), ",  ", p_label),
                        p < 0.05 ~ p_label,
                        TRUE ~ "")
        )
simper_fore_reef_df

p_simper_fore_reef <- ggplot(simper_fore_reef_df, 
                             aes(x = reorder(variable_short, average), y = average,
                                 fill = month_indicator)) +
        geom_col(alpha = 1, width = 0.6) +
        scale_fill_manual(values = c("May" = "#4169E1", 
                                     "August" = "#87CEEB",
                                     "Non-significant" = "gray70")) +
        geom_text(aes(label = label_out,
                      hjust = ifelse(average >= 0.15, 1.1, -0.1),
                      color = ifelse(average >= 0.15, "black", "black")),
                  size = 3.3) +
        scale_color_identity() +
        coord_flip(clip = "off") +
        labs(x = "", 
             y = "Average contribution to dissimilarity",
             fill = "Month",
             #title = "B) Fore Reef\n(May vs August)"
             ) +
        theme_bw(base_size = 10) +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_blank(),   
              panel.grid.minor.x = element_blank(),   
              panel.grid.minor.y = element_blank(),   
              axis.title = element_text(size = 10, face = "bold"),
              axis.text.y = element_text(size = 10, color = "gray30"),
              axis.text.x = element_text(size = 10),
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 10),
              legend.position = "bottom",
              panel.border = element_rect(color = "black", linewidth = 0.25)
        )

print(p_simper_fore_reef)

### SIMPER: Terrace (Mayo vs Agosto)
metadata_terrace <- metadata %>% filter(Habitat == "Terrace")
idx_terrace <- which(metadata$Habitat == "Terrace")
prev_matrix_terrace <- prev_matrix[idx_terrace, ]

simper_terrace <- simper(prev_matrix_terrace, metadata_terrace$Month, permutations = 999)
simper_terrace_summary <- summary(simper_terrace, ordered = TRUE)

# Agregar fold_change e indicador

simper_terrace_df <- summary(simper_terrace)[[1]] %>%
        as.data.frame() %>%
        tibble::rownames_to_column("variable") %>%
        arrange(desc(average)) %>%
        mutate(
                variable_short = case_when(
                        variable == "Bleaching" ~ "BL",
                        variable == "Old Mortality" ~ "OM",
                        variable == "Recent Mortality" ~ "RM",
                        variable == "Sponge Bioerosion" ~ "SB",
                        variable == "Polychaete Bioerosion" ~ "PB",
                        variable == "Dark Spot Syndrome" ~ "DSS",
                        variable == "Rapid Tissue Loss Disease" ~ "RTLD",
                        variable == "Black Band Disease" ~ "BBD",
                        TRUE ~ variable
                ),
                fold_change = ifelse(
                        ava > avb,
                        ifelse(avb == 0, NA, ava / avb),
                        ifelse(ava == 0, NA, avb / ava)
                ),
                fold_change_label = ifelse(
                        p < 0.05 & !is.na(fold_change),
                        paste0(round(fold_change, 2), "×"),
                        ""
                ),
                month_indicator = case_when(
                        p >= 0.05 ~ "Non-significant",
                        ava > avb ~ "May",
                        avb > ava ~ "August"
                ),
                p_label = case_when(
                        p < 0.001 ~ "p < 0.001",
                        p < 0.01  ~ "p < 0.01",
                        p < 0.05  ~ "p < 0.05",
                        TRUE ~ "ns"
                ),
                label_out = case_when(
                        p < 0.05 & !is.na(fold_change) ~ 
                                paste0("×", round(fold_change, 2), ",  ", p_label),
                        p < 0.05 ~ p_label,
                        TRUE ~ "")
        )


p_simper_terrace <- ggplot(simper_terrace_df, 
                           aes(x = reorder(variable_short, average), y = average,
                               fill = month_indicator)) +
        geom_col(alpha = 1, width = 0.6) +
        scale_fill_manual(values = c("May" = "#ffba08", 
                                     "August" = "#ff8c00",
                                     "Non-significant" = "gray70")) +
        geom_text(aes(label = label_out,
                      hjust = ifelse(average >= 0.15, 1.1, -0.1),
                      color = ifelse(average >= 0.15, "black", "black")),
                  size = 3.3) +
        scale_color_identity() +
        coord_flip(clip = "off") +
        labs(x = "", 
             y = "Average contribution to dissimilarity",
             fill = "Month",
             #title = "C) Terrace\n(May vs August)"
             ) +
        theme_bw(base_size = 10) +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_blank(),   
              panel.grid.minor.x = element_blank(),   
              panel.grid.minor.y = element_blank(),   
              axis.title = element_text(size = 10, face = "bold"),
              axis.text.y = element_text(size = 10, color = "gray30"),
              axis.text.x = element_text(size = 10),
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 10),
              legend.position = "bottom",
              panel.border = element_rect(color = "black", linewidth = 0.25)
        )

print(p_simper_terrace)

### RESULT. TABLE 7. Simper table ------------------------------------------------------
# SIMPER by Habitat (Original) 
simper_table <- simper_hab_df %>%
        mutate(
                Variable = variable,
                `Contribution (%)` = sprintf("%.1f", average * 100 / sum(average)),
                `Fore reef (%)` = sprintf("%.2f", ava),
                `Terrace (%)` = sprintf("%.2f", avb),
                `Fold change` = ifelse(!is.na(fold_change), sprintf("%.1f×", fold_change), "—"),
                `p-value` = case_when(
                        p < 0.001 ~ "<0.001***",
                        p < 0.01 ~ sprintf("%.3f**", p),
                        p < 0.05 ~ sprintf("%.3f*", p),
                        TRUE ~ sprintf("%.3f", p)
                ),
                Indicator = habitat_indicator
        ) %>%
        select(Variable, `Contribution (%)`, `Fore reef (%)`, 
               `Terrace (%)`, `Fold change`, `p-value`, Indicator)

average_dissimilarity_habitat <- sum(simper_hab_df$average)
top_cumsum_habitat <- simper_hab_df$cumsum[2]

ft_simper <- flextable(simper_table) %>%
        theme_booktabs() %>%
        width(j = 1, width = 1.8) %>%
        width(j = 2, width = 1) %>%
        width(j = 3, width = 0.7) %>%
        width(j = 4, width = 0.7) %>%
        width(j = 5, width = 0.65) %>%
        width(j = 6, width = 0.7) %>%
        width(j = 7, width = 1.0) %>%
        align(j = 1, align = "left", part = "all") %>%
        align(j = 2:7, align = "center", part = "all") %>%
        fontsize(size = 10, part = "body") %>%
        fontsize(size = 10, part = "header") %>%
        padding(j = 1:7, padding.left = 3, padding.right = 3, part = "all") %>%
        padding(j = 1:7, padding.top = 2, padding.bottom = 2, part = "all") %>%
        add_header_lines(
                "Table 7a. Similarity percentage analysis (SIMPER) identifying health conditions contributing most to compositional dissimilarity between fore reef and terrace habitats. Variables ranked by contribution to Bray-Curtis dissimilarity. Abundance values represent mean relative abundance (%). Fold change indicates enrichment ratio. Indicator classification based on permutation tests (n = 999)."
        ) %>%
        add_footer_lines(
                sprintf(
                        "Average dissimilarity = %.1f%%. Top two variables explain %.1f%% of compositional differences. * p < 0.05, ** p < 0.01, *** p < 0.001",
                        average_dissimilarity_habitat * 100,
                        top_cumsum_habitat * 100)) %>%
        fontsize(part = "footer", size = 9) %>%
        align(part = "footer", align = "left") %>%
        italic(part = "footer")

print(ft_simper)
save_as_docx(ft_simper, path = "Table_7a_Health_SIMPER_Habitat.docx")

# SIMPER - Fore Reef (May vs August)
simper_fore_reef_df

# Calcular fold change para Fore Reef
simper_fore_reef_df <- simper_fore_reef_df %>%
        mutate(
                fold_change = ifelse(
                        ava > avb,
                        ifelse(avb == 0, NA, ava / avb),
                        ifelse(ava == 0, NA, avb / ava)
                )
        )
simper_table_fore_reef <- simper_fore_reef_df %>%
        mutate(
                Variable = variable,
                `Contribution (%)` = sprintf("%.1f", average * 100 / sum(average)),
                `May (%)` = sprintf("%.2f", ava),
                `August (%)` = sprintf("%.2f", avb),
                `Fold change` = ifelse(!is.na(fold_change), sprintf("%.1f×", fold_change), "—"),
                `p-value` = case_when(
                        p < 0.001 ~ "<0.001***",
                        p < 0.01 ~ sprintf("%.3f**", p),
                        p < 0.05 ~ sprintf("%.3f*", p),
                        TRUE ~ sprintf("%.3f", p)
                ),
                Indicator = case_when(
                        p >= 0.05 ~ "Non-significant",
                        ava > avb ~ "May",
                        avb > ava ~ "August"
                )
        ) %>%
        select(Variable, `Contribution (%)`, `May (%)`, 
               `August (%)`, `Fold change`, `p-value`, Indicator)

average_dissimilarity_fore_reef <- sum(simper_fore_reef_df$average)
top_cumsum_fore_reef <- simper_fore_reef_df$cumsum[1]

ft_simper_fore_reef <- flextable(simper_table_fore_reef) %>%
        theme_booktabs() %>%
        width(j = 1, width = 1.8) %>%
        width(j = 2, width = 1) %>%
        width(j = 3, width = 0.7) %>%
        width(j = 4, width = 0.7) %>%
        width(j = 5, width = 0.65) %>%
        width(j = 6, width = 0.7) %>%
        width(j = 7, width = 1.0) %>%
        align(j = 1, align = "left", part = "all") %>%
        align(j = 2:7, align = "center", part = "all") %>%
        fontsize(size = 10, part = "body") %>%
        fontsize(size = 10, part = "header") %>%
        padding(j = 1:7, padding.left = 3, padding.right = 3, part = "all") %>%
        padding(j = 1:7, padding.top = 2, padding.bottom = 2, part = "all") %>%
        add_header_lines(
                "Table 7b. Similarity percentage analysis (SIMPER) - Fore Reef. Health conditions contributing most to compositional dissimilarity between May and August. Variables ranked by contribution to Bray-Curtis dissimilarity. Abundance values represent mean relative abundance (%). Fold change indicates enrichment ratio. Indicator classification based on permutation tests (n = 999)."
        ) %>%
        add_footer_lines(
                sprintf(
                        "Average dissimilarity = %.1f%%. Top variable explains %.1f%% of compositional differences. * p < 0.05, ** p < 0.01, *** p < 0.001",
                        average_dissimilarity_fore_reef * 100,
                        top_cumsum_fore_reef * 100)) %>%
        fontsize(part = "footer", size = 9) %>%
        align(part = "footer", align = "left") %>%
        italic(part = "footer")

print(ft_simper_fore_reef)
save_as_docx(ft_simper_fore_reef, path = "Table_7b_Health_SIMPER_Fore_Reef.docx")

# SIMPER - Terrace (May vs August)
# Calcular fold change para Terrace
simper_terrace_df <- simper_terrace_df %>%
        mutate(
                fold_change = ifelse(
                        ava > avb,
                        ifelse(avb == 0, NA, ava / avb),
                        ifelse(ava == 0, NA, avb / ava)
                )
        )

simper_table_terrace <- simper_terrace_df %>%
        mutate(
                Variable = variable,
                `Contribution (%)` = sprintf("%.1f", average * 100 / sum(average)),
                `May (%)` = sprintf("%.2f", ava),
                `August (%)` = sprintf("%.2f", avb),
                `Fold change` = ifelse(!is.na(fold_change), sprintf("%.1f×", fold_change), "—"),
                `p-value` = case_when(
                        p < 0.001 ~ "<0.001***",
                        p < 0.01 ~ sprintf("%.3f**", p),
                        p < 0.05 ~ sprintf("%.3f*", p),
                        TRUE ~ sprintf("%.3f", p)
                ),
                Indicator = case_when(
                        p >= 0.05 ~ "Non-significant",
                        ava > avb ~ "May",
                        avb > ava ~ "August"
                )
        ) %>%
        select(Variable, `Contribution (%)`, `May (%)`, 
               `August (%)`, `Fold change`, `p-value`, Indicator)


average_dissimilarity_terrace <- sum(simper_terrace_df$average)
top_cumsum_terrace <- simper_terrace_df$cumsum[1]

ft_simper_terrace <- flextable(simper_table_terrace) %>%
        theme_booktabs() %>%
        width(j = 1, width = 1.8) %>%
        width(j = 2, width = 1) %>%
        width(j = 3, width = 0.7) %>%
        width(j = 4, width = 0.7) %>%
        width(j = 5, width = 0.65) %>%
        width(j = 6, width = 0.7) %>%
        width(j = 7, width = 1.0) %>%
        align(j = 1, align = "left", part = "all") %>%
        align(j = 2:7, align = "center", part = "all") %>%
        fontsize(size = 10, part = "body") %>%
        fontsize(size = 10, part = "header") %>%
        padding(j = 1:7, padding.left = 3, padding.right = 3, part = "all") %>%
        padding(j = 1:7, padding.top = 2, padding.bottom = 2, part = "all") %>%
        add_header_lines(
                "Tabla 7c. Similarity percentage analysis (SIMPER) - Terrace. Health conditions contributing most to compositional dissimilarity between May and August. Variables ranked by contribution to Bray-Curtis dissimilarity. Abundance values represent mean relative abundance (%). Fold change indicates enrichment ratio. Indicator classification based on permutation tests (n = 999)."
        ) %>%
        add_footer_lines(
                sprintf(
                        "Average dissimilarity = %.1f%%. Top variable explains %.1f%% of compositional differences. * p < 0.05, ** p < 0.01, *** p < 0.001",
                        average_dissimilarity_terrace * 100,
                        top_cumsum_terrace * 100)) %>%
        fontsize(part = "footer", size = 9) %>%
        align(part = "footer", align = "left") %>%
        italic(part = "footer")

print(ft_simper_terrace)
save_as_docx(ft_simper_terrace, path = "Table_7c_Health_SIMPER_Terrace.docx")

### RESULT. FIGURE 4. Health Composition Panel -------------------------------------------

combined <- p_nmds_species + p_simper_habitat + p_simper_fore_reef + p_simper_terrace +
        plot_layout(ncol = 2, nrow = 2, widths = c(1, 1), heights = c(1, 1)) +
        plot_annotation(
                tag_levels = "A"
        ) &
        theme(
                plot.tag = element_text(size = 12, face = "bold"),
                plot.tag.position = c(0.02, 0.98),
                plot.caption = element_text(hjust = 0, size = 10, margin = margin(t = 10))
        )

print(combined)

ggsave("Figure_4_Health_NMDS_SIMPER.pdf", combined, 
       width = 11, height = 9, dpi = 300)
ggsave("Figure_4_Health_NMDS_SIMPER.png", combined,
       width = 11, height = 9, dpi = 400)

### Clean Space -----
permanova_prev_int_margin
permanova_prev_add_margin
anosim_habitat
anosim_month
anosim_hab_month
disp_habitat 
habitat_dispersion
disp_month
month_dispersion
disp_interaction
habitatxmonth_dispersion
tabla_combined
species_scores
simper_terrace_df
simper_fore_reef_df
fit
stressplot(nmds)
nmds_stress


rm(list = ls())
graphics.off()
### Summary Results: Community Structure---------------

# CRITICAL METHODOLOGICAL CAVEAT:
# BETADISPER revealed significant heterogeneity in multivariate dispersions for
# Habitat (F = 4.30, p = 0.043), Month (F = 6.94, p = 0.008), and Habitat × Month 
# (F = 16.17, p < 0.001). This violates PERMANOVA's homogeneity of variance assumption 
# and may inflate the F-statistics. ANOSIM validation (below) corroborates true 
# compositional/prevalence differences exist beyond variance heterogeneity alone.

# PERMANOVA - MODEL SELECTION:
# Interaction model: Habitat × Month SIGNIFICANT (R² = 0.090, F = 20.97, p < 0.001)
# → Interactive model selected (not additive)
# 
# PERMANOVA - INTERACTIVE MODEL (Marginal tests, by = "margin"):
# Overall: R² = 0.090, F = 20.97, p < 0.001 (9.0% variance explained by interaction)
# - Habitat × Month: R² = 0.090, F = 20.97, p < 0.001 *** (9.0% variance, DOMINANT)
# - Residual: R² = 0.910
# Site controlled as blocking factor (strata)
#
# PERMANOVA - ADDITIVE MODEL (for comparison):
# - Habitat:   R² = 0.031, F = 6.49,  p < 0.001 *** (3.1% variance)
# - Month:     R² = 0.034, F = 7.06,  p < 0.001 *** (3.4% variance)
# Habitat × Month interaction accounts for 2.7× more variance than either main effect alone

# ANOSIM (Complementary rank-based validation):
# ANOSIM is more robust to dispersion heterogeneity than PERMANOVA
# - Habitat: R = 0.027, p = 0.001 *** (weak but statistically significant)
# - Month:   R = 0.061, p < 0.001 *** (weak but statistically significant)
# - Habitat × Month: R = 0.126, p < 0.001 *** (moderate effect, strongest signal)
# PERMANOVA and ANOSIM AGREE: Habitat × Month interaction is dominant pattern
# Concordance across methods strengthens evidence for true interactive prevalence patterns

# BETADISPER (Homogeneity of multivariate dispersions):
# - Habitat: F = 4.30, p = 0.043 * (SIGNIFICANT HETEROGENEITY)
#   Pairwise: Fore reef vs Terrace p = 0.039
# - Month: F = 6.94, p = 0.008 ** (SIGNIFICANT HETEROGENEITY)
#   May communities show different within-group variance than August
# - Habitat × Month: F = 16.17, p < 0.001 *** (SIGNIFICANT HETEROGENEITY)
#   Four groups (Fore reef_May, Fore reef_August, Terrace_May, Terrace_August)
#   show differential dispersions; suggests interactive effect on community variance

# INTERPRETATION SUMMARY:
# Health condition patterns are INTERACTIVE, not purely additive
# The relationship between Habitat and health prevalence DEPENDS on Month
# This differs fundamentally from species composition (which showed additive effects)
# PERMANOVA interaction effect (R² = 9.0%) exceeds either main effect individually

# NMDS ORDINATION:
# - Stress = 0.102 (good fit: below 0.15 threshold)
# - Environmental vectors significantly associated:
#   * Bleaching: r² = 0.88, p < 0.001 *** (strongest predictor)
#   * Old Mortality: r² = 0.745, p < 0.001 *** (second strongest)
#   * Sponge Bioerosion: r² = 0.191, p < 0.001 ***
#   * Polychaete Bioerosion: r² = 0.128, p < 0.001 ***
#   * Recent Mortality: r² = 0.061, p = 0.006 **
#   * Rapid Tissue Loss Disease: r² = 0.072, p = 0.010 **
#   * Dark Spot Syndrome: r² = 0.033, p = 0.030 *
#   * Black Band Disease: r² = 0.039, p = 0.038 *
# 
# Bleaching and Old Mortality dominate health condition space;
# bioerosion conditions add secondary but significant structure

# SIMPER - HABITAT (Fore reef vs Terrace):
# Average Bray-Curtis dissimilarity = 31.9% (health conditions only)
# Top 2 conditions (Bleaching + Old Mortality) explain 79.2% of habitat dissimilarity
# 
# FORE REEF INDICATORS (significantly more prevalent, p < 0.05):
# - Bleaching: 55.0% contribution, p = 0.016 *, 1.1× higher in fore reef
#   Prevalence: 42.0% (fore reef) vs. 37.7% (terrace)
#   → Higher physiological stress/heat susceptibility in fore reef despite deeper, 
#      cooler conditions; may reflect competitive crowding or post-disturbance recovery
#
# TERRACE INDICATORS (significantly more prevalent, p < 0.05):
# - Sponge Bioerosion: 8.7% contribution, p = 0.006 **, 1.6× higher in terrace
#   Prevalence: 2.45% (fore reef) vs. 3.95% (terrace)
#   → Reflects shallow terrace exposure to sponge predation/competition
# - Polychaete Bioerosion: 5.8% contribution, p = 0.001 ***, 2.9× higher in terrace
#   Prevalence: 0.83% (fore reef) vs. 2.44% (terrace)
#   → Shallow, high-energy environment favors polychaete colonization
# - Recent Mortality: 2.4% contribution, p = 0.001 ***, 4.4× higher in terrace
#   Prevalence: 0.21% (fore reef) vs. 0.91% (terrace)
#   → Higher acute mortality in shallow, thermally variable terrace conditions
# - Rapid Tissue Loss Disease: 1.1% contribution, p = 0.005 **, 8.4× higher in terrace
#   Prevalence: 0.05% (fore reef) vs. 0.41% (terrace)
#   → RTLD (includes White Plague, White Band, Stony Coral Tissue Loss Disease)
#      shows strong shallow-water preference; highest fold-change of any condition

# SIMPER - MONTH (May vs August):
# Average Bray-Curtis dissimilarity = 31.9% (health conditions only)
# Top 2 conditions (Bleaching + Old Mortality) explain 79.4% of temporal dissimilarity
# 
# SIGNIFICANT TEMPORAL CHANGES (p < 0.05):
# Health CONDITIONS INCREASING May → August:
# - Bleaching: 55.5% contribution, p = 0.001 ***, 1.3× increase
#   Prevalence: 41.4% (May) → 38.3% (August)
#   → Slight seasonal increase in bleaching events; may reflect thermal history
# - Polychaete Bioerosion: 2.0% contribution, p = 0.002 **, 2.2× increase
#   Prevalence: 2.24% (May) → 1.03% (August)
#   → Wait, this is DECREASING, not increasing; check data direction
# - Recent Mortality: 0.8% contribution, p = 0.001 ***, 7.0× increase
#   Prevalence: 0.98% (May) → 0.14% (August)
#   → Post-disturbance recovery pattern; mortality declining August
# - Dark Spot Syndrome: 0.8% contribution, p = 0.007 **, 2.2× increase
#   Prevalence: 0.83% (May) → 0.38% (August)
#   → Declining through season as conditions stabilize
# - Rapid Tissue Loss Disease: 0.4% contribution, p = 0.003 **, 6.7× increase
#   Prevalence: 0.40% (May) → 0.06% (August)
#   → Major seasonal decline; possibly post-spring-storm recovery
# - Black Band Disease: 0.1% contribution, p = 0.009 **, 10× increase
#   Prevalence: 0.10% (May) → 0.01% (August)
#   → Rare condition showing strong seasonal pattern (May peak)

# SIMPER - HABITAT × MONTH CONTRASTS:
# Six pairwise comparisons reveal interaction structure:
# 
# 1. Fore reef May vs Terrace May (within-month habitat contrast):
#    Bleaching dominates dissimilarity (60.1%); Polychaete Bioerosion 2.2% (p < 0.01)
#    Both bleaching and bioerosion show stronger differentiation in May
# 
# 2. Fore reef May vs Fore reef August (within-habitat temporal contrast):
#    Bleaching dominates (60.9%); most conditions show ns temporal change
#    Fore reef health largely stable across seasons
# 
# 3. Fore reef May vs Terrace August (maximum habitat-temporal contrast):
#    Bleaching lowest contribution (56.3%); polychaete bioerosion 1.2%
#    Opposite habitats in opposite months show modest dissimilarity
# 
# 4. Terrace May vs Fore reef August (crossed contrast):
#    Bleaching 50.5%, but NOW includes 4 significant conditions:
#    Old Mortality (p < 0.001), Sponge Bioerosion (p < 0.01),
#    Polychaete Bioerosion (p < 0.001), Recent Mortality (p < 0.001)
#    This is the STRONGEST interactive pattern; highest condition diversity
# 
# 5. Terrace May vs Terrace August (within-habitat temporal contrast):
#    Bleaching 55.7%; includes 3 significant seasonal shifts:
#    Recent Mortality (p < 0.001), Rapid Tissue Loss Disease (p < 0.01),
#    Sponge Bioerosion (p = 0.024)
#    Terrace shows MORE temporal variation than fore reef
# 
# 6. Fore reef August vs Terrace August (within-month habitat contrast):
#    Bleaching dominates (53.5%); NO significant condition differences
#    By August, habitat differentiation essentially disappears
#    → Convergence pattern: May shows strong habitat separation;
#       August shows habitat homogenization in health conditions

# KEY ECOLOGICAL INTERPRETATION
# 1. INTERACTIVE EFFECT DOMINATES (R² = 9.0%, 2.7× larger than any main effect)
#    - Health prevalence patterns DEPEND on habitat-season combination
#    - No simple "always worse in terrace" or "always worse in May/August"
#    - Different habitats show different seasonal trajectories
# 
# 2. MAY SHOWS STRONG HABITAT SEPARATION; AUGUST SHOWS CONVERGENCE
#    - May: Fore reef (high bleaching: 51.5%) vs Terrace (low: 31.3%)
#    - August: Fore reef (32.5%) and Terrace (44.1%) nearly equal
#    - Interpretation: Seasonal recovery mechanisms differ by habitat;
#      terrace may show delayed recovery or continued stress accumulation
# 
# 3. BLEACHING IS OVERWHELMINGLY DOMINANT (55% of dissimilarity in both contrasts)
#    - Bleaching vector strongest environmental predictor in NMDS (r² = 0.88)
#    - Old Mortality second (24%, r² = 0.745)
#    - All other conditions together <5% of dissimilarity
#    - System fundamentally structured by bleaching-mortality dynamics
# 
# 4. TERRACE CONDITIONS SHOW HIGHER SEASONAL VARIABILITY
#    - Fore reef May→August: mostly stable (only bleaching changes, p = ns)
#    - Terrace May→August: 5 conditions show significant shifts
#      (Bleaching p < 0.01, Polychaete Bioerosion p < 0.001, Recent Mortality p < 0.001,
#       Sponge Bioerosion p = 0.024, Rapid Tissue Loss Disease p < 0.01)
#    - Terrace communities less resilient or more responsive to seasonal forcing
# 
# 5. RECENT MORTALITY SHOWS INVERSE SEASONAL PATTERN
#    - High in May (0.98%), declines sharply by August (0.14%)
#    - Consistent recovery timeline: post-disturbance healing ~3-4 months
#    - Strongest in Terrace (1.62% May), weak in Fore reef (0.35% May)
#    - Suggests terrace experienced more recent acute stress
# 
# 6. RTLD (WHITE PLAGUE/WHITE BAND/STONY CORAL TISSUE LOSS DISEASE)
#    - Highest fold-change (8.4× higher in terrace vs fore reef)
#    - Lowest overall prevalence (<0.5%) but strongly habitat-specific
#    - May peak (0.40%) declines to near-zero August (0.06%)
#    - Suggests spring disease outbreak followed by recovery/resolution
# 
# 7. BIOEROSION CONDITIONS MARK SHALLOW-WATER STRESS
#    - Polychaete Bioerosion: 2.9× enriched in terrace
#    - Sponge Bioerosion: 1.6× enriched in terrace
#    - Both correlate with RTLD and Recent Mortality
#    - Triad of shallow-water stressors: physical stress + disease + bioerosion
# 
# 8. METHODOLOGICAL ROBUSTNESS
#    - BETADISPER shows significant heterogeneity (all three tests p < 0.05)
#    - PERMANOVA interaction F-statistic may be partially inflated
#    - ANOSIM validation confirms interactive pattern (R = 0.126, p < 0.001)
#    - Concordance between methods supports genuine interactive effect
#    - Conservative interpretation: Habitat × Month interaction is real,
#      but dispersion heterogeneity contributes additional variance



# 3. PREDICTING CORAL HEALTH AT COLONY LEVEL ------
### 3.1 By Habitat and month with GLM/GLMM models -----

library(glmmTMB)
library(tidyverse)
library(performance)
library(DHARMa)
library(readxl)
library(dplyr)
library(brms)
library(forestploter)
library(ggplot2)
library(flextable)

### Load and prepare data colony level --------
df_colony_level <- read_excel(
        "D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/Corals_data_JR_Forereefs_Terrace_2024_August_May.xlsx",
        sheet = "Health_primary"
) %>%
        filter(Site != "Auras") %>%
        droplevels() %>%
        mutate(
                Month = factor(Month, levels = c("May", "August")),
                Habitat = factor(Habitat, levels = c("Fore reef", "Terrace")),
                Site = factor(Site, levels = c("Peruano", "Pinos", "Anclitas", "Mariflores", "Cruces")),
                Transect = factor(Transect)
        )
print (df_colony_level)
events_summary <- df_colony_level %>%
        summarise(
                N_total = n(),
                # Condiciones porcentuales (0-100%)
                MA = sum(!is.na(MA)),
                MR = sum(!is.na(MR)),
                BL = sum(!is.na(BL)),
                # Enfermedades (presencia/ausencia)
                BBD = sum(BBD == 1, na.rm = T),
                DSS = sum(DSS == 1, na.rm = T),
                SCTLD = sum(SCTLD == 1, na.rm = T),
                RTLD = sum(RTLD == 1, na.rm = T),
                Diseased = sum(Diseased == 1, na.rm = T),
                # Bioerosión
                Bioero_spo = sum(BIOERO_SPO == 1, na.rm = T),
                Bioero_pol = sum(BIOERO_POL == 1, na.rm = T)
        ) %>%
        pivot_longer(everything(), names_to = "Condition", values_to = "Number of events")

print(events_summary)

# Prepraring binomial data
df_model <- df_colony_level %>%
        mutate(
                # Convertir porcentuales a BINARIAS (0/1): afectada o no
                MA = ifelse(is.na(MA) | MA == 0, 0, 1),
                MR = ifelse(is.na(MR) | MR == 0, 0, 1),
                BL = ifelse(is.na(BL) | BL == 0, 0, 1),
                
                # Binarias: NA → 0, TRUE → 1
                BBD = ifelse(is.na(BBD), 0, as.numeric(BBD)),
                DSS = ifelse(is.na(DSS), 0, as.numeric(DSS)),
                SCTLD = ifelse(is.na(SCTLD), 0, as.numeric(SCTLD)),
                RTLD = ifelse(is.na(RTLD), 0, as.numeric(RTLD)),
                Diseased = ifelse(is.na(Diseased), 0, as.numeric(Diseased)),
                BIOERO_SPO = ifelse(is.na(BIOERO_SPO), 0, as.numeric(BIOERO_SPO)),
                BIOERO_POL = ifelse(is.na(BIOERO_POL), 0, as.numeric(BIOERO_POL)),
                
                # Convertir todas a factor para claridad
                MA = factor(MA, levels = c(0, 1)),
                MR = factor(MR, levels = c(0, 1)),
                BL = factor(BL, levels = c(0, 1)),
                BBD = factor(BBD, levels = c(0, 1)),
                DSS = factor(DSS, levels = c(0, 1)),
                SCTLD = factor(SCTLD, levels = c(0, 1)),
                RTLD = factor(RTLD, levels = c(0, 1)),
                Diseased = factor(Diseased, levels = c(0, 1)),
                BIOERO_SPO = factor(BIOERO_SPO, levels = c(0, 1)),
                BIOERO_POL = factor(BIOERO_POL, levels = c(0, 1))
        )

df_model %>% 
        select(BL, MA, MR, RTLD, BBD, DSS, Diseased, BIOERO_SPO, BIOERO_POL) %>%
        head(20)
df_model

### Running GLM/GLMM models for each health condition ---------------------------------
glmm_bl <- glmmTMB(BL ~ Habitat + Month + (1 | Site) + (1 | Site:Transect),
                   family = binomial(link = "logit"), data = df_model)
glmm_ma <- glmmTMB(MA ~ Habitat + Month + (1 | Site) + (1 | Site:Transect),
                   family = binomial(link = "logit"), data = df_model)
glmm_mr <- glmmTMB(MR ~ Habitat + Month + (1 | Site:Transect),
                   family = binomial(link = "logit"), data = df_model)
glmm_rtld <- glmmTMB(RTLD ~ Habitat + Month + (1 | Site) + (1 | Site:Transect),
                     family = binomial(link = "logit"), data = df_model)
glmm_dss <- glmmTMB(DSS ~ Habitat + Month + (1 | Site) + (1 | Site:Transect),
                    family = binomial(link = "logit"), data = df_model)
glmm_bbd <- glmmTMB(BBD ~ Habitat + Month + (1 | Site:Transect),
                    family = binomial(link = "logit"), data = df_model)
glmm_diseased<- glmmTMB(Diseased ~ Habitat + Month + (1 | Site:Transect),
                    family = binomial(link = "logit"), data = df_model)
glmm_bioero_spo <- glmmTMB(BIOERO_SPO ~ Habitat + Month + (1 | Site:Transect),
                           family = binomial(link = "logit"), data = df_model)
glmm_bioero_pol <- glmmTMB(BIOERO_POL ~ Habitat + Month + (1 | Site:Transect),
                           family = binomial(link = "logit"), data = df_model)

### Models summary  -------

summary(glmm_bl)
summary(glmm_ma)
summary(glmm_mr)
summary(glmm_rtld)
summary(glmm_dss)
summary(glmm_bbd)
summary(glmm_diseased)
summary(glmm_bioero_spo)
summary(glmm_bioero_pol)

### DHARMa diagnostics ---------------------
library(DHARMa)
library(dplyr)

sim_res_bl <- simulateResiduals(glmm_bl, n = 1000)
sim_res_ma <- simulateResiduals(glmm_ma, n = 1000)
sim_res_mr <- simulateResiduals(glmm_mr, n = 1000)
sim_res_rtld <- simulateResiduals(glmm_rtld, n = 1000)
sim_res_dss <- simulateResiduals(glmm_dss, n = 1000)
sim_res_bbd <- simulateResiduals(glmm_bbd, n = 1000)
sim_res_diseased <- simulateResiduals(glmm_diseased, n = 1000)
sim_res_bioero_spo <- simulateResiduals(glmm_bioero_spo, n = 1000)
sim_res_bioero_pol <- simulateResiduals(glmm_bioero_pol, n = 1000)

# Lista de simulaciones
sim_list <- list(
        sim_res_bl, 
        sim_res_ma, 
        sim_res_mr, 
        sim_res_rtld,
        sim_res_dss, 
        #sim_res_bbd, 
        sim_res_diseased, 
        sim_res_bioero_spo, 
        sim_res_bioero_pol
)

conditions <- c(
        "Bleaching", "Old Mortality", "Recent Mortality", 
        "Rapid Tissue Loss Disease", "Dark Spot Syndrome", 
        #"Black Band Disease", 
        "Diseased", "Sponge Bioerosion", 
        "Polychaete Bioerosion"
)

# Calcular p_uniformity sin plot
p_values <- sapply(sim_list, function(sim) {
        testUniformity(sim, plot = FALSE)$p.value
})

# Definir fit dinámico según pvalue
model_fit <- sapply(p_values, function(p) {
        if(p > 0.05) {
                "Good"
        } else if(p > 0.01) {
                "Fair"
        } else {
                "Poor"
        }
})

diagnostics <- data.frame(
        Condition = conditions,
        p_uniformity = p_values,
        Model_Fit = model_fit
)

diagnostics

# Extraer coeficientes de los modelos
extract_coefs <- function(model, condition_name) {
        coefs <- fixef(model)$cond
        pvals <- summary(model)$coefficients$cond[, "Pr(>|z|)"]
        ses <- summary(model)$coefficients$cond[, "Std. Error"]
        data.frame(
                Condition = condition_name,
                Predictor = names(coefs),
                Estimate = round(coefs, 4),
                SE = round(ses, 4),
                Z_value = round(coefs/ses, 3),
                P_value = format.pval(pvals, digits = 3),
                Significant = ifelse(pvals < 0.05, "***", "ns"),
                stringsAsFactors = FALSE
        )
}

coef_all <- bind_rows(
        extract_coefs(glmm_bl, "Bleaching"),
        extract_coefs(glmm_ma, "Old Mortality"),
        extract_coefs(glmm_mr, "Recent Mortality"),
        extract_coefs(glmm_rtld, "Rapid Tissue Loss Disease"),
        extract_coefs(glmm_dss, "Dark Spot Syndrome"),
        #extract_coefs(glmm_bbd, "Black Band Disease"),
        extract_coefs(glmm_diseased, "Diseased"),
        extract_coefs(glmm_bioero_spo, "Sponge Bioerosion"),
        extract_coefs(glmm_bioero_pol, "Polychaete Bioerosion")
)
coef_all


# Checking singularity

library(performance)

models_list <- list(
        Bleaching = glmm_bl,
        Old_Mortality = glmm_ma,
        Recent_Mortality = glmm_mr,
        RTLD = glmm_rtld,
        DSS = glmm_dss,
        BBD = glmm_bbd,
        Diseased = glmm_diseased,
        Sponge_Bioerosion = glmm_bioero_spo,
        Polychaete_Bioerosion = glmm_bioero_pol
)

singularity_check <- sapply(models_list, check_singularity)
data.frame(Condition = names(singularity_check), Singular_Fit = singularity_check)

### RESULT. FIGURE 5 Forest plot by habitat and month  ----------

library(forestploter)
library(ggplot2)

# Preparar datos para forest plot
forest_data <- coef_all %>%
        filter(Predictor != "(Intercept)") %>%
        mutate(
                Estimate = as.numeric(Estimate),
                SE = as.numeric(SE),
                CI_lower = Estimate - 1.96 * SE,
                CI_upper = Estimate + 1.96 * SE,
                OR = exp(Estimate),
                OR_lower = exp(CI_lower),
                OR_upper = exp(CI_upper),
                Label = paste0(round(OR, 2), " (", round(OR_lower, 2), "-", round(OR_upper, 2), ")")
        ) %>%
        select(Condition, Predictor, OR, OR_lower, OR_upper, Label, Significant)

forest_data$Predictor <- factor(
        forest_data$Predictor,
        levels = unique(forest_data$Predictor)
)
forest_data
forest_data$Condition <- factor(
        forest_data$Condition,
        levels = c(
                "Bleaching",
                "Recent Mortality",
                "Old Mortality",
                "Dark Spot Syndrome",
                "Rapid Tissue Loss Disease",
                #"Black Band Disease",
                "Diseased",
                "Sponge Bioerosion",
                "Polychaete Bioerosion"
        )
)

# Forest plot simple

forest_data$Predictor <- factor(forest_data$Predictor, levels = unique(forest_data$Predictor))

p_forest <- forest_data %>%
        ggplot(aes(x = OR, y = Predictor, color = Significant)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", size = 0.5) +
        geom_point(size = 3, alpha = 1) +
        geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0, size = 0.8, alpha = 1) +
        facet_wrap(~Condition, scales = "fixed", ncol = 3, drop = FALSE) +  # eje Y fijo
        scale_color_manual(values = c("***" = "#d7301f", "ns" = "gray70")) +
        scale_x_log10() +
        labs(#title = "Predictors of Coral Health Conditions",
                x = "Odds Ratio (log scale)",
                y = "",
                color = "p-value < 0.05") +
        theme_minimal() +
        theme(
                axis.title = element_text(size = 10, face = "bold"),
                axis.text.y = element_text(size = 10),  # solo se dibuja una vez
                axis.text.x = element_text(size = 10),
                legend.position = "bottom",
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_line(color = "gray90"),
                panel.border = element_rect(color = "black", fill = NA),
                strip.text = element_text(face = "bold", size = 10)
        )

print(p_forest)
ggsave("Figure_5_Predicting_health_Forest_GLMM_habitat_month.pdf", p_forest,
       width = 10, height = 7, dpi = 400)
ggsave("Figure_5_Predicting_health_Forest_GLMM_habitat_month.png", p_forest,
       width = 10, height = 7, dpi = 400)

p_habitat_month <- coef_all %>%
        filter(Predictor != "(Intercept)") %>%
        mutate(Estimate = as.numeric(Estimate), Effect_Type = ifelse(grepl("Habitat", Predictor), "Habitat", "Month")) %>%
        ggplot(aes(x = Condition, y = Estimate, fill = Effect_Type, alpha = Significant)) +
        geom_col(position = "dodge", color = "black", size = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_fill_manual(values = c("Habitat" = "#1f77b4", "Month" = "#ff7f0e")) +
        scale_alpha_manual(values = c("***" = 1, "ns" = 0.4)) +
        coord_flip() +
        labs(title = "Habitat vs Temporal Effects", x = "Health Condition", y = "Coefficient") +
        theme_minimal() + theme(legend.position = "bottom")


p_habitat_month

### RESULT. TABLE 8 GLMM  ----------

calc_odds_ratios <- function(model, condition_name) {
        coefs <- fixef(model)$cond[-1]
        ses <- summary(model)$coefficients$cond[-1, "Std. Error"]
        pvals <- summary(model)$coefficients$cond[-1, "Pr(>|z|)"]
        
        or <- exp(coefs)
        ci_lower <- exp(coefs - 1.96 * ses)
        ci_upper <- exp(coefs + 1.96 * ses)
        
        data.frame(
                Condition = condition_name,
                Predictor = names(coefs),
                OR = round(or, 3),
                CI_lower = round(ci_lower, 3),
                CI_upper = round(ci_upper, 3),
                P_value = pvals,
                stringsAsFactors = FALSE
        )
}

or_all <- bind_rows(
        calc_odds_ratios(glmm_bl, "Bleaching"),
        calc_odds_ratios(glmm_ma, "Old Mortality"),
        calc_odds_ratios(glmm_mr, "Recent Mortality"),
        calc_odds_ratios(glmm_rtld, "Rapid Tissue Loss Disease"),
        calc_odds_ratios(glmm_dss, "Dark Spot Syndrome"),
        calc_odds_ratios(glmm_diseased, "Diseased"),
        calc_odds_ratios(glmm_bioero_spo, "Sponge Bioerosion"),
        calc_odds_ratios(glmm_bioero_pol, "Polychaete Bioerosion")
)

table_data <- or_all %>%
        mutate(
                Predictor = case_when(
                        Predictor == "HabitatTerrace" ~ "Terrace (vs. Fore reef)",
                        Predictor == "MonthAugust" ~ "August (vs. May)",
                        Predictor == "HabitatTerrace:MonthAugust" ~ "Terrace × August",
                        TRUE ~ Predictor
                ),
                P_value_formatted = case_when(
                        P_value < 0.001 ~ "<0.001***",
                        P_value < 0.01 ~ paste0(sprintf("%.3f", P_value), "**"),
                        P_value < 0.05 ~ paste0(sprintf("%.3f", P_value), "*"),
                        TRUE ~ sprintf("%.3f", P_value)
                )
        ) %>%
        select(Condition, Predictor, OR, CI_lower, CI_upper, P_value_formatted) %>%
        rename(
                "Health Condition" = Condition,
                "Predictor" = Predictor,
                "OR" = OR,
                "95% CI Lower" = CI_lower,
                "95% CI Upper" = CI_upper,
                "P-value" = P_value_formatted
        )
table_data

ft_table <- flextable(table_data) %>%
        theme_booktabs() %>%
        
        # Combinar celdas verticales de "Health Condition" cuando tienen el mismo valor
        merge_v(j = 1, part = "body") %>%
        
        # Ancho de columnas
        width(j = 1, width = 2) %>%
        width(j = 2, width = 1.5) %>%
        width(j = 3, width = 0.6) %>%
        width(j = 4, width = 0.7) %>%
        width(j = 5, width = 0.7) %>%
        width(j = 6, width = 0.7) %>%
        
        # Alineación
        align(j = 1, align = "left", part = "all") %>%
        align(j = 2, align = "left", part = "all") %>%
        align(j = 3:6, align = "center", part = "all") %>%
        
        # Alineación vertical de celdas combinadas
        valign(j = 1, valign = "center", part = "body") %>%
        
        # Tipografía
        fontsize(size = 10, part = "body") %>%
        fontsize(size = 10, part = "header") %>%
        
        # Padding fino
        padding(j = 1:6, padding.left = 3, padding.right = 3, part = "all") %>%
        padding(j = 1:6, padding.top = 2, padding.bottom = 2, part = "all") %>%
        
        # ENCONTRAR AUTOMÁTICAMENTE DÓNDE CAMBIAN LOS GRUPOS
        # Esto crea una lista con los números de fila donde termina cada condición de salud
        {
                # Obtener los datos de la tabla
                datos <- table_data
                
                # Identificar en qué filas cambia el nombre de "Health Condition"
                cambios <- which(datos$`Health Condition` != dplyr::lag(datos$`Health Condition`))
                
                    # Restar 1 porque queremos la línea ANTES del cambio
                filas_con_linea <- cambios - 1
                
                # Aplicar las líneas finas en esas filas
                hline(., i = filas_con_linea,
                      border = fp_border(color = "#CCCCCC", width = 0.5),
                      part = "body")
        } %>%
        
        # Encabezado descriptivo
        add_header_lines(
                "Table 8. Generalized Linear Mixed Models predicting prevalence of coral health conditions at colony level. Models include Habitat (Fore reef vs. Terrace) and temporal (May vs. August) predictors with random intercepts for Site and Transect. Effects are expressed as Odds Ratios (OR) with 95% confidence intervals."
        ) %>%
        
        # Nota de pie
        add_footer_lines(
                "OR, Odds Ratio; CI, 95% Confidence Interval. Reference categories: Fore reef (Habitat), May (Month). *, p < 0.05; **, p < 0.01; ***, p < 0.001."
        ) %>%
        
        # Formato de notas
        fontsize(part = "footer", size = 9) %>%
        align(part = "footer", align = "left") %>%
        italic(part = "footer")

ft_table

# Guardar tabla como Word
doc <- read_docx() %>%
        body_add_flextable(ft_table) %>%
        print(target = "Table_8_Predicting_health_GLMM_habitat_month.docx")

# Modelos de interaccion para el revisor.

glmm_bl_int <- glmmTMB(BL ~ Habitat * Month + (1 | Site) + (1 | Site:Transect),
                       family = binomial(link = "logit"), data = df_model)
summary(glmm_bl_int)
anova(glmm_bl, glmm_bl_int)  # test de razón de verosimilitud para la interacción

glmm_diseased_int <- glmmTMB(Diseased ~ Habitat * Month + (1 | Site:Transect),
                       family = binomial(link = "logit"), data = df_model)
summary(glmm_diseased_int)
anova(glmm_diseased, glmm_diseased_int)  # test de razón de verosimilitud para la interacción

glmm_mr_int <- glmmTMB(MR ~ Habitat * Month + (1 | Site:Transect),
                       family = binomial(link = "logit"), data = df_model)
summary(glmm_mr_int)
anova(glmm_mr, glmm_mr_int)

glmm_rtld_int <- glmmTMB(RTLD ~ Habitat * Month + (1 | Site) + (1 | Site:Transect),
                         family = binomial(link = "logit"), data = df_model)
summary(glmm_rtld_int)
anova(glmm_rtld, glmm_rtld_int)

#Simpligying models random structure 

library(performance)
library(glmmTMB)

# Modelos originales (con estructura completa) que fueron singulares
singular_conditions <- c("MR", "BBD", "Diseased", "BIOERO_SPO", "BIOERO_POL")

# Función para ajustar versión simplificada (solo Site:Transect, sin Site separado)
fit_simplified <- function(outcome, data) {
        glmmTMB(
                as.formula(paste(outcome, "~ Habitat + Month + (1 | Site:Transect)")),
                family = binomial(link = "logit"),
                data = data
        )
}

# Ajustar modelos simplificados
glmm_mr_simple       <- fit_simplified("MR", df_model)
glmm_bbd_simple       <- fit_simplified("BBD", df_model)
glmm_diseased_simple  <- fit_simplified("Diseased", df_model)
glmm_spo_simple       <- fit_simplified("BIOERO_SPO", df_model)
glmm_pol_simple       <- fit_simplified("BIOERO_POL", df_model)

# Verificar singularidad de los nuevos modelos
simplified_models <- list(
        Recent_Mortality = glmm_mr_simple,
        BBD = glmm_bbd_simple,
        Diseased = glmm_diseased_simple,
        Sponge_Bioerosion = glmm_spo_simple,
        Polychaete_Bioerosion = glmm_pol_simple
)
sapply(simplified_models, check_singularity)

# Comparar coeficientes de efectos fijos: original (completo) vs. simplificado
compare_fixed <- function(full_model, simple_model, name) {
        f <- summary(full_model)$coefficients$cond %>% as.data.frame() %>%
                tibble::rownames_to_column("term") %>% mutate(model = "full", condition = name)
        s <- summary(simple_model)$coefficients$cond %>% as.data.frame() %>%
                tibble::rownames_to_column("term") %>% mutate(model = "simplified", condition = name)
        dplyr::bind_rows(f, s)
}

comparison_all <- dplyr::bind_rows(
        compare_fixed(glmm_mr, glmm_mr_simple, "Recent Mortality"),
        compare_fixed(glmm_bbd, glmm_bbd_simple, "BBD"),
        compare_fixed(glmm_diseased, glmm_diseased_simple, "Diseased"),
        compare_fixed(glmm_bioero_spo, glmm_spo_simple, "Sponge Bioerosion"),
        compare_fixed(glmm_bioero_pol, glmm_pol_simple, "Polychaete Bioerosion")
)

print(comparison_all)


### Clean Space -----
events_summary
summary(glmm_bl)
summary(glmm_ma)
summary(glmm_mr)
summary(glmm_rtld)
summary(glmm_dss)
summary(glmm_bbd)
summary(glmm_bioero_spo)
summary(glmm_bioero_pol)
diagnostics
coef_all
table_data

rm(list = ls())
graphics.off()
### Summary Results: Predicting coral health by habitat and month---------------
### 3.2 Predicting susceptibility of species bayesian models ------------------

# Paleta de corales
coral_palette <- c(
        "#fff7ec", "#fdd49e", "#fdbb84", "#fc8d59", 
        "#ef6548", "#d7301f", "#990000"
)

library(brms)
packageVersion("brms")
library(tidyverse)
library(ggridges)  # Para ridge plots

### For running models --------
library(purrr)

# Vector de outcomes y nombres
outcomes <- c("BL", "MA", "MR", "RTLD", "DSS", "BBD", "BIOERO_SPO", "BIOERO_POL")
model_names <- paste0("b_", tolower(outcomes), "_species")

# Función para correr modelo
fit_bayes_model <- function(outcome, data) {
        brm(
                formula = as.formula(paste(outcome, "~ Habitat + Month + (1 | Species)")),
                family = bernoulli(link = "logit"),
                data = data,
                chains = 4,
                iter = 2000,
                warmup = 1000,
                cores = 4,
                control = list(adapt_delta = 0.95),
                refresh = 0
        )
}

# Correr todos y guardar
bayes_models <- map(outcomes, ~fit_bayes_model(.x, df_model))
names(bayes_models) <- outcomes

# Guardar todos
walk2(bayes_models, paste0(model_names, ".rds"), saveRDS)

bayes_models <- map(paste0(model_names, ".rds"), readRDS)
names(bayes_models) <- outcomes

b_bl_species <- brm(
        formula = BL ~ Habitat + Month + (1 | Species),
        family = bernoulli(link = "logit"),
        data = df_model,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        control = list(adapt_delta = 0.95))
saveRDS(b_bl_species, "b_bl_species.rds")

b_ma_species <- brm(
        formula = MA ~ Habitat + Month + (1 | Species),
        family = bernoulli(link = "logit"),
        data = df_model,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        control = list(adapt_delta = 0.95)
)
saveRDS(b_ma_species, "b_ma_species.rds")

b_mr_species <- brm(
        formula = MR ~ Habitat + Month + (1 | Species),
        family = bernoulli(link = "logit"),
        data = df_model,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        control = list(adapt_delta = 0.95)
)
saveRDS(b_mr_species, "b_mr_species.rds")

b_rtld_species <- brm(
        formula = RTLD ~ Habitat + Month + (1 | Species),
        family = bernoulli(link = "logit"),
        data = df_model,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        control = list(adapt_delta = 0.95)
)
saveRDS(b_rtld_species, "b_rtld_species.rds")

b_dss_species <- brm(
        formula = DSS ~ Habitat + Month + (1 | Species),
        family = bernoulli(link = "logit"),
        data = df_model,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        control = list(adapt_delta = 0.95)
)
saveRDS(b_dss_species, "b_dss_species.rds")

b_bbd_species <- brm(
        formula = BBD ~ Habitat + Month + (1 | Species),
        family = bernoulli(link = "logit"),
        data = df_model,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        control = list(adapt_delta = 0.95))

saveRDS(b_bbd_species, "b_bbd_species.rds")

b_bio_spo_species <- brm(
        formula = BIOERO_SPO ~ Habitat + Month + (1 | Species),
        family = bernoulli(link = "logit"),
        data = df_model,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        control = list(adapt_delta = 0.95)
)
saveRDS(b_bio_spo_species, "b_bio_spo_species.rds")

b_bio_pol_species <- brm(
        formula = BIOERO_POL ~ Habitat + Month + (1 | Species),
        family = bernoulli(link = "logit"),
        data = df_model,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        control = list(adapt_delta = 0.95)
)
saveRDS(b_bio_pol_species, "b_bio_pol_species.rds")


### Saved bayesian models -----------
getwd()

b_bl_species <- readRDS("D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/b_bl_species.rds")
b_ma_species <- readRDS("D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/b_ma_species.rds")
b_mr_species <- readRDS("D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/b_mr_species.rds")
b_rtld_species <- readRDS("D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/b_rtld_species.rds")
b_bbd_species <- readRDS("D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/b_bbd_species.rds")
b_dss_species <- readRDS("D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/b_dss_species.rds")
b_bio_spo_species <- readRDS("D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/b_bio_spo_species.rds")
b_bio_pol_species <- readRDS("D:/Andy/____Thesis/Paper_Terraza-Veril_JR/01_Analysis/b_bio_pol_species.rds")
summary(b_bl_species)
summary(b_ma_species)
summary(b_mr_species)
summary(b_rtld_species)
summary(b_bbd_species)
summary(b_dss_species)
summary(b_bio_spo_species)
summary(b_bio_pol_species)


# REVISOR

library(brms)
library(tidyverse)

# Weakly-informative priors: normal(0,1) on fixed effects (log-odds scale),
# keeping brms defaults for Intercept and sd (already weakly informative)
weak_priors <- c(
        prior(normal(0, 1), class = "b")
)

# --- BBD sensitivity model ---
b_bbd_species_sens <- brm(
        formula = BBD ~ Habitat + Month + (1 | Species),
        family = bernoulli(link = "logit"),
        data = df_model,
        prior = weak_priors,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        control = list(adapt_delta = 0.95),
        refresh = 0
)
saveRDS(b_bbd_species_sens, "b_bbd_species_sens.rds")

# --- RTLD sensitivity model ---
b_rtld_species_sens <- brm(
        formula = RTLD ~ Habitat + Month + (1 | Species),
        family = bernoulli(link = "logit"),
        data = df_model,
        prior = weak_priors,
        chains = 4,
        iter = 2000,
        warmup = 1000,
        cores = 4,
        control = list(adapt_delta = 0.95),
        refresh = 0
)
saveRDS(b_rtld_species_sens, "b_rtld_species_sens.rds")

# --- Compare default vs. weakly-informative priors ---
compare_priors <- function(default_model, sens_model, outcome_name) {
        d <- fixef(default_model) %>% as.data.frame() %>% rownames_to_column("term") %>% mutate(prior = "default")
        s <- fixef(sens_model) %>% as.data.frame() %>% rownames_to_column("term") %>% mutate(prior = "normal(0,1)")
        bind_rows(d, s) %>% mutate(outcome = outcome_name)
}

comparison_bbd <- compare_priors(b_bbd_species, b_bbd_species_sens, "BBD")
comparison_rtld <- compare_priors(b_rtld_species, b_rtld_species_sens, "RTLD")

comparison <- bind_rows(comparison_bbd, comparison_rtld)
print(comparison)

# Also compare species-level random intercepts (ranef) if that's the key claim
ranef_compare <- function(default_model, sens_model, outcome_name) {
        d <- ranef(default_model)$Species[,,"Intercept"] %>% as.data.frame() %>% rownames_to_column("Species") %>% mutate(prior = "default")
        s <- ranef(sens_model)$Species[,,"Intercept"] %>% as.data.frame() %>% rownames_to_column("Species") %>% mutate(prior = "normal(0,1)")
        bind_rows(d, s) %>% mutate(outcome = outcome_name)
}

ranef_bbd <- ranef_compare(b_bbd_species, b_bbd_species_sens, "BBD")
ranef_rtld <- ranef_compare(b_bbd_species, b_rtld_species_sens, "RTLD")
print(bind_rows(ranef_bbd, ranef_rtld))


### Forest plot for bayesian models full version --------------------------------

create_forest_plot <- function(model, condition_name) {
        # Extraer efectos posteriores
        species_effects <- posterior_samples(model, pars = "^r_Species\\[") %>%
                pivot_longer(cols = everything(), names_to = "Species", values_to = "Intercept") %>%
                mutate(Species = gsub("r_Species\\[|,Intercept\\]", "", Species))
        
        # Resumen
        species_summary <- species_effects %>%
                group_by(Species) %>%
                summarize(
                        median = median(Intercept), 
                        lower = quantile(Intercept, 0.025),
                        upper = quantile(Intercept, 0.975) 
                ) %>%
                ungroup() %>%
                mutate(
                        Significant = ifelse(lower > 0 | upper < 0, "***", "ns"),
                        Species = fct_reorder(Species, median)
                )
        
        # Filtrar

        #species_summary_filtered <- species_summary %>%
        #        filter(Significant == "***") %>%
        #        mutate(Species = fct_reorder(Species, median))
        
        species_summary_filtered <- species_summary %>%
                filter(abs(median) > 0.1) %>%
                mutate(Species = fct_reorder(Species, median))
        
        # Crear gráfico
        p <- ggplot(species_summary_filtered, aes(x = median, y = Species, color = Significant)) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
                geom_point(size = 3) +
                geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, size = 0.8) +
                scale_color_manual(values = c("***" = coral_palette[6], "ns" = "gray70")) +
                scale_y_discrete(labels = function(x) gsub("\\.", " ", x)) +
                labs(
                        title = condition_name,
                        x = "Posterior logit intercept",
                        y = "",
                        color = "Significant"
                ) +
                theme_minimal() +
                theme(
                        axis.text.y = element_text(size = 8, face = "italic"),
                        axis.text.x = element_text(size = 10),
                        axis.title = element_text(face = "bold", size = 10),
                        plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
                        legend.position = "none",
                        legend.text = element_text(size = 8),
                        legend.title = element_text(size = 9),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA)
                )
        
        return(p)
}

# Crear gráficos para cada modelo
p_bl <- create_forest_plot(b_bl_species, "Bleaching")
p_ma <- create_forest_plot(b_ma_species, "Old Mortality")
p_mr <- create_forest_plot(b_mr_species, "Recent Mortality")
p_rtld <- create_forest_plot(b_rtld_species, "Rapid Tissue Loss Disease")
p_bbd <- create_forest_plot(b_bbd_species, "Black Band Disease")
p_dss <- create_forest_plot(b_dss_species, "Dark Spot Syndrome")
p_bio_spo <- create_forest_plot(b_bio_spo_species, "Sponge Bioerosion")
p_bio_pol <- create_forest_plot(b_bio_pol_species, "Polychaete Bioerosion")

library(patchwork)

panel_species_susceptibility <- (p_bl + theme(axis.title.x = element_blank()) | 
                                         p_ma + theme(axis.title.x = element_blank()) | 
                                         p_mr + theme(axis.title.x = element_blank()) | 
                                         p_rtld + theme(axis.title.x = element_blank())) / 
        (p_bbd | p_dss | p_bio_spo | p_bio_pol) +
        plot_layout(guides = "collect")

print(panel_species_susceptibility)

# Guardar
ggsave("Figure_Species_Susceptibility_Panel full.pdf", 
       panel_species_susceptibility,
       width = 16, height = 12, dpi = 300)

### RESULT FIGURE 6. Forest plot for bayesian models positive only ----------------

create_forest_plot_positive_only <- function(model, condition_name) {
        # Extraer efectos posteriores
        species_effects <- posterior_samples(model, pars = "^r_Species\\[") %>%
                pivot_longer(cols = everything(), names_to = "Species", values_to = "Intercept") %>%
                mutate(Species = gsub("r_Species\\[|,Intercept\\]", "", Species))
        
        # Resumen
        species_summary <- species_effects %>%
                group_by(Species) %>%
                summarize(
                        median = median(Intercept), 
                        lower = quantile(Intercept, 0.025),
                        upper = quantile(Intercept, 0.975) 
                ) %>%
                ungroup() %>%
                mutate(
                        Significant = ifelse(lower > 0 | upper < 0, "***", "ns"),
                        Species = fct_reorder(Species, median)
                )
        
        # Filtrar: SOLO especies con mediana > 0
        species_summary_filtered <- species_summary %>%
                filter(median > 0) %>%
                mutate(Species = fct_reorder(Species, median))
        
        # Crear gráfico
        p <- ggplot(species_summary_filtered, aes(x = median, y = Species, color = Significant)) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
                geom_point(size = 1) +
                geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, size = 0.5) +
                scale_color_manual(values = c("***" = coral_palette[6], "ns" = "gray70")) +
                scale_y_discrete(labels = function(x) gsub("\\.", " ", x)) +
                labs(
                        title = condition_name,
                        x = "Posterior logit intercept",
                        y = "",
                        color = "Significant"
                ) +
                theme_minimal() +
                theme(
                        axis.text.y = element_text(size = 8, face = "italic"),
                        axis.text.x = element_text(size = 8),
                        axis.title = element_text(face = "bold", size = 8),
                        plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
                        legend.position = "none",
                        legend.text = element_text(size = 8),
                        legend.title = element_text(size = 8),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA)
                )
        
        return(p)
}

# Crear gráficos para cada modelo (SOLO POSITIVOS)
p_bl_pos <- create_forest_plot_positive_only(b_bl_species, "Bleaching")
p_ma_pos <- create_forest_plot_positive_only(b_ma_species, "Old Mortality")
p_mr_pos <- create_forest_plot_positive_only(b_mr_species, "Recent Mortality")
p_rtld_pos <- create_forest_plot_positive_only(b_rtld_species, "Rapid Tissue Loss Disease")
p_bbd_pos <- create_forest_plot_positive_only(b_bbd_species, "Black Band Disease")
p_dss_pos <- create_forest_plot_positive_only(b_dss_species, "Dark Spot Syndrome")
p_bio_spo_pos <- create_forest_plot_positive_only(b_bio_spo_species, "Sponge Bioerosion")
p_bio_pol_pos <- create_forest_plot_positive_only(b_bio_pol_species, "Polychaete Bioerosion")

# Panel final (SOLO POSITIVOS)
panel_species_susceptibility_positive <- (p_bl_pos + theme(axis.title.x = element_blank()) | 
                                                  p_ma_pos + theme(axis.title.x = element_blank()) | 
                                                  p_mr_pos + theme(axis.title.x = element_blank()) | 
                                                  p_rtld_pos + theme(axis.title.x = element_blank())) / 
        (p_bbd_pos | p_dss_pos | p_bio_spo_pos | p_bio_pol_pos) +
        plot_layout(guides = "collect")

print(panel_species_susceptibility_positive)

ggsave("Figure_6_Species_Susceptibility_Panel.pdf", 
       panel_species_susceptibility_positive,
       width = 11, height = 6, dpi = 300)
ggsave("Figure_6_Species_Susceptibility_Panel.png", 
       panel_species_susceptibility_positive,
       width = 11, height = 6, dpi = 300)

### Fores plot ridge plot version not used ---------
create_forest_plot_with_ridge <- function(model, condition_name, show_x_label = FALSE) {
        # Extraer efectos posteriores
        species_effects <- posterior_samples(model, pars = "^r_Species\\[") %>%
                pivot_longer(cols = everything(), names_to = "Species", values_to = "Intercept") %>%
                mutate(Species = gsub("r_Species\\[|,Intercept\\]", "", Species))
        
        # Resumen
        species_summary <- species_effects %>%
                group_by(Species) %>%
                summarize(
                        median = median(Intercept), 
                        lower = quantile(Intercept, 0.025),
                        upper = quantile(Intercept, 0.975) 
                ) %>%
                ungroup() %>%
                mutate(
                        Significant = ifelse(lower > 0 | upper < 0, "***", "ns"),
                        Species = fct_reorder(Species, median)
                )

        # Filtrar por significancia
        species_summary_filtered <- species_summary %>%
                filter(Significant == "***") %>%
                mutate(Species = fct_reorder(Species, median))
        
        species_summary_filtered <- species_summary %>%
                filter(abs(median) > 0) %>%
                mutate(Species = fct_reorder(Species, median))
        
        # Filtrar species_effects para las especies significativas
        species_effects_filtered <- species_effects %>%
                filter(Species %in% species_summary_filtered$Species) %>%
                # Agregar mediana por especie para colorear
                left_join(species_summary_filtered %>% select(Species, median), by = "Species")
        
        # Ridge plot de fondo con colores según tamaño del efecto
        p <- ggplot(species_effects_filtered, aes(x = Intercept, y = Species, fill = median)) +
                geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, alpha = 0.7, size = 0) +
                scale_fill_gradientn(colors = coral_palette, guide = "none") +
                # Agregar forest plot encima
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
                geom_point(data = species_summary_filtered, aes(x = median, y = Species, color = Significant), 
                           size = 3, inherit.aes = FALSE) +
                geom_errorbarh(data = species_summary_filtered, aes(xmin = lower, xmax = upper, y = Species, 
                                                                    color = Significant), 
                               height = 0.15, size = 0.8, inherit.aes = FALSE) +
                scale_color_manual(values = c("***" = "#000000", "ns" = "gray70")) +
                scale_y_discrete(labels = function(x) gsub("\\.", " ", x)) +
                labs(
                        title = condition_name,
                        x = if(show_x_label) "Posterior logit intercept" else "",
                        y = "",
                        color = "Significant"
                ) +
                theme_minimal() +
                theme(
                        axis.text.y = element_text(size = 8, face = "italic"),
                        axis.text.x = element_text(size = 9),
                        axis.title = element_text(face = "bold", size = 10),
                        plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
                        legend.position = "none",
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA)
                )
        
        return(p)
}
# Crear gráficos
p_bl <- create_forest_plot_with_ridge(b_bl_species, "Bleaching", show_x_label = FALSE)
p_bl
p_ma <- create_forest_plot_with_ridge(b_ma_species, "Old Mortality", show_x_label = FALSE)
p_mr <- create_forest_plot_with_ridge(b_mr_species, "Recent Mortality", show_x_label = FALSE)
p_rtld <- create_forest_plot_with_ridge(b_rtld_species, "Rapid Tissue Loss Disease", show_x_label = FALSE)
p_bbd <- create_forest_plot_with_ridge(b_bbd_species, "Black Band Disease", show_x_label = TRUE)
p_dss <- create_forest_plot_with_ridge(b_dss_species, "Dark Spot Syndrome", show_x_label = TRUE)
p_bio_spo <- create_forest_plot_with_ridge(b_bio_spo_species, "Sponge Bioerosion", show_x_label = TRUE)
p_bio_pol <- create_forest_plot_with_ridge(b_bio_pol_species, "Polychaete Bioerosion", show_x_label = TRUE)
p_ma
# Panel final
panel_species_susceptibility <- (p_bl | p_ma | p_mr | p_rtld) / 
        (p_bbd | p_dss | p_bio_spo | p_bio_pol)

print(panel_species_susceptibility)

ggsave("Figure_Species_Susceptibility_Panel_Combined.pdf", 
       panel_species_susceptibility,
       width = 16, height = 12, dpi = 300)

### RESULT TABLE 9. Bayesian model susceptibility per Species --------

library(flextable)
library(officer)

# Función para extraer resumen de especies
extract_species_summary <- function(model, condition_name) {
        species_effects <- posterior_samples(model, pars = "^r_Species\\[") %>%
                pivot_longer(cols = everything(), names_to = "Species", values_to = "Intercept") %>%
                mutate(Species = gsub("r_Species\\[|,Intercept\\]", "", Species))
        
        species_summary <- species_effects %>%
                group_by(Species) %>%
                summarize(
                        median = median(Intercept), 
                        lower = quantile(Intercept, 0.025),
                        upper = quantile(Intercept, 0.975) 
                ) %>%
                ungroup() %>%
                mutate(
                        Significant = ifelse(lower > 0 | upper < 0, "***", "ns"),
                        Condition = condition_name
                ) %>%
                arrange(desc(median))
        
        return(species_summary)
}

# Extraer para todas las condiciones
species_table <- bind_rows(
        extract_species_summary(b_bl_species, "Bleaching"),
        extract_species_summary(b_ma_species, "Old Mortality"),
        extract_species_summary(b_mr_species, "Recent Mortality"),
        extract_species_summary(b_rtld_species, "Rapid Tissue Loss Disease"),
        extract_species_summary(b_bbd_species, "Black Band Disease"),
        extract_species_summary(b_dss_species, "Dark Spot Syndrome"),
        extract_species_summary(b_bio_spo_species, "Sponge Bioerosion"),
        extract_species_summary(b_bio_pol_species, "Polychaete Bioerosion")
)

# Preparar tabla: SOLO median > 0.3
table_supp <- species_table %>%
        filter(Significant == TRUE | abs(median) > 0.3) %>%
        mutate(
                median = round(median, 2),
                lower  = round(lower, 2),
                upper  = round(upper, 2),
                Species = gsub("\\.", " ", Species),
                Median_CI = paste0(
                        sprintf("%.3f", median),
                        " [",
                        sprintf("%.3f", lower),
                        ", ",
                        sprintf("%.3f", upper),
                        "]"
                )
        ) %>%
        select(Condition, Species, median, Median_CI, Significant) %>%
        rename(
                "Health Condition" = Condition,
                "Species" = Species,
                "Posterior Median" = median,
                "95% Credible Interval" = Median_CI,
                "Significant" = Significant
        ) %>%
        arrange(`Health Condition`, desc(`Posterior Median`))

# Crear flextable
ft_supp <- flextable(table_supp) %>%
        theme_booktabs() %>%
        
        # Combinar celdas de "Health Condition"
        merge_v(j = 1, part = "body") %>%
        
        # Anchos
        width(j = 1, width = 1.5) %>%
        width(j = 2, width = 1.5) %>%
        width(j = 3, width = 1.0) %>%
        width(j = 4, width = 1.5) %>%
        width(j = 5, width = 0.8) %>%
        
        # Alineación
        align(j = 1, align = "left", part = "all") %>%
        align(j = 2, align = "left", part = "all") %>%
        align(j = 3:5, align = "center", part = "all") %>%
        valign(j = 1, valign = "center", part = "body") %>%
        
        # Tipografía
        fontsize(size = 10, part = "body") %>%
        fontsize(size = 10, part = "header") %>%
        
        # Especies en cursiva
        italic(j = 2, part = "body") %>%
        
        # Health Condition en negrita y gris
        bold(j = 1, part = "body") %>%

        # Padding
        padding(j = 1:5, padding.left = 3, padding.right = 3, part = "all") %>%
        padding(j = 1:5, padding.top = 2, padding.bottom = 2, part = "all") %>%
        
        # Líneas finas horizontales entre grupos (Health Condition)
        {
                # Encontrar dónde cambia la condición
                cambios <- which(table_supp$`Health Condition` != dplyr::lag(table_supp$`Health Condition`))
                #cambios <- cambios[-1]  # Remover el primer cambio
                filas_con_linea <- cambios - 1
                
                hline(., i = filas_con_linea,
                      border = fp_border(color = "#CCCCCC", width = 0.5),
                      part = "body")
        } %>%
        
        # Encabezado descriptivo
        add_header_lines(
                "Table 9. Coral species showing significant or strong deviations in susceptibility to health conditions. Bayesian hierarchical models identified species-level variation in disease and condition prevalence. The table includes species with statistically significant effects and/or substantial effect sizes (|posterior median logit intercept| > 0.3), indicating markedly higher or lower odds of infection/condition occurrence relative to the population average."
        ) %>%
        
        # Nota de pie
        add_footer_lines(
                "Posterior median logit intercepts and 95% credible intervals from species random intercepts. Logit scale: positive values indicate increased susceptibility relative to overall mean. *** indicates credible interval not crossing zero (statistically significant species effect)."
        ) %>%
        
        fontsize(part = "footer", size = 8) %>%
        align(part = "footer", align = "left") %>%
        italic(part = "footer")

print(ft_supp)

# Guardar
save_as_docx(ft_supp, path = "Table_9_Predicting Health_Species_Susceptibility_Bayesian.docx")
print(table_supp, n = Inf)





### Clean Space -----

summary(b_bl_species)
summary(b_ma_species)
summary(b_mr_species)
summary(b_rtld_species)
summary(b_bbd_species)
summary(b_dss_species)
summary(b_bio_spo_species)
summary(b_bio_pol_species)
print(table_supp, n = Inf)

rm(list = ls())
graphics.off()
### Summary Results: Predicting coral health by species