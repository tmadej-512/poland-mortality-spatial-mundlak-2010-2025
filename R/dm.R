# ============================================================
# Dataset powiat × rok ze standaryzowanym wiekiem i płcią
# współczynnikiem zgonów
#
# Standard: średnia struktura ludności Polski 2010–2019
# Grupy wieku: 5-letnie + 85+
# Bez ESP
# ============================================================

rm(list = ls())

detachAllPackages <- function() {
  basic.packages <- c(
    "package:stats",
    "package:graphics",
    "package:grDevices",
    "package:utils",
    "package:datasets",
    "package:methods",
    "package:base"
  )
  
  package.list <- search()[
    ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)
  ]
  
  package.list <- setdiff(package.list, basic.packages)
  
  if (length(package.list) > 0) {
    for (package in package.list) {
      detach(package, character.only = TRUE)
    }
  }
}

detachAllPackages()

options(dplyr.legacy_locale = TRUE)

library(readxl)
library(tidyr)
library(dplyr)
library(mgcv)
library(ggplot2)
library(ISOweek)
library(tsibble)
library(fable)
library(ggpubr)
library(sf)
library(cols4all)
library(Cairo)
library(stringr)
library(gratia)
library(patchwork)
library(tibble)

# ============================================================
# 0. Ścieżki
# ============================================================

setwd("C:/Users/panzc/Documents/")

base_dir <- "C:/Users/panzc/Documents/Instytut/covid"

path_deaths_m <- "~/Instytut/covid/zgony_gus/zgony_powiaty/zgony_1999_2023_m.csv"
path_deaths_f <- "~/Instytut/covid/zgony_gus/zgony_powiaty/zgony_1999_2023_f.csv"

path_population <- "C:/Users/panzc/Documents/Instytut/covid/analiza_przestrzenna/population_vars/populacja_powiaty_2002_2023.csv"

path_powiaty_shp <- "C:/Users/panzc/Documents/Instytut/covid/miasta/Powiaty.shp"

path_unmp <- "~/Instytut/covid/analiza_przestrzenna/exploratory_vars/bezrobocie_w_aktywnych.csv"
path_inc <- "~/Instytut/covid/analiza_przestrzenna/exploratory_vars/przecietne_wyn_msc_brutto.csv"
path_docs <- "~/Instytut/covid/analiza_przestrzenna/exploratory_vars/lekarze_pracujacy.csv"
path_beds <- "~/Instytut/covid/analiza_przestrzenna/exploratory_vars/liczba_lozek.csv"
path_med_exp <- "~/Instytut/covid/analiza_przestrzenna/exploratory_vars/wydatki_na_ochrone_zdrowia_851.csv"
path_edu <- "~/Instytut/covid/analiza_przestrzenna/exploratory_vars/wyksztalcenie.csv"
path_inc_bud <- "~/Instytut/covid/analiza_przestrzenna/exploratory_vars/dochody_wlasne_budzetow_powiatow.csv"
path_mig <- "~/Instytut/covid/analiza_przestrzenna/exploratory_vars/migracje.csv"
path_temp_poll <- "~/Instytut/covid/analiza_przestrzenna/exploratory_vars/pollution_temp.csv"

path_out <- "C:/Users/panzc/Documents/Instytut/covid/analiza_przestrzenna/dataset/zbior_asmr_5yr_age_sex_2010_2019std.rds"

# ============================================================
# 1. Funkcje pomocnicze
# ============================================================

age_levels <- c(
  "0-4", "5-9", "10-14", "15-19",
  "20-24", "25-29", "30-34", "35-39",
  "40-44", "45-49", "50-54", "55-59",
  "60-64", "65-69", "70-74", "75-79",
  "80-84", "85+"
)

fix_powiat_names_codes <- function(df, code_col, name_col) {
  code_col <- rlang::ensym(code_col)
  name_col <- rlang::ensym(name_col)
  
  df %>%
    mutate(
      !!name_col := if_else(
        as.character(!!name_col) == "Powiat warszawski",
        "Powiat m. st. Warszawa",
        as.character(!!name_col)
      ),
      !!name_col := if_else(
        as.character(!!name_col) %in% c(
          "Powiat m. Wałbrzych do 2002",
          "Powiat m. Wałbrzych od 2013"
        ),
        "Powiat wałbrzyski",
        as.character(!!name_col)
      ),
      !!code_col := as.numeric(as.character(!!code_col)),
      !!code_col := if_else(!!code_col == 1431000, 1465000, !!code_col),
      !!code_col := if_else(!!code_col %in% c(263000, 265000), 221000, !!code_col)
    )
}

make_code_pow <- function(code) {
  code_chr <- as.character(code)
  
  code_chr <- if_else(
    nchar(code_chr) == 6,
    paste0("0", code_chr),
    code_chr
  )
  
  substring(code_chr, 1, 4)
}

recode_age_population <- function(age_group) {
  age_group <- as.character(age_group)
  
  case_when(
    age_group %in% c("0-4") ~ "0-4",
    age_group %in% c("5-9") ~ "5-9",
    age_group %in% c("10-14") ~ "10-14",
    age_group %in% c("15-19") ~ "15-19",
    age_group %in% c("20-24") ~ "20-24",
    age_group %in% c("25-29") ~ "25-29",
    age_group %in% c("30-34") ~ "30-34",
    age_group %in% c("35-39") ~ "35-39",
    age_group %in% c("40-44") ~ "40-44",
    age_group %in% c("45-49") ~ "45-49",
    age_group %in% c("50-54") ~ "50-54",
    age_group %in% c("55-59") ~ "55-59",
    age_group %in% c("60-64") ~ "60-64",
    age_group %in% c("65-69") ~ "65-69",
    age_group %in% c("70-74") ~ "70-74",
    age_group %in% c("75-79") ~ "75-79",
    age_group %in% c("80-84") ~ "80-84",
    age_group %in% c("85 i więcej", "85 i wiecej", "85+", "85 lat i więcej", "85 lat i wiecej") ~ "85+",
    TRUE ~ NA_character_
  )
}

recode_age_deaths <- function(age_group) {
  age_group <- as.character(age_group)
  
  case_when(
    age_group %in% c("0", "1-4") ~ "0-4",
    age_group %in% c("5-9") ~ "5-9",
    age_group %in% c("10-14") ~ "10-14",
    age_group %in% c("15-19") ~ "15-19",
    age_group %in% c("20-24") ~ "20-24",
    age_group %in% c("25-29") ~ "25-29",
    age_group %in% c("30-34") ~ "30-34",
    age_group %in% c("35-39") ~ "35-39",
    age_group %in% c("40-44") ~ "40-44",
    age_group %in% c("45-49") ~ "45-49",
    age_group %in% c("50-54") ~ "50-54",
    age_group %in% c("55-59") ~ "55-59",
    age_group %in% c("60-64") ~ "60-64",
    age_group %in% c("65-69") ~ "65-69",
    age_group %in% c("70-74") ~ "70-74",
    age_group %in% c("75-79") ~ "75-79",
    age_group %in% c("80-84") ~ "80-84",
    age_group %in% c("85 i więcej", "85 i wiecej", "85+", "85 lat i więcej", "85 lat i wiecej") ~ "85+",
    TRUE ~ NA_character_
  )
}

safe_num <- function(x) {
  as.numeric(gsub(",", ".", as.character(x)))
}

# ============================================================
# 2. Wczytanie zgonów i geometrii powiatów
# ============================================================

zgony_1999_2023_m <- read.csv(path_deaths_m, sep = ";")
zgony_1999_2023_f <- read.csv(path_deaths_f, sep = ";")

zgony_1999_2023 <- bind_rows(
  zgony_1999_2023_m,
  zgony_1999_2023_f
)

powiaty_spdf <- st_read(path_powiaty_shp) %>%
  mutate(
    JPT_NAZWA_ = if_else(
      JPT_NAZWA_ %in% c("powiat Wałbrzych", "powiat wałbrzyski"),
      "powiat wałbrzyski",
      JPT_NAZWA_
    ),
    JPT_KOD_JE = if_else(
      JPT_KOD_JE %in% c("0265"),
      "0221",
      JPT_KOD_JE
    )
  ) %>%
  group_by(JPT_KOD_JE, JPT_NAZWA_) %>%
  summarise(
    geometry = st_union(geometry),
    .groups = "drop"
  )

# ============================================================
# 3. Wczytanie i przygotowanie ludności
# ============================================================

ludnosc_powiaty <- read.csv(path_population)[, -1]

ludnosc <- ludnosc_powiaty %>%
  select(
    code = id_teryt,
    year,
    name,
    age_group,
    pop,
    sex
  ) %>%
  fix_powiat_names_codes(code_col = code, name_col = name) %>%
  mutate(
    year = as.integer(year),
    pop = as.numeric(pop),
    sex = as.character(sex),
    age_group_raw = as.character(age_group),
    age_group = recode_age_population(age_group_raw)
  ) %>%
  filter(
    year >= 2002,
    !is.na(age_group),
    !is.na(sex),
    !is.na(pop),
    pop > 0
  ) %>%
  group_by(code, name, year, sex, age_group) %>%
  summarise(
    pop = sum(pop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    code = as.character(code),
    name = as.character(name),
    sex = as.character(sex),
    age_group = factor(age_group, levels = age_levels)
  )

# ============================================================
# 4. Przygotowanie zgonów
# ============================================================
# Wyrzucamy agregaty, żeby nie dublować zgonów:
# - "poniżej 5" = 0 + 1-4
# - "65 i więcej", "70 i więcej" itd. to agregaty wieku
# Zostają szczegółowe grupy: 0, 1-4, 5-9, ..., 80-84, 85+
# Następnie 0 + 1-4 agregujemy do 0-4.

DEATHS <- zgony_1999_2023 %>%
  filter(
    str_detect(Nazwa, "Powiat"),
    Wiek != "ogółem",
    Rok >= 2002
  ) %>%
  filter(
    !Wiek %in% c(
      "poniżej 5",
      "65 i więcej",
      "70 i więcej",
      "75 i więcej",
      "80 i więcej"
    )
  ) %>%
  fix_powiat_names_codes(code_col = Kod, name_col = Nazwa) %>%
  mutate(
    year = as.integer(Rok),
    code = as.character(Kod),
    name = as.character(Nazwa),
    sex = as.character(Płeć),
    age_group_raw = as.character(Wiek),
    age_group = recode_age_deaths(age_group_raw),
    deaths = as.numeric(Wartosc)
  ) %>%
  filter(
    !is.na(age_group),
    !is.na(sex),
    !is.na(deaths)
  ) %>%
  group_by(code, name, sex, year, age_group) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    age_group = factor(age_group, levels = age_levels),
    code_pow = make_code_pow(code)
  ) %>%
  left_join(
    powiaty_spdf,
    by = c("code_pow" = "JPT_KOD_JE")
  ) %>%
  select(
    code,
    code_pow,
    name,
    sex,
    year,
    age_group,
    deaths,
    geometry
  )

# ============================================================
# 5. Kontrole grup wieku i płci
# ============================================================

cat("\nGrupy wieku w populacji:\n")
print(
  ludnosc %>%
    count(age_group) %>%
    arrange(age_group),
  n = 100
)

cat("\nGrupy wieku w zgonach:\n")
print(
  DEATHS %>%
    st_drop_geometry() %>%
    count(age_group) %>%
    arrange(age_group),
  n = 100
)

cat("\nPłeć w populacji:\n")
print(
  ludnosc %>%
    count(sex)
)

cat("\nPłeć w zgonach:\n")
print(
  DEATHS %>%
    st_drop_geometry() %>%
    count(sex)
)

# ============================================================
# 6. Dołączenie populacji do komórek zgony × wiek × płeć
# ============================================================

mort_cells <- DEATHS %>%
  st_drop_geometry() %>%
  as_tibble() %>%
  select(code, name, year, sex, age_group, deaths) %>%
  inner_join(
    ludnosc %>%
      select(code, name, year, sex, age_group, pop),
    by = c("code", "name", "year", "sex", "age_group")
  ) %>%
  mutate(
    mx = deaths / pop
  )

missing_pop_check <- DEATHS %>%
  st_drop_geometry() %>%
  as_tibble() %>%
  select(code, name, year, sex, age_group, deaths) %>%
  anti_join(
    ludnosc %>%
      select(code, name, year, sex, age_group, pop),
    by = c("code", "name", "year", "sex", "age_group")
  )

if (nrow(missing_pop_check) > 0) {
  warning("Są komórki zgonów bez odpowiadającej populacji.")
  print(missing_pop_check, n = 200)
}

cat("\nLiczba komórek zgonów po połączeniu z populacją:\n")
print(nrow(mort_cells))

# ============================================================
# 7. Populacja standardowa:
#    średnia struktura Polski z lat 2010–2019
#    po płci i 5-letnich grupach wieku
# ============================================================

std_pop_2010_2019 <- ludnosc %>%
  filter(year >= 2010, year <= 2019) %>%
  group_by(year, sex, age_group) %>%
  summarise(
    pop_cell_year = sum(pop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sex, age_group) %>%
  summarise(
    std_pop = mean(pop_cell_year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    std_weight = std_pop / sum(std_pop, na.rm = TRUE)
  ) %>%
  arrange(sex, age_group) %>%
  select(sex, age_group, std_pop, std_weight)

cat("\nStandard 2010–2019, sex × age_group:\n")
print(std_pop_2010_2019, n = 100)

std_weight_sum <- sum(std_pop_2010_2019$std_weight, na.rm = TRUE)

cat("\nSuma wag standardowych:\n")
print(std_weight_sum)

if (abs(std_weight_sum - 1) > 1e-8) {
  warning("Suma wag standardowych nie wynosi 1.")
}

cat("\nLiczba komórek standardu:\n")
print(nrow(std_pop_2010_2019))

# Przy 2 płciach i 18 grupach wieku powinno być 36 komórek.
if (nrow(std_pop_2010_2019) != 36) {
  warning("Standard nie ma 36 komórek. Sprawdź grupy wieku/płeć.")
}

# ============================================================
# 8. ASMR: standaryzowany wiekiem i płcią współczynnik zgonów
# ============================================================

asmr_powiat <- mort_cells %>%
  inner_join(
    std_pop_2010_2019 %>%
      select(sex, age_group, std_weight),
    by = c("sex", "age_group")
  ) %>%
  mutate(
    weighted_mx = mx * std_weight
  ) %>%
  group_by(code, name, year) %>%
  summarise(
    asmr_2010_2019std = sum(weighted_mx, na.rm = TRUE) * 100000,
    deaths_total = sum(deaths, na.rm = TRUE),
    pop_total = sum(pop, na.rm = TRUE),
    crude_mort_rate = deaths_total / pop_total * 100000,
    .groups = "drop"
  ) %>%
  arrange(code, year)

coverage_check <- mort_cells %>%
  distinct(code, name, year, sex, age_group) %>%
  anti_join(
    std_pop_2010_2019 %>%
      distinct(sex, age_group),
    by = c("sex", "age_group")
  )

if (nrow(coverage_check) > 0) {
  warning("Są komórki wieku/płci bez wagi standardowej.")
  print(coverage_check, n = 200)
} else {
  message("OK: wszystkie komórki wieku/płci mają wagę standardową.")
}

cat("\nPodgląd ASMR:\n")
print(asmr_powiat, n = 100)

# ============================================================
# 9. Centroidy powiatów i dystans do najbliższego miasta
# ============================================================

toxy <- DEATHS %>%
  group_by(geometry) %>%
  distinct(code) %>%
  ungroup()

toxy$x <- do.call(rbind, st_centroid(toxy$geometry))[, 1]
toxy$y <- do.call(rbind, st_centroid(toxy$geometry))[, 2]

toxy <- toxy %>%
  st_drop_geometry() %>%
  select(code, x, y) %>%
  distinct()

DEATHS <- DEATHS %>%
  left_join(toxy, by = "code")

df_map <- DEATHS %>%
  ungroup() %>%
  mutate(
    city = if_else(grepl("Powiat m\\.", name), "city", "rural")
  ) %>%
  st_as_sf(coords = c("x", "y"), crs = 2180)

cities <- df_map %>%
  filter(city == "city")

others <- df_map %>%
  filter(city != "city")

nearest_idx <- st_nearest_feature(others, cities)

dist_km <- as.numeric(
  st_distance(
    others,
    cities[nearest_idx, ],
    by_element = TRUE
  )
) / 1000

others2 <- others %>%
  mutate(
    nearest_city = cities$name[nearest_idx],
    distance_km = dist_km
  )

cities2 <- cities %>%
  mutate(
    nearest_city = name,
    distance_km = 0
  )

df_distance_city <- bind_rows(others2, cities2) %>%
  bind_cols(
    st_coordinates(.) |>
      as.data.frame() |>
      dplyr::rename(x = X, y = Y)
  ) %>%
  st_drop_geometry() %>%
  select(
    code,
    distance_to_city = distance_km
  ) %>%
  distinct()

# ============================================================
# 10. Zmienne pomocnicze z populacji
# ============================================================

POP_MF <- ludnosc %>%
  group_by(code, name, year, sex) %>%
  summarise(
    pop = sum(pop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sex,
    values_from = pop
  ) %>%
  mutate(
    pop_mf = `mężczyźni` / kobiety * 100
  ) %>%
  select(code, name, year, pop_mf)

ludnosc_total <- ludnosc %>%
  group_by(code, name, year) %>%
  summarise(
    pop_tot = sum(pop, na.rm = TRUE),
    .groups = "drop"
  )

POP_65 <- ludnosc %>%
  filter(age_group %in% c("65-69", "70-74", "75-79", "80-84", "85+")) %>%
  group_by(code, name, year) %>%
  summarise(
    pop_65_n = sum(pop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    ludnosc_total,
    by = c("code", "name", "year")
  ) %>%
  mutate(
    pop_65 = pop_65_n / pop_tot * 100
  ) %>%
  select(code, name, year, pop_65, pop_tot)

ludnosc_mianownik <- ludnosc_total %>%
  select(code, year, pop = pop_tot)

# ============================================================
# 11. Wczytanie zmiennych objaśniających
# ============================================================

UNMP <- read.csv(path_unmp, sep = ";") %>%
  mutate(
    Wartosc = safe_num(Wartosc)
  ) %>%
  select(
    code = Kod,
    year = Rok,
    unmp = Wartosc
  ) %>%
  mutate(
    code = as.character(code),
    year = as.integer(year)
  )

INC <- read.csv(path_inc, sep = ";") %>%
  mutate(
    Wartosc = safe_num(Wartosc)
  ) %>%
  select(
    code = Kod,
    year = Rok,
    inc = Wartosc
  ) %>%
  mutate(
    code = as.character(code),
    year = as.integer(year)
  )

DOCS <- read.csv(path_docs, sep = ";") %>%
  select(
    code = Kod,
    year = Rok,
    docs = Wartosc
  ) %>%
  mutate(
    code = as.character(code),
    year = as.integer(year),
    docs = as.numeric(docs)
  ) %>%
  arrange(code, year) %>%
  group_by(code) %>%
  fill(docs, .direction = "down") %>%
  ungroup() %>%
  left_join(
    ludnosc_mianownik,
    by = c("code", "year")
  ) %>%
  mutate(
    docs = round(docs / pop * 1e5, 2)
  ) %>%
  select(code, year, docs)

BEDS <- read.csv(path_beds, sep = ";") %>%
  select(
    code = Kod,
    year = Rok,
    beds = Wartosc
  ) %>%
  mutate(
    code = as.character(code),
    year = as.integer(year),
    beds = as.numeric(beds)
  ) %>%
  arrange(code, year) %>%
  group_by(code) %>%
  fill(beds, .direction = "down") %>%
  ungroup() %>%
  left_join(
    ludnosc_mianownik,
    by = c("code", "year")
  ) %>%
  mutate(
    beds = round(beds / pop * 1e5, 2)
  ) %>%
  select(code, year, beds)

MED_EXP <- read.csv(path_med_exp, sep = ";") %>%
  mutate(
    Wartosc = safe_num(Wartosc)
  ) %>%
  filter(
    Jednostki.terytorialne == "gminy łącznie z miastami na prawach powiatu"
  ) %>%
  select(
    code = Kod,
    year = Rok,
    med_exp = Wartosc
  ) %>%
  mutate(
    code = as.character(code),
    year = as.integer(year)
  ) %>%
  arrange(code, year) %>%
  group_by(code) %>%
  fill(med_exp, .direction = "down") %>%
  ungroup() %>%
  left_join(
    ludnosc_mianownik,
    by = c("code", "year")
  ) %>%
  mutate(
    med_exp = round(med_exp / pop * 1e5, 2)
  ) %>%
  select(code, year, med_exp)

SEC <- read.csv(path_edu, sep = ";") %>%
  filter(
    Poziomy.wykształcenia %in% c(
      "ogółem",
      "gimnazjalne i podstawowe ukończone"
    )
  ) %>%
  select(
    code = Kod,
    edu = Poziomy.wykształcenia,
    value = Wartosc
  ) %>%
  mutate(
    code = as.character(code),
    value = as.numeric(value)
  ) %>%
  pivot_wider(
    names_from = edu,
    values_from = value
  ) %>%
  mutate(
    sec = `gimnazjalne i podstawowe ukończone` / ogółem * 100
  ) %>%
  select(code, sec)

UNI <- read.csv(path_edu, sep = ";") %>%
  filter(
    Poziomy.wykształcenia %in% c("ogółem", "wyższe")
  ) %>%
  select(
    code = Kod,
    edu = Poziomy.wykształcenia,
    value = Wartosc
  ) %>%
  mutate(
    code = as.character(code),
    value = as.numeric(value)
  ) %>%
  pivot_wider(
    names_from = edu,
    values_from = value
  ) %>%
  mutate(
    uni = wyższe / ogółem * 100
  ) %>%
  select(code, uni)

INC_BUD <- read.csv(path_inc_bud, sep = ";") %>%
  mutate(
    Wartosc = safe_num(Wartosc)
  ) %>%
  filter(
    Jednostki.terytorialne == "gminy łącznie z miastami na prawach powiatu"
  ) %>%
  select(
    code = Kod,
    year = Rok,
    inc_bud = Wartosc
  ) %>%
  mutate(
    code = as.character(code),
    year = as.integer(year)
  ) %>%
  left_join(
    ludnosc_mianownik,
    by = c("code", "year")
  ) %>%
  mutate(
    inc_bud = round(inc_bud / pop * 1e5, 2)
  ) %>%
  filter(complete.cases(.)) %>%
  select(code, year, inc_bud)

MIG <- read.csv(path_mig, sep = ";") %>%
  select(
    code = Kod,
    year = Rok,
    mig = Wartosc
  ) %>%
  mutate(
    code = as.character(code),
    year = as.integer(year),
    mig = as.numeric(mig)
  )

TEMP_POLL <- read.csv(path_temp_poll) %>%
  pivot_wider(
    names_from = indicator,
    values_from = val
  ) %>%
  rename(
    code_pow = JPT_KOD_JE
  ) %>%
  mutate(
    code_pow = as.character(code_pow),
    code_pow = if_else(
      nchar(code_pow) == 3,
      paste0("0", code_pow),
      code_pow
    ),
    year = as.integer(year)
  ) %>%
  select(-JPT_NAZWA_)

# ============================================================
# 12. Miasta, powierzchnia, gęstość
# ============================================================

CITIES <- powiaty_spdf %>%
  mutate(
    city = str_detect(JPT_NAZWA_, "powiat [A-ZĄĆĘŁŃÓŚŹŻ]")
  ) %>%
  select(
    code_pow = JPT_KOD_JE,
    city
  ) %>%
  st_drop_geometry()

AREA <- data.frame(
  code_pow = powiaty_spdf$JPT_KOD_JE,
  area = as.numeric(st_area(powiaty_spdf) / 1e6)
)

# ============================================================
# 13. Finalny dataset powiat × rok
# ============================================================

df <- asmr_powiat %>%
  mutate(
    code = as.character(code),
    year = as.integer(year),
    code_pow = make_code_pow(code)
  ) %>%
  left_join(
    POP_65 %>%
      select(code, year, pop_65, pop_tot),
    by = c("code", "year")
  ) %>%
  left_join(
    POP_MF %>%
      select(code, year, pop_mf),
    by = c("code", "year")
  ) %>%
  left_join(UNMP, by = c("code", "year")) %>%
  left_join(INC, by = c("code", "year")) %>%
  left_join(DOCS, by = c("code", "year")) %>%
  left_join(BEDS, by = c("code", "year")) %>%
  left_join(MED_EXP, by = c("code", "year")) %>%
  left_join(UNI, by = "code") %>%
  left_join(SEC, by = "code") %>%
  left_join(INC_BUD, by = c("code", "year")) %>%
  left_join(MIG, by = c("code", "year")) %>%
  left_join(AREA, by = "code_pow") %>%
  left_join(CITIES, by = "code_pow") %>%
  left_join(TEMP_POLL, by = c("code_pow", "year")) %>%
  left_join(df_distance_city, by = "code") %>%
  mutate(
    dens = pop_tot / area,
    beds = if_else(is.na(beds), 0, beds),
    docs = if_else(is.na(docs), 0, docs),
    Intercept = 1
  ) %>%
  filter(year >= 2010) %>%
  select(
    code,
    code_pow,
    name,
    year,
    asmr_2010_2019std,
    crude_mort_rate,
    deaths_total,
    pop_total,
    pop_tot,
    pop_65,
    pop_mf,
    dens,
    area,
    city,
    distance_to_city,
    unmp,
    inc,
    docs,
    beds,
    med_exp,
    uni,
    sec,
    inc_bud,
    mig,
    everything()
  )

# ============================================================
# 14. Kontrole końcowe
# ============================================================

cat("\nLiczba wierszy finalnego datasetu:\n")
print(nrow(df))

cat("\nZakres lat:\n")
print(range(df$year, na.rm = TRUE))

cat("\nLiczba powiatów:\n")
print(n_distinct(df$code))

cat("\nBraki w kluczowych zmiennych:\n")
print(
  df %>%
    summarise(
      n = n(),
      missing_asmr = sum(is.na(asmr_2010_2019std)),
      missing_crude = sum(is.na(crude_mort_rate)),
      missing_pop_tot = sum(is.na(pop_tot)),
      missing_area = sum(is.na(area)),
      missing_distance = sum(is.na(distance_to_city))
    )
)

cat("\nPodgląd finalnego datasetu:\n")
print(df, n = 100)

# ============================================================
# 15. Zapis
# ============================================================

saveRDS(df, path_out)

message("Zapisano plik: ", path_out)


# ============================================================
# 16. Sanity check standaryzacji
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ------------------------------------------------------------
# 16.1. Podstawowe porównanie crude vs ASMR
# ------------------------------------------------------------

sanity_summary <- df %>%
  summarise(
    n = n(),
    crude_min = min(crude_mort_rate, na.rm = TRUE),
    crude_p25 = quantile(crude_mort_rate, 0.25, na.rm = TRUE),
    crude_median = median(crude_mort_rate, na.rm = TRUE),
    crude_p75 = quantile(crude_mort_rate, 0.75, na.rm = TRUE),
    crude_max = max(crude_mort_rate, na.rm = TRUE),
    asmr_min = min(asmr_2010_2019std, na.rm = TRUE),
    asmr_p25 = quantile(asmr_2010_2019std, 0.25, na.rm = TRUE),
    asmr_median = median(asmr_2010_2019std, na.rm = TRUE),
    asmr_p75 = quantile(asmr_2010_2019std, 0.75, na.rm = TRUE),
    asmr_max = max(asmr_2010_2019std, na.rm = TRUE),
    corr_crude_asmr = cor(crude_mort_rate, asmr_2010_2019std, use = "complete.obs"),
    corr_crude_pop65 = cor(crude_mort_rate, pop_65, use = "complete.obs"),
    corr_asmr_pop65 = cor(asmr_2010_2019std, pop_65, use = "complete.obs")
  )

cat("\nSanity summary crude vs ASMR:\n")
print(sanity_summary)

# Interpretacja szybka:
# - corr_crude_asmr powinno być dodatnie i raczej wysokie.
# - corr_asmr_pop65 powinno być niższe niż corr_crude_pop65.
#   To znaczy, że standaryzacja ograniczyła wpływ struktury wieku.

# ------------------------------------------------------------
# 16.2. Trend krajowy: crude vs ASMR
# ------------------------------------------------------------
# Liczymy średnią ważoną powiatowych współczynników populacją.
# Dla ASMR ważenie populacją powiatów służy tylko do trendu krajowego.

trend_country <- df %>%
  group_by(year) %>%
  summarise(
    crude_country = sum(deaths_total, na.rm = TRUE) / sum(pop_total, na.rm = TRUE) * 100000,
    asmr_country = weighted.mean(asmr_2010_2019std, w = pop_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(crude_country, asmr_country),
    names_to = "rate_type",
    values_to = "rate"
  ) %>%
  mutate(
    rate_type = recode(
      rate_type,
      crude_country = "Surowy współczynnik",
      asmr_country = "Standaryzowany wiekiem i płcią"
    )
  )

p_trend <- ggplot(
  trend_country,
  aes(x = year, y = rate, color = rate_type)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Trend krajowy: surowy vs standaryzowany współczynnik zgonów",
    subtitle = "ASMR liczony na standardzie: średnia struktura ludności Polski 2010–2019",
    x = "Rok",
    y = "Zgony na 100 000",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom"
  )

print(p_trend)

# ------------------------------------------------------------
# 16.3. Scatter: crude vs ASMR
# ------------------------------------------------------------

p_scatter <- ggplot(
  df,
  aes(x = crude_mort_rate, y = asmr_2010_2019std)
) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Powiat-rok: surowy vs standaryzowany współczynnik zgonów",
    subtitle = "Linia przerywana oznacza y = x",
    x = "Surowy współczynnik zgonów na 100 000",
    y = "ASMR na 100 000"
  ) +
  theme_minimal(base_size = 13)

print(p_scatter)

# ------------------------------------------------------------
# 16.4. Czy ASMR jest mniej zależny od udziału osób 65+?
# ------------------------------------------------------------

corr_by_year <- df %>%
  group_by(year) %>%
  summarise(
    corr_crude_pop65 = cor(crude_mort_rate, pop_65, use = "complete.obs"),
    corr_asmr_pop65 = cor(asmr_2010_2019std, pop_65, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(corr_crude_pop65, corr_asmr_pop65),
    names_to = "metric",
    values_to = "correlation"
  ) %>%
  mutate(
    metric = recode(
      metric,
      corr_crude_pop65 = "Surowy vs udział 65+",
      corr_asmr_pop65 = "ASMR vs udział 65+"
    )
  )

p_corr <- ggplot(
  corr_by_year,
  aes(x = year, y = correlation, color = metric)
) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  labs(
    title = "Korelacja z udziałem ludności 65+",
    subtitle = "Po standaryzacji korelacja z pop_65 powinna zwykle spaść",
    x = "Rok",
    y = "Korelacja Pearsona",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom"
  )

print(p_corr)

# ------------------------------------------------------------
# 16.5. Różnica ASMR - crude
# ------------------------------------------------------------

df_sanity <- df %>%
  mutate(
    asmr_minus_crude = asmr_2010_2019std - crude_mort_rate,
    asmr_ratio_crude = asmr_2010_2019std / crude_mort_rate
  )

diff_summary <- df_sanity %>%
  summarise(
    diff_min = min(asmr_minus_crude, na.rm = TRUE),
    diff_p25 = quantile(asmr_minus_crude, 0.25, na.rm = TRUE),
    diff_median = median(asmr_minus_crude, na.rm = TRUE),
    diff_p75 = quantile(asmr_minus_crude, 0.75, na.rm = TRUE),
    diff_max = max(asmr_minus_crude, na.rm = TRUE),
    ratio_min = min(asmr_ratio_crude, na.rm = TRUE),
    ratio_median = median(asmr_ratio_crude, na.rm = TRUE),
    ratio_max = max(asmr_ratio_crude, na.rm = TRUE)
  )

cat("\nRóżnica ASMR - crude:\n")
print(diff_summary)

p_diff_hist <- ggplot(
  df_sanity,
  aes(x = asmr_minus_crude)
) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Rozkład różnicy: ASMR minus surowy współczynnik",
    subtitle = "Wartości dodatnie: po standaryzacji współczynnik jest wyższy niż surowy",
    x = "ASMR - crude, na 100 000",
    y = "Liczba obserwacji powiat-rok"
  ) +
  theme_minimal(base_size = 13)

print(p_diff_hist)

# ------------------------------------------------------------
# 16.6. Ranking powiatów, gdzie standaryzacja zmienia najwięcej
# ------------------------------------------------------------

rank_diff <- df_sanity %>%
  group_by(code, name) %>%
  summarise(
    mean_crude = mean(crude_mort_rate, na.rm = TRUE),
    mean_asmr = mean(asmr_2010_2019std, na.rm = TRUE),
    mean_diff = mean(asmr_minus_crude, na.rm = TRUE),
    mean_abs_diff = mean(abs(asmr_minus_crude), na.rm = TRUE),
    mean_pop65 = mean(pop_65, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_diff))

cat("\nPowiaty z największą średnią bezwzględną różnicą ASMR vs crude:\n")
print(rank_diff, n = 30)

p_rank <- rank_diff %>%
  slice_max(mean_abs_diff, n = 25) %>%
  mutate(
    name = reorder(name, mean_abs_diff)
  ) %>%
  ggplot(
    aes(x = name, y = mean_abs_diff)
  ) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Powiaty, w których standaryzacja zmienia współczynnik najbardziej",
    subtitle = "Średnia bezwzględna różnica ASMR vs crude",
    x = NULL,
    y = "Średnia |ASMR - crude| na 100 000"
  ) +
  theme_minimal(base_size = 13)

print(p_rank)

# ------------------------------------------------------------
# 16.7. Mapa różnicy ASMR - crude dla wybranego roku
# ------------------------------------------------------------

map_year <- 2023

df_map_diff <- df_sanity %>%
  filter(year == map_year) %>%
  mutate(
    code_pow = make_code_pow(code)
  ) %>%
  left_join(
    powiaty_spdf,
    by = c("code_pow" = "JPT_KOD_JE")
  ) %>%
  st_as_sf()

p_map_diff <- ggplot(df_map_diff) +
  geom_sf(aes(fill = asmr_minus_crude), color = NA) +
  scale_fill_gradient2(
    labels = comma,
    midpoint = 0
  ) +
  labs(
    title = paste0("Różnica ASMR - crude w ", map_year),
    subtitle = "Dodatnie wartości oznaczają, że standaryzacja podnosi współczynnik",
    fill = "ASMR - crude"
  ) +
  theme_void(base_size = 13)

print(p_map_diff)