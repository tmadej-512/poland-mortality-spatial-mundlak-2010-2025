# Data

## Outcome

Age-sex standardised all-cause mortality rates (ASMR) per 100,000 population for 379 Polish counties (powiaty) over 2010–2023. Death counts at the county level were provided directly by Statistics Poland (GUS) and are **not publicly available** at the required level of disaggregation for all study years. An anonymised version may be requested from Statistics Poland.

Population denominators were obtained from the Statistics Poland Local Data Bank (Bank Danych Lokalnych): https://bdl.stat.gov.pl

## Predictors

All eight time-varying determinants and three between-only determinants were derived from publicly available sources:

| Variable | Source |
|---|---|
| Gross wages (PLN) | Statistics Poland BDL |
| Unemployment rate (%) | Statistics Poland BDL |
| Physicians per 100,000 | Statistics Poland BDL |
| Hospital beds per 100,000 | Statistics Poland BDL |
| Health expenditure per capita (PLN) | Statistics Poland BDL |
| Population density (persons/km²) | Statistics Poland BDL |
| Net migration rate (per 1,000) | Statistics Poland BDL |
| Annual mean PM2.5 (µg/m³) | Chief Inspectorate of Environmental Protection (GIOŚ) |
| Higher education share (%) | 2021 National Population and Housing Census (GUS) |
| Low education share (%) | 2021 National Population and Housing Census (GUS) |
| Distance to nearest city with county rights (km) | Computed from TERYT administrative boundaries (GUS) |

## Spatial data

County (powiat) boundary shapefiles are available from the Head Office of Geodesy and Cartography (GUGiK): https://www.geoportal.gov.pl

The analysis used the 2019 NUTS-4 county boundaries harmonised to a consistent 379-unit panel (Wałbrzych county was aggregated due to an administrative status change during the study period).

## Expected file structure

```
data/
  asmr_panel.csv          # NOT included — request from GUS
  determinants_panel.csv  # Compiled from BDL and GIOŚ (publicly available)
  powiaty_shapefile/      # Download from GUGiK
```

Scripts expect data to be loaded via `dm.R` (data management script). See `R/dm.R` for variable construction and standardisation details.
