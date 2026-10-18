# Zenodo dataset deposit — PREPARATION ONLY

**Status: draft metadata; no Zenodo deposit has been created. No DOI has been assigned.**

## Suggested title
County-level mortality and determinant panel for Poland, 2010–2025 (379 harmonised counties): analysis data supporting a spatial Mundlak GAM

## Recommended deposit structure

```text
study_data/
  finalny_panel_przestrzenny_379_powiatow.csv.gz (if redistribution permitted)
  finalny_panel_przestrzenny_379_powiatow.rds    (optional if permitted)
  county_boundaries.gpkg                        (only if redistribution permitted)
  variable_dictionary.csv                       (needs producing from actual data)
  SOURCES_AND_LICENSES.md                       (verify original terms)
  README_DATA.md                                (scope, units, construction, missingness)
```

The existing `svc_reporting` CSVs are reporting outputs, **not** the complete analysis dataset. Do not label them as a complete replication package.

## Publication sequence

1. Find the actual final panel locally under `data/`.
2. Check permissions for mortality, census, CAMS and geospatial redistribution.
3. Produce a variable dictionary and sources/licensing statement; inspect data for personal information.
4. Create a **draft** Zenodo record. Keep it unsubmitted while reviewing anonymity and rights; adjust access restrictions as appropriate.
5. Upload permitted materials and review content/checksums.
6. Publish when authorised, then record the **actual** DOI in GitHub and the SSM Data and Code Availability section.

Researcher identity and ORCID should be populated by the account owner in Zenodo; this draft does not assume or upload unverified author metadata.
