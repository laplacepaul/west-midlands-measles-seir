# Data provenance and limits

**Source:** UK Health Security Agency dashboard API, metric `measles_cases_casesByOnsetWeek`, region West Midlands (`E45000005`), all age groups and sexes, retrieved October 2026. Source: https://ukhsa-dashboard.data.gov.uk/metrics-documentation/measles-cases-by-week-of-symptom-onset

The extracted source counts are laboratory-confirmed measles cases grouped by symptom-onset week; they are not total infections, infection dates, or an independent measure of transmission.

**Fitting subset:** 6 November 2023 through 29 April 2024 (26 consecutive weeks). A row for 23 October 2023 is absent from the API extract and is never treated as a known zero. The entire studied October 2023–December 2024 case series has 63 recorded rows across 64 expected weekly slots; the observed records sum to 714. UKHSA revisions mean other retrieval dates may differ.

**Additional context, not used as direct parameter estimates:** ONS/Nomis mid-year West Midlands population estimates (2023: 6,105,177; 2024: 6,181,633), and NHS England MMR coverage figures for childhood cohorts. Those vaccination percentages cannot be converted directly into an all-age immune proportion.

For geography comparisons distinguish UKHSA region `E45000005` and ONS statistical region `E12000005`; correspondence was examined by administrative coverage, not assumed from identical geography codes.

## Licence and attribution

The UKHSA data dashboard states in its site footer: "All content is available under the Open Government Licence v3.0, except where otherwise stated", with the notice "© Crown copyright" (checked 25 November 2026 at https://ukhsa-dashboard.data.gov.uk/about). No dataset-specific exception was found for the measles onset-week metric used here.

The two files in `data/` derived from that source are redistributed under the same licence, with this attribution:

> Contains public sector information licensed under the Open Government Licence v3.0. Source: UK Health Security Agency, UKHSA data dashboard, metric `measles_cases_casesByOnsetWeek`, West Midlands (`E45000005`), retrieved October 2026. © Crown copyright.

Licence text: https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/

Conditions to keep when reusing these files: acknowledge the source as above, link to the licence where possible, and do not imply endorsement by UKHSA. The processed file `observed_weekly_cases_with_missingness.csv` is an adaptation (subset, explicit missing-row flag, week index); its values are not altered.

The ONS/Nomis population figures and NHS England MMR coverage figures are quoted as context in the report only; their data files are not redistributed here. Check the licence statement on each source page before redistributing those data.
