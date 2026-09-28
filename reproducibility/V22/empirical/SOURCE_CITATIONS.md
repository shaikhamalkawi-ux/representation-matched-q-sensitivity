# Source and encoding map

All data tables here are author-created adaptations of the cited aggregate data;
they are not endorsed or official rankings of any source agency. Original-source
fingerprints and exact derived-table hashes appear in `PROVENANCE.json`.

| ID | Tables | Source and exact analyst encoding |
|---|---|---|
| EC_QOL2023 | qol_primary (83x2), qol_secondary (83x10) | European Commission, DG REGIO, 2023 Quality of Life in European Cities, workbook revised September 2024. Close the five reported response shares by their total T; mu=(rather+very satisfied)/T, nu=(rather+very unsatisfied)/T; nonresponse remains residual. All 83 cities retained. Closure changes are below 5.4e-15. |
| CMS_HCAHPS2026 | cms_primary (3168x2), cms_strict (3134x2) | CMS HCAHPS, August 13, 2026 release, October 1, 2024–September 30, 2025 reporting period. Nurse/doctor communication: Always/100 and Sometimes-or-Never/100; Usually remains residual. 4790 units screened, six unique measures, common dates, numeric rates/counts, blank exclusion footnotes; no 300-response threshold. Strict removes 34 boundary cases. Public adjusted survey summaries, not binomial patient data. |
| OECD_PISA2022 | pisa_primary (77x2), pisa_secondary (65x2) | OECD PISA 2022 Results Volume I, Tables I.B1.3.1/2, mathematics/reading. mu=Levels5+6/100, nu=BelowLevel2/100, no renormalization; total within 0.1 percentage point. Exclude absent positive tails and Viet Nam reading linkage; secondary also excludes 12 sampling-warning units. |
| CAA_PUNCTUALITY2025 | caa_primary (23x2) | UK CAA 2025 flight punctuality, arrivals/departures. For each stratum M=matched,U=unmatched,C=cancelled,D=M+U+C; denominator N=sum(M+C). mu=sum(D*g)/(100N), where g is early through 15 minutes late; nu=(sum(D*b)/100+sum C)/N, b more than 60 minutes late. Unmatched not treated as neutral. Four gross category-closure failures exclude whole Bournemouth/Isle of Man airports. Rounded bin shares are not flight counts. |
| EPA_AQI2024_2025 | epa_primary (741x2), epa_secondary (333x2) | US EPA annual county AQI, years 2024/2025; mu=GoodDays/reportedDays, nu=days AQI>=101/reportedDays, Moderate residual. Exact State|County matching; US states/DC, coverage >=330/366 and >=329/365; secondary full-calendar. County daily maximum over reporting sites/pollutants is not mean personal exposure. Download September 28, 2026, catalog June 25 versus HTTP modified July 13 both retained; not proven June25 bytes. |
| ACS_HOUSING2024 | housing_primary (51x2), housing_secondary (51x2) | US Census ACS 2024 one-year, B25070/B25091, renters/owners for states/DC. F=below30% income,A=50%or more,M=30–49.9%,U=notcomputed,T=F+M+A+U,C=T-U. Primary(F/C,A/C); secondary(F/T,A/T). Pool mortgage subgroups before division. Integer point estimates, not exact census counts; computability is not response rate; no joint margin-of-error uncertainty box. |

Official landing pages and methods:

- EC: https://ec.europa.eu/regional_policy/information-sources/maps/quality-of-life_en
- CMS: https://data.cms.gov/provider-data/dataset/dgck-syfz ; https://www.cms.gov/medicare/quality/initiatives/hospital-quality-initiative/hcahps-patients-perspectives-care-survey
- OECD: https://www.oecd.org/en/publications/pisa-2022-results-volume-i_53f23881-en/full-report/results-for-countries-and-economies_360c8f67.html
- CAA: https://www.caa.co.uk/data-and-analysis/uk-aviation-market/flight-punctuality/uk-flight-punctuality-statistics/2025/ ; https://www.caa.co.uk/data-and-analysis/uk-aviation-market/flight-punctuality/uk-flight-punctuality-statistics-notes/
- EPA: https://aqs.epa.gov/aqsweb/airdata/download_files.html ; https://www.epa.gov/outdoor-air-quality-data/about-airdata-reports
- ACS: https://www.census.gov/programs-surveys/acs/data/summary-file.2024.html ; https://www2.census.gov/programs-surveys/acs/tech_docs/subject_definitions/2024_ACSSubjectDefinitions.pdf

Operator literature: Liu and Wang (2018), https://doi.org/10.1002/int.21927
(WA/WG; author-posted full-text formula access); Khan et al. (2021),
https://doi.org/10.1007/s40747-021-00425-7 (cosine similarity Eq.20 and Algorithm1,
published PDF). WG correction: https://doi.org/10.1155/int/9864340.
No publisher or subscription PDF is redistributed here.
