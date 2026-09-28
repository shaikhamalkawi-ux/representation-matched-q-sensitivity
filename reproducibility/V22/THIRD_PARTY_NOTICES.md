# Third-party sources, attribution and licence boundaries

Source-policy pages checked on 28 September 2026. These notices identify
applicable source conditions; they are not legal advice, a rights warranty or
an assertion that the authors own third-party data.

The repository's existing CC0 dedication applies only to the authors' original
material, to the extent they hold the relevant rights. It does not relicense
the source-derived input tables or override the conditions below. Keep these
notices and the source/provenance record with redistributed copies. Public
access to data is not, by itself, an unrestricted licence.

The six empirical domains are European city satisfaction (`qol`), hospital
communication (`cms`), educational proficiency (`pisa`), airport punctuality
(`aviation`), county air quality (`air_quality`), and housing burden (`housing`).
The release uses selected aggregate values and derived rational proportions,
not a new release of respondent-level data or complete official databases.
Cohort exclusions, category grouping and mathematical encodings are the
research authors' adaptations, not methods endorsed by the source agencies.

## European Commission: quality of life in European cities, 2023

Source: European Commission, Directorate-General for Regional and Urban
Policy, *Report on the quality of life in European cities, 2023*, accompanying
[aggregate city data](https://ec.europa.eu/regional_policy/information-sources/maps/quality-of-life_en).
The retained data are transformed into selected favourable/adverse shares;
this adaptation and its results are the research authors' responsibility.

The [Commission copyright notice](https://commission.europa.eu/legal-notice_en)
places EU-owned website content under CC BY 4.0 unless otherwise indicated,
requiring credit and identification of changes. Third-party content and logos
are not covered by that default. For material obtained from Eurostat, its
[reuse notice](https://ec.europa.eu/eurostat/help/copyright-notice) also requires
source acknowledgement, disclosure of modifications and a disclaimer of
Eurostat responsibility; specific-source exceptions remain applicable. This
artifact grants no additional rights over the source statistics. Neither the
European Commission nor Eurostat is responsible for or endorses this analysis.

## CMS: HCAHPS hospital-level aggregate data

Source: U.S. Centers for Medicare & Medicaid Services, [HCAHPS Hospital data,
Provider Data Catalog, dataset dgck-syfz](https://data.cms.gov/provider-data/dataset/dgck-syfz).
The mathematical input uses publicly reported hospital-level response
percentages, with the source qualifications and exclusion rules recorded in
the provenance files; it contains no patient-level observations.

The Provider Data Catalog's [government-data notice](https://data.cms.gov/provider-data/topics/long-term-care-hospitals/about-using-government-data)
states that U.S. Government works are public domain, appreciates agency
attribution and warns against implying government endorsement. This notice is
not a claim that every third-party item available through CMS shares the same
status. CMS has not approved or endorsed the encoding, rankings or conclusions.

## OECD: PISA 2022 proficiency distributions

Source: OECD (2023), *PISA 2022 Results (Volume I): The State of Learning and
Equity in Education*, [country/economy results and associated tables](https://www.oecd.org/en/publications/pisa-2022-results-volume-i_53f23881-en/full-report/results-for-countries-and-economies_360c8f67.html),
including the [StatLink workbook](https://stat.link/znxau8), accessed for this
analysis on 28 September 2026. Selected mathematics/reading proficiency levels
are grouped into favourable/adverse/residual categories; these are adaptations,
not official OECD scores or recommended decision rankings.

Section 3 of the [OECD terms and conditions](https://www.oecd.org/en/about/terms-conditions.html)
permits data extraction, adaptation and distribution subject to attribution
and any dataset-specific or third-party restrictions. It requires the source
acknowledgement requirement to be passed on when work using the data is shared
or sublicensed, including to further sublicensees. Keep the OECD citation and
this acknowledgement condition with this derived input and its redistribution.
The source data are not dedicated to CC0 by this repository. No OECD logo is
provided, and the encoding and conclusions do not represent OECD or member
country views or endorsement. Terms for OECD written publications are not
silently substituted for the separate data terms.

## UK Civil Aviation Authority: 2025 punctuality data

Source: UK Civil Aviation Authority, [2025 Annual Punctuality Statistics, Full
Analysis Arrival Departure](https://www.caa.co.uk/data-and-analysis/uk-aviation-market/flight-punctuality/uk-flight-punctuality-statistics/2025/).
The released numerical input consists of the research authors' selected
airport-level aggregates, not the complete airline/route dataset. Preserve the
source definition, conditional denominator, excluded airports and coverage
qualifications documented in its provenance.

The source page requires CAA acknowledgement and states: “No statistical data
provided by CAA may be sold on to a third party.” Preserve that condition; this
repository does not authorize resale of the underlying CAA data. No CC0 or
Open Government Licence is asserted for the CAA-derived input. CAA disclaims
accuracy/reliability warranties and liability for reliance on its statistics.
The [methodological notes](https://www.caa.co.uk/data-and-analysis/uk-aviation-market/flight-punctuality/uk-flight-punctuality-statistics-notes/)
identify cooperation with airport operators and Airport Coordination Ltd and
the omission of operators lacking publication consent. Do not infer rights in
unpublished underlying records or imply CAA endorsement of this model.

## U.S. EPA: AirData/AQS county AQI aggregates

Source: U.S. Environmental Protection Agency, AirData/Air Quality System,
[annual county AQI files for 2024 and 2025](https://aqs.epa.gov/aqsweb/airdata/download_files.html).
The source counts are grouped and filtered for the declared coverage cohorts;
the resulting mathematical scores are not EPA health advice or official
county rankings.

EPA's [specific AirData permission FAQ](https://www.epa.gov/outdoor-air-quality-data/do-i-need-request-permission-use-monitoring-data-and-graphics-airdata)
identifies ambient AQS monitoring data as public domain and permits their use
without a permission request. Credit EPA and preserve dates, coverage and
methodological qualifications. No EPA endorsement or warranty of the
adaptation is implied. This specific data notice does not relicense unrelated
third-party documents or imagery hosted on EPA websites.

## U.S. Census Bureau: 2024 ACS housing tables

Source: U.S. Census Bureau, [2024 ACS one-year table-based Summary File](https://www.census.gov/programs-surveys/acs/data/summary-file.2024.html),
Detailed Tables B25070 and B25091. The inputs contain derived state/DC aggregate
ratios, not confidential household records. Their grouping and alternative
denominators are research-author choices; no Census approval or endorsement is
claimed.

The Bureau's [Research Transparency and Public Access policy, DS027, page 9](https://www2.census.gov/foia/ds_policies/ds027.pdf)
explains that employee-created data/works generally lack U.S. copyright
protection, while foreign protection and nonemployee rights can differ. This
repository does not enlarge that statement into a universal worldwide CC0
dedication of Census material. Preserve source attribution and the published
estimates' statistical limitations; do not use the derived ratios as validated
uncertainty intervals or identify individual respondents.

## Publications, software dependencies and other material

Published mathematical definitions are cited in the claim/source map. This
release does not redistribute subscription/reference PDFs, journal figures,
publisher page images, logos or manuscript PDFs. The source provenance of a
numerical table is not a claim to copyright in the corresponding publication.

Third-party Python libraries are separately installed dependencies and retain
their own licences. No copied binary dependency or installed environment is
relicensed by this repository. Earlier public versions and their original
licensing records remain preserved. No author-approval statement, peer review,
journal acceptance or journal submission is implied by publication of this
reproducibility artifact.
