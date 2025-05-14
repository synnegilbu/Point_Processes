# AVOTREX: a global dataset of extinct birds and their traits (v. 1.0)

**Authors:** Ferran Sayol, Joseph P. Wayman, Paul Dufour, Thomas E. Martin, Julian P. Hume, Maria Wagner Jørgensen, Natàlia Martínez-Rubio, Ariadna Sanglas, Filipa C. Soares, Rob Cooke, Chase D. Mendenhall, Jay R. Margolis, Juan Carlos Illera, Rhys Lemoine, Eva Benavides, Oriol Lapiedra, Kostas A. Triantis, Alex L. Pigot, Joseph A. Tobias, Søren Faurby & Thomas J. Matthews.

Comments and requests should be addressed to Ferran Sayol: [fsayol@gmail.com](mailto:fsayol@gmail.com).

### Description of the data and file structure:

Avotrex_metadata.xlsx : Spreadsheet file containing the definitions of all variables (with their corresponding units), abbreviations, and acronyms used within the dataset files. The document comprises the following sheets:

* Metadata_Avotrex: Description of variables from AVOTREX (v. 1.0) files
* Metadata_Avonet: Description of variables from AVONET database
* Metadata_RawData_D1-D8: Description of variables from RawData_D1-D8
* Ref_IDs: List of reference IDs used within the AVOTREX dataset files and their corresponding full references from which the information for the species in the database was extracted. Short codes were created using the following formula: 'first author surname' + underscore + 'year of publication'.
* Collection_IDs: List of collection ID acronyms and their corresponding full name, city, and location, from which the specimens measured for the AVOTREX database belong.
* Measurer_IDs: List of unique identifiers (initials or other abbreviation) and their corresponding full names, used to identify each individual measurer within the AVOTREX database.

Avotrex_taxonomy_v1.csv : Taxonomy of the 610 species that are part of the AVOTREX database.

Avotrex_traits_v1.csv : Status, distribution and external morphological measurements, either measured or imputed, of the 610 species that are part of the AVOTREX database.

DB_Avonet_BirdLife.csv : Contains the external morphological measurements of 11009 extant bird species from the Avonet database. For more information see Tobias et al 2022 Ecology Letters DOI: [https://doi.org/10.111/ele.13898](https://doi.org/10.111/ele.13898)

RawData_D1.csv : External morphological measurements for extinct species taken from preserved skins of specimens in natural history collections.

RawData_D2.csv : External morphological measurements for extinct species from preserved skins of specimens gathered from literature sources.

RawData_D3.csv : Skeletal measurements for extinct species taken from museum specimens.

RawData_D4.csv : Skeletal measurements for extinct species gathered from literature sources.

RawData_D5.csv : Body mass for extinct species gathered from literature sources.

RawData_D6.csv : Inferred body mass for extinct species based either on qualitative descriptions in the literature (i.e., reported to be a similar size to a specific close relative) or skeletal proportions (i.e., 10% larger than a relative based on femur measurements).

RawData_D7.csv : Skeletal measurements for extant species taken from museum specimens.

RawData_D8.csv : Skeletal measurements for extant species gathered from literature sources.

### Version changes:

**26-oct-2024:** First publication of the dataset.

**Jan-2025:** Removed `specimen_id` from several items in the CM collection in RawData_D3 and RawData_D7 due to errors.

### User notes:

Across all the spreadsheet files, any cells containing "NA" indicate that data is not available for that particular trait or variable.

The datasets included in this repository are part of the AVOTREX (v. 1.0) database.

If you find any errors, you can report them here to help update future versions of AVOTREX: [https://forms.gle/HWPoAri1gyNzBENt8](https://forms.gle/HWPoAri1gyNzBENt8)

### Sharing/Access information:

Full reference: Sayol, F., Wayman, J., Dufour, P., Martin, T.E., Hume, J.P., Jørgensen, M., Martínez-Rubio, N., Sanglas, A., Soares, F.C., Cooke, R., Mendenhall, C.D., Margolis, J.R.,Illera, J.C., Lemoine, R., Benavides, E., Lapiedra, O., Triantis, K.A., Pigot, A.L., Tobias, J.A., Faurby, S. & Matthews, T. J. (2024). AVOTREX: a global dataset of extinct birds and their traits. *Global Ecology and Biogeography,* e13927. [https://doi.org/10.1111/geb.13927](https://doi.org/10.1111/geb.13927)
