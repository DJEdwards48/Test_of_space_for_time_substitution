# Test_of_space_for_time_substitution
Code used for conducting tests of space-for-time substitution which uses Bayesian mixed effect models to compare between space-for-time substitution and time-for-time substitution across zooplankton diversity (alpha and local contribution to beta diversity) calculated under the Hill's diversity framework and 3 environmental variables (temperature, total nitrogen, and total phosphorus)

## Data files

### Copy.of.lake.coords.landformregions.csv
Contains the land region data for lakes. (Site rows and variable columns), including: Unique Lakes (Lake Identifiers); lat (latitude); long (longitude); LANDFORM_R (Land region of Iowa that site is located within).

### LakeAreas.csv
Contains information related to each lakes' size. (site rows and variable columns), including: LakeID (unique lake identifiers); TMDL(Y/N) (Whether documentation on lake with a 'total maximum daily load' data file was present or not (Y = present, N = absent)); Surface Area (Acre) (surface area of each lake); Watershed_LakeArea_Ratio (watershed to lake area ratio (acres)); Residence time(days) (Lake residence time in days)

### WQ.taxa25.csv
Contains information related to each lake's community composition and water quality variables. (site rows and variable columns), including: lake.date (unique Identifier of Lakes and sample date, lake unique ID_Year sampled_day of the year sampled); the following columns are zooplankton taxa biomass estimates (microgram/liter): Keratella, Cyclopoida, Nauplii, Daphnia,	Calanoida,	Brachionus,	Bosmina,	Ceriodaphnia,	Filinia,	Polyarthra,	Macrochaetus,	Hexarthra,	Kellicottia,	Lecane,	Trichocerca,	Asplanchna,	Monostyla,	Conochilus,	Pompholyx,	Ascomorpha,	Euchlanis,	Chydorus,	Leptodora,	Synchaeta,	Notholca,	Testudinella,	Lepadella,	Disparalona,	Anuraeopsis,	Proales,	Moina,	Gastropus,	Diaphanosoma,	Mytilina,	Platyias,	Alona,	Pseudochydorus,	Asplanchnopus,	Coulrella,	Collotheca,	Pleuroxus,	Alonella,	Trichotria,	Graptoleberis,	Camptocercus,	Ploesoma,	Simocephalus,	Scapholeberis; Temp.water.degC (surface water temperatures in celcius);	Kjeldahl_nitrogen.as_N. (total Kjeldahl nitrogen in mg/l); Total_Phorphorus_Mixed_Forms_mg.l.as.P.	(total phosphorus in mg/l); Orthophosphate_.as_P. mg.l (Orthophosphate in mg/l);	Nutrient.nitrogen_mg.l (alternative name for Kjeldahl Nitrogen, resolved in code in mg/l)	Nitrate.Nitrite_.as_N._mg.l (Combined values of Nitrate and Nitrite in mg/l) 

### Sft-tft-wholecomparison.R
R code used to conduct test space-for-time substitution assumptions with all data 

### Sft-tft-bylakesize.R

### Sft-tft-testbylandregion.R
 
