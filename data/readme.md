### Code and data from Kandlikar et al. "Functional traits predict species responses to environmental variation in a California grassland annual plant community"

**all_environmental_data.csv**: This file contains all the environmental data collected for our study.  

`ele` refers to Elevation (meters ASL);

`Tmax` and `Tmin` refer to the average daily maximum and minimum temperature during the growing period as recorded using iButtons (degrees C);

The columns `organic_matter_ENR`, `pH`, `CEC_meq_100g`, `K_ppm`, `Mg_ppm`, `Ca_ppm`, `NH4_N_ppm`, `Nitrate_ppm`,  `sand`, `clay`,  `P1_weak_bray_ppm`, and `NaHCO3P_olsen_ppm` refer to soil chemical and physical characteristics measured by A&L agricultural labs (http://www.al-labs-plains.com/soil/2511974). `ppm` = Parts per million; `meq` = Milliequivalents per 100 grams; the `sand` and `clay` columns are in units of `%`.

`microsite` refers to whether a plot was situated on a serpentine hummock or in the grassland matrix.

`depth` is reported in units if `cm`

`soil_moisture` is the average gravimetric water content ((weight of fresh soil - weight of dry soil)/weight of dry soil) in early- and mid- growing season (March and April, respectively)


**all_trait_data.csv**  

`USDA_code` is the abbreviated USDA code for each species
`leaf size` is reported in cm^2; 
`SLA` is specific leaf area, reported in (cm^2/g); 
`LDMC` is leaf dry matter contet, reported in (mg/g); 
`seed_mass` is reported in mg; 
`max_height` is reported in cm; 
`SRL` is specific root length, reported in (m/g);
`relative_spread` is the ratio of a plant's lateral spread vs. height (dimensionless); 
`projected_area` is the average projected canopy area of a plant, as area of oval with each plant's longest axis and axis perpendicular to that providing the orthogonal dimensions, reported as cm^2; 
`phenology` is day of year at which approximately 50\% of individual plants begin fruiting in the field (doy);
`foliar_N` is the leaf Nitrogen as a percent of total leaf mass (%),;  
`CN_ratio` is the ratio of leaf carbon to leaf nitrogen (unitless);
`d13C` is the relative enrichment of leaf carbon with carbon-13 isotopes;
`d15N` is the relative enrichment of leaf nitrogen with nitrogen-15 isotopes

**germination_data.csv** 

`plot_num` is the plot number (740-755), used to keep track of different sites; 
`rep` is the replicate planting at each plot (1-5);  
`sp_code` is the species code (see Table 1 of manuscript for full names);  
`num_germinants` is the number of germinants recorded across three germination surveys in mid-Jan to mid-February;   
`expected_germinants_l` is the number of viable seeds planted at each site (expected number of germinants);  
`num_germinants_corrected` is identical to `num_germinants`, except in cases where `num_germinants` > `expected_germinants_l` (i.s. where we saw more germinants than the number of seeds we planted, presumably due to germination from the background seedbank.) In these cases, we set `num_germinants` = `expected_germinants_l`, such that germination fraction = 1.   

**sedgwick_cover_2017.csv**  

`year`: year of cover estimate
`site`: site number, 740 to 755
`plot`: plot number, 1 to 5
`USDA_symbol`: USDA symbol used to identify species
`area`: Sample area in cm^2
`cover`: estimated cover of each species in each plot
