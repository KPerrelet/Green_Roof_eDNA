# Green_Roof_eDNA
Data and code for the article titled "Green roofs harbor different and non-substituting invertebrate communities than surrounding ground-level habitats"


The information contained in this repository is referenced in the manuscript. 


*****************************************************


skyline_controlmat.csv - OTU table for the control samples (2022 and 2023). 

Rows = sites, columns = OTU. Values indicate number of reads per OTU (after filtering, contamination removal, and rarefaction). 


skyline_covmat_greenroof.csv - Covariate values extracted for all green roof sites. As many sites are private, their location cannot be disclosed publically.

skyline_covmat_groundlevel.csv - Covariate values extracted for all ground-level sites. As many sites are private, their location cannot be disclosed publically.

skyline_shareddist.csv - Vertical (building height) and horizontal (distance to ground-level site) isolation metrics for each pairs of site


skyline_sharedmat.csv - Taxa table for the shared taxa. 

Rows = sites, columns = species. Values indicate number of reads per OTU (after filtering, contamination removal, and rarefaction). 

/!\ This table also contains results extracted from null models and can be recreated using the null_models.R script /!\


skyline_spmat.csv - Taxa table for all sites. 

Rows = sites, columns = taxa. Values indicate number of reads per OTU (after filtering, contamination removal, and rarefaction). 


skyline_spmat_per_sample.csv - Taxa table for all samples (3 samples per site). 

Rows = samples, columns = taxa. Values indicate number of reads per OTU (after filtering, contamination removal, and rarefaction). 


skyline_taxo.csv - Taxonomy table, also containing mobility mode values, for all taxa. 

Rows = taxa, columns = phylogeny and mobility mode. 


Raw sequences for the data generated in this study will be uploaded on the European Nucleotide Archive upon acceptance.


*****************************************************


fig2.R - generate Figure 2. 

fig3_figS4.R - generate Figure 3 and Figure S4 (which is an alternative to figure 2, but separated per site). 

fig4.R - generate Figure 4. 

fig5.R - generate Figure 5. 

figS1.R - generate Figure S1. 

figS2.R - generate Figure S2. 

figS3.R - generate Figure S3. 

figS5.R - generate Figure S5. 

figS6.R - generate Figure S6. 

figS7.R - generate Figure S7. 

figS8.R - generate Figure S8. 

figS9.R - generate Figure S9. 

figS10.R - generate Figure S10. 

null_models.R - generate null models for figure 5, as well as compute the shared species table (i.e., skyline_sharedmat.csv)

*****************************************************

Data dictionary: 

site_no          |  paired site ID. A green roof site has the same ID as the ground-level site only if they are paired. 

type             |  site type ("Roof" = green roof, "Ground" = ground-level)

vegetation       |  vegetation coverage (%)

substrate type   |  most dominant substrate type (Soil/Stones)

height           |  building height

age_roof         |  Age of the roof

depth            |  substrate depth

grey_frac_50     |  fraction of grey surface in a 50 m buffer

grey_frac_500    |  fraction of grey surface in a 500 m buffer

green_frac_50    |  fraction of green surface in a 50 m buffer

green_frac_500   |  fraction of green surface in a 50 m buffer

pv               |  presence of solar panels (0 = FALSE, 1 = TRUE)

distance_gr      |  distance to the closest green roof

paired_distance  |  distance between paired sites (one pair = one green roof + one ground-level site)

year             |  year in which the samples were taken
