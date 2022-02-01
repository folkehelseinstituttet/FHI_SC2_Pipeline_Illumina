metadata <- read_tsv(file = "/home/jonr/Downloads/metadata.tsv")

metadata <- metadata %>% 
  filter(str_detect(Location, "Norway")) %>% 
  select(`Virus name`, `Accession ID`)

to_BN <- metadata %>% 
  mutate("GISAID_ISOLAT_NAVN" = str_remove(`Virus name`, "hCoV-19/")) %>% 
  select("gisaid_epi_isl" = `Accession ID`,
         "GISAID_ISOLAT_NAVN")

write_csv(to_BN, 
          file = "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/4-GISAIDsubmisjon/2022.02.01_BN_import_Gisaid_isolate.csv",
          quote = "all")
