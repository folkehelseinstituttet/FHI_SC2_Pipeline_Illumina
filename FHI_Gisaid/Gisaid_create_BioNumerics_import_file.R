pacman::p_load(tidyverse)
# Download from Gisaid Sequencing Technology metadata

# Read the metadata file
gisaid_md <- read_tsv(file = "/home/jonr/Downloads/gisaid_hcov-19_2022_02_01_09.tsv")

to_BN_OUS <- gisaid_md %>% 
  select(`Virus name`, `Accession ID`, `Sequencing technology`) %>% 
  filter(str_detect(`Virus name`, "OUS")) %>% 
  separate(`Virus name`, into = c(NA, NA, "Key", NA), sep = "/") %>% 
  select(Key, 
         "gisaid_epi_isl" = `Accession ID`,
         "Platform" = `Sequencing technology`)

write.csv(to_BN_OUS, 
          file = "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/4-GISAIDsubmisjon/BN-import-filer/2022.02.01_MIK_batch_import.csv",
          quote = TRUE,
          row.names = FALSE)

to_BN_rest <- gisaid_md %>% 
  # Remove Ahus samples
  filter(str_detect(`Virus name`, "Ahus", negate = TRUE)) %>% 
  select(`Virus name`, `Accession ID`, `Sequencing technology`) %>% 
  # Remove OUS samples
  filter(str_detect(`Virus name`, "OUS", negate = TRUE)) %>% 
  separate(`Virus name`, into = c(NA, NA, "Key", "year"), sep = "/") %>%
  # Left pad the Key with zeroes to a total of 5 digits
  mutate("Key" = str_pad(Key, width = 5, side = c("left"), pad = "0")) %>% 
  add_column("nr" = 25) %>% 
  mutate("year" = str_sub(year, 3, 4)) %>% 
  unite("tmp", c(nr, year, Key), sep = "") %>% 
  select("Key" = tmp,
         "gisaid_epi_isl" = `Accession ID`,
         "Platform" = `Sequencing technology`)

write.csv(to_BN_rest, 
          file = "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/4-GISAIDsubmisjon/BN-import-filer/2022.02.01_BN_batch_import.csv",
          quote = TRUE,
          row.names = FALSE)
