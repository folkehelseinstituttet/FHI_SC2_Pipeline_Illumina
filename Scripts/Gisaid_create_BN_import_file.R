# Define clean-up function ------------------------------------------------
clean_up <- function(submission, HF){
  if (HF == "FHI") {
    #Not for OUS-data
    to_BN <- submission %>% 
      separate(V1, into = c(NA, NA, "unique", "year"), sep = "/") %>% 
      add_column("prefix" = 25) %>% 
      mutate("year" = str_sub(year, start = 3, end = -1)) %>% 
      # Create the Key for BN matching
      unite(Key, c(prefix, year, unique), sep = "", remove = F) %>% 
      mutate("GISAID_EPI_ISL" = str_sub(V2, start = 1, end = -3)) %>% 
      # Re-create the virus name to match with the Gisaid metadata
      add_column("Gisaid_prefix" = "hCoV-19/Norway",
                 "year_full" = 2021) %>% 
      unite("covv_virus_name", c(Gisaid_prefix, unique, year_full), sep = "/") %>% 
      select(Key, GISAID_EPI_ISL, covv_virus_name)
    return(to_BN)
  }else if (HF == "OUS") {
    # Only for OUS data
    to_BN <- submission %>% 
      separate(V1, into = c(NA, NA, "Key", NA), sep = "/") %>% 
      add_column("year" = 2021) %>% 
      mutate("GISAID_EPI_ISL" = str_sub(V2, start = 1, end = -3)) %>% 
      # Re-create the virus name to match with the Gisaid metadata
      add_column("Gisaid_prefix" = "hCoV-19/Norway") %>% 
      unite("covv_virus_name", c(Gisaid_prefix, Key, year), sep = "/", remove = F) %>% 
      select(Key, GISAID_EPI_ISL, covv_virus_name)
    return(to_BN)
  }else {
    print("HF needs to be either FHI or OUS")
  }
}


# Start script ------------------------------------------------------------

library(tidyverse)

try(setwd("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/4-GISAIDsubmisjon/"))
try(setwd("N:/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/4-GISAIDsubmisjon/"))

# Set the path for the submission directory
path <- "2021-12-06_Run681/"

# Oppsett name
oppsett <- "2021-12-10_Run681"

# Read the gisaid uploader log file
submission <- read_lines(list.files(path = path, 
                      pattern = ".log$",
                      full.names = TRUE))
# keep only sucessfully uploaded Keys
submission <- submission[str_detect(submission, "error", negate = TRUE)]
submission <- submission[str_detect(submission, "EPI_ISL")]
# Get Gisaid acc. and Virus name on different columns
submission <- as_tibble(as.data.frame(str_split(submission, pattern = ";", simplify = TRUE)))


# Clean up data and create columns for import -----------------------------
# HF is either "FHI" or "OUS"
to_BN <- clean_up(submission, HF = "FHI")


# Get the Gisaid metadata -------------------------------------------------
files <- list.files(path = path, 
                    pattern = "*.csv",
                    full.names = TRUE)
files <- files[-grep("failed", files)]

if (length(files) == 1) {
  gisaid_metadata <- read_csv(files,
                              col_types = cols(covv_seq_technology = col_character())) # Add this because comma in the text for Nanopore
} else {
  print("There can only be one csv file")
}



# Merge seq. tech. info to BN import data ---------------------------------
to_BN <- left_join(to_BN, gisaid_metadata, by = "covv_virus_name") %>% 
  select(Key, 
         "gisaid_epi_isl" = GISAID_EPI_ISL, 
         "Platform" = covv_seq_technology)


# Write file --------------------------------------------------------------
write_csv(to_BN, file = paste0(path, oppsett, "_BN_import.csv"))
