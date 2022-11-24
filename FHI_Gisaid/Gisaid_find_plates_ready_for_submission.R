#############################################
## Maintaned by Jon Bråte - jon.brate@fhi.no
## This script checks which SARS-CoV-2 samples at NIPH have not been submitted to Gisaid.
#############################################

# Load packages
pacman::p_load(tidyverse, readxl, stringr, lubridate)

# Read data from BioNumerics ----------------------------------------------
# NB: Husk å kjøre scriptet refresh_data_from_BN.R først
try(load(file = "/mnt/N/Virologi/JonBrate/Prosjekter/BN.RData"))
try(load(file = "N:/Virologi/JonBrate/Prosjekter/BN.RData"))


# Check Swift FHI
BN %>%
  # Convert empty strings to NA
  mutate_all(list(~na_if(.,""))) %>%
  # Get the Swift FHI plates
  filter(str_detect(SEKV_OPPSETT_SWIFT7, "FHI")) %>%
  # Filtrer på coverage >= 94%
  filter(Dekning_Swift >=94) %>%
  # Format the date
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>%
  # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
  filter(!is.na(MELDT_SMITTESPORING)) %>%
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  filter(is.na(GISAID_EPI_ISL)) %>% 
  # Get the various plates)
  select(SEKV_OPPSETT_SWIFT7) %>%
  distinct() %>%
  # Create numeric column for sorting
  mutate("tmp" = as.numeric(str_remove(SEKV_OPPSETT_SWIFT7, "FHI"))) %>% 
  arrange(desc(tmp)) %>% 
  select(SEKV_OPPSETT_SWIFT7) %>% 
  View("Swift_FHI")

# Check Artic Illumina
BN %>%
  # Convert empty strings to NA
  mutate_all(list(~na_if(.,""))) %>% 
  # Filtrer på coverage >= 94%
  filter(Dekning_Artic >=94) %>%
  # Format the date
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>%
  # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
  filter(!is.na(MELDT_SMITTESPORING)) %>%
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  filter(is.na(GISAID_EPI_ISL)) %>%
  # Get the various plates
  select(SAMPLE_CATEGORY) %>%
  # Remove some wrong entries
  filter(str_detect(SAMPLE_CATEGORY, "naekstrakt", negate = TRUE)) %>% 
  filter(str_detect(SAMPLE_CATEGORY, "Nano", negate = TRUE)) %>% 
  # Create column for sorting
  separate(SAMPLE_CATEGORY, into = c(NA, "tmp"), sep = "/", remove = FALSE) %>% 
  mutate("tmp" = as.numeric(tmp)) %>% 
  arrange(desc(tmp)) %>% 
  select(SAMPLE_CATEGORY) %>% 
  distinct() %>%
  View("Artic_Illumina")

# Check Nanopore
BN %>%
  # Convert empty strings to NA
  mutate_all(list(~na_if(.,""))) %>%
  # Filtrer på coverage >= 94%
  filter(Dekning_Nano >=94) %>%
  # Format the date
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>%
  # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
  filter(!is.na(MELDT_SMITTESPORING)) %>%
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  filter(is.na(GISAID_EPI_ISL)) %>%
  # Get the various plates
  select(SEKV_OPPSETT_NANOPORE) %>%
  distinct() %>% 
  # Remove a wrong entry
  filter(SEKV_OPPSETT_NANOPORE != "Nr227/705") %>% 
  # Create columns for sorting
  separate(SEKV_OPPSETT_NANOPORE, into = c("tmp"), sep = "/", remove = FALSE) %>% 
  mutate("tmp" = str_remove(tmp, "Nr")) %>% 
  mutate("tmp" = str_remove(tmp, "A|B")) %>% 
  mutate("tmp" = str_remove(tmp, "/Midnight")) %>% 
  mutate(tmp = as.numeric(tmp)) %>% 
  arrange(desc(tmp)) %>% 
  select(SEKV_OPPSETT_NANOPORE) %>% 
  View("Nanopore")

# Check MIK
BN %>%
  # Convert empty strings to NA
  mutate_all(list(~na_if(.,""))) %>% 
  # Filtrer på coverage >= 94%
  filter(Dekning_Swift >=94) %>%
  # Trekk ut OUS-prøver
  filter(str_detect(SEKV_OPPSETT_SWIFT7, "MIK")) %>%
  # Format the date
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>%
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  filter(is.na(GISAID_EPI_ISL)) %>%
  # Get the various plates
  select(SEKV_OPPSETT_SWIFT7) %>%
  distinct() %>%
  # Create numeric column for sorting
  mutate("tmp" = as.numeric(str_remove(SEKV_OPPSETT_SWIFT7, "MIK"))) %>% 
  arrange(desc(tmp)) %>% 
  select(SEKV_OPPSETT_SWIFT7) %>% 
  View("Swift_MIK")
