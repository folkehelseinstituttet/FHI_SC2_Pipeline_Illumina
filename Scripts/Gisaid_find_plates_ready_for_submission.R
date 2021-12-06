# Load packages
pacman::p_load(tidyverse, readxl, stringr, lubridate)

# Read data from BioNumerics ----------------------------------------------
# Les inn BN spørring. Husk å Refreshe og lagre den originale excel-fila først (N:/Virologi/Influensa/2021/Spørringsfiler BN/SQLSERVER_TestBN_Spørring_Entrytable.xlsx)
BN <- read_excel("/mnt/N/Virologi/Influensa/2021/Spørringsfiler BN/SQLSERVER_TestBN_Spørring_Entrytable.xlsx", sheet = "Sporring BN") %>% 
  select(KEY, REKVNR, PROVE_TATT, FYLKENAVN, MATERIALE, PROSENTDEKNING_GENOM, DEKNING_NANOPORE, SEKV_OPPSETT_NANOPORE, DEKNING_NANOPORE, 
         SEKV_OPPSETT_SWIFT7, SEQUENCEID_NANO29, SEQUENCEID_SWIFT, COVERAGE_BREADTH_SWIFT, GISAID_PLATFORM, GISAID_EPI_ISL, 
         GENOTYPE_SVART_I_LABWARE, COVERAGE_BREATH_EKSTERNE, SAMPLE_CATEGORY, INNSENDER, COVERAGE_DEPTH_SWIFT, COVARAGE_DEPTH_NANO, 
         RES_CDC_INFA_RX, RES_CDC_INFB_CT, MELDT_SMITTESPORING) %>% 
  rename("Dekning_Artic" = PROSENTDEKNING_GENOM, "Dekning_Swift" = COVERAGE_BREADTH_SWIFT, "Dekning_Nano" = DEKNING_NANOPORE)

# Check everything except MIK
BN %>%
  # Format the date 
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
  # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
  filter(!is.na(MELDT_SMITTESPORING)) %>%  
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  # Using the is.na() filter because there could be other strings than the EPI_ISL accession written
  filter(is.na(GISAID_EPI_ISL)) %>% 
  # Get the various plates
  select(KEY, PROVE_TATT, SEKV_OPPSETT_NANOPORE, SEKV_OPPSETT_SWIFT7, SAMPLE_CATEGORY) %>%
  # Select one sample per oppsett
  group_by(SEKV_OPPSETT_NANOPORE, SEKV_OPPSETT_SWIFT7, SAMPLE_CATEGORY) %>% 
  slice_head(n = 1) %>% 
  # Sort by date
  arrange(desc(PROVE_TATT)) %>% View()

# Check MIK
BN %>%
  filter(str_detect(SEKV_OPPSETT_SWIFT7, "MIK")) %>% 
  # Format the date 
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  # Using the is.na() filter because there could be other strings than the EPI_ISL accession written
  filter(is.na(GISAID_EPI_ISL)) %>% 
  # Get the various plates
  select(KEY, PROVE_TATT, SEKV_OPPSETT_SWIFT7) %>%
  # Select one sample per oppsett
  group_by(SEKV_OPPSETT_SWIFT7) %>% 
  slice_head(n = 1) %>% 
  # Sort by date
  arrange(desc(PROVE_TATT)) %>% View()

# Check already submitted MIK samples
BN %>%
  filter(str_detect(SEKV_OPPSETT_SWIFT7, "MIK")) %>% 
  # Format the date 
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  filter(str_detect(GISAID_EPI_ISL, "^EPI")) %>% 
  select(KEY, PROVE_TATT, GISAID_EPI_ISL, SEKV_OPPSETT_SWIFT7) %>% 
  arrange(desc(PROVE_TATT)) %>% View()
