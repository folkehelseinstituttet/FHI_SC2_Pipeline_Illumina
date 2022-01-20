# Load packages
pacman::p_load(tidyverse, readxl, stringr, lubridate)

# Read data from BioNumerics ----------------------------------------------
try(load(file = "/mnt/N/Virologi/JonBrate/Prosjekter/BN.RData"))
try(load(file = "N:/Virologi/JonBrate/Prosjekter/BN.RData"))
# Les inn BN spørring. Husk å Refreshe og lagre den originale excel-fila først (N:/Virologi/Influensa/2021/Spørringsfiler BN/SQLSERVER_TestBN_Spørring_Entrytable.xlsx)
#BN <- read_excel("/mnt/N/Virologi/Influensa/2021/Spørringsfiler BN/SQLSERVER_TestBN_Spørring_Entrytable.xlsx", sheet = "Sporring BN") %>%
#  select(KEY, REKVNR, PROVE_TATT, FYLKENAVN, MATERIALE, PROSENTDEKNING_GENOM, DEKNING_NANOPORE, SEKV_OPPSETT_NANOPORE, DEKNING_NANOPORE,
#         SEKV_OPPSETT_SWIFT7, SEQUENCEID_NANO29, SEQUENCEID_SWIFT, COVERAGE_BREADTH_SWIFT, GISAID_PLATFORM, GISAID_EPI_ISL,
#         GENOTYPE_SVART_I_LABWARE, COVERAGE_BREATH_EKSTERNE, SAMPLE_CATEGORY, INNSENDER, COVERAGE_DEPTH_SWIFT, COVARAGE_DEPTH_NANO,
#         RES_CDC_INFA_RX, RES_CDC_INFB_CT, MELDT_SMITTESPORING) %>%
#  rename("Dekning_Artic" = PROSENTDEKNING_GENOM, "Dekning_Swift" = COVERAGE_BREADTH_SWIFT, "Dekning_Nano" = DEKNING_NANOPORE)

# Check Swift FHI
BN %>%
  # Filtrer på coverage >= 94%
  filter(Dekning_Swift >=94) %>%
  # Format the date
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>%
  # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
  filter(!is.na(MELDT_SMITTESPORING)) %>%
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  # Using the is.na() filter because there could be other strings than the EPI_ISL accession written
  filter(GISAID_EPI_ISL == "") %>%
  #filter(is.na(GISAID_EPI_ISL)) %>% pull(GISAI
  # Get the various platesD_EPI_ISL)
  select(KEY, PROVE_TATT, SEKV_OPPSETT_SWIFT7) %>%
  filter(str_detect(SEKV_OPPSETT_SWIFT7, "FHI")) %>%
  # Select one sample per oppsett
  group_by(SEKV_OPPSETT_SWIFT7) %>%
  slice_head(n = 1) %>%
  # Sort by date
  arrange(desc(PROVE_TATT)) %>% View()

# Check Artic Illumina
BN %>%
  # Filtrer på coverage >= 94%
  filter(Dekning_Artic >=94) %>%
  # Format the date
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>%
  # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
  filter(!is.na(MELDT_SMITTESPORING)) %>%
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  # Using the is.na() filter because there could be other strings than the EPI_ISL accession written
  filter(GISAID_EPI_ISL == "") %>%
  #filter(is.na(GISAID_EPI_ISL)) %>%
  # Get the various plates
  select(KEY, PROVE_TATT, SAMPLE_CATEGORY) %>%
  filter(str_detect(SAMPLE_CATEGORY, "naekstrakt", negate = TRUE)) %>%
  # Select one sample per oppsett
  group_by(SAMPLE_CATEGORY) %>%
  slice_head(n = 1) %>%
  # Sort by date
  arrange(desc(PROVE_TATT)) %>% View()

# Check Nanopore
BN %>%
  # Filtrer på coverage >= 94%
  filter(Dekning_Nano >=94) %>%
  # Format the date
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>%
  # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
  filter(!is.na(MELDT_SMITTESPORING)) %>%
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  # Using the is.na() filter because there could be other strings than the EPI_ISL accession written
  filter(GISAID_EPI_ISL == "") %>%
  #filter(is.na(GISAID_EPI_ISL)) %>%
  # Get the various plates
  select(KEY, PROVE_TATT, SEKV_OPPSETT_NANOPORE) %>%
  # Select one sample per oppsett
  group_by(SEKV_OPPSETT_NANOPORE) %>%
  slice_head(n = 1) %>%
  # Sort by date
  arrange(desc(PROVE_TATT)) %>% View()

# Check MIK
BN %>%
  # Filtrer på coverage >= 94%
  filter(Dekning_Swift >=94) %>%
  filter(str_detect(SEKV_OPPSETT_SWIFT7, "MIK")) %>%
  # Format the date
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>%
  # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
  # filter(!is.na(MELDT_SMITTESPORING)) %>%
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  # Using the is.na() filter because there could be other strings than the EPI_ISL accession written
  filter(GISAID_EPI_ISL == "") %>%
  #filter(is.na(GISAID_EPI_ISL)) %>%
  # Get the various plates
  select(KEY, PROVE_TATT, SEKV_OPPSETT_SWIFT7) %>%
  # Select one sample per oppsett
  group_by(SEKV_OPPSETT_SWIFT7) %>%
  slice_head(n = 1) %>%
  # Sort by date
  arrange(desc(PROVE_TATT)) %>% View()
