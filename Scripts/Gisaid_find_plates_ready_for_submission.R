# Load packages
pacman::p_load(tidyverse, readxl, stringr, lubridate)

# Read data from BioNumerics ----------------------------------------------
load(file = "/mnt/N/Virologi/JonBrate/Prosjekter/BN.RData")
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

# Check Omikron
tmp2 <- BN %>%
  filter(str_detect(PANGOLIN_NOM, "^BA")) %>%
  # Format the date
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>%
  # Fjerne evt positive kontroller
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fjerne de som allerede er submittet til Gisaid. NB - husk å importere submisjonsresultater først.
  # Using the is.na() filter because there could be other strings than the EPI_ISL accession written
  filter(is.na(GISAID_EPI_ISL)) %>%
  filter(Dekning_Swift >=97)
  filter(Dekning_Nano >=97) %>%
  filter(str_detect(SEKV_OPPSETT_NANOPORE, "198", negate = TRUE))

  filter(grepl(NA, Dekning_Artic))
  !grepl("str", NA) returns TRUE, so is kept.


  # Filtrer på coverage >= 97%
  filter(Dekning_Artic >97) %>%
  filter(Dekning_Nano >=97) %>%
  filter(Dekning_Swift >=97) %>%
  # Get the various plates
  select(KEY, PROVE_TATT, SEKV_OPPSETT_NANOPORE, SEKV_OPPSETT_SWIFT7, SAMPLE_CATEGORY) %>%
  # Select one sample per oppsett
  group_by(SEKV_OPPSETT_NANOPORE, SEKV_OPPSETT_SWIFT7, SAMPLE_CATEGORY) %>%
  slice_head(n = 1) %>%
  # Sort by date
  arrange(desc(PROVE_TATT)) %>% View()

if (platform == "Artic_Illumina") {
  oppsett_details <- tmp %>%
    filter(str_detect(SAMPLE_CATEGORY, oppsett)) %>%

    mutate(SEARCH_COLUMN = KEY) %>%
    rename("COVERAGE" = RES_CDC_INFA_RX)
} else if (platform == "Artic_Nanopore") {
  oppsett_details <- tmp %>% 
    filter(str_detect(SEKV_OPPSETT_NANOPORE, oppsett)) %>%
    # Filtrer på coverage >= 97%

    mutate(SEARCH_COLUMN = KEY) %>%
    rename("COVERAGE" = COVARAGE_DEPTH_NANO)
} else if (platform == "Swift_FHI") {
  oppsett_details <- tmp %>%
    filter(SEKV_OPPSETT_SWIFT7 == oppsett) %>%
    # Filtrer på coverage >= 97%

    # Create column for looping through later
    mutate(SEARCH_COLUMN = KEY) %>%
    rename("COVERAGE" = COVERAGE_DEPTH_SWIFT)
} else if (platform == "Swift_MIK") {
  oppsett_details <- tmp %>%
    filter(SEKV_OPPSETT_SWIFT7 == oppsett) %>%
    # Filtrer på coverage >= 97%
    filter(Dekning_Swift >=97) %>%
    # Remove "OUS-" from Sequence ID
    mutate(SEQUENCE_ID_TRIMMED = str_remove(SEQUENCEID_SWIFT, "OUS-")) %>%
    # Create column for looping through later
    mutate(SEARCH_COLUMN = SEQUENCE_ID_TRIMMED) %>%
    rename("COVERAGE" = COVERAGE_DEPTH_SWIFT)
