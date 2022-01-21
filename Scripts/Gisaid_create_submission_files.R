#!/usr/bin/env Rscript

# Load packages
pacman::p_load(optparse, phylotools, tidyverse, readxl, stringr, lubridate)

# Parse options -----------------------------------------------------------
option_list <- list(
  make_option(c("-p", "--platform"), type="character", default=NULL,
              help="Type platform. Må være enten: Swift_FHI, Swift_MIK, Artic_Illumina, Artic_Nanopore", metavar="character"),
  make_option(c("-o", "--oppsett"), type = "character", default = NULL,
              help = "Navn på plate/oppsett. Må tilsvare mappenavn. F.eks. FHI133 for Swift_FHI, MIK172 for Swift_MIK, 646 for Artic_Illumina eller Nr134A/Nano for Artic_Nanopore", metavar = "character"),
  make_option(c("-s", "--specimen"), type="character", default="Unknown",
              help="Type prøvemateriale. E.g. \"Upper respiratory swab\" [default= %default]", metavar="character"),
  make_option(c("-f", "--fasta"), type="character", default="out_fasta.fasta",
              help="Navn på fasta-fil som skal submittes. F.eks. FHI133.fasta [default= %default]", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default="FHI133.csv",
              help="Navn på metadata-fil som skal submittes. F.eks. FHI133.csv [default= %default]", metavar="character"),
  make_option(c("-S", "--submitter"), type="character", default=NULL,
              help="Gisaid brukernavn til den som submitter. F.eks. jonbra", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$platform)){
  print_help(opt_parser)
  stop("At least one argument must be supplied", call.=FALSE)
}
fasta_filename <- opt$fasta
# Read data from BioNumerics ----------------------------------------------
# Les inn BN spørring. Husk å Refreshe og lagre den originale excel-fila først (N:/Virologi/Influensa/2021/Spørringsfiler BN/SQLSERVER_TestBN_Spørring_Entrytable.xlsx)
#BN <- suppressWarnings(read_excel("/home/docker/N/Influensa/2021/Spørringsfiler BN/SQLSERVER_TestBN_Spørring_Entrytable.xlsx", sheet = "Sporring BN") %>%
#  select(KEY, REKVNR, PROVE_TATT, FYLKENAVN, MATERIALE, PROSENTDEKNING_GENOM, DEKNING_NANOPORE, SEKV_OPPSETT_NANOPORE, DEKNING_NANOPORE, SEKV_OPPSETT_SWIFT7,
#         SEQUENCEID_NANO29, SEQUENCEID_SWIFT, COVERAGE_BREADTH_SWIFT, GISAID_PLATFORM, GISAID_EPI_ISL, GENOTYPE_SVART_I_LABWARE, COVERAGE_BREATH_EKSTERNE,
#         SAMPLE_CATEGORY, INNSENDER, COVERAGE_DEPTH_SWIFT, COVARAGE_DEPTH_NANO, RES_CDC_INFA_RX, RES_CDC_INFB_CT, MELDT_SMITTESPORING) %>%
#  rename("Dekning_Artic" = PROSENTDEKNING_GENOM,
#         "Dekning_Swift" = COVERAGE_BREADTH_SWIFT,
#         "Dekning_Nano" = DEKNING_NANOPORE))
try(load(file = "/mnt/N/Virologi/JonBrate/Prosjekter/BN.RData"))
try(load(file = "/home/docker/N/JonBrate/Prosjekter/BN.RData"))
# BN <- read_excel("/mnt/N/Virologi/Influensa/2021/Spørringsfiler BN/SQLSERVER_TestBN_Spørring_Entrytable.xlsx", sheet = "Sporring BN") %>% select(KEY, REKVNR, PROVE_TATT, FYLKENAVN, MATERIALE, PROSENTDEKNING_GENOM, DEKNING_NANOPORE, SEKV_OPPSETT_NANOPORE, DEKNING_NANOPORE, SEKV_OPPSETT_SWIFT7, SEQUENCEID_NANO29, SEQUENCEID_SWIFT, COVERAGE_BREADTH_SWIFT, GISAID_PLATFORM, GISAID_EPI_ISL, GENOTYPE_SVART_I_LABWARE, COVERAGE_BREATH_EKSTERNE, SAMPLE_CATEGORY, INNSENDER, COVERAGE_DEPTH_SWIFT, COVARAGE_DEPTH_NANO, RES_CDC_INFA_RX, RES_CDC_INFB_CT, MELDT_SMITTESPORING) %>% rename("Dekning_Artic" = PROSENTDEKNING_GENOM, "Dekning_Swift" = COVERAGE_BREADTH_SWIFT, "Dekning_Nano" = DEKNING_NANOPORE)

# Set parameters ----------------------------------------------------------

# Set platform
platform <- opt$platform

# Plate / oppsett
oppsett <- opt$oppsett

# Add fixed columns for metadata
submitter <- opt$submitter
type <- "betacoronavirus"
passage <- "original"
host <- "Human"
gender <- "Unknown"
age <- "Unknown"
status <- "Unknown"
covv_subm_sample_id <- "Unknown"
covv_outbreak <- "Unknown"
covv_add_host_info <- "Unknown"
covv_add_location <- "Unknown"
covv_provider_sample_id <- "Unknown"
covv_last_vaccinated <- "Unknown"
covv_treatment <- "Unknown"
orig_lab <- NA
orig_adr <- NA

# Add specimen
specimen <- opt$specimen


# Set originating labs and adresses ---------------------------------------
lab_lookup_table <- tribble(
  ~`Lab code`, ~`Lab`, ~`Lab address`,
  0,	"Norwegian Institute of Public Health, Department of Virology",	"P.O.Box 222 Skoyen, 0213 Oslo, Norway",
  1,	"Ostfold Hospital Trust - Kalnes, Centre for Laboratory Medicine, Section for gene technology and infection serology", "P.O.Box 300, N-1714 Graalum, Norway",
  2,	"Akershus University Hospital, Department for Microbiology and Infectious Disease Control",	"P.O.Box 1000, N-1478 Loerenskog, Norway",
  3,	"Oslo University Hospital, Department of Microbiology",	"P.O.Box 4956 Nydalen, N-0424 Oslo, Norway",
  4,	"Furst Medical Laboratory",	"Soeren Bulls vei 25, N-1051 Oslo, Norway",
  5,	"Innlandet Hospital Trust, Division Lillehammer, Department for Medical Microbiology", "P.O.Box 990, N-2629 Lillehammer, Norway",
  6,	"Medical Microbiology Unit, Department for Laboratory Medicine, Drammen Hospital, Vestre Viken Health Trust", "P.O.Box 800, N-3004 Drammen, Norway",
  7,	"Vestfold Hospital, Toensberg, Department of Microbiology",	"P.O.Box 2168, N-3103 Toensberg, Norway",
  8,	"Unilabs Laboratory Medicine", "Leirvollen 19, N-3736 Skien, Norway",
  9, NA, NA,
  10,	"Hospital of Southern Norway - Kristiansand, Department of Medical Microbiology",	"P.O.Box 416 Lundsiden, N-4604 Kristiansand, Norway",
  11,	"Dept. of Medical Microbiology, Stavanger University Hospital, Helse Stavanger HF", "P.O.Box 8100, N-4068 Stavanger, Norway",
  12,	"Haukeland University Hospital, Dept. of Microbiology",	"P.O.Box 1400, N-5021 Bergen, Norway",
  13,	"Haugesund Hospital, laboratory for Medical Microbiology", "P.O.Box 2170, N-5504 Haugesund, Norway",
  14,	"Foerde Hospital, Department of Microbiology", "P.O.Box 1000, N-6807 Foerde, Norway",
  15,	"Department of Medical Microbiology - section Molde, Molde Hospital",	"Parkveien 84, N-6407 Molde, Norway",
  16,	"Department of Medical Microbiology, St. Olavs hospital",	"P.O.box 3250 Torgarden, N-7006 Trondheim, Norway",
  17,	"Haugesund Hospital, laboratory for Medical Microbiology", "P.O.Box 2170, N-5504 Haugesund, Norway",
  18,	"Nordland Hospital - Bodo, Laboratory Department, Molecular Biology Unit", "P.O.Box 1480, N-8092 Bodo, Norway",
  19,	"University Hospital of Northern Norway, Department for Microbiology and Infectious Disease Control",	"P.O.Box 56, N-9038 Tromsoe, Norway",
  20, NA, NA,
  21,	NA, NA,
  22,	"Department of medical microbiology, section Aalesund, Aalesund Hospital", "N-6026 Aalesund, Norway",
  23,	NA, NA,
  24, "Department Medical Microbiology, Baerum Hospital, Vestre Viken Health Trust", "P.O.Box 800, N-3004 Drammen, Norway",
  25,	"Telemark Hospital Trust – Skien, Dept. of Medical Microbiology",	"P.O.Box 2900 Kjørbekk, N-3710 Skien",
  26,	"Unilabs Laboratory Medicine", "Silurveien 2 B, N-0380 Oslo, Norway",
  27,	"Oslo Helse", "Hegdehaugsveien 36, 0352 Oslo"
) %>%
  mutate(`Lab code` = as.character(`Lab code`))


# Define lookup function to decide originating lab and andress ------------
lookup_function <- function(metadata) {
  for (row in seq_along(metadata$INNSENDER)) {
    if (metadata[row,]$INNSENDER == 0){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[1,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[1,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 1){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[2,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[2,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 2){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[3,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[3,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 3){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[4,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[4,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 4){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[5,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[5,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 5){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[6,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[6,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 6){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[7,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[7,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 7){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[8,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[8,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 8){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[9,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[9,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 9){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[10,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[10,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 10){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[11,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[11,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 11){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[12,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[12,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 12){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[13,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[13,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 13){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[14,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[14,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 14){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[15,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[15,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 15){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[16,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[16,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 16){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[17,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[17,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 17){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[18,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[18,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 18){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[19,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[19,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 19){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[20,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[20,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 20){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[21,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[21,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 21){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[22,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[22,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 22){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[23,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[23,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 23){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[24,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[24,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 24){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[25,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[25,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 25){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[26,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[26,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 26){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[27,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[27,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 27){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[28,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[28,]$`Lab address`
    }
  }
  return(metadata)
}





# Define filter function --------------------------------------------------
filter_BN <- function(BN) {
  tmp <- BN %>%
    # Remove previously submitted samples
    filter(!is.na(GISAID_EPI_ISL)) %>% 
    filter(GISAID_EPI_ISL == "") %>%
    # Fjerne evt positiv controll
    filter(str_detect(KEY, "pos", negate = TRUE)) %>%
    # Fjerne prøver som mangler Fylkenavn
    filter(!is.na(FYLKENAVN)) %>%
    # Det kan også stå "Ukjent" som Fylkenavn - ta bort
    filter(str_detect(FYLKENAVN, "kjent", negate = TRUE)) %>%
    # Endre Trøndelag til Trondelag
    mutate("FYLKENAVN" = str_replace(FYLKENAVN, "Tr\xf8ndelag", "Trondelag")) %>%
    # Endre Møre og Romsdal
    mutate("FYLKENAVN" = str_replace(FYLKENAVN, "M\xf8re", "More")) %>%
    # Endre Sør
    mutate("FYLKENAVN" = str_replace(FYLKENAVN, "S\xf8r", "Sor")) %>%
    # Fix date format
    mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
    # Drop samples witout collection date
    filter(!is.na(PROVE_TATT))

    if (platform == "Artic_Illumina") {
    oppsett_details <- tmp %>%
      filter(str_detect(SAMPLE_CATEGORY, oppsett)) %>%
      # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
      filter(!is.na(MELDT_SMITTESPORING)) %>%
      # Filtrer på coverage >= 94%
      filter(Dekning_Artic >=94) %>%
      mutate(SEARCH_COLUMN = RES_CDC_INFB_CT) %>%
      rename("COVERAGE" = RES_CDC_INFA_RX)
  } else if (platform == "Artic_Nanopore") {
    oppsett_details <- tmp %>%
      filter(str_detect(SEKV_OPPSETT_NANOPORE, oppsett)) %>%
      # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
      filter(!is.na(MELDT_SMITTESPORING)) %>%
      # Filtrer på coverage >= 94%
      filter(Dekning_Nano >=94) %>%
      mutate(SEARCH_COLUMN = SEQUENCEID_NANO29) %>%
      rename("COVERAGE" = COVARAGE_DEPTH_NANO)
  } else if (platform == "Swift_FHI") {
    oppsett_details <- tmp %>%
      filter(str_detect(SEKV_OPPSETT_SWIFT7, oppsett)) %>% 
      # Behold bare de som er meldt smittesporing. Disse skal da være godkjent.
      filter(!is.na(MELDT_SMITTESPORING)) %>%
      # Filtrer på coverage >= 94%
      filter(Dekning_Swift >=94) %>%
      # Create column for looping through later
      mutate(SEARCH_COLUMN = SEQUENCEID_SWIFT) %>%
      rename("COVERAGE" = COVERAGE_DEPTH_SWIFT)
  } else if (platform == "Swift_MIK") {
    oppsett_details <- tmp %>%
      filter(str_detect(SEKV_OPPSETT_SWIFT7, oppsett)) %>% 
      # Filtrer på coverage >= 94%
      filter(Dekning_Swift >=94) %>%
      # Remove "OUS-" from Sequence ID
      mutate(SEQUENCE_ID_TRIMMED = str_remove(SEQUENCEID_SWIFT, "OUS-")) %>%
      # Create column for looping through later
      mutate(SEARCH_COLUMN = SEQUENCEID_SWIFT) %>%
      rename("COVERAGE" = COVERAGE_DEPTH_SWIFT)
  }
    return(oppsett_details)
}

# Find sequences on N: and create fasta object ----------------------------
find_sequences <- function(platform, oppsett) {
  if (platform == "Swift_FHI"){
    # Search the N: disk for consensus sequences
    try(dirs_fhi <- c(list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2021/", recursive = FALSE), 
                      list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2022/", recursive = FALSE)))
    try(dirs_fhi <- c(list.dirs("/home/docker/N/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2021/", recursive = FALSE),
                      list.dirs("/home/docker/N/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2022/", recursive = FALSE)))

    # Pick our the relevant oppsett
    dir <- dirs_fhi[grep(paste0(oppsett, "\\b"), dirs_fhi)]

    # List the files
    filepaths <- list.files(path = dir,
                            pattern = "ivar\\.consensus\\.masked_Nremoved\\.fa$",
                            full.names = TRUE,
                            recursive = TRUE)
    samples <- gsub("_.*", "", basename(filepaths))

  } else if (platform == "Swift_MIK") {
    # Search the N: disk for consensus sequences.
    try(dirs_fhi <- list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK", recursive = FALSE))
    try(dirs_fhi <- list.dirs("/home/docker/N/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK", recursive = FALSE))
    # Pick our the relevant oppsett
    dir <- dirs_fhi[grep(paste0(oppsett, "\\b"), dirs_fhi)]

    # List the files
    filepaths <- list.files(path = dir,
                            pattern = "ivar\\.consensus\\.masked_Nremoved\\.fa$",
                            full.names = TRUE,
                            recursive = TRUE)

    samples <- gsub("_.*","", gsub(".*/","", filepaths))
  } else if (platform == "Artic_Illumina") {
    # Search the N: disk for consensus sequences.
    try(dirs_fhi <- c(list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina/2021", recursive = FALSE), 
                      list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina/2022", recursive = FALSE)))
    try(dirs_fhi <- c(list.dirs("/home/docker/N/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina/2021", recursive = FALSE), 
                      list.dirs("/home/docker/N/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina/2022", recursive = FALSE)))
    
    # Pick our the relevant oppsett
    dir <- dirs_fhi[grep(oppsett, dirs_fhi)]

    # List the files
    filepaths <- list.files(path = dir,
                            pattern = "consensus\\.fa$",
                            full.names = TRUE,
                            recursive = TRUE)

    # Dropper det siste tallet.
    samples <- gsub("_.*", "", basename(filepaths))
    #samples <- str_sub(gsub("Artic", "", gsub("_.*","", gsub(".*/","", filepaths))), start = 1, end = -2)
  } else if (platform == "Artic_Nanopore") {
    # Search the N: disk for consensus sequences.
    try(dirs_fhi <- c(list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/2021", recursive = FALSE),
                      list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/2022", recursive = FALSE)))
    try(dirs_fhi <- c(list.dirs("/home/docker/N/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/2021", recursive = FALSE),
                      list.dirs("/home/docker/N/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/2022", recursive = FALSE)))
    
    # Pick our the relevant oppsett
    oppsett <- gsub("Nr", "", (gsub("/Nano", "", oppsett)))
    dir <- dirs_fhi[grep(paste0(oppsett), dirs_fhi)]
    
    # List the files
    filepaths <- list.files(path = dir,
                            pattern = "consensus\\.fasta$",
                            full.names = TRUE,
                            recursive = TRUE)
  }

  # Find which filepaths to keep
  keep <- vector("character", length = length(oppsett_details$SEARCH_COLUMN))
  for (i in seq_along(oppsett_details$SEARCH_COLUMN)){
    keep[i] <- filepaths[grep(oppsett_details$SEARCH_COLUMN[i], filepaths)]
  }

  # Read each fasta file and combine them to create one file
  # First create empty data frame to fill
  fastas <- data.frame(seq.name = character(),
                       seq.text = character())

  for (f in seq_along(keep)){
    tmp <- read.fasta(keep[f])      # read the file
    fastas <- rbind(fastas, tmp)    # append the current file
  }

  # Convert to tibble for easier manipulation
  fastas <- as_tibble(fastas)

  # Fix names to match KEY
  if (platform == "Swift_FHI") {
    # Fix names to match SEQUENCEID_SWIFT
    fastas <- fastas %>%
      mutate(SEQUENCEID_SWIFT = str_remove(seq.name, "_ivar_masked"))
      #mutate(tmp = str_remove(seq.name, "_ivar_masked")) %>%
      #mutate(tmp = str_remove(tmp, "B$")) %>%
      #mutate(KEY = str_remove(tmp, "SWIFT"))
  } else if (platform == "Swift_MIK") {
    # Fix names to match SEQUENCEID_SWIFT
    fastas <- fastas %>%
      mutate("tmp" = str_remove(seq.name, "_ivar_masked")) %>%
      # Some of the MIK-samples have a leading number in front of OUS
      mutate(SEQUENCE_ID_TRIMMED = gsub(".*OUS-", "", .$tmp))
  } else if (platform == "Artic_Illumina") {
    fastas <- fastas %>%
      rename("RES_CDC_INFB_CT" = `seq.name`)
      #mutate(tmp = str_remove(seq.name, "Artic")) %>%
      #mutate(KEY = str_sub(tmp, start = 1, end = -2))
  } else if (platform == "Artic_Nanopore") {
    # Fix names to match SEQUENCEID_NANO29
    fastas <- fastas %>%
      separate("seq.name", into = c("SEQUENCEID_NANO29", NA, NA), sep = "/", remove = F)

  }

  # Sett Virus name som fasta header
  # Først lage en mapping mellom KEY og virus name
  if (platform == "Swift_FHI") {
    SEQUENCEID_virus_mapping <- oppsett_details %>%
      # Trenger også å lage Virus name
      # Lage kolonne for "year"
      separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
      # Trekke ut sifrene fra 5 og til det siste fra BN KEY
      mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>%
      # Fjerne ledende nuller fra stammenavnet
      mutate("Uniq_nr" = str_remove(Uniq_nr, "^0+")) %>% 
      # Legge til kolonner med fast informasjon for å lage "Virus name" senere
      add_column("Separator" = "/",
                 "GISAID_prefix" = "hCoV-19/",
                 "Country" = "Norway/",
                 "Continent" = "Europe/") %>%
      # Make "Virus name" column
      unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
      select(KEY, SEQUENCEID_SWIFT, covv_virus_name)
    
    fastas <- left_join(fastas, SEQUENCEID_virus_mapping, by = "SEQUENCEID_SWIFT") %>%
      select(`seq.name` = covv_virus_name,
             seq.text)
    
  } else if (platform == "Swift_MIK") {
    KEY_virus_mapping <- oppsett_details %>%
      # Lage kolonne for "year"
      separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
      # Trekke ut sifrene fra 5 og til det siste fra BN KEY
      mutate("Uniq_nr" = str_sub(KEY, start = 1, end = -1)) %>%
      # Legge til kolonner med fast informasjon for å lage "Virus name" senere
      add_column("Separator" = "/",
                 "GISAID_prefix" = "hCoV-19/",
                 "Country" = "Norway/",
                 "Continent" = "Europe/") %>%
      # Make "Virus name" column
      unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
      select(SEQUENCEID_SWIFT, KEY, covv_virus_name, SEQUENCE_ID_TRIMMED)


    fastas <- left_join(fastas, KEY_virus_mapping, by = "SEQUENCE_ID_TRIMMED") %>%
      select(`seq.name` = covv_virus_name,
             seq.text)
  } else if (platform == "Artic_Illumina") {
    SEQUENCEID_virus_mapping <- oppsett_details %>%
      # Trenger også å lage Virus name
      # Lage kolonne for "year"
      separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
      # Trekke ut sifrene fra 5 og til det siste fra BN KEY
      mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>%
      # Fjerne ledende nuller fra stammenavnet
      mutate("Uniq_nr" = str_remove(Uniq_nr, "^0+")) %>% 
      # Legge til kolonner med fast informasjon for å lage "Virus name" senere
      add_column("Separator" = "/",
                 "GISAID_prefix" = "hCoV-19/",
                 "Country" = "Norway/",
                 "Continent" = "Europe/") %>%
      # Make "Virus name" column
      unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
      select(KEY, RES_CDC_INFB_CT, covv_virus_name)
    
    fastas <- left_join(fastas, SEQUENCEID_virus_mapping, by = "RES_CDC_INFB_CT") %>%
      select(`seq.name` = covv_virus_name,
             seq.text)

  } else if (platform == "Artic_Nanopore") {
    SEQUENCEID_virus_mapping <- oppsett_details %>%
      # Trenger også å lage Virus name
      # Lage kolonne for "year"
      separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
      # Trekke ut sifrene fra 5 og til det siste fra BN KEY
      mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>%
      # Fjerne ledende nuller fra stammenavnet
      mutate("Uniq_nr" = str_remove(Uniq_nr, "^0+")) %>% 
      # Legge til kolonner med fast informasjon for å lage "Virus name" senere
      add_column("Separator" = "/",
                 "GISAID_prefix" = "hCoV-19/",
                 "Country" = "Norway/",
                 "Continent" = "Europe/") %>%
      # Make "Virus name" column
      unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
      select(KEY, SEQUENCEID_NANO29, covv_virus_name)

    fastas <- left_join(fastas, SEQUENCEID_virus_mapping, by = "SEQUENCEID_NANO29") %>%
      select(`seq.name` = covv_virus_name,
             seq.text)
  }
  return(fastas)
}


# Define metadata function ------------------------------------------------
create_metadata <- function(oppsett_details) {

  metadata <- oppsett_details %>%
    # Lage kolonne for "year"
    separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE)

  if (platform == "Swift_MIK") {
    metadata <- metadata %>%
      # For OUS bruke hele KEY som stammenr
      mutate("Uniq_nr" = str_sub(KEY, start = 1, end = -1))
  } else {
    metadata <- metadata %>%
      # Trekke ut sifrene fra 5 og til det siste fra BN KEY
      mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>% 
      # Fjerne ledende nuller fra stammenavnet
      mutate("Uniq_nr" = str_remove(Uniq_nr, "^0+"))
  }

    metadata <- metadata %>%
      # Legge til kolonner med fast informasjon for å lage "Virus name" senere
      add_column("Separator" = "/",
                 "GISAID_prefix" = "hCoV-19/",
                 "Country" = "Norway/",
                 "Continent" = "Europe/") %>%
      # Make "Virus name" column
      unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>%
      # Lage Location-kolonne
      unite("covv_location", c(Continent, Country, FYLKENAVN), sep = "", remove = FALSE) %>%
      # Legge til faste kolonner
      add_column("submitter" = submitter,
                 "fn" = fasta_filename,
                 "covv_type" = type,
                 "covv_passage" = passage,
                 "covv_host" = host,
                 "covv_gender" = gender,
                 "covv_patient_age" = age,
                 "covv_patient_status" = status,
                 "covv_specimen" = specimen,
                 "covv_seq_technology" = seq_tech,
                 "covv_assembly_method" = ass_method,
                 "covv_orig_lab" = orig_lab,
                 "covv_orig_lab_addr" = orig_adr,
                 "covv_subm_lab" = sub_lab,
                 "covv_subm_lab_addr" = address,
                 "covv_authors" = authors,
                 "covv_subm_sample_id" = covv_subm_sample_id,
                 "covv_outbreak" = covv_outbreak,
                 "covv_add_host_info" = covv_add_host_info,
                 "covv_add_location" = covv_add_location,
                 "covv_provider_sample_id" = covv_provider_sample_id,
                 "covv_last_vaccinated" = covv_last_vaccinated,
                 "covv_treatment" = covv_treatment) %>%
      # Beholde endelige kolonner og rekkefølge
      select("submitter",
             "fn",
             "covv_virus_name",
             "covv_type",
             "covv_passage",
             "covv_collection_date" = PROVE_TATT,
             "covv_location",
             "covv_host",
             "covv_gender",
             "covv_patient_age",
             "covv_patient_status",
             "covv_specimen",
             "covv_seq_technology",
             "covv_assembly_method",
             "covv_orig_lab",
             "covv_orig_lab_addr",
             "covv_subm_lab",
             "covv_subm_lab_addr",
             "covv_authors",
             "covv_subm_sample_id",
             "covv_outbreak",
             "covv_add_host_info",
             "covv_add_location",
             "covv_provider_sample_id",
             "covv_last_vaccinated",
             "covv_treatment",
             "covv_coverage" = COVERAGE,
             "INNSENDER")

  if (platform == "Swift_MIK") {
    # Remove column INNSENDER
    metadata <- metadata %>% select(-INNSENDER)
  } else {
    # Legge inn orig lab og adresse
    metadata <- lookup_function(metadata)

    # Remove column INNSENDER
    metadata <- metadata %>% select(-INNSENDER)
  }
  return(metadata)
}
# Define Frameshift function ----------------------------------------------
FS <- function(fastas){
  #### Run Frameshift analysis ####
  # write temporary fasta file to Frameshift folder
  suppressMessages(try(setwd("/home/jonr/tmp_gisaid/Frameshift/")))
  suppressMessages(try(setwd("/home/docker/Fastq/Frameshift")))
  # dat2fasta(fastas, outfile = "/home/jonr/FHI_SC2_Pipeline_Illumina/Frameshift/tmp.fasta")

  dat2fasta(fastas, outfile = "tmp.fasta")

  # Run the frameshift script
  system("Rscript --vanilla /home/docker/Scripts/CSAK_Frameshift_Finder_docker.R")
  # system("docker run --rm -v $(pwd):/home/docker/Fastq garcianacho/fhisc2:Illumina Rscript --vanilla /home/docker/Scripts/CSAK_Frameshift_Finder_docker.R")

}

# Remove bad FS from fasta ------------------------------------------------
remove_FS_fasta <- function(fastas){
  #### Remove any samples with bad FS ####
  suppressMessages(try(setwd("/home/jonr/tmp_gisaid/Frameshift/")))
  suppressMessages(try(setwd("/home/docker/Fastq/Frameshift")))
  FS_OK <- read_excel("FrameShift_tmp.xlsx") %>%
    filter(Ready == "YES") %>%
    rename(`seq.name` = "Sample")

  # Rename navn til å matche navn i fastas
  # Join fastas with FS to keep
  fastas_clean <- left_join(FS_OK, fastas, by = "seq.name") %>%
    select(`seq.name`, `seq.text`)
  return(fastas_clean)
}

# Remove bad FS from metadata ---------------------------------------------
remove_FS_metadata <- function(metadata){
  #### Drop the same samples from the metadata file #####
  suppressMessages(try(setwd("/home/jonr/tmp_gisaid/Frameshift/")))
  suppressMessages(try(setwd("/home/docker/Fastq/Frameshift")))
  # Define samples to keep (i.e. with OK FS)
  FS_OK <- read_excel("FrameShift_tmp.xlsx") %>%
    filter(Ready == "YES") %>%
    rename("covv_virus_name" = "Sample")

  metadata_clean <- left_join(FS_OK, metadata, by = "covv_virus_name") %>%
    select(-Deletions, -Frameshift, -Insertions, -Ready, -Comments)

  return(metadata_clean)
}


# Define clean up and write function --------------------------------------
clean_up_and_write <- function(fastas_clean, metadata_clean) {

  suppressMessages(try(setwd("/home/jonr/tmp_gisaid/")))
  suppressMessages(try(setwd("/home/docker/Fastq/")))

  if (platform == "Artic_Nanopore") {
    oppsett_stripped <- gsub("/", "", oppsett)
    file.remove(dir("Frameshift/", pattern = "csv|fasta", full.names = T))
    file.rename("Frameshift/FrameShift_tmp.xlsx", paste0("FrameShift_", oppsett_stripped, ".xlsx"))
    write_csv(metadata_clean, file = paste0(oppsett_stripped, ".csv"))
    dat2fasta(fastas_clean, outfile = paste0(oppsett_stripped, ".fasta"))
  } else {
    file.remove(dir("Frameshift/", pattern = "csv|fasta", full.names = T))
    file.rename("Frameshift/FrameShift_tmp.xlsx", paste0("FrameShift_", oppsett, ".xlsx"))
    write_csv(metadata_clean, file = opt$metadata)
    dat2fasta(fastas_clean, outfile = fasta_filename)
  }
}

# Start script ------------------------------------------------------------

if (platform == "Swift_FHI"){
  print ("Platform is Swift FHI")

  # Add platform-specific columns.
  seq_tech <- "Illumina Swift Amplicon SARS-CoV-2 protocol at Norwegian Sequencing Centre"
  ass_method <- "Assembly by reference based mapping using Bowtie2 with iVar majority rules consensus"
  sub_lab <- "Norwegian Institute of Public Health, Department of Virology"
  address <- "P.O.Box 222 Skoyen, 0213 Oslo, Norway"
  authors <- "Kathrine Stene-Johansen, Kamilla Heddeland Instefjord, Hilde Elshaug, Garcia Llorente Ignacio, Jon Bråte, Engebretsen Serina Beate, Pedersen Benedikte Nevjen, Line Victoria Moen, Hilde Nordby Falkenhaug, Debech Nadia, Atiya R Ali, Marie Paulsen Madsen, Rasmus Riis Kopperud, Hilde Vollan, Karoline Bragstad, Olav Hungnes"

  #### Trekke ut prøver ####
  oppsett_details <- filter_BN(BN)

  #### Lage metadata ####
  metadata <- create_metadata(oppsett_details)

  #### Find sequences on N: ####
  fastas <- find_sequences(platform, oppsett)

  #### Run Frameshift analysis ####
  FS(fastas)
  fastas_clean <- remove_FS_fasta(fastas)
  metadata_clean <- remove_FS_metadata(metadata)

  #### Clean up and write files
  clean_up_and_write(fastas = fastas_clean, metadata = metadata_clean)

} else if (platform == "Swift_MIK") {
  print ("Platform is Swift MIK")

  # Add platform-specific columns.
  seq_tech <- "Illumina Swift Amplicon SARS-CoV-2 protocol at Norwegian Sequencing Centre"
  ass_method <- "Assembly by reference based mapping using Bowtie2 with iVar majority rules consensus"
  sub_lab <- "Norwegian Institute of Public Health, Department of Virology"
  address <- "P.O.Box 222 Skoyen, 0213 Oslo, Norway"
  authors <- "Mona Holberg-Petersen, Lise Andresen, Cathrine Fladeby, Mariann Nilsen, Teodora Plamenova Ribarska, Pål Marius Bjørnstad, Gregor D. Gilfillan, Arvind Yegambaram Meenakshi Sundaram, Kathrine Stene-Johansen, Kamilla Heddeland Instefjord, Hilde Elshaug, Garcia Llorente Ignacio, Jon Bråte, Pedersen Benedikte Nevjen, Line Victoria Moen, Rasmus Riis Kopperud, Hilde Vollan, Olav Hungnes, Karoline Bragstad"
  orig_lab <- "Oslo University Hospital, Department of Microbiology"
  orig_adr <- "P.O.Box 4956 Nydalen, N-0424 Oslo, Norway"

  #### Trekke ut prøver ####
  oppsett_details <- filter_BN(BN)

  #### Lage metadata ####
  metadata <- create_metadata(oppsett_details)

  #### Find sequences on N: ####
  fastas <- find_sequences(platform, oppsett)

  #### Run Frameshift analysis ####
  FS(fastas)
  fastas_clean <- remove_FS_fasta(fastas)
  metadata_clean <- remove_FS_metadata(metadata)

  #### Clean up and write files
  clean_up_and_write(fastas = fastas_clean, metadata = metadata_clean)

} else if (platform == "Artic_Illumina") {
  print ("Platform is Artic Illumina")

  # Add platform-specific columns.
  seq_tech <- "Illumina MiSeq, modified ARTIC protocol with V4.1 primers"
  ass_method <- "Assembly by reference based mapping using Tanoti with iVar majority rules consensus"
  sub_lab <- "Norwegian Institute of Public Health, Department of Virology"
  address <- "P.O.Box 222 Skoyen, 0213 Oslo, Norway"
  authors <- "Kathrine Stene-Johansen, Kamilla Heddeland Instefjord, Hilde Elshaug, Garcia Llorente Ignacio, Jon Bråte, Engebretsen Serina Beate, Pedersen Benedikte Nevjen, Line Victoria Moen, Hilde Nordby Falkenhaug, Debech Nadia, Atiya R Ali, Marie Paulsen Madsen, Rasmus Riis Kopperud, Hilde Vollan, Karoline Bragstad, Olav Hungnes"

  ##### Trekke ut prøver ####
  oppsett_details <- filter_BN(BN)

  #### Lage metadata ####
  metadata <- create_metadata(oppsett_details)

  #### Find sequences on N: ####
  fastas <- find_sequences(platform, oppsett)

  #### Run Frameshift analysis ####
  FS(fastas)
  fastas_clean <- remove_FS_fasta(fastas)
  metadata_clean <- remove_FS_metadata(metadata)

  #### Clean up and write files
  clean_up_and_write(fastas = fastas_clean, metadata = metadata_clean)

} else if (platform == "Artic_Nanopore") {
  print ("Platform is Artic Nanopore")

  # Add platform-specific columns.
  seq_tech <- "Nanopore GridIon, Artic V4.1 protocol modified"
  ass_method <- "Assembly by reference based mapping using the Artic Nanopore protocol with medaka"
  sub_lab <- "Norwegian Institute of Public Health, Department of Virology"
  address <- "P.O.Box 222 Skoyen, 0213 Oslo, Norway"
  authors <- "Kathrine Stene-Johansen, Kamilla Heddeland Instefjord, Hilde Elshaug, Garcia Llorente Ignacio, Jon Bråte, Engebretsen Serina Beate, Pedersen Benedikte Nevjen, Line Victoria Moen, Hilde Nordby Falkenhaug, Debech Nadia, Atiya R Ali, Marie Paulsen Madsen, Rasmus Riis Kopperud, Hilde Vollan, Karoline Bragstad, Olav Hungnes"

  ##### Trekke ut prøver ####
  oppsett_details <- filter_BN(BN)

  #### Lage metadata ####
  metadata <- create_metadata(oppsett_details)

  #### Find sequences on N: ####
  fastas <- find_sequences(platform, oppsett)

  #### Run Frameshift analysis ####
  FS(fastas)
  fastas_clean <- remove_FS_fasta(fastas)
  metadata_clean <- remove_FS_metadata(metadata)

  #### Clean up and write files
  clean_up_and_write(fastas = fastas_clean, metadata = metadata_clean)
}
