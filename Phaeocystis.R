library(maps) # Draw Geographical Maps
library(rgbif) # Interface to the Global Biodiversity Information Facility API
library(robis) # Ocean Biodiversity Information System (OBIS) Client
library(tidyverse) # Easily Install and Load the 'Tidyverse'
library(worrms) # World Register of Marine Species (WoRMS) Client
library(dggridR) # Discrete Global Grids
library(ggthemes) # Extra Themes, Scales and Geoms for 'ggplot2'
library(RColorBrewer) # ColorBrewer Palettes
library(pangaear) # Client for the 'Pangaea' Database
library(actverse) # Tools for Actigraphy Data Analysis

# Define the taxonomic name -----------------------------------------------

taxon_name <- c(
  "Phaeocystis globosa",
  "Phaeocystis pouchetii",
  "Phaeocystis antarctica"
)

# Search for occurrences of the taxon using the OBIS API ------------------

obis_occ_raw <- occurrence(taxon_name)

# Search for occurrences of the taxon using the GBIF API ------------------

# Extract the usageKey for each taxonomic name using name_backbone()
usageKeys <- lapply(taxon_name, function(x) name_backbone(x)$usageKey)

source(file = "~/Library/CloudStorage/OneDrive-UniversitéLaval/Script_R/gbifcode.R")

# Download the occurrence data
gbif_download <- occ_download(pred_in("taxonKey", usageKeys), user = username, pwd = password, email = courriel)

occ_download_wait(gbif_download)

# Convert into dataframe
gbif_occ_raw <-
  occ_download_get(gbif_download) %>%
  occ_download_import()

# Column selection --------------------------------------------------------

col <- c(
  "basisOfRecord",
  "bibliographicCitation",
  "catalogNumber",
  "collectionCode",
  "datasetID",
  "datasetName",
  "day",
  "decimalLatitude",
  "decimalLongitude",
  "depth",
  "eventDate",
  "eventID",
  "eventRemarks",
  # "fieldNotes",
  "identificationQualifier",
  "identificationReferences",
  "identificationRemarks",
  "individualCount",
  "institutionCode",
  "materialSampleID",
  "month",
  "nomenclaturalCode",
  "occurrenceID",
  "occurrenceRemarks",
  "occurrenceStatus",
  # "organismID",
  "organismQuantity",
  "organismQuantityType",
  "references",
  "sampleSizeValue",
  "sampleSizeValue",
  "samplingEffort",
  "samplingProtocol",
  "scientificName",
  "species",
  "year"
)

# Filter out occurrences that do not contain GPS coordinates and/or year information

gbif_occ <-
  gbif_occ_raw %>%
  mutate_all(~ paste(., "  /n  ")) %>%
  mutate_all(~ str_replace_all(., "/n", "")) %>%
  mutate_all(str_trim) %>%
  mutate_all(~ ifelse(. == "NA", NA, .)) %>%
  mutate_all(~ ifelse(. == "", NA, .)) %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude), !is.na(eventDate), !is.na(depth)) %>%
  mutate(
    sourceArchive = "GBIF",
    year = as.numeric(year),
    individualCount = as.numeric(individualCount),
    modified = as.Date(modified),
    day = as.numeric(day),
    month = as.numeric(month),
    eventDate = as.Date(eventDate),
    depth = as.numeric(depth),
    organismQuantity = as.numeric(organismQuantity),
    sampleSizeValue = as.double(sampleSizeValue),
    decimalLatitude = as.numeric(decimalLatitude),
    decimalLongitude = as.numeric(decimalLongitude)
  ) %>%
  select(all_of(col)) %>%
  rename("verbatimscientificName" = "scientificName", "scientificName" = "species")


obis_occ <-
  obis_occ_raw %>%
  mutate_all(~ paste(., "  /n  ")) %>%
  mutate_all(~ str_replace_all(., "/n", "")) %>%
  mutate_all(str_trim) %>%
  mutate_all(~ ifelse(. == "NA", NA, .)) %>%
  mutate_all(~ ifelse(. == "", NA, .)) %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude), !is.na(eventDate), !is.na(depth)) %>%
  mutate(
    sourceArchive = "OBIS",
    year = as.numeric(year),
    individualCount = as.numeric(individualCount),
    modified = as.Date(modified),
    day = as.numeric(day),
    month = as.numeric(month),
    eventDate = as.Date(eventDate),
    depth = as.numeric(depth),
    organismQuantity = as.numeric(organismQuantity),
    sampleSizeValue = as.double(sampleSizeValue),
    decimalLatitude = as.numeric(decimalLatitude),
    decimalLongitude = as.numeric(decimalLongitude)
  ) %>%
  select(all_of(col)) %>%
  rename("verbatimscientificName" = "scientificName", "scientificName" = "species")

# Add data from PhytoBase -------------------------------------------------

dl_Righetti2020 <- pg_data("10.1594/PANGAEA.904397")

raw_Righetti2020 <- read_csv(dl_Righetti2020[[1]]$path)

Righetti2020 <-
  raw_Righetti2020 %>%
  mutate_all(~ ifelse(. == "", NA, .)) %>%
  mutate_all(~ ifelse(. == "NA", NA, .)) %>%
  filter(scientificName %in% taxon_name) %>%
  mutate(
    day = as.double(day),
    individualCount = as.double(individualCount),
    organismQuantity = as.double(organismQuantity)
  ) %>%
  select(-c(phylum, class))


# Add data from dataBaseline Arctic ---------------------------------------

file <- "dataBaseline_Arctic.csv"

get_from_zenodo(doi = "10.5281/zenodo.13376814", path = tempdir(), file = file)

# Define the path to the temporary directory
temp_directory <- tempdir()

# Specify the full path to your CSV file
file_path <- file.path(temp_directory, "dataBaseline_Arctic.csv")

# Load the CSV file into a data frame
dataBaseline_Arctic <- read.csv(file_path, dec = ",", sep = ";")

Phaeocystis_Arctic <-
  dataBaseline_Arctic %>%
  filter(scientificName %in% taxon_name)

# Taxonomy verification with worms  ---------------------------------------

verif_sp <-
  bind_rows(
    gbif_occ %>%
      distinct(scientificName) %>%
      pull(scientificName) %>%
      wm_records_names(., fuzzy = F, marine_only = F) %>%
      enframe() %>%
      unnest(cols = c(value)),
    obis_occ %>%
      distinct(scientificName) %>%
      pull(scientificName) %>%
      wm_records_names(., fuzzy = F, marine_only = F) %>%
      enframe() %>%
      unnest(cols = c(value))
  ) %>%
  distinct(scientificname, .keep_all = T) %>%
  select(-name)

# Merging the different dataset -------------------------------------------

data <-
  bind_rows(
    obis_occ,
    gbif_occ,
    Righetti2020
  ) %>%
  left_join(., verif_sp, by = c("scientificName" = "scientificname")) %>%
  mutate(
    decimalLatitude = as.numeric(decimalLatitude),
    decimalLongitude = as.numeric(decimalLongitude),
    depth = gsub("-", "", depth), # Records with negative recording depths were flagged and changed to positive, assuming that their original sign was mistaken.
    depth = as.numeric(depth),
    year = ifelse(is.na(year), year(as.POSIXct(eventDate, format = "%Y-%m-%d")), as.numeric(year)),
    month = ifelse(is.na(month), month(as.POSIXct(eventDate, format = "%Y-%m-%d")), as.numeric(month)),
    day = ifelse(is.na(day), day(as.POSIXct(eventDate, format = "%Y-%m-%d")), as.numeric(day)),
    occurrenceStatus = ifelse(is.na(occurrenceStatus), "PRESENT", toupper(occurrenceStatus)),
    basisOfRecord = case_when(
      basisOfRecord == "HumanObservation" ~ "HUMAN_OBSERVATION",
      basisOfRecord == "MaterialSample" ~ "MATERIAL_SAMPLE",
      basisOfRecord == "Occurrence" ~ "OCCURRENCE",
      basisOfRecord == "PreservedSpecimen" ~ "PRESERVED_SPECIMEN",
      TRUE ~ basisOfRecord
    )
  ) %>%
  filter(
    !depth > 100,
    !is.na(year),
    !occurrenceStatus %in% "ABSENT",
    !basisOfRecord %in% "FOSSIL_SPECIMEN"
  ) %>%
  distinct(decimalLatitude, decimalLongitude, depth, day, month, year, verbatimscientificName, scientificName, basisOfRecord, .keep_all = T)

data_sf_point <-
  data %>%
  distinct(valid_name, decimalLatitude, decimalLongitude) %>%
  mutate_if(is.double, as.character) %>% # convert double columns to character to avoid errors
  sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat")

world_map <- rnaturalearth::ne_countries(scale = "small", returnclass = c("sf"))

Map <-
  ggplot() +
  geom_sf(data = data_sf_point, aes(shape = valid_name, fill = valid_name, color = valid_name), alpha = 0.5, size = 2.5) +
  geom_sf(data = world_map, fill = "grey50") +
  ggtitle("Occurrence data") +
  labs(colour = "Species", shape = "Species", fill = "Species") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(21:24)) +
  coord_sf(crs = "+proj=moll") +
  ggthemes::theme_map(base_size = 15) +
  theme(legend.text = element_text(face = "italic"))

print(Map)

# ggsave("Distribution_Phaeocystis.tiff", Map, height = 10, width = 12.5, dpi = 300, units = "in")
#
# Lat_Distribution_Phaeocystis <-
#   data %>%
#   distinct(valid_name, decimalLatitude, decimalLongitude) %>%
#   ggplot(., aes(valid_name, decimalLatitude)) +
#   geom_violin(aes(fill = valid_name), alpha = 0.5) +
#   # geom_jitter(alpha = 0.5, width = 0.075, size = 0.5) +
#   geom_boxplot(width = 0.05) +
#   labs(y = "Latitute (°N)", x = "") +
#   scale_fill_brewer(palette = "Dark2") +
#   # scale_x_discrete(name = "", breaks = lab$`break`, labels = lab$label) +
#   scale_y_continuous(breaks = seq(-90, 90, by = 25)) +
#   # ylim(-90, 90) +
#   ggtitle("Latitudinal repartition") +
#   theme_base() +
#   theme(legend.position = "none", axis.text.x = element_text(face = "italic"))
#
# print(Lat_Distribution_Phaeocystis)
#
# ggsave("Lat_Distribution_Phaeocystis.tiff", Lat_Distribution_Phaeocystis, height = 5, width = 7.5, dpi = 300, units = "in")
#
# data %>%
#   count(valid_name, basisOfRecord) %>%
#   group_by(valid_name) %>%
#   arrange(desc(n), .by_group = TRUE)
#
# Lat_Distribution_basisOfRecord <-
#   data %>%
#   mutate(lat_cut = cut(decimalLatitude, breaks = seq(-90, 90, by = 15))) %>%
#   count(valid_name, lat_cut, basisOfRecord, sort = T) %>%
#   ggplot() +
#   geom_col(aes(x = n, y = lat_cut, fill = basisOfRecord)) +
#   scale_fill_brewer(palette = "Set2", na.value = "grey50") +
#   facet_wrap(~valid_name) +
#   labs(y = "Latitute (°N)", x = "Number of observation") +
#   ggtitle("Latitudinal distribution of the number of the specific nature of the record") +
#   theme_base(base_size = 15) +
#   theme(strip.text = element_text(face = "italic"), strip.background = element_rect(fill = "grey", color = "black"))
#
# print(Lat_Distribution_basisOfRecord)
#
# ggsave("Lat_Distribution_basisOfRecord.tiff", Lat_Distribution_basisOfRecord, height = 5, width = 7.5, dpi = 300, units = "in")
#
# data %>%
#   filter(valid_name %in% "Phaeocystis antarctica", between(decimalLatitude, 0, 90)) %>%
#   count(basisOfRecord)
