#!/usr/bin/env Rscript

library(sf)
library(dplyr)
library(ggplot2)
library(geodata)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

# -----------------------------
# Paths
# -----------------------------
base_dir <- "map"

input_shp <- file.path(
  base_dir,
  "redlist_species_data_fd293426-0bea-4397-8f6b-85447aba5418",
  "data_0.shp"
)

cache_dir <- file.path(base_dir, "map_cache")

output_pdf <- file.path(
  base_dir,
  "panda_iucn_mountain_river_context.pdf"
)

output_shp <- file.path(
  base_dir,
  "panda_iucn_extant_zoom.shp"
)

dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Map extent
# -----------------------------
xmin <- 101.5
xmax <- 109.5
ymin <- 27.5
ymax <- 34.5

bbox_sf <- st_as_sfc(
  st_bbox(
    c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    crs = st_crs(4326)
  )
)

bbox_vect <- terra::vect(
  terra::ext(xmin, xmax, ymin, ymax),
  crs = "EPSG:4326"
)

# -----------------------------
# Administrative boundaries for main map
# -----------------------------
china0 <- geodata::gadm(
  country = "CHN",
  level = 0,
  path = cache_dir
)
china1 <- geodata::gadm(
  country = "CHN",
  level = 1,
  path = cache_dir
)
china2 <- geodata::gadm(
  country = "CHN",
  level = 3,
  path = cache_dir
)

china0 <- st_as_sf(china0) |>
  st_make_valid() |>
  st_transform(4326)
china1 <- st_as_sf(china1) |>
  st_make_valid() |>
  st_transform(4326)
china2 <- st_as_sf(china2) |>
  st_make_valid() |>
  st_transform(4326)

target_provinces <- c(
  "Sichuan",
  "Shaanxi",
  "Gansu",
  "Chongqing",
  "Yunnan",
  "Guizhou"
)

china1_main <- china1 %>%
  filter(NAME_1 %in% target_provinces)
china2_main <- china2 %>%
  filter(NAME_1 %in% target_provinces)
china0_crop <- st_crop(china0, bbox_sf)
china1_crop <- st_crop(china1_main, bbox_sf)
china2_crop <- st_crop(china2_main, bbox_sf)

# -----------------------------
# Rivers
# -----------------------------
# rivers <- ne_download(
#   scale = 50,
#   type = "rivers_lake_centerlines",
#   category = "physical",
#   returnclass = "sf"
# )
# 
# rivers <- rivers |>
#   st_make_valid() |>
#   st_transform(4326) |>
#   st_crop(bbox_sf)

# -----------------------------
# Elevation background
# -----------------------------
# This downloads DEM data for China and caches it locally.
dem <- geodata::elevation_30s(
  country = "CHN",
  path = cache_dir
)

dem_crop <- terra::crop(dem, bbox_vect)

# Convert elevation raster to data frame
elev_df <- as.data.frame(
  dem_crop,
  xy = TRUE,
  na.rm = FALSE
)

names(elev_df) <- c("lon", "lat", "elev")

# -----------------------------
# IUCN panda range
# -----------------------------
panda_range <- st_read(input_shp, quiet = TRUE)

panda_range <- panda_range |>
  st_make_valid() |>
  st_transform(4326)

presence_field <- intersect(
  c("presence", "PRESENCE", "presence_c", "PRESENCE_C"),
  names(panda_range)
)[1]

if (is.na(presence_field)) {
  stop("No presence field found in the IUCN shapefile. Please check names(panda_range).")
}

# Keep extant and probably extant; remove possibly extant
if (is.numeric(panda_range[[presence_field]]) || is.integer(panda_range[[presence_field]])) {
  panda_range <- panda_range %>%
    filter(.data[[presence_field]] %in% c(1, 2))
} else {
  pres <- tolower(trimws(as.character(panda_range[[presence_field]])))
  panda_range <- panda_range %>%
    filter(pres %in% c("extant", "probably extant"))
}

panda_range_crop <- st_crop(panda_range, bbox_sf)

st_write(
  panda_range_crop,
  output_shp,
  delete_dsn = TRUE,
  quiet = TRUE
)

# -----------------------------
# Main map
# -----------------------------
p <- ggplot() +
  # Soft elevation background
  geom_raster(
    data = elev_df,
    aes(x = lon, y = lat, fill = elev),
    alpha = 0.58
  ) +
  scale_fill_gradientn(
    colors = c(
      "#EEF4E3",
      "#C9D9A3",
      "#E3D29A",
      "#D5B47B",
      "#F4F1EA"
    ),
    na.value = "white",
    guide = "none"
  ) +
  
  # National boundary
  geom_sf(
    data = china0_crop,
    fill = NA,
    color = "grey50",
    linewidth = 0.30
  ) +
  # Prefecture-level boundaries
  geom_sf(
    data = china2_crop,
    fill = NA,
    color = "grey88",
    linewidth = 0.12
  ) +
  # Province boundaries
  geom_sf(
    data = china1_crop,
    fill = NA,
    color = "grey70",
    linewidth = 0.25
  ) +
  
  # Rivers
  # geom_sf(
  #   data = rivers,
  #   color = "#9BCBEA",
  #   linewidth = 0.18,
  #   alpha = 0.45
  # ) +
  
  # Panda IUCN range
  geom_sf(
    data = panda_range_crop,
    fill = "#D9D9D9",
    color = "#3F3F3F",
    linewidth = 0.30,
    alpha = 0.88
  ) +
  
  coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    expand = FALSE
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

p

ggsave(
  filename = output_pdf,
  plot = p,
  device = cairo_pdf,
  width = 9,
  height = 7
)