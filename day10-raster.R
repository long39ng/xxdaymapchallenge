library(tidyverse)
library(httr)
library(sf)
library(stars)
library(gstat)
library(patchwork)
library(magick)


# Data --------------------------------------------------------------------

get_list <- function(url, query = NULL) {
  GET(url, query = query) |>
    content("text", encoding = "UTF-8") |>
    jsonlite::fromJSON(simplifyVector = FALSE)
}

extract_no2_tbl <- function(list) {
  list$data |>
    map_dfr(\(x) set_names(x[1:2], c("station", "value")))
}

extract_station_tbl <- function(list) {
  list$stations |>
    map_dfr(\(x) {
      set_names(
        x[c(1, 8, 9, 15, 16, 17)],
        c("station", "long", "lat", "setting", "setting_short", "type")
      )
    })
}

get_list("https://www.umweltbundesamt.de/api/air_data/v2/components/json") |>
  listviewer::jsonedit() # Component ID for NO2 is 5

annual_no2 <- 2001:2020 |>
  map_dfr(\(x) {
    Sys.sleep(1)
    get_list(
      "https://www.umweltbundesamt.de/api/air_data/v2/annualbalances/json",
      list(component = 5, year = x)
    ) |>
      extract_no2_tbl() |>
      mutate(year = x)
  }) |>
  mutate(value = as.numeric(value))

stations <- get_list(
  "https://www.umweltbundesamt.de/api/air_data/v2/meta/json",
  list(use = "measure")
) |>
  extract_station_tbl()

de_utm <- giscoR::gisco_get_nuts(country = "Germany", nuts_level = 1, resolution = "10") |>
  st_transform(st_crs(32632))

de_utm_rast <- de_utm |>
  st_cast("MULTILINESTRING") |>
  st_rasterize(dx = 1e4) |>
  st_as_sf(as_points = TRUE)

de_grid <- de_utm |>
  st_bbox() |>
  st_as_stars(dx = 1e4) |>
  st_crop(de_utm)

annual_no2_utm <- annual_no2 |>
  inner_join(stations, by = "station") |>
  drop_na(value) |>
  # Drop traffic and industry stations
  filter(type == "background") |>
  st_as_sf(coords = c("long", "lat"), crs = st_crs(4326)) |>
  st_transform(st_crs(32632)) |>
  nest(data = -year)


# Interpolation -----------------------------------------------------------

annual_no2_krige <- annual_no2_utm |>
  mutate(
    vgram = map(data, \(.x) {
      variogram(value ~ 1, .x)
    }),
    vgram_fit = map(vgram, \(.x) {
      fit.variogram(.x, vgm(1, "Ste", 5e4, 1))
    }),
    results = map2(data, vgram_fit, \(.x, .y) {
      krige(value ~ 1, .x, de_grid, .y)
    })
  )

map2(annual_no2_krige$vgram, annual_no2_krige$vgram_fit, plot)


# Plotting ----------------------------------------------------------------

combined_range <- function(list, var) {
  ranges <- map_dfr(list, \(.x) {
    range(.x[[var]], na.rm = TRUE) |>
      set_names(c("min", "max"))
  })
  c(min(ranges$min), max(ranges$max))
}

dots_stations_by_year <- function(list) {
  map(list, \(.x) {
    .x |>
      ggplot(aes(value, 1, fill = value)) +
      ggfx::with_shadow(
        ggbeeswarm::geom_quasirandom(
          groupOnX = FALSE,
          shape = 21,
          colour = alpha("white", .5),
          size = 2
        ),
        x_offset = 3, y_offset = 3, sigma = 3
      ) +
      scale_x_continuous(labels = \(x) {
        if_else(
          x == max(x, na.rm = TRUE),
          paste(strrep(" ", 10), x, "µg/m³"),
          as.character(x)
        )
      }) +
      scale_fill_gradientn(colours = pal, limit = combined_range(list, "value")) +
      coord_cartesian(
        xlim = combined_range(list, "value"),
        clip = "off"
      ) +
      labs(
        subtitle = "Kriging interpolation using NO₂ measurements from\nbackground monitoring stations",
        caption = "Data retrieved via the Umweltbundesamt Air Data API"
      ) +
      theme_void(base_family = "Red Hat Text", base_size = 20) +
      theme(
        axis.text.x = element_text(colour = text_col, margin = margin(3)),
        axis.ticks.x = element_line(text_col, .6),
        axis.ticks.length.x = unit(4, "mm"),
        plot.subtitle = element_text(
          colour = text_col,
          hjust = .5,
          lineheight = 1.4,
          margin = margin(b = 20)
        ),
        plot.caption = element_text(
          face = "italic",
          colour = text_col,
          margin = margin(20),
        ),
        plot.margin = margin(0, 30, 20, 30)
      )
  })
}

map_kriging_by_year <- function(list) {
  map(list, \(.x) {
    .x <- .x |>
      mutate(visible = !is.na(var1.pred))

    ggplot() +
      ggfx::with_shadow(
        geom_stars(data = .x, aes(x, y, fill = var1.pred, alpha = visible)),
        x_offset = 8, y_offset = 8, sigma = 5
      ) +
      geom_sf(
        data = de_utm_rast,
        colour = alpha("white", .5),
        size = .25
      ) +
      scale_fill_gradientn(
        limit = combined_range(list, "var1.pred"),
        colours = pal
      ) +
      # Make NA raster cells transparent
      scale_alpha_discrete(range = c(0, 1)) +
      theme_void() +
      theme(plot.margin = margin())
  })
}

text_col <- alpha("white", .9)
bg_col <- "#1c1c1c"
pal <- c("#0E1214", "#405258", "#759095", "#ACCACD", "#E3FEFF")

krige_plots <- map_kriging_by_year(annual_no2_krige$results)
station_plots <- dots_stations_by_year(annual_no2_krige$data)

gif_frames <- image_graph(800, 1080, bg = bg_col, res = 80)

imap(
  2001:2020,
  \(.x, .y) {
    krige_plots[[.y]] / station_plots[[.y]] +
      plot_layout(heights = c(7, 1)) +
      plot_annotation(
        title = paste("Air pollution in Germany –", .x)
      ) &
      theme(
        legend.position = "none",
        plot.title = element_text(
          family = "Red Hat Display",
          face = "bold",
          size = 36,
          colour = "white",
          margin = margin(t = 40, l = 30)
        ),
        plot.background = element_rect(fill = bg_col, size = 0)
      )
  }
)

gif <- gif_frames |>
  image_animate(delay = 60, optimize = TRUE)

image_write(gif, "day10-de-no2.gif")
