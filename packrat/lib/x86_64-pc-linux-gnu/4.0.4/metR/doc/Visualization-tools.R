## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  dev = "jpeg",
collapse = TRUE,
message = FALSE,
fig.width = 7,
comment = "#>")

## ---- message = FALSE---------------------------------------------------------
# Packages and data use throught 
library(metR)
library(ggplot2)
library(data.table)
temperature <- copy(temperature)
temperature[, air.z := Anomaly(air), by = .(lat, lev)]

## -----------------------------------------------------------------------------
ggplot(temperature[lon %~% 180], aes(lat, lev, z = air.z)) +
  geom_contour_fill(breaks = MakeBreaks(binwidth = 2, exclude = 0)) +
  scale_fill_divergent(breaks = MakeBreaks(binwidth = 2, exclude = 0))

## -----------------------------------------------------------------------------
data(volcano)
volcano <- setDT(reshape2::melt(volcano))
volcano[, value.gap := value]
volcano[(Var1 - 40)^2 + (Var2 - 35)^2 < 50, value.gap := NA]

ggplot(volcano, aes(Var1, Var2, z = value.gap)) +
    geom_contour_fill() 

## -----------------------------------------------------------------------------
ggplot(volcano, aes(Var1, Var2, z = value.gap)) +
    geom_contour_fill(na.fill = TRUE) 

## -----------------------------------------------------------------------------
data(surface)
ggplot(surface, aes(lon, lat)) +
  geom_point(aes(color = height))

ggplot(surface, aes(x, y)) +
  geom_point(aes(color = height))

## -----------------------------------------------------------------------------
proj_string <- "+proj=lcc +lat_1=-30.9659996032715 +lat_2=-30.9659996032715 +lat_0=-30.9660034179688 +lon_0=-63.5670013427734 +a=6370000 +b=6370000"

ggplot(surface, aes(x, y)) +
  geom_contour_fill(aes(z = height))

ggplot(surface, aes(x, y)) +
  geom_contour_fill(aes(z = height), proj = proj_string)

## -----------------------------------------------------------------------------
set.seed(42)

some_volcano <- volcano[sample(.N, .N/7)]   # select 70% of the points
some_volcano[, Var1 := Var1 + rnorm(.N)]    # add some random noise
some_volcano[, Var2 := Var2 + rnorm(.N)]

ggplot(some_volcano, aes(Var1, Var2)) +
  geom_point(aes(color = value))

ggplot(some_volcano, aes(Var1, Var2)) +
  geom_contour_fill(aes(z = value), kriging = TRUE) +
  geom_point(size = 0.2)

## -----------------------------------------------------------------------------
ggplot(temperature[lev %in% c(1000, 300)], aes(lon, lat, z = air.z)) +
  geom_contour_fill() +
  scale_fill_divergent() +
  facet_grid(~lev)

## -----------------------------------------------------------------------------
ggplot(temperature[lev %in% c(1000, 300)], aes(lon, lat, z = air.z)) +
  geom_contour_fill(global.breaks = FALSE) +
  scale_fill_divergent() +
  facet_grid(~lev)

## -----------------------------------------------------------------------------
ggplot(temperature[lev == 300], aes(lon, lat, z = air.z)) +
    geom_contour_fill()

ggplot(temperature[lev == 300], aes(lon, lat, z = air.z)) +
    geom_contour_fill(aes(fill = stat(level)))

## -----------------------------------------------------------------------------
ggplot(temperature[lev == 300], aes(lon, lat, z = air.z)) +
    geom_contour_fill(aes(fill = stat(level_d)))

## -----------------------------------------------------------------------------
ggplot(temperature[lev == 300], aes(lon, lat, z = air.z)) +
    geom_contour_fill() +
    geom_contour(color = "black") +
    geom_text_contour() +
    scale_fill_divergent() 

## -----------------------------------------------------------------------------
ggplot(temperature[lev == 300], aes(lon, lat, z = air.z)) +
    geom_contour_fill() +
    geom_contour2(color = "black") +
    geom_text_contour(stroke = 0.2) +
    scale_fill_divergent() 

## -----------------------------------------------------------------------------
ggplot(temperature[lev == 300], aes(lon, lat, z = air.z)) +
    geom_contour_fill() +
    geom_contour2(color = "black") +
    geom_text_contour(stroke = 0.2, label.placement = label_placement_random()) +
    scale_fill_divergent() 

## -----------------------------------------------------------------------------
ggplot(temperature[lev == 300], aes(lon, lat, z = air.z)) +
     geom_contour_fill() +
     geom_contour_tanaka() +
     scale_fill_divergent()

## -----------------------------------------------------------------------------
data(geopotential)    # geopotential height at 700hPa for the Southern Hemisphere. 

ggplot(geopotential[, gh.base := gh[lon == 120 & lat == -50], by = date][
    , .(correlation = cor(gh.base, gh)), 
    by = .(lon, lat)],
    aes(lon, lat, z = correlation)) +
    geom_contour_fill(breaks = MakeBreaks(0.1)) +
    stat_subset(aes(subset = correlation > 0.5),
                geom = "point", size = 0.1) +
    scale_fill_divergent() 

## -----------------------------------------------------------------------------
ggplot(volcano, aes(Var1, Var2, z = value.gap)) +
    geom_contour_fill(na.fill = TRUE) +
    stat_subset(aes(subset = is.na(value.gap)), geom = "raster", 
                fill = "#a56de2")

## -----------------------------------------------------------------------------
temperature[, c("t.dx", "t.dy") := Derivate(air.z ~ lon + lat,
                                            cyclical = c(TRUE, FALSE)), 
            by = lev]

(g <- ggplot(temperature[lev == 500], aes(lon, lat)) +
    geom_contour_fill(aes(z = air.z)) +
    geom_vector(aes(dx = t.dx, dy = t.dy), skip.x = 2, 
                skip.y = 1) +
    scale_mag())

## -----------------------------------------------------------------------------
g + coord_polar()

## -----------------------------------------------------------------------------
ggplot(temperature[lon %between% c(100, 200) & lat == -50], aes(lon, lev)) + 
    geom_arrow(aes(dx = dx(t.dx, lat), dy = dy(t.dy)), skip = 1) +
    scale_y_level() +
    scale_mag()

## -----------------------------------------------------------------------------
(g <- ggplot(temperature[lev == 500], aes(lon, lat)) +
     geom_contour_fill(aes(z = air.z)) +
     geom_streamline(aes(dx = t.dy, dy = -t.dx), L = 10, res = 2,   
                     arrow.length = 0.3, xwrap = c(0, 360)))

## -----------------------------------------------------------------------------
g + coord_polar()

## -----------------------------------------------------------------------------
ggplot(temperature[lev == 500], aes(lon, lat)) +
    geom_streamline(aes(dx = t.dy, dy = -t.dx, size = ..step.., alpha = ..step..,
                        color = sqrt(..dx..^2 + ..dy..^2)), arrow = NULL,
                    L = 10, res = 2, xwrap = c(0, 360), lineend = "round") + 
    scale_color_viridis_c(guide = "none") +
    scale_size(range = c(0, 1), guide = "none") +
    scale_alpha(guide = "none")

## -----------------------------------------------------------------------------
ggplot(temperature[lev == 300], aes(lon, lat, z = air.z)) +
    geom_contour_fill(aes(fill = stat(level)), breaks = c(-10, -8, -6, -2, -1, 0, 6, 8, 10)) +
  guides(fill = guide_colorsteps())

## -----------------------------------------------------------------------------
ggplot(temperature[lev == 300], aes(lon, lat, z = air.z)) +
  geom_contour_fill(aes(fill = stat(level)), breaks = c(-10, -8, -6, -2, -1, 0, 6, 8, 10)) +
  scale_fill_discretised()

## -----------------------------------------------------------------------------
ggplot(temperature[lev == 300], aes(lon, lat, z = air.z)) +
  geom_contour_fill(aes(fill = stat(level)), breaks = c(-10, -8, -6, -2, -1, 0, 6, 8, 10)) +
  scale_fill_divergent_discretised(midpoint = 3)

## ---- message = FALSE---------------------------------------------------------
# Plot made with base ggplot
(g <- ggplot(temperature[lon %~% 180], aes(lat, lev, z = air.z)) +
        geom_contour2(aes(color = ..level..)))

## -----------------------------------------------------------------------------
g + 
    scale_y_level() +
    scale_x_latitude(ticks = 15, limits = c(-90, 90)) +
    scale_color_divergent()

## -----------------------------------------------------------------------------
ggplot(volcano, aes(Var1, Var2)) +
    geom_relief(aes(z = value))

