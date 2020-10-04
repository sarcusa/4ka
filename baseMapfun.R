baseMAP  <- function (lon, lat, map.type = "line", f = 0.1, 
                      restrict.map.range = TRUE,
                      projection = "mollweide", 
                      bound.circ = FALSE, global = FALSE, 
                      extend.range = 5, country.boundaries = T){
  if (length(lat) == 1 & length(lon) == 1) {
    lon = lon + c(-extend.range, extend.range)
    lat = lat + c(-extend.range, extend.range)
  }
  if (bound.circ & restrict.map.range) {
    warning("bound.circ and restrict.map.range can't both be true - setting restrict.map.range to FALSE")
    restrict.map.range = FALSE
  }
  if (global & restrict.map.range) {
    warning("global and restrict.map.range can't both be true - setting restrict.map.range to FALSE")
    restrict.map.range = FALSE
  }
  if (map.type == "google") {
    warning("google maps are no longer supported, using `stamen` instead")
    map.type <- "stamen"
  }
  if (map.type == "stamen") {
    try(load(file.path(tempdir(), "bound.box.Rdata")), silent = TRUE)
    if (!exists("bb")) {
      bb = 1000
    }
  }
  else {
    bb = 1000
  }
  bbnew <- ggmap::make_bbox(lon, lat, f = f)
  if (all(bbnew == bb)) {
    lnp = TRUE
  }
  else {
    lnp = FALSE
    bb = bbnew
    if (bb[4] > 90) {
      bb[4] = 90
    }
    save(bb, file = file.path(tempdir(), "bound.box.Rdata"))
  }
  if (map.type == "stamen") {
    if (global) {
      bb <- c(-180, -70, 180, 70)
      names(bb) <- c("left", "bottom", "right", "top")
      lnp <- FALSE
      warning("Stamen globally doesn't work well at high latitudes. Truncating at 70.")
    }
    if (lnp) {
      suppressWarnings(try(load(file.path(tempdir(), "newmap.Rdata")), 
                           silent = TRUE))
      if (!exists("newmap")) {
        newmap <- get_map(location = bb, maptype = "terrain", 
                          source = "stamen")
        save(newmap, file = file.path(tempdir(), "newmap.Rdata"))
      }
    }
    else {
      newmap <- get_map(location = bb, maptype = "terrain", 
                        source = "stamen")
      save(newmap, file = file.path(tempdir(), "newmap.Rdata"))
    }
    baseMap = ggmap(newmap, maprange = TRUE)
  }
  else if (map.type == "line") {
    res = 1
    low.lat <- min(lat) - 5
    if (restrict.map.range) {
      x_lim = bb[c(1, 3)]
      y_lim = bb[c(2, 4)]
    }
    else {
      x_lim = c(-190, 190)
      y_lim = c(-91, 91)
    }
    x_cell_lim <- x_lim + c(1, -1) * res/2
    y_cell_lim <- y_lim + c(1, -1) * res/2
    if (global) {
      if(country.boundaries){
        dum = maps::map(plot = FALSE)
    }else{
      dum = maps::map(plot = FALSE, interior = FALSE)
    }
      x_lim = c(-190, 190)
      y_lim = c(-91, 91)
    }
    else {
      dum = maps::map(xlim = x_lim, ylim = y_lim, plot = FALSE, 
                      wrap = TRUE)
    }
    ant_ggplot = data.frame(dum[c("x", "y")])
    if (min(x_lim) <= -180 & max(x_lim) >= 180) {
      badLines = which(ant_ggplot$y <= min(y_lim) | ant_ggplot$y >= 
                         max(y_lim))
    }
    else {
      badLines = which(ant_ggplot$x <= min(x_lim) | ant_ggplot$x >= 
                         max(x_lim) | ant_ggplot$y <= min(y_lim) | ant_ggplot$y >= 
                         max(y_lim))
    }
    ant_ggplot[badLines, ] = NA
    bbdf = data.frame(bb)
    baseMap <- ggplot() + geom_path(aes(x = x, y = y), data = ant_ggplot, 
                                    size = 0.1)
    if (restrict.map.range) {
      baseMap = baseMap + geom_rect(aes(xmax = bb[3] - 
                                          0.1, xmin = bb[1] + 0.1,
                                        ymax = bb[4] - 0.1, 
                                        ymin = bb[2] + 0.1), fill = NA, 
                                    colour = "black", 
                                    data = bbdf)
    }
    else if (bound.circ) {
      wbb = which(min(abs(bb[c(2, 4)])) == abs(bb[c(2, 
                                                    4)]))
      wbb = c(2, 4)[wbb]
      bound.circ = data.frame(x = seq(-180, 180), y = rep(bb[wbb] + 0.1, length.out = length(seq(-180, 180))))
      baseMap = baseMap + geom_path(aes(x = x, y = y), 
                                    data = bound.circ)
    }
    baseMap = baseMap + theme(panel.background = element_rect(fill = "white"), 
                              axis.ticks = element_blank(), 
                              axis.text.y = element_blank(), 
                              axis.text.x = element_blank(), 
                              axis.title.x = element_blank(), 
                              axis.title.y = element_blank(), 
                              panel.border = element_blank())
    if (global) {
      baseMap = baseMap + coord_map(projection, xlim = c(-180,180)) + 
        geom_rect(aes(xmax = 180.1, xmin = -180.1,ymax = 90.1, ymin = -90.1),
                  fill = NA, colour = "black")
    }
    else {
      baseMap = baseMap + coord_map(projection)
    }
  }
  else {
    stop(paste("Dont recognize map.type =", map.type))
  }
  return(baseMap)
}