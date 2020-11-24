histogram_net  <- function(EX, MS, BS, param){
  
  data_EX = EX
  data_MS = MS
  data_BS = BS  
  
  datDir = file.path(createPaths(), 'RData')
  list2env(loading(file.path(datDir, 'histogram_M_EX.RData')),envir=.GlobalEnv)
  list2env(loading(file.path(datDir, 'histogram_T_EX.RData')), envir=.GlobalEnv)
  list2env(loading(file.path(datDir, 'histogram_M_MS.RData')),envir=.GlobalEnv)
  list2env(loading(file.path(datDir, 'histogram_T_MS.RData')),envir=.GlobalEnv)
  list2env(loading(file.path(datDir, 'histogram_M_BS.RData')),envir=.GlobalEnv)
  list2env(loading(file.path(datDir, 'histogram_T_BS.RData')),envir=.GlobalEnv)
  
  list2env(loading(file.path(datDir, 'histogram_ALL_M_EX.RData')),envir=.GlobalEnv)
  list2env(loading(file.path(datDir, 'histogram_ALL_T_EX.RData')),envir=.GlobalEnv)
  list2env(loading(file.path(datDir, 'histogram_ALL_M_MS.RData')),envir=.GlobalEnv)
  list2env(loading(file.path(datDir, 'histogram_ALL_T_MS.RData')),envir=.GlobalEnv)
  list2env(loading(file.path(datDir, 'histogram_ALL_M_BS.RData')),envir=.GlobalEnv)
  list2env(loading(file.path(datDir, 'histogram_ALL_T_BS.RData')),envir=.GlobalEnv)
  
  posCol_M = '#003c30'
  posCol_T = '#67001f'
  negCol_M = '#543005'
  negCol_T = '#053061'
  quantCol = c('#fed976', '#fd8d3c', '#fc4e2a')
  eventYrs = param$eventYrs 
  figDir = file.path(createPaths(), 'histograms')
  
  ## Excursion
  
  p1 = ggplot() + 
    geom_segment(aes(x = eventYrs,
                     y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = posEvents_T_EX), 
                 color = posCol_T) +
    geom_segment(aes(x = eventYrs, y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = -negEvents_T_EX),
                 color = negCol_T) +
    geom_col(aes(x = eventYrs, y = posDiff_T_EX), fill = posCol_T) +
    geom_col(aes(x = eventYrs, y = negDiff_T_EX), fill = negCol_T) +
    geom_line(aes(x = eventYrs, y = quants_T_EX[1,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_T_EX[2,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_T_EX[3,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_T_EX[1,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_T_EX[2,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_T_EX[3,]), 
    #           color = quantCol[3]) +
    geom_line(aes(x = eventYrs, y = quants_T_EX[4,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_T_EX[5,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_T_EX[6,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_T_EX[4,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_T_EX[5,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_T_EX[6,]), 
    #           color = quantCol[3]) +
    scale_x_continuous(name = '', breaks = eventYrs[seq(1,25,by=2)], 
                       labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('') + #ylab('Excursion:\nFraction of events') +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 14), 
                       axis.text = element_text(size = 12), 
                       axis.title = element_text(size = 12)) +
    ggtitle('Net temperature events')
  
  p2 = ggplot() + 
    geom_segment(aes(x = eventYrs, 
                     y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = posEvents_M_EX), 
                 color = posCol_M) +
    geom_segment(aes(x = eventYrs, y = rep(0,length(eventYrs)),
                     xend = eventYrs, yend = -negEvents_M_EX), 
                 color = negCol_M) +
    geom_col(aes(x = eventYrs, y = posDiff_M_EX), fill = posCol_M) +
    geom_col(aes(x = eventYrs, y = negDiff_M_EX), fill = negCol_M) +
    geom_line(aes(x = eventYrs, y = quants_M_EX[1,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_M_EX[2,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_M_EX[3,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_M_EX[1,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_M_EX[2,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_M_EX[3,]), 
    #           color = quantCol[3]) +
    geom_line(aes(x = eventYrs, y = quants_M_EX[4,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_M_EX[5,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_M_EX[6,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_M_EX[4,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_M_EX[5,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_M_EX[6,]), 
    #           color = quantCol[3]) +
    scale_x_continuous(name = '', breaks = eventYrs[seq(1,25,by=2)], 
                       labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('') +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 14), 
                       axis.text = element_text(size = 12)) +
    ggtitle('Net moisture events')
  
  # Mean Shift
  
  p3 = ggplot() + 
    geom_segment(aes(x = eventYrs, y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = posEvents_T_MS), 
                 color = posCol_T) +
    geom_segment(aes(x = eventYrs, y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = -negEvents_T_MS), 
                 color = negCol_T) +
    geom_col(aes(x = eventYrs, y = posDiff_T_MS), fill = posCol_T) +
    geom_col(aes(x = eventYrs, y = negDiff_T_MS), fill = negCol_T) +
    geom_line(aes(x = eventYrs, y = quants_T_MS[1,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_T_MS[2,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_T_MS[3,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_T_MS[1,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_T_MS[2,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_T_MS[3,]), 
    #           color = quantCol[3]) +
    geom_line(aes(x = eventYrs, y = quants_T_MS[4,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_T_MS[5,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_T_MS[6,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_T_MS[4,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_T_MS[5,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_T_MS[6,]), 
    #           color = quantCol[3]) +
    scale_x_continuous(name = '', breaks = eventYrs[seq(1,25,by=2)], 
                       labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('') + #ylab('Regime change:\nFraction of events') +
    theme_bw() + theme(axis.text = element_text(size = 12), 
                       axis.title = element_text(size = 12))
  
  p4 = ggplot() + 
    geom_segment(aes(x = eventYrs, y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = posEvents_M_MS), 
                 color = posCol_M) +
    geom_segment(aes(x = eventYrs, y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = -negEvents_M_MS), 
                 color = negCol_M) +
    geom_col(aes(x = eventYrs, y = posDiff_M_MS), fill = posCol_M) +
    geom_col(aes(x = eventYrs, y = negDiff_M_MS), fill = negCol_M) +
    geom_line(aes(x = eventYrs, y = quants_M_MS[1,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_M_MS[2,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_M_MS[3,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_M_MS[1,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_M_MS[2,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_M_MS[3,]), 
    #           color = quantCol[3]) +
    geom_line(aes(x = eventYrs, y = quants_M_MS[4,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_M_MS[5,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_M_MS[6,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_M_MS[4,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_M_MS[5,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_M_MS[6,]), 
    #           color = quantCol[3]) +
    scale_x_continuous(name = '', breaks = eventYrs[seq(1,25,by=2)], 
                       labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('') +
    theme_bw() + theme(axis.text = element_text(size = 12))
  
  # Trend Changes
  
  p5 = ggplot() + 
    geom_segment(aes(x = eventYrs, y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = posEvents_T_BS), 
                 color = posCol_T) +
    geom_segment(aes(x = eventYrs, y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = -negEvents_T_BS), 
                 color = negCol_T) +
    geom_col(aes(x = eventYrs, y = posDiff_T_BS), fill = posCol_T) +
    geom_col(aes(x = eventYrs, y = negDiff_T_BS), fill = negCol_T) +
    geom_line(aes(x = eventYrs, y = quants_T_BS[1,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_T_BS[2,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_T_BS[3,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_T_BS[1,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_T_BS[2,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_T_BS[3,]), 
    #           color = quantCol[3]) +
    geom_line(aes(x = eventYrs, y = quants_T_BS[4,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_T_BS[5,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_T_BS[6,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_T_BS[4,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_T_BS[5,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_T_BS[6,]), 
    #           color = quantCol[3]) +
    scale_x_continuous(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], 
                       labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('') + #ylab('Trend change:\nFraction of events') +
    theme_bw() + theme(axis.text = element_text(size = 12), 
                       axis.title = element_text(size = 12))
  
  p6 = ggplot() + 
    geom_segment(aes(x = eventYrs, y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = posEvents_M_BS), 
                 color = posCol_M) +
    geom_segment(aes(x = eventYrs, y = rep(0,length(eventYrs)), 
                     xend = eventYrs, yend = -negEvents_M_BS), 
                 color = negCol_M) +
    geom_col(aes(x = eventYrs, y = posDiff_M_BS), fill = posCol_M) +
    geom_col(aes(x = eventYrs, y = negDiff_M_BS), fill = negCol_M) +
    geom_line(aes(x = eventYrs, y = quants_M_BS[1,]), color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_M_BS[2,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_M_BS[3,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_M_BS[1,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_M_BS[2,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_M_BS[3,]), 
    #           color = quantCol[3]) +
    geom_line(aes(x = eventYrs, y = quants_M_BS[4,]), 
              color = quantCol[1]) +
    geom_line(aes(x = eventYrs, y = quants_M_BS[5,]), color = quantCol[2]) +
    geom_line(aes(x = eventYrs, y = quants_M_BS[6,]), color = quantCol[3]) +
    #geom_point(aes(x = eventYrs, y = quants_M_BS[4,]), 
    #           color = quantCol[1]) +
    #geom_point(aes(x = eventYrs, y = quants_M_BS[5,]), 
    #           color = quantCol[2]) +
    #geom_point(aes(x = eventYrs, y = quants_M_BS[6,]), 
    #           color = quantCol[3]) +
    scale_x_continuous(name = 'ky BP', breaks = eventYrs[seq(1,25,by=2)], 
                       labels = eventYrs[seq(1,25,by=2)]/1000) +
    ylab('') +
    theme_bw() + theme(axis.text = element_text(size = 12))
  
  g1 = ggplotGrob(p1)
  g2 = ggplotGrob(p2)
  g3 = ggplotGrob(p3)
  g4 = ggplotGrob(p4)
  g5 = ggplotGrob(p5)
  g6 = ggplotGrob(p6)
  g_c1 = rbind(g1, g3, g5, size = "first")
  g_c2 = rbind(g2, g4, g6, size = "first")
  g_c1$widths <- unit.pmax(g1$widths, g3$widths, g5$widths)
  g_c2$widths <- unit.pmax(g2$widths, g4$widths, g6$widths)
  g = cbind(g_c1, g_c2, size = 'first')
    
  pdf(file.path(figDir, 'netT_doubleAxis.pdf'))
  grid.newpage()
  grid.draw(g_c1)
  dev.off()
  
  pdf(file.path(figDir, 'netM_doubleAxis.pdf'))
  grid.newpage()
  grid.draw(g_c2)
  dev.off()
  
  plots  <- list(g_c1,g_c2)
  
  return(plots)
}
