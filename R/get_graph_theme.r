get_graph_theme <- function(){
  theTheme = theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(size = .5, color = "grey"),
        axis.line = element_line(size=.7, color = "black"),legend.position = c(0.85, 0.7), text = element_text(size=14),
        axis.title.x = element_text(vjust = -0.35), axis.title.y = element_text(vjust = 0.6)) + theme(legend.position="none")
  return(theTheme)
}