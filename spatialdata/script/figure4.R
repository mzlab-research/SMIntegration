suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(showtext))
suppressMessages(library(ggtext))
suppressMessages(library(VennDiagram))
suppressMessages(library(Cairo))
suppressMessages(library(dplyr))
suppressMessages(library(gridExtra))
suppressMessages(library(magick))
suppressMessages(library(stringr))
suppressMessages(library(readr)) # for read_delim
suppressMessages(library(magrittr)) # for %<>% pipe
suppressMessages(library(Matrix))
#4a---------------------------------------------------------------------------------------
#' @title Figure 4a: Spatial Group Assignment
#' @description Visualizes the spatial distribution of Treatment (CA) and Control (PAG) groups.
#' @return Generates a PDF '4a.pdf'.

# font_add("sans", regular = "arial.ttf", italic = "ariali.ttf")
showtext_auto()

data=readRDS("../data/fig4a_group_assignment.rds")
data <- data |>
  group_by(groups) 
group_color<-c("treatment"="red","control"="blue","0"="grey")
p1 <-  ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(color = groups), size = 1) +
  guides(colour = guide_legend(title = "Group",override.aes = list(size=3), nrow = 10))+
  scale_color_manual(
    values = group_color,
    labels = c("0" = "Other", 
               "treatment" = "CA", 
               "control" = "PAG")
  ) +
  xlab("") + ylab("") +
  xlim(min(data$x),max(data$x))+
  ylim(min(data$y),max(data$y))+
  coord_equal()+
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    legend.title = element_text(family = "Arial"),
    legend.text = element_markdown(family = "Arial"),
    axis.title = element_text(family = "Arial"),
    axis.text = element_text(family = "Arial"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )


pdf("4a.pdf", width = 4, height = 6)
print(p1)
dev.off()
#4b---------------------------------------------------------------------------------------
#' @title Figure 4b: Differential Feature Counts
#' @description Bar plot showing the number of Up- and Down-regulated features for each sample type.
#' @return Generates a PDF '4b.pdf'.
data=readRDS("../data/fig4b_feature_counts.rds")
UD_Class <- unique(data$State)
UD_len <- length(UD_Class)

if(UD_len == 1){
  if(UD_Class == 'Up'){
    scale_fill_values <- "HotPink"
    scale_fill_labels <- "Up"
  }else if(UD_Class == 'Down'){
    scale_fill_values <- "LightSkyBlue"
    scale_fill_labels <- "Down"
  }
}else{
  scale_fill_values <- c("HotPink", "LightSkyBlue")
  scale_fill_labels <- c("Up","Down")
}
p2 <- ggplot(data,mapping = aes(x=Sample_type, y=value, fill=State)) +
  theme_bw()+ theme(panel.grid=element_blank())+
  geom_bar(stat="identity",width = 0.5, position = position_dodge(0.7)) +
  geom_text(aes(label=value),hjust = 0.3,vjust = 1,position = position_dodge(0.7))+
  scale_fill_manual(values = scale_fill_values, labels = scale_fill_labels)+
  labs(x="",y="")+
  ylim(0, 1000) +
  theme(
    text = element_text(family = "Arial"),
    legend.title = element_text(family = "Arial", size = 12),
    legend.text = element_markdown(family = "Arial", size = 12),
    axis.title = element_text(family = "Arial", size = 12, face = "bold"),
    axis.text = element_text(family = "Arial"),
    axis.text.y = element_text(size = 10, angle = 0),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  )
# theme(
#       text = element_text(size=15))

pdf("4b.pdf", width = 6, height = 6)
print(p2)
dev.off()

#4c---------------------------------------------------------------------------------------
#' @title Figure 4c: Pathway Overlap Venn Diagram
#' @description Visualizes the intersection of pathways annotated by metabolites and genes.
#' @return Generates a PDF '4c.pdf'.
vennplot_data=readRDS("../data/fig4c_pathway_venn.rds")

p3 <- venn.diagram(vennplot_data,
                   col = "transparent",
                   fill = c('yellow', 'skyblue'),
                   alpha = 0.5, cex = 2,
                   fontfamily = "Arial",
                   cat.cex = 2, cat.pos = 0,
                   cat.fontfamily = "Arial",
                   rotation.degree = 0,
                   resolution = 300, filename = NULL)

CairoPDF("4c.pdf", width = 6, height = 6 ,family = "Arial")
grid.draw(p3)
dev.off()



#4d---------------------------------------------------------------------------------------
#' @title Figure 4d: Pathway annotation Dot Plot
#' @description Visualizes top annotated pathways with bubble size representing count.
#' @return Generates a PDF '4d.pdf'.
df=readRDS("../data/fig4d_pathway_dotplot.rds")
plot_title = ""
plot_title_size = 15
text_size =10 
axis_title_size = 12
legend_title_size = 12
legend_text_size =12 
top_num=20
# Calculate Total Count for each Pathway to determine Top 20
pathway_counts <- df %>%
  dplyr::group_by(Pathway) %>%
  dplyr::summarise(TotalCount = sum(Count)) %>%
  dplyr::arrange(desc(TotalCount)) %>%
  dplyr::slice_head(n = top_num)

top_pathways <- pathway_counts$Pathway

# Filter original df to keep only Top 20 pathways
df <- df %>% 
  dplyr::filter(Pathway %in% top_pathways) %>%
  # Join with total counts for sorting
  dplyr::left_join(pathway_counts, by = "Pathway") %>%
  dplyr::mutate(Pathway = factor(Pathway, levels = rev(top_pathways))) # Sort Y-axis by Total Count
# Plot: Annotation Count Dot Plot
# Color and Shape both map to Types (Gene/Metabolite)
g <- ggplot(df)+
  geom_point(aes(x= Count, y= Pathway, colour= Types, shape = Types), size = 4)+ 
  theme_minimal() +
  labs(title="", x="Annotation Count", y="Pathway", colour="Omics Type", shape="Omics Type")+
  scale_shape_manual(values = c("Gene" = 16, "Metabolite" = 17, "Co-annotation" = 15)) + # Circle, Triangle
  scale_x_continuous(breaks = function(x) unique(floor(pretty(x)))) +
  scale_colour_manual(values = c("Gene" = "#377EB8", "Metabolite" = "#E41A1C", "Co-annotation" = "#4DAF4A")) + # Blue, Red, Green
  theme(
    text = element_text(family = "Arial"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 10, angle = 0),
    plot.title = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", size = 1), 
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    panel.border = element_blank()
  )


# Title styling
if(plot_title == ""){
  g <- g + theme(plot.title = element_blank())
}else{
  g <- g + theme(plot.title = element_text(size = plot_title_size, hjust = 0.5, face="bold"))
}

# Axis title styling
if(axis_title_size == 0){
  g <- g + theme(axis.title = element_blank())
}else{
  g <- g + theme(axis.title = element_text(size = axis_title_size, face="bold"))
}
pdf("4d.pdf", width = 8, height = 6) 
print(g)
dev.off()
#4e---------------------------------------------------------------------------------------
#' @title Figure 4e: KEGG Pathway Map
#' @description Overlays differential feature data onto a KEGG pathway map image.
#' Draws colored circles for metabolites and rectangles for genes based on regulation status.
#' @return Generates a PDF '4e.pdf'.
annotation_plotdata=readRDS("../data/fig4ef_kegg_map_data.rds")
m_fplot <- annotation_plotdata[[1]]
t_fplot <- annotation_plotdata[[2]]
pathway_id <- annotation_plotdata[[5]]
species <- "mmu"
devs <- dev.list()
if (!is.null(devs)) {
  for (d in devs) {
    tryCatch({
      dev.off(d)
    }, error = function(e) {})
  }
}

mappath <- "../source/Database/KEGG/map"
img <- image_read(file.path(mappath, paste0(pathway_id, ".png")))
img <- image_convert(img, colorspace = 'gray')

lines <- suppressWarnings(readLines(file.path(mappath, paste0(pathway_id, ".conf"))))
dim <- image_info(img)
width <- dim$width
height <- dim$height

canvas <- image_blank(width, height, "none")

img_draw <- image_draw(canvas)
on.exit(dev.off(), add = TRUE)

circle <- function(x, y, r, border = "black", lty = "solid", lwd = 1) {
  angle <- seq(0, 2 * pi, length.out = 100)
  x_circle <- x + r * cos(angle)
  y_circle <- y + r * sin(angle)
  lines(x_circle, y_circle, col = border, lty = lty, lwd = lwd)
}

for (line in lines) {
  line <- trimws(line)
  
  if (grepl("^circ", line) && any(!is.na(m_fplot))) {
    coords <- as.numeric(str_extract_all(line, "\\d+")[[1]])
    id_match <- str_match(line, "(C\\d{5})")
    matched_ids <- intersect(id_match, names(m_fplot))
    if (length(matched_ids) > 0 && length(coords) >= 3) {
      x <- coords[1]
      y <- coords[2]
      r <- coords[3]
      
      values <- unlist(m_fplot[matched_ids])
      median_value <- median(values)
      
      if (!is.na(median_value)) {
        color <- ifelse(median_value > 0, "orange",
                        ifelse(median_value < 0, "blue", "gray"))
        circle(x, y, r, border = color, lty = "solid", lwd = 3)
      }
    }
  } else if (grepl("^rect", line) && any(!is.na(t_fplot))) {      
    coords <- as.numeric(str_extract_all(line, "\\d+")[[1]])
    ids <- str_extract_all(line, "K\\d{5}")[[1]]
    matched_ids <- intersect(ids, names(t_fplot))
    
    if (length(matched_ids) > 0 && length(coords) >= 4) {
      xleft <- coords[1]
      ybottom <- coords[2]
      xright <- coords[3]
      ytop <- coords[4]
      
      values <- unlist(t_fplot[matched_ids])
      median_value <- median(values)
      if (!is.na(median_value)) {
        color <- ifelse(median_value > 0, "red",
                        ifelse(median_value < 0, "green", "gray"))
        rect(xleft, ybottom, xright, ytop, border = color, lty = "solid", lwd = 3)
      }
    }
  }
}


img <- image_composite(img, img_draw, offset = "+0+0")

output_file <-  file.path( "4e.pdf")
target_width <- 5  
target_height <- 6  
dpi <- 300  
width_px <- target_width * dpi
height_px <- target_height * dpi

img_adjusted <- img %>%
  image_scale(geometry = sprintf("%dx%d", width_px, height_px)) %>% 
  image_extent(geometry = sprintf("%dx%d", width_px, height_px),   
               gravity = "center",                                  
               color = "white")                                  

image_write(
  img_adjusted,
  path = output_file,
  format = "pdf",
  density = dpi  
)
#4f---------------------------------------------------------------------------------------
#' @title Figure 4f: Feature Spatial Distribution
#' @description Visualizes the spatial expression/abundance of key pathway genes and metabolites.
#' @return Generates a PDF '4f.pdf'.
devs <- dev.list()
if (!is.null(devs)) {
  for (d in devs) {
    tryCatch({
      dev.off(d)
    }, error = function(e) {})
  }
}

annotation_plotdata=readRDS("../data/fig4ef_kegg_map_data.rds")
data_rds=readRDS("../data/fig4f_data_rds.rds")
m<-annotation_plotdata[[3]]
g<-annotation_plotdata[[4]]
# Prkacb PKA
# Slc17a6 vGLUT
# Slc32a1 VGAT
# Gng4 Gi/o
m1=data_rds[[1]]
m2=data_rds[[2]]
i1=data_rds[[3]]
i2=data_rds[[4]]
#' @title Generate Spatial Feature Plots
#' @description Creates spatial scatter plots for a list of features (genes/metabolites).
#' Colors points by normalized intensity.
#' @param feature Vector of feature IDs to plot.
#' @param meta.data Data frame containing spatial coordinates (x, y).
#' @param counts Expression/Intensity matrix.
#' @param size Point size.
#' @param featurename Optional custom titles for features.
#' @param anno_name Optional annotation to append to titles.
#' @return A list of ggplot objects.
runplot <- function(feature, meta.data, counts, size = 1, featurename = NULL, anno_name = NULL) {
  if (all(is.null(featurename))) {
    featurename <- feature
  }
  
  plots <- list()
  
  for (i in 1:length(feature)) {
    plotdata <- meta.data
    if(!is.null(nrow(counts))){
      plotdata$intensity <- counts[feature[i], ]
    }else{
      plotdata$intensity <- as.numeric(counts)
    }
    
    plotdata$norm_intensity <- 100 * (plotdata$intensity) / max(plotdata$intensity)
    heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
    
    p <- ggplot(plotdata, aes(x = x, y = y)) +
      geom_point(aes(color = norm_intensity), size = size) +
      scale_color_gradientn(
        colours = heatmap_Palette(100),
        name = "Relative Abundance",
        guide = guide_colourbar(
          title.theme = element_text(size = 10, family = "Arial"),  
          label.theme = element_text(family = "Arial") 
        )
      ) +
      xlim(min(plotdata$x), max(plotdata$x)) +
      ylim(min(plotdata$y), max(plotdata$y)) +
      coord_equal() +
      xlab("") + ylab("") +
      theme_minimal() +
      theme(
        text = element_text(family = "Arial"),  
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
      )
    
    if (all(is.null(anno_name))) {
      plots[[i]] <- p + 
        labs(title = featurename[i]) +
        theme(plot.title = element_text(hjust = 0.5, size = 10, family = "Arial")) 
    } else {
      title_expr <- bquote(italic(.(featurename[i])) ~ "(" * .(anno_name[i]) * ")")
      plots[[i]] <- p + 
        labs(title = title_expr) +
        theme(plot.title = element_text(hjust = 0.5, size = 10, family = "Arial")) 
    }
  }
  return(plots)
}

mp<-runplot(feature=m$metabolite,meta.data=m1,counts=i1,size=0.1,featurename="GABA")
tp<-runplot(feature=g$gene,meta.data=m2,counts=i2,size=0.1,featurename=c("Prkacb","Slc17a6","Slc32a1","Gng4"),anno_name=c("PKA","vGLUT","VGAT","Gi/o"))
pdf("4f.pdf", width = 10, height = 6)
do.call(gridExtra::grid.arrange, c(c(mp,tp), ncol = 3))
dev.off()


# mp<-runplot(feature="gamma-Aminobutyric acid",meta.data=m1,counts=i1,size=0.5)
# tp<-runplot(feature="Slc32a1",meta.data=m2,counts=i2,size=0.5)
# pdf("4f-2.pdf", width = 8, height = 6)
# do.call(gridExtra::grid.arrange, c(c(mp,tp), ncol = 2))
# dev.off()

#4G---------------------------------------------------------------------------------------
#' @title Figure 4g: Multi-Modal Co-Visualization
#' @description Creates an RGB overlay plot (p5) and a categorical co-expression plot (p6) 
#' for selected metabolite (GABA) and gene (Slc32a1).
#' @return Generates a PDF '4g.pdf'.
spatial_coordlist=readRDS("data/fig4g_co_visualization.rds")
spatial_coord=spatial_coordlist[[1]]
pair=spatial_coordlist[[2]]
pt.size=0.5
pairname=c("GABA",NA,"Slc32a1")
LRpair =pair[!is.na(pair)]
data = spatial_coord[,c('x','y',LRpair)] 
if(!is.na(pair[1])){
  f1=data[,pair[1]]
}else{
  f1=0
}
if(!is.na(pair[2])){
  f2=data[,pair[2]]
}else{
  f2=0
}
if(!is.na(pair[3])){
  f3=data[,pair[3]]
}else{
  f3=0
}


data$color <- rgb(f1,f2,f3, maxColorValue = 1)


p5<-ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(color = color), size = pt.size) +
  xlim(min(data$x), max(data$x)*1.3) +
  ylim(min(data$y), max(data$y)) +
  coord_equal() +
  scale_color_identity() +  
  xlab("")+ylab("")+
  theme_minimal() +
  theme(legend.text = element_text(size = 4 ) ,
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) 
if(!is.na(pair[1])){
  
  p5<-p5 + annotate("rect", xmin = max(data$x)+1, xmax =max(data$x)+ 3, ymin = max(data$y)*1.4/2, ymax = max(data$y)*1.4/2+4, 
                    fill = "red", color = "black") +  
    annotate("text", x = max(data$x)+ 5, y = max(data$y)*1.4/2+2, label = pairname[1], 
             hjust = 0, vjust = 0.5, size = 3, color = "black", family = "Arial")   
  
}
if(!is.na(pair[2])){
  
  p5<-p5 +annotate("rect", xmin = max(data$x)+1, xmax =max(data$x)+ 3, ymin = max(data$y)/2, ymax = (max(data$y)/2)+4, 
                   fill = "green", color = "black") +  
    annotate("text", x = max(data$x)+ 5, y = (max(data$y)/2)+2, label = pairname[2], 
             hjust = 0, vjust = 0.5, size = 3, color = "black", family = "Arial")   
}
if(!is.na(pair[3])){
  
  p5<-p5 +annotate("rect", xmin = max(data$x)+1, xmax =max(data$x)+ 3, ymin = max(data$y)*1.2/2, ymax = (max(data$y)*1.2/2)+4, 
                   fill = "blue", color = "black") +  
    annotate("text", x = max(data$x)+ 5, y =  (max(data$y)*1.2/2)+2, label =pairname[3], 
             hjust = 0, vjust = 0.5, size = 3, color = "black", family = "Arial")  #, fontface = "italic"
}
#
pairname=c("GABA",NA,"Slc32a1")
f1_col=NULL
f2_col=NULL
f3_col=NULL
pt.size=0.5
Fpair =pair[!is.na(pair)]

location = spatial_coord[,c('x','y')]
topn=floor(0.2*dim(location)[1])
expr =spatial_coord[,Fpair,drop=FALSE]
ncell<-dim(expr)[1]

expstatus<-function(f,fname){
  f_q25 = quantile(f, probs = 0.25)
  f_q75 = quantile(f, probs = 0.75)
  n1 = which(f > f_q75)
  n2 = which(f < f_q25)
  f_col<-rep(paste0(fname,"_medium"),ncell)
  f_col[n1]<-paste0(fname,"_high")
  f_col[n2]<-paste0(fname,"_low")
  return(f_col)
}
if(!is.na(pair[1])){
  f1<-expr[,pair[1]]
  f1_col<-expstatus(f1,fname=pairname[1])
}
if(!is.na(pair[2])){
  f2<-expr[,pair[2]]
  f2_col<-expstatus(f2,fname=pairname[2])
}
if(!is.na(pair[3])){
  f3<-expr[,pair[3]]
  f3_col<-expstatus(f3,fname=pairname[3])
}


expcol <- paste0(
  ifelse(nzchar(f1_col), paste0(f1_col, "_"), ""),
  ifelse(nzchar(f2_col), paste0(f2_col, "_"), ""),
  ifelse(nzchar(f3_col), f3_col, "")
)

expcol <- gsub("_$", "", expcol)

tmp<-data.frame(x=location[,1],y=location[,2],Exp=as.factor(expcol))

legend_labels <- levels(tmp$Exp)
gene_name <- "Slc32a1"
italic_labels <- gsub(
  gene_name, 
  paste0("<i>", gene_name, "</i>"), 
  legend_labels
)

tmp$Exp <- factor(tmp$Exp, labels = italic_labels)

p6 <- ggplot(tmp, aes(x = x, y = y, col = Exp)) +
  geom_point(size = pt.size) +
  labs(color = "") + 
  xlim(min(tmp$x), max(tmp$x)) +
  ylim(min(tmp$y), max(tmp$y)) +
  coord_equal() +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(
    legend.text = element_markdown(size = 4, family = "Arial"),  
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    text = element_text(family = "Arial")  
  )
pdf("4g.pdf", width = 10, height = 6)
gridExtra::grid.arrange(p5,p6, ncol = 2)
dev.off()
