# plotColorPalettes is the function for generating the pallettes, 
# the colors in the palettes depends on the number we put in 
# in our case, the palettes is called plotcolor

plotColorPalettes <- function(g){
  d <- 360/g
  h <- cumsum(c(15, rep(d,g - 1)))
  hcl(h = h, c = 100, l = 65)
}

grid_arrange_shared_legend <- 
  function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
    plots <- list(...) 
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] +theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, 
                                                               lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, 
                                                             lwidth)))
    grid.newpage()
    grid.draw(combined)
    # return gtable invisibly
    invisible(combined)
    
  }

# functions to compute the max and min value 
##  Only useful for those has different percentiles
### n and m represent 2 different data frames that you want to compare with

find_max<-function(n, m)
{
  max_n <- max(max(subset(n, Percentile <= 95,select=c(value))), 
               max(subset(m, Percentile <= 95,select=c(value))))
  return(max_n)
}

#find_max(QCmetricsLongPrecMW_ALL, QCmetricsLongPrecMW_ID)

find_min <- function(n,m)
{
  min_n <- min(min(subset(n, Percentile <= 95,select=c(value))), 
               min(subset(m, Percentile <= 95,select=c(value))))
  return(min_n)
}

# Build a function for generating the a table with 
## we have 2 parameter one is the kind we want for rows and the other is the dataframe we are sing to generate the table 
##remember to put the "" mark, it only has 2 choices one is fraction and the other one is analytical samples, each kind will generate a pre-determined table
# when kind ==  fraction meaning we are geneating plots whose rows are fractions, like the table of 0 for precursor intensity 
# when kind == analytical sample the rows will be sample names like the table for the percentage of samples under 20 in MS2 per MS2 value

percentage_table <- function(kind, n)
{
  
  B=c(5,25,50,75,95)
  
  # generate a dataframe
  num_row <-  as.numeric(count(unique(select(n, contains(kind)))))
  
  if(kind=="fraction")
  {
    perc_table = matrix(,nrow = length(fractions), ncol = length(B))
  } else 
  {
    if(kind=="analytical")
    {
      perc_table = matrix(,nrow = length(analyticalSamples), ncol = length(B))
    }
  }
  
  if(kind == "fraction"){
    for (i in 1:nrow(perc_table)){
      for (j in 1:length(B)){
        if(sum(n$Fraction == i) > 0)
        {
          perc_table[i,j]<-
            round(length(which(n$Percentile == B[j] & 
                                 n$Fraction == i))/
                    length(which(QCmetricsLongPrecInt_ALL$Fraction == i &
                                   QCmetricsLongPrecInt_ALL$Percentile == B[j]))*100, 
                  3)
        }
        else
        {
          perc_table[i,j]<- 0
        }
        
      }
    }
    colnames(perc_table) <- c("5th Perc", "25th Perc", "50th Perc", 
                              "75th Perc","95th Perc")
    rownames(perc_table) <- c(paste0("fraction", 1:24))##c(1:25)))
    kable(perc_table)# %>% kable_styling(bootstrap_options = "striped", full_width = F)
    
  }else 
  {
    if(kind == "analytical"){
      for (i in 1:length(analyticalSamples)){
        for(j in 1:length(B)){
          perc_table[i,j] <- 
            round(length(which(n$Percentile == B[j] & 
                                 n$value <20 & 
                                 n$analyticalBasename == analyticalSamples[i] ))/ 
                    length(which(n$Percentile == B[j] & 
                                   n$analyticalBasename == analyticalSamples[i]))*100,
                  3)
        }
      }
      colnames(perc_table) <- c(paste0(B, rep("th Perc", 5)))
      rownames(perc_table) <- c(as.character(analyticalSamples))
      kable(perc_table, caption = "Percent of samples' value less than 20") #%>%
      # kable_styling(bootstrap_options = "striped") 
    }
  }
}

# Functions to make specific plots (Jiahao)
gg_points <- function(data, y_string, x_string,
                      col_string, label_string, label2_string="",
                      min_y=NULL, max_y=NULL, 
                      sparse_x_labels=0)
{
  if (label2_string == "") {
    label2_string = y_string;
  }
  data.original.level=as.vector(unique(data[, x_string]))
  break.points=as.numeric(unique(data[, x_string]))
  new.break.points=break.points
  if(length(break.points)>(sparse_x_labels+1) && sparse_x_labels>0){
    divisor=floor(length(break.points)/sparse_x_labels)
    max.x=divisor*sparse_x_labels
    new.break.points=break.points[seq(from=1,to=max.x, by=divisor)]
    if(new.break.points[length(new.break.points)]<max(break.points)) 
      new.break.points=c(new.break.points, max(break.points))
  }
  
  b <- ggplot(data, aes_string(y=y_string, 
                               x=x_string, 
                               col=col_string, 
                               label=label_string, 
                               label2=label2_string)) +
    geom_point() +
    scale_x_discrete(breaks = data.original.level[new.break.points], 
                     labels = data.original.level[new.break.points]) + 
    theme(# axis.title.x=element_blank(),
      plot.margin = margin(0.2, 0.5, 0.2, 0, "cm"),
      # axis.text.x=element_blank(),
      legend.position="none")
  if (!is.null(min_y) && !is.null(max_y)) {
    b <- b + ylim(min_y,max_y)
  } else if (!is.null(min_y)) {
    b <- b + ylim(min_y,NA)
  } else if (!is.null(max_y)) {
    b <- b + ylim(NA,max_y)
  }
  b
}

gg_points_plotly <- function(df, y_string, title, as_title, min_y=0)
{
  gg <- gg_points(df, 
                  y_string="ProteinCounts", x_string="analyticalBasename",
                  col_string="analyticalBasename", 
                  label_string="`Analytical sample`",min_y=min_y) 
  pg <- ggplotly(gg, tooltip=c("label","label2")) %>% 
    layout(annotations=title) %>% layout(annotations=as_title)
  
  subplot(pg, margin=0.05)##, titleX=TRUE, titleY=TRUE)
}

gg_box <- function(data, y_string, x_string,
                   col_string, label_string, label2_string="Fraction",
                   same_scale=FALSE, min_y=NULL, max_y=NULL,sparse_x_labels=0)
{
  data.original.level=as.vector(unique(data[, x_string]))
  break.points=as.numeric(unique(data[, x_string]))
  new.break.points=break.points
  if(length(break.points)>(sparse_x_labels+1) && sparse_x_labels>0){
    divisor=floor(length(break.points)/sparse_x_labels)
    max.x=divisor*sparse_x_labels
    new.break.points=break.points[seq(from=1,to=max.x, by=divisor)]
    if(new.break.points[length(new.break.points)]<max(break.points)) 
      new.break.points=c(new.break.points, max(break.points))
  }
  
  b <- ggplot(data, aes_string(y=y_string, x=x_string, 
                               col=col_string, label=label_string,
                               label2=label2_string)) +
    geom_boxplot()+
    geom_point(alpha=0.6) +
    scale_x_discrete(breaks = data.original.level[new.break.points], 
                     labels = data.original.level[new.break.points]) + 
    theme(axis.title.x=element_blank(),
          plot.margin = margin(0.2, 1.0, 0.2, 0, "cm"),
          # axis.text.x=element_blank(),
          legend.position="none")
  if(same_scale)
  {
    b <- b+ylim(min_y, max_y)
  } else if (!is.null(min_y) && !is.null(max_y)) {
    b <- b + ylim(min_y,max_y)
  } else if (!is.null(min_y)) {
    b <- b + ylim(min_y,NA)
  } else if (!is.null(max_y)) {
    b <- b + ylim(NA,max_y)
  }
  b
}

gg_points_frac <- function(data, y_string, x_string,
                           col_string, label_string, label2_string="Fraction",
                           smooth=FALSE,boxplot=FALSE,                   
                           same_scale=FALSE, min_y=NULL, max_y=NULL, sparse_x_labels=0)
{
  x_string_fact <- paste("as.factor(",x_string,")",sep="")
  x_string_used = ""
  if(boxplot)
  {
    x_string_used <- x_string_fact
  } else {
    x_string_used <- x_string
  }
  
  
  fractions.unique.orig=unique(data[,x_string])
  fractions.relevel=factor(data[,x_string], levels=ordered(fractions.unique.orig))
  data[,x_string]=fractions.relevel
  
  data.original.level=as.vector(unique(data[, x_string]))
  data[, x_string] = as.numeric(data[, x_string])
  
  
  break.points=as.numeric(unique(data[, x_string]))
  new.break.points=break.points
  if(length(break.points)>(sparse_x_labels+1) && sparse_x_labels>0){
    divisor=floor(length(break.points)/sparse_x_labels)
    max.x=divisor*sparse_x_labels
    new.break.points=break.points[seq(from=1,to=max.x, by=divisor)]
    if(new.break.points[length(new.break.points)]<max(break.points)) 
      new.break.points=c(new.break.points, max(break.points))
  }
}

gg_points_frac_perc <- function(data, y_string, x_string,
                                col_string, group_string,
                                label_string, 
                                label2_string="Fraction",
                                boxplot=FALSE,
                                same_scale=TRUE, min_y=NULL, max_y=NULL,sparse_x_labels=0)
{
  
  
  fractions.unique.orig=unique(data[,x_string])
  fractions.relevel=factor(data[,x_string], levels=ordered(fractions.unique.orig))
  data[,x_string]=fractions.relevel
  
  data.original.level=as.vector(unique(data[, x_string]))
  data[, x_string] = as.numeric(data[, x_string])
  
  
  break.points=as.numeric(unique(data[, x_string]))
  new.break.points=break.points
  if(length(break.points)>(sparse_x_labels+1)&&sparse_x_labels>0){
    divisor=floor(length(break.points)/sparse_x_labels)
    max.x=divisor*sparse_x_labels
    new.break.points=break.points[seq(from=1,to=max.x, by=divisor)]
    if(new.break.points[length(new.break.points)]<max(break.points)) new.break.points=c(new.break.points, max(break.points))
  }
  
  
  
  pfc <- 
    ggplot(data, 
           aes_string(y=y_string, x=x_string, 
                      col=col_string, group=group_string,
                      label=label_string, label2=label2_string)) + 
    facet_grid(~ Percentile)   +
    scale_x_continuous(breaks = new.break.points, labels=data.original.level[new.break.points]) +
    theme(plot.margin = margin(0.2, 1.0, 0.2, 0, "cm")) #change 
  
  
  if(same_scale){
    pfc <- pfc +  ylim(min_y,max_y)
  }  
  if(boxplot){
    pfc <- pfc + geom_boxplot() +
      theme(axis.title.x=element_blank(),
            #      axis.text.x=element_blank(),
            legend.position="none")
  } else {
    pfc <- pfc +
      theme(legend.position = "none")
  }
  pfc + geom_point(aes_string(col=col_string), alpha=0.6, size=0.9)
}

gg_pfp_plotly_all_id <- function(qc_all, qc_id, y_string,
                                 min_y_all, max_y_all,
                                 min_y_id, min_x_id,
                                 annot_all_1, annot_all_2,
                                 annot_id_1, annot_id_2)
{
  g_all <- gg_points_frac_perc(filter(qc_all, 
                                      Percentile %in% c(5,25,50,75,95)),
                               y_string = y_string,
                               x_string = "Fraction",
                               col_string = "analyticalBasename",
                               group_string="analyticalBasename",
                               label_string="`Analytical sample`",
                               min_y=min_y_all,
                               max_y=max_y_all,
                               sparse_x_labels=3)
  
  p_all <- ggplotly(g_all, tooltip=c("label","label2")) %>%
    layout(annotations=annot_all_1)  %>% layout(annotations=annot_all_2)
  
  # Or, we might not need to change the variable name
  g_id <- gg_points_frac_perc(filter(qc_id, 
                                     Percentile %in% c(5,25,50,75,95)),
                              y_string = y_string,
                              x_string = "Fraction",
                              col_string = "analyticalBasename",
                              group_string="analyticalBasename",
                              label_string="`Analytical sample`",
                              min_y=min_y_id,
                              max_y=max_y_id,
                              sparse_x_labels=3)
  p_id <- ggplotly(g_id, tooltip=c("label","label2")) %>%
    layout(annotations=annot_id_1) %>% layout(annotations=annot_id_2)
  
  subplot(p_all, p_id, titleY=TRUE, margin=0.05) #change
  
}

# Add ggbars
gg_bars <- function(data, y_string, x_string, col_string, label_string)
{
  ggplot(data, aes_string(y=y_string, x=x_string, fill=col_string, label=label_string)) +
    geom_bar(stat = "identity") +
    theme(axis.title.x=element_blank(),
          plot.margin = margin(0.2, 0.5, 0.2, 0, "cm"),
          # axis.text.x=element_blank(),
          legend.title = element_blank()
    )
}
  
# Function to melt variables by fraction
melt_by_fraction <- function(qc_df, variable.name)
{  
  qc_df_subset <- qc_df[,c(2,3,grep(paste(variable.name,".",sep=""),
                                    colnames(qc_df)))]
  qc_long <- melt(qc_df_subset,
                  variable.name=variable.name,
                  value.name="value",
                  id.vars=c("analyticalBasename","Fraction"))
  
  qc_long <- cbind(qc_long[,1:2],
                   colsplit(qc_long[,3], "\\.",
                            c(variable.name,"Percentile")),
                   QCmetricsLongPrecInt_ALL[,4])
  
  colnames(qc_long)[5] <- "value"
  
  qc_long$`Analytical sample` <- qc_long$analyticalBasename
  
  qc_long
}
  
