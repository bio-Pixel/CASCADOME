meta <- readRDS("/Users/panxu/Desktop/Project/MacroEnv/Merge_CSS_osp_530581_metadata.rds")
meta$Ann_Level0 = as.factor(meta$Ann_Level0)
meta$Ann_Level1 = as.factor(meta$Ann_Level1)
ml = split(meta, meta$batch)
dat = do.call(rbind, lapply(ml, function(dat){
  prop = table(dat$Ann_Level0)/dim(dat)[1]
  prop[levels(meta$Ann_Level0)]
}))
colnames(dat) = levels(meta$Ann_Level0)
dat = data.frame(state = gsub(".*_","",rownames(dat)),
                 organ = gsub("_Healthy|_Pancreatitis|_PDAC|_Cachexia","",rownames(dat)),
                 dat)
dat$state = factor(dat$state,
                     levels = c("Healthy","Pancreatitis","PDAC","Cachexia"))

library(ggplot2)

dat2 <- as.data.frame(dat)

stage_order <- c(
  "Healthy",
  "Pancreatitis",
  "PDAC",
  "Cachexia"
)

celltypes <- setdiff(
  colnames(dat2),
  c("state", "organ")
)

# numeric
dat2[celltypes] <- lapply(
  dat2[celltypes],
  function(x) as.numeric(as.character(x))
)

score_df <- lapply(celltypes, function(ct){
  
  pri_vec <- sapply(unique(dat2$organ), function(org){
    
    sub <- dat2[dat2$organ == org, ]
    
    vals <- sapply(stage_order, function(st){
      
      mean(
        sub[sub$state == st, ct],
        na.rm = TRUE
      )
      
    })
    
    H = vals[1]
    P = vals[2]
    D = vals[3]
    C = vals[4]
    
    d1 = P - H
    d2 = D - P
    d3 = C - D
    
    delta = c(d1,d2,d3)
    
    # 主方向
    dir = sign(sum(delta))
    
    consistency =
      mean(sign(delta) == dir)
    
    pri =
      sum(delta) * consistency
    
    pri
    
  })
  
  data.frame(
    celltype = ct,
    score = mean(pri_vec, na.rm = TRUE),
    consistency = mean(pri_vec > 0, na.rm = TRUE)
  )
  
})

score_df <- do.call(rbind, score_df)

score_df <- score_df[
  order(score_df$score, decreasing = TRUE),
]

score_df$celltype <- factor(
  score_df$celltype,
  levels = score_df$celltype
)

ggplot(score_df,
       aes(x = score, y = celltype , color = score)) +
  
  geom_point(size = 4) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  
  theme_classic(base_size = 14)

library(ggplot2)
library(ggrepel)

plot_df <- score_df
plot_df <- plot_df[order(plot_df$score, decreasing = TRUE), ]
plot_df$rank <- seq_len(nrow(plot_df))

top10 <- head(plot_df, 10)
bottom10 <- tail(plot_df, 10)

top10$type <- "Top10"
bottom10$type <- "Bottom10"

label_df <- rbind(top10, bottom10)

ggplot(plot_df, aes(x = rank, y = score)) + 
  geom_area(fill = "grey80", alpha = 0.75) + 
  geom_point( color = "grey45", size = 1.8, alpha = 0.75 ) + 
  geom_point( data = top10, aes(x = rank, y = score), color = "#B22222", size = 1.8 ) + 
  geom_point( data = bottom10, aes(x = rank, y = score), color = "#2166AC", size = 1.8 ) +    
  geom_text_repel(
    data = label_df,
    aes(label = celltype),
    size = 3.8, 
    segment.color = "grey60",
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    nudge_x = ifelse(label_df$type == "Top10", 6, -6),
    direction = "y"
  ) +
  geom_hline( yintercept = 0, linetype = 2, color = "grey50" ) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0(nrow(plot_df), " cell types"),
    y = "Systemic progression score"
  ) + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(aspect.ratio = 1,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave("../../CellPropDiff.pdf", height = 4, width = 4)
