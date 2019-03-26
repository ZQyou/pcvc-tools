library(tidyr)
library(ggplot2)

num_panel = 4
df <- read.csv("test.csv")
df <- df[1:num_panel,]

# round up the values
df[,2:ncol(df)] <- round(df[,2:ncol(df)],2)
# transform wide data to long
df_long <- gather(df, contig, value, M1:Y2B, factor_key=TRUE)
# plot
p <- ggplot(data = df_long[which(df_long$value>0.05),], aes(x = "", y = value, fill = contig)) + 
     geom_bar(stat = "identity", position = position_fill()) +
     geom_text(aes(label = value), position = position_fill(vjust = 0.5)) +
     coord_polar(theta = "y") +
     facet_grid(facets=. ~ subclusters) + 
     theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
           axis.ticks = element_blank(), axis.text = element_blank(), panel.grid = element_blank(), # panel control
           panel.background = element_blank()) + 
     theme(strip.background = element_blank()) + 	# facet control
     theme(legend.position='bottom') +
     guides(fill=guide_legend(nrow=1,byrow=TRUE)) + 
     scale_fill_manual(values=c("#99ffcc","#66ffcc","#00ffcc","#009999","#006666","#66ff66","#00cc33",
                                "#009900","#006600","#ff6666","#ff00cc","#666666","#33ffcc","#99ff99",
                                "#00ff66","#669900","#666600","#ffcc33","#ff9933","#ff6666","#ff00cc"))

p 

# 
# Below are the RGB color codes for samples. 
#
# M1 (153, 255, 204)		#99ffcc
# M15 (102, 255, 204)           #66ffcc
# M2 (0, 255, 204)              #00ffcc
# M25 (0, 153, 153)             #009999
# M35 (0, 102, 102)             #006666
# M5 (102, 255, 102)            #66ff66
# M6 (0, 204, 51)               #00cc33
# M7 (0, 153, 0)                #009900
# M75 (0, 102, 0)               #006600
# Y1 (255, 102, 102)            #ff6666
# Y2 (255, 0, 204)              #ff00cc
# RefSeq (102, 102, 102)        #666666
# M1B (51, 255, 204)            #33ffcc
# M4B (153, 255,153)            #99ff99
# M5B (0, 255, 102)             #00ff66
# M6B (102, 153, 0)             #669900
# M7B (102, 102, 0)             #666600
# SB1 (255, 204, 51)            #ffcc33
# SB2 (255, 153, 51)            #ff9933
# Y1B (255, 102, 102)           #ff6666
# Y2B (255, 0, 204)             #ff00cc
#
