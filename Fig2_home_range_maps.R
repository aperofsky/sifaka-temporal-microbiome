# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

### this script is sourced within the PCOA script to make Figure 2B
library(dplyr)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggpubfigs) #devtools::install_github("JLSteenwyk/ggpubfigs")
library(cowplot)
####################################################################
## Home range maps: Groups I to VI
####################################################################
load("Rdata/home_range_gps_df_2012.Rdata")
g_2011_2012 <- ggplot(df1_2012, aes(x = long, y = lat, fill = Group, group = group)) +
  geom_polygon(alpha = 0.7) +
  scale_fill_manual(values=friendly_pal("bright_seven"))+
  coord_equal() + 
  theme_bw(base_size=15) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(fill="Group")
g_2011_2012

### 2014-2015
load("Rdata/home_range_gps_df_2015.Rdata")
head(df1_2015)
g_2014_2015 <- ggplot(df1_2015, aes(x = long, y = lat, fill = Group, group = group)) +
  geom_polygon(alpha = 0.7) +
  scale_fill_manual(values=friendly_pal("bright_seven"))+
  coord_equal() + 
  theme_bw(base_size=15) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(fill="Group")
g_2014_2015

### 2015-2016
load("Rdata/home_range_gps_df_2016.Rdata")
g_2015_2016 <- ggplot(df1_2016, aes(x = long, y = lat, fill = Group, group = group)) +
  geom_polygon(alpha = 0.7) +
  scale_fill_manual(values=friendly_pal("bright_seven"))+
  coord_equal() + 
  theme_bw(base_size=15) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(fill="Group")
g_2015_2016

plot_grid(g_2011_2012,g_2014_2015,g_2015_2016,nrow=1)

prow <- plot_grid(
  g_2011_2012+ theme(legend.position="none"),
  g_2014_2015 + theme(legend.position="none"),
  g_2015_2016 + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)
prow
legend <- get_legend(
  # create some space to the left of the legend
  g_2011_2012 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

multi = plot_grid(prow, legend, rel_widths = c(3, .4))
multi

labels_12 = df1_2012 %>%
  dplyr::group_by(Group)%>%
  dplyr::summarise(long=mean(long),
                   lat = mean(lat))%>%
  ungroup()%>%
  mutate()%>%
  as.data.frame()

g_2011_2012_full = g_2011_2012 + scale_x_continuous(labels = scales::number_format(accuracy = 0.001,
                                                                                   decimal.mark = '.'))+ theme(legend.position = "none")
year1 = g_2011_2012_full+ geom_text(data=labels_12, aes( x=long, y=lat,group=Group,label=Group), size=4)

labels_15 = df1_2015 %>%
  dplyr::group_by(Group)%>%
  dplyr::summarise(long=mean(long),
                   lat = mean(lat))%>%
  ungroup()%>%
  mutate()%>%
  as.data.frame()

g_2014_2015_full = g_2014_2015 + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.001,
                                 decimal.mark = '.'))+  
  scale_y_continuous(breaks=c(-20.787,-20.784,-20.781,-20.778))+
  theme(legend.position = "none")

year2 = g_2014_2015_full + geom_text(data=labels_15, aes( x=long, y=lat,group=Group,label=Group), size=4)

labels_16 = df1_2016 %>%
  dplyr::group_by(Group)%>%
  dplyr::summarise(long=mean(long),
                   lat = mean(lat))%>%
  ungroup()%>%
  mutate()%>%
  as.data.frame()

g_2015_2016_full = g_2015_2016 + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.001,
                                 decimal.mark = '.'))+ theme(legend.position = "none")

year3 = g_2015_2016_full + geom_text(data=labels_16, aes( x=long, y=lat,group=Group,label=Group), size=4)

gps = plot_grid(year1,year2,year3,hjust = -1,nrow=1)
gps
group_fig = plot_grid(p+theme(legend.position = "none"),gps,nrow=2,labels=c("a","b"))## p is pcoa figure 2

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  p + theme(legend.box.margin = ggplot2::margin(0, 0, 0, 12))
)

multi = plot_grid(group_fig, legend, rel_widths = c(3, .4))
multi
