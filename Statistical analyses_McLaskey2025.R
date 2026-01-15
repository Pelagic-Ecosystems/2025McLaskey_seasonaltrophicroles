
require(tidyverse)
require(vegan)
require(Polychrome)
require(gmodels)
require(adespatial)
require(ggfortify)
require(clustsig)
require(cluster)
require(labdsv)
require(dendextend)
require(cowplot)
library(car)



# Load full FA dataset 
zoopData <- read.csv("processed_data/zoopsprintData all 2025-03-17.csv")
str(zoopData)
zoopData$Date <- as.Date(zoopData$Date, "%Y-%m-%d")

clrs <- c( "#55405C", "#FD1600", "#1C22FE", "#35FE1C", "#FE00D2", "#FEA60D", "#901C45",
                    "#00BBFC", "#1CFFCD", "#D1E265")
POMclrs <- c("#3EBA34", "#55405C", "#FD1600", "#1C22FE", "#35FE1C", "#FE00D2", "#FEA60D", "#901C45",
                                          "#00BBFC", "#1CFFCD", "#D1E265")
POMsizesclrs <- c("#c2e699","#006837", "#78c679", "#55405C", "#FD1600", "#1C22FE", "#35FE1C", "#FE00D2", "#FEA60D", "#901C45",
                                                                     "#00BBFC", "#1CFFCD", "#D1E265")



# POM size fractions
POM.sizes <- zoopData %>% filter(Type == "POM" & Pore.Size !="bulk")

# Bulk POM
bulk.POM <- zoopData %>% filter(Pore.Size =="bulk")

# Zooplankton species
zoopData.sm <- zoopData  %>% filter(Type=="Zoop") %>%
  filter(Species!="250")


taxa.summary.small <- zoopData.sm %>% select(Species, SumFA_mg.g, DHA.EPA,C20.4n.6_PERCENT)
zoopData.sm.grouped <- taxa.summary.small %>% group_by(Species)
proportions.summary <- zoopData.sm.grouped %>% dplyr::summarise(across(SumFA_mg.g:C20.4n.6_PERCENT, ~ mean(.x, na.rm = TRUE)))

colnames(proportions.summary)[2:5] <- paste0("Mean.", colnames(proportions.summary)[2:5])



zoopData.sm.means <- zoopData.sm %>% dplyr::select(Species, Month, SumFA_mg.g, 
                                                   DHA.EPA, percent.SFA, percent.PUFA, percent.MUFA,
                                                   C18.3andC18.4, C16PUFA, Bacteria_no15, Copes.C20C.1C22.1,
                                                   Diatom.Flag, Ratio16.1, C16.3n.4_PERCENT, 
                                                   Date) %>% 
  group_by(Species, Month) %>% 
  summarise_all( mean)


colnames(zoopData.sm.means)[3:14] <- paste0("mean.", colnames(zoopData.sm.means)[3:14])

zoopData.sm.wMeans <- zoopData.sm.means %>% dplyr::select(-Date) %>% 
  full_join(.,zoopData.sm, by = join_by(Species,Month))



# *Fig. S5C ----------------------------------------------------------------


TFA.plot <- ggplot(zoopData.sm.wMeans, aes(x=Date, y=mean.SumFA_mg.g, color=Species, linetype=Species)) + 
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank())  +  
  scale_y_continuous(expand=c(0,0), limits = c(0,235)) + 
  geom_point(aes(x=Date, y=SumFA_mg.g)) + 
  geom_line(size=1.2) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")) +
  scale_color_manual(values = clrs) + 
  labs(x=element_blank(), y=expression(paste("Total FA (mg g"^-1*" DW)")))  + 
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May\n2021", "Jul\n2021", "Sep\n2021", 
                          "Nov\n2021", "Jan\n2022", "Mar\n2022")) 
  TFA.plot

# tiff("TFA timeseries 20250717.tiff", width=100, height=70, units="mm", res=500)
# TFA.plot
# dev.off()




# *Fig. S6 ----------------------------------------------------------------


zoopsprint.bulk.means <- zoopData.sm %>% dplyr::select(Species, Month, SumFA_mg.g,
                                                   DHA.EPA, percent.SFA, percent.PUFA, percent.MUFA,
                                                   C18.3andC18.4, C16PUFA, Bacteria_no15, Copes.C20C.1C22.1,
                                                   Diatom.Flag, Ratio16.1, C16.3n.4_PERCENT, 
                                                   Date) %>% 
  group_by(Species, Month) %>% 
  summarise_all( mean)

colnames(zoopsprint.bulk.means)[3:14] <- paste0("mean.", colnames(zoopsprint.bulk.means)[3:14])

zoopsprint.bulk.wMeans <- zoopsprint.bulk.means %>% dplyr::select(-Date) %>% 
  full_join(.,zoopData.sm, by = join_by(Species,Month))


# POM size means independently
POM.sizes$Species <- POM.sizes$Pore.Size

zoopsprint.POM.means <- POM.sizes %>% dplyr::select(Species, Month, SumFA_mg.g,
                                                       DHA.EPA, percent.SFA, percent.PUFA, percent.MUFA,
                                                       C18.3andC18.4, C16PUFA, Bacteria_no15, Copes.C20C.1C22.1,
                                                       Diatom.Flag, Ratio16.1, C16.3n.4_PERCENT, 
                                                       Date) %>% 
  group_by(Species, Month) %>% 
  summarise_all( funs(mean), na.rm = T )

colnames(zoopsprint.POM.means)[3:14] <- paste0("mean.", colnames(zoopsprint.POM.means)[3:14])

zoopsprint.POM.wMeans <- zoopsprint.POM.means %>% dplyr::select(-Date) %>% 
  full_join(.,POM.sizes, by = join_by(Species,Month))
zoopsprint.POM.wMeans<- zoopsprint.POM.wMeans %>%  mutate(Species = factor(Species, 
                                                                           levels = c("0.7", "3","20")))


pDHA.EPA <- ggplot(zoopsprint.bulk.wMeans, aes(x=Date, y=mean.DHA.EPA, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=DHA.EPA)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")  +  
  geom_line(aes(linetype = Species), size=1.2) +
  # scale_y_continuous(limits = c(0, 2.5), expand = c(0,0)) + 
  scale_linetype_manual(values = c( "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")) + 
  scale_color_manual(values = clrs) +
  scale_fill_manual(values =  clrs)  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs(x="", y="DHA:EPA") 
pDHA.EPA
# tiff("DHA.EPA 20250717.tiff", width=110, height=70, units="mm", res=500)
# pDHA.EPA
# dev.off() 

pDHA.EPA2 <- ggplot(zoopsprint.POM.wMeans, aes(x=Date, y=mean.DHA.EPA, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=DHA.EPA)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(),
        axis.title.x = element_blank())  +  
  geom_line(linetype = "twodash", size=1.2) +
  # scale_y_continuous(limits = c(0, 2.5), expand = c(0,0)) + 
  scale_color_manual(values = c("#b8e166", "#0073e6","#eb59a1")) +
  scale_fill_manual(values =  c("#b8e166", "#0073e6","#eb59a1"))  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs(y="DHA:EPA",
       color="Species" , fill="Species") 
pDHA.EPA2

prow1 <- plot_grid(pDHA.EPA2, 
                  pDHA.EPA,
                  rel_heights = c(45, 70),
                  nrow = 2,
                  align = "v")
prow1
# tiff("DHA.EPA Zoop w POM sizes 20250717.tiff", width=120, height=120, units="mm", res=500)
# prow1
# dev.off()



pSFA <- ggplot(zoopsprint.bulk.wMeans, aes(x=Date, y=mean.percent.SFA*100, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=percent.SFA*100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")  +  
  geom_line(aes(linetype = Species), size=1.2) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values =  clrs)  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs(x="", y="% SFA",
       color="Species" , fill="Species") 
pSFA
# tiff("SFA 20250717.tiff", width=110, height=70, units="mm", res=500)
# pSFA
# dev.off()

pSFA2 <- ggplot(zoopsprint.POM.wMeans, aes(x=Date, y=mean.percent.SFA*100, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=percent.SFA*100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(),
        axis.title.x = element_blank())  +  
  geom_line(linetype = "twodash", size=1.2) +
  scale_color_manual(values = c("#b8e166", "#0073e6","#eb59a1")) + 
  scale_fill_manual(values =  c("#b8e166", "#0073e6","#eb59a1"))  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs( y="% SFA",
       color="Species" , fill="Species") 
pSFA2
prow2 <- plot_grid(pSFA2, 
                  pSFA,
                  rel_heights = c(45, 70),
                  nrow = 2,
                  align = "v")
prow2
# tiff("SFA Zoop w POM 20250717.tiff", width=120, height=120, units="mm", res=500)
# prow2
# dev.off()


pPUFA <- ggplot(zoopsprint.bulk.wMeans, aes(x=Date, y=mean.percent.PUFA*100, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=percent.PUFA*100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")  +  
  geom_line(aes(linetype = Species), size=1.2) +
  scale_linetype_manual(values = c( "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values =  clrs)  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs(y="% PUFA",
       color="Species" , fill="Species") 

pPUFA
# tiff("PUFA 20250424.tiff", width=110, height=70, units="mm", res=500)
# pPUFA
# dev.off()

pPUFA2 <- ggplot(zoopsprint.POM.wMeans, aes(x=Date, y=mean.percent.PUFA*100, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=percent.PUFA*100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(),
        axis.title.x = element_blank())  +  
  geom_line(linetype = "twodash", size=1.2) +
  scale_color_manual(values = c("#b8e166", "#0073e6","#eb59a1")) + 
  scale_fill_manual(values =  c("#b8e166", "#0073e6","#eb59a1"))  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs(y="% PUFA",
       color="Species" , fill="Species") 
pPUFA2
prow3 <- plot_grid(pPUFA2, 
                  pPUFA,
                  rel_heights = c(45, 70),
                  nrow = 2,
                  align = "v")
prow3
# tiff("PUFA Zoop w POM 20250717.tiff", width=120, height=120, units="mm", res=500)
# prow3
# dev.off()


pBAFA <- ggplot(zoopsprint.bulk.wMeans, aes(x=Date, y=mean.Bacteria_no15*100, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=Bacteria_no15*100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")  +  
  geom_line(aes(linetype = Species), size=1.2) +
  scale_linetype_manual(values = c( "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values =  clrs)  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs(x="", y="% Bacterial FA",
       color="Species" , fill="Species") 

pBAFA
# tiff("Bacterial FA 20250717.tiff", width=110, height=70, units="mm", res=500)
# pBAFA
# dev.off()

pBAFA2 <- ggplot(zoopsprint.POM.wMeans, aes(x=Date, y=mean.Bacteria_no15*100, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=Bacteria_no15*100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(),
        axis.title.x = element_blank())  +  
  geom_line(linetype = "twodash", size=1.2) +
  scale_color_manual(values = c("#b8e166", "#0073e6","#eb59a1")) + 
  scale_fill_manual(values =  c("#b8e166", "#0073e6","#eb59a1"))  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs( y="% Bacterial FA",
       color="Species" , fill="Species") 
pBAFA2
prow4 <- plot_grid(pBAFA2, 
                  pBAFA,
                  rel_heights = c(45, 70),
                  nrow = 2,
                  align = "v")
prow4
# tiff("Bacterial FA Zoops w POM sizes 20250717.tiff", width=120, height=120, units="mm", res=500)
# prow4
# dev.off()


pC16PUFA <- ggplot(zoopsprint.bulk.wMeans, aes(x=Date, y=mean.C16PUFA*100, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=C16PUFA*100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")  +  
  geom_line(aes(linetype = Species), size=1.2) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values =  clrs)  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May\n2021", "Jul\n2021", "Sep\n2021", 
                          "Nov\n2021", "Jan\n2022", "Mar\n2022")) + 
  labs(x="", y="% C16-PUFAs") 

pC16PUFA
# tiff("C16PUFA 20250717.tiff", width=110, height=70, units="mm", res=500)
# pC16PUFA
# dev.off()

pC16PUFA2 <- ggplot(zoopsprint.POM.wMeans, aes(x=Date, y=mean.C16PUFA*100, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=C16PUFA*100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(),
        axis.title.x = element_blank())  +  
  geom_line(linetype = "twodash", size=1.2) +
  scale_color_manual(values = c("#b8e166", "#0073e6","#eb59a1")) + 
  scale_fill_manual(values =  c("#b8e166", "#0073e6","#eb59a1"))  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs( y="% C16-PUFAs",
       color="Species" , fill="Species") 
pC16PUFA2
prow5 <- plot_grid(pC16PUFA2, 
                  pC16PUFA,
                  rel_heights = c(45, 70),
                  nrow = 2,
                  align = "v")
prow5
# tiff("C16PUFA Zoop w POM sizes 20250717.tiff", width=120, height=120, units="mm", res=500)
# prow5
# dev.off()


pC18PUFA <- ggplot(zoopsprint.bulk.wMeans, aes(x=Date, y=mean.C18.3andC18.4*100, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=C18.3andC18.4*100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")  +  
  geom_line(aes(linetype = Species), size=1.2) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values =  clrs)  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May\n2021", "Jul\n2021", "Sep\n2021", 
                          "Nov\n2021", "Jan\n2022", "Mar\n2022")) + 
  labs(x="", y="% C18-PUFAs")  
pC18PUFA 
# tiff("C18PUFA 20250717.tiff", width=110, height=70, units="mm", res=500)
# pC18PUFA
# dev.off()

pC18PUFA2 <- ggplot(zoopsprint.POM.wMeans, aes(x=Date, y=mean.C18.3andC18.4*100, color=Species, fill=Species)) + 
  geom_point(aes(x=Date, y=C18.3andC18.4*100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(),
        axis.title.x = element_blank())  +  
  geom_line(linetype = "twodash", size=1.2) +
  scale_color_manual(values = c("#b8e166", "#0073e6","#eb59a1")) + 
  scale_fill_manual(values =  c("#b8e166", "#0073e6","#eb59a1"))  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs( y="% C18-PUFAs",
       color="Species" , fill="Species") 
pC18PUFA2
prow6 <- plot_grid(pC18PUFA2, 
                  pC18PUFA,
                  rel_heights = c(45, 70),
                  nrow = 2,
                  align = "v")
prow6
# tiff("C18PUFA Zoop w POM sizes 20250717.tiff", width=120, height=120, units="mm", res=500)
# prow6
# dev.off()





# Chlorophyll -------------------------------------------------------------


chl_all <- read.csv("processed_data/zoopsprint Chl data 2025-03-17.csv")
chl_all$Date <- as.Date(chl_all$Date, "%Y-%m-%d")

# Chl size class concentrations
chlSumm.concs.small <- chl_all %>% select(Date, chl_GF.F, chl_20um, chl_3um)
chlSumm.concs.small.long <- chlSumm.concs.small %>% pivot_longer(-Date, names_to = "Size.class", values_to = "Conc")
# chlSumm.concs.small.long$Size.class <- as.factor(chlSumm.concs.small.long$Size.class)
chlSumm.concs.small.long<- chlSumm.concs.small.long %>%  mutate(Size.class = factor(Size.class, 
                                            levels = c("chl_GF.F", "chl_3um","chl_20um")))



# Chlorophyll 
chlSumm.props.small <- chl_all %>%   mutate(chla_prop_GF.F = chl_GF.F / SumChl)  %>% 
  mutate(chla_prop_20um = chl_20um / SumChl)  %>% 
  mutate(chla_prop_3um = chl_3um / SumChl)  
chlSumm.props.small <- chlSumm.props.small %>% select(Date, chla_prop_GF.F, chla_prop_20um, chla_prop_3um)

chlSumm.props.small.long <- chlSumm.props.small %>% pivot_longer(-Date, names_to = "Size.class", values_to = "Prop")
chlSumm.props.small.long<- chlSumm.props.small.long %>%  mutate(Size.class = factor(Size.class, 
                                                                                    levels = c("chla_prop_GF.F", "chla_prop_3um","chla_prop_20um")))




  

# Fig. 5A -----------------------------------------------------------------


pChl.props <- ggplot(chlSumm.props.small.long, aes(Date, Prop, fill=Size.class)) + geom_area() + 
  theme_bw()+ 
  # scale_fill_manual(values=c("#c2e699","#78c679", "#006837"),
  scale_fill_manual(values = c("#b8e166", "#0073e6","#eb59a1"),
                    labels=c(expression(paste("0.7 ",mu,"m")), 
                             expression(paste("3 ",mu,"m")),
                             expression(paste("20 ",mu,"m")))) + 
  labs(y="Chl proportions", fill="Size Class")   + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 0.3),
        panel.grid = element_blank()) + 
  scale_x_date(limits = as.Date(c("2021-01-01", "2022-05-01")), expand = c(0,0),
               breaks = as.Date(c("2021-01-01", "2021-03-01",
                                  "2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("Jan","Mar", "May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  scale_y_continuous(expand = c(0,0))
pChl.props

# tiff("Chl a props 20250715.tiff", width=90, height=70, units="mm", res=500)
# pChl.props
# dev.off()


pChl.legend <- ggplot(chlSumm.concs.small.long, aes(Date, Conc, fill=Size.class)) + geom_area() + 
  theme_bw()+ 
  scale_fill_manual(values=c("#89ce00", "#0073e6","#e6308a"),
                    labels=c(expression(paste("0.7 ",mu,"m")), 
                             expression(paste("3 ",mu,"m")),
                             expression(paste("20 ",mu,"m")))) + 
  labs(y="Chl proportions", fill="Size Class")    


# tiff("Chl a lengend 20250703.tiff", width=100, height=60, units="mm", res=500)
# pChl.legend
# dev.off() 


# Overlay of Chl-a concentration for Fig. 5A 

ggp_transparent1<- ggplot(chl_all, aes(x=Date, y=SumChl)) + 
  geom_line(color = "black", size=1) + 
  labs(x="", color="Species" , fill="Species",
       y=expression(paste("Chl-a (",mu,"g L"^-1*")"))) + 
  theme(rect = element_rect(fill = "transparent"),
        legend.position= "none",
        panel.background = element_rect(fill = "transparent"),
        panel.grid = element_blank()) + 
  scale_y_continuous(expand=c(0,0)) + 
  scale_x_date(limits = as.Date(c("2021-01-01", "2022-05-01")), expand = c(0,0),
               breaks = as.Date(c("2021-01-01", "2021-03-01",
                                  "2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")))
               
ggp_transparent1

# ggsave(ggp_transparent1,            # Save transparent png file
#        filename = "transparent Sum Chl black.tiff",
#        bg = "transparent",
#        width=90, height=70, unit="mm", dpi=500)







# Stable isotopes -------------------------------------------------------


SI.all <- read.csv("processed_data/full_zoopsprint SI 20250912.csv")
SI.all$Date <- as.Date(SI.all$Date, "%Y-%m-%d")


QU39POM.ZSperiod <- SI.all %>% filter(Pore.Size == "bulk")


SI.POM <- SI.all %>% filter(Species=="POM") %>% 
  filter(Pore.Size != "bulk") %>%   select(Date, Month, Pore.Size, 
                                         C_N, delta13c, delta15n)


QU39POM.ZSperiod$Month <- format(as.Date(QU39POM.ZSperiod$Date), "%m")

SI.POM <- SI.POM %>% filter(Date != "2022-03-30")
SI.POM.sm <- SI.POM 
ZSdates <- c("2021-05-20", "2021-07-22", "2021-09-02", "2021-11-04", "2022-01-15", "2022-03-30")

SI.POM.sm$Date2 <- c(rep(ZSdates[1],9), 
                     rep(ZSdates[2],9),
                     rep(ZSdates[3],6),
                     rep(ZSdates[5],12),
                     rep(ZSdates[4],6))


SI.POM.sm$Date.num <- as.numeric(as.Date(SI.POM.sm$Date2))


SI.POM.sm$delta15n <- as.numeric((SI.POM.sm$delta15n))
SI.POM.sm$Date <- as.Date(SI.POM.sm$Date, "%Y-%m-%d")
SI.POM$Month <- format(as.Date(SI.POM$Date), "%m")




# *Table S2 and S3 --------------------------------------------------------




# # # # # 

SI.POM$Month <- factor(SI.POM$Month, levels = c( "09", "07", "05", "11", "01"))
SI.POM$Pore.Size <- factor(SI.POM$Pore.Size, levels = c( "3", "20", "0.7"))

mod <- lm(SI.POM$delta13c ~ SI.POM$Pore.Size + as.factor(SI.POM$Month))
summary(mod)
# July different from Sept, Nov, and Jan but not May
# Jan different from Sept, May, and July but not Nov
# May different from Nov, Sept, Jan but not July
# Sept different from all

plot(mod)

leveneTest(lm(delta13c ~ Pore.Size*Month, data=SI.POM))
# no evidence of unequal variances
leveneTest(lm(delta15n ~ Pore.Size*as.factor(SI.POM$Month), data=SI.POM))
# no evidence of unequal variances

del13.corr.AOV <- aov(SI.POM$delta13c ~ SI.POM$Pore.Size + as.factor(SI.POM$Month))
summary(del13.corr.AOV)
TukeyHSD(del13.corr.AOV)


del15AOV <- aov(SI.POM$delta15n ~ SI.POM$Pore.Size + as.factor(SI.POM$Month))
summary(del15AOV)
TukeyHSD(del15AOV)
# 20 is different from pico, but not other combos



# *Fig.2 A,B -------------------------------------------------------------------

SI.POM$Month <- as.character(SI.POM$Month)
SI.POM$Month[SI.POM$Month=="01"] <- "13"
SI.POM$Month[SI.POM$Month=="03"] <- "15"

SI.POM<- SI.POM %>%  mutate(Pore.Size = factor(Pore.Size, 
                                            levels = c("0.7", "3","20")))


sizes.C <- ggplot(SI.POM, aes(x=as.factor(Month), y=delta13c, fill=Pore.Size)) + 
  geom_boxplot(width = 0.5) + labs(x="") + 
  scale_y_continuous(limits=c(-30, -19)) +
  scale_fill_manual(values = c("#89ce00", "#0073e6","#e6308a")) +
  # scale_fill_manual(values=c("#c2e699","#78c679", "#006837")) + 
  labs(y=expression(paste(delta^{13}, "C (\u2030)")), x="")+ 
  theme_classic() + theme_bw() + 
  theme(legend.position = "none",
        panel.grid = element_blank())
sizes.C

sizes.N <- ggplot(SI.POM, aes(x=as.factor(Month), y=delta15n, fill=Pore.Size)) + 
  geom_boxplot(width = 0.5) + labs(x="") + 
  scale_y_continuous(limits=c(2.5, 11.5)) +
  scale_x_discrete(labels = c("May\n2021", "Jul\n2021", "Sep\n2021", "Nov\n2021", "Jan\n2022")) +
  scale_fill_manual(values = c("#89ce00", "#0073e6","#e6308a")) +
  # scale_fill_manual(values=c("#c2e699","#78c679", "#006837")) + 
  labs(y=expression(paste(delta^{15}, "N (\u2030)")), x="")+ 
  theme_classic() + theme_bw() + 
  theme(legend.position = "none",
        panel.grid = element_blank())
sizes.N

prow <- plot_grid(sizes.C, 
                  sizes.N,
                  align = "v",
                  # rel_heights = c(0.9, 1),
                  nrow = 2)
prow

# png("20250913 Figure 2 a,b.png", width=80, height=130, units="mm", res=300)
# prow
# dev.off()




# Fig. 5B,C ---------------------------------------------------------------



bulk.N <- ggplot(QU39POM.ZSperiod, aes(x=Date, y= delta15n)) + 
  geom_vline(xintercept = c(as.Date(c("2021-05-20", "2021-07-22", "2021-09-02", 
                                      "2021-11-04", "2022-01-15","2022-03-30"))),
             color = "grey60") + 
  geom_point(size=2) + theme_bw() + 
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(2.5, 11.5)) +
  scale_x_date(limits = as.Date(c("2021-01-01", "2022-05-01")), expand = c(0,0),
               breaks = as.Date(c("2021-01-01", "2021-03-01",
                                  "2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("Jan\n2021","Mar\n2021", "May\n2021", "Jul\n2021", "Sep\n2021", 
                          "Nov\n2021", "Jan\n2022", "Mar\n2022"))  +
  labs(y=expression(paste(delta^{15}, "N (\u2030)")), x="") 
bulk.N

bulk.C <- ggplot(QU39POM.ZSperiod, aes(x=Date, y= delta13c)) + 
  geom_vline(xintercept = c(as.Date(c("2021-05-20", "2021-07-22", "2021-09-02", 
                                      "2021-11-04", "2022-01-15","2022-03-30"))),
             color = "grey60") + 
  geom_point(size=2) + theme_bw() + 
  theme(legend.position = "none",
        panel.grid = element_blank()) + 
  scale_y_continuous(limits=c(-30, -19)) +
  geom_point(size=2) + 
  scale_x_date(limits = as.Date(c("2021-01-01", "2022-05-01")), expand = c(0,0),
               breaks = as.Date(c("2021-01-01", "2021-03-01",
                                  "2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("Jan","Mar", "May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar")) + 
  labs(y=expression(paste(delta^{13}, "C (\u2030)")), x="")  
bulk.C


prow <- plot_grid(bulk.C, 
                  bulk.N, 
                  align = "v",
                  # rel_heights = c(0.9, 1),
                  nrow = 2)
prow

# png("20250804 Figure 5 b,c.png", width=90, height=130, units="mm", res=300)
# prow
# dev.off()




# Chl : C ratio ----------------------------------------------------------

QU39POM.ZSperiod <- QU39POM.ZSperiod %>% select(-Depth)

chl_all$Sum.Chl <- chl_all$chl_GF.F + chl_all$chl_20um + chl_all$chl_3um

Chl.carbon <- full_join(chl_all, QU39POM.ZSperiod)

Chl.carbon$C.Chl.acid <- Chl.carbon$C_ug.L.acidified / Chl.carbon$Sum.Chl



# * Fig. S10 ---------------------------------------------------------------

pPOC.Chl <- ggplot(Chl.carbon, aes(x=Date, y=log10(C.Chl.acid))) + geom_point() +
  scale_x_date(limits = as.Date(c("2021-01-01", "2022-05-01")), expand = c(0,0),
               breaks = as.Date(c("2021-01-01", "2021-03-01",
                                  "2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("Jan","Mar", "May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar"))  +
  theme_bw() + geom_hline(yintercept = log10(200)) + 
  theme(panel.grid.minor = element_blank()) +
  labs(y="log10(POC:Chl)", x=element_blank()) 
pPOC.Chl

# tiff("POC.Chl ratio 20250218.tiff", width=100, height=70, units="mm", res=500)
# pPOC.Chl
# dev.off()


pPOC <- ggplot(Chl.carbon, aes(x=Date, y=C_ug.L.acidified)) + geom_point() +
  scale_x_date(limits = as.Date(c("2021-01-01", "2022-05-01")), expand = c(0,0),
               breaks = as.Date(c("2021-01-01", "2021-03-01",
                                  "2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("Jan","Mar", "May", "Jul", "Sep", 
                          "Nov", "Jan", "Mar"))  +
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  labs(x=element_blank(), y=expression(paste("POC (",mu,"g L"^-1*")")))  

pPOC

# tiff("POC timeseries.tiff", width=100, height=70, units="mm", res=500)
# pPOC
# dev.off()





# POM FA  -------------------------------------------------------------
POM.sizes$Month
POM.sizes$Month[POM.sizes$Month==1] <- 13

POM.sizes<- POM.sizes %>%  mutate(Pore.Size = factor(Pore.Size, 
                                            levels = c("0.7", "3","20")))
POM.sizes<- POM.sizes %>%  mutate(Month = factor(Month, 
                                                     levels = c("5", "7","9","11","13")))


# Figs. S1 + S2 -----------------------------------------------------------


            
# Sign. IndVals of micro-POM
ggplot(POM.sizes, aes(x=Pore.Size, y=C20.5n.3_PERCENT, fill = as.factor(Month))) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = c( "#A6D96A",  "#FDAE61", 
                                         "#D53E4F",  "#67001F", "#053061"),
                                         labels = c("May 2021", "July", "Sept", "Nov", "Jan 2022")) +
  scale_x_discrete(labels = c("pico-POM", "nano-POM", "micro-POM")) + 
  labs(x=element_blank(), y = expression(paste("20:5",omega,"3")),
       fill = element_blank()) + 
  theme_bw()

ggplot(POM.sizes, aes(x=Pore.Size, y=C22.6n.3_PERCENT, fill = as.factor(Month))) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = c( "#A6D96A",  "#FDAE61", 
                                         "#D53E4F",  "#67001F", "#053061"),
                                         labels = c("May 2021", "July", "Sept", "Nov", "Jan 2022")) +
  scale_x_discrete(labels = c("pico-POM", "nano-POM", "micro-POM")) + 
  labs(x=element_blank(), y = expression(paste("22:6",omega,"3")),
       fill = element_blank()) + 
  theme_bw()

ggplot(POM.sizes, aes(x=Pore.Size, y=C22.1n.9_PERCENT, fill = as.factor(Month))) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = c( "#A6D96A",  "#FDAE61", 
                                         "#D53E4F",  "#67001F", "#053061"),
                                         labels = c("May 2021", "July", "Sept", "Nov", "Jan 2022")) +
  scale_x_discrete(labels = c("pico-POM", "nano-POM", "micro-POM")) + 
  labs(x=element_blank(), y = expression(paste("22:1",omega,"9")),
       fill = element_blank()) + 
  theme_bw()


# Sign. IndVals of pico-POM
ggplot(POM.sizes, aes(x=Pore.Size, y=C18.3n.3_PERCENT, fill = as.factor(Month))) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = c( "#A6D96A",  "#FDAE61", 
                                         "#D53E4F",  "#67001F", "#053061"),
                                         labels = c("May 2021", "July", "Sept", "Nov", "Jan 2022")) +
  scale_x_discrete(labels = c("pico-POM", "nano-POM", "micro-POM")) + 
  labs(x=element_blank(), y = expression(paste("18:3",omega,"3")),
       fill = element_blank()) + 
  theme_bw()

ggplot(POM.sizes, aes(x=Pore.Size, y=C18.1n.7_PERCENT, fill = as.factor(Month))) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = c( "#A6D96A",  "#FDAE61", 
                                         "#D53E4F",  "#67001F", "#053061"),
                    labels = c("May 2021", "July", "Sept", "Nov", "Jan 2022")) +
  scale_x_discrete(labels = c("pico-POM", "nano-POM", "micro-POM")) + 
  labs(x=element_blank(), y = expression(paste("18:1",omega,"7")),
       fill = element_blank()) +  
    theme_bw()


# Nutritional quality


# EFA = 18C PUFA and longer; HUFA = 20C PUFA and longer
POM.sizes <- POM.sizes %>% mutate(EFA.mg.g = ((C18.2n.6c_PERCENT + C18.3n.6_PERCENT + C18.3n.3_PERCENT + C18.4n.3_PERCENT + 
                                                 C20.2n.6_PERCENT + C20.3n.3_PERCENT + C20.3n.6_PERCENT + C20.4n.6_PERCENT + 
                                                 C20.4n.3_PERCENT + C20.5n.3_PERCENT + C22.6n.3_PERCENT +C22.5n.6_PERCENT + 
                                                 C22.4n.6_PERCENT + C22.5n.3_PERCENT)* SumFA_mg.g))
POM.sizes <- POM.sizes %>% mutate(HUFA.mg.g = ((  C20.3n.3_PERCENT + C20.3n.6_PERCENT + C20.4n.6_PERCENT + 
                                                  C20.4n.3_PERCENT + C20.5n.3_PERCENT + C22.6n.3_PERCENT + C22.5n.6_PERCENT + 
                                                  C22.4n.6_PERCENT + C22.5n.3_PERCENT)* SumFA_mg.g))

POM.sizes <- POM.sizes %>% mutate(EFA.percent = ((C18.2n.6c_PERCENT + C18.3n.6_PERCENT + C18.3n.3_PERCENT + C18.4n.3_PERCENT + 
                                                    C20.2n.6_PERCENT + C20.3n.3_PERCENT + C20.3n.6_PERCENT + C20.4n.6_PERCENT + 
                                                    C20.4n.3_PERCENT + C20.5n.3_PERCENT + C22.6n.3_PERCENT + C22.5n.6_PERCENT + 
                                                    C22.4n.6_PERCENT + C22.5n.3_PERCENT)))
POM.sizes <- POM.sizes %>% mutate(HUFA.percent = ((C20.3n.3_PERCENT + C20.3n.6_PERCENT + C20.4n.6_PERCENT + C20.4n.3_PERCENT + 
                                                     C20.5n.3_PERCENT + C22.6n.3_PERCENT + C22.5n.6_PERCENT + C22.4n.6_PERCENT + 
                                                     C22.5n.3_PERCENT)))


# Fig. S3 -----------------------------------------------------------------

ggplot(POM.sizes, aes(x=as.factor(Month), y=SumFA_mg.g, fill = Pore.Size)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#b8e166", "#0073e6","#eb59a1"),
                    labels = c("pico-POM", "nano-POM", "micro-POM")) +
  scale_x_discrete(labels = c("May\n2021", "Jul\n2021", "Sep\n2021", 
                              "Nov\n2021", "Jan\n2022")) +
  theme_bw() + 
  labs(x=element_blank(), y=expression(paste("Total FA (mg g"^-1*")")), fill="Size Class")   


ggplot(POM.sizes, aes(x=as.factor(Month), y=EFA.mg.g, fill = Pore.Size)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#b8e166", "#0073e6","#eb59a1"),
                    labels = c("pico-POM", "nano-POM", "micro-POM")) +
  scale_x_discrete(labels = c("May\n2021", "Jul\n2021", "Sep\n2021", 
                              "Nov\n2021", "Jan\n2022")) +
  theme_bw() + 
  labs(x=element_blank(), y=expression(paste("EFA (mg g"^-1*")")), fill="Size Class")   


ggplot(POM.sizes, aes(x=as.factor(Month), y=HUFA.mg.g, fill = Pore.Size)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#b8e166", "#0073e6","#eb59a1"),
                    labels = c("pico-POM", "nano-POM", "micro-POM")) +
  scale_x_discrete(labels = c("May\n2021", "Jul\n2021", "Sep\n2021", 
                              "Nov\n2021", "Jan\n2022")) +
  theme_bw() + 
  labs(x=element_blank(), y=expression(paste("HUFA (mg g"^-1*")")), fill="Size Class")   




mod <- lm(POM.sizes$SumFA_ug.L ~ POM.sizes$Pore.Size + as.factor(POM.sizes$Month))
summary(mod)
# July different from Sept, Nov, and Jan but not May
# Jan different from Sept, May, and July but not Nov
# May different from Nov, Sept, Jan but not July
# Sept different from all

plot(mod)


leveneTest(lm(SumFA_ug.L ~ Pore.Size*as.factor(POM.sizes$Month), data=POM.sizes))
leveneTest(lm(EFA.mg.g ~ Pore.Size*as.factor(POM.sizes$Month), data=POM.sizes))
leveneTest(lm(HUFA.mg.g ~ Pore.Size*as.factor(POM.sizes$Month), data=POM.sizes))
# no evidence of unequal variances


SumFA_ug.L.corr.AOV <- aov(POM.sizes$SumFA_ug.L ~ POM.sizes$Pore.Size * as.factor(POM.sizes$Month))
summary(SumFA_ug.L.corr.AOV)
# All significant
TukeyHSD(SumFA_ug.L.corr.AOV)
#  All three size fractions distinct
# 20-0.7  0.015
# 3-0.7   0.00002
# 3-20    0.033


EFA.mg.g.corr.AOV <- aov(POM.sizes$EFA.mg.g ~ POM.sizes$Pore.Size * as.factor(POM.sizes$Month))
summary(EFA.mg.g.corr.AOV)
# All significant
TukeyHSD(EFA.mg.g.corr.AOV)
# 20-0.7  0.0240095
# 3-0.7   0.0000006
# 3-20    0.0006521


HUFA.mg.g.corr.AOV <- aov(POM.sizes$HUFA.mg.g ~ POM.sizes$Pore.Size * as.factor(POM.sizes$Month))
summary(HUFA.mg.g.corr.AOV)
# Micro and pico not significantly different
TukeyHSD(HUFA.mg.g.corr.AOV)
# 20-0.7  0.1915106
# 3-0.7   0.0001905
# 3-20    0.0000023




# CTD data ----------------------------------------------------------------



# # CTD data summary 
CTD_all <- read.csv("processed_data/CTD summary stats for each day_20251108.csv")
CTD_all$Date <- as.Date(CTD_all$Date, "%Y-%m-%d")


# *Fig. S9 ----------------------------------------------------------------


ggplot(CTD_all, aes(x=Date, y=SST)) + geom_point() + theme_bw() + 
  geom_vline(xintercept = c(as.Date(c("2021-05-20", "2021-07-22", "2021-09-02", 
                                      "2021-11-04", "2022-01-15","2022-03-30"))),
             color = "grey60") + 
  scale_x_date(limits = as.Date(c("2021-01-01", "2022-05-01")), expand = c(0,0),
               breaks = as.Date(c("2021-01-01", "2021-03-01",
                                  "2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("Jan\n2021","Mar\n2021", "May\n2021", "Jul\n2021", "Sep\n2021", 
                          "Nov\n2021", "Jan\n2022", "Mar\n2022"))  + 
  theme(axis.title.x = element_blank()) + 
  labs(y="temperature (°C)")


ggplot(CTD_all, aes(x=Date, y=SSS)) + geom_point() + theme_bw() + 
  geom_vline(xintercept = c(as.Date(c("2021-05-20", "2021-07-22", "2021-09-02", 
                                      "2021-11-04", "2022-01-15","2022-03-30"))),
             color = "grey60") + 
  scale_x_date(limits = as.Date(c("2021-01-01", "2022-05-01")), expand = c(0,0),
               breaks = as.Date(c("2021-01-01", "2021-03-01",
                                  "2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("Jan\n2021","Mar\n2021", "May\n2021", "Jul\n2021", "Sep\n2021", 
                          "Nov\n2021", "Jan\n2022", "Mar\n2022"))  + 
  theme(axis.title.x = element_blank()) + 
  labs(y="salinity (psu)")


ggplot(CTD_all, aes(x=Date, y=Mean.Chl)) + geom_point() + theme_bw() + 
  geom_vline(xintercept = c(as.Date(c("2021-05-20", "2021-07-22", "2021-09-02", 
                                      "2021-11-04", "2022-01-15","2022-03-30"))),
             color = "grey60") + 
  scale_x_date(limits = as.Date(c("2021-01-01", "2022-05-01")), expand = c(0,0),
               breaks = as.Date(c("2021-01-01", "2021-03-01",
                                  "2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("Jan\n2021","Mar\n2021", "May\n2021", "Jul\n2021", "Sep\n2021", 
                          "Nov\n2021", "Jan\n2022", "Mar\n2022"))  + 
  theme(axis.title.x = element_blank()) + 
  labs(y="Chl fluorescence")




# POM FA NMDS -----------------------------------------------------------------

# Combine POM FA and SI
SI.POM$Month <- as.factor(SI.POM$Month)
POM.sizes$Month

POM.sizes.wSI_2 <- SI.POM %>% select(-Month) %>% right_join( POM.sizes)

# Need to average by survey prior to joining 
# # CTD data summary 
# CTD.data <- read.csv("/Users/amclaskey/My Drive/UBC Work/Work/Research/2021 Zoopsprint/Data/CTD summary stats for each day_20251108.csv")

# For Zoopsprint dates, calculate 
ZSdates <- c("2021-05-20", "2021-07-22", "2021-09-02", "2021-11-04", "2022-01-15", "2022-03-30")
ZSdates <- as.Date(ZSdates, "%Y-%m-%d")

# Calculate 14 day moving average of POM DelN15
CTD.dataframe <- data.frame(Date = ZSdates)
for(i in 1:length(ZSdates)){
  CTD.dataframe$SST[i] <- mean(CTD_all$SST[CTD_all$Date >= (ZSdates[i]-14) & 
                                              CTD_all$Date<=ZSdates[i] ], na.rm = T)
  CTD.dataframe$SSS[i] <- mean(CTD_all$SSS[CTD_all$Date >= (ZSdates[i]-14) & 
                                              CTD_all$Date<=ZSdates[i] ])  
  CTD.dataframe$Chl[i] <- mean(CTD_all$Mean.Chl[CTD_all$Date >= (ZSdates[i]-14) & 
                                                   CTD_all$Date<=ZSdates[i] ]) } 
CTD.dataframe$Month <- format(as.Date(CTD.dataframe$Date), "%m")
CTD.dataframe$Month <- as.integer(CTD.dataframe$Month)

CTD.data.monthlymeans <- CTD.dataframe %>% select(Month, SST, SSS, Chl)


CTD.data.monthlymeans$Month <- (as.character(CTD.data.monthlymeans$Month))
CTD.data.monthlymeans$Month[CTD.data.monthlymeans$Month=="1"] <- "13"

POM.sizes.wSI <- full_join(POM.sizes.wSI_2, CTD.data.monthlymeans) 

# remove CTD data from March - no size frac POM samples
POM.sizes.wSI <- POM.sizes.wSI[-40,]



covariates <- POM.sizes.wSI %>% select(delta13c, delta15n, SST, SSS, Chl)
colnames(covariates)[c(1:2)] <- c("del13C", "del15N")
taxa.names <- as.factor(POM.sizes.wSI$Pore.Size)



# Make reduced data matrix
POM.matrix <- POM.sizes.wSI %>%  dplyr::select(C14.0_PERCENT:iso.17.0_PERCENT, -C20.0_PERCENT, -C24.1n.9_PERCENT)

# What are the average contributions of each FA
contributions <- data.frame(matrix(nrow=ncol(POM.matrix), ncol=3))
for(i in 1:ncol(POM.matrix)){
  contributions[i,1] <- names(POM.matrix[i])
  contributions[i,2] <- mean(POM.matrix[,i], na.rm = T)*100
  temp <- (POM.matrix[i]==0)
  contributions[i,3] <- length(temp[temp==TRUE])
}
names(contributions) <- c("FA", "Mean.Percentage", "Num zeros")

# separate out peaks that are >1%
abundant.pom <- contributions$FA[contributions$Mean.Percentage>=0.5]
abundant.pom

POM.matrix.abund <- POM.matrix[,abundant.pom]
colnames(POM.matrix.abund) <- paste(str_remove(colnames(POM.matrix.abund), "_PERCENT"), "", sep = "")

POM.sizes.wSI$Month <- format(as.Date(POM.sizes.wSI$Date), "%m")
POM.sizes.wSI$Month <- as.factor(POM.sizes.wSI$Month)
Months <- as.factor(POM.sizes.wSI$Month)
sizes <- as.factor(POM.sizes.wSI$Pore.Size)


#create matrix for NMDS calculation:
species.matrix.sm <- as.matrix(POM.matrix.abund)
#change diet dataframe into a matrix
class(species.matrix.sm) <- "numeric"
#make sure your numbers are treated as numbers



# *Table S3 ---------------------------------------------------------------
# PERMANOVA of POM FAs by size class and collection month

permanova_eco.bc<-adonis2(species.matrix.sm ~ Months*sizes, permutations = 999, method="bray")
permanova_eco.bc #if significant, then plot it

set.seed(2)
eco.nmds.bc <- metaMDS(species.matrix.sm, distance="bray",labels=Months, trymax = 100, autotransform = FALSE)
# eco.nmds.bc <- metaMDS(transformed_matrix, distance="bray",labels=taxa.names, trymax = 100, autotransform = FALSE, k=3)
eco.nmds.bc
stressplot(eco.nmds.bc)


species.scores <- as.data.frame(scores(eco.nmds.bc, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- paste(str_remove(row.names(species.scores), "_PERCENT"), "", sep = "")

head(species.scores)  #look at the data



# This code is from Vanessa
NMDS.bc<-data.frame(NMDS1.bc=eco.nmds.bc$points[,1],NMDS2.bc=eco.nmds.bc$points[,2], size=sizes)
# dataframe for plotting NMDS

# Working on putting species (FAs) onto plot
species.dataframe <- data.frame(FA = rownames(eco.nmds.bc$species), 
                                NMDS.1 = eco.nmds.bc$species[,1],
                                NMDS.2 = eco.nmds.bc$species[,2])


NMDS.bc<- NMDS.bc %>%  mutate(size = factor(size, 
                                            levels = c("0.7", "3","20")))




# Project SI POM over FA ordination

(fit <- envfit(eco.nmds.bc, covariates, perm = 999, na.rm = T))
# Del15N not significantly correlated p = 0.099
# Del13C  significant p = 0.001
# SST is significant p = 0.03
summary(fit)
plot(eco.nmds.bc)
plot(fit)
#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(eco.nmds.bc)$sites)
#add 'season' column as before
data.scores$taxa.names = sizes

en_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit)


# This code is from Vanessa
NMDS.bc<-data.frame(NMDS1.bc=eco.nmds.bc$points[,1],NMDS2.bc=eco.nmds.bc$points[,2],
                    group=sizes)
NMDS.bc<- NMDS.bc %>%  mutate(group = factor(group, 
                                            levels = c("0.7", "3","20")))


####

# https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Ellipses are standard deviation, no scaling of data (can use standard error, scaling, and confidence limit options)

plot.new()
ord.bc<-ordiellipse(eco.nmds.bc,sizes,display="sites",kind="sd", conf = 0.95, label=T)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


df_ell.bc <- data.frame()
for(g in levels(NMDS.bc$group)){
  df_ell.bc <- rbind(df_ell.bc, cbind(as.data.frame(with(NMDS.bc[NMDS.bc$group==g,],
                                                         veganCovEllipse(ord.bc[[g]]$cov,ord.bc[[g]]$center))),group=g))
}
#https://www.rpubs.com/RGrieger/545184


df_ell.bc<- df_ell.bc %>%  mutate(group = factor(group, 
                                             levels = c("0.7", "3","20")))


# Working on putting species (FAs) onto plot
species.dataframe <- data.frame(FA = rownames(eco.nmds.bc$species), 
                                NMDS.1 = eco.nmds.bc$species[,1],
                                NMDS.2 = eco.nmds.bc$species[,2])

####

#  *Fig. 2C -----------------------------------------------------------------

POM.NMDS <- ggplot(NMDS.bc, aes(NMDS1.bc, NMDS2.bc))+ theme_bw() + 
  geom_segment(data = species.scores, aes(x = 0, xend=NMDS1*1, y=0, yend=NMDS2*1), arrow = arrow(length = unit(0.15, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  # geom_text(data = species.dataframe, aes(NMDS.1, NMDS.2, label=FA)) + 
  geom_path(data=df_ell.bc, aes(x=NMDS1, y=NMDS2,colour=group), size=1, linetype=1) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont[c(1,3),], lwd=0.6, colour = "black") +
  geom_point(stat = "identity", aes(color= group, shape=Months, fill=group), size=1.5, stroke = 0.8) +
  scale_shape_manual(values = c(16, 17,15,3,7), 
    labels = c("Jan 2022", "May 2021", "July", "Sept", "Nov"))  +
  scale_color_manual(values = c("#89ce00", "#0073e6","#e6308a"),
                    labels=c(expression(paste("0.7 ",mu,"m")), 
                             expression(paste("3 ",mu,"m")),
                             expression(paste("20 ",mu,"m")))) + 
  #  14:0
  annotate("text", x=species.dataframe[c(1),2]*1.15, y=(species.dataframe[c(1),3]*1.2), 
           label=c("14:0"),  size=3) +
  #  15:0
  annotate("text", x=species.dataframe[c(2),2]*1.2, y=(species.dataframe[c(2),3]*1.25), 
           label=c("15:0"),  size=3) +
  #  16:0
  annotate("text", x=species.dataframe[c(3),2]*1.8, y=(species.dataframe[c(3),3]*1.2), 
           label=c("16:0"),  size=3) +
  # 16:1n-7
  annotate("text", x=species.dataframe[c(4),2]*1.2, y=species.dataframe[c(4),3]*1, 
           label=expression(paste("16:1",omega,"7")), fontface =2, size=3) + 
  # 16:2n-4
  annotate("text", x=species.dataframe[c(5),2]*1.1, y=(species.dataframe[c(5),3]*1.05), 
           label=expression(paste("16:2",omega,"4")), fontface =2, size=3) +
  # 16:4n-1
  annotate("text", x=species.dataframe[c(6),2]*1.25, y=(species.dataframe[c(6),3])*1, 
           label=expression(paste("16:4",omega,"1")), fontface =2, size=3) + 
  # 17:0
  annotate("text", x=species.dataframe[c(7),2]*1.1, y=(species.dataframe[c(7),3]*1.2), 
           label="17:0",  size=3)  + 
  # 18:0
  annotate("text", x=species.dataframe[c(8),2]*1.1, y=(species.dataframe[c(8),3]*1.2), 
           label="18:0", size=3) +
  # 18:1n-7
  annotate("text", x=species.dataframe[c(9),2]*1.15, y=(species.dataframe[c(9),3]*1.11), 
           label=expression(paste("18:1",omega,"7")), size=3)  + 
  # 18:1n-9
  annotate("text", x=species.dataframe[c(10),2]*1.1, y=(species.dataframe[c(10),3]*1.1), 
           label=expression(paste("18:1",omega,"9")), size=3)  + 
  # 18:2n-6
  annotate("text", x=species.dataframe[c(11),2]*1.1, y=(species.dataframe[c(11),3]*1.1), 
           label=expression(paste("18:2",omega,"6")), size=3)  + 
  # ALA
  annotate("text", x=species.dataframe[c(12),2]*1.1, y=(species.dataframe[c(12),3]*1.1), 
           label=expression(paste("18:3",omega,"3")), size=3) + 
  # SDA
  annotate("text", x=species.dataframe[c(13),2]*1.34, y=(species.dataframe[c(13),3]*1.25),
           label=expression(paste("18:4",omega,"3")),size=3) +
  # EPA
  annotate("text", x=species.dataframe[c(14),2]*1.17, y=(species.dataframe[c(14),3]*2.2), 
           label=expression(paste("20:5",omega,"3")), size=3) +
  # 22:1n-9
  annotate("text", x=species.dataframe[c(15),2]*1.1, y=(species.dataframe[c(15),3]*1.1), 
           label=expression(paste("22:1",omega,"9")), size=3) +
  # DHA
  annotate("text", x=species.dataframe[c(16),2]*1.3, y=(species.dataframe[c(16),3]*1.1), 
           label=expression(paste("22:6",omega,"3")), size=3) +
  labs(x="NMDS1",
       y="NMDS2",
       shape = "Survey",
       color = "Size Class") + 
  annotate("text", x=en_coord_cont[c(1),1]*1.08, y=en_coord_cont[c(1),2]*1.05,
           label=expression(paste(delta^{13}, "C")),  size=3) +
  # annotate("text", x=en_coord_cont[c(2),1]*1.1, y=en_coord_cont[c(2),2]*1.1, 
  #          label=expression(paste(delta^{15}, "N")),  size=3) +
  annotate("text", x=en_coord_cont[c(3),1]*1.15, y=en_coord_cont[c(3),2], 
           label="SST",  size=3) +
  theme(panel.grid = element_blank(),
        legend.key.size = unit(0.8, 'mm'), #change legend key size
        legend.spacing.y = unit(2.0, "mm"),
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8),
        axis.title = element_text(size=9),
        axis.text = element_text(size=7),
        plot.margin = margin(1, 0, 0, 0), 
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(3, "pt"),
        legend.title.align = 0,
        axis.title.y = element_text(margin = margin(0, 0, 0, 0),
                                    vjust = 0))

POM.NMDS 

# tiff("POMFA NMDS w SI 20251125.tiff", width=110, height=100, units="mm", res=500)
# POM.NMDS
# dev.off()








# Zoop FA NMDS ------------------------------------------------------------


SI.zoop <- SI.all %>% filter(Species!="POM") %>% filter(Species!="z.250") %>% filter(Species!="March paraeuchaeta eggs") %>% 
   select(Month, Species, 
          C_N, delta13c, delta15n)
SI.zoop.monthlymeans <- SI.zoop %>% group_by(Species, Month) %>% 
  summarise_all( mean)


# Combine with survey CTD data
SI.zoop.monthlymeans$Month <- as.factor(SI.zoop.monthlymeans$Month)
zoopSI.CTD <- full_join(SI.zoop.monthlymeans, CTD.data.monthlymeans)

zoopData.sm$Month <- as.character(zoopData.sm$Month)
zoopData.sm$Month[zoopData.sm$Month=="01"] <- "13"
zoopData.sm$Month[zoopData.sm$Month=="03"] <- "15"
# Join Zoop species FA and SI data
zoopData.wSI.CTD <- left_join(zoopData.sm, zoopSI.CTD)




covariates <- zoopData.wSI.CTD %>% select(delta13c, delta15n, SST, SSS, Chl)
# colnames(covariates)[c(1:2)] <- c("del13C", "del15N")
taxa.names <- as.factor(zoopData.wSI.CTD$Species)
# zoopData.wSI.CTD$Month <- format(as.Date(zoopData.wSI.CTD$Date), "%m")
# zoopData.wSI.CTD$Month <- as.factor(zoopData.wSI.CTD$Month)
Months <- as.factor(zoopData.wSI.CTD$Month)
Seasons <- as.factor(zoopData.wSI.CTD$season)

# Make reduced data matrix
Q20.matrix <- zoopData.wSI.CTD %>%  dplyr::select(C14.0_PERCENT:iso.17.0_PERCENT)


# What are the average contributions of each FA
contributions <- data.frame(matrix(nrow=ncol(Q20.matrix), ncol=3))
for(i in 1:ncol(Q20.matrix)){
  contributions[i,1] <- names(Q20.matrix[i])
  contributions[i,2] <- mean(Q20.matrix[,i], na.rm = T)*100
  temp <- (Q20.matrix[i]==0)
  contributions[i,3] <- length(temp[temp==TRUE])
}
names(contributions) <- c("FA", "Mean.Percentage", "Num zeros")
# write.csv(contributions, "Q20 zoop contributions.csv", row.names = F)

# separate out peaks that are >1%
abundant.z <- contributions$FA[contributions$Mean.Percentage>=0.5]
# take out C15:0
abundant.z
abundant.z <- abundant.z[-2]


Q20.matrix.abund <- Q20.matrix[,abundant.z]

colnames(Q20.matrix.abund) <- paste(str_remove(colnames(Q20.matrix.abund), "_PERCENT"), "", sep = "")





# *Table S7 ---------------------------------------------------------------

# make a summary table of FAs by taxa

fullsummary <- zoopData.sm %>% select(Species, all_of(abundant.z), Bacteria_all, percent.SFA, percent.MUFA, percent.PUFA, DHA.EPA, SumFA_mg.g, SumFA_mg.gWW) %>% 
  group_by(Species) %>% summarise_all(mean)

# fullsummary <- zsAllw.POM %>% select(Species, C12.0_PERCENT:iso.17.0_PERCENT, 
#                                      Bacteria_all, percent.SFA, percent.MUFA, percent.PUFA, 
#                                      DHA.EPA, C16PUFA, C18.3andC18.4, Copes.C20C.1C22.1,
#                                      SumFA_mg.g, SumFA_mg.gWW, SumFA_ug.L) %>% 
#   group_by(Species) %>% summarise_all(mean)


 # write.csv(fullsummary, "Taxa FA mean 2025-07-17.csv")




# taxa.names <- as.factor(zoopData.sm$Species)
# zoopData.sm$Month <- format(as.Date(zoopData.sm$Date), "%m")
# zoopData.sm$Month <- as.factor(zoopData.sm$Month)
# Months <- as.factor(zoopData.sm$Month)
# Seasons <- as.factor(zoopData.sm$season)


Q20.matrix.abund <- Q20.matrix[,abundant.z]
colnames(Q20.matrix.abund) <- paste(str_remove(colnames(Q20.matrix.abund), "_PERCENT"), "", sep = "")

Q20.matrix.abund[is.na(Q20.matrix.abund)] <- 0

#create matrix for NMDS calculation:
species.matrix.sm <- as.matrix(Q20.matrix.abund)
#change diet dataframe into a matrix
class(species.matrix.sm) <- "numeric"
#make sure your numbers are treated as numbers
proportions_matrix <- decostand(species.matrix.sm, "total")
# transformed_matrix <- asin(sqrt(species.matrix.sm))


# start here if all loaded from above 

# 
bray_curtis_matrix <- vegdist(species.matrix.sm, "bray")

perm.eg.betadisper <- betadisper(bray_curtis_matrix, 
                                 taxa.names)

anova(perm.eg.betadisper)
# Months is non-significant (p=0.32) 
# species are significant (p<0.0001)

zoopData.sm$Species.Month <- paste(zoopData.sm$Species, zoopData.sm$Month, sep=".")
Species.Months <- as.factor(zoopData.sm$Species.Month)

perm.eg.betadisper <- betadisper(bray_curtis_matrix, 
                                 Species.Months)
anova(perm.eg.betadisper)

# when interaction is considered, species*months have homogeneous dispersions



# *Table S5 ---------------------------------------------------------------

permanova_eco.bc<-adonis2(species.matrix.sm ~ taxa.names*Months, permutations = 999, method="bray")
permanova_eco.bc #if significant, then plot it


# set.seed(91203)
set.seed(91203)

eco.nmds.bc <- metaMDS(species.matrix.sm, distance="bray",labels=taxa.names, trymax = 100, autotransform = FALSE)
# eco.nmds.bc <- metaMDS(transformed_matrix, distance="bray",labels=taxa.names, trymax = 100, autotransform = FALSE, k=3)
eco.nmds.bc
stressplot(eco.nmds.bc)

eco.nmds.bc$species

species.scores <- as.data.frame(scores(eco.nmds.bc, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- paste(str_remove(row.names(species.scores), "_PERCENT"), "", sep = "")

head(species.scores)  #look at the data
sort(species.scores$NMDS1)
sort(species.scores$NMDS2)

# write.csv(species.scores, "species.scores.csv")


# This code is from Vanessa
NMDS.bc<-data.frame(NMDS1.bc=eco.nmds.bc$points[,1],NMDS2.bc=eco.nmds.bc$points[,2], month=Months, 
                    sizes=taxa.names)
# dataframe for plotting NMDS


ggplot(NMDS.bc, aes(NMDS1.bc, NMDS2.bc))+ theme_bw() + 
  geom_point(stat = "identity", aes(color=taxa.names), size=2) + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.01) + 
  scale_color_manual(values = clrs) 




#dataframe for plotting NMDS
NMDS.bc<-data.frame(NMDS1.bc=eco.nmds.bc$points[,1],NMDS2.bc=eco.nmds.bc$points[,2],group=taxa.names, month=Months)

zoopData.sm$NMDS1 <- eco.nmds.bc$points[,1]
zoopData.sm$NMDS2 <- eco.nmds.bc$points[,2]

# https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Ellipses are standard deviation, no scaling of data (can use standard error, scaling, and confidence limit options)

plot.new()
ord.bc<-ordiellipse(eco.nmds.bc,taxa.names,display="sites",kind="sd", conf = 0.95, label=T)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


df_ell.bc <- data.frame()
for(g in levels(NMDS.bc$group)){
  df_ell.bc <- rbind(df_ell.bc, cbind(as.data.frame(with(NMDS.bc[NMDS.bc$group==g,],
                                                         veganCovEllipse(ord.bc[[g]]$cov,ord.bc[[g]]$center))),group=g))
}
#https://www.rpubs.com/RGrieger/545184




# Working on putting species (FAs) onto plot
species.dataframe <- data.frame(FA = rownames(eco.nmds.bc$species), 
                                NMDS.1 = eco.nmds.bc$species[,1],
                                NMDS.2 = eco.nmds.bc$species[,2])


(fit <- envfit(eco.nmds.bc, covariates, perm = 999, na.rm = T))
# Del15N significantly correlated p = 0.001
# Del13C nsignificantly correlated p = 0.002
# No enviro variables are significant 
summary(fit)
plot(eco.nmds.bc)
plot(fit)
#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(eco.nmds.bc)$sites)
#add 'season' column as before
data.scores$taxa.names = taxa.names

en_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit)






# *Fig. S7 ----------------------------------------------------------------


NMDS1<- ggplot(NMDS.bc, aes(NMDS1.bc, NMDS2.bc))+
  #to change color: fill = other_ids
  geom_path(data=df_ell.bc, aes(x=NMDS1, y=NMDS2,colour=group), size=1, linetype=1) +
  scale_fill_manual(values=clrs, name="Size", guide="legend") +
  guides(fill= guide_legend(override.aes = list(shape=21)))+
  scale_color_manual(values=c(clrs), guide=FALSE) + 
  theme_bw()+
  scale_x_continuous(limits = c(-0.6,1.1)) +
  theme(axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(angle=90,size=9),
        axis.text.y=element_text(size=9),
        panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        legend.position = "none")+
  geom_point(stat = "identity", aes(fill=taxa.names, color=taxa.names), size=2) +
  geom_segment(data = species.scores, aes(x = 0, xend=NMDS1*0.9, y=0, yend=NMDS2*0.9), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  # segments for bioenv analysis
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              data = en_coord_cont, lwd=0.6, colour = "red", alpha = 0.4) +
  # annotate("text", x=en_coord_cont[,1]*1, y=en_coord_cont[,2]*1, 
  #          label=row.names(en_coord_cont),  size=4, colour = "red") +
  # scale_shape_manual(values = c(15:18,20,23,25,8,11,12)) +
  # geom_text(data = species.dataframe, aes(NMDS.1, NMDS.2, label=FA)) + 
  # 14:0 
  annotate("text", x=species.dataframe[c(1),2]*1, y=species.dataframe[c(1),3]*0.95, 
           label="14:0",  size=4) + 
  #  16:0
  annotate("text", x=species.dataframe[c(2),2]*1.6, y=(species.dataframe[c(2),3]*1.3), 
           label=c("16:0"),  size=4) +
  # 16:1n-7
  annotate("text", x=species.dataframe[c(3),2]*0.5, y=species.dataframe[c(3),3]*1.05, 
           label=expression(paste("16:1",omega,"7")), fontface =2, size=4) + 
  # 16:2n-4
  annotate("text", x=species.dataframe[c(4),2]*0.8, y=(species.dataframe[c(4),3]*1.05), 
           label=expression(paste("16:2",omega,"4")), fontface =2, size=4) +
  # 16:3n-4
  annotate("text", x=species.dataframe[c(5),2]*1.1, y=(species.dataframe[c(5),3]*0.95), 
           label=expression(paste("16:3",omega,"4")), fontface =2, size=4) +
  # 16:4n-1
  annotate("text", x=species.dataframe[c(6),2]*0.95, y=(species.dataframe[c(6),3])*1.0, 
           label=expression(paste("16:4",omega,"1")), fontface =2, size=4) + 
  # 17:0
  annotate("text", x=species.dataframe[c(7),2], y=(species.dataframe[c(7),3]*1), 
           label="17:0",  size=4)  + 
  # 18:0
  annotate("text", x=species.dataframe[c(8),2]*2, y=(species.dataframe[c(8),3]*1.05), 
           label="18:0", size=4) +
  # 18:1n-7
  annotate("text", x=species.dataframe[c(9),2]*1.6, y=(species.dataframe[c(9),3]*0.95), 
           label=expression(paste("18:1",omega,"7")), size=4)  + 
  # 18:1n-9
  annotate("text", x=species.dataframe[c(10),2]*0.8, y=(species.dataframe[c(10),3]*1.05), 
           label=expression(paste("18:1",omega,"9")), size=4)  + 
  # ALA
  annotate("text", x=species.dataframe[c(12),2]*1.4, y=(species.dataframe[c(12),3]*0.9), 
           label=expression(paste("18:3",omega,"3")), size=4) + 
  # 20:1n-11
  annotate("text", x=species.dataframe[c(14),2]*0.9, y=(species.dataframe[c(14),3]*1.25),
           label=expression(paste("20:1",omega,"11")),size=4) +
  # 20:1n-9
  annotate("text", x=species.dataframe[c(15),2]*0.85, y=(species.dataframe[c(15),3]), 
           label=expression(paste("20:1",omega,"9")), size=4) +
  # EPA
  annotate("text", x=species.dataframe[c(17),2]*1.5, y=(species.dataframe[c(17),3]*1.2), 
           label=expression(paste("20:5",omega,"3")), size=4) +
  # 22:1n-11
  annotate("text", x=species.dataframe[c(18),2]*0.82, y=(species.dataframe[c(18),3]*0.95), 
           label=expression(paste("22:1",omega,"11")), size=4) +
  # DHA
  annotate("text", x=species.dataframe[c(20),2]*1.1, y=(species.dataframe[c(20),3]*1.05), 
           label=expression(paste("22:6",omega,"3")), size=4) +
  labs(x="NMDS1",
       y="NMDS2")

NMDS1

# png("ZS NMDS 20250610.png", width=120, height=120, units="mm", res=300)
# NMDS1
# dev.off()






inset <- ggplot(NMDS.bc, aes(NMDS1.bc, NMDS2.bc))+ theme_bw() + 
  geom_segment(data = species.scores, aes(x = 0, xend=NMDS1*0.9, y=0, yend=NMDS2*0.9), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  #  16:0
  annotate("text", x=species.dataframe[c(2),2]*1.3, y=(species.dataframe[c(2),3]*0.8), 
           label=c("16:0"),  size=4) +
  # 18:0
  annotate("text", x=species.dataframe[c(8),2]*2.8, y=(species.dataframe[c(8),3]*0.95), 
           label="18:0", size=4) +
  # 18:1n-7
  annotate("text", x=species.dataframe[c(9),2]*1.3, y=(species.dataframe[c(9),3]*0.8), 
           label=expression(paste("18:1",omega,"7")), size=4)  + 
  # 18:2n-6
  annotate("text", x=species.dataframe[c(11),2]*1.6, y=(species.dataframe[c(11),3]*1.05), 
           label=expression(paste("18:2",omega,"6")), size=4)  + 
  # ALA
  annotate("text", x=species.dataframe[c(12),2]*1.4, y=(species.dataframe[c(12),3]*0.8), 
           label=expression(paste("18:3",omega,"3")), size=4) + 
  # SDA 
  annotate("text", x=species.dataframe[c(13),2]*-1.2, y=(species.dataframe[c(13),3]*1.2),
           label=expression(paste("18:4",omega,"3")),size=4) +
  # ARA
  annotate("text", x=species.dataframe[c(16),2]*5, y=(species.dataframe[c(16),3])*0.85, 
           label=expression(paste("20:4",omega,"6")), size=4) +
  # EPA
  annotate("text", x=species.dataframe[c(17),2]*1.3, y=(species.dataframe[c(17),3]*1.05), 
           label=expression(paste("20:5",omega,"3")), size=4) +
  # DPA
  annotate("text", x=species.dataframe[c(19),2]*0.8, y=(species.dataframe[c(19),3]*0.5), 
           label=expression(paste("22:5",omega,"3")), size=4) +
  # DHA
  annotate("text", x=species.dataframe[c(20),2]*1.6, y=(species.dataframe[c(20),3]*0.85), 
           label=expression(paste("22:6",omega,"3")), size=4) +
  # 24:1n-9
  annotate("text", x=species.dataframe[c(21),2]*1.9, y=(species.dataframe[c(21),3]*1.05), 
           label=expression(paste("24:1",omega,"9")), size=4) +
  labs(x="NMDS1",
       y="NMDS2") + 
  theme(panel.grid = element_blank())
inset

# png("ZS NMDS inset 20240709.png", width=170, height=170, units="mm", res=300)
# inset
# dev.off()




# *Fig S5 A,B zoop SI time series ----------------------------------------------------

Zoop.SI <- SI.all %>% filter(Species!= "POM" & Species!= "z.250" & Species!= "March paraeuchaeta eggs")

zoopsprint.SI.means <- Zoop.SI %>% dplyr::select(Species, Month,
                                                 delta13c, delta15n, delta13c.notlipidcorrected, 
                                                 Date) %>% 
  group_by(Species, Month) %>% 
  summarise_all( mean)


colnames(zoopsprint.SI.means)[3:5] <- paste0("mean.", colnames(zoopsprint.SI.means)[3:5])

Zoop.SI.wMeans <- zoopsprint.SI.means %>% dplyr::select(-Date) %>% 
  full_join(.,Zoop.SI, by = join_by(Species,Month))


pDel13C <- ggplot(Zoop.SI.wMeans, aes(x=Date, y=mean.delta13c, color=Species, fill=Species)) + 
  geom_jitter(aes(x=Date, y=delta13c), width = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none")  +  
  geom_line(aes(linetype = Species), size=1.2) +
  # scale_y_continuous(limits = c(0, 2.5), expand = c(0,0)) + 
  scale_linetype_manual(values = c( "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values =  clrs)  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May\n2021", "Jul\n2021", "Sep\n2021", 
                          "Nov\n2021", "Jan\n2022", "Mar\n2022")) +
  labs(x=element_blank(), y=expression(paste(delta^{13}, "C (\u2030)")),
       color="Species" , fill="Species") 
pDel13C

# tiff("Zoop d13C 20250619.tiff", width=110, height=100, units="mm", res=500)
# pDel13C
# dev.off() 

pDel15N <- ggplot(Zoop.SI.wMeans, aes(x=Date, y=mean.delta15n, color=Species, fill=Species)) + 
  geom_jitter(aes(x=Date, y=delta15n), width = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none")  +  
  geom_line(aes(linetype = Species), size=1.2) +
  # scale_y_continuous(limits = c(0, 2.5), expand = c(0,0)) + 
  scale_linetype_manual(values = c( "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid")) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values =  clrs)  +
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May\n2021", "Jul\n2021", "Sep\n2021", 
                          "Nov\n2021", "Jan\n2022", "Mar\n2022"))  + 
  labs(x=element_blank(), y=expression(paste(delta^{15}, "N (\u2030)")),
       color="Species" , fill="Species") 
pDel15N

# tiff("Zoop d15N 20250619.tiff", width=110, height=100, units="mm", res=500)
# pDel15N
# dev.off() 

prow <- plot_grid(pDel13C, 
                  pDel15N,
                  TFA.plot,
                  rel_widths = c(1.1, 1, 1),
                  # rel_heights = c(0.9, 1),
                  nrow = 3)
prow

# tiff("Fig S4 time series 20250804.tiff", width=100, height=250, units="mm", res=500)
# prow
# dev.off()







         


# *Table S6 ---------------------------------------------------------------

allSIdata.sm <- Zoop.SI %>% select(Species, Month, delta13c.notlipidcorrected, delta15n)
colnames(allSIdata.sm)[4] <- c("delta15n")

taxa.names <- as.factor(allSIdata.sm$Species)
Months <- as.factor(allSIdata.sm$Month)

SI.coretaxa.sm <- allSIdata.sm[,3:4]
SI.coretaxa.sm$delta13c <- SI.coretaxa.sm$delta13c *-1


permanova_eco.bc<-adonis2(SI.coretaxa.sm ~ taxa.names*Months, permutations = 999, method="bray")
permanova_eco.bc





# *Fig. S4 ---------------------------------------------------------------


Zoop.sizes.SI <- SI.all %>% filter( Species!= "z.250" & Species!= "March paraeuchaeta eggs")

Zoop.sizes.SI$delta13c.notlipidcorrected[Zoop.sizes.SI$Species=="POM"] <- Zoop.sizes.SI$delta13c[Zoop.sizes.SI$Species=="POM"]
Zoop.sizes.SI$Species[Zoop.sizes.SI$Species=="POM"] <- Zoop.sizes.SI$Pore.Size[Zoop.sizes.SI$Species=="POM"]
Zoop.sizes.SI <- Zoop.sizes.SI %>% filter(Species!="bulk")

Zoop.sizes.SI.sm <- Zoop.sizes.SI %>% select(Species, Month, delta13c.notlipidcorrected, delta15n)

taxa.names <- as.factor(Zoop.sizes.SI.sm$Species)
Months <- as.factor(Zoop.sizes.SI.sm$Month)

#dataframe for plotting NMDS
NMDS.bc.SI<-data.frame(NMDS1.bc=Zoop.sizes.SI.sm[,3],NMDS2.bc=Zoop.sizes.SI.sm[,4], group=taxa.names, Month=Zoop.sizes.SI.sm[,1])

plot.new()
ord.bc.SI <-ordiellipse(Zoop.sizes.SI.sm[,3:4], taxa.names, display="sites",kind="sd", conf = 0.95, label=T)


veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell.bc.SI <- data.frame()

for(g in levels(NMDS.bc.SI$group)){
  df_ell.bc.SI <- rbind(df_ell.bc.SI, cbind(as.data.frame(with(NMDS.bc.SI[NMDS.bc.SI$group==g,],
                                                               veganCovEllipse(ord.bc.SI[[g]]$cov,ord.bc.SI[[g]]$center))),group=g))
}

#https://www.rpubs.com/RGrieger/545184


NMDS.bc.SI$Month <- as.factor(NMDS.bc.SI$Month)



pSIbiplot <- ggplot(NMDS.bc.SI, aes(NMDS1.bc, NMDS2.bc))+
  geom_point(stat = "identity", aes(color=group), size=1)+
  geom_path(data=df_ell.bc.SI, aes(x=delta13c.notlipidcorrected, y=delta15n, colour=group), size=1.5) +
  scale_fill_manual(values=POMsizesclrs) +
  labs(x=expression(paste(delta^{13}, "C (\u2030)")), 
       y=expression(paste(delta^{15}, "N (\u2030)")))+ 
  scale_color_manual(values=POMsizesclrs) + theme_bw() + 
  scale_y_continuous(breaks = c(6,9,12,15)) + 
  theme(legend.position = "none", 
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title.y = element_text(vjust=-1),
        axis.title.x = element_text(vjust=1))
pSIbiplot

 # png("ZS SI biplot 20251202.png", width=120, height=120, units="mm", res=300)
 # pSIbiplot
 # dev.off()





# Indicator value FAs -----------------------------------------------------

# Zooplankton taxa 
taxa.names <- as.factor(zoopData.sm$Species)

# Indicator species
(iva <- indval(Q20.matrix.abund, taxa.names, numitr = 10000))

# Correct the p-values for multiple testing:
pval.adj <- p.adjust(iva$pval, method = "holm")

gr <- iva$maxcls[pval.adj <= 0.05]
iv <- iva$indcls[pval.adj <= 0.05]
pv <- iva$pval[pval.adj <= 0.05]
fr <- apply(Q20.matrix.abund > 0, 2, sum)[pval.adj <= 0.05]
fidg <- data.frame(
  group = gr,
  indval = iv,
  pvalue = pv,
  freq = fr
)
fidg <- fidg[order(fidg$group, -fidg$indval), ]
fidg



# POM size classes
sizes <- as.factor(POM.sizes.wSI$Pore.Size)

# Indicator species
(iva <- indval(POM.matrix.abund, sizes, numitr = 10000))

# Correct the p-values for multiple testing:
pval.adj <- p.adjust(iva$pval, method = "holm")

gr <- iva$maxcls[pval.adj <= 0.05]
iv <- iva$indcls[pval.adj <= 0.05]
pv <- iva$pval[pval.adj <= 0.05]
fr <- apply(POM.matrix.abund > 0, 2, sum)[pval.adj <= 0.05]
fidg <- data.frame(
  group = gr,
  indval = iv,
  pvalue = pv,
  freq = fr
)
fidg <- fidg[order(fidg$group, -fidg$indval), ]
fidg




# Trophic Position --------------------------------------------------------

Zoop.SI <- SI.all %>% filter(Species!= "POM" & Species!= "z.250" & Species!= "March paraeuchaeta eggs")

# Calculate Trophic Position from Bulk POM baseline

# Need to calculate 14 day running average of POM Del15N 
QU39POM.sm <- QU39POM.ZSperiod %>% dplyr::select(Date, delta15n)
QU39POM.sm <- QU39POM.sm[complete.cases(QU39POM.sm$delta15n),]


# For Zoopsprint dates, calculate 
ZSdates <- c("2021-05-20", "2021-07-22", "2021-09-02", "2021-11-04", "2022-01-15", "2022-03-30")
ZSdates <- as.Date(ZSdates, "%Y-%m-%d")


length(QU39POM.sm[QU39POM.sm$Date >= (ZSdates[i]-14) & 
                    QU39POM.sm$Date<=ZSdates[i], ])

# Calculate 14 day moving average of POM DelN15
POM.dataframe <- data.frame(Date = ZSdates)

for(i in 1:length(ZSdates)){
  POM.dataframe$D15NmeanPOM14[i] <- mean(QU39POM.sm$delta15n[QU39POM.sm$Date >= (ZSdates[i]-7) & 
                                                                       QU39POM.sm$Date<=ZSdates[i] ], na.rm = T)
  POM.dataframe$n.D15NPOM14[i] <- length(QU39POM.sm$delta15n[QU39POM.sm$Date >= (ZSdates[i]-7) & 
                                                                       QU39POM.sm$Date<=ZSdates[i] ])
}

POM.dataframe$Month <- format(as.Date(POM.dataframe$Date), "%m")
POM.dataframe$Month <- as.numeric(POM.dataframe$Month)
POM.dataframe$Month[5:6] <- c(13, 15)

colnames(POM.dataframe)[2] <- "delta15n.POM"


SI.nonacid.forTP <- full_join(Zoop.SI, POM.dataframe, by="Month")


SI.nonacid.forTP$TrophicPosition <- ((SI.nonacid.forTP$delta15n - SI.nonacid.forTP$delta15n.POM) / 3.4) +1

SI.nonacid.forTP <- SI.nonacid.forTP %>% select(-Date.x) %>% 
  rename(Date = Date.y)

# Mean TP each Species 

TP.means <- SI.nonacid.forTP %>% group_by(Month, Date, Species) %>% 
  summarise(TrophicPosition = mean(TrophicPosition))






#  Using >20 POM only 
SI.POM$Month[SI.POM$Month==1] <- 13
SI.POM$Month <- as.numeric(SI.POM$Month)

SI.POM.micro <- SI.POM %>% filter(Pore.Size == "20" )

# Calculate means for each taxa in each month
SI.POM.micro.means <- SI.POM.micro %>% dplyr::select(Month, Date, delta15n) %>% 
  group_by(Month) %>% summarise_all(mean)

SI.POM.micro.means <- SI.POM.micro.means %>% rename(delta15n.POM = delta15n)
SI.POM.micro.means <- SI.POM.micro.means %>%  dplyr::select(Month, delta15n.POM) 


SI.nonacid.forTP <- full_join(Zoop.SI, SI.POM.micro.means, by="Month")
SI.nonacid.forTP$TrophicPosition <- ((SI.nonacid.forTP$delta15n - SI.nonacid.forTP$delta15n.POM) / 3.4) +1

TP.Micro.means <- SI.nonacid.forTP %>% group_by(Month, Species) %>% 
  summarise(TrophicPosition = mean(TrophicPosition))

TP.Micro.means <- TP.Micro.means %>%  mutate(Date = case_when(Month == 5 ~ "2021-05-20", 
                                                              Month == 7 ~ "2021-07-22", 
                                                              Month == 9 ~ "2021-09-02", 
                                                              Month == 11 ~ "2021-11-04", 
                                                              Month == 13 ~ "2022-01-15", 
                                                              Month == 15 ~ "2022-03-30"))
TP.Micro.means$Date <- as.Date(TP.Micro.means$Date)





#  Using >0.7 POM only 
SI.POM.pico <- SI.POM %>% filter(Pore.Size == "0.7" )

# Calculate means for each taxa in each month
SI.POM.pico.grouped <- SI.POM.pico %>% dplyr::select(Month, Date, delta15n) %>% 
  group_by(Month)
SI.POM.pico.means <- summarise_all(SI.POM.pico.grouped, mean)

SI.POM.pico.means <- SI.POM.pico.means %>% rename(delta15n.POM = delta15n)
SI.POM.pico.means <- SI.POM.pico.means %>%  dplyr::select(Month, delta15n.POM) 


SI.nonacid.forTP <- full_join(Zoop.SI, SI.POM.pico.means, by="Month")
SI.nonacid.forTP$TrophicPosition <- ((SI.nonacid.forTP$delta15n - SI.nonacid.forTP$delta15n.POM) / 3.4) +1

TP.pico.means <- SI.nonacid.forTP %>% group_by(Month, Species) %>% 
  summarise(TrophicPosition = mean(TrophicPosition))

TP.pico.means <- TP.pico.means %>%  mutate(Date = case_when(Month == 5 ~ "2021-05-20", 
                                                              Month == 7 ~ "2021-07-22", 
                                                              Month == 9 ~ "2021-09-02", 
                                                              Month == 11 ~ "2021-11-04", 
                                                              Month == 13 ~ "2022-01-15", 
                                                              Month == 15 ~ "2022-03-30"))
TP.pico.means$Date <- as.Date(TP.pico.means$Date)







#  Using nano POM only 
SI.POM.nano <- SI.POM %>% filter(Pore.Size == "3" )

# Calculate means for each taxa in each month
SI.POM.nano.grouped <- SI.POM.nano %>% dplyr::select(Month, Date, delta15n) %>% 
  group_by(Month)
SI.POM.nano.means <- summarise_all(SI.POM.nano.grouped, mean)

SI.POM.nano.means <- SI.POM.nano.means %>% rename(delta15n.POM = delta15n)
SI.POM.nano.means <- SI.POM.nano.means %>%  dplyr::select(Month, delta15n.POM) 


SI.nonacid.forTP <- full_join(Zoop.SI, SI.POM.nano.means, by="Month")
SI.nonacid.forTP$TrophicPosition <- ((SI.nonacid.forTP$delta15n - SI.nonacid.forTP$delta15n.POM) / 3.4) +1

TP.nano.means <- SI.nonacid.forTP %>% group_by(Month, Species) %>% 
  summarise(TrophicPosition = mean(TrophicPosition))

TP.nano.means <- TP.nano.means %>%  mutate(Date = case_when(Month == 5 ~ "2021-05-20", 
                                                            Month == 7 ~ "2021-07-22", 
                                                            Month == 9 ~ "2021-09-02", 
                                                            Month == 11 ~ "2021-11-04", 
                                                            Month == 13 ~ "2022-01-15", 
                                                            Month == 15 ~ "2022-03-30"))
TP.nano.means$Date <- as.Date(TP.nano.means$Date)


ggplot(TP.nano.means, aes(x= Month, y=TrophicPosition, color=Species, fill = Species)) + geom_point(size=3) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values = clrs) + 
  geom_line(linewidth=1)  + 
  theme_bw() + labs(title="nano-POM baseline")



#  Using 250-um zoop baseline 

SI.sizefracs.ZS <- SI.all %>% filter(Species == "z.250")

# Calculate 28 day moving average of zoops DelN15
z250.dataframe <- data.frame(Date = ZSdates)

for(i in 1:length(ZSdates)){
  z250.dataframe$D15NmeanPOM14[i] <- mean(SI.sizefracs.ZS$delta15n[SI.sizefracs.ZS$Date >= (ZSdates[i]-28) & 
                                                                     SI.sizefracs.ZS$Date<=ZSdates[i] ], na.rm = T)
  z250.dataframe$n.D15NPOM14[i] <- length(SI.sizefracs.ZS$delta15n[SI.sizefracs.ZS$Date >= (ZSdates[i]-28) & 
                                                                     SI.sizefracs.ZS$Date<=ZSdates[i] ])
}

z250.dataframe$Month <- format(as.Date(z250.dataframe$Date), "%m")
z250.dataframe$Month <- as.numeric(z250.dataframe$Month)
z250.dataframe$Month[5:6] <- c(13, 15)


colnames(z250.dataframe)[2] <- "delta15n.z250"



SI.nonacid.forTP.z250 <- full_join(Zoop.SI, z250.dataframe, by="Month")

SI.nonacid.forTP.z250$TrophicPosition <- ((SI.nonacid.forTP.z250$delta15n - SI.nonacid.forTP.z250$delta15n.z250) / 3.4) + 2

SI.nonacid.forTP.z250 <- SI.nonacid.forTP.z250 %>% select(-Date.x) %>% 
  rename(Date = Date.y)

# Mean TP each Species 
TP.means.z250 <- SI.nonacid.forTP.z250 %>% group_by(Month, Species) %>% 
  summarise(TrophicPosition = mean(TrophicPosition))

TP.means.z250 <- TP.means.z250 %>%  mutate(Date = case_when(Month == 5 ~ "2021-05-20", 
                                                            Month == 7 ~ "2021-07-22", 
                                                            Month == 9 ~ "2021-09-02", 
                                                            Month == 11 ~ "2021-11-04", 
                                                            Month == 13 ~ "2022-01-15", 
                                                            Month == 15 ~ "2022-03-30"))
TP.means.z250$Date <- as.Date(TP.means.z250$Date)







# *Fig. 3 -------------------------------------------------


  
pTP.bulk <- ggplot(TP.means, aes(x= Date, y=TrophicPosition, color=Species, fill = Species)) + geom_point(size=3) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values = clrs) + 
  # scale_shape_manual(values = species.shapes) + 
  geom_line(linewidth=1) +  theme_bw() + 
  scale_y_continuous(limits = c(0.75,3.5), expand = c(0,0),
                     breaks = c(1,1.5,2,2.5,3)) + 
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May\n2021", "Jul", "Sep", 
                          "Nov", "Jan\n2022", "Mar")) +
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        plot.margin = margin(5.5, 1, 5.5, 1),
        axis.title.x = element_blank())+ 
  geom_hline(yintercept=2) + 
  annotate("text",x=as.Date("2021-10-25"), y=3.4, label="bulk POM baseline")

pTP.bulk

pTP.micro <- ggplot(TP.Micro.means, aes(x= Date, y=TrophicPosition, color=Species, fill = Species)) + geom_point(size=3) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values = clrs) + 
  # scale_shape_manual(values = species.shapes) + 
  geom_line(linewidth=1)  + 
  scale_y_continuous(limits = c(0.75,3.5), expand = c(0,0),
                     breaks = c(1,2,3)) + 
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May\n2021", "Jul", "Sep", 
                          "Nov", "Jan\n2022", "Mar")) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin = margin(5.5, 1, 5.5, 1),
        axis.title.x = element_blank())+ 
  geom_hline(yintercept=2) + 
  annotate("text",x=as.Date("2021-10-25"), y=3.4, label="micro-POM baseline")
pTP.micro

pTP.pico <- ggplot(TP.pico.means, aes(x= Date, y=TrophicPosition, color=Species, fill = Species)) + geom_point(size=3) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values = clrs) + 
  # scale_shape_manual(values = species.shapes) + 
  geom_line(linewidth=1)  + 
  scale_y_continuous(limits = c(0.75,3.5), expand = c(0,0),
                     breaks = c(1,2,3)) + 
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May\n2021", "Jul", "Sep", 
                          "Nov", "Jan\n2022", "Mar")) +
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin = margin(5.5, 1, 5.5, 1),
        axis.title.x = element_blank()) + 
  geom_hline(yintercept=2) + 
  annotate("text",x=as.Date("2021-10-25"), y=3.4, label="pico-POM baseline")
pTP.pico

pTP.z250 <- ggplot(TP.means.z250, aes(x= Date, y=TrophicPosition, color=Species)) + geom_point(size=3) + 
  scale_color_manual(values = clrs) + 
  scale_fill_manual(values = clrs) + 
  # scale_shape_manual(values = species.shapes) + 
  scale_y_continuous(limits = c(0.75,3.5), expand = c(0,0),
                     breaks = c(1,2,3)) + 
  geom_line(linewidth=1) + theme_bw() + 
  scale_x_date(limits = as.Date(c("2021-05-01", "2022-04-10")), expand = c(0,0),
               breaks = as.Date(c("2021-05-01", "2021-07-01",
                                  "2021-09-01", "2021-11-01", "2022-01-01", 
                                  "2022-03-01")),
               labels = c("May\n2021", "Jul", "Sep", 
                          "Nov", "Jan\n2022", "Mar")) +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin = margin(5.5, 1, 5.5, 1),
        axis.title.x = element_blank()) + 
  geom_hline(yintercept=2) + 
  annotate("text",x=as.Date("2021-10-25"), y=3.4, label=expression(paste("250-",mu,"m zoop baseline")))
pTP.z250


prow <- plot_grid(pTP.bulk, pTP.micro, 
                  pTP.pico, pTP.z250,
                  rel_widths = c(1.1, 1, 1, 1),
                  # rel_heights = c(0.9, 1),
                  nrow = 1)
prow



# png("20250718 Trophic position plot.png", width=240, height=90, units="mm", res=300)
# prow
# dev.off()





# Cluster analysis --------------------------------------------------------


Q20.abund.brayFA <- vegdist(Q20.matrix.abund, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
# Ward's clustering method 
Q20.abundant.clust.ward<-hclust(Q20.abund.brayFA, method="ward.D2") 
k=6
# k=4
clusters <- cutree(Q20.abundant.clust.ward, k = k)

plot(Q20.abundant.clust.ward,labels = FALSE, hang= -1)


zoopData.sm$cluster <- as.factor(clusters)


NMDS.bc<-data.frame(NMDS1.bc=eco.nmds.bc$points[,1],NMDS2.bc=eco.nmds.bc$points[,2], 
                    group=taxa.names,
                    Cluster=zoopData.sm$cluster)


ggplot(NMDS.bc, aes(NMDS1.bc, NMDS2.bc))+ theme_bw() + 
  geom_point(stat = "identity", aes(color=Cluster), size=2) + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.01) +  # add the species labels
  scale_color_manual(values = c("#DF536B", "black","#CD0BBC" ,  "#28E2E5","#61D04F","#2297E6"))



palette(clrs) 
ordiplot(eco.nmds.bc,type="n", xlim=c(-1,1))
ordihull(eco.nmds.bc, groups=clusters, draw="polygon",label=F)
points(NMDS.bc$NMDS2.bc ~ NMDS.bc$NMDS1.bc, col=taxa.names, pch=16, cex=1.2)

arrows(0,0, y1=species.dataframe$NMDS.2, x1=species.dataframe$NMDS.1)
text(y=species.dataframe$NMDS.2, x=species.dataframe$NMDS.1, labels=species.dataframe$FA)




# Thinking more about the optimal number of clusters 



# simprof test on FA data

simprofTest<- simprof(Q20.matrix.abund, num.expected=1000, num.simulated=999,
                      method.cluster="ward.D2", method.distance="braycurtis", 
                      method.transform="identity", alpha=0.01,
                      sample.orientation="row", const=0,
                      silent=F, increment=100,
                      undef.zero=TRUE, warn.braycurtis=TRUE)

# Graph the result
pl.color <- simprof.plot(simprofTest)
summary(simprofTest)
# suggests seven groups



# 1) Silhouette plots
# Choose and rename the dendrogram ("hclust" object)
hc <- Q20.abundant.clust.ward

# Plot average silhouette widths (using Ward clustering) for all
# partitions except for the trivial partitions (k = 1 or k = n)
Si <- numeric(nrow(Q20.matrix.abund))
for (k in 2:(nrow(Q20.matrix.abund) - 1)){
  sil <- silhouette(cutree(hc, k = k), Q20.abund.brayFA)
  Si[k] <- summary(sil)$avg.width
}
k.best <- which.max(Si)
plot(
  1:nrow(Q20.matrix.abund),
  Si,
  type = "h",
  main = "Silhouette-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "Average silhouette width"
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best,
       max(Si),
       pch = 16,
       col = "red",
       cex = 1.5)

# Ward: optimum=5



# *Fig. S8 ------------------------------------------------

substring(zoopData.sm$Species, 1, 6)
format(as.Date(zoopData.sm$Date), "%m")
zoopData.sm$Month <- format(as.Date(zoopData.sm$Date), "%m")

r.names <-  paste(substring(zoopData.sm$Species, 1, 6), format(as.Date(zoopData.sm$Date), "%m"),  sep = " ")
r.nums <- sample(100:300, 149, replace = F)
r.names2 <- paste(r.names, r.nums, sep = " ")

row.names(Q20.matrix.abund) <- r.names2
 
Q20.abund.brayFA <- vegdist(Q20.matrix.abund, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
 
# Ward's clustering method 
Q20.abundant.clust.ward<-hclust(Q20.abund.brayFA, method="ward.D2")


dend <- as.dendrogram(Q20.abundant.clust.ward)


taxa.names <- as.factor(zoopData.sm$Species)
taxa.nums <- as.factor(taxa.names)
taxa.nums <- as.numeric(taxa.nums)


colors_to_use <- taxa.nums[order.dendrogram(dend)]
# 
labels_colors(dend) <- colors_to_use



# png("Dendrogram 20250707.png", width=125, height=500, units="mm", res=300)
par(mar=c(0,1,0,6))
dend %>%
  set("labels_cex", 0.9 )%>% set("branches_k_color", k=6, value = c("#CD0BBC", "#2297E6","#61D04F", "black", "#28E2E5", "#DF536B"),alpha=0.30) %>%
  plot(horiz=T)
# dev.off()






# *Fig. 4A --------------------------------------------------------


# Calculate mean SI for each taxa in each month
colnames(Zoop.SI)
Zoop.SI.grouped <- Zoop.SI %>% dplyr::select(Species, Month, delta13c, delta15n, C_N) %>% 
  group_by(Species, Month)

SI.nonacid.means <- summarise_all(Zoop.SI.grouped, mean)

# Project SI vectors onto fatty acid ordination 
# SI.nonacid.means from stable isotopes 20210127.R 
SI.nonacid.means$Month <- as.character(SI.nonacid.means$Month)
SI.nonacid.means$Month[SI.nonacid.means$Month=="5"] <- "05"
SI.nonacid.means$Month[SI.nonacid.means$Month=="7"] <- "07"
SI.nonacid.means$Month[SI.nonacid.means$Month=="9"] <- "09"
SI.nonacid.means$Month <- as.factor(SI.nonacid.means$Month)


zoopData.wSI <- right_join(SI.nonacid.means, zoopData.sm)


# C18.3andC18.4
isotopes <- zoopData.wSI %>% ungroup() %>% 
  select(delta13c, delta15n, DHA.EPA, Bacteria_no15, C16PUFA,
         percent.MUFA, percent.PUFA, percent.SFA, Copes.C20C.1C22.1, C18.3andC18.4)
colnames(isotopes)[c(1:2,4)] <- c("del13C", "del15N", "BAFA")
taxa.names <- as.factor(zoopData.wSI$Species)



Q20.matrix.abund <- zoopData.wSI[,abundant.z]
#create matrix for NMDS calculation:
species.matrix.sm <- as.matrix(Q20.matrix.abund)
#change diet dataframe into a matrix
class(species.matrix.sm) <- "numeric"
#make sure your numbers are treated as numbers

# set.seed(91203)
set.seed(91203)

eco.nmds.bc <- metaMDS(species.matrix.sm, distance="bray",labels=taxa.names, trymax = 100, autotransform = FALSE)

# BioEnv
(fit <- envfit(eco.nmds.bc, isotopes, perm = 999, na.rm = T))
# all significant except C18PUFAs 
summary(fit)
plot(eco.nmds.bc)
plot(fit)





NMDS.bc<-data.frame(NMDS1.bc=eco.nmds.bc$points[,1],NMDS2.bc=eco.nmds.bc$points[,2], 
                    group=taxa.names,  month=zoopData.wSI$Month)

plot(-NMDS.bc$NMDS2.bc ~ NMDS.bc$NMDS1.bc, col=taxa.names, pch=16, cex=1)


#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(eco.nmds.bc)$sites)
#add 'season' column as before
data.scores$taxa.names = taxa.names

en_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit)

data.scores$cluster <- zoopData.wSI$cluster

grp.a <- data.scores[data.scores$cluster == 1, ][chull(data.scores[data.scores$cluster == 
                                                                     1, c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$cluster == 2, ][chull(data.scores[data.scores$cluster == 
                                                                     2, c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.c <- data.scores[data.scores$cluster == 3, ][chull(data.scores[data.scores$cluster == 
                                                                     3, c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.d <- data.scores[data.scores$cluster == 4, ][chull(data.scores[data.scores$cluster == 
                                                                     4, c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.e <- data.scores[data.scores$cluster == 5, ][chull(data.scores[data.scores$cluster == 
                                                                     5, c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.f <- data.scores[data.scores$cluster == 6, ][chull(data.scores[data.scores$cluster == 
                                                                     6, c("NMDS1", "NMDS2")]), ]  # hull values for grp B

hull.data <- rbind(grp.a, grp.b, grp.c, grp.d, grp.e, grp.f)  #combine grp.a and grp.b
hull.data


data.scores$months <- (NMDS.bc$month)



gg = ggplot(data = data.scores, aes(x = NMDS1, y = -NMDS2)) +
  # add the convex hulls
  geom_polygon(data=hull.data,aes(x=NMDS1,y=-NMDS2,fill=as.factor(cluster),group=cluster),alpha=0.30) +
  scale_fill_manual(values = c("black","#DF536B", "#CD0BBC" ,  "#28E2E5","#61D04F","#2297E6")) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = -NMDS2), 
               data = en_coord_cont[-10,], size =1, alpha = 0.5, colour = "black") +
  # geom_text(data = en_coord_cont, aes(x = NMDS1*1.1, y = -NMDS2*1.1), colour = "grey30", 
  #           fontface = "bold", label = row.names(en_coord_cont)) + 
  geom_point(data = data.scores, aes(colour = taxa.names, shape=as.factor(months)), size = 2) + 
  # scale_y_continuous(limits = c(-0.6,0.9)) +
  # scale_x_continuous(limits = c(-0.65,0.6)) +
  scale_shape_manual(values = c(16,17,15,3,7,8)) +
  scale_colour_manual(values = clrs) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "none") + 
  labs(y="NMDS2") +
  annotate("text", x=en_coord_cont[c(1),1]*1.3, y=en_coord_cont[c(1),2]*-1.3, 
           label=expression(paste(delta^{13}, "C")),  size=3.5) + 
  annotate("text", x=en_coord_cont[c(2),1]*1.2, y=en_coord_cont[c(2),2]*-1, 
           label=expression(paste(delta^{15}, "N")),  size=3.5) + 
  annotate("text", x=en_coord_cont[c(3),1]*1.15, y=en_coord_cont[c(3),2]*-1.15, 
           label="DHA:EPA",  size=3.5) + 
  annotate("text", x=en_coord_cont[c(4),1]*2, y=en_coord_cont[c(4),2]*-1.2, 
           label="Bacterial FA",  size=3.5) + 
  annotate("text", x=en_coord_cont[c(5),1]*1.1, y=en_coord_cont[c(5),2]*-1.13, 
           label="C16-PUFAs",  size=3.5) + 
  annotate("text", x=en_coord_cont[c(6),1]*1.2, y=en_coord_cont[c(6),2]*-1.1, 
           label="% MUFA",  size=3.5) + 
  annotate("text", x=en_coord_cont[c(7),1]*1.2, y=en_coord_cont[c(7),2]*-1.11, 
           label="% PUFA",  size=3.5) +
  annotate("text", x=en_coord_cont[c(8),1]*1.17, y=en_coord_cont[c(8),2]*-1.1, 
           label="% SFA",  size=3.5) +
  annotate("text", x=en_coord_cont[c(9),1]*0.9, y=en_coord_cont[c(9),2]*-1.07, 
           label="C20,22-MUFAs",  size=3.5) 


gg

# png("NMDS w projected vars 20250504.png", width=112, height=112, units="mm", res=300)
# gg
# dev.off()




ggLegend = ggplot(data = data.scores, aes(x = NMDS1, y = -NMDS2)) +
  # add the convex hulls
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = -NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "black") +
  geom_text(data = en_coord_cont, aes(x = NMDS1*1.1, y = -NMDS2*1.1), colour = "grey30",
            fontface = "bold", label = row.names(en_coord_cont)) +
  geom_point(data = data.scores, aes(colour = taxa.names, shape=as.factor(months)), size = 2) + 
  scale_shape_manual(values = c(16,17,15,3,7,8), 
                     labels = c("May", "July", "Sept", "Nov", "Jan", "Mar")) +
  scale_colour_manual(values = clrs, 
                      labels = c("Calanus spp.", "C. limacina", "C. challengeri", 
                                 "E. bungii", "E. pacifica", "L. helicina",
                                 "M. pacifica", "P. elongata", "P. abyssalis", "T. pacifica")) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(4, 'mm'), #change legend key size
        legend.spacing.y = unit(2.0, "mm"),
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8, face = "italic"),
        legend.key = element_rect(fill = NA)) 

ggLegend

# png("Cluster plot legend 20250204.png", width=112, height=112, units="mm", res=300)
# ggLegend
# dev.off()





# SIMPER analysis ---------------------------------------------------------




# This is looking at differences in the clustering, which is based on the untransformed data

sim <- simper(Q20.matrix.abund, clusters,permutations=1000)

sim
simper.summ <- summary(sim)
simper.summ






# * Table S1 --------------------------------------------------------------

n.POM.SI <- CrossTable(SI.POM$Pore.Size, SI.POM$Month, prop.t=F, prop.r=F, prop.c=F, prop.chisq=F)

n.Zoop.SI <- CrossTable(Zoop.SI$Species, Zoop.SI$Month, prop.t=F, prop.r=F, prop.c=F, prop.chisq=F)

n.POM.FA <- CrossTable(POM.sizes$Pore.Size, POM.sizes$Month, prop.t=F, prop.r=F, prop.c=F, prop.chisq=F)

n.Zoop.FA <- CrossTable(zoopData.sm$Species, zoopData.sm$Month, prop.t=F, prop.r=F, prop.c=F, prop.chisq=F)




