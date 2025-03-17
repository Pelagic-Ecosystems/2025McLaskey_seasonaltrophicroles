
library(tidyverse)
library(vegan)
library(Polychrome)
library(gmodels)
library(adespatial)


# Load full FA dataset 
Q21data <- read.csv("raw_data/Q21 joined w location 2023-11-21.csv")
str(Q21data)
Q21data$Date <- as.Date(Q21data$Date, "%Y-%m-%d")
Q21data <- Q21data %>% rename(Species = Species.AKM)


# Q22 data 
Q22data <- read.csv("raw_data/Q22 joined w location 2023-05-10.csv")
str(Q22data)
Q22data$Date <- as.Date(Q22data$Date, "%Y-%m-%d")
# One tube was run twice and not removed in Q22 compiling yet
Q22data <- Q22data[-(grep("Q22-024-2", Q22data$file)),]
# Clione sample ZFA-338 / tube Q22-107 is an obvious outlier and very low chromatogram. 
Q22data <- Q22data[Q22data$Sample!="ZFA-338",]


# Q23 data 
# Q23data <- Q23.joined.wLocation
Q23data <- read.csv("raw_data/Q23 joined w location 2024-10-03.csv")
str(Q23data)
Q23data$Date <- as.Date(Q23data$Date, "%Y-%m-%d")


allFAdata.1 <- full_join(Q21data, Q22data)
allFAdata <- full_join(allFAdata.1, (Q23data))


zoopsprintData <- allFAdata %>% filter(Project =="Zoopsprint")

table(zoopsprintData$Date)

zoopsprintData$Month <- format(as.Date(zoopsprintData$Date), "%m")
zoopsprintData$Month <- as.factor(zoopsprintData$Month)

zoopsprintData$Species[zoopsprintData$Species==""] <- "POM"



##

# Only most abundant species
zoopsprintData <- zoopsprintData %>% filter(Species !="Aequorea victoria" & 
                                     Species !="Larvaceans" & 
                                     Species !="chaetegnath" & 
                                     Species !="Rhynconerella" & 
                                     Species !="Parasagitta elegans" & 
                                     Species !="Scina spp." & 
                                     Species !="Hyperia" & 
                                     Species !="Beroe" & 
                                     Species !="Clytia gregaria" &
                                     Species !="Thysanoessa spinifera"&
                                     Species !="Thysanoessa longipes"&
                                     Species !="Tomopteris pacificus"&
                                     Species !="Tomopteris"&
                                     Species !="Tomopteris septentrionalis"&
                                     Species !="Neocalanus plumchrus"&
                                     Species !="Munida quadraspina"&
                                     Species !="Disconchoecia large"&
                                     Species !="Disconchoecia elegans"&
                                     Species !="Disconchoecia small"&
                                     Species !="Ostracods"&
                                     Species !="Pasiphaea pacifica"&
                                     Species !="500" &
                                     # Species !="250" &
                                     Species !="Pleurobrachia")



# remove the couple of Calanus samples in March that I separated by sex/stage
zoopsprintData$Stage

zoopsprintData <- zoopsprintData[-grep("all male", zoopsprintData$Stage),]
zoopsprintData <- zoopsprintData[-grep("all female", zoopsprintData$Stage),]
# Also need to remove Nov Calanus that contained C3-C4 (stage lost when combined) 
# Also removing Nov Paraeuchaeta - it was a very small sample and a strong outlier 
zoopsprintData <- zoopsprintData[zoopsprintData$Sample != "ZFA-219" & zoopsprintData$Sample != "ZFA-231" & zoopsprintData$Sample != "ZFA-234",]




# Calculate Prop ID for peak area 
zoopsprintData$PropID <- zoopsprintData$SumFA / zoopsprintData$Total.FA_G_AREA 

# Make factor for season 
zoopsprintData <- zoopsprintData %>% mutate(season = ifelse( Month %in% c("03","04","05"), "spring" ,
                                                       ifelse(Month %in% c("06", "07", "08"), "summer", 
                                                              ifelse(Month %in% c("09", "10", "11"), "fall",
                                                                     "winter"))))

zoopsprintData <- zoopsprintData %>% select(Sample:Site.ID, Type:Species, 
                                            SumFA_mg.g, SumFA_mg.gWW, 
                                            SumFA_ug.L:iso.17.0_PERCENT,
                                            Size.Fraction, Month:season, 
                                            file, Tube.ID)



# Calculate Markers -------------------------------------------------------

# 16:1n7/16:0
zoopsprintData$Ratio16.1 <- (zoopsprintData$C16.1n.7_PERCENT/ zoopsprintData$C16.0_PERCENT)

# Diatom FAs / Flagellate FAs
zoopsprintData$Diatom.Flag <- ((zoopsprintData$C16.1n.7_PERCENT + zoopsprintData$C20.5n.3_PERCENT +  zoopsprintData$C16.2n.4_PERCENT + zoopsprintData$C16.3n.4_PERCENT + zoopsprintData$C16.4n.1_PERCENT) 
                            / (zoopsprintData$C22.6n.3_PERCENT + zoopsprintData$C18.3n.3_PERCENT + zoopsprintData$C18.3n.6_PERCENT + zoopsprintData$C18.4n.3_PERCENT))

zoopsprintData$C16PUFA <- rowSums(zoopsprintData[,c("C16.2n.4_PERCENT", "C16.3n.4_PERCENT", "C16.4n.1_PERCENT")], na.rm=TRUE)


# DHA/EPA : Dinos/Diatoms
zoopsprintData$DHA.EPA <- (zoopsprintData$C22.6n.3_PERCENT/ zoopsprintData$C20.5n.3_PERCENT)


# 15:0 and 17:0 : Bacteria
zoopsprintData$Bacteria_all <- (zoopsprintData$C15.0_PERCENT + zoopsprintData$C17.0_PERCENT + 
                               zoopsprintData$iso.15.0_PERCENT + zoopsprintData$iso.17.0_PERCENT + zoopsprintData$ant.15.0_PERCENT)

zoopsprintData$Bacteria_no15 <- (zoopsprintData$C17.0_PERCENT + 
                                  zoopsprintData$iso.15.0_PERCENT + zoopsprintData$iso.17.0_PERCENT + zoopsprintData$ant.15.0_PERCENT)

# 18:1n-9/18:1n-7 Carnivory or omnivory
zoopsprintData$Carn18.1N9_n7 <- (zoopsprintData$C18.1n.9c_PERCENT/ zoopsprintData$C18.1n.7_PERCENT)


# Make 18:3n-3_18:4n-4 (pico-Chl)
zoopsprintData$C18.3andC18.4 <- zoopsprintData$C18.3n.3_PERCENT + zoopsprintData$C18.4n.3_PERCENT


# Make Copepood marker C20, C22 MUFAs
zoopsprintData$Copes.C20C.1C22.1 <- zoopsprintData$C20.1n.11_PERCENT + zoopsprintData$C20.1n.9_PERCENT +
  zoopsprintData$C22.1n.11_PERCENT + zoopsprintData$C22.1n.9_PERCENT



# Calculate proportions for each main class

percent.SFA <- vector(length=nrow(zoopsprintData))
for(i in 1:nrow(zoopsprintData)){
  percent.SFA[i]=sum(
    zoopsprintData$C10.0_PERCENT[i], 
    zoopsprintData$C11.0_PERCENT[i], 
    zoopsprintData$C12.0_PERCENT[i], 
    zoopsprintData$C13.0_PERCENT[i], 
    zoopsprintData$C14.0_PERCENT[i], 
    zoopsprintData$C15.0_PERCENT[i], 
    zoopsprintData$C16.0_PERCENT[i], 
    zoopsprintData$C17.0_PERCENT[i], 
    zoopsprintData$C18.0_PERCENT[i], 
    zoopsprintData$C20.0_PERCENT[i],
    zoopsprintData$C22.0_PERCENT[i],
    zoopsprintData$C23.0_PERCENT[i],
    zoopsprintData$C24.0_PERCENT[i],
    na.rm=TRUE)
}
zoopsprintData$percent.SFA <- percent.SFA


percent.MUFA <- vector(length=nrow(zoopsprintData))
for(i in 1:nrow(zoopsprintData)){
  percent.MUFA[i]=sum(zoopsprintData$C14.1_PERCENT[i], 
                      zoopsprintData$C16.1n.7_PERCENT[i], 
                      zoopsprintData$C18.1n.7_PERCENT[i], 
                      zoopsprintData$C18.1n.9c_PERCENT[i], 
                      zoopsprintData$C20.1n.9_PERCENT[i],
                      zoopsprintData$C22.1n.9_PERCENT[i],
                      zoopsprintData$C24.1n.9_PERCENT[i],na.rm=TRUE)
}
zoopsprintData$percent.MUFA <- percent.MUFA


percent.PUFA <- vector(length=nrow(zoopsprintData))
for(i in 1:nrow(zoopsprintData)){
  percent.PUFA[i]=sum(zoopsprintData$C18.2n.6c_PERCENT[i], 
                      zoopsprintData$C18.3n.3_PERCENT[i], 
                      zoopsprintData$C18.3n.6_PERCENT[i], 
                      zoopsprintData$C18.4n.3_PERCENT[i],
                      zoopsprintData$Diatom.III_PERCENT[i],
                      zoopsprintData$C18.1n.12_PERCENT[i], 
                      zoopsprintData$C20.3n.3_PERCENT[i], 
                      zoopsprintData$C20.4n.6_PERCENT[i], 
                      zoopsprintData$C20.5n.3_PERCENT[i],
                      zoopsprintData$C22.2n.6_PERCENT[i], 
                      zoopsprintData$C22.4n.6_PERCENT[i], 
                      zoopsprintData$C22.5n.3_PERCENT[i], 
                      zoopsprintData$C22.5n.6._PERCENT[i],
                      zoopsprintData$C22.6n.3_PERCENT[i], na.rm=TRUE)
} 
zoopsprintData$percent.PUFA <- percent.PUFA



# PUFA/SFA carnivory index
zoopsprintData$PUFA.SFA <- zoopsprintData$percent.PUFA/zoopsprintData$percent.SFA



percent.n3_PUFA <- vector(length=nrow(zoopsprintData))
for(i in 1:nrow(zoopsprintData)){
  percent.n3_PUFA[i]=sum(zoopsprintData$C16.2n.4_PERCENT[i],
                         zoopsprintData$C16.3n.4_PERCENT[i],
                         zoopsprintData$C16.4n.1_PERCENT[i],
                         zoopsprintData$C18.3n.3_PERCENT[i], 
                         zoopsprintData$C18.4n.3_PERCENT[i],
                         zoopsprintData$C20.3n.3_PERCENT[i], 
                         zoopsprintData$C20.4n.3_PERCENT[i], 
                         zoopsprintData$C20.5n.3_PERCENT[i],
                         zoopsprintData$C22.5n.3_PERCENT[i], 
                         zoopsprintData$C22.6n.3_PERCENT[i], na.rm=TRUE)
} 
zoopsprintData$percent.n3_PUFA <- percent.n3_PUFA

percent.n6_PUFA <- vector(length=nrow(zoopsprintData))
for(i in 1:nrow(zoopsprintData)){
  percent.n6_PUFA[i]=sum(zoopsprintData$C18.2n.6_PERCENT[i], 
                         zoopsprintData$C18.3n.6_PERCENT[i], 
                         zoopsprintData$C20.3n.6_PERCENT[i], 
                         zoopsprintData$C20.4n.6_PERCENT[i], 
                         zoopsprintData$C22.2n.6_PERCENT[i], 
                         zoopsprintData$C22.4n.6_PERCENT[i], 
                         zoopsprintData$C22.5n.6_PERCENT[i], na.rm=TRUE)
} 
zoopsprintData$percent.n6_PUFA <- percent.n6_PUFA


# Bacterial FATM without 15:0
zoopsprintData$Bacteria_no15 <- (zoopsprintData$C17.0_PERCENT + 
                                   zoopsprintData$iso.15.0_PERCENT + 
                                  zoopsprintData$iso.17.0_PERCENT + 
                                   zoopsprintData$ant.15.0_PERCENT)



zoopsprintData$SumFA_ug.L[is.infinite(zoopsprintData$SumFA_ug.L)] <- NA





# POM data assembly -------------------------------------------------------

# POM data ----------------------------------------------------------------


POMData <- zoopsprintData %>% filter(Type =="POM")


# Load POM metadata
POMmetadata <- read.csv("metadata/2021_seasonal_POM samples.csv")
POMmetadata$Date <- as.Date(POMmetadata$Date, "%m/%d/%y")
colnames(POMmetadata)[2] <- "Sample"
colnames(POMmetadata)[8] <- "Site.ID"

POMData.full <- full_join(POMmetadata, POMData)
# one tube broke, not on data. 
# one sample from McKertcher creek not on metadata sheet
POMData.full <- POMData.full[!is.na(POMData.full$file),]
POMData.full <- POMData.full[!is.na(POMData.full$Size),]


# Take out replicates from days I sampled both pooled depths and 5 m only.
POMData.full <- POMData.full %>% filter(!(Depth =="5" & Date == as.Date("2021-07-20") & Size == "bulk") &
                                  !(Depth =="5" & Date == as.Date("2021-07-21")& Size == "bulk") &
                                  !(Depth =="5" & Date == as.Date("2021-09-01")& Size == "bulk") &
                                  !(Depth =="5" & Date == as.Date("2021-09-02")& Size == "bulk") &
                                  !(Depth =="5" & Date == as.Date("2022-03-30")& Size == "bulk"))

POMData.full$Species <- "POM"
POMData.full <- POMData.full %>% filter(Date < as.Date("2022-03-30"))

# # 2024-04-17 removing Q22-004 and Q22-055 
# # lots of strange peaks in chromatograms w very low sample signal - don't trust 
POMData.full <- POMData.full %>% filter(Tube.ID != "Q22-004" & Tube.ID != "Q22-055")



# full FA dataset ---------------------------------------------------------

ZoopData <- zoopsprintData %>% filter(Type =="Zoop")

zsAll <- full_join(POMData.full, ZoopData)

zsAll <- zsAll %>% select(-file, -Tube.ID, -Vol..mL., -sample.type, -Notes, 
                          -Run.Date, -Sample, -PropID)
zsAll <- zsAll %>% rename(Line.Out.Depth = Depth, 
                          Pore.Size = Size)
 # write.csv(zsAll, "zoopsprintData all 2025-03-17.csv", row.names = F)





# Stable Isotopes ---------------------------------------------------------


# Load collection metadata
zoopMetadata <- read.csv("metadata/Zoopsprint_FAsamples.csv", stringsAsFactors = F)
colnames(zoopMetadata)[2] <- "Date"
zoopMetadata$Date <- as.Date(zoopMetadata$Date, "%m/%d/%Y")
zoopMetadata$Species

POMmetadata <- read.csv("raw_data/2024-01-18_131849_HakaiData_poms.csv")
colnames(POMmetadata)
POMmetadata <- POMmetadata %>% dplyr::select(Date, Hakai.ID, Site.ID, Line.Out.Depth, Volume..ml., Weight.Before..mg.)
colnames(POMmetadata)[2:3] <- c("Sample", "Site")

POMsizes <- read.csv("metadata/Zoopsprint POM SI sizes.csv")
POMsizes <- POMsizes %>% rename(Site = Site.ID)
POMmetadata <- full_join(POMmetadata, POMsizes)
POMmetadata$Date <- as.Date(POMmetadata$Date)
POMmetadata$Species <- "POM"

str(zoopMetadata)
colnames(zoopMetadata)[1:2] <- c("Sample", "Date")


hakaizoop <- read.csv("raw_data/2023-04-25_145549_HakaiData_zooplankton_isotope.csv")
hakaizoop$Date <- as.Date(hakaizoop$Date, "%Y-%m-%d")
hakaizoop <- hakaizoop %>% dplyr::select(Date, Site.ID, Hakai.ID, Size.Fraction..um.)
colnames(hakaizoop)[2:3] <- c("Site", "Sample")

zoopMetadata <- full_join(zoopMetadata, hakaizoop)


# Load SI sent by Windsor on 2022-11-29

SIdata1 <- read.csv("raw_data/SI data MASTER Anna McLaskey 2022_20230327.csv", na="")

colnames(SIdata1)
colnames(SIdata1)[1:6] <- c("Sample.SI", "delta13c", "percent.C", 
                            "delta15n", "percent.N", "C_N")

SIdata1 <- SIdata1 %>% dplyr::select(-starts_with('X'))
SIdata1 <- SIdata1 %>% filter(Sample.SI != "")


# Only ZFA samples 
SIdata.zoops <- SIdata1 %>% filter(grepl("ZFA",Sample.SI))
SIdata.zoops$Sample <- substring(SIdata.zoops$Sample.SI, 1,7)

# Only POM / zoop size frac samples 
SIdata.POM <- SIdata1 %>% filter(grepl("QF",Sample.SI))
SIdata.POM$Sample <- substring(SIdata.POM$Sample.SI, 1,6)

SIdata2 <- full_join(SIdata.zoops, SIdata.POM)



# Load SI sent by Windsor on 2023-10-29
SIdata3 <- read.csv("raw_data/SI data MASTER Brian Hunt 2023_20231019.csv", na="")

colnames(SIdata3)
colnames(SIdata3)[1:6] <- c("Sample", "delta13c", "percent.C", 
                            "delta15n", "percent.N", "C_N")

# Only POM / zoop size frac samples 
SIdata3.QF <- SIdata3 %>% filter(grepl("QF",Sample))

SIdata4 <- full_join(SIdata2, SIdata3.QF)

# QF8700 was accidentally assigned QF9700 in results back from Windsor
SIdata4$Sample[SIdata4$Sample=="QF9700"] <- "QF8700"

# Fix ZFA105. sample name
SIdata4$Sample[SIdata4$Sample=="ZFA105."] <- "ZFA-105"



metadataAll <- full_join(zoopMetadata, POMmetadata)
# ZFA-292 was assigned ZFA-292.1 in metadata
metadataAll$Sample[metadataAll$Sample=="ZFA-292.1"] <- "ZFA-292"


metadataAll$Month <- format(as.Date(metadataAll$Date), "%m")
metadataAll$Month <- as.factor((metadataAll$Month))

SI.combined <- right_join(metadataAll, SIdata4)

SI.combined <- SI.combined %>% filter(Date > as.Date("2021-01-01"))
# Separate out 250 size fraction only
SI.combined <- SI.combined %>% filter(is.na(Size.Fraction..um.) | Size.Fraction..um. == "250")

SI.combined$Species[SI.combined$Size.Fraction..um.=="250"] <- "z.250"
SI.combined$Species



# acidified or not? 
# Some of the non acidified were left the same, not marked with N 
SI.combined$Sample.number <- substring(SI.combined$Sample.SI, 4,15)
SI.combined$Acidified <- ifelse(grepl("A", SI.combined$Sample.number), "Yes",
                                "Non")
SI.combined$Acidified  <- as.factor(SI.combined$Acidified )


SI.combined$Month <- as.numeric(as.character(SI.combined$Month))
SI.combined$Month[SI.combined$Month==1] <- 13
SI.combined$Month[SI.combined$Month==3] <- 15


SI.combined.nonacid <- SI.combined %>% filter(Acidified == "Non") %>% 
  rename(delta13c.not = delta13c, 
         delta15n.not = delta15n)
colnames(SI.combined.nonacid)



# Correction factor for Limacina percentC and delta13C 
# SI.combined.Limacina <- SI.combined %>% filter(Species == "Limacina helicina") %>% 
#   select(Date, delta13c, Acidified) 
# 
# SI.combined.Limacina.agg <- SI.combined.Limacina %>%   aggregate(by = list(SI.combined.Limacina$Date,
#                                                                            SI.combined.Limacina$Acidified),
#             FUN = mean, na.rm=TRUE)
#   
# SI.combined.Limacina.wide <- SI.combined.Limacina.agg %>% pivot_wider(values_from = delta13c, names_from = Group.2)
# 
# SI.combined.Limacina.wide$C13.offset <- SI.combined.Limacina.wide$Yes - SI.combined.Limacina.wide$Non
# mean(SI.combined.Limacina.wide$C13.offset, na.rm = T)
# 
# SI.combined.Limacina.wide$percent.C.offset <- SI.combined.Limacina.wide$Yes - SI.combined.Limacina.wide$Non
# mean(SI.combined.Limacina.wide$percent.C.offset, na.rm = T)


# Have not made offset for C:N ratio yet 
offset.delta13C <- 3.870198 # 3.870198 before=3.882431
offset.percent.C <- 11.8205 # 11.8205 before=10.78028

SI.combined.nonacid.Limacina <- SI.combined.nonacid %>% filter(Species == "Limacina helicina")  
SI.combined.nonacid.Limacina$delta13c.not <- SI.combined.nonacid.Limacina$delta13c.not - offset.delta13C
SI.combined.nonacid.Limacina$percent.C <- SI.combined.nonacid.Limacina$percent.C - offset.percent.C


SI.combined.nonacid.NOLimacina <- SI.combined.nonacid %>% filter(Species != "Limacina helicina")  
SI.combined.nonacid.corrected <- full_join(SI.combined.nonacid.Limacina, SI.combined.nonacid.NOLimacina)

# Calculate C:N (from non acidified samples, need to use corrected Limacina)
SI.combined.nonacid.corrected$C_N <- (SI.combined.nonacid.corrected$percent.C * 12.01) / (SI.combined.nonacid.corrected$percent.N * 14.007)


# Only lipid correct Zooplankton NOT POMs

SI.combined.nonacid.corrected.noPOM <- SI.combined.nonacid.corrected %>% filter(Species != "POM")
SI.combined.nonacid.corrected.POM <- SI.combined.nonacid.corrected %>% filter(Species == "POM")
colnames(SI.combined.nonacid.corrected.POM)
SI.combined.nonacid.corrected.POM <- SI.combined.nonacid.corrected.POM %>% mutate(delta13c =  delta13c.not)



# lipid correction on Del13C following regression for SoG copepods
# δ13C_corrected = δ13C_bulk + (0.38 * C:N_bulk) - 1.85  (El-Sabaawi et al. 2009)
SI.combined.nonacid.corrected.noPOM <- SI.combined.nonacid.corrected.noPOM %>% mutate(delta13c =  delta13c.not + (0.38 * C_N) - 1.85 )
SI.combined.nonacid.corrected.noPOM <- SI.combined.nonacid.corrected.noPOM %>% rename(delta13c.notlipidcorrected = delta13c.not)

SI.combined.nonacid.corrected.final <- full_join(SI.combined.nonacid.corrected.noPOM, SI.combined.nonacid.corrected.POM)




# QU39 POMs isotope data compiling and filtering 

# Load POM Isotope Data 
POMisotopeDataQU39 <- read.csv("raw_data/2024-01-18_131849_HakaiData_poms.csv", stringsAsFactors = FALSE, na=c("", "NA"))

# POMisotopeDataQU39 <- read.csv("Covariates raw data/2020-12-03_110551_HakaiData_poms.csv", stringsAsFactors = FALSE, na="")
POMisotopeDataQU39$Date <- as.Date(POMisotopeDataQU39$Date, "%Y-%m-%d")

# use carbon isotopes from acidified samples, nitrogen and C:N from not acidified samples
# Need acidified and non-acidified ug.C, ug.N, corr_delta15n, corr_delta13c
# Also need Volume filtered, C.Flag, N.Flag, date, Line.Out.Depth
smallPOMisoData <- POMisotopeDataQU39 %>% select(Date, Line.Out.Depth, Acidified, Volume..ml.,  corr_delta15n, corr_delta13c, ug.C, ug.N, C.Flag, N.Flag) 

# Do all samples have acidified or not?
table(is.na(smallPOMisoData$Acidified))
# Four are missing Acidified. Remove them for now. 
smallPOMisoData <- smallPOMisoData %>%filter(!is.na(Acidified))

# Check for sample dates+depths where both samples were not Acidified, or both were, so they don't have unique combo of keys.
smallPOMisoData[duplicated(smallPOMisoData[,1:3]),]
smallPOMisoData[duplicated(smallPOMisoData[,1:3], fromLast = T),]

# There are four dates+depths with duplicate samples
# As well as the zoopsprint samples - will need to pull in separately 

# There is only one sample from 0 or 5 m, looks like both samples really weren't Acidified. Just going to remove duplicates for now
smallPOMisoData.new <- smallPOMisoData %>% distinct(Date, Line.Out.Depth, Acidified, .keep_all = T)


# Check quality flags
table(smallPOMisoData.new$C.Flag)
table(is.na(smallPOMisoData.new$C.Flag))
table(smallPOMisoData.new$N.Flag)
table(is.na(smallPOMisoData.new$N.Flag))


# Look at measurements flagged as too low biomass

# Calculate elemental concentrations based on volume filtered
smallPOMisoData.new <- smallPOMisoData.new %>% 
  mutate(C_ug.L =  (ug.C / (Volume..ml./1000))) %>% 
  mutate(N_ug.L =  (ug.N / (Volume..ml./1000))) 


# spread the acidified/non-acidified samples 
# Remove Vol filtered, ug.C, and ug.N which aren't meaningful to spread
smallPOMisoData.new <- smallPOMisoData.new %>% select(-c(Volume..ml., ug.C, ug.N))
smallPOMisoData.new_long <- pivot_wider(smallPOMisoData.new, names_from = Acidified, values_from = c(corr_delta15n, corr_delta13c, C_ug.L, N_ug.L))

# I don't find true=acidified and false=not acidified to be very helpful. Change col names 
# CHECK TO MAKE SURE COLS ARE STILL IN SAME ORDER
colnames(smallPOMisoData.new_long)[5:12]<- c("delta15n.notacid", "delta15n.acidified","delta13c.notacid",
                                             "delta13c.acidified",  "C_ug.L.notacid", "C_ug.L.acidified", "N_ug.L.notacid",
                                             "N_ug.L.acidified")


# Acidified samples are either good for carbon or aren't good (i.e. flags on nitrogen don't matter)
# Non-Acidified samples  could be good for just nitrogen, or also C_N (i.e. flags on C and N both matter)

# Split the dataset in two to QC acidfied and not-Acidified separately
smallPOMisoData.acid <- smallPOMisoData.new %>% filter(Acidified=="TRUE")
smallPOMisoData.notacid <- smallPOMisoData.new %>% filter(Acidified=="FALSE")

# Only the carbon matters out of the Acidified sample
smallPOMisoData.acid.noCflag <- smallPOMisoData.acid %>% filter(is.na(C.Flag))
# both N and C matter out of the unAcidified sample, although I may want to separate these later to maximize the amount of data I have 
smallPOMisoData.notacid.noflags <- smallPOMisoData.notacid %>% filter(is.na(N.Flag) & is.na(C.Flag))
# I need them both for unAcidified samples because I'm calculated C:N, but only need the C.Flag for Acidified samples
smallPOMisoData.QC <- rbind(smallPOMisoData.acid.noCflag, smallPOMisoData.notacid.noflags)

# Before spreading the samples I need to remove N.Flag because it is different between the rows now
smallPOMisoData.QC <- smallPOMisoData.QC %>% select(-N.Flag)

# Calculate C:N 
smallPOMisoData.QC <- smallPOMisoData.QC %>% mutate(C_N =  ((C_ug.L*12.01) / (N_ug.L*14.007)) )


# spread the acidified/non-acidified samples 
# Remove Vol filtered, ug.C, and ug.N which aren't meaningful to spread
smallPOMisoData.QC <- smallPOMisoData.QC %>% select(-c(C.Flag))
# smallPOMisoData.QC <- smallPOMisoData.QC %>% select(-c(Volume..ml., ug.C, ug.N, C.Flag))
smallPOMisoData_long <- pivot_wider(smallPOMisoData.QC, names_from = Acidified, values_from = c(C_N, corr_delta15n, corr_delta13c, C_ug.L, N_ug.L))

# I don't find true=acidified and false=not acidified to be very helpful. Change col names 
# CHECK TO MAKE SURE COLS ARE STILL IN SAME ORDER
colnames(smallPOMisoData_long)[3:12]<- c("C_N_acidified", "C_N_notacid",  "delta15n.acidified","delta15n.notacid",
                                         "delta13c.acidified", "delta13c.notacid", "C_ug.L.acidified", "C_ug.L.notacid", 
                                         "N_ug.L.acidified", "N_ug.L.notacid")

# Look at distribution of sampling depths
table(smallPOMisoData_long$Line.Out.Depth)

## filter out 0 and 5 m only
# Average the zero and 5 m samples 
POMisoData_0and5m <- smallPOMisoData_long %>% filter(Line.Out.Depth==5 | Line.Out.Depth==0)
POMisoData_0and5m.agg = aggregate(POMisoData_0and5m,
                                  by = list(POMisoData_0and5m$Date),
                                  FUN = mean, na.rm=TRUE)

# Isotope data all dates
POMisoData_all <- POMisoData_0and5m.agg %>% select(-Group.1) %>% 
  filter(Date >= "2021-01-01" & Date < "2022-03-30")

POMisoData_all$Species <- "POM"
POMisoData_all$Pore.Size <- "bulk"
POMisoData_all <- POMisoData_all %>% rename(C_N = C_N_notacid) %>% 
  rename(delta15n.not = delta15n.notacid) %>% 
  rename(delta13c = delta13c.acidified) %>% 
  select(Date, Species, Pore.Size, Line.Out.Depth, C_N, delta15n.not, delta13c, delta13c.notacid, C_ug.L.acidified, C_ug.L.notacid, N_ug.L.notacid) 

SI.combined.nonacid.corrected.final$Pore.Size <- as.character(SI.combined.nonacid.corrected.final$Pore.Size)

# remove 3-30-2022 POM samples that were filtered >8 hrs after collection and were very different than samples from the day before
SI.combined.nonacid.corrected.final <- SI.combined.nonacid.corrected.final %>% filter(Species !="POM" | Date != as.Date("2022-03-30"))

ZS.SIdata <- full_join(SI.combined.nonacid.corrected.final, POMisoData_all)

ZS.SIdata <- ZS.SIdata %>% select(Date:Species, Site, 
                                  Line.Out.Depth, Pore.Size, Month, C_N, 
                                  delta13c, delta15n.not, C_ug.L.acidified)
ZS.SIdata <- ZS.SIdata %>% rename(Depth = Line.Out.Depth,
                                  Site.ID = Site, 
                                  delta15n = delta15n.not)
  
 # write.csv(ZS.SIdata, "full_zoopsprint SI 20250317.csv", row.names = F)









# Chlorophyll -------------------------------------------------------------

# Load full dataset 
QU39Chl <- read.csv("raw_data/2023-01-18_141158_HakaiData_chlorophyll.csv")

# QU39Chl <- read.csv("G:/My Drive/UBC Work/Work/Research/2021 Zoopsprint/Analysis/2023-07-28_115109_HakaiData_chlorophyll_2018toCurrent.csv", stringsAsFactors = F, na="")
str(QU39Chl)
QU39Chl$Date <- as.Date(QU39Chl$Date, "%Y-%m-%d")

QU39Chl <- QU39Chl %>% filter(Site.ID == "QU39")
QU39Chl <- QU39Chl %>% filter(Survey != "ZOOPSPRINT")
QU39Chl <- QU39Chl %>% filter(Replicate.Number < 2)
# QU39Chl <- QU39Chl %>% filter(Survey == "ZOOPSPRINT")

colnames(QU39Chl)

chlSumm.FA2 <- QU39Chl %>% select(Date, Line.Out.Depth, 
                                  Hakai.ID, Filter.Type, Before.Acid, After.Acid, Acid.Flag,
                                  Chla, Chla.Flag, Phaeo, Phaeo.Flag,
                                  Quality.Level, Quality.Log)


# Use Chla_final to avoid negative values. Those below detection limit have Chla.Flag= "DL"

# Looking at the acid ratio is an important QC step
plot(chlSumm.FA2$Before.Acid ~ chlSumm.FA2$After.Acid)
abline(0,1)
# There are a couple outliers

table(chlSumm.FA2$Chla.Flag)
# no flags


# Although the units are not actually indicated, they are already ug/L
chl.noflags <- chlSumm.FA2 %>% mutate(chla_ug.L = Chla)


# Will need to spread filter type. Make reduced dataset
colnames(chl.noflags)
chl.noflags <- select(chl.noflags, Date, Line.Out.Depth, Filter.Type, chla_ug.L)
# Check for duplicates
dupes <-chl.noflags[duplicated(chl.noflags[1:3]),]
chl.noflags[duplicated(chl.noflags[1:3], fromLast = T),]
# No issues

# Spread
chl.noflags.wide <- chl.noflags %>% pivot_wider(names_from = Filter.Type, values_from = chla_ug.L) 
colnames(chl.noflags.wide)[3:6] <- paste("chl", colnames(chl.noflags.wide)[3:6], sep = "_")

colnames(chl.noflags.wide)[c(5,6)] <- c("chl_GF.F", "chl_Bulk.GF.F")


ggplot(chl.noflags.wide, aes(x=Date, y=chl_Bulk.GF.F, color=as.factor(Line.Out.Depth))) + 
  geom_line(size=1.2) + theme_bw() +
  scale_color_manual(values=clrs)


# Find mean of 0 and 5 m 
chl.noflags.0_5m <- chl.noflags.wide %>% 
  filter(Line.Out.Depth<=5)

# Find mean by filter size
chl.noflags.agg.mean = aggregate(chl.noflags.0_5m,
                                       by = list(chl.noflags.0_5m$Date),
                                       FUN = mean, na.rm=TRUE)
# remove extra cols
chl_all <- chl.noflags.agg.mean[-c(1)] 
# Now I have my file of all dates 


# Calculate SumChl for days that have all three size fracs
chl_all <- chl_all %>%  mutate(SumChl = (chl_20um + chl_3um + chl_GF.F))
chl_all <- chl_all %>% select(-chl_Bulk.GF.F)


# write.csv(chl_all, "zoopsprint Chl data 2025-03-17.csv", row.names = F)








