




####REQUIRED PACKAGES###########################################################
packages <- c(
  "readxl",     #reads excel files
  "ggplot2",    #makes better plots
  "dplyr",      #helps make data tables more editable
  "tidyr",      
  "GGally",     #additional plots etc for ggplot
  "corrplot",   #makes correlation plots
  "ggforce",    #for convex hulls
  "ggridges",   #to make ridge line plots
  "forcats",    #reorder factors to make plots more meaningful
  "ggrepel" ,   #avoids overlapping text
  "rstatix",    #genereate the W statistic for pairwise wilcox tests
  "viridis",    #color blind friendly color palettes
  "viridisLite",
  "ggnewscale", #so I can have two scale fills 
  "patchwork",  #allows you to put two figures side by side
  "stringr"     #str_extract
)

#installs missing packages from this list
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)

if(length(to_install) > 0){
  install.packages(to_install, dependencies = TRUE)
  }

#load the packages
lapply(packages, library, character.only = TRUE)

################################################################################

####IMPORT DATA#################################################################
M1_data <- read.csv("https://raw.githubusercontent.com/jstamer95/Phaleron-Serial-Sampling/refs/heads/main/M1_data.csv")

BI_comparative_data <-read.csv("https://raw.githubusercontent.com/jstamer95/Phaleron-Serial-Sampling/refs/heads/main/isoarch_data.csv")

model_diets <- read.csv("https://raw.githubusercontent.com/jstamer95/Phaleron-Serial-Sampling/refs/heads/main/model_diets.csv") #already offset to diet

M1_data <- M1_data %>%
  rename(
    `wt.%.C` = `wt...C`,
    `wt.%.N` = `wt...N`,
    `C:N` = `C.N`
  )

BI_comparative_data <- BI_comparative_data %>%
  rename(
    `δ13C (‰vpdb)` = "δ13C...vpdb.",
    `δ15N (‰air)` = "δ15N...air.",
  )

##Phaleron adult comparative data
##see Hannigan et al (forthcoming) for these data
Phaleron_adult_diet <- data.frame(
  burial.number = c("5_252", "IV_646", "IV_692", "IV_694"),
  d15N.Air = c(10.4, 9.4, 10.6, 9.7),
  d13C = c(-19.5, -19.7, -19.8, -19.3),
  group = "adult")

Phaleron_adult_summary <- Phaleron_adult_diet %>%
  summarise(
    across(c(d15N.Air, d13C), 
           list(
             mean = ~mean(. , na.rm = TRUE), 
             sd = ~ifelse(n() > 1, sd(., na.rm = TRUE), NA), 
             max = ~max(., na.rm = TRUE), 
             min = ~min(., na.rm = TRUE), 
             range_low = ~min(., na.rm = TRUE), 
             range_high = ~max(., na.rm = TRUE)
           )
    ), 
    .groups = "drop"
  )

adult.Phaleron.N <- Phaleron_adult_summary$d15N.Air_mean
adult.Phaleron.C <- Phaleron_adult_summary$d13C_mean
################################################################################

####COLLAGEN QUALITY CRITERIA:##################################################

##data clean up using standard CQC measures:

M1_data <- M1_data %>%
  mutate(stage = factor(stage)) %>%
  rename(`d13C` = "d13.vs.VPDB") %>%
  filter(`wt.%.C` > 5.19) %>% #going to remove low wt % C and N per Ambrose 1990 
  filter(`wt.%.N` > 1.79) %>%
  filter(`C:N` > 2.89 & `C:N` < 3.61) #Trimming C:N values according to DeNiro (1985) (2.9-3.6)
  
  
##Guiry and Szpak (2021) suggest plotting C:N ratio to d13N and d15N to look for correlation, indicator of humic contaminates
  
humic.con.d13C.regression <- lm(`d13C`~`C:N`, data = M1_data) #regression of d13C and C:N ratio
summary_humic_c <- summary(humic.con.d13C.regression) #storing the summary of the model
humic_C_r2 <- round(summary_humic_c$r.squared, 4) #getting the R squared value
pearson_test.C <- cor.test(M1_data$`C:N`, M1_data$`d13C`, use = "complete.obs", method = "pearson") #running a Pearson's R for d13C and C:N
pearson_r.C <- pearson_test.C$estimate  #storing Pearson's R 
pearson_p.C <- pearson_test.C$p.value  #storing Pearson's R p-value
  
r2_label.C <- paste0("R² = ", round(humic_C_r2, 3), #creating a data label for the graph with the R squared value, Pearson's R, and p-value
                       "\nPearson's r = ", round(pearson_r.C, 3),
                       "\np = ", round(pearson_p.C, 6))
  
humic.con.d13C <- ggplot(M1_data, aes(x= `C:N`, y= `d13C`)) + #graphing the relationship between d13C and C:N
  geom_point()+
  geom_smooth(method=lm, color= "red") +
  labs(
    x = "atomic C:N", 
    y = "δ13C vs VPDB",
    title = "Collagen Quality Criteria: C:N v δ13C", 
    subtitle = "C:N Trimmed to DeNiro (1985) (2.9-3.6)"
   ) +
  theme(
    plot.title = element_text(hjust = 0.5), 
    plot.subtitle = element_text(hjust = 0.5)
  )+
  annotate("text", x = Inf, y = Inf, label = r2_label.C , size = 5, hjust=1.5, vjust = 2.5, color = "green", parse = FALSE) #adding the data label 
plot(humic.con.d13C) #plotting the graph
  
humic.con.d15N.regression <- lm(`d15N.Air`~`C:N`, data= M1_data) #regression of d15N and C:N ratio
summary_humic_n <-summary(humic.con.d15N.regression) #storing the summary of that model
humic_N_r2 <- round(summary_humic_n$r.squared, 4) #storing the r squared value
pearson_test.N <- cor.test(M1_data$`C:N`, M1_data$`d15N.Air`, use = "complete.obs", method = "pearson") #running a Pearson's R 
pearson_r.N <- pearson_test.N$estimate  #storing the Pearson's r value
pearson_p.N <- pearson_test.N$p.value #storing the p-value
  
r2_label.N <- paste0("R² = ", round(humic_N_r2, 8), #creating a data label with the R squared, Pearson's R, and p-value
                      "\nPearson's r = ", round(pearson_r.N, 3),
                      "\np = ", round(pearson_p.N, 3))
humic.con.d15N <- ggplot(M1_data, aes(x= `C:N`, y= `d15N.Air`))+
  geom_point()+
  geom_smooth(method=lm, color= "red")+
  labs(
    x = "atomic C:N", 
    y = "δ15N v Air",
    title = "Collagen Quality Criteria: C:N v δ15N", 
    subtitle = "C:N Trimmed to DeNiro (1985) (2.9-3.6)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), 
    plot.subtitle = element_text(hjust = 0.5) 
  )+
  annotate("text", x = Inf, y = Inf, label = r2_label.N , size = 5, hjust=1.5, vjust = 2.5, color = "green", parse = FALSE)
plot(humic.con.d15N)
  
  #These plots indicate that there is likely some remaining humic contamination 
  #as Guiry and Szpak 2021 noted in their cod samples. 
  #however, the expected d13C values from other sites in Athens, including other 
  #Phaleron data, suggest that -18% to -20% values are possible. Table 4 of Guiry 
  #and Szpak 2021 suggests a conservative upper limit of 3.45 to 3.50, and a liberal
  #tolerance of 3.70 to 3.90 for these d13C values.
  #We will stick with the DeNiro 1985 range (2.9-3.6), noting that for samples could 
  #be altered by as much as -1.71 d13C and -0.08 N, but likley are under that threashhold.
  #We will try to avoid interpreting significant dietary difference ~+/-2% for d13c 
  #and ~+/-0.4% for d15N to account for both humic contamination as well as total 
  #analytical uncertainty (Szpak et al 2017). 
  
################################################################################

####WEANING CURVES##############################################################

#this is the set up for generating the weaning curves
#isolate the individuals who have more than 4 serial samples
M1_data_serial_only <- M1_data %>%
  filter(!(`burial.number` == "IV_202")) %>%
  filter(!(`burial.number` == "IV_220")) %>%
  filter(!(`burial.number` == "IV_380")) %>%
  filter(!(`burial.number` == "IV_400")) %>%
  filter(!(`burial.number` == "IX_826"))

#establish weaning age
M1_data_serial_only_weaning <- M1_data_serial_only %>%
  arrange(burial.number, MIDPOINT.serial.age) %>%  
  group_by(burial.number) %>%
  mutate(d15N.Air_diff = d15N.Air - lag(d15N.Air)) %>% #create a new column called "d15N.Air_diff" that is the difference between the first d15N value and next d15N value for that individual by age
  mutate(d15N.Air_diff = replace_na(d15N.Air_diff, 0)) %>%  #Replace NA with 0
  mutate(d15N.Air_total = cumsum(d15N.Air_diff)) %>% #also creates a new column that sums the total change over time for that individual called "d15N.Air_total"
  ungroup()

#WEANING AND POST-WEANING DIET
summary(M1_data$d15N.Air)
#filtered out the post weaning data
post_weaning <- M1_data %>%
  filter(stage == "post-weaning")
weaning <- M1_data %>%
  filter(stage == "weaning")


#M1 DATA SUMMARY GROUPED
M1_data_summary_grouped <- M1_data %>%
  group_by(burial.number) %>%
  summarise(
    across(c(d15N.Air, d13C, `C:N`), 
           list(
             mean = ~mean(. , na.rm = TRUE), 
             sd = ~ifelse(n() > 1, sd(., na.rm = TRUE), NA), 
             max = ~max(., na.rm = TRUE), 
             min = ~min(., na.rm = TRUE), 
             range_low = ~min(., na.rm = TRUE), 
             range_high = ~max(., na.rm = TRUE)
           )
    ), 
    .groups = "drop"
  )
M1_data_summary_grouped <-M1_data_summary_grouped %>%
  mutate(d13C.range = abs(d13C_max-d13C_min)) %>%
  mutate(d15N.range = abs(d15N.Air_max-d15N.Air_min)) %>%
  select(-c(d15N.Air_range_high, d15N.Air_range_low, d13C_range_high, d13C_range_low, d13C_sd, d15N.Air_sd, `C:N_sd`, `C:N_max`, `C:N_min`, `C:N_range_low`, `C:N_range_high`))



#WEANING AGE
weaning_ages <- data.frame(
  burial.number = c("5_102", "5_135", "IV_356", "IV_372", "IV_488", "IV_743"),
  weaning.age = c(2.9, 2.9, 2.55, 3.0, 2.7, 2.9)
)
weaning_age_summary<-summary(weaning_ages)

  ##SUPPLEMENTAL FIGURE 1: WEANING CURVES####
  #actually create the graphs
  #Find global axis limits
  x_limits <- range(M1_data_serial_only$MIDPOINT.serial.age, na.rm = TRUE)
  y_limits_N <- range(M1_data_serial_only$d15N.Air, na.rm = TRUE)
  y_limits_C <- range(M1_data_serial_only$d13C, na.rm = TRUE)
  
  #This is supplemental figure 1: it shows the serial sampling data for all individuals with more than 1 serial sample
  #it is not possible to get ggplot2 to graph two axes on the same plot, so what I did was generate each of these and overlay them in illustrator
  #there is a base R example below that generates the plot all together but it is not as pretty and I did not include it in the manuscript
  big_C_plot_scaled <- ggplot(M1_data_serial_only, aes(x = MIDPOINT.serial.age, y = d13C)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 24, color= "black", fill="orange", size=2) +
    labs(x =  "",
         y = expression(paste(delta^13, "C"[collagen], " (\u2030, VPDB)"))) +
    facet_wrap(~ burial.number) +
    theme_minimal()+
    ggtitle("Early Childhood Diet at Phaleron") +
    geom_hline(yintercept = adult.Phaleron.C, #the adult carbon values from Hannigan et al (forthcoming)
               linetype = "dashed", 
               color = "orange", linewidth = 2)  +
    scale_x_continuous(limits = x_limits) +  
    scale_y_continuous(limits = y_limits_C) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 16, color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),
      axis.title = element_text(size = 18, color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text.x = element_blank(),
      plot.title.position = "plot",
      axis.ticks = element_blank(),
      axis.line.y = element_line(color = "grey50"),
      axis.line.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
    )    
  print(big_C_plot_scaled)
  
  #same y scale but for N
  big_N_plot_scaled <- ggplot(M1_data_serial_only, aes(x = MIDPOINT.serial.age, y = d15N.Air)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 21, color= "black", fill="#69b3a2", size=2) + 
    labs(x = "", y = expression(paste(delta^15, "N"[collagen], " (\u2030, AIR)"))) +
    facet_wrap(~ burial.number) +
    theme_minimal() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = adult.Phaleron.N,
               linetype = c("dotted"), 
               color = c("#69b3a2"), linewidth = 1.2) +
    scale_x_continuous(limits = x_limits) +  
    scale_y_continuous(limits = y_limits_N, position = "right")  +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 16, color = "grey50"),
      axis.title = element_text(size = 18, color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text.x = element_blank(),
      plot.title.position = "plot",
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "grey50"),
      panel.grid = element_line(color = "#b4aea9"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linetype = "solid"),
      panel.grid.major.y = element_line(linetype = "solid")
    ) 
  

  print(big_N_plot_scaled)
  

  ##PLOT OF C AND N TOGETHER (THANKS SHANNON)
    #list of individuals with good serial sampling data
    good_data = c("5_102", "5_135", "5_218", "IV_063", "IV_225", "IV_281", "IV_356", "IV_372", "IV_488", "IV_560", "IV_734", "IV_743", "IV_782")
    par(mfrow = c(4, 4)) #parameters of a grid, 4X4
    
    #create a for loop that loops from 1 to n of good_data
    for (i in 1:length(good_data)) { #i is number of loop
      
      plot_data = M1_data_serial_only %>% 
        filter(burial.number == good_data[i]) %>% 
        arrange(MIDPOINT.serial.age) #filtering to the list from those with good data
      
      par(mar = c(5, 4, 3, 3) +0.5) #figure size
      
      plot(plot_data$MIDPOINT.serial.age, plot_data$d15N.Air, pch = 15, ces = 2, type = "b", 
           col = "#69b3a2", xlab = "Age", ylab = "", main = as.character(good_data[i]), cex.main = 2, 
           ylim = c(8, 15), xlim = c(0, 9), axes = FALSE, lwd = 2.5) #this is the nitrogen plot
      axis(1)
      axis(2, las= 1, col.axis = "#69b3a2") #setting conditions for right side y axis
      mtext("δ15N vs Air", side = 2, line = 3, col = "#69b3a2") #right side y axis data label
      abline(h = adult.Phaleron.N, lty = 2, col = "#69b3a2", lwd = 2.5) #this is the nitrogen adult data as a line
      
      par(new = TRUE) #adding new data to the same plot
      
      plot(plot_data$MIDPOINT.serial.age, plot_data$d13C, pch = 17, ces = 2, type = "b",  col = "orange", 
           xlab = "", ylab = "", axes = FALSE, bty= "n", 
           ylim = c(-21, -14), xlim = c(0, 9), lwd = 1.5) #carbon plots
      abline(h = adult.Phaleron.C, lty = 3, col = "orange", lwd = 2) 
      
      axis(4, las = 1, col.axis = "orange") #setting conditions for right side y axis
      mtext("δ13C vs VPDB", side = 4, line = 3, col = "orange") #right side y axis data label
      
    }
  
  ##FIGURE 3: INDIVIDUAL WEANING CURVES####
  #C and N are plotted individually for each individual. This has to be combined into figure 3 in illustrator because ggplot wont plot two axes
  
    ##5_102####
  x5_102 <- M1_data_serial_only_weaning %>% #weaning age 2.9, positive change is 0.10, less than our uncertainty in N
    filter(burial.number == "5_102")

  x5_102_N <- ggplot(x5_102, aes(x = MIDPOINT.serial.age, y = d15N.Air)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 21, color= "black", fill="#69b3a2", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    ggtitle("5_102") +
    geom_hline(yintercept = adult.Phaleron.N,
               linetype = c("dotted"), 
               color = c("#69b3a2"), linewidth = 2) +
    geom_vline(xintercept = 2.9, 
               linetype = "dashed",  
               color = "blue", linewidth = 1.2) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 16, color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),
      axis.title = element_text(size = 18, color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text.x = element_blank(),
      plot.title.position = "plot",
      axis.ticks = element_blank(),
      axis.line.y = element_line(color = "grey50"),
      axis.line.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
    ) +
    scale_x_continuous(limits = x_limits) +  
    scale_y_continuous(limits = y_limits_N)          
  print(x5_102_N)
  
  x5_102_C <- ggplot(x5_102, aes(x = MIDPOINT.serial.age, y = d13C)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 24, color= "black", fill="orange", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      legend.position = "none", 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.line.y = element_line(color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),)+
    ggtitle("") +
    geom_hline(yintercept = adult.Phaleron.C,
               linetype = c("dashed"), 
               color = c("orange"), linewidth = 2) +
    scale_x_continuous(limits = x_limits) + 
    scale_y_continuous(limits = y_limits_C, position = "right")         
  print(x5_102_C)
  
  

    ##5_135####
  x5_135 <-M1_data_serial_only_weaning %>% #weaning age 2.9, positive change is 0.37
    filter(burial.number == "5_135")
  
  x5_135_N <- ggplot(x5_135, aes(x = MIDPOINT.serial.age, y = d15N.Air)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 21, color= "black", fill="#69b3a2", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none")+
    ggtitle("5_135") +
    geom_hline(yintercept = adult.Phaleron.N,
               linetype = c("dotted"), 
               color = c("#69b3a2"), linewidth = 2) +
    geom_vline(xintercept = 2.9, 
               linetype = "dashed",  
               color = "blue", linewidth = 1.2) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 16, color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),
      axis.title = element_text(size = 18, color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text.x = element_blank(),
      plot.title.position = "plot",
      axis.ticks = element_blank(),
      axis.line.y = element_line(color = "grey50"),
      axis.line.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
    ) +
    scale_x_continuous(limits = x_limits) +  
    scale_y_continuous(limits = y_limits_N)    
  print(x5_135_N)
  
  x5_135_C <- ggplot(x5_135, aes(x = MIDPOINT.serial.age, y = d13C)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 24, color= "black", fill="orange", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      legend.position = "none", 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.line.y = element_line(color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),)+
    ggtitle("") +
    geom_hline(yintercept = adult.Phaleron.C,
               linetype = c("dashed"), 
               color = c("orange"), linewidth = 2) +
    scale_x_continuous(limits = x_limits) + 
    scale_y_continuous(limits = y_limits_C, position = "right")        
  print(x5_135_C)
  

    ##IV_356####
  IV_356 <-M1_data_serial_only_weaning %>% #weaning age 2.55, positive change is 1.06 which is greater than uncertainty, but they do continue to have a curve afterwards, maybe their age is actually, 2.55 with a change of 0.10, less than uncertainty
    filter(burial.number == "IV_356")

  IV_356_N <- ggplot(IV_356, aes(x = MIDPOINT.serial.age, y = d15N.Air)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 21, color= "black", fill="#69b3a2", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none")+
    ggtitle("IV_356") +
    geom_hline(yintercept = adult.Phaleron.N,
               linetype = c("dotted"), 
               color = c("#69b3a2"), linewidth = 2) +
    geom_vline(xintercept = 2.9, 
               linetype = "dashed",  
               color = "blue", linewidth = 1.2) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 16, color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),
      axis.title = element_text(size = 18, color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text.x = element_blank(),
      plot.title.position = "plot",
      axis.ticks = element_blank(),
      axis.line.y = element_line(color = "grey50"),
      axis.line.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
    ) +
    scale_x_continuous(limits = x_limits) +  
    scale_y_continuous(limits = y_limits_N)    
  print(IV_356_N)
  
  IV_356_C <- ggplot(IV_356, aes(x = MIDPOINT.serial.age, y = d13C)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 24, color= "black", fill="orange", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      legend.position = "none", 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.line.y = element_line(color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),)+
    ggtitle("") +
    geom_hline(yintercept = adult.Phaleron.C,
               linetype = c("dashed"), 
               color = c("orange"), linewidth = 2) +
    scale_x_continuous(limits = x_limits) + 
    scale_y_continuous(limits = y_limits_C, position = "right")         
  print(IV_356_C)
  
  
    ##IV_372####
  IV_372 <-M1_data_serial_only_weaning %>% #can weaning be estimated for this person? I think it should go at 3.0, but there is not 2% dip or an increase!
    filter(burial.number == "IV_372")

  IV_372_N <- ggplot(IV_372, aes(x = MIDPOINT.serial.age, y = d15N.Air)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 21, color= "black", fill="#69b3a2", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none")+
    ggtitle("IV_372") +
    geom_hline(yintercept = adult.Phaleron.N,
               linetype = c("dotted"), 
               color = c("#69b3a2"), linewidth = 2) +
    geom_vline(xintercept = 2.9, 
               linetype = "dashed",  
               color = "blue", linewidth = 1.2) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 16, color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),
      axis.title = element_text(size = 18, color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text.x = element_blank(),
      plot.title.position = "plot",
      axis.ticks = element_blank(),
      axis.line.y = element_line(color = "grey50"),
      axis.line.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
    ) +
    scale_x_continuous(limits = x_limits) +  
    scale_y_continuous(limits = y_limits_N)   
  print(IV_372_N)
  
  IV_372_C <- ggplot(IV_372, aes(x = MIDPOINT.serial.age, y = d13C)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 24, color= "black", fill="orange", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      legend.position = "none", 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.line.y = element_line(color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),)+
    ggtitle("") +
    geom_hline(yintercept = adult.Phaleron.C,
               linetype = c("dashed"), 
               color = c("orange"), linewidth = 2) +
    scale_x_continuous(limits = x_limits) + 
    scale_y_continuous(limits = y_limits_C, position = "right")         
  print(IV_372_C)
  
  
    ##IV_488####
  IV_488 <-M1_data_serial_only_weaning %>% #weaning age 2.7, 2% decrease and positive change greater than uncertainty, actually 2nd increase, but idk
    filter(burial.number == "IV_488")

  IV_488_N <- ggplot(IV_488, aes(x = MIDPOINT.serial.age, y = d15N.Air)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 21, color= "black", fill="#69b3a2", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none")+
    ggtitle("IV_488") +
    geom_hline(yintercept = adult.Phaleron.N,
               linetype = c("dotted"), 
               color = c("#69b3a2"), linewidth = 2) +
    geom_vline(xintercept = 2.9, 
               linetype = "dashed",  
               color = "blue", linewidth = 1.2) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 16, color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),
      axis.title = element_text(size = 18, color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text.x = element_blank(),
      plot.title.position = "plot",
      axis.ticks = element_blank(),
      axis.line.y = element_line(color = "grey50"),
      axis.line.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
    ) +
    scale_x_continuous(limits = x_limits) +  
    scale_y_continuous(limits = y_limits_N)   
  print(IV_488_N)
  
  IV_488_C <- ggplot(IV_488, aes(x = MIDPOINT.serial.age, y = d13C)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 24, color= "black", fill="orange", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      legend.position = "none", 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.line.y = element_line(color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),)+
    ggtitle("") +
    geom_hline(yintercept = adult.Phaleron.C,
               linetype = c("dashed"), 
               color = c("orange"), linewidth = 2) +
    scale_x_continuous(limits = x_limits) + 
    scale_y_continuous(limits = y_limits_C, position = "right")        
  print(IV_488_C)
  
  
    ##IV_743####
  IV_743 <-M1_data_serial_only_weaning %>% #weaning age 2.9, never goes below 2% total decrease, but increase is greater than uncertainty
    filter(burial.number == "IV_743")

  IV_743_N <- ggplot(IV_743, aes(x = MIDPOINT.serial.age, y = d15N.Air)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 21, color= "black", fill="#69b3a2", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none")+
    ggtitle("IV_743") +
    geom_hline(yintercept = adult.Phaleron.N,
               linetype = c("dotted"), 
               color = c("#69b3a2"), linewidth = 2) +
    geom_vline(xintercept = 2.9, 
               linetype = "dashed",  
               color = "blue", linewidth = 1.2) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 16, color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),
      axis.title = element_text(size = 18, color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text.x = element_blank(),
      plot.title.position = "plot",
      axis.ticks = element_blank(),
      axis.line.y = element_line(color = "grey50"),
      axis.line.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
    ) +
    scale_x_continuous(limits = x_limits) +  
    scale_y_continuous(limits = y_limits_N)   
  print(IV_743_N)
  
  IV_743_C <- ggplot(IV_743, aes(x = MIDPOINT.serial.age, y = d13C)) +
    geom_line(color = "grey", linewidth = 1.2) + 
    geom_point(shape = 24, color= "black", fill="orange", size=2) + 
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      legend.position = "none", 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.line.y = element_line(color = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "grey50"),
      axis.text = element_text(size = 16, color = "grey50"),)+
    ggtitle("") +
    geom_hline(yintercept = adult.Phaleron.C,
               linetype = c("dashed"), 
               color = c("orange"), linewidth = 2) +
    scale_x_continuous(limits = x_limits) + 
    scale_y_continuous(limits = y_limits_C, position = "right")         
  print(IV_743_C)

################################################################################

####SOCIAL AGE GROUP COMPARSONS#################################################

#INFANT DIET: all serial samples under 3 years, this number was chosen after weaning age was established for Phaleron
infant_diet <- M1_data %>%
  filter(`MIDPOINT.serial.age` < 2.99) %>%
  mutate(group = "infant")
  
infant_diet_summary_grouped <- infant_diet %>%
  group_by(burial.number) %>%
  summarise(
    across(c(d15N.Air, d13C), 
            list(
              mean = ~mean(. , na.rm = TRUE), 
              sd = ~ifelse(n() > 1, sd(., na.rm = TRUE), NA), 
              max = ~max(., na.rm = TRUE), 
              min = ~min(., na.rm = TRUE), 
              range_low = ~min(., na.rm = TRUE), 
              range_high = ~max(., na.rm = TRUE),
              n = ~sum(!is.na(.)) 
            )
    ), 
    .groups = "drop"
  )

infant_diet_summary_grouped <-infant_diet_summary_grouped %>%
  mutate(d13C.range = abs(d13C_max-d13C_min)) %>%
  mutate(d15N.range = abs(d15N.Air_max-d15N.Air_min)) %>%
  select(-c(d15N.Air_range_high, d15N.Air_range_low, d13C_range_high, d13C_range_low, d13C_sd, d15N.Air_sd))
  
#CHILD DIET
child_diet <- M1_data %>%
  filter(`MIDPOINT.serial.age` > 3) %>%
  mutate(group = "child")
  
child_diet_summary_grouped <- child_diet %>%
  group_by(burial.number) %>%
  summarise(
    across(c(d15N.Air, d13C), 
            list(
              mean = ~mean(. , na.rm = TRUE), 
              sd = ~ifelse(n() > 1, sd(., na.rm = TRUE), NA), 
              max = ~max(., na.rm = TRUE), 
              min = ~min(., na.rm = TRUE), 
              range_low = ~min(., na.rm = TRUE), 
              range_high = ~max(., na.rm = TRUE),
              n = ~sum(!is.na(.)) 
            )), .groups = "drop")

child_diet_summary_grouped <-child_diet_summary_grouped %>%
  mutate(d13C.range = abs(d13C_max-d13C_min)) %>%
  mutate(d15N.range = abs(d15N.Air_max-d15N.Air_min)) %>%
  select(-c(d15N.Air_range_high, d15N.Air_range_low, d13C_range_high, d13C_range_low, d13C_sd, d15N.Air_sd))
  
  
##PAIRWISE WILCOX TEST:
  
age_comparison <- infant_diet %>%
  full_join(Phaleron_adult_diet)
age_comparison <- age_comparison %>%
  full_join(child_diet)

age.comparison.N <- age_comparison %>%
  pairwise_wilcox_test(
    formula = d15N.Air ~ group,
    p.adjust.method = "holm"  
    )
age.comparison.N

age.comparison.C <- age_comparison %>%
  pairwise_wilcox_test(
    formula = d13C ~ group,   # replace d15N with your variable of interest
    p.adjust.method = "holm"  # you can also use "bonferroni", "fdr", etc.
  )
age.comparison.C
  ##only child and infant N is significantly different
  
  
  ##FIGURE 2: SOCIAL AGE COMPARISON####
convex_hulls_age <- age_comparison %>% #this determines the shape of the filled in area in the graph
  group_by(group) %>%
  slice(chull(d13C, d15N.Air))  
  
  shape_values <- c("child" = 24, "infant" = 21, "adult" = 22)
  
  age.comparison <- ggplot(age_comparison, aes(x = d13C, y = d15N.Air)) +
    geom_polygon(
      data = convex_hulls_age, 
      aes(fill = group, group = group), 
      alpha = 0.7, 
      color = NA
    ) +
    geom_point(
      aes(shape = group, fill = group), 
      size = 4,
      color = "black",   
      stroke = 0.8       
    ) +
    scale_shape_manual(values = shape_values) +
    scale_fill_viridis_d(option = "mako", alpha = 0.8) + 
    guides(fill = "none") +
    theme_minimal() +
    ggtitle("Infant, Child, and Adult Diet at Phaleron") +
    labs(x =  expression(paste(delta^13, "C"[collagen], " (\u2030, VPDB)")),
         y = expression(paste(delta^15, "N"[collagen], " (\u2030, AIR)")), 
         shape = "Age Group", 
         fill = "Age Group") +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 16, colour = "grey50"),
      legend.title = element_text(size = 18, colour = "grey50"),
      axis.text = element_text(size = 16, colour = "grey50"),
      axis.title = element_text(size = 18, colour = "grey50"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", colour = "grey50"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, colour = "grey50"),
      plot.title.position = "plot",
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "grey50"),
      panel.grid = element_line(color = "#b4aea9"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linetype = "solid"),
      panel.grid.major.y = element_line(linetype = "solid"),
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
    )
  
  age.comparison

################################################################################
  
####MODELLING DIET OF INFANTS###################################################
  #filter the individual data for all data points under 3 years old, trying to see what the weaning diet would be
  IV_281 <- M1_data %>%
    filter(burial.number == "IV_281")
  IV_281_weaning <- IV_281 %>%
    filter(MIDPOINT.serial.age < 2.99) %>% # n = 4 
    mutate(group = "IV_281")
  IV_734 <- M1_data %>%
    filter(burial.number == "IV_734")
  IV_734_weaning <- IV_734%>%
    filter(MIDPOINT.serial.age < 2.99)%>% # n = 1
    mutate(group = "IV_734")
  weaning <- infant_diet %>%
    filter(!burial.number %in% c("IV_734", "IV_281")) %>%
    mutate(group = "Other Individuals")

    under_3_diets <- bind_rows(weaning, IV_281_weaning, IV_734_weaning)
    
    under_3_diets <- under_3_diets %>%
      mutate(age.bin = cut(MIDPOINT.serial.age, 
                           breaks = 4, 
                           labels = c("0–0.75", "0.75–1.5", "1.5–2.25", "2.25–3")))
    
    # accounting for the trophic enrichment effect: Ambrose and Norr, 1993; Koch, 1998; Krueger and Sullivan, 1984; Schoeninger, 1985
    ## -3.4% for d13C -- originally had these as added, but they really should be subtracted from human data 7.15.25
    ## -1.8% for d15N 
    ###7/22/25: going to use the diet-to-collagen enrichment from Fernandes et al 2015 instead:
    ## -4.8 +/- 0.5% for d13C
    ## -5.5 +/- 0.5% for d15N 
    ##I think I will *add* these to the dietary sources instead of subtracting from the humans
    
    seuss.effect <- 1.5 #value added to modern plant sources to correct for carbon enrichment in the atmosphere
    d13C.diet.to.collagen <- 4.8 #diet-to-collagen enrichment factor from Fernandes et al 2015 *added to dietary sources OR subtracted from human data
    d15N.diet.to.collagen <- 5.5 #diet-to-collagen enrichment factor from Fernandes et al 2015 *added to dietary sources OR subtracted from human data
    
    
    #these data were downloaded from IsoArcH on Feb 26, 2025, I want only the 'barley' data
    #reference information is available in the data
    barley_data <- BI_comparative_data %>% 
      filter(`Reported.biological.identification` == "barley")
    barley_data <- barley_data %>%
      mutate(d13C.corrected = case_when(
        Context == "archaeological" ~ `δ13C (‰vpdb)` + d13C.diet.to.collagen, #for archaeological samples, will only add diet-to-collagen
        Context == "modern" ~ `δ13C (‰vpdb)` + seuss.effect + d13C.diet.to.collagen, #for modern samples, will account for seuss effect and diet-to-collagen offset
        TRUE ~ `δ13C (‰vpdb)`)) #the barley data is from Nitsch et al 2017 and is archaeological
    barley_data <- barley_data %>%
      select(c(`Location.name`, 
               `Reported.biological.identification`, 
               `Context`, 
               `δ13C (‰vpdb)`, 
               `δ15N (‰air)`,
               `Short.references`,
               `d13C.corrected`))
    barley_data <- barley_data %>%
      mutate(type = `Reported.biological.identification`,
             identifier = `Location.name`,
             d13C = `δ13C (‰vpdb)`,
             d15N.corrected = `δ15N (‰air)` + d15N.diet.to.collagen, #diet to collagen
             d15N = `δ15N (‰air)`, 
             source = `Short.references`)
    barley_data <- barley_data %>%
      select(c(`type`, 
               `identifier`, 
               `Context`, 
               `d13C`,
               `d13C.corrected`,
               `d15N`,
               `d15N.corrected`, 
               `source`)) %>%
      mutate(type = case_when(
        type == "barley" ~ "Barley Mash",  # Check if 'type' is 'broomcorn millet'
        TRUE ~ type )) # Keep the existing value if the condition is not met
    

    #these data are from mulitple sources, you can find the reference in 'source'
    model_diets <- model_diets %>%
      mutate(d15N = as.numeric(d15N)) %>%
      filter(!type == "Millet Mash") #decided to remove millet here
    model_diets <- model_diets %>%
      mutate(d13C.corrected = case_when(
        Context == "archaeological" ~ d13C + d13C.diet.to.collagen, #for archaeological samples, will only add diet-to-collagen
        Context == "modern" ~ d13C + seuss.effect + d13C.diet.to.collagen, #for modern samples, will account for seuss effect and diet-to-collagen offset
        Context == "modeled" ~ d13C + d13C.diet.to.collagen, #the modeled foods here come from archaeological data and do not need to be corrected for the Seuss effect
        TRUE ~ d13C),
        d15N.corrected = d15N + d15N.diet.to.collagen) #diet to collagen
    
    
    model_diets <- model_diets %>%
      full_join(barley_data)

    
    convex_hulls.diet <- model_diets %>%
      group_by(type) %>%
      slice(chull(d13C.corrected, d15N.corrected))  # Get boundary points for convex hull
    
    dietary_contributions <- ggplot() + 
      geom_polygon(
        data = convex_hulls.diet,   # Convex hulls for modeled diets
        aes(x = d13C.corrected, 
            y = d15N.corrected, 
            fill = type, 
            group = type), 
            alpha = 0.6, color = NA, inherit.aes = FALSE
      ) +
      scale_fill_viridis_d(name = "Dietary Source", option = "rocket") + 
      new_scale_fill() + 
      geom_point(
        data = under_3_diets,   # Individual diet points
        aes(x = d13C, 
            y = d15N.Air,
            size = MIDPOINT.serial.age,
            fill = "grey50",
            shape = group),
            color = "grey20",
            stroke = 1,  
            alpha = 0.9) +
      scale_shape_manual(
        name = "Individual",
        values = c(
          "IV_281" = 21,
          "IV_734" = 24,
          "Other Individuals" = 22)) +
      scale_size_continuous(
        name = "Age",
        breaks = c(0.5, 1, 1.5, 2, 2.5, 3),
        label = c("0-0.5", "0.5-1", "1-1.5", "1.5-2", "2-2.5", "2.5-3")) +
      guides(
        fill = "none",
        color = guide_legend(title = "Individual"),
        size = guide_legend(title = "Age of Serial Sample")) +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 24)) +
      ggtitle("Potential Contributions to the Diet of Individuals Under 3 At Phaleron") +
      labs(x = expression(paste(delta^13, "C"[collagen], " (\u2030, VPDB) Dietary Component")), 
           y = expression(paste(delta^15, "N"[collagen], " (\u2030, AIR) Dietary Component"))) +
      theme(
        legend.position = "right",
        legend.text = element_text(size = 16, colour = "grey50"),
        legend.title = element_text(size = 18, colour = "grey50"),
        axis.text = element_text(size = 16, colour = "grey50"),
        axis.title = element_text(size = 18, colour = "grey50"),
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold", colour = "grey50"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16, colour = "grey50"),
        plot.title.position = "plot",
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linetype = "solid"),
        panel.grid.major.y = element_line(linetype = "solid"),
        panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
      )
    
    dietary_contributions
    
    
################################################################################
     
####GENDER######################################################################

     M1_data <- M1_data %>% #condensing sex estimation into male, female, or ambiguous only
       mutate(`cond.sex.est` = case_when(
         `sex.estimation`== "possible female" ~ "female",
         `sex.estimation` == "possible male" ~ "male",
         TRUE ~ `sex.estimation`
       ))

     infant_diet <- infant_diet %>%
       mutate(`cond.sex.est` = case_when(
         `sex.estimation`== "possible female" ~ "female",
         `sex.estimation` == "possible male" ~ "male",
         TRUE ~ `sex.estimation`
       ))
     
     child_diet <- child_diet %>%
       mutate(`cond.sex.est` = case_when(
         `sex.estimation`== "possible female" ~ "female",
         `sex.estimation` == "possible male" ~ "male",
         TRUE ~ `sex.estimation`
       ))
     
  ##KRUSKALL-WALLIS####
  gender.kruskal.N <- kruskal.test(d15N.Air ~ cond.sex.est, data = M1_data)
  gender.kruskal.C <- kruskal.test(d13C ~ cond.sex.est, data = M1_data)
       
  gender.infant.kruskal.N <- kruskal.test(d15N.Air ~ cond.sex.est, data = infant_diet)
  gender.infant.kruskal.C <- kruskal.test(d13C ~ cond.sex.est, data = infant_diet)
       
  gender.children.kruskal.N <- kruskal.test(d15N.Air ~ cond.sex.est, data = child_diet)
  gender.children.kruskal.C <- kruskal.test(d13C ~ cond.sex.est, data = child_diet)
       
      kruskal.table.sex <- tibble(
         test = c(
           "Overall d15N",
           "Overall d13C",
           "Infant d15N",
           "Infant d13C",
           "Child d15N",
           "Child d13C"),
         H_statistic = c(
           gender.kruskal.N$statistic,
           gender.kruskal.C$statistic,
           gender.infant.kruskal.N$statistic,
           gender.infant.kruskal.C$statistic,
           gender.children.kruskal.N$statistic,
           gender.children.kruskal.C$statistic),
         df = c(
           gender.kruskal.N$parameter,
           gender.kruskal.C$parameter,
           gender.infant.kruskal.N$parameter,
           gender.infant.kruskal.C$parameter,
           gender.children.kruskal.N$parameter,
           gender.children.kruskal.C$parameter),
         p_value = c(
           gender.kruskal.N$p.value,
           gender.kruskal.C$p.value,
           gender.infant.kruskal.N$p.value,
           gender.infant.kruskal.C$p.value,
           gender.children.kruskal.N$p.value,
           gender.children.kruskal.C$p.value))
      
       kruskal.table.sex
       #write.csv(kruskal.table.sex, "kruskal.table.sex.csv")
       
  ##PAIRWISE WILCOX####
       #Overall diet
       cond.sex.est.pairwise.wilcox.N <- M1_data %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ cond.sex.est,
           p.adjust.method = "holm")
       
       cond.sex.est.pairwise.wilcox.C <- M1_data %>%
         pairwise_wilcox_test(
           formula = d13C ~ cond.sex.est,
           p.adjust.method = "holm")
       
       #Infants
       infant.cond.sex.est.pairwise.wilcox.N <- infant_diet %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ cond.sex.est,
           p.adjust.method = "holm")
       
       infant.cond.sex.est.pairwise.wilcox.C <- infant_diet %>%
         pairwise_wilcox_test(
           formula = d13C ~ cond.sex.est,
           p.adjust.method = "holm")
       
       #Children
       child.cond.sex.est.pairwise.wilcox.N <- child_diet %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ cond.sex.est,
           p.adjust.method = "holm")
       
       child.cond.sex.est.pairwise.wilcox.C <- child_diet %>%
         pairwise_wilcox_test(
           formula = d13C ~ cond.sex.est,
           p.adjust.method = "holm")
       
       #Results
       pairwise.wilcox.table.sex <- bind_rows(
         cond.sex.est.pairwise.wilcox.N %>% mutate(test = "Overall d15N"),
         cond.sex.est.pairwise.wilcox.C %>% mutate(test = "Overall d13C"),
         infant.cond.sex.est.pairwise.wilcox.N %>% mutate(test = "Infant d15N"),
         infant.cond.sex.est.pairwise.wilcox.C %>% mutate(test = "Infant d13C"),
         child.cond.sex.est.pairwise.wilcox.N %>% mutate(test = "Child d15N"),
         child.cond.sex.est.pairwise.wilcox.C %>% mutate(test = "Child d13C")) %>%
         select(test, group1, group2, W = statistic, p.adj, p.adj.signif)
       
       pairwise.wilcox.table.sex
       #write.csv(pairwise.wilcox.table.sex, "pairwise.wilcox.table.sex.csv")
       
  ##FIGURE 4: GENDERED DIFFERENCES####     

       age_comparison <- age_comparison %>%
         mutate(`cond.sex.est` = case_when(
           `sex.estimation`== "possible female" ~ "female",
           `sex.estimation` == "possible male" ~ "male",
           TRUE ~ `sex.estimation`
         ))
       
       age_comparison_summary <- age_comparison %>%
         group_by(`burial.number`) %>%
         summarise(n = sum(!is.na(specimen.number))) %>%
         filter(!(n <= 2)) %>%
         ungroup()
       
       age_comparison_filtered <- age_comparison %>%
         filter(`burial.number` %in% age_comparison_summary$burial.number)
       
       age_comparison_filtered <- age_comparison_filtered %>%
         left_join(M1_data_summary_grouped)
       
       age_comparison_filtered <- age_comparison_filtered %>%
         mutate(sex_order = case_when(
           cond.sex.est == "male" ~ 1,
           cond.sex.est == "female" ~ 3,
           cond.sex.est == "ambiguous" ~ 2,
           TRUE ~ 0
         )) %>%
         mutate(burial.number = forcats::fct_reorder2(burial.number,  sex_order, d15N.Air_mean)) %>%
         mutate(burial.number = forcats::fct_reorder2(burial.number, d15N.Air_mean, sex_order)) ##idk how to not have to reorder it twice, this is just how it is
       
       

       diet.age.C <- ggplot(age_comparison_filtered, 
                            aes(x = d13C, y = burial.number, fill = cond.sex.est)) +
         geom_density_ridges(alpha = 0.3, color = "grey20") +
         geom_point(aes(shape = group, color = group),
                    position = position_jitter(height = 0.1, width = 0),
                    color = "grey20",
                    size = 3,
                    alpha = 0.8) +
         ggtitle("") +
         labs(
           x = expression(paste(delta^13, "C"[collagen], " (\u2030, VPDB)")),
           y = "Burial Number",
           fill = "Sex Estimation",
           shape = "Social Age",
           color = "Social Age"
         ) +
         scale_fill_viridis_d(option = "mako") +   # ridge fill
         scale_color_viridis_d(option = "mako") +  # points color
         scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
         geom_vline(
           xintercept = adult.Phaleron.C, 
           linetype = "dashed", 
           color = "#FFA600", 
           linewidth = 2
         ) +
         theme_minimal() +
         theme(
           legend.position = "right",
           legend.text = element_text(size = 16, colour = "grey50"),
           legend.title = element_text(size = 18, colour = "grey50"),
           legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
           axis.text = element_text(size = 16, colour = "grey50"),
           axis.title = element_text(size = 18, colour = "grey50"),
           plot.title = element_text(hjust = 0.5, size = 24, face = "bold", colour = "grey50"),
           axis.text.x = element_text(angle = 45, hjust = 1, size = 16, colour = "grey50"),
           plot.title.position = "plot",
           axis.ticks = element_blank(),
           axis.line = element_line(colour = "grey50"),
           panel.grid = element_line(color = "#b4aea9"),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_line(linetype = "solid"),
           panel.grid.major.y = element_line(linetype = "solid"),
           panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
           plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
         )
       diet.age.C
       
       diet.age.N <- ggplot(age_comparison_filtered, 
                            aes(x = d15N.Air, y = burial.number, fill = cond.sex.est)) +
         geom_density_ridges(alpha = 0.3, color = "grey20") +
         geom_point(aes(shape = group, color = group),
                    position = position_jitter(height = 0.1, width = 0),
                    color = "grey20",
                    size = 3,
                    alpha = 0.8) +
         ggtitle("") +
         labs(
           x = expression(paste(delta^15, "N"[collagen], " (\u2030, AIR)")),
           y = "Burial Number",
           fill = "Sex Estimation",
           shape = "Social Age",
           color = "Social Age"
         ) +
         scale_fill_viridis_d(option = "mako") +   # ridge fill
         scale_color_viridis_d(option = "mako") +  # points color
         scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
         geom_vline(
           xintercept = adult.Phaleron.N, 
           linetype = "dotted", 
           color = "#69b3a2", 
           linewidth = 3
         ) +
         theme_minimal() +
         theme(
           legend.position = "none",
           axis.text = element_text(size = 16, colour = "grey50"),
           axis.title = element_text(size = 18, colour = "grey50"),
           plot.title = element_text(hjust = 0.5, size = 24, face = "bold", colour = "grey50"),
           axis.text.x = element_text(angle = 45, hjust = 1, size = 16, colour = "grey50"),
           plot.title.position = "plot",
           axis.ticks = element_blank(),
           axis.line = element_line(colour = "grey50"),
           panel.grid = element_line(color = "#b4aea9"),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_line(linetype = "solid"),
           panel.grid.major.y = element_line(linetype = "solid"),
           panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
           plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
         )
       diet.age.N
       
       gender.figure <- diet.age.N + diet.age.C + 
         plot_layout(ncol = 2) + 
         plot_annotation(
           title = "Early Childhood Diet by Sex and Age",
           theme = theme(
             plot.title = element_text(hjust = 0.5, size = 24, face = "bold", colour = "grey50"), 
             plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
           )
         )
       
       gender.figure
       

  
################################################################################
       
####MORTUARY ANALYSIS###########################################################
  ##TORSO POSITION (SUPINE, PRONE, LEFT SIDE, RIGHT SIDE)####
    ##KRUSKALL-WALLIS####
       torso.position.kruskal.N <- kruskal.test(d15N.Air ~ torso.position, data = M1_data)
       torso.position.kruskal.C <- kruskal.test(d13C ~ torso.position, data = M1_data)
       
       torso.position.infant.kruskal.N <- kruskal.test(d15N.Air ~ torso.position, data = infant_diet)
       torso.position.infant.kruskal.C <- kruskal.test(d13C ~ torso.position, data = infant_diet)
       
       torso.position.children.kruskal.N <- kruskal.test(d15N.Air ~ torso.position, data = child_diet)
       torso.position.children.kruskal.C <- kruskal.test(d13C ~ torso.position, data = child_diet)
       
       kruskal.table.torso <- tibble(
         test = c(
           "Overall d15N",
           "Overall d13C",
           "Infant d15N",
           "Infant d13C",
           "Child d15N",
           "Child d13C"
         ),
         H_statistic = c(
           torso.position.kruskal.N$statistic,
           torso.position.kruskal.C$statistic,
           torso.position.infant.kruskal.N$statistic,
           torso.position.infant.kruskal.C$statistic,
           torso.position.children.kruskal.N$statistic,
           torso.position.children.kruskal.C$statistic
         ),
         df = c(
           torso.position.kruskal.N$parameter,
           torso.position.kruskal.C$parameter,
           torso.position.infant.kruskal.N$parameter,
           torso.position.infant.kruskal.C$parameter,
           torso.position.children.kruskal.N$parameter,
           torso.position.children.kruskal.C$parameter
         ),
         p_value = c(
           torso.position.kruskal.N$p.value,
           torso.position.kruskal.C$p.value,
           torso.position.infant.kruskal.N$p.value,
           torso.position.infant.kruskal.C$p.value,
           torso.position.children.kruskal.N$p.value,
           torso.position.children.kruskal.C$p.value
         )
       )
       
       kruskal.table.torso
       #write.csv(kruskal.table.torso, "kruskal.table.torso.csv")
       
    ##PAIRWISE WILCOX####
       #Overall diet
       torso.position.pairwise.wilcox.N <- M1_data %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ torso.position,
           p.adjust.method = "holm")
       
       torso.position.pairwise.wilcox.C <- M1_data %>%
         pairwise_wilcox_test(
           formula = d13C ~ torso.position,
           p.adjust.method = "holm")
       
       #Infants
       infant.torso.position.pairwise.wilcox.N <- infant_diet %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ torso.position,
           p.adjust.method = "holm")
       
       infant.torso.position.pairwise.wilcox.C <- infant_diet %>%
         pairwise_wilcox_test(
           formula = d13C ~ torso.position,
           p.adjust.method = "holm")
       
       #Children
       child.torso.position.pairwise.wilcox.N <- child_diet %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ torso.position,
           p.adjust.method = "holm")
       
       child.torso.position.pairwise.wilcox.C <- child_diet %>%
         pairwise_wilcox_test(
           formula = d13C ~ torso.position,
           p.adjust.method = "holm")
       
       #Results
       pairwise.wilcox.table.torso.position <- bind_rows(
         torso.position.pairwise.wilcox.N %>% mutate(test = "Overall d15N"),
         torso.position.pairwise.wilcox.C %>% mutate(test = "Overall d13C"),
         infant.torso.position.pairwise.wilcox.N %>% mutate(test = "Infant d15N"),
         infant.torso.position.pairwise.wilcox.C %>% mutate(test = "Infant d13C"),
         child.torso.position.pairwise.wilcox.N %>% mutate(test = "Child d15N"),
         child.torso.position.pairwise.wilcox.C %>% mutate(test = "Child d13C")) %>%
         select(test, group1, group2, W = statistic, p.adj, p.adj.signif)
       
       pairwise.wilcox.table.torso.position
       #write.csv(pairwise.wilcox.table.torso.position, "pairwise.wilcox.table.torso.position.csv")
       
       
  ##BODY POSITIONING (EXTENDED, SEMI-FLEXED, FLEXED, HYPERFLEXED)####
    ##condensing semi-flexed and hyper-flexed into flexed
       M1_data <- M1_data %>%
         mutate(body.pos.redux = case_when(
           body.position.2 %in% c("semi-flexed", "hyper-flexed") ~ "flexed",
           TRUE ~ body.position.2
         ))
       infant_diet <- infant_diet %>%
         mutate(body.pos.redux = case_when(
           body.position.2 %in% c("semi-flexed", "hyper-flexed") ~ "flexed",
           TRUE ~ body.position.2
         ))
       child_diet <- child_diet %>%
         mutate(body.pos.redux = case_when(
           body.position.2 %in% c("semi-flexed", "hyper-flexed") ~ "flexed",
           TRUE ~ body.position.2
         )) 
       
    ##KRUSKALL-WALLIS#### 
       body.pos.redux.kruskal.N <- kruskal.test(d15N.Air ~ body.pos.redux, data = M1_data)
       body.pos.redux.kruskal.C <- kruskal.test(d13C ~ body.pos.redux, data = M1_data)
       
       body.pos.redux.infant.kruskal.N <- kruskal.test(d15N.Air ~ body.pos.redux, data = infant_diet)
       body.pos.redux.infant.kruskal.C <- kruskal.test(d13C ~ body.pos.redux, data = infant_diet)
       
       body.pos.redux.children.kruskal.N <- kruskal.test(d15N.Air ~ body.pos.redux, data = child_diet)
       body.pos.redux.children.kruskal.C <- kruskal.test(d13C ~ body.pos.redux, data = child_diet)
       
       kruskal.table.body.pos.redux <- tibble(
         test = c(
           "Overall d15N",
           "Overall d13C",
           "Infant d15N",
           "Infant d13C",
           "Child d15N",
           "Child d13C"
         ),
         H_statistic = c(
           body.pos.redux.kruskal.N$statistic,
           body.pos.redux.kruskal.C$statistic,
           body.pos.redux.infant.kruskal.N$statistic,
           body.pos.redux.infant.kruskal.C$statistic,
           body.pos.redux.children.kruskal.N$statistic,
           body.pos.redux.children.kruskal.C$statistic
         ),
         df = c(
           body.pos.redux.kruskal.N$parameter,
           body.pos.redux.kruskal.C$parameter,
           body.pos.redux.infant.kruskal.N$parameter,
           body.pos.redux.infant.kruskal.C$parameter,
           body.pos.redux.children.kruskal.N$parameter,
           body.pos.redux.children.kruskal.C$parameter
         ),
         p_value = c(
           body.pos.redux.kruskal.N$p.value,
           body.pos.redux.kruskal.C$p.value,
           body.pos.redux.infant.kruskal.N$p.value,
           body.pos.redux.infant.kruskal.C$p.value,
           body.pos.redux.children.kruskal.N$p.value,
           body.pos.redux.children.kruskal.C$p.value
         )
       )
       
       kruskal.table.body.pos.redux
       #write.csv(kruskal.table.body.pos.redux, "kruskal.table.body.pos.redux.csv")
       
    ##PAIRWISE WILCOX####
       #Overall diet
       body.pos.redux.pairwise.wilcox.N <- M1_data %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ body.pos.redux,
           p.adjust.method = "holm")
       
       body.pos.redux.pairwise.wilcox.C <- M1_data %>%
         pairwise_wilcox_test(
           formula = d13C ~ body.pos.redux,
           p.adjust.method = "holm")
       
       #Infants
       infant.body.pos.redux.pairwise.wilcox.N <- infant_diet %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ body.pos.redux,
           p.adjust.method = "holm")
       
       infant.body.pos.redux.pairwise.wilcox.C <- infant_diet %>%
         pairwise_wilcox_test(
           formula = d13C ~ body.pos.redux,
           p.adjust.method = "holm")
       
       #Children
       child.body.pos.redux.pairwise.wilcox.N <- child_diet %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ body.pos.redux,
           p.adjust.method = "holm")
       
       child.body.pos.redux.pairwise.wilcox.C <- child_diet %>%
         pairwise_wilcox_test(
           formula = d13C ~ body.pos.redux,
           p.adjust.method = "holm")
       
       #Results
       pairwise.wilcox.table.body.pos.redux <- bind_rows(
         body.pos.redux.pairwise.wilcox.N %>% mutate(test = "Overall d15N"),
         body.pos.redux.pairwise.wilcox.C %>% mutate(test = "Overall d13C"),
         infant.body.pos.redux.pairwise.wilcox.N %>% mutate(test = "Infant d15N"),
         infant.body.pos.redux.pairwise.wilcox.C %>% mutate(test = "Infant d13C"),
         child.body.pos.redux.pairwise.wilcox.N %>% mutate(test = "Child d15N"),
         child.body.pos.redux.pairwise.wilcox.C %>% mutate(test = "Child d13C")) %>%
         select(test, group1, group2, W = statistic, p.adj, p.adj.signif)
       
       pairwise.wilcox.table.body.pos.redux
       #write.csv(pairwise.wilcox.table.body.pos.redux, "pairwise.wilcox.table.body.pos.redux.csv")
       
       
  ##RESTRAINED (LIKELY UNRESTRAINED, POSSIBLE RESTRAINT, METAL RESTRAINT, UNOBSERVABLE)####
       ##going to recode hands and feet so that it reflects "metal restraint" or "possible restraint"  
       ##if either hands or feet are bound, and "likely unrestrained" if either hands or feet are 
       ##observable.

       M1_data <- M1_data %>%
         mutate(restrained = case_when(
           restrained.hands == "metal restraint" | restrained.feet == "metal restraint" ~ "metal restraint",
           restrained.hands == "possible restraint" | restrained.feet == "possible restraint" ~ "possible restraint",
           restrained.hands == "likely unrestrained" | restrained.feet == "likely unrestrained" ~ "likely unrestrained",
           TRUE ~ "unobservable"
         ))

       infant_diet <- infant_diet %>%
         mutate(restrained = case_when(
           restrained.hands == "metal restraint" | restrained.feet == "metal restraint" ~ "metal restraint",
           restrained.hands == "possible restraint" | restrained.feet == "possible restraint" ~ "possible restraint",
           restrained.hands == "likely unrestrained" | restrained.feet == "likely unrestrained" ~ "likely unrestrained",
           TRUE ~ "unobservable"
         ))
       
       child_diet <- child_diet %>%
         mutate(restrained = case_when(
           restrained.hands == "metal restraint" | restrained.feet == "metal restraint" ~ "metal restraint",
           restrained.hands == "possible restraint" | restrained.feet == "possible restraint" ~ "possible restraint",
           restrained.hands == "likely unrestrained" | restrained.feet == "likely unrestrained" ~ "likely unrestrained",
           TRUE ~ "unobservable"
         ))
       
    ##KRUSKALL-WALLIS####      
       restrained.kruskal.N <- kruskal.test(d15N.Air ~ restrained, data = M1_data)
       restrained.kruskal.C <- kruskal.test(d13C ~ restrained, data = M1_data)
       
       restrained.infant.kruskal.N <- kruskal.test(d15N.Air ~ restrained, data = infant_diet)
       restrained.infant.kruskal.C <- kruskal.test(d13C ~ restrained, data = infant_diet)
       
       restrained.children.kruskal.N <- kruskal.test(d15N.Air ~ restrained, data = child_diet)
       restrained.children.kruskal.C <- kruskal.test(d13C ~ restrained, data = child_diet)
       
       kruskal.table.restrained <- tibble(
         test = c(
           "Overall d15N",
           "Overall d13C",
           "Infant d15N",
           "Infant d13C",
           "Child d15N",
           "Child d13C"
         ),
         H_statistic = c(
           restrained.kruskal.N$statistic,
           restrained.kruskal.C$statistic,
           restrained.infant.kruskal.N$statistic,
           restrained.infant.kruskal.C$statistic,
           restrained.children.kruskal.N$statistic,
           restrained.children.kruskal.C$statistic
         ),
         df = c(
           restrained.kruskal.N$parameter,
           restrained.kruskal.C$parameter,
           restrained.infant.kruskal.N$parameter,
           restrained.infant.kruskal.C$parameter,
           restrained.children.kruskal.N$parameter,
           restrained.children.kruskal.C$parameter
         ),
         p_value = c(
           restrained.kruskal.N$p.value,
           restrained.kruskal.C$p.value,
           restrained.infant.kruskal.N$p.value,
           restrained.infant.kruskal.C$p.value,
           restrained.children.kruskal.N$p.value,
           restrained.children.kruskal.C$p.value
         )
       )
       
       kruskal.table.restrained
       #write.csv(kruskal.table.restrained, "kruskal.table.restrained.csv")
       
    ##PAIRWISE WILCOX####
       #Overall diet
       restrained.pairwise.wilcox.N <- M1_data %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ restrained,
           p.adjust.method = "holm")
       
       restrained.pairwise.wilcox.C <- M1_data %>%
         pairwise_wilcox_test(
           formula = d13C ~ restrained,
           p.adjust.method = "holm")
       
       #Infants
       infant.restrained.pairwise.wilcox.N <- infant_diet %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ restrained,
           p.adjust.method = "holm")
       
       infant.restrained.pairwise.wilcox.C <- infant_diet %>%
         pairwise_wilcox_test(
           formula = d13C ~ restrained,
           p.adjust.method = "holm")
       
       #Children
       child.restrained.pairwise.wilcox.N <- child_diet %>%
         pairwise_wilcox_test(
           formula = d15N.Air ~ restrained,
           p.adjust.method = "holm")
       
       child.restrained.pairwise.wilcox.C <- child_diet %>%
         pairwise_wilcox_test(
           formula = d13C ~ restrained,
           p.adjust.method = "holm")
       
       #Results
       pairwise.wilcox.table.restrained <- bind_rows(
         restrained.pairwise.wilcox.N %>% mutate(test = "Overall d15N"),
         restrained.pairwise.wilcox.C %>% mutate(test = "Overall d13C"),
         infant.restrained.pairwise.wilcox.N %>% mutate(test = "Infant d15N"),
         infant.restrained.pairwise.wilcox.C %>% mutate(test = "Infant d13C"),
         child.restrained.pairwise.wilcox.N %>% mutate(test = "Child d15N"),
         child.restrained.pairwise.wilcox.C %>% mutate(test = "Child d13C")) %>%
         select(test, group1, group2, W = statistic, p.adj, p.adj.signif)
       
       pairwise.wilcox.table.restrained
       #write.csv(pairwise.wilcox.table.restrained, "pairwise.wilcox.table.restrained.csv")
       
     
       
  ##MULTIPLE BURIAL####
   ##don't need Kruskall-Wallis since there are only two levels (yes or no)
       ##WILCOX####
       #Overall diet
       multiple.burial.wilcox.N <- M1_data %>%
         wilcox_test(
           formula = d15N.Air ~ multiple.burial,
           p.adjust.method = "holm")
       
       multiple.burial.wilcox.C <- M1_data %>%
         wilcox_test(
           formula = d13C ~ multiple.burial,
           p.adjust.method = "holm")
       
       #Infants
       infant.multiple.burial.wilcox.N <- infant_diet %>%
         wilcox_test(
           formula = d15N.Air ~ multiple.burial,
           p.adjust.method = "holm")
       
       infant.multiple.burial.wilcox.C <- infant_diet %>%
         wilcox_test(
           formula = d13C ~ multiple.burial,
           p.adjust.method = "holm")
       
       #Children
       child.multiple.burial.wilcox.N <- child_diet %>%
         wilcox_test(
           formula = d15N.Air ~ multiple.burial,
           p.adjust.method = "holm")
       
       child.multiple.burial.wilcox.C <- child_diet %>%
         wilcox_test(
           formula = d13C ~ multiple.burial,
           p.adjust.method = "holm")
       
       #Results
       wilcox.table.multiple.burial <- bind_rows(
         multiple.burial.wilcox.N %>% mutate(test = "Overall d15N"),
         multiple.burial.wilcox.C %>% mutate(test = "Overall d13C"),
         infant.multiple.burial.wilcox.N %>% mutate(test = "Infant d15N"),
         infant.multiple.burial.wilcox.C %>% mutate(test = "Infant d13C"),
         child.multiple.burial.wilcox.N %>% mutate(test = "Child d15N"),
         child.multiple.burial.wilcox.C %>% mutate(test = "Child d13C")
       ) %>%
         adjust_pvalue(method = "holm") %>%  # creates column p.adj
         add_significance() %>%
         select(test, group1, group2, W = statistic, p.adj, p.adj.signif)
       
       wilcox.table.multiple.burial
       #write.csv(wilcox.table.multiple.burial, "wilcox.table.multiple.burial.csv")
       
  ##FIGURE 6####
       ##need to re-sort the data so it plots nicely
       burial_type_filtered <- age_comparison_filtered %>%
         mutate(burial.type.sup = factor(burial.type.sup,
                                         levels = c("extended supine", 
                                                    "flexed supine", 
                                                    "flexed side", 
                                                    "CMB3", 
                                                    "single shackled", 
                                                    "prone")),

                burial.num = as.numeric(str_extract(burial.number, "\\d+$")),
                burial.type.num = as.numeric(burial.type.sup),
                burial.sort.string = sprintf("%02d_%04d_%s", 
                                             burial.type.num, 
                                             burial.num, 
                                             burial.number)) %>%
         arrange(burial.sort.string) %>%
         mutate(burial.y = factor(burial.number, levels = unique(burial.number)))
       
  diet.grave.type.C <- ggplot(burial_type_filtered,
                              aes(x = d13C, y = burial.y, fill = burial.type.sup)) +
         geom_density_ridges(alpha = 0.3, 
                             color = "grey20",
                             scale = 0.9) +
         geom_point(aes(shape = group, color = group),
                    position = position_jitter(height = 0.1, width = 0),
                    color = "grey20",
                    size = 3,
                    alpha = 0.8)+
         labs(
           x = expression(delta^13*C[collagen] * " (\u2030, VPDB)"),
           y = "Burial Number",
           fill = "Burial Type",
           shape = "Social Age",
           color = "Social Age"
         ) +
         scale_fill_viridis_d(option = "mako") +
         scale_color_viridis_d(option = "mako") +
         scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
         geom_vline(xintercept = adult.Phaleron.C, 
                    linetype = "dashed", 
                    color = "#FFA600", 
                    linewidth = 2) +
         theme_minimal() +
         theme(
           legend.position = "right",
           legend.text = element_text(size = 16, colour = "grey50"),
           legend.title = element_text(size = 18, colour = "grey50"),
           legend.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
           axis.text = element_text(size = 16, colour = "grey50"),
           axis.title = element_text(size = 18, colour = "grey50"),
           plot.title = element_text(hjust = 0.5, size = 24, face = "bold", colour = "grey50"),
           axis.text.x = element_text(angle = 45, hjust = 1, size = 16, colour = "grey50"),
           plot.title.position = "plot",
           axis.ticks = element_blank(),
           axis.line = element_line(colour = "grey50"),
           panel.grid = element_line(color = "#b4aea9"),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_line(linetype = "solid"),
           panel.grid.major.y = element_line(linetype = "solid"),
           panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
           plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
         )
       
       diet.grave.type.C
       
  diet.grave.type.N <- ggplot(burial_type_filtered, 
                            aes(x = d15N.Air, y = burial.y, fill = burial.type.sup)) +
         geom_density_ridges(alpha = 0.3, 
                             color = "grey20",
                             scale = 0.9) +
         geom_point(aes(shape = group, color = group),
                    position = position_jitter(height = 0.1, width = 0),
                    color = "grey20",
                    size = 3,
                    alpha = 0.8) +
         ggtitle("") +
         labs(
           x = expression(paste(delta^15, "N"[collagen], " (\u2030, AIR)")),
           y = "Burial Number",
           fill = "Burial Type",
           shape = "Social Age",
           color = "Social Age"
         ) +
         scale_fill_viridis_d(option = "mako") +   # ridge fill
         scale_color_viridis_d(option = "mako") +  # points color
         scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
         geom_vline(
           xintercept = adult.Phaleron.N, 
           linetype = "dotted", 
           color = "#69b3a2", 
           linewidth = 3
         ) +
         theme_minimal() +
         theme(
           legend.position = "none",
           axis.text = element_text(size = 16, colour = "grey50"),
           axis.title = element_text(size = 18, colour = "grey50"),
           plot.title = element_text(hjust = 0.5, size = 24, face = "bold", colour = "grey50"),
           axis.text.x = element_text(angle = 45, hjust = 1, size = 16, colour = "grey50"),
           plot.title.position = "plot",
           axis.ticks = element_blank(),
           axis.line = element_line(colour = "grey50"),
           panel.grid = element_line(color = "#b4aea9"),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_line(linetype = "solid"),
           panel.grid.major.y = element_line(linetype = "solid"),
           panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
           plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
         )
       diet.grave.type.N
       
       burial.figure <- diet.grave.type.N + diet.grave.type.C + 
         plot_layout(ncol = 2, guides = "collect") + 
         plot_annotation(
           title = "Early Childhood Diet by Burial Type",
           theme = theme(
             plot.title = element_text(hjust = 0.5, size = 24, face = "bold", colour = "grey50"), 
             plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
           )
         )
       
       burial.figure
       
       
       
################################################################################



