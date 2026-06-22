library(tidyverse) # Load the tidyverse package
Infection_chain  <-  read_csv("MMF - Final-3.csv") # read in the "data"
Drs_Visit  <-  read_csv("MMF - DoctorVisits.csv") # read in the "data"

dim(Infection_chain)			# Determine the number of rows and columns
Infection_chain		        # Look at the beginning of the dataset

#fixing date format
Infection_chain$date <- as.Date(Infection_chain$Date, format="%d/%m/%Y")

dim(Drs_Visit)			# Determine the number of rows and columns
Drs_Visit		        # Look at the beginning of the dataset

Drs_Visit$date <- as.Date(Drs_Visit$date, format="%d/%m/%Y")
