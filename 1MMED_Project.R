############################################################
# MMF Outbreak Investigation
# Group members: Rabi and Nurudeen
# MMED 2026
############################################################

# Load required packages
# tidyverse is used for data cleaning, analysis, and plotting
# lubridate is used for working with dates and times
library(tidyverse)
library(lubridate)

############################################################
# 1. Import the datasets
############################################################

# Read the main MMF infection dataset
mmf <- read_csv("MMF - Final-3.csv")

# Read the doctor visit dataset, which contains symptomatic cases
doctor <- read_csv("MMF - DoctorVisits.csv")

############################################################
# 2. Explore the datasets
############################################################

# Check the number of rows and columns
dim(mmf)
mmf
dim(doctor)
doctor

# Check column names
names(mmf)
names(doctor)

# View the first few rows
head(mmf)
head(doctor)

############################################################
# 3. Rename columns for easier coding
############################################################

# Rename columns so we avoid using spaces in variable names
mmf <- mmf |>
  rename(
    person_infected = `Person Infected`,
    infected_by = `Infected by`
  )

############################################################
# 4. Clean text variables
############################################################

# Remove extra spaces from names
mmf <- mmf |>
  mutate(
    person_infected = str_trim(person_infected),
    infected_by = str_trim(infected_by)
  )

doctor <- doctor |>
  mutate(
    Name = str_trim(Name)
  )

############################################################
# 5. Correct inconsistent name spellings
############################################################

# First check all unique names in both datasets
sort(unique(mmf$person_infected))
sort(unique(mmf$infected_by))
sort(unique(doctor$Name))

# Correct known spelling differences
mmf <- mmf |>
  mutate(
    person_infected = case_when(
      person_infected == "Thuli" ~ "Thulisile",
      TRUE ~ person_infected
    ),
    infected_by = case_when(
      infected_by == "Thuli" ~ "Thulisile",
      TRUE ~ infected_by
    )
  )

doctor <- doctor |>
  mutate(
    Name = case_when(
      Name == "Upendo" ~ "Pendo",
      TRUE ~ Name
    )
  )

# Check names again
sort(unique(mmf$person_infected))
sort(unique(mmf$infected_by))
sort(unique(doctor$Name))


############################################################
# 6. Create date-time variables
############################################################

# 1. Correct incomplete date first
mmf <- mmf |>
  mutate(
    Date = case_when(
      Date == "18/06" ~ "18/06/2026",
      TRUE ~ Date
    )
  )

# 2. Convert Date from character to date format
mmf <- mmf |>
  mutate(
    Date = dmy(Date)
  )

# 3. Create infection datetime
mmf <- mmf |>
  mutate(
    infection_datetime = as.POSIXct(
      paste(Date, Time),
      format = "%Y-%m-%d %H:%M:%S"
    )
  )

# 4. Convert doctor Date from character to date format
doctor <- doctor |>
  mutate(
    Date = dmy(Date)
  )

# 5. Create doctor visit datetime
doctor <- doctor |>
  mutate(
    doctor_datetime = as.POSIXct(
      paste(Date, Time),
      format = "%Y-%m-%d %H:%M:%S"
    )
  )
###confirm change to date format
head(mmf)
head(doctor)

############################################################
# 7. Arrange infections by time
############################################################

# Sort the dataset from first infection to last infection
# Then create infection_order to show the sequence of infections
mmf <- mmf |>
  arrange(infection_datetime) |>
  mutate(infection_order = row_number())

# View the infection order
mmf |>
  select(infection_order, infection_datetime, person_infected, infected_by)

############################################################
# 8. Plot epidemic curve
############################################################

# This shows how infections occurred over time
# binwidth = 3600 means infections are grouped by 1 hour
ggplot(mmf |> filter(!is.na(infection_datetime)), 
       aes(x = infection_datetime)) +
  geom_histogram(binwidth = 3600) +
  theme_bw() +
  labs(
    title = "Epidemic Curve of MMF Outbreak",
    x = "Time of infection",
    y = "Number of infections"
  )

####better graph..Daily Epidemic Curve

ggplot(mmf, aes(x = Date)) +
  geom_bar() +
  theme_bw() +
  labs(
    title = "Daily Epidemic Curve of MMF Outbreak",
    x = "Date",
    y = "Number of infections"
  )

####Cumulative epidemic curve
cum_cases <- mmf |>
  arrange(infection_datetime) |>
  mutate(cumulative_cases = row_number())

ggplot(cum_cases,
       aes(x = infection_datetime,
           y = cumulative_cases)) +
  geom_step() +
  theme_bw() +
  labs(
    title = "Cumulative MMF Cases",
    x = "Date",
    y = "Cumulative infections"
  )

###plot by Infection Timeline
ggplot(
  mmf |> arrange(infection_datetime),
  aes(
    x = infection_datetime,
    y = reorder(person_infected, infection_datetime),
    colour = infected_by
  )
) +
  geom_point(size = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Timeline of MMF Infections",
    x = "Time of infection",
    y = "Participant",
    colour = "Infected by"
  )

############################################################
# 9. Estimate secondary infections
############################################################

# Count how many people each infected person infected
secondary_cases <- mmf |>
  count(infected_by, name = "secondary_infections") |>
  arrange(desc(secondary_infections))

# View secondary infections
secondary_cases

############################################################
# 8. Estimate effective reproduction number
############################################################

# Re is estimated as the average number of secondary infections
# caused by each infected person
Re_estimate <- mean(secondary_cases$secondary_infections)

# View Re estimate
Re_estimate

############################################################
# 9. Merge infection data with doctor visit data
############################################################

# Merge doctor visit data with the main MMF data
# Participants found in doctor dataset are classified as symptomatic
mmf_doctor <- mmf |>
  left_join(
    doctor |> select(Name, doctor_datetime, symptomatic),
    by = c("person_infected" = "Name")
  ) |>
  mutate(
    symptomatic = replace_na(symptomatic, 0)
  )

mmf_doctor |>
  select(person_infected, infected_by, infection_datetime, symptomatic, doctor_datetime) |>
  print(n = 33)

############################################################
# 10. Calculate symptomatic proportion
############################################################

# Summarise total infected, symptomatic cases, and proportion symptomatic
mmf_doctor |>
  summarise(
    total_infected = n(),
    symptomatic_cases = sum(symptomatic, na.rm = TRUE),
    symptomatic_proportion = mean(symptomatic, na.rm = TRUE)
  )

############################################################
# Deterministic SIR base case for MMF
############################################################

library(deSolve)
library(tidyr)

# Define SIR model
sir <- function(t, y, parms){
  with(c(as.list(y), parms), {
    dSdt <- -beta * S * I / N
    dIdt <- beta * S * I / N - gamma * I
    dRdt <- gamma * I
    
    return(list(c(dSdt, dIdt, dRdt)))
  })
}

# Set initial conditions
N <- 33        # total participants, change if your true class size differs
I0 <- 1        # initial infected
R0_init <- 0   # no one immune at start
S0 <- N - I0 - R0_init

initial_state <- c(S = S0, I = I0, R = R0_init)

# Set starting parameter values
values <- c(
  beta = 4,     # transmission rate;
  gamma = 1       # recovery rate; 1 means average infectious period = 1 day
)

# Time in days
time.out <- seq(0, 5, by = 0.04)

# Run deterministic SIR model
sir_output <- data.frame(
  lsoda(
    y = initial_state,
    times = time.out,
    func = sir,
    parms = values
  )
)

# Plot deterministic SIR output(I only)
ggplot(sir_output, aes(x = time, y = I)) +
  geom_line() +
  theme_bw() +
  labs(
    title = "Deterministic SIR Model for MMF",
    x = "Time since first infection (days)",
    y = "Number infectious"
  )


#######PLOT  SIR model
sir_long <- sir_output |>
  pivot_longer(
    cols = c(S, I, R),
    names_to = "Compartment",
    values_to = "Count"
  )

ggplot(sir_long,
       aes(x = time,
           y = Count,
           colour = Compartment)) +
  geom_line(linewidth = 1.2) +
  theme_bw() +
  labs(
    title = "Deterministic SIR Model for MMF",
    x = "Time since first infection (days)",
    y = "Number of participants"
  )

############################################################
# Stochastic SIR using Gillespie algorithm
############################################################

event_sir <- function(time, S, I, R, params) {
  with(as.list(params), {
    
    # Event rates
    trans.rate <- beta * S * I / N
    recov.rate <- gamma * I
    tot.rate <- trans.rate + recov.rate
    
    # Stop if no event can occur
    if (tot.rate == 0) {
      return(tibble(time = time, S = S, I = I, R = R))
    }
    
    # Time to next event
    event.time <- time + rexp(1, tot.rate)
    
    # Choose event type
    dd <- runif(1)
    
    if (dd < trans.rate / tot.rate) {
      # Infection event: S -> I
      S <- S - 1
      I <- I + 1
    } else {
      # Recovery event: I -> R
      I <- I - 1
      R <- R + 1
    }
    
    tibble(time = event.time, S = S, I = I, R = R)
  })
}

run_sir_gillespie <- function(run_id, params, MAXTIME) {
  
  ts <- list(
    tibble(time = 0, S = params$N - 1, I = 1, R = 0)
  )
  
  current <- ts[[1]]
  i <- 1
  
  while (current$time < MAXTIME && current$I > 0) {
    current <- event_sir(
      current$time,
      current$S,
      current$I,
      current$R,
      params
    )
    
    i <- i + 1
    ts[[i]] <- current
  }
  
  bind_rows(ts) |>
    mutate(run = run_id)
}

# Parameters
params <- list(
  N = 33,
  beta = 4,
  gamma = 1
)

MAXTIME <- 5
NUMRUNS <- 20

# Run several stochastic simulations
gillespie_runs <- map_dfr(
  1:NUMRUNS,
  run_sir_gillespie,
  params = params,
  MAXTIME = MAXTIME
)

# Plot stochastic simulations
ggplot(gillespie_runs, aes(x = time, y = I, group = run)) +
  geom_step(alpha = 0.4) +
  theme_bw() +
  labs(
    title = "Stochastic SIR Gillespie Simulations for MMF",
    x = "Time since first infection (days)",
    y = "Number infectious"
  )



