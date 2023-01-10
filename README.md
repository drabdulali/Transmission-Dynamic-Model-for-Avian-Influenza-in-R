# Transmission-Dynamic-Model-for-Avian-Influenza-in-R
Avian influenza, also known as bird flu, is a viral infection that primarily affects domestic poultry, such as chickens, ducks, and turkeys, as well as wild birds. 
OPTIMAL CONTROL FOR INFECTED IMMIGRANTS IN AN AVIAN INFLUENZA TRANSMISSION DYNAMICS – A codec Transmission Dynamic Model in R for analysis with definite parameters

Introduction:

Avian influenza, also known as bird flu, is a viral infection that primarily affects domestic poultry, such as chickens, ducks, and turkeys, as well as wild birds. The virus is transmitted through the saliva, nasal secretions, and feces of infected birds, and it can be spread through direct contact with infected birds or through contact with contaminated feed, water, or equipment. In severe cases, avian influenza can lead to high mortality rates in infected birds, causing significant economic losses for farmers and the poultry industry.
There are different subtypes of avian influenza virus, with varying degrees of pathogenicity. The H5N1 subtype is considered to be the most virulent, as it can cause severe illness and death in birds and has also been known to cause infections in humans. Other subtypes, such as H7N9 and H9N2, can also cause disease in poultry, but they are not known to be as deadly as H5N1.
One of the key challenges in controlling avian influenza outbreaks is the potential for infected birds to migrate to new locations, spreading the virus to new areas and potentially leading to new outbreaks. In order to better understand and control the transmission dynamics of avian influenza in these situations, mathematical modeling can be used to study the spread of the virus and to identify potential intervention strategies.
One approach that can be used to model avian influenza transmission dynamics is called "optimal control for infected immigrants." In this model, the total number of birds in the location of interest is represented by the variable N, and the average birth rate in birds is represented by the variable β. The proportion of infected birds among the migrated birds is represented by the variable ρ, and the total number of migrated birds is represented by the variable M. The total number of susceptible birds is represented by the variable S, and the total number of infected birds is represented by the variable I. The natural death rate in birds is represented by the variable μ, and the infection transmission rate from bird to bird is represented by the variable α. The flu-induced death rate for birds is represented by the variable ν.
The model represents the transmission dynamics of avian influenza as a set of differential equations. These equations describe the rate at which susceptible birds become infected, the rate at which infected birds recover or die, and the rate at which new birds are born or migrate into the population. The model also takes into account the potential for interventions, such as vaccination or culling of infected birds, to control the spread of the virus.
Using this model, researchers can estimate the values of the different parameters based on data from avian influenza outbreaks. For example, it is estimated that in the location of interest there are total of 20000 birds, average birth rate in birds is 0.05R.Proportion of infected in migrated birds is 0.01M and the total number of migrated birds is 2000. The total number of susceptible birds is 40000, total number of infected birds is 3. The natural death rate in birds is 0.0255, infection transmission rate from bird to bird is 0.8, and flu-induced death rate for birds is 0.653. These estimates can be used to make predictions about the potential impact of different intervention strategies on the spread of the virus and to identify the most effective ways to control avian influenza outbreaks.
However, It's important to note that the model must be calibrated with real data and validated with real-world outbreaks and experiments in order to have valid results. Also, the parameter values can vary depending on the location and the specific strain of avian influenza virus.

Analysis in R a Transmission Dynamic Model for Avian Influenza with following parameters

Estimated value N Total number of birds in the location of interest 20000
ß Average birth rate in birds 0.05R 
Proportion of infected in migrated birds 0.01M 
Total number of migrated birds 2000 S Total number of susceptible birds 40000 
I Total number of infected birds 3µ 
Natural death rate in birds 0.0255 α 
Infection transmission rate from bird to bird 0.8ν 
Flu-induced death rate for birds 0.653

Analysis in R

library(deSolve)
library(ggplot2)

# Define the model
avian_influenza_model <- function(t, y, p) {
    with(as.list(c(p, y)), {
        dS <- (beta*N - mu - alpha*I*S - S*R*u)*S
        dI <- (alpha*I*S + S*R*u - (mu + nu)*I)*I
        dM <- (beta*N - mu - alpha*I*S - S*R*u)*M
        dU <- -k*u
        return(list(c(dS, dI, dM, dU)))
    })
}

# Define the parameters
p <- c(beta = 0.05, mu = 0.0255, alpha = 0.8, nu = 0.653, R = 0.01, N = 20000, k = 0.01)

# Define the initial conditions
y0 <- c(S = 40000, I = 3, M = 2000, u = 0.05)

# Define the time points for the simulation
times <- seq(from = 0, to = 365, by = 1)

# Run the simulation
out <- as.data.frame(ode(y = y0, times = times, func = avian_influenza_model, parms = p))

# Set custom theme
custom_theme <- theme(
    axis.title.x = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.background = element_rect(fill = "gray"),
    plot.title = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "bottom"
)

# Plot the results
ggplot(data = out, aes(x = time)) +
    geom_line(aes(y = S, color = "Susceptible")) +
    geom_line(aes(y = I, color = "Infected")) +
    geom_line(aes(y = M, color = "Migrated")) +
    geom_line(aes(y = u, color = "Control variable")) +
    ggtitle("Avian influenza transmission dynamic model with optimal control") +
    ylab("Number of individuals") +
    xlab("Time (days)") +
    scale_color_manual(name = 'Compartment', values = c("Susceptible" = "forestgreen", "Infected" = "red", "Migrated" = "skyblue", "Control variable" = "purple")) +
    scale_y_continuous(limits = c(0, max(out$S, out$I, out$M, out$u))) +
    custom_theme

 
