#Deliverable D5.4: Spatial Stochastic Modelling of African Swine Fever Spread and Control in Denmark

This repository contains the modelling framework and analyses developed under WiLiMan Deliverable D5.4, led by the University of Copenhagen (UCPH). The project evaluates the potential spread of African Swine Fever (ASF) among Danish pig farms and assesses the effectiveness of current EU-recommended control measures.

#Overview

A spatially explicit stochastic farm-level transmission model was developed using the Gillespie simulation framework to simulate ASF spread across Danish pig farms connected through a spatial contact network. Farms transition through five epidemiological states: Susceptible (S), Exposed (E), Infectious (I), Detected (D), and Culled (C).

#The model incorporates:

Farm-to-farm transmission through contact networks
Residual transmission from detected farms
Farm-level environmental contamination
Background risk of external virus introduction
Detection, movement restrictions, zoning, surveillance, culling, and biosecurity interventions
Scenarios Evaluated
EFSA-like control measures: Immediate implementation of movement restrictions, protection and surveillance zones, enhanced surveillance, culling, and biosecurity following detection.
No control measures: Unrestricted disease spread after introduction.

Multiple stochastic simulations were performed to quantify outbreak uncertainty and epidemic trajectories.

#Key Findings
ASF can spread rapidly through the Danish pig farm network in the absence of intervention.
EFSA-like control measures substantially reduce outbreak size, transmission, and geographic spread.
The mean farm-level reproduction number decreased from 20.68 (R₀) to 5.15 (Rₑff), representing an estimated 75% reduction in transmission potential.
Spatial spread remained largely localized under control measures but expanded extensively without intervention.
Sensitivity analyses identified post-detection transmission reduction as one of the most important determinants of outbreak severity.
Highly connected farms contribute disproportionately to transmission, highlighting the value of risk-based surveillance and targeted interventions.
Implications for Preparedness

The results demonstrate that rapid detection, movement restrictions, zoning, surveillance, culling, and strengthened biosecurity can significantly reduce ASF spread. The framework supports contingency planning and preparedness efforts in Denmark and other ASF-free regions.

#Author
Ofosuhene Okofrobour Apenteng

#Deliverable Output

The project delivered a comprehensive spatial stochastic modelling framework for ASF preparedness, outbreak simulation, reproduction number estimation, intervention assessment, and sensitivity analysis. The work forms the basis of a scientific manuscript describing the methodology and findings.

