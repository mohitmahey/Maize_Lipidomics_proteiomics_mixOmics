# 🌱 Maize Lipidomics & Proteomics – Multi-Omics Integration with mixOmics

This repository contains scripts and workflows for correlating **lipidomics** and **proteomics** datasets derived from **plastoglobule (PG)** and **thylakoid (Thy)** fractions of maize under **drought** and **heat stress** conditions.  

The analyses are performed using the [mixOmics](http://mixomics.org/) R package, focusing on **Regularized Canonical Correlation Analysis (rCCA)**, **correlation networks**, and **PQN normalization** to integrate lipid–protein interactions across stresses.  

---

## 📂 Repository Structure

- **`normalised_pqn.R`**  
  Script for **Probabilistic Quotient Normalization (PQN)** of proteomics and lipidomics data.  

- **`correlation_among_matrix_drought.R`**  
  Calculates correlation matrices among proteomics/lipidomics datasets under **drought** stress.  

- **`correlation_among_matrix_heat.R`**  
  Calculates correlation matrices among proteomics/lipidomics datasets under **heat** stress.  

- **`correlation_drought_networks.R`**  
  Builds **network graphs** from drought correlation matrices for visualization (e.g., in Cytoscape).  

- **`rCCA_drought_PG.R`**  
  rCCA analysis integrating **plastoglobule proteome** and **lipidome** under **drought**.  

- **`rCCA_drought_thy.R`**  
  rCCA analysis integrating **thylakoid proteome** and **lipidome** under **drought**.  

- **`rCCA_heat_PG.R`**  
  rCCA analysis integrating **plastoglobule proteome** and **lipidome** under **heat** stress.  

- **`README.md`**  
  Documentation for the project.  

- **`LICENSE`**  
  Licensed under **Apache 2.0**.  

---

## 🚀 Methods Overview

1. **Data Preprocessing**
   - PQN normalization applied to proteomics & lipidomics datasets.  
   - Filtering steps to remove low-information features.  

2. **Correlation Analysis**
   - Pairwise **Pearson correlations** among proteins and lipids.  
   - Separate analyses for **drought** and **heat** conditions.  

3. **Canonical Correlation (rCCA)**
   - Multi-block integration using mixOmics.  
   - Separate models for plastoglobule vs thylakoid datasets.  

4. **Network Construction**
   - Edges derived from correlation matrices (above threshold).  
   - Exportable to **Cytoscape** for visualization.  

---

## 📖 Usage

Clone the repo:

```bash
git clone https://github.com/mohitmahey/Maize_Lipidomics_proteiomics_mixOmics.git
cd Maize_Lipidomics_proteiomics_mixOmics
