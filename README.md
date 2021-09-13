# TFM: Hourly Prediction of Air Temperature in Spain with Deep Learning-based Strategies

This repository includes the main code required to replicate the experiments performed as part of the Thesis 'Hourly Prediction of Air Temperature in Spain with Deep Learning-based Strategies' for the Master in Big Data Analytics (Universidad Carlos III de Madrid).


## Description

This project’s experimental part focuses on forecasting the hourly mean air temperature of different regions throughout Spain. First, diverse datasets of climate features from observational stations in Spain were analyzed to obtain reliable data to perform the analysis. Second, preprocessing and cleaning methodologies were carried out to prepare the data for deep learning-based forecasting. Then, two novel networks were proposed for this task: a *Graph Convolutional Long Short Term Memory Network* (GCN-LSTM) that considersthe spatial location of the sites, and a *Multifactor Spatio-Temporal Correlation model based on aHybrid Convolutional Long Short Term Memory Network* (MFSTC-CNN-LSTM) that accounts for therelationship between air temperature and other climate parameters. Lastly, the performance of the different models was evaluated and compared with other proposals in the literature.

The code is divided in two folders. The folder `Missing values imputation` includes an R code that loads NOAA data and imputes missing values using a user-defined function that computes the *Inverse Distance Weighting* (IDW) method, followed by either *Autoregressive Integrated Moving Average with Kalman Smoothing* (ARIMA) or *Moving Average* (MA) imputation. A statistical analysis is performed to select the optimal method: deterministic statistics and kernel density functions of the imputed data are compared with those of the raw data. The folder `Models` includes five Python codes for each of the studied models (LSTM, Stacked-LSTM, CNN-LSTM, GCN-LSTM, MFSTC-CNN-LSTM) that preprocesses the imputed data, trains the model and evaluates its performance. They hyperparameters to be tunned (`seq_len`, `units`) and modified (`pre_len`) are indictated in the notebooks.

## Getting Started

### Dependencies

* R (4.0.3)
* R libraries: `dplyr` (1.0.4), `tidyverse` (1.3.0), `plyr` (1.8.6), `naniar` (0.6.0), `xts` (0.12.1), `fpp2` (2.4), `cowplot` (1.1.1), `reshape2` (1.4.4), `data.table` (1.13.6), `imputeTS` (3.2), `worldmet` (0.9.2), `moments` (0.14), `ggplot2` (3.3.3)
* Python (3.8.5)
* Python libraries: `numpy` (1.19.2), `pandas` (1.1.3), `matplotlib` (3.3.2), `tensorflow` (2.5.0), `stellargraph` (1.2.1), `sklearn` (0.23.2)
* Anaconda (2020.11)

### Installing

* Notice that there are files to paths in your computer that must be changed to your desired destination.


## Authors

Lucía García-Duarte Sáenz

## Version History

* 0.1
    * Initial Release
