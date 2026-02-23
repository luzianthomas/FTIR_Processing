# FT-IR Data Processing

## Overview
This is a small project to demonstrate basic data processing and analysis of time-domain
interferograms recorded on an FT-IR spectrometer. The script extracts the
zero-path-difference (ZPD) regions, performs Fourier transformation, computes
transmittance and absorbance spectra, plots the spectra in the typically relevant region (400 to 4000 cm^-1)
and automatically detects absorption peaks.

I intend this code to be a compact example of scientific data handling and signal
processing in Python.

## Input Data
The script expects comma-separated text files containing interferogram data
with two columns of which only the first column is used:

1. acquisition index
2. Signal intensity
