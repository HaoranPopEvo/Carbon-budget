# Carbon-budget

This repository provides modelling code for tree carbon budget at a species level. 

Code files: <br>
**FE_weather.R**: Code for 'weather' module. It simulates dirnual variation in air temperature, vapour pressure decifit, and net shortwave radiation.<br>
**FE_photo.R**: Code for 'photosynthesis' module. It simulates the net CO2 assimilation rate.<br>
**FE_transpir.R**: Code for 'transpiration' module. It simulates the transpiration demand.<br>
**FE_allocation.R**: Code for 'respiration' and 'allocation' module. It simulates maintainence respiration for each plant tissue, and the carbon allocation process.<br>
**Benchmark_NPP.R**: Code for benchmark the carbon budget model with real measurements on European ash tree (*Fraxinus excelsior* L.).<br><br>

Data files: <br>
**allocation-data.csv**: measured carbon allocation data.<br>
**Amance-climate-ECMWF.csv**: daily climatic time series.<br>
**Majewski-2024-data.csv**: measured relationship between net CO2 assimilation and stomatal conductance.<br>
**benchmark_result.csv**: result for benchmarking the model with tree carbon budget data.<br>

