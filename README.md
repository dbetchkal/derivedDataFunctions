# derivedDataFunctions
Functions to calculate acoustical metrics from standard NPS derived data files

## Functions contained in this module:


#### AMPLITUDE METRICS FROM SRCID
```python
quantile_amplitude(srcid, q, metric = "Lmax", weight = "A", source = "all")

mad_amplitude(srcid, metric = "Lmax", weight = "A", source = "all")

iqr_amplitude(srcid, metric = "Lmax", weight = "A", source = "all")

mean_amplitude(srcid, metric = "Lmax", weight = "A", source = "all")

stdev_amplitude(srcid, metric = "Lmax", weight = "A", source = "all")

stderr_amplitude(srcid, metric = "Lmax", weight = "A", source = "all")
```
______

#### DURATION METRICS FROM SRCID
```python
total_event_duration(srcid, source = "all")

quantile_event_duration(srcid, q, source = "all")

mad_event_duration(srcid, source = "all")

iqr_event_duration(srcid, source = "all")

mean_event_duration(srcid, source = "all")

stdev_event_duration(srcid, source = "all")

stderr_event_duration(srcid, source = "all")

total_audible_dur_hourly(dailypa, hour, source = "all")

mean_audible_duration_hourly(dailypa, hour, source = "all")
```
______

#### COUNT METRICS FROM SRCID
```python
total_count(srcid, source = "all")

percentageOfAll_bySource(srcid, id_code)

percentageOfAir_bySource(srcid, id_code)

propJetRatio(srcid)

DENABCMP_SPL_exceedance(srcid, zone, source = "all")

DENABCMP_SPL_exceedanceRate(srcid, zone, source = "all")
```
______

#### DATASET DESCRIPTION METRICS FROM SRCID
```python
number_of_days_splatted(srcid)

days_splatted(srcid)

SPLAT_center_date(srcid)
```
______

#### NOISE FREE INTERVAL
```python
mean_NFI(srcid, source = "all", unit="hours")

quantile_NFI(srcid, q, source = "all", unit="hours")

NFI_list(srcid, source = "all", unit="hours")
```
______

#### PERCENT TIME AUDIBLE METRICS FROM DAILYPA
```python
quantile_hourlyPA(dailypa, q, hour, source = "all")

quantile_dailyPA(dailypa, q, source = "all")

DENABCMP_PA_exceedance_all(dailypa, zone, start_hour = 0, end_hour = 23, source = "all")

overall_PA(dailypa, source = "all")
```
______

#### EVENT RATE, COUNT, AND SATURATION METRICS FROM DAILYPA
```python
quantile_eventsPerDay(dailypa, q, source = "all")

total_events(dailypa, source = "all")

event_saturation(dailypa, start_hour = 0, end_hour = 23, source = "all")
```
______

#### EVENT RATES ABOVE AMBIENT FROM LOUDEVENTS
```python
DENABCMP_events_exceedance(loudevents, zone)

quantile_eventRate_overAmbient(loudevents, q)

mean_eventRate_overAmbient(loudevents)

stdev_eventRate_overAmbient(loudevents)

stderr_eventRate_overAmbient(loudevents)
```
______

#### STANDARD ACOUSTIC EXCEEDANCE METRICS
```python
L90(metrics, season="Summer", weight = "A")

Lnat(metrics, season="Summer", weight = "A")

L50(metrics, season="Summer", weight = "A")

L10(metrics, season="Summer", weight = "A")
```
