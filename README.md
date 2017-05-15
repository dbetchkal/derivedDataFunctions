# derivedDataFunctions
Functions to calculate acoustical metrics from standard NPS derived data files

<!DOCTYPE html>
<html>
<head>
</head>

<body>
<h2> Functions contained in this module:</h2>


<h3> AMPLITUDE METRICS FROM SRCID</h3>
<p> quantile_Lmax(srcid, q, weight = "A", source = "all") </p>
<p> mad_Lmax(srcid, weight = "A", source = "all")</p>
<p> iqr_Lmax(srcid, weight = "A", source = "all")</p>
<p> mean_Lmax(srcid, weight = "A", source = "all")</p>
<p> stdev_Lmax(srcid, weight = "A", source = "all")</p>
<p> stderr_Lmax(srcid, weight = "A", source = "all")</p>
<p> quantile_SEL(srcid, q, weight = "A", source = "all") </p>


<hr>
<h3> DURATION METRICS FROM SRCID</h3>
<p> total_event_duration(srcid, source = "all")</p>
<p> quantile_event_duration(srcid, q, source = "all")</p>
<p> mad_event_duration(srcid, source = "all")</p>
<p> iqr_event_duration(srcid, source = "all")</p>
<p> mean_event_duration(srcid, source = "all")</p>
<p> stdev_event_duration(srcid, source = "all")</p>
<p> stderr_event_duration(srcid, source = "all")</p>
<p> total_audible_dur_hourly(dailypa, hour, source = "all")</p>
<p> mean_audible_duration_hourly(dailypa, hour, source = "all")</p>


<hr>
<h3> COUNT METRICS FROM SRCID</h3>
<p> total_count(srcid, source = "all")</p>
<p> percentageOfAll_bySource(srcid, id_code)</p>
<p> percentageOfAir_bySource(srcid, id_code)</p>
<p> propJetRatio(srcid)</p>
<p> DENABCMP_SPL_exceedance(srcid, zone, source = "all")</p>
<p> DENABCMP_SPL_exceedanceRate(srcid, zone, source = "all")</p>


<hr>
<h3> DATASET DESCRIPTION METRICS FROM SRCID</h3>
<p> number_of_days_splatted(srcid)</p>
<p> days_splatted(srcid)</p>
<p> SPLAT_center_date(srcid)</p>

<hr>
<h3> NOISE FREE INTERVAL</h3>
<p> mean_NFI(srcid, source = "all", unit="hours")</p>
<p> quantile_NFI(srcid, q, source = "all", unit="hours")</p>
<p> NFI_list(srcid, source = "all", unit="hours") </p>


<hr>
<h3> PERCENT TIME AUDIBLE METRICS FROM DAILYPA</h3>
<p> quantile_hourlyPA(dailypa, q, hour, source = "all")</p>
<p> quantile_dailyPA(dailypa, q, source = "all")</p>
<p> DENABCMP_PA_exceedance_all(dailypa, zone, start_hour = 0, end_hour = 23, source = "all")</p>
<p> overall_PA(dailypa, source = "all")</p>
<p> </p>
<p> </p>
<p> </p>
<p> </p>

<hr>
<h3> EVENT RATE, COUNT, AND SATURATION METRICS FROM DAILYPA</h3>
<p> quantile_eventsPerDay(dailypa, q, source = "all")</p>
<p> total_events(dailypa, source = "all")</p>
<p> event_saturation(dailypa, start_hour = 0, end_hour = 23, source = "all")</p>

<hr>
<h3> EVENT RATES ABOVE AMBIENT FROM LOUDEVENTS </h3>
<p> DENABCMP_events_exceedance(loudevents, zone)</p>
<p> quantile_eventRate_overAmbient(loudevents, q)</p>
<p> mean_eventRate_overAmbient(loudevents)</p>
<p> stdev_eventRate_overAmbient(loudevents)</p>
<p> stderr_eventRate_overAmbient(loudevents)</p>

<hr>
<h3> STANDARD ACOUSTIC EXCEEDANCE METRICS </h3>
<p> L90(metrics, season="Summer", weight = "A")</p>
<p> Lnat(metrics, season="Summer", weight = "A")</p>
<p> L50(metrics, season="Summer", weight = "A")</p>
<p> L10(metrics, season="Summer", weight = "A")</p>

</body>
</html>
