import pandas as pd
import numpy as np

#------------------------------------------------------------------------------------------------------------------
# ### AMPLITUDE METRICS FROM SRCID


def quantile_amplitude(srcid, q, metric = "Lmax", weight = "A", source = "all"):
    """
    Calculate a quantile for amplitude values of SPLAT-annotated sources at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    q: float, a value that indicates the quantile desired, from 0.0 (minimum) to 1.0 (maximum.) 
    metric: str, optional.  The amplitude metric to use when preforming the calculation, either "Lmax" or "SEL". Defaults to "Lmax" if unspecified.
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    formatted float
    """

    # allow the user to enter weighting networks either way, but convert to upper case
    w = weight.upper()

    # intialize a weighting lookup function based on the metric used
    if(metric == "Lmax"):
        lookup = {"A":"MaxSPL", "T":"MaxSPLt"}
    elif(metric == "SEL"):
        lookup = {"A":"SEL", "T":"SELt"}
    else:
        raise ValueError('metric must be either "Lmax" or "SEL"')

    
    if(type(source) == str):
        if(source.lower() == "all"):
            return float("{0:.1f}".format(srcid.loc[:,lookup[w]].quantile(q)))
            
        elif(source.lower() == "air"):
            return float("{0:.1f}".format(srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),lookup[w]].quantile(q)))
    
    else:
        return float("{0:.1f}".format(srcid.loc[srcid.srcID.isin(source), lookup[w]].quantile(q)))




def mad_amplitude(srcid, metric = "Lmax", weight = "A", source = "all"):
    """
    Calculate the median absolute deviation for amplitude values of SPLAT-annotated sources at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    metric: str, optional.  The amplitude metric to use when preforming the calculation, either "Lmax" or "SEL". Defaults to "Lmax" if unspecified.
    weight: str, specifies the acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.

    Returns
    -------
    formatted float
    """

    # allow the user to enter weighting networks either way, but convert to upper case
    w = weight.upper()

    # intialize a weighting lookup function based on the metric used
    if(metric == "Lmax"):
        lookup = {"A":"MaxSPL", "T":"MaxSPLt"}
    elif(metric == "SEL"):
        lookup = {"A":"SEL", "T":"SELt"}
    else:
        raise ValueError('metric must be either "Lmax" or "SEL"')
    

    if(type(source) == str):
        if(source.lower() == "all"):
            return float("{0:.1f}".format(pd.Series(abs(srcid.loc[:,lookup[w]].median() - srcid.loc[:,lookup[w]])).median()))
            
        elif(source.lower() == "air"):
            return float("{0:.1f}".format(pd.Series(abs(srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),lookup[w]].median() - srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),lookup[w]])).median()))
    
    else:
        return float("{0:.1f}".format(pd.Series(abs(srcid.loc[srcid.srcID.isin(source),lookup[w]].median() - srcid.loc[srcid.srcID.isin(source),lookup[w]])).median()))




def iqr_amplitude(srcid, metric = "Lmax", weight = "A", source = "all"):
    """
    Calculate the interquartile range for amplitude values of SPLAT-annotated sources at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    metric: str, optional.  The amplitude metric to use when preforming the calculation, either "Lmax" or "SEL". Defaults to "Lmax" if unspecified.
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    formatted float
    """

    # allow the user to enter weighting networks either way, but convert to upper case
    w = weight.upper()

    # intialize a weighting lookup function based on the metric used
    if(metric == "Lmax"):
        lookup = {"A":"MaxSPL", "T":"MaxSPLt"}
    elif(metric == "SEL"):
        lookup = {"A":"SEL", "T":"SELt"}
    else:
        raise ValueError('metric must be either "Lmax" or "SEL"')
    

    if(type(source) == str):
        if(source.lower() == "all"):
            return float("{0:.1f}".format(srcid.loc[:,lookup[w]].quantile(0.75) - srcid.loc[:,lookup[w]].quantile(0.25)))
            
        elif(source.lower() == "air"):
            return float("{0:.1f}".format(srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),lookup[w]].quantile(0.75) - srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),lookup[w]].quantile(0.25)))
    
    else:
        return float("{0:.1f}".format(srcid.loc[srcid.srcID.isin(source), lookup[w]].quantile(0.75) - srcid.loc[srcid.srcID.isin(source), lookup[w]].quantile(0.25)))        




def mean_amplitude(srcid, metric = "Lmax", weight = "A", source = "all"):
    """
    Calculate a mean amplitude value for SPLAT-annotated sources at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    metric: str, optional.  The amplitude metric to use when preforming the calculation, either "Lmax" or "SEL". Defaults to "Lmax" if unspecified.
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    formatted float
    """
    # allow the user to enter weighting networks either way, but convert to upper case
    w = weight.upper()

    # intialize a weighting lookup function based on the metric used
    if(metric == "Lmax"):
        lookup = {"A":"MaxSPL", "T":"MaxSPLt"}
    elif(metric == "SEL"):
        lookup = {"A":"SEL", "T":"SELt"}
    else:
        raise ValueError('metric must be either "Lmax" or "SEL"')


    if(type(source) == str):
        if(source.lower() == "all"):
            return float("{0:.1f}".format(srcid.loc[:,lookup[w]].mean()))
            
        elif(source.lower() == "air"):
            return float("{0:.1f}".format(srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),lookup[w]].mean()))
    
    else:
        return float("{0:.1f}".format(srcid.loc[srcid.srcID.isin(source), lookup[w]].mean()))




def stdev_amplitude(srcid, metric = "Lmax", weight = "A", source = "all"):
    """
    Calculate the standard deviation for amplitude values of SPLAT-annotated sources at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    metric: str, optional.  The amplitude metric to use when preforming the calculation, either "Lmax" or "SEL". Defaults to "Lmax" if unspecified.
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    formatted float
    """
    # allow the user to enter weighting networks either way, but convert to upper case
    w = weight.upper()

    # intialize a weighting lookup function based on the metric used
    if(metric == "Lmax"):
        lookup = {"A":"MaxSPL", "T":"MaxSPLt"}
    elif(metric == "SEL"):
        lookup = {"A":"SEL", "T":"SELt"}
    else:
        raise ValueError('metric must be either "Lmax" or "SEL"')

    
    if(type(source) == str):
        if(source.lower() == "all"):
            return float("{0:.1f}".format(srcid.loc[:,lookup[w]].std()))
            
        elif(source.lower() == "air"):
            return float("{0:.1f}".format(srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),lookup[w]].std()))
    
    else:
        return float("{0:.1f}".format(srcid.loc[srcid.srcID.isin(source), lookup[w]].std()))




def stderr_amplitude(srcid, metric = "Lmax", weight = "A", source = "all"):
    """
    Calculate the standard error of the mean for amplitude values of SPLAT-annotated sources at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    metric: str, optional.  The amplitude metric to use when preforming the calculation, either "Lmax" or "SEL". Defaults to "Lmax" if unspecified.
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    formatted float
    """
    # allow the user to enter weighting networks either way, but convert to upper case
    w = weight.upper()

    # intialize a weighting lookup function based on the metric used
    if(metric == "Lmax"):
        lookup = {"A":"MaxSPL", "T":"MaxSPLt"}
    elif(metric == "SEL"):
        lookup = {"A":"SEL", "T":"SELt"}
    else:
        raise ValueError('metric must be either "Lmax" or "SEL"')
    

    if(type(source) == str):
        if(source.lower() == "all"):
            return float("{0:.1f}".format(srcid.loc[:,lookup[w]].std()/np.sqrt(srcid.loc[:,lookup[w]].count())))
            
        elif(source.lower() == "air"):
            return float("{0:.1f}".format(srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),lookup[w]].std()/np.sqrt(srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),lookup[w]].count())))
    
    else:
        return float("{0:.1f}".format(srcid.loc[srcid.srcID.isin(source), lookup[w]].std()/np.sqrt(srcid.loc[srcid.srcID.isin(source), lookup[w]].count())))  


#------------------------------------------------------------------------------------------------------------------
# ### DURATION METRICS FROM SRCID

  
def total_event_duration(srcid, source = "all"):
    """
    Sum the duration of SPLAT-annotated sources at a site.  Note that this is NOT EQUAL to percent time audible.
    However, the following is always true:  total_event_duration >= total_audible_duration.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    formatted string (from timedelta)
    """
    
    if(type(source) == str):
        if(source.lower() == "all"):
            return srcid.loc[:,"len"].sum()
            
        elif(source.lower() == "air"):
            return srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),"len"].sum()
    
    else:
        return srcid.loc[srcid.srcID.isin(source),"len"].sum()



def quantile_event_duration(srcid, q, source = "all"):  
    """
    Return a quantile of durations from SPLAT-annotated sources at a site.  
    Note that this is NOT EQUAL to percent time audible.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    q: float, a value that indicates the quantile desired, from 0.0 (minimum) to 1.0 (maximum.)   
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    timedelta
    """
    
    if(type(source) == str):
        if(source.lower() == "all"):
            import datetime
            return datetime.timedelta(seconds = srcid.loc[:,"len"].quantile(q).total_seconds())
            
        elif(source.lower() == "air"):
            import datetime
            return datetime.timedelta(seconds = srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),"len"].quantile(q).total_seconds())
    
    else:
        import datetime
        return datetime.timedelta(seconds = srcid.loc[srcid.srcID.isin(source),"len"].quantile(q).total_seconds())




def mad_event_duration(srcid, source = "all"):  
    """
    Calculate the median absolute deviation for duration values of SPLAT-annotated sources at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library. 
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    timedelta
    """
    
    if(type(source) == str):
        if(source.lower() == "all"):
            import datetime
            mad = pd.Series(abs(srcid.loc[:,"len"].median() - srcid.loc[:,"len"])).median()
            return datetime.timedelta(seconds = mad.total_seconds())
            
        elif(source.lower() == "air"):
            import datetime
            mad = pd.Series(abs(srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),"len"].median() - srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),"len"])).median()
            return datetime.timedelta(seconds = mad.total_seconds())
    
    else:
        import datetime
        mad = pd.Series(abs(srcid.loc[srcid.srcID.isin(source),"len"].median() - srcid.loc[srcid.srcID.isin(source),"len"])).median()
        return datetime.timedelta(seconds = mad.total_seconds())




def iqr_event_duration(srcid, source = "all"):  
    """
    Calculate the interquartile range for durations of SPLAT-annotated sources at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library. 
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    timedelta
    """
    
    if(type(source) == str):
        if(source.lower() == "all"):
            import datetime
            iqr = srcid.loc[:,"len"].quantile(0.75) - srcid.loc[:,"len"].quantile(0.25)
            return datetime.timedelta(seconds = iqr.total_seconds())
            
        elif(source.lower() == "air"):
            import datetime
            iqr = srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),"len"].quantile(0.75) - srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),"len"].quantile(0.25)
            return datetime.timedelta(seconds = iqr.total_seconds())
    
    else:
        import datetime
        iqr = srcid.loc[srcid.srcID.isin(source),"len"].quantile(0.75) - srcid.loc[srcid.srcID.isin(source),"len"].quantile(0.25)
        return datetime.timedelta(seconds = iqr.total_seconds())


        
  
def mean_event_duration(srcid, source = "all"):
    """
    The mean duration of SPLAT-annotated sources at a site.  

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    timedelta
    """
    
    if(type(source) == str):
        if(source.lower() == "all"):
            import datetime
            return datetime.timedelta(seconds = srcid.loc[:,"len"].mean().total_seconds())
            
        elif(source.lower() == "air"):
            import datetime
            return datetime.timedelta(seconds = srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),"len"].mean().total_seconds())
    
    else:
        import datetime
        return datetime.timedelta(seconds = srcid.loc[srcid.srcID.isin(source),"len"].mean().total_seconds())


  
def stdev_event_duration(srcid, source = "all"):
    """
    The standard deviation of durations for SPLAT-annotated sources at a site.  

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    timedelta
    """
    
    if(type(source) == str):
        if(source.lower() == "all"):
            import datetime
            return datetime.timedelta(seconds = srcid.loc[:,"len"].std().total_seconds())
            
        elif(source.lower() == "air"):
            import datetime
            return datetime.timedelta(seconds = srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),"len"].std().total_seconds())
    
    else:
        import datetime
        return datetime.timedelta(seconds = srcid.loc[srcid.srcID.isin(source),"len"].std().total_seconds())


  
def stderr_event_duration(srcid, source = "all"):
    """
    The standard error of the mean for durations of SPLAT-annotated sources at a site.  

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    timedelta
    """
    
    if(type(source) == str):
        if(source.lower() == "all"):
            import datetime
            return datetime.timedelta(seconds = srcid.loc[:,"len"].std().total_seconds()/np.sqrt(srcid.loc[:, "len"].count()))
            
        elif(source.lower() == "air"):
            import datetime
            return datetime.timedelta(seconds = srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),"len"].std().total_seconds()/np.sqrt(srcid.loc[:, "len"].count()))
    
    else:
        import datetime
        return datetime.timedelta(seconds = srcid.loc[srcid.srcID.isin(source),"len"].std().total_seconds()/np.sqrt(srcid.loc[srcid.srcID.isin(source), "len"].count()))



def total_audible_dur_hourly(dailypa, hour, source = "all"): 
    """
    The total percent time audible for SPLAT-annotated sources at a site.  

    Parameters
    ----------
    dailypa: pandas dataframe representing NPS NSNSD dailypa file, formatted by soundDB library.
    hour: int, hour of the day to be summarized from 0 to 23.  To summarize the entire day iterate this function with "for hour in range(0, 23)"
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    timedelta
    """
    h = str(hour).zfill(2) + "h"
    
    if(type(source) == str):
        if(source.lower() == "all"):
            import datetime
            secs = (dailypa.loc[(slice(None), "Total_All"), "00h":"23h"]*3600)/100
            return datetime.timedelta(seconds = secs.sum()[h])
            
        elif(source.lower() == "air"):
            import datetime
            secs = (dailypa.loc[(slice(None), "Total_1"), "00h":"23h"]*3600)/100
            return datetime.timedelta(seconds = secs.sum()[h])
    
    else:
        import datetime
        secs = (dailypa.loc[(slice(None), str(source)), "00h":"23h"]*3600)/100
        return datetime.timedelta(seconds = secs.sum()[h])



def mean_audible_duration_hourly(dailypa, hour, source = "all"): 
    """
    The average percent time audible for SPLAT-annotated sources at a site.  

    Parameters
    ----------
    dailypa: pandas dataframe representing NPS NSNSD dailypa file, formatted by soundDB library.
    hour: int, hour of the day to be summarized from 0 to 23.  To summarize the entire day iterate this function with "for hour in range(0, 23)"
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    timedelta
    """
    h = str(hour).zfill(2) + "h"
    
    if(type(source) == str):
        if(source.lower() == "all"):
            import datetime
            tot_count = dailypa.loc[(slice(None), "Total_All"), "nEvents_24Hr"].sum()
            secs = (dailypa.loc[(slice(None), "Total_All"), "00h":"23h"]*3600)/100
            return datetime.timedelta(seconds = secs.sum()[h]/tot_count)
            
        elif(source.lower() == "air"):
            import datetime
            tot_count = dailypa.loc[(slice(None), "Total_1"), "nEvents_24Hr"].sum()
            secs = (dailypa.loc[(slice(None), "Total_1"), "00h":"23h"]*3600)/100
            return datetime.timedelta(seconds = secs.sum()[h]/tot_count)
    
    else:
        import datetime
        tot_count = dailypa.loc[(slice(None), str(source)), "nEvents_24Hr"].sum()
        secs = (dailypa.loc[(slice(None), str(source)), "00h":"23h"]*3600)/100
        return datetime.timedelta(seconds = secs.sum()[h]/tot_count)



#------------------------------------------------------------------------------------------------------------------
# ### COUNT METRICS FROM SRCID

  
def total_count(srcid, source = "all"):  # will give total props, total jets, total helicopters
    """
    The total number of noise events by source type.  

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", or "air".  Defaults to "all" if unspecified.
    
    Returns
    -------
    int
    """
    if(type(source) == str):
        if(source.lower() == "all"):
            return srcid.loc[:,"MaxSPL"].count()
        
        elif(source.lower() == "air"):
            return srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.),"MaxSPL"].count()

    else:       
        return srcid.loc[srcid.srcID.isin(source),"MaxSPL"].count()



  
def percentageOfAll_bySource(srcid, id_code):  
    """
    Counts by srcID code expressed as a percentage of all annotated sources at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    id_code: float, an srcID code used for annotating noise sources in SPLAT.
    
    Returns
    -------
    float, as a percentage
    """
    return 100*(srcid.loc[srcid.srcID == id_code,"MaxSPL"].count()/srcid.loc[:, "MaxSPL"].count())



  
def percentageOfAir_bySource(srcid, id_code):  
    """
    Counts by srcID code expressed as a percentage of aviation sources annotated at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    id_code: float, an srcID code used for annotating noise sources in SPLAT.
    
    Returns
    -------
    float, as a percentage
    """
    return 100*(srcid.loc[srcid.srcID == id_code,"MaxSPL"].count()/srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.), "MaxSPL"].count())



  
def propJetRatio(srcid):  
    """
    Returns the count of props at a site by the count of jets at a site.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    
    Returns
    -------
    float
    """
    return srcid.loc[srcid.srcID == 1.2,"MaxSPL"].count()/srcid.loc[srcid.srcID == 1.1,"MaxSPL"].count()



  
def DENABCMP_SPL_exceedance(srcid, zone, source = "all"):  # report the percentage of events exceeding DENA BCMP SPL standard
    """
    Reports the percentage of noise events exceeding the Denali Backcountry Managment Plan SPL standard.  

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    zone: str, the Denali backcountry managment plan zone:  "low", "medium", "high", or "very high".  case insensitive.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", or "air".  Defaults to "all" if unspecified, which is the truest expression of the plan as well.
    
    Returns
    -------
    float, as a percentage
    """

    key = zone.lower()
    lookup = {"low":40.0, "med":40.0, "medium":40.0, "high":60., "very high":60., "v. high":60., "veryhigh":60., "v high":60.}
    

    if(source.lower() == "all"):
        return 100*(srcid.loc[srcid.MaxSPL > lookup[key], "MaxSPL"].count()/srcid.loc[:,"MaxSPL"].count())

    elif(source.lower() == "air"):
        subset = srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.), "MaxSPL"]
        return 100*(subset.loc[subset > lookup[key]].count()/subset.count())



  
def DENABCMP_SPL_exceedanceRate(srcid, zone, source = "all"):  # the number of events exceeding DENA BCMP SPL standard per day
    """
    Reports the number of noise events per day exceeding the Denali Backcountry Managment Plan SPL standard.  

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    zone: str, the Denali backcountry managment plan zone:  "low", "medium", "high", or "very high".  case insensitive.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", or "air".  Defaults to "all" if unspecified, which is the truest expression of the plan as well.
    
    Returns
    -------
    float, as a percentage
    """

    key = zone.lower()
    lookup = {"low":40.0, "med":40.0, "medium":40.0, "high":60., "very high":60., "v. high":60., "veryhigh":60., "v high":60.}
    days = len(pd.Series(srcid.index.values).dt.date.unique())

    if(source.lower() == "all"):
        return srcid.loc[srcid.MaxSPL > lookup[key], "MaxSPL"].count()/days

    elif(source.lower() == "air"):
        subset = srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.), "MaxSPL"]
        return subset.loc[subset > lookup[key]].count()/days



#------------------------------------------------------------------------------------------------------------------
# ### DATASET DESCRIPTION METRICS FROM SRCID

  
def number_of_days_splatted(srcid):  
    """
    Returns the number of days of SPLAT analysis conducted for a site.  

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    
    Returns
    -------
    int
    """
    import datetime as dt
    return len(pd.Series(srcid.index.values).dt.date.unique())



  
def days_splatted(srcid): 
    """
    Returns a list of the days that were analyzed in SPLAT.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    
    Returns
    -------
    list of datetimes in chronological order
    """
    dates = pd.DataFrame(pd.Series(srcid.index).dt.date.unique())
    dates.columns = ["date"]
    return dates.sort(['date'])  



  
def SPLAT_center_date(srcid): 
    """
    Returns the center date of the analyzed period (useful for typical day lengths, etc...)

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    
    Returns
    -------
    datetime
    """

    dList = pd.Series(srcid.index).dt.date.unique()
    # remember, .strftime("%Y-%m-%d") will turn this result into a string if you need it in that format
    # for MM/DD format, .strftime("%Y-%m-%d")[5:12].replace("-", "/")  will do the trick
    return dList[len(dList)/2]



#------------------------------------------------------------------------------------------------------------------
# ### NOISE FREE INTERVAL


def mean_NFI(srcid, source = "all", unit="hours"): 
    """
    Returns the average NFI for selected source type.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    unit: str, a value that indicates the units desired for the output value.  Defaults to "hours".
    
    Returns
    -------
    numpy.float64
    """

    unitDict = {"seconds":1, "minutes":60, "hours":3600, "days":86400}

    if(type(source) == str):
        if(source.lower() == "all"):  
 
            dt = srcid.sort_index()

            NFIlst = pd.Series([(dt.index[i+1] - (dt.index[i] + dt['len'][i])).total_seconds()/unitDict[unit] for i in range(len(dt.index)+1) 
                if(i < len(dt.index)-1)])
    
            return NFIlst.mean()
        
        elif(source.lower() == "air"):
            dt = srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.), :]
            dt.sort_index()

            NFIlst = pd.Series([(dt.index[i+1] - (dt.index[i] + dt['len'][i])).total_seconds()/unitDict[unit] for i in range(len(dt.index)+1) 
                if(i < len(dt.index)-1)])
    
            return NFIlst.mean()
    
    else: 

        dt = srcid.loc[srcid.srcID.isin(source), :]
        dt.sort_index()

        NFIlst = pd.Series([(dt.index[i+1] - (dt.index[i] + dt['len'][i])).total_seconds()/unitDict[unit] for i in range(len(dt.index)+1) 
            if(i < len(dt.index)-1)])
    
        return NFIlst.mean()




def quantile_NFI(srcid, q, source = "all", unit="hours"): 
    """
    Returns the average NFI for selected source type.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    q: float, a value that indicates the quantile desired, from 0.0 (minimum) to 1.0 (maximum.)   
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    unit: str, a value that indicates the units desired for the output value.  Defaults to "hours".

    Returns
    -------
    numpy.float64
    """
    unitDict = {"seconds":1, "minutes":60, "hours":3600, "days":86400}

    if(type(source) == str):
        if(source.lower() == "all"):  
 
            dt = srcid.sort_index()

            NFIlst = pd.Series([(dt.index[i+1] - (dt.index[i] + dt['len'][i])).total_seconds()/unitDict[unit] for i in range(len(dt.index)+1) 
                if(i < len(dt.index)-1)])
    
            return NFIlst.quantile(q)
        
        elif(source.lower() == "air"):
            dt = srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.), :]
            dt.sort_index()

            NFIlst = pd.Series([(dt.index[i+1] - (dt.index[i] + dt['len'][i])).total_seconds()/unitDict[unit] for i in range(len(dt.index)+1) 
                if(i < len(dt.index)-1)])
    
            return NFIlst.quantile(q)
    
    else: 

        dt = srcid.loc[srcid.srcID.isin(source), :]
        dt.sort_index()

        NFIlst = pd.Series([(dt.index[i+1] - (dt.index[i] + dt['len'][i])).total_seconds()/unitDict[unit] for i in range(len(dt.index)+1) 
            if(i < len(dt.index)-1)])
    
        return NFIlst.quantile(q)
      

def NFI_list(srcid, source = "all", unit="hours"): 
    """
    Returns a list of all Noise Free Intervals for selected source type.

    Parameters
    ----------
    srcid: pandas dataframe representing NPS NSNSD srcid file, formatted by soundDB library.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    unit: str, a value that indicates the units desired for the output value.  Defaults to "hours".

    Returns
    -------
    pandas Series
    """
    unitDict = {"seconds":1, "minutes":60, "hours":3600, "days":86400}

    if(type(source) == str):
        if(source.lower() == "all"):  
 
            dt = srcid.sort_index()

            NFIlst = pd.Series([(dt.index[i+1] - (dt.index[i] + dt['len'][i])).total_seconds()/unitDict[unit] for i in range(len(dt.index)+1) 
                if(i < len(dt.index)-1)])
    
            return NFIlst
        
        elif(source.lower() == "air"):
            dt = srcid.loc[(srcid.srcID > 0) & (srcid.srcID < 2.), :]
            dt.sort_index()

            NFIlst = pd.Series([(dt.index[i+1] - (dt.index[i] + dt['len'][i])).total_seconds()/unitDict[unit] for i in range(len(dt.index)+1) 
                if(i < len(dt.index)-1)])
    
            return NFIlst
    
    else: 

        dt = srcid.loc[srcid.srcID.isin(source), :]
        dt.sort_index()

        NFIlst = pd.Series([(dt.index[i+1] - (dt.index[i] + dt['len'][i])).total_seconds()/unitDict[unit] for i in range(len(dt.index)+1) 
            if(i < len(dt.index)-1)])
    
        return NFIlst   


#------------------------------------------------------------------------------------------------------------------
# ###PERCENT TIME AUDIBLE METRICS FROM DAILYPA

def quantile_hourlyPA(dailypa, q, hour, source = "all"): # PA quantiles by hour for all sources 
    """
    Returns a quantile of percent time audible for a particular hour and source type.

    Parameters
    ----------
    dailypa: pandas dataframe representing NPS NSNSD dailypa file, formatted by soundDB library.
    q: float, a value that indicates the quantile desired, from 0.0 (minimum) to 1.0 (maximum.) 
    hour: int, hour of the day to be summarized from 0 to 23.  To summarize the entire day iterate this function with "for hour in range(0, 23):"  
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    float
    """
    
    h = str(hour).zfill(2) + "h"
    
    if(type(source) == str):
        if(source.lower() == "all"):
            return dailypa.loc[(slice(None), "Total_All"), h].quantile(q)
            
        elif(source.lower() == "air"):
            return dailypa.loc[(slice(None), "Total_1"), h].quantile(q)
    else:
        return dailypa.loc[(slice(None), [str(s) for s in source]), h].quantile(q)



def quantile_dailyPA(dailypa, q, source = "all", hour_range = [0, 23]): # PA quantiles by hour for all sources 
    """
    Returns a pandas Series containing percent time audible quantiles 
    for each hour of the day for a user-defined source type.

    Parameters
    ----------
    dailypa: pandas dataframe representing NPS NSNSD dailypa file, formatted by soundDB library.
    q: float, a value that indicates the quantile desired, from 0.0 (minimum) to 1.0 (maximum.) 
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    hour_range: list of integers. Define a time range (inclusive) for which percent time audible should be summarized.  Values are expected in 24-hour time. 
    
    
    Returns
    -------
    pandas Series of floats
    """
    
    start = str(hour_range[0]).zfill(2) + "h"
    end = str(hour_range[1]).zfill(2) + "h"

    if(type(source) == str):
        if(source.lower() == "all"):
            return dailypa.loc[(slice(None), "Total_All"), start:end].quantile(q)
            
        elif(source.lower() == "air"):
            return dailypa.loc[(slice(None), "Total_1"), start:end].quantile(q)
    else:
        if(len(source) > 1):
            raise ValueError("A single source code must be defined.")

        elif(len(source) == 1):
            return dailypa.loc[(slice(None), [str(s) for s in source]), start:end].quantile(q)



def DENABCMP_PA_exceedance(dailypa, zone, start_hour = 0, end_hour = 23, source = "all"): 
    """
    Returns the percentage of sampled hours exceeding the Denali Backcountry Management Plan PA standard.  
    Exceedance can be modified to include only specific source types or sequences of hours.

    Parameters
    ----------
    dailypa: pandas dataframe representing NPS NSNSD dailypa file, formatted by soundDB library.
    zone: str, the Denali backcountry managment plan zone:  "low", "medium", "high", or "very high".  case insensitive.
    start_hour: int, the first hour in the range to be evaluated
    end_hour: int, the last hour in the range to be evaluated
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    float, as a percentage
    """
    key = zone.lower()
    lookup = {"low": 5., "med": 15., "medium":15., "high":25., "very high":50., "v. high":50., "veryhigh":50., "v high":50.}
    
    hs = str(start_hour).zfill(2) + "h"
    hf = str(end_hour).zfill(2) + "h"
    
    if(type(source) == str):
        if(source.lower() == "all"):
            data = dailypa.loc[(slice(None), "Total_All"), hs:hf]

            p = []
            for index, values in data.iterrows():
                date = pd.Timestamp(index[0]) #convert string dates to timestamps
                hours = pd.to_timedelta(values.index.str.slice(0,2), unit='h') #this contains hourly offsets in time 
                values.index = date + hours
                p.append(values)

            d = pd.concat(p)
            return 100*(d.loc[d > lookup[key]].count()/d.count())
            
        elif(source.lower() == "air"):
            data = dailypa.loc[(slice(None), "Total_1"), hs:hf]

            p = []
            for index, values in data.iterrows():
                date = pd.Timestamp(index[0]) #convert string dates to timestamps
                hours = pd.to_timedelta(values.index.str.slice(0,2), unit='h') #this contains hourly offsets in time 
                values.index = date + hours
                p.append(values)

            d = pd.concat(p)
            return 100*(d.loc[d > lookup[key]].count()/d.count())
            
    else:
        data = dailypa.loc[(slice(None), str(source)), hs:hf]

        p = []
        for index, values in data.iterrows():
            date = pd.Timestamp(index[0]) #convert string dates to timestamps
            hours = pd.to_timedelta(values.index.str.slice(0,2), unit='h') #this contains hourly offsets in time 
            values.index = date + hours
            p.append(values)

        d = pd.concat(p)
        return 100*(d.loc[d > lookup[key]].count()/d.count())


def overall_PA(dailypa, source = "all"): 
    """
    Returns the percentage of time a source type was audible over the entire sampling period.

    Parameters
    ----------
    dailypa: pandas dataframe representing NPS NSNSD dailypa file, formatted by soundDB library.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    float, as a percentage

    """

    
    if(type(source) == str):
        if(source.lower() == "all"):
            data = dailypa.loc[(slice(None), "Total_All"), "00h":"23h"]

            p = []
            for index, values in data.iterrows():
                date = pd.Timestamp(index[0]) #convert string dates to timestamps
                hours = pd.to_timedelta(values.index.str.slice(0,2), unit='h') #this contains hourly offsets in time 
                values.index = date + hours
                p.append(values)

            d = pd.concat(p)
    
            return 100*((d/100)*3600).sum()/(len(d)*3600)
            
        elif(source.lower() == "air"):
            data = dailypa.loc[(slice(None), "Total_1"), "00h":"23h"]

            p = []
            for index, values in data.iterrows():
                date = pd.Timestamp(index[0]) #convert string dates to timestamps
                hours = pd.to_timedelta(values.index.str.slice(0,2), unit='h') #this contains hourly offsets in time 
                values.index = date + hours
                p.append(values)

            d = pd.concat(p)

            data1 = dailypa.loc[(slice(None), "Total_All"), "00h":"23h"]

            q = []
            for index, values in data1.iterrows():
                date = pd.Timestamp(index[0]) #convert string dates to timestamps
                hours = pd.to_timedelta(values.index.str.slice(0,2), unit='h') #this contains hourly offsets in time 
                values.index = date + hours
                q.append(values)

            tot = pd.concat(q)

            return 100*((d/100)*3600).sum()/(len(tot)*3600)
            
    else:

        query = [str(s) for s in source]
        data = dailypa.loc[(slice(None), query), "00h":"23h"] # if a whole number, expects source as "2" instead of "2.0"

        p = []
        for index, values in data.iterrows():
            date = pd.Timestamp(index[0]) #convert string dates to timestamps
            hours = pd.to_timedelta(values.index.str.slice(0,2), unit='h') #this contains hourly offsets in time 
            values.index = date + hours
            p.append(values)

        d = pd.concat(p)

        data1 = dailypa.loc[(slice(None), "Total_All"), "00h":"23h"]

        q = []
        for index, values in data1.iterrows():
            date = pd.Timestamp(index[0]) #convert string dates to timestamps
            hours = pd.to_timedelta(values.index.str.slice(0,2), unit='h') #this contains hourly offsets in time 
            values.index = date + hours
            q.append(values)

        tot = pd.concat(q)

        return 100*((d/100)*3600).sum()/(len(tot)*3600)



#------------------------------------------------------------------------------------------------------------------
# ### EVENT RATE, COUNT, & SATURATION METRICS FROM DAILYPA

def quantile_eventsPerDay(dailypa, q, source = "all"):
    """
    Returns a quantile of daily event rates.  Rates are calculated by source type.

    Parameters
    ----------
    dailypa: pandas dataframe representing NPS NSNSD dailypa file, formatted by soundDB library.
    q: float, a value that indicates the quantile desired, from 0.0 (minimum) to 1.0 (maximum.) 
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    float
    """
    
    
    if(type(source) == str):
        if(source.lower() == "all"):
            return dailypa.loc[(slice(None), "Total_All"), "nEvents_24Hr"].quantile(q)
            
        elif(source.lower() == "air"):
            return dailypa.loc[(slice(None), "Total_1"), "nEvents_24Hr"].quantile(q)
            
        elif(source.lower() == "low"):
            #there's some error handling for the cases where either props or helicopters 
            #were not detected during the sampling period
            try:
                props = dailypa.loc[(slice(None), str(1.2)), "nEvents_24Hr"]
            except KeyError:
                props = pd.Series()

            try:
                helos = dailypa.loc[(slice(None), str(1.3)), "nEvents_24Hr"]
            except KeyError:
                helos = pd.Series()

            d = pd.concat([props, helos], axis=1).sum(axis = 1).groupby(level = 0).sum()  #this collapses the hierarchical indexing
            return d.quantile(q)
    else:
            query = [str(s) for s in source]
            return dailypa.loc[(slice(None), query), "nEvents_24Hr"].quantile(q)


 
def total_events(dailypa, source = "all"): 
    """
    Returns a total count of events by source type.

    Parameters
    ----------
    dailypa: pandas dataframe representing NPS NSNSD dailypa file, formatted by soundDB library.
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    float
    """   
    if(type(source) == str):
        if(source.lower() == "all"):
            return dailypa.loc[(slice(None), "Total_All"), "nEvents_24Hr"].sum()
            
        elif(source.lower() == "air"):
            return dailypa.loc[(slice(None), "Total_1"), "nEvents_24Hr"].sum()
    else:
        query = [str(s) for s in source]
        return dailypa.loc[(slice(None), query), "nEvents_24Hr"].sum()



def event_saturation(dailypa, start_hour = 0, end_hour = 23, source = "all"): 
    """
    Returns a tuple of floats containing information on how many hours have events of a certain source type.
    The first value is the total number of hours with events of a certain source type.
    The second value is the percentage of all hours containing events of a certain source type.

    Parameters
    ----------
    dailypa: pandas dataframe representing NPS NSNSD dailypa file, formatted by soundDB library.
    start_hour: int, the first hour in the range to be evaluated
    end_hour: int, the last hour in the range to be evaluated
    source: str or list of floats, optional.  Which subset of srcid codes to summarize - choose either "all", "air", or specify a list of srcID codes as float.  Defaults to "all" if unspecified.
    
    Returns
    -------
    tuple: (number_of_hours_with_events_of_source_type, percentage_of_hours_with_events_of_source_type)

    """
    hs = str(start_hour).zfill(2) + "h"
    hf = str(end_hour).zfill(2) + "h"
    
    totHours = len(dailypa.index.levels[0])*24
    
    def sat(dailypa, query, hs, hf, tot=totHours):
        data = dailypa.loc[(slice(None), query), hs:hf] 

        # create a pandas series that contains hourly PA values and is indexed by datetime of each hour
        indexer = []
        hourlyPAs = []
        for index, day in data.iterrows():
            for i, hour in enumerate(day):
                indexer.append(pd.Timestamp(index[0]) + pd.to_timedelta(pd.to_numeric(day.index[i][:2]), unit='h'))
                hourlyPAs.append(float(hour))

        d = pd.Series(hourlyPAs)
        d.index = indexer

        # it doesn't matter what the value is for this calculation: it's just presence/absence
        # so overwrite the series with only the unique indices
        d = d.groupby(d.index).first() 
        
        return (len(d), 100*(d.loc[d > 0].count()/tot))
    
    
    # now use the source parameter to generate a query and calculate saturation
    if(type(source) == str):
        if(source.lower() == "all"):
            query = "Total_All"
            return sat(dailypa, query, hs, hf)
            
        elif(source.lower() == "air"):
            query = "Total_1"
            return sat(dailypa, query, hs, hf)
  
    elif(type(source) == list):    
        query = [str(s) for s in source]
        return sat(dailypa, query, hs, hf)



#------------------------------------------------------------------------------------------------------------------
# ### EVENT RATES ABOVE AMBIENT FROM LOUDEVENTS

 
def DENABCMP_events_exceedance(loudevents, zone): 
    """
    Returns the percentage of days in the sampling period that exceed the Denali Backcountry Management Plan 'events per day' standard.

    Parameters
    ----------
    loudevents: pandas dataframe representing NPS NSNSD 'loudevents' file, formatted by soundDB library.
    zone: str, the Denali backcountry managment plan zone:  "low", "medium", "high", or "very high".  case insensitive.
    
    Returns
    -------
    float, as percentage
    """
    # does loud events deal with all event types or just air?  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    key = zone.lower()
    lookup = {"low": 1., "med": 10., "medium":10., "high":25., "very high":50., "v. high":50., "veryhigh":50., "v high":50.}
    day_counts = loudevents.above.sum(axis=1)
    return 100*(len(day_counts.loc[day_counts > lookup[key]])/len(day_counts))

 
def quantile_eventRate_overAmbient(loudevents, q): #quantiles of the number of events per day over the natural ambient level
    """
    """
    return loudevents.above.sum(axis=1).quantile(q)

 
def mean_eventRate_overAmbient(loudevents): #average number of events per day over the natural ambient level
    """
    """
    return loudevents.above.sum(axis=1).mean()

 
def stdev_eventRate_overAmbient(loudevents): #standard deviation of the number of events per day over the natural ambient level
    """
    """
    return loudevents.above.sum(axis=1).std()

 
def stderr_eventRate_overAmbient(loudevents): #standard deviation of the number of events per day over the natural ambient level
    """
    """
    return loudevents.above.sum(axis=1).std()/np.sqrt(loudevents.above.sum(axis=1).count())



#------------------------------------------------------------------------------------------------------------------
# ### STANDARD ACOUSTIC EXCEEDANCE METRICS

def L90(metrics, season="Summer", weight = "A"): # the sound pressure level exceeded 10% of the time
    """
    Returns the sound pressure level that was exceeded during 90% of the sampling period.  The 10th percentile SPL.

    Parameters
    ----------
    metrics: pandas dataframe representing NPS NSNSD metrics file, formatted by soundDB library.
    season: the season for which the metric is desired:  "Summer", "Fall", "Winter", "Spring" all case-sensitive.  Defaults to "Summer".
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  

    Returns
    -------
    formatted float
    """
    w = weight.upper()
    lookup = {"A":"dBA", "T":"dBT"}
    return float("{0:.1f}".format(metrics.ambient.data.loc[season, lookup[w], "overall", "L090"]))



def Lnat(metrics, season="Summer", weight = "A"): 
    """
    Estimates the median sound pressure level that would exist in the absence of human noise.  
    Calculated as the exceedance value Lx, where x = 100*((1 + PA) / 2). 
    Thus, when the percent time audible is zero, x = 100*((1 + 0)/2) = 50, and Lnat = L50.

    Parameters
    ----------
    metrics: pandas dataframe representing NPS NSNSD metrics file, formatted by soundDB library.
    season: the season for which the metric is desired:  "Summer", "Fall", "Winter", "Spring" all case-sensitive.  Defaults to "Summer".
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  

    Returns
    -------
    formatted float
    """
    w = weight.upper()
    lookup = {"A":"dBA", "T":"dBT"}
    return float("{0:.1f}".format(metrics.ambient.data.loc[season, lookup[w], "overall", "Lnat"]))



def L50(metrics, season="Summer", weight = "A"): # the sound pressure level exceeded 50% of the time (the median SPL)
    """
    Returns the sound pressure level that was exceeded during 50% of the sampling period.  The 50th percentile SPL.

    Parameters
    ----------
    metrics: pandas dataframe representing NPS NSNSD metrics file, formatted by soundDB library.
    season: the season for which the metric is desired:  "Summer", "Fall", "Winter", "Spring" all case-sensitive.  Defaults to "Summer".
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  

    Returns
    -------
    formatted float
    """
    w = weight.upper()
    lookup = {"A":"dBA", "T":"dBT"}
    return float("{0:.1f}".format(metrics.ambient.data.loc[season, lookup[w], "overall", "L050"]))



def L10(metrics, season="Summer", weight = "A"): 
    """
    Returns the sound pressure level that was exceeded during 10% of the sampling period.  The 90th percentile SPL.

    Parameters
    ----------
    metrics: pandas dataframe representing NPS NSNSD metrics file, formatted by soundDB library.
    season: the season for which the metric is desired:  "Summer", "Fall", "Winter", "Spring" all case-sensitive.  Defaults to "Summer".
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  

    Returns
    -------
    formatted float
    """
    w = weight.upper()
    lookup = {"A":"dBA", "T":"dBT"}
    return float("{0:.1f}".format(metrics.ambient.data.loc[season, lookup[w], "overall", "L010"]))



def Leq(metrics, season="Summer", weight = "A"): 
    """
    Returns the equivalent sound pressure level, the logarithmic average of all sample values.

    Parameters
    ----------
    metrics: pandas dataframe representing NPS NSNSD metrics file, formatted by soundDB library.
    season: the season for which the metric is desired:  "Summer", "Fall", "Winter", "Spring" all case-sensitive.  Defaults to "Summer".
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  

    Returns
    -------
    formatted float
    """
    w = weight.upper()
    lookup = {"A":"dBA", "T":"dBT"}
    return float("{0:.1f}".format(metrics.ambient.data.loc[season, lookup[w], "overall", "Leq"]))



def Ldn(metrics, season="Summer", weight = "A"): 
    """
    Returns the day-night level (Ldn or DNL), which is the equivalent sound pressure level with a 10dB penalty for the hours between 22:00 and 07:00.

    Parameters
    ----------
    metrics: pandas dataframe representing NPS NSNSD metrics file, formatted by soundDB library.
    season: the season for which the metric is desired:  "Summer", "Fall", "Winter", "Spring" all case-sensitive.  Defaults to "Summer".
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  

    Returns
    -------
    formatted float
    """
    w = weight.upper()
    lookup = {"A":"dBA", "T":"dBT"}
    
    #
    hourlyLeq = metrics.hourlyMedian.data.loc[season, lookup[w], "Leq"]
    increaseNight10dB = hourlyLeq.iloc[0:7].append(hourlyLeq.iloc[22:])-10
    artificialIncrease = increaseNight10dB.append(hourlyLeq.iloc[7:22]).sort_index()
    Ldn = 10*np.log10(artificialIncrease.apply(lambda x: pow(10, x/10)).sum())

    return float("{0:.1f}".format(Ldn))


def percentTimeAbove(metrics, threshold, season="Summer", weight = "A", timeRange = "overall"): 
    """
    Returns the percentage of a record above a select SPL threshold.

    Parameters
    ----------
    metrics: pandas dataframe representing NPS NSNSD metrics file, formatted by soundDB library.
    threshold: int. The SPL threshold upon which the percentage of time above calculation is based.  Choose from 35, 45, 52, or 60 dB.
    season: the season for which the metric is desired:  "Summer", "Fall", "Winter", "Spring" all case-sensitive.  Defaults to "Summer".
    weight: str, optional.  The acoustic weighting used to calculate Lmax, either "A" or "T". Defaults to "A" if unspecified.  
    timeRange: str, optional. The summary period, either "Day", "Night", or "overall". Defaults to "overall" if unspecified.
    

    Returns
    -------
    formatted float
    """
    w = weight.upper()
    lookup = {"A":"dBA", "T":"dBT"}
    
    if threshold not in [35, 45, 52, 60]: 
        raise ValueError("Percent time above threshold must be either 35, 45, 52, or 60 dB.  Please choose another value.")
        
    thresholdStr = str(threshold)+"dB"
    
    pTA = metrics.percentTimeAbove.data.loc[season, lookup[w], timeRange, thresholdStr]

    return float("{0:.2f}".format(pTA))


#------------------------------------------------------------------------------------------------------------------
# ### STANDARD ACOUSTIC EXCEEDANCE METRICS

def Lx(nvspl, x, dBA_only=True):
    """
    Returns the exceedance percentile (Lx) for bands passed from an NVSPL file.

    Parameters
    ----------
    nvspl: pandas dataframe representing NPS NSNSD NVSPL file, formatted by soundDB library.
    x: float, the exceedance level = (100 - percentile), such that x = 10 is the 90th percentile.
    dBA_only: boolean, optional.  Whether to return a single broadband A-weighted value or all the bands passed in the NVSPL DataFrame. Defaults to a single A-weighted value if unspecified.  

    Returns
    -------
    pandas Series (default) or DataFrame

    """

    # tidy up the NVSPL file for processing by dropping the outer index
    nvspl.reset_index(level=1, drop=True)
    
    # which columns in the NVSPL file contain SPL values?
    SPLColumns = ['12.5', '15.8', '20', '25', '31.5', '40', '50', '63', '80',
                   '100', '125', '160', '200', '250', '315', '400', '500', '630', '800',
                   '1000', '1250', '1600', '2000', '2500', '3150', '4000', '5000', '6300',
                   '8000', '10000', '12500', '16000', '20000', 'dbA']

    # convert x to a quantile value
    q = (100 - x)/100

    # calculate the quantile for each band passed into the function
    bands = pd.DataFrame.from_dict({column: nvspl[column].astype(float).quantile(q) for column in nvspl if column in SPLColumns}, orient='index')
    bands.columns = ["L" + str(x)]
    bands["BANDS"] = bands.index

    # set up a template with all the columns in order
    # then merge the columns that were passed
    templateFrame = pd.DataFrame({'TODROP' : SPLColumns})
    out = pd.merge(templateFrame, bands, left_on='TODROP', right_on='BANDS', how='outer')
    
    # reset the index to the band names
    out.set_index(out['BANDS'], drop=True, inplace=True)
    out.index.name = None
    
    # then drop the template columns and unused bands
    out.drop(['TODROP', 'BANDS'], axis=1, inplace=True)
    out.dropna(inplace=True)
    
    # the default will be to return a dBA value
    if(dBA_only):
        return out.loc["dbA",:]
    
    # else return a dataframe that has all the bands!
    else:
        return out