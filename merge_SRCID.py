import pandas as pd
import numpy as np
import datetime as dt

def join_srcID_rows(df):

    '''
    Join a sequence of spectrogram annotations into a single annotation.
    '''

    new_begin = df.head(1).index[0]
    new_length = df['len'].sum()
    new_srcID = df["srcID"].values[0]
    new_L = df['Hz_L'].min()
    new_U = df['Hz_U'].max()
    new_MaxA = "{0:.01f}".format(df["MaxSPL"].max())

    # because SEL values are already normalized you just logarithmically add them
    # this gives the total energy dose
    new_SELA = "{0:.01f}".format(10*np.log10(np.power(10, df["SEL"]/10).sum()))

    # repeat for truncated values
    new_MaxT = "{0:.01f}".format(df["MaxSPLt"].max())
    new_SELT = "{0:.01f}".format(10*np.log10(np.power(10, df["SELt"]/10).sum()))

    new_user = df.head(1)["userName"].values[0]
    new_tagdate = df.tail(1)["tagDate"][0]


    joined = pd.DataFrame([new_length, new_srcID, new_L, new_U, new_MaxA, new_SELA, new_MaxT, new_SELT, new_user, new_tagdate], 
                 columns=[new_begin], index=df.columns[:-1]).T

    return joined


def merge_SRCID(src):

    import pandas as pd
    import datetime as dt

    '''
    Find SRCID annotations that break across hours, and join them to create a
    final, merged SRCID for more accurate calculations.
    '''

    # these are only the events that end during the last second of the hour
    end_at_hour = src.loc[((src.index + src["len"]).dt.minute==59)&
                          ((src.index + src["len"]).dt.second==59)]

    start_at_hour = src.loc[(src.index.minute==0)&(src.index.second==0)]

    # there are only the events that start during the first second of the hour
    starts_at_hour_break = src.index[(src.index.minute==0)&(src.index.second==0)].to_series()

    # IMPORTANT: these are the starts where there is definitely an ending one second before
    # conveniently the next two lines handle both (1) matching srcid values, 
    # and (2) lots of back-to-back-to-back hour long annotations!
    matches_end = start_at_hour.copy()
    matches_end.index = (start_at_hour.index + start_at_hour["len"] + dt.timedelta(seconds=1))


    # the .isin method is extremely helpful for working backwards to get the matched starts
    matches_start = end_at_hour[(end_at_hour.index + end_at_hour["len"] + dt.timedelta(seconds=1)).isin(start_at_hour.index)].copy()

    # these are the real break-point annotations
    cons = pd.concat([matches_start, matches_end]).sort_index().dropna()
    cons = cons.drop_duplicates(keep="first") # frankly, it doesn't matter


    # assign groups to each set of consecutive annotations 
    # it's a huge benefit to group by source type first!
    group = 1
    for srcID, source_group in cons.groupby("srcID"):

        for ts, annotation in source_group.sort_index().iterrows():

            ends = ts + annotation["len"]

            if((ts.minute!=0)&(ts.second!=0)&
               (ends.minute==59)&(ends.second==59)):

                group = group + 1
                cons.loc[ts, "group"] = group

            elif((ts.minute==0)&(ts.second==0)&
               (ends.minute!=59)&(ends.second!=59)):

                cons.loc[ts, "group"] = group
                group = group + 1

            else:
                cons.loc[ts, "group"] = group

        # at the end of the current source type, 
        # "flush" the current group
        group = group + 1


    frames = []

    # # # now we can actually perform the joining operations
    for group_number, pieces in cons.groupby('group'):

        if(len(pieces) < 2):
            # this removes single values - they don't actually need to be merged
            cons.drop(pieces.index, inplace=True)

        else:
            frames.append(join_srcID_rows(pieces))

    # here are all the joined data
    merged_breaks = pd.concat(frames)

    # # now that everything is neat and tidy, we can get the lines
    # # not representing true breaks
    no_breaks = src.loc[~src.index.isin(cons.index)]

    # final SRCID file with events across hour breaks merged
    final_src = pd.concat([merged_breaks, no_breaks])
    final_src = final_src.sort_index()

    return final_src