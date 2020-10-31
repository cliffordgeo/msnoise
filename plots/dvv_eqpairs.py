#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 11:49:06 2020

@author: tclifford



Custom station averages
"""

for i, mov_stack in enumerate(mov_stacks):
        current = start
        first = True
        alldf = []
        while current <= end:
            for comp in components:
                day = os.path.join('DTT', "%02i" % filterid, "%03i_DAYS" %
                                   mov_stack, comp, '%s.txt' % current)
                if os.path.isfile(day):
                    df = pd.read_csv(day, header=0, index_col=0,
                                     parse_dates=True)
                    alldf.append(df)
            current += datetime.timedelta(days=1)
        if len(alldf) == 0:
            print("No Data for %s m%i f%i" % (components, mov_stack, filterid))
            continue


#get above code to just grab the station pairs I want

#sample dtt /Users/tclifford/msnoise/2017/0.2-0.5/DTT/01/001_DAYS/ZZ/2017-01-02.txt

import os
import pandas as pd

df = pd.read_csv('/Users/tclifford/msnoise/2017/0.2-0.5/DTT/01/001_DAYS/ZZ/2017-01-02.txt')

df.Date
df.Pairs

if df.Pairs matches eq_pair, append to all:
    for i in df.Pairs:
        if i in eq_pairs:
            eq_pairs[append]
        else:
            non_eq_pairs.append(df.Pairs[i])
            
            
#https://www.geeksforgeeks.org/selecting-rows-in-pandas-dataframe-based-on-conditions/

eq_pair = ['XU_BR01_XU_LS02']  #expand to the stations I want

eq_data = df[df['Pairs'].isin(eq_pair)]
