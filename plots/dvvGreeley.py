"""
This plot shows the final output of MSNoise.


.. include:: clickhelp/msnoise-plot-dvv.rst


Example:

``msnoise plot dvv`` will plot all defaults:

.. image:: .static/dvv.png



MODIFIED BY TOM CLIFFORD STARTING JUNE 25TH 2020

ADDING GREELEY EARTHQUAKE CATALOGUE AND LOCAL INJECTION INFO TO DVV PLOT FOR EASY COMPARISON
Make plotting each an optional tag
remember to add to the main msnoise plot python file

"""
import datetime
import pandas as pd
import numpy as np

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from matplotlib.dates import DateFormatter
from matplotlib import rcParams


from ..api import *

import pandas as pd
import numpy as np
import csv


#------------------------- Yearly Seasonal Curve import -----------------
import pandas as pd
import matplotlib.pyplot as plt
import obspy
import numpy as np
from scipy import signal
import datetime
from obspy import read, UTCDateTime

from matplotlib import rcParams

import matplotlib 
# matplotlib.rc('xtick', labelsize=20) 
# matplotlib.rc('ytick', labelsize=20) 

# plt.style.use('ggplot')
# plt.rcParams['figure.figsize'] = 15, 7
# rcParams["figure.subplot.hspace"] = (0.5)

curve = pd.read_csv('/Users/tclifford/msnoise/allyears/dvv_curve')

# plt.figure(figsize=(20,3))

st = obspy.core.stream.Stream()
tr = obspy.core.trace.Trace()
stats = obspy.core.trace.Stats()
tr.stats = stats
data = np.asarray(curve.M.tolist())
def convert(seconds): 
    min, sec = divmod(seconds, 60) 
    hour, min = divmod(min, 60) 
    return "%02d:%02d:%02d" % (hour, min, sec)  

times = []
for i in range(len(data)):
    # print(str(i+1))
    times.append(str(convert(i+1)))
times = [datetime.datetime.strptime(x, '%H:%M:%S') for x in times ]

tr.data = data
tr.stats.starttime = '2014-06-09T00:00:00.000000Z'
tr.stats.npts = len(data)
tr.stats.sampling_rate=1
#made interval 1 second, scale up to 1 day after filtering
#instead of freq = 1 year, do 365seconds
#period = 365s, so frequency = 1/365

# tr.plot()

tr_low = tr.copy()
tr_low.filter('lowpass', freq= 1/365)
# tr_low.plot()

tr_high = tr.copy()
tr_high.filter('highpass', freq = 1/365)
# tr_high.plot()

#now scale up to year ------------------
tr_year = tr.copy()
tr_year.stats.sampling_rate = 1/(60*60*24)
# tr_year.plot()

tr_low = tr_year.copy()
tr_low.filter('lowpass', freq= 1/(60*60*24*365))
# tr_low.plot()

tr_high = tr_year.copy()
tr_high.filter('highpass', freq= 1/(60*60*24*365))
# tr_high.plot()

#------------------------


plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 15, 7
rcParams["figure.subplot.hspace"] = (0.5)

#import and define station pairs that intersect greeley seismicity
justgreeley = pd.read_csv('/Users/tclifford/Documents/Greeley_Data/stations_justgreeley.txt')
justgreeley = justgreeley.values.tolist()
justgreeley = [val for sublist in justgreeley for val in sublist]

# then stack selected pairs to GROUPS:
groups = {
    'greeley': justgreeley
    # 'test': ['XU_GRCO_XU_GRRO']
    }

with open('/Users/tclifford/Documents/Greeley_Data/greeley_pairs.txt', newline='') as f:
    reader = csv.reader(f)
    greeley = list(reader)

print(greeley)


#make this more general later
GR_3 = greeley[0][1:-1]
nGR_3 = greeley[1][1:-1]


GR = []
for i in GR_3:
    # print(i[11:])
    if i[3:7] == i[11:]:
        GR.append(i)
        
NGR = []
for i in nGR_3:
    # print(i[11:])
    if i[3:7] == i[11:]:
        NGR.append(i)
        
#redefine GR_3, nGR as only AC
        
pd.set_option("display.max_rows", None, "display.max_columns", None)

# grl = pd.DataFrame(greeley[0][1:-1],columns=[greeley[0][0]])
# grl[greeley[1][0]]=greeley[1][1:-1]

#read in relevant dtt file, grab those that correspond to greeley group, then feed into below function

def wavg(group, dttname, errname):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    return wavg


def wstd(group, dttname, errname):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    N = len(np.nonzero(w)[0])
    wstd = np.sqrt(np.sum(w * (d - wavg) ** 2) / ((N - 1) * np.sum(w) / N))
    return wstd


def get_wavgwstd(data, dttname, errname):
    grouped = data.groupby(level=0)
    g = grouped.apply(wavg, dttname=dttname, errname=errname)
    h = grouped.apply(wstd, dttname=dttname, errname=errname)
    return g, h

# 
def main(mov_stack=None, dttname="M0", components='ZZ', filterid=1,            #modified dttname to M0
         pairs=[], showALL=False, show=False, outfile=None, pair_type=None):
    db = connect()

    start, end, datelist = build_movstack_datelist(db)

    if mov_stack != 0:
        mov_stacks = [mov_stack, ]
    else:
        mov_stack = get_config(db, "mov_stack")
        if mov_stack.count(',') == 0:
            mov_stacks = [int(mov_stack), ]
        else:
            mov_stacks = [int(mi) for mi in mov_stack.split(',')]

    if components.count(","):
        components = components.split(",")
    else:
        components = [components, ]

    low = high = 0.0
    for filterdb in get_filters(db, all=True):
        if filterid == filterdb.ref:
            low = float(filterdb.low)
            high = float(filterdb.high)
            break

    dvv_series = []
    # gs = gridspec.GridSpec(len(mov_stacks)+3, 1)                   #######################################
    gs = gridspec.GridSpec(5, 1)                                #adding two more rows

                            #adding two more rows
    fig = plt.figure(figsize=(16, 16))

    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    first_plot = True
    mov_stacks = [90, 90]
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

        alldf = pd.concat(alldf)
        print(mov_stack, alldf.head())
        if 'alldf' in locals():
            errname = "E" + dttname

            alldf[dttname] *= -100                                              #where dtt is converted to dvv? Why x100?
            alldf[errname] *= -100

            ALL = alldf[alldf['Pairs'] == 'ALL'].copy()
            allbut = alldf[alldf['Pairs'] != 'ALL'].copy()
            
            '''Plot only AC - comment out to undo'''
            if pair_type == 'AConly':
                
                allbut = allbut[allbut.Pairs.str[:7] == allbut.Pairs.str[-7:]]  #only AC
                allbut = allbut[allbut.Pairs != 'ALL'] 
            
            
            # '''Plot only stations that intersect Greeley'''
            # if pair_type == 'greeley':
            #     # allbut = allbut[allbut.Pairs.isin(justgreeley)]
            #     allbut = allbut[allbut.Pairs.isin(GR_3)]
            
            # '''Plot non-Greeley stations'''
            # if pair_type == 'nongreeley':
            #     # allbut = allbut[~allbut.Pairs.isin(justgreeley)]
            #     allbut = allbut[allbut.Pairs.isin(nGR_3)]
            #     # allbut = allbut[~allbut.Pairs.isin(nGR_3)]


                # print(allbut['Pairs' ])
            
            
            ''''testing DTT dataframe'''
            
            testdf = []
            testdf.append(pd.read_csv('/Users/tclifford/msnoise/AC2018_3/DTT/01/001_DAYS/ZZ/2018-01-05.txt'))
            testdf.append(pd.read_csv('/Users/tclifford/msnoise/2017/1-2/DTT/01/001_DAYS/ZZ/2017-01-07.txt'))
            
            testdf = pd.concat(testdf) #this combines all df into one


            # testALL = testdf[testdf['Pairs'] == 'ALL'].copy()  #Weighted mean of all slopes, calc'd in compute_dtt.py
            
            #refine testdf by whether matches justgreeley
            
            testdf = testdf[testdf.Pairs.isin(justgreeley)]  #works!
            
            

            if first_plot == 1:
                ax = plt.subplot(gs[i], )
            else:
                plt.subplot(gs[i], sharex=ax, )
            # x = {}
            # for group in groups.keys():
            #     pairindex = []
            #     for j, pair in enumerate(allbut['Pairs']):
            #         net1, sta1, net2, sta2 = pair.split('_')
            #         if sta1 in groups[group] and sta2 in groups[group]:
            #             pairindex.append(j)
            #     tmp = allbut.iloc[np.array(pairindex)]
            #     tmp = tmp.resample('D', how='mean')
            #     #~ plt.plot(tmp.index, tmp[dttname], label=group)
            #     x[group] = tmp
            #
            # tmp = x["CRATER"] - x["VOLCAN"]
            # plt.plot(tmp.index, tmp[dttname], label="Crater - Volcan")

            #modify below to plot all pairs
            #pairs of interest
            poi = [] #random pairs for now
            
                
            # for pair in pairs:
            #     print(pair)
            #     pair1 = alldf[alldf['Pairs'] == pair].copy()
            #     # testpair = testdf[testdf['Pairs'] == 'XU_LS02_XU_ORC1'].copy()

            #     print(pair1.head())
            #     plt.plot(pair1.index, pair1[dttname], label=pair)
            #     plt.fill_between(pair1.index, pair1[dttname]-pair1[errname],
            #                      pair1[dttname]+pair1[errname], zorder=-1,
            #                      alpha=0.5)
            #     pair1.to_csv('%s-m%i-f%i.csv'%(pair, mov_stack, filterid))
                
            '''below modified by Tom to plot all dvv pairs'''
            if pair_type == 'all':
            # # grab every unique pair for dftdf
                allp = alldf.Pairs.unique()
                for pair in allp:
                    print(pair)
                    pair1 = alldf[alldf['Pairs'] == pair].copy()
                    # testpair = testdf[testdf['Pairs'] == 'XU_LS02_XU_ORC1'].copy()
    
                    print(pair1.head())
                    plt.plot(pair1.index, pair1[dttname], label=pair)
                    plt.fill_between(pair1.index, pair1[dttname]-pair1[errname],
                                      pair1[dttname]+pair1[errname], zorder=-1,
                                      alpha=0.5)
                    pair1.to_csv('%s-m%i-f%i.csv'%(pair, mov_stack, filterid))
            
            
            #only AC df
            ACdf = testdf.copy()
            ACdf = ACdf[ACdf.Pairs.str[:7] == ACdf.Pairs.str[-7:]]  #only AC
            ACdf = ACdf[ACdf.Pairs != 'ALL']  #not 'All' rows
            # ACdf.mean()
            # ACdf.median()       
            
            ''' all AC pairs '''
            
            if pair_type == 'allAC':
                # '''below plots all AC pairs'''
                allp = alldf.Pairs.unique()
                for pair in allp:
                    # p = 'XU_LS01_XU_LS02'
                    if pair[:7] == pair [-7:]:    #checks if AC
                        # if 'GRRO' not in pair:  #exclude pairs if needed
                        # if pair != 'ALL':
                            print(pair)
                            pair1 = alldf[alldf['Pairs'] == pair].copy()
                            # testpair = testdf[testdf['Pairs'] == 'XU_LS02_XU_ORC1'].copy()
            
                            print(pair1.head())
                            plt.plot(pair1.index, pair1[dttname], label=pair)
                            plt.fill_between(pair1.index, pair1[dttname]-pair1[errname],
                                              pair1[dttname]+pair1[errname], zorder=-1,
                                              alpha=0.5)
                            pair1.to_csv('%s-m%i-f%i.csv'%(pair, mov_stack, filterid))
                
      
                
                
            '''below mod by Tom to plot only AC'''
            
            if pair_type == 'aconly':
                print("aconly option works")
          
                AC = []
                for pair in alldf.Pairs:
                    # print(pair)
                    if alldf.Pairs[0][0:7] == alldf.Pairs[0][-7:]:
                        AC.append(pair)
                    
                    
                ACdf = alldf[alldf['Pairs'].isin(AC)]
                for pair in ACdf:
                    print(pair)
                    pair1 = alldf[alldf['Pairs'] == pair].copy()
                    # testpair = testdf[testdf['Pairs'] == 'XU_LS02_XU_ORC1'].copy()
    
                    print(pair1.head())
                    plt.plot(pair1.index, pair1[dttname], label=pair)
                    plt.fill_between(pair1.index, pair1[dttname]-pair1[errname],
                                      pair1[dttname]+pair1[errname], zorder=-1,
                                      alpha=0.5)
                    pair1.to_csv('%s-m%i-f%i.csv'%(pair, mov_stack, filterid))
                                                       #does this append each day??
                    #append the whole series, not just the pair?
                

            if showALL:
                plt.plot(ALL.index, ALL[dttname], c='r',
                         label='ALL: $\delta v/v$ of the mean network')
                
            
            #AC or CC only
                
            AC = []
            CC = []
            for pair in alldf.Pairs:
                # print(pair)
                # print(pair)
                # if alldf.Pairs[0] == 'ALL':
                #     pass
                # if alldf.Pairs[0][0:7] == alldf.Pairs[0][-7:]:
                if pair[0:7] == pair[-7:]:
                    AC.append(pair)
                  
                    # print(alldf.Pairs[0][0:7])
                    # print(alldf.Pairs[0][-7:])
                    
                    # print(pair)
                # if alldf.Pairs[0][0:7] != alldf.Pairs[0][-7:]:
                #     CC.append(pair)
                else:
                    CC.append(pair)
                    # print('CC')
                    # print(pair)
            
                
                        
            if pair_type == 'AC':
                allbut = allbut[allbut.Pairs.isin(AC)]
                print('AUTO-CORRELATIONS ONLY')
            
            if pair_type == 'CC':
                allbut = allbut[allbut.Pairs.isin(CC)]
                print('CROSS-CORRELATIONS ONLY')

            tmp2 = allbut[dttname].resample('D').mean()
            # plt.plot(tmp2.index, tmp2, label="mean")
            # tmp2.plot(label='mean',)
            # dvv_series.append([tmp2.index, tmp2])
            # tmp2('/Users/tclifford/msnoise/allyears.dvvcurve.csv')
            # tmp2.to_csv('/Users/tclifford/msnoise/allyears/dvv_curve_5')
            # f = open('/Users/tclifford/msnoise/allyears/dvvcurve_5.txt', 'a')
            # f.write(str(tmp2))
            # f.close()
            
            tmp3 = allbut[dttname].resample('D').median()
            # tmp3.plot(label='median')
            # plt.plot(tmp3.index, tmp3, label="median")
            # plt.ylabel('dv/v (%)')
            
            
            
             #'''Plot only stations that intersect Greeley'''
            # if pair_type == 'greeley':
                # allbut = allbut[allbut.Pairs.isin(justgreeley)]
            
            if i == 0:
                GR = allbut[allbut.Pairs.isin(GR_3)]
                GRmean = GR[dttname].resample('D').mean()
                # fig.set_figheight(8)
                plt.title('%i Days Moving Window' % mov_stack + ' Greeley Pairs')
                plt.ylim([-0.5, 0.5])
                # plt.figure(figsize=(20,8)
                plt.plot(GRmean.index, GRmean, )#label="Greeley Mean")
                # plt.plot(tr_low.times('matplotlib'), tr_low.data, label='Annual', linewidth=2, alpha=0.7)

                # GRmed = GR[dttname].resample('D').median()
            # plt.plot(GRmed.index, GRmed, label="Greeley Med")

                
            #'''Plot non-Greeley stations'''
            # if pair_type == 'nongreeley':
                # allbut = allbut[~allbut.Pairs.isin(justgreeley)]
            if i == 1:
                nonGR = allbut[allbut.Pairs.isin(nGR_3)]
                nonGRmean = nonGR[dttname].resample('D').mean()
                # fig.set_figheight(8)

                plt.ylim([-0.5, 0.5])
                plt.title('%i Days Moving Window' % mov_stack + ' Non-Greeley Pairs')

                plt.plot(nonGRmean.index, nonGRmean, color='blue')#label="Non-Greeley Mean", )#alpha=0.7)
            # nonGRmed = nonGR[dttname].resample('D').median()
            # plt.plot(nonGRmed.index, nonGRmean, label="Non-Greeley Median", alpha=0.7)
            
            # sp_pair = ['']
            # spair = []
            # for i in sp_pair:
            #     spair.append('XU_'+i+'_XU_'+i)
            # single_pair = allbut[allbut.Pairs.isin(spair)]
            # special = single_pair[dttname].resample('D').mean()
            # # plt.plot(special.index, special, label=str(sp_pair), alpha=0.7)
            
            
            
            plt.ylabel('dv/v (%)')
            
            #------------------------------------------ Tom's Plotting Options --------------------------------------

            # #plot greeley intersect
            # if pair_type == 'greeley':
            # # if pair_type == 'nongreeley':
            #     # allbut = allbut[~allbut.Pairs.isin(justgreeley)]
            #     nongr = allbut[~allbut.Pairs.isin(nGR_3)]
            #     nongr2 = nongr[dttname].resample('D').mean()
            #     plt.plot(nongr2.index, nongr2, label="nongreeley mean")
            #     nongr3 = nongr[dttname].resample('D').median()
            #     plt.plot(nongr3.index, nongr3, label="nongreeley median")
            #     # plt.ylabel('dv/v (%)')
                
          
            #     # allbut = allbut[allbut.Pairs.isin(justgreeley)]
            #     gree = allbut[allbut.Pairs.isin(GR_3)]
            #     gr2 = gree[dttname].resample('D').mean()
            #     plt.plot(gr2.index, gr2, label="greeley mean", alpha = 0.3)
            #     gr3 = gree[dttname].resample('D').median()
            #     plt.plot(gr3.index, gr2, label="greeley median", alpha = 0.3)
            #     plt.ylabel('dv/v (%)')
                
            

            if first_plot == 1:
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.0))
                # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                #             ncol=2, borderaxespad=0.)
                # plt.legend(loc=1)
                left, right = tmp2.index[0], tmp2.index[-1]
                if mov_stack == 1:
                    plt.title('1 Day')
                # else:
                    # plt.title('%i Days Moving Window' % mov_stack)
                first_plot = False
            else:
                plt.xlim(left, right)
                # plt.title('%i Days Moving Window' % mov_stack)

            plt.grid(True)
            plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d %H:%M"))
            fig.autofmt_xdate()
            title = '%s, Filter %d (%.2f - %.2f Hz)' % \
                    (",".join(components), filterid, low, high)
            plt.suptitle(title)
            del alldf
            
    # print(pwd)
    # print(dvv_series)
    # f = open('/Users/tclifford/msnoise/allyears/dvvcurve.txt', 'a')
    # f.write(str(dvv_series))
    # f.close()
    #--------------------------------------------------------------------------------------------------------------------------
    #Injection Data
    inject = pd.read_csv('/Users/tclifford/Documents/Greeley_Data/Injection Volumes.csv')
    inject.set_index('Well', inplace=True)
    inject = inject.loc[:, ~inject.columns.str.contains('^Unnamed')]

    gt = inject.xs('greeley_total')
    gt.index = pd.to_datetime(gt.index)
    gtt = pd.DataFrame(gt)
    start = pd.to_datetime(allbut.index[0])
    end = pd.to_datetime(allbut.index[-1])
    # start = pd.to_datetime(tmp3.index[0])
    # end = pd.to_datetime(tmp3.index[-1])

    gt_curr = gtt.loc[(gtt.index >= start) & (gtt.index < end)]  #should only be 12 data points if it's a year
    # test = pd.concat([gt_curr]*30, ignore_index=True)
    # test = pd.DataFrame(np.repeat(gt_curr.values,12,axis=0)) #this extends index into previosu years, not what I want
    gt_daily = gt_curr.resample('D', convention='end').asfreq() #only goes up to 12/01 

    #inject.xs('greeley_total').plot()
    plt.subplot(gs[2], sharex=ax)  #, sharex=ax)

    # plt.plot(pair1.index, pair1[dttname], label=pair)
    # plt.plot(tmp3.index,gt_curr.greeley_total)
    # gt_curr.plot()

    plt.title('Injection Volumes')
    plt.ylabel('Volume (bbls)')
    # plt.ylim([0,700000])
    plt.plot(gt_curr.index, gt_curr.greeley_total)
    # gt['2019-07-01T00:00:00.000000000','2019-01-01T00:00:00.000000000']
    # gtt.between_times('2019-07-01T00:00:00.000000000','2019-01-01T00:00:00.000000000')
    # print('TEST TEST TEST')
    # print(start)
    # print(end)
    # print(tmp3.index[0])
    # print(tmp3.index[-1])
    
    #--------------------------------------------------------------------------------------------------------------------------
    #EARTHQUAKE DATA
    eq = pd.read_csv('/Users/tclifford/Documents/Greeley_Data/greeley_cat.csv')
    eq['Time'] = pd.to_datetime(eq[['Year','Month','Day','Hour','Minute','Second']])
    eq['Time'] = pd.to_datetime(eq['Time'])
    # eq['Date'] = pd.to_datetime(eq[['Year','Month','Day']])
    eq['Date'] = eq['Time'].apply( lambda eq : datetime.datetime(year=eq.year, month=eq.month, day=eq.day))
    # eq.set_index(eq["Date"],inplace=True)
   
    count = eq.groupby('Date')['Time'].count()    #count eq/day
    count = pd.DataFrame(count)

    eq_curr = count.loc[(count.index >= start) & (count.index < end)] 
    
    

    plt.subplot((gs[3]), sharex=ax)
    plt.title('Earthquake Frequency')
    plt.ylabel('Earthquakes Per Day')
    plt.bar(eq_curr.index, eq_curr.Time, width=20)
    print('EARTH EARTH EARTH')
    
    print(start)
    print(end)
    
    #--------------------------------------------------------------------------------------------------------------------------
    #EQ MOMENT MAGNITUDE replaceing previous eq plot 
    # from pandas.plotting import register_matplotlib_converters
    # register_matplotlib_converters()
    # import datetime as dt
    # import matplotlib.dates as mdates
    # import matplotlib.dates as dates
    # import matplotlib.ticker as ticker
    
    # gr = pd.read_csv('/Users/tclifford/Documents/Greeley_Data/greeley_cat.csv')

    # gr['time'] = (gr['Year'].astype(str) + '-' + 
    #               gr['Month'].astype(str) + '-' + 
    #               gr['Day'].astype(str) + 'T' + 
    #               gr['Hour'].astype(str) + ':' + 
    #               gr['Minute'].astype(str) + ':' + 
    #               gr['Second'].astype(str))
    
    
    # df['time'] = df['time'].astype('datetime64')
    # gr['time'] = gr['time'].astype('datetime64')
    
    # df = gr

    # start_date = '2014-06-01'
    # end_date = '2020-01-01'
    
    # mask = (df['time'] > start_date) & (df['time'] <= end_date)
    
    # df = df.loc[mask]

    # df3 = df#[df.Ml >= 3]

    
    # a = df3['Ml'].tolist() #df to an array/list

    # print(df3['Ml'])
    
    # moment_array = []
    
    # for i in a:
    #         #moment_array.append((pow(10,(i + 10.73)*1.5))*1e-7) #S&W formula calc moment then convert Dyn-cm to Nm
    #         moment_array.append((pow(10,((1.5*i + 9.05))))) #Kanamori formula in Nm

    
    
    # cum_array = np.nancumsum(moment_array) #plot cumulative moment instead of individual moments
    # #cum array only 544 long - why?
    
    # print (cum_array)
    
  
    # y1 = df.time.tolist()
    

    # # plt.subplot((gs[4]), sharex=ax)

    # plt.figure(figsize=(10,8))
    # plt.plot(y1, cum_array, marker='x') #changed plt.scatter to plt.plot to get lines between points
    # plt.title('Cumulative Moment Release', fontsize = 25)
    # plt.xlabel('Year', fontsize = 20)
    # plt.ylabel('Cumulative Moment (Newton-Meters)', fontsize = 20)
    # # plt.yscale('log') #Tom, also make your y-axis logarithmic, not really sure what log and sci together do considering the error message
    # plt.show()

    
     #--------------------------------------------------------------------------------------------------------------------------
    #PRECIPITATION DATA
    rain = pd.read_csv('/Users/tclifford/Documents/Greeley_Data/Greeley_precip.csv')
    
    rain = rain[rain.STATION == 'US1COWE0003']
    rain.DATE = pd.to_datetime(rain.DATE)
    
    
    # start = '2014-06-01'
    start = str(start)[:10]
    start_dt = datetime.datetime.strptime(start, '%Y-%m-%d')

    # end = '2019-12-31'
    end = str(end)[:10]
    end_dt = datetime.datetime.strptime(end, '%Y-%m-%d')
    
    
    
    
    
    date_mask = (rain['DATE'] > start_dt) & (rain['DATE'] <= end_dt)
    
    rain = rain.loc[date_mask]



    plt.subplot((gs[4]), sharex=ax)

    # plt.ylabel('Daily Precipitation in mm')
    # plt.xlabel('Time')
    plt.title('Rainfall')
    plt.ylabel('Daily Precipitation (mm)')
    plt.plot(rain.DATE, rain.PRCP)
    print('rain rain rain')
    
    
    if outfile:
        if outfile.startswith("?"):
            if len(mov_stacks) == 1:
                outfile = outfile.replace('?', '%s-f%i-m%i-M%s' % (components,
                                                                   filterid,
                                                                   mov_stack,
                                                                   dttname))
            else:
                outfile = outfile.replace('?', '%s-f%i-M%s' % (components,
                                                               filterid,
                                                               dttname))
        outfile = "dvv " + outfile
        print("output to:", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()


if __name__ == "__main__":
    main()
