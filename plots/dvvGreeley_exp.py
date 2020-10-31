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

from ..api import *

import pandas as pd
import numpy as np
import csv

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
    print(i[11:])
    if i[3:7] == i[11:]:
        GR.append(i)
        
NGR = []
for i in nGR_3:
    print(i[11:])
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
def main(mov_stack=None, dttname="M", components='ZZ', filterid=1,
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
    gs = gridspec.GridSpec(len(mov_stacks)+3, 1)                                #adding two more rows
    fig = plt.figure(figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    first_plot = True
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
            
            
            '''Plot only stations that intersect Greeley'''
            if pair_type == 'greeley':
                # allbut = allbut[allbut.Pairs.isin(justgreeley)]
                allbut = allbut[allbut.Pairs.isin(GR_3)]
            
            '''Plot non-Greeley stations'''
            if pair_type == 'nongreeley':
                # allbut = allbut[~allbut.Pairs.isin(justgreeley)]
                allbut = allbut[~allbut.Pairs.isin(nGR_3)]

                print(allbut['Pairs' ])
            
            
            ''''testing DTT dataframe'''
            
            testdf = []
            testdf.append(pd.read_csv('/Users/tclifford/msnoise/AC2018_3/DTT/01/001_DAYS/ZZ/2018-01-05.txt'))
            testdf.append(pd.read_csv('/Users/tclifford/msnoise/2017/1-2/DTT/01/001_DAYS/ZZ/2017-01-07.txt'))
            
            testdf = pd.concat(testdf) #this combines all df into one


            # testALL = testdf[testdf['Pairs'] == 'ALL'].copy()  #Weighted mean of all slopes, calc'd in compute_dtt.py
            
            #refine testdf by whether matches justgreeley
            
            testdf = testdf[testdf.Pairs.isin(justgreeley)]  #works!
            
            

            if first_plot == 1:
                ax = plt.subplot(gs[i])
            else:
                plt.subplot(gs[i], sharex=ax)
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
                
            
                

            tmp2 = allbut[dttname].resample('D').mean()
            plt.plot(tmp2.index, tmp2, label="mean")
            # tmp2.plot(label='mean',)
            # dvv_series.append([tmp2.index, tmp2])
            # tmp2('/Users/tclifford/msnoise/allyears.dvvcurve.csv')
            tmp2.to_csv('/Users/tclifford/msnoise/allyears/dvv_curve')
            f = open('/Users/tclifford/msnoise/allyears/dvvcurve.txt', 'a')
            f.write(str(tmp2))
            f.close()
            
            tmp3 = allbut[dttname].resample('D').median()
            # tmp3.plot(label='median')
            plt.plot(tmp3.index, tmp3, label="median")
            plt.ylabel('dv/v (%)')

            if first_plot == 1:
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.0))
                # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                #             ncol=2, borderaxespad=0.)
                # plt.legend(loc=1)
                left, right = tmp2.index[0], tmp2.index[-1]
                if mov_stack == 1:
                    plt.title('1 Day')
                else:
                    plt.title('%i Days Moving Window' % mov_stack)
                first_plot = False
            else:
                plt.xlim(left, right)
                plt.title('%i Days Moving Window' % mov_stack)

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
    plt.subplot(gs[3], sharex=ax)  #, sharex=ax)

    # plt.plot(pair1.index, pair1[dttname], label=pair)
    # plt.plot(tmp3.index,gt_curr.greeley_total)
    # gt_curr.plot()

    plt.title('Injection Volumes')
    plt.ylabel('Volume (bbls)')

    plt.plot(gt_curr.index, gt_curr.greeley_total)
    # gt['2019-07-01T00:00:00.000000000','2019-01-01T00:00:00.000000000']
    # gtt.between_times('2019-07-01T00:00:00.000000000','2019-01-01T00:00:00.000000000')
    print('TEST TEST TEST')
    print(start)
    print(end)
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
    
    plt.subplot((gs[4]), sharex=ax)
    plt.bar(eq_curr.index, eq_curr.Time, width=20)
    
    # print(start)
    # print(end)
    
    #--------------------------------------------------------------------------------------------------------------------------
    #PRECIPITATION DATA
    rain = pd.read_csv('/Users/tclifford/Documents/Greeley_Data/Greeley_precip.csv')

    rain = rain[rain.STATION == 'US1COWE0003']
    rain.DATE = pd.to_datetime(rain.DATE)
    
    
    start = '2014-06-01'
    start_dt = datetime.datetime.strptime(start, '%Y-%m-%d')
    
    end = '2019-12-31'
    end_dt = datetime.datetime.strptime(end, '%Y-%m-%d')
    
    
    # date_mask = (rain['DATE'] > start_dt) & (rain['DATE'] <= end_dt)
    rain_curr = rain.loc[(rain.index >= start) & (rain.index < end)] 
    
    # rain = rain.loc[date_mask]
    rain = rain.loc[rain_curr]

    plt.subplot((gs[5]), sharex=ax)

    # plt.ylabel('Daily Precipitation in mm')
    # plt.xlabel('Time')
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
