from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import os
import spotpy
import spotpy_setup_MONICA
import csv
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

font = {'family' : 'calibri',
    'weight' : 'normal',
    'size'   : 18}

def produce_plot(experiments, variable, ylabel='Best model simulation', xlabel='Date'):    
    #cnames = list(colors.cnames)

    plt.rc('font', **font)
    colors = ['grey', 'black', 'brown', 'red', 'orange', 'yellow', 'green', 'blue']
    n_subplots = max(2, len(experiments))
    # N subplots, sharing x axis
    width = 20
    height = n_subplots * 10
    f, axarr = plt.subplots(n_subplots, sharex=True, figsize=(width, height))
    i=0
    for exp in experiments:
        RMSE = spotpy.objectivefunctions.rmse(experiments[exp]["obs"], experiments[exp]["sims"])
        axarr[i].plot(experiments[exp]["dates"], experiments[exp]["obs"], 'ro', markersize=8, label='obs data')
        #axarr[i].plot(experiments[exp]["dates"], experiments[exp]["sims"],'-', color=colors[7], linewidth=2, label='exp ' + exp + ': RMSE=' + str(round(RMSE, 2)))
        axarr[i].plot(experiments[exp]["all_dates"], experiments[exp]["daily"],'-', color=colors[7], linewidth=2, label='exp ' + exp + ': RMSE=' + str(round(RMSE, 3)))
        axarr[i].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
        #axarr[i].set_title(str(exp))
        i +=1
    filename = variable + '.png'
    f.savefig(filename)
    text = 'A figure has been saved as ' + filename
    print(text)

#read general settings
exp_maps = []
basepath = os.path.dirname(os.path.abspath(__file__))
with open('crop_sim_site_MAP.csv') as exp_mapfile:
    reader = csv.reader(exp_mapfile)
    next(reader, None)  # skip the header
    for row in reader:
        exp_map = {}
        exp_map["exp_ID"] = row[0]
        exp_map["sim_file"] = basepath+"\\sim_files\\"+row[1]
        exp_map["crop_file"] = basepath+"\\crop_files\\"+row[2]
        exp_map["site_file"] = basepath+"\\site_files\\"+row[3]
        exp_map["climate_file"] = basepath+"\\climate_files\\"+row[4]
        exp_map["species_file"] = basepath+"\\param_files\\"+row[5]
        exp_map["cultivar_file"] = basepath+"\\param_files\\"+row[6]
        exp_map["where_in_rotation"] = row[7].split("-")
        exp_map["crop_ID"] = row[8]
        exp_maps.append(exp_map)

#read observations
obslist = [] #for envs (outputs)
with open('observations.csv') as obsfile:
    reader = csv.reader(obsfile)
    next(reader, None)  # skip the header
    for row in reader:
        if row[6].upper() == "Y":
            record = {}
            record["exp_ID"] = row[0]
            record["date"] = date(int(row[1]), int(row[2]), int(row[3])).isoformat()
            record["variable"] = row[4]
            record["value"] = float(row[5])
            if row[8] != "":
                aggregation = []
                nested_arr=[]
                nested_arr.append(int(row[8]))
                nested_arr.append(int(row[9]))
                nested_arr.append(unicode(row[10].upper()))
                aggregation.append(unicode(row[4]))               
                aggregation.append(nested_arr)
                record["aggregation"] = aggregation
            obslist.append(record) #TODO:Add weight here?

#order obslist by exp_id to avoid mismatch between observation/evaluation lists
def getKey(record):
    return int(record["exp_ID"])
obslist = sorted(obslist, key=getKey)

#read params to be calibrated
params = []
with open('calibratethese.csv') as paramscsv:
    reader = csv.reader(paramscsv)
    next(reader, None)  # skip the header
    for row in reader:
        p={}
        p["name"] = row[0]
        p["array"] = row[1]
        p["low"] = row[2]
        p["high"] = row[3]
        p["stepsize"] = row[4]
        p["optguess"] = row[5]
        p["minbound"] = row[6]
        p["maxbound"] = row[7]
        params.append(p)

spot_setup = spotpy_setup_MONICA.spot_setup(params, exp_maps, obslist)
rep = 5
results = []

sampler = spotpy.algorithms.sceua(spot_setup, dbname='SCEUA', dbformat='ram')
sampler.sample(rep, ngs=3, kstop=10)
#sampler.sample(rep, ngs=3, kstop=50, pcento=0.01, peps=0.05)
results.append(sampler.getdata())

best = sampler.status.params

with open('optimizedparams.csv', 'wb') as outcsvfile:
    writer = csv.writer(outcsvfile)        
    for i in range(len(best)):
        outrow=[]
        arr_pos = ""
        if params[i]["array"].upper() != "FALSE":
            arr_pos = params[i]["array"]        
        outrow.append(params[i]["name"]+arr_pos)
        outrow.append(best[i])
        writer.writerow(outrow)
    text='optimized parameters saved!'
    print(text)

#PLOTTING
#get the best model run
for i in range(len(results)):
    index,maximum=spotpy.analyser.get_maxlikeindex(results[i])

bestmodelrun=list(spotpy.analyser.get_modelruns(results[i])[index][0]) #Transform values into list to ensure plotting

obs_dates = spot_setup.evaluation(get_dates_dict=True)
obs_values = spot_setup.evaluation(get_values_dict=True)

#Run with optimized params
print("running simulations with optimized params")
spot_setup = spotpy_setup_MONICA.spot_setup(params, exp_maps, obslist, True)
daily_out = spot_setup.simulation(best, True)

#retrieve info for plots
print("preparing charts...")
for variable in obs_dates:
    exps = {}
    for experiment in obs_dates[variable]:
        sims = []
        obs = []
        dates = []
        for k,v in obs_dates[variable][experiment]:
            sims.append(bestmodelrun[k])
            dates.append(v)
        for k,v in obs_values[variable][experiment]:
            obs.append(v)
        exps[experiment] = {}
        exps[experiment]["dates"] = dates
        exps[experiment]["sims"] = sims
        exps[experiment]["obs"] = obs
        exps[experiment]["daily"] = daily_out[int(experiment)][variable]
        exps[experiment]["all_dates"] = daily_out[int(experiment)]["Date"]
    produce_plot(exps,variable)

    print("finished!")



