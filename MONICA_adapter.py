import json
import sys
#sys.path.insert(0, "C:\\Users\\stella\\Documents\\GitHub\\monica\\project-files\\Win32\\Release")
#sys.path.insert(0, "C:\\Users\\stella\\Documents\\GitHub\\monica\\src\\python")
import monica_io
import zmq
import csv
import os
from datetime import date
import collections
import threading
from threading import Thread
import copy
import sympy
from sympy import symbols, Eq, solve
from collections import defaultdict


class monica_adapter(object):
    def __init__(self, exp_maps, obslist, finalrun):

        #for multi-experiment: create a M-2 relationship between exp_IDs and param files
        self.IDs_paramspaths = {}
        for exp_map in exp_maps:
            self.IDs_paramspaths[exp_map["exp_ID"]] = {}
            self.IDs_paramspaths[exp_map["exp_ID"]]["species"] = exp_map["species_file"]
            self.IDs_paramspaths[exp_map["exp_ID"]]["cultivar"] = exp_map["cultivar_file"]

        #observations data structures
        self.observations = [] #for spotpy
        self.evaluationdates = {} #for plotting outputs
        self.obsdict = {} #for plotting outputs
        i = 0 #i is a reference to match the element in result array (from spotpy)
        for record in obslist:
            var_name = record["variable"]
            if "aggregation" in record.keys():
                #for plotting purposes the variable name must be unique (and match variable in collect results method)
                fromL = str(record["aggregation"][1][0])
                toL = str(record["aggregation"][1][1])
                var_name += " " + fromL + " to " + toL 
            self.observations.append(record["value"])
            if var_name not in self.evaluationdates: #add the variable as a key
                self.evaluationdates[var_name] = {}
                self.obsdict[var_name] = {}
            if record["exp_ID"] not in self.evaluationdates[var_name]: #add the experiment as a key
                self.evaluationdates[var_name][record["exp_ID"]] = []
                self.obsdict[var_name][record["exp_ID"]] = []
            thisdate = record["date"].split("-")#self.evaluationdates needs a date type (not isoformat)
            self.evaluationdates[var_name][record["exp_ID"]].append([i, date(int(thisdate[0]), int(thisdate[1]), int(thisdate[2]))])
            self.obsdict[var_name][record["exp_ID"]].append([i, record["value"]])
            i += 1

        #normalization data structure
        self.normalize = defaultdict(lambda: defaultdict(dict))
        for var in self.obsdict:
            var_indexes = []
            self.normalize[var]["max_obs_value"] = 0.01
            for exp_id in self.obsdict[var]:
                inds_values = self.obsdict[var][exp_id]
                for index_value in inds_values:
                    var_indexes.append(index_value[0]) #index in evallist
                    self.normalize[var]["max_obs_value"] = max(self.normalize[var]["max_obs_value"], index_value[1]) #max observed value
            self.normalize[var]["where"] = var_indexes
        
        #apply normalization to obs
        for anorm in self.normalize:
            n_factor = self.normalize[anorm]["max_obs_value"] / 100
            for jjj in self.normalize[anorm]["where"]:
                self.observations[jjj] /= n_factor

        self.species_params={} #map to store different species params sets avoiding repetition
        self.cultivar_params={} #map to store different cultivar params sets avoiding repetition

        #create envs
        self.envs = []
        for exp_map in exp_maps:
            with open(exp_map["sim_file"]) as simfile:
                sim = json.load(simfile)
                sim["crop.json"] = exp_map["crop_file"]
                sim["site.json"] = exp_map["site_file"]
                sim["climate.csv"] = exp_map["climate_file"]

            with open(exp_map["site_file"]) as sitefile:
                site = json.load(sitefile)

            with open(exp_map["crop_file"]) as cropfile:
                crop = json.load(cropfile)
                mycrop = exp_map["crop_ID"]
                crop["crops"][mycrop]["cropParams"]["species"][1] = exp_map["species_file"]
                crop["crops"][mycrop]["cropParams"]["cultivar"][1] = exp_map["cultivar_file"]

            env = monica_io.create_env_json_from_json_config({
                "crop": crop,
                "site": site,
                "sim": sim
            })

            #add required outputs
            for record in obslist:
                if record["exp_ID"] == exp_map["exp_ID"]:
                    if not finalrun: #output only info required by spotpy
                        env["events"].append(unicode(record["date"]))                        
                        var = [unicode(record["variable"])]
                        if "aggregation" in record.keys():
                            var = []
                            var.append(record["aggregation"])                            
                        env["events"].append(var) #TODO:Add weight here?

                    elif finalrun: #daily output for plots
                        if "daily" not in env["events"]:
                            env["events"].append(unicode("daily"))
                            env["events"].append([]) #empty list of daily variables
                            env["events"][1].append(unicode("Date"))
                        var = [unicode(record["variable"])]
                        if "aggregation" in record.keys():
                            var = record["aggregation"]
                        if var not in env["events"][1]: #avoid to ask twice the same var as out
                            env["events"][1].append(var)
                            
            position = int(exp_map["where_in_rotation"][0])

            for ws in env["cropRotation"][position]["worksteps"]:
                if ws["type"] == "Seed" or ws["type"] == "Sowing":
                    self.species_params[exp_map["species_file"]] = ws["crop"]["cropParams"]["species"]
                    self.cultivar_params[exp_map["cultivar_file"]] = ws["crop"]["cropParams"]["cultivar"]
                    break

            monica_io.add_climate_data_to_env(env, sim)
            env["customId"] = exp_map["exp_ID"]
            env["where_in_rotation"] = exp_map["where_in_rotation"]
            self.envs.append(env)

        self.context = zmq.Context()
        self.socket_producer = self.context.socket(zmq.PUSH)
        #self.socket_producer.connect("tcp://cluster2:6666")
        self.socket_producer.connect("tcp://localhost:6666")

    def run(self,args):
        return self._run(*args)

    def _run(self,vector, user_params, finalrun):

        evallist = []
        self.out = {}

        def seek_set_param(par, p_value, model_params):
            p_name = par["name"]
            array = par["array"]
            add_index = False
            if len(model_params[p_name]) > 1 and isinstance(model_params[p_name][1], basestring):
                add_index = True #the param contains text (e.g., units)
            if array.upper() == "FALSE":
                if add_index:
                    model_params[p_name][0] = p_value
                else:
                    model_params[p_name] = p_value
            else: #param is in an array (possibly nested)
                array = array.split("_") #nested array
                if add_index:
                    array = [0] + array
                if len(array) == 1:
                    model_params[p_name][int(array[0])] = p_value
                elif len(array) == 2:
                    model_params[p_name][int(array[0])][int(array[1])] = p_value
                elif len(array) == 3:
                    model_params[p_name][int(array[0])][int(array[1])][int(array[2])] = p_value
                else:
                    print "param array too nested, contact developers"
            

        #set params according to spotpy sampling. Update all the species/cultivar available
        for i in range(len(user_params)):                        #loop on the user params
            for s in self.species_params:               #loop on the species
                if user_params[i]["name"] in self.species_params[s]:
                    seek_set_param(user_params[i],
                    user_params[i]["derive_function"](vector, self.species_params[s]) if "derive_function" in user_params[i] else vector[i],
                    self.species_params[s])
                else:
                    break                                   #break loop on species if the param is not there
            for cv in self.cultivar_params:                 #loop on the cultivars
                if user_params[i]["name"] in self.cultivar_params[cv]:
                    seek_set_param(user_params[i],
                    user_params[i]["derive_function"](vector, self.cultivar_params[cv]) if "derive_function" in user_params[i] else vector[i],
                    self.cultivar_params[cv])
                else:
                    break

        #launch parallel thread for the collector
        collector = Thread(target=self.collect_results, kwargs={'finalrun': finalrun})
        collector.daemon = True
        collector.start()


        #send jobs to the cluster
        for env in self.envs:
            species = self.species_params[self.IDs_paramspaths[env["customId"]]["species"]]
            cultivar = self.cultivar_params[self.IDs_paramspaths[env["customId"]]["cultivar"]]
            for crop_to_cal in env["where_in_rotation"]:
            #if the crop appears more than once in the rotation, the same params will be set
                for ws in env["cropRotation"][int(crop_to_cal)]["worksteps"]:
                    if ws["type"] == "Seed" or ws["type"] == "Sowing":
                        ws["crop"]["cropParams"]["species"] = species
                        ws["crop"]["cropParams"]["cultivar"] = cultivar
                        break

            self.socket_producer.send_json(env)


        #wait until the collector finishes
        collector.join()
        
        #build the evaluation list for spotpy        
        if not finalrun:
            ordered_out = collections.OrderedDict(sorted(self.out.items()))
            for k, v in ordered_out.iteritems():
                for value in v:
                    evallist.append(float(value))
                        
            #apply normalization to sims
            for anorm in self.normalize:
                n_factor = self.normalize[anorm]["max_obs_value"] / 100
                for jjj in self.normalize[anorm]["where"]:
                    evallist[jjj] /= n_factor
                
            return evallist

        #return daily outputs
        elif finalrun:
            return self.out

    def collect_results(self, finalrun):
        socket_collector = self.context.socket(zmq.PULL)
        #socket_collector.connect("tcp://cluster2:7777")
        socket_collector.connect("tcp://localhost:7777")
        received_results = 0
        leave = False
        while not leave:
            try:
                rec_msg = socket_collector.recv_json()
            except:
                continue
            
            if not finalrun:
                results_rec = []
                for res in rec_msg["data"]:
                    results_rec.append(res["results"][0][0])
                self.out[int(rec_msg["customId"])] = results_rec
                #print (rec_msg["customId"], results_rec)

            elif finalrun:
                #print rec_msg["customId"]
                self.out[int(rec_msg["customId"])] = {}
                indexes_variables = []
                indexes_layeraggr =[]
                outputIds = rec_msg["data"][0]["outputIds"]
                for index in range(len(outputIds)):
                    indexes_variables.append(outputIds[index]["name"])
                    fromL_toL = [] #store info about out aggregation
                    fromL_toL.append(outputIds[index]["fromLayer"] + 1)
                    fromL_toL.append(outputIds[index]["toLayer"] + 1)
                    indexes_layeraggr.append(fromL_toL)
                results = rec_msg["data"][0]["results"]
                for res in range(len(results)):
                    variable = indexes_variables[res]
                    if indexes_layeraggr[res][0] != 0: #for variables related to soil layers
                        variable += " " + str(indexes_layeraggr[res][0]) + " to " + str(indexes_layeraggr[res][1])
                    daily_out = results[res]
                    if variable == "Date":
                        for t in range(len(daily_out)):
                            day = daily_out[t].split("-")#charts need a date type (not isoformat)
                            daily_out[t] = date(int(day[0]), int(day[1]), int(day[2]))
                    self.out[int(rec_msg["customId"])][variable] = daily_out
            
            received_results += 1

            #print("total received: " + str(received_results))

            if received_results == len(self.envs):
                leave = True
