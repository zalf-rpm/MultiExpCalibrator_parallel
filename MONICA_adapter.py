import json
import sys
sys.path.insert(0, "C:\\Users\\stella\\Documents\\GitHub\\monica\\project-files\\Win32\\Release")
sys.path.insert(0, "C:\\Users\\stella\\Documents\\GitHub\\monica\\src\\python")
import monica_io
import zmq
import csv
import os
from datetime import date
import collections
import threading
from threading import Thread


class monica_adapter(object):
    def __init__(self, exp_maps, obslist):
        
        #for multi-experiment: create a M-2 relationship between exp_IDs and param files
        self.IDs_paramspaths = {}
        for exp_map in exp_maps:
            self.IDs_paramspaths[exp_map["exp_ID"]] = {}
            self.IDs_paramspaths[exp_map["exp_ID"]]["species"] = exp_map["species_file"]
            self.IDs_paramspaths[exp_map["exp_ID"]]["cultivar"] = exp_map["cultivar_file"]
        
        #observations data structures
        self.observations  = [] #for spotpy
        self.evaluationdates = {} #for plotting outputs
        self.obsdict = {} #for plotting outputs
        i=0 #i is a reference to match the element in result array (from spotpy)
        for record in obslist:
            self.observations.append(record["value"])
            if record["variable"] not in self.evaluationdates: #add the variable as a key
                self.evaluationdates[record["variable"]] = {}
                self.obsdict[record["variable"]] = {}
            if record["exp_ID"] not in self.evaluationdates[record["variable"]]: #add the experiment as a key
                self.evaluationdates[record["variable"]][record["exp_ID"]] = [] 
                self.obsdict[record["variable"]][record["exp_ID"]] = []
            thisdate = record["date"].split("-")#self.evaluationdates needs a date type (not isoformat)
            self.evaluationdates[record["variable"]][record["exp_ID"]].append([i, date(int(thisdate[0]), int(thisdate[1]), int(thisdate[2]))])
            self.obsdict[record["variable"]][record["exp_ID"]].append([i, record["value"]])
            i +=1

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
                    env["events"].append(unicode(record["date"]))
                    var = [unicode(record["variable"])]
                    env["events"].append(var) #TODO:Add weight here?

            position = int(exp_map["where_in_rotation"][0])

            for ws in env["cropRotation"][position]["worksteps"]:
                if ws["type"] == "Seed":
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

    def _run(self,vector, user_params):
        
        evallist=[]
        self.out={}

        #set params according to spotpy sampling. Update all the species/cultivar available
        for i in range(len(vector)):                        #loop on the vector
            par_name = user_params[i]["name"]
            for s in self.species_params:                   #loop on the species
                species = self.species_params[s]
                if par_name in species.keys():              #check for parameter existence in the dict
                    if user_params[i]["array"] == "FALSE":  #check the parameter is not part of an array   
                        species[par_name] = vector[i]
                    else:
                        arr_index= user_params[i]["array"]
                        species[par_name][int(arr_index)] = vector[i]
                else:
                    break                                   #break loop on species if the param is not there
            for cv in self.cultivar_params:                 #loop on the cultivars
                cultivar = self.cultivar_params[cv]
                if par_name in cultivar.keys():             
                    if user_params[i]["array"] == "FALSE":     
                        cultivar[par_name] = vector[i]
                    else:
                        arr_index= user_params[i]["array"]
                        if isinstance(cultivar[par_name][1], basestring): 
                            #q&d way to understand if parameters' values are in a nested array (e.g., StageTemperatureSum)
                            #a better -and more generic- solution must be found!
                            cultivar[par_name][0][int(arr_index)] = vector[i]
                            #for maize only (while calibrating TSUM1):
                            #cultivar[par_name][0][2] = vector[i] * 0.290666667
                            #cultivar[par_name][0][3] = vector[i] * 0.333333333
                            #cultivar[par_name][0][5] = vector[i] * 0.533333333
                            #for WW only (while calibrating TSUM1):
                            cultivar[par_name][0][2] = vector[i] * 1.338028169
                            cultivar[par_name][0][4] = vector[i] * 1.478873239
                        else:
                            cultivar[par_name][int(arr_index)] = vector[i]
                else:
                    break

        #launch parallel thread for the collector
        collector = Thread(target=self.collect_results)
        collector.daemon = True
        collector.start()


        #send jobs to the cluster
        for env in self.envs:
            species = self.species_params[self.IDs_paramspaths[env["customId"]]["species"]]
            cultivar = self.cultivar_params[self.IDs_paramspaths[env["customId"]]["cultivar"]]
            for crop_to_cal in env["where_in_rotation"]:
            #if the crop appears more than once in the rotation, the same params will be set
                for ws in env["cropRotation"][int(crop_to_cal)]["worksteps"]:
                    if ws["type"] == "Seed":
                        ws["crop"]["cropParams"]["species"] = species
                        ws["crop"]["cropParams"]["cultivar"] = cultivar
                        break
            
            self.socket_producer.send_json(env)
            
        
        #wait until the collector finishes
        #print "wait for collector to finish"
        collector.join()
        #print "collector finished"
        
        #build the evaluation list
        ordered_out = collections.OrderedDict(sorted(self.out.items()))
        for k, v in ordered_out.iteritems():
            for value in v:
                evallist.append(float(value))
        
        return evallist

    def collect_results(self):
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
            
            results_rec = []
            for res in rec_msg["data"]:
                results_rec.append(res["results"][0][0])
            
            self.out[int(rec_msg["customId"])] = results_rec
            received_results += 1

            #print("total received: " + str(received_results))

            if received_results == len(self.envs):
                leave = True

                