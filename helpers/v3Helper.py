import json
import os
from pickle import FALSE
import multiprocessing as mp
import random
import string
import math
import pandas as pd
import decimal
from django.db import connection
from .mysql_connection import close_and_init_ssh_tunnel, close_ssh_tunnel,init_ssh_tunnel
import atexit
import time
import send2trash
from .models import ProfitInformation
from .serializers import ProfitInformationSerializer
from .storeprocedurecall import insert_data
import threading

def id_generator(size=30, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def execute_sql(sql_template, values_list):
    insert_data(sql_template, values_list)


#cnx_a = psycopg2.connect(host="localhost", user="postgres", password="CitrusABM21", database="netbenefits")
#cursor = cnx_a.cursor()
    
cluster_name = "TESTER"
#cluster_id = createCluster(cnx_a,cursor,cluster_name)
global_bioid = 1


#Open default config file
config_file = None 

import datetime
import copy 
configList = []
with open("Agentbasedmodel/configs/bioConfig.json", "r") as read_file:
        config_file = json.load(read_file)

econConfig = None 
with open("Agentbasedmodel/configs/econConfig.json", "r") as read_file:
    econConfig = json.load(read_file)

#List to hold all of our triples
masterConfigList = []

setid_dict = {}
def appendMasterTriples(econ, biocons, folder, master, numSets=0, nameSuffix=""):
    print("biocons,econ : {} {} " .format(biocons,econ))
    if numSets==0:
        rndBio = random.randint(0, len(biocons)-1)
        bio_copy = copy.deepcopy(biocons[rndBio])
        econ_copy = copy.deepcopy(econ) 
        rndmid = id_generator()
        if len(nameSuffix) == 0:
           rndmid = nameSuffix + rndmid
        econ_copy["outputFilename"] = f"{rndmid}_econ.csv"
        bio_copy["fileName"] = econ_copy["outputFilename"].replace("econ", "bio")
        print(f"bio_copy : {bio_copy[0]}")
        master.append([econ_copy, bio_copy, folder])
    else: 
        #for ns in range(numSets):
            #for bio in biocons:
        econ_copy = copy.deepcopy(econ)
        bio_copy = copy.deepcopy(biocons[0])
        rndmid = id_generator()
        if len(nameSuffix) == 0:
            rndmid = nameSuffix + rndmid
        bio_copy["fileName"] = f"{rndmid}_bio.csv"
        econ_copy["outputFilename"] = bio_copy["fileName"].replace("bio", "econ")
        master.append([econ_copy, bio_copy, folder])
                


def createeconconfig(configtype):
    econ_copy = copy.deepcopy(econConfig)
    for wl in range(-0.0002,-0.2,0,0.0002,0.2,0.3):
        econ_copy["wl"]=wl
        for wh in range(-0.0002,-0.2,0,0.0002,0.2,0.3):
            econ_copy["wh"]=wh
            for tyn in range(1,2,3):
                econ_copy["typeofnetwork"]=tyn


                if configtype == 1:
                    # #NO ACTION BASELINE
                    for x in range(0,1):
                        stf = ""
                        stp = ""
                        for k in range(0,8):
                            stf = 0,0,0,0,0 + ";"
                            stp = f"0,0,0,0,0,0,0,0" + ";" 
                        econ_copy["strategyFlags"] = stf
                        econ_copy["strategyParameters"] = stp
                        appendMasterTriples(econ_copy, configList, f"noAction_baseCase",  masterConfigList, numSets=1)
                else if configtype == 2:
                    #SPRAY TEST
                    for x in range(0,1):
                        for e in [600,700,800,900]:
                            stf = ""
                            stp = ""
                            for k in range(0,8):
                                stf = 0,1,0,0,0 + ";"
                                stp = f"5,5,0,0,0.1,{e/1000},0.246,5,0" + ";" 
                            econ_copy["strategyFlags"] = stf
                            econ_copy["strategyParameters"] = stp
                            appendMasterTriples(econ_copy, configList, f"sprayVariations",  masterConfigList, numSets=1)
                else if configtype == 3:
                    # # # #ROGUE TEST
                    for x in range(0,1):
                        for radius in [0]: #20,40,60
                            for frequency in [45,70,105,135]: #1,6,12,18,24
                                for cost in [0]: #3,6,9
                                    for thershold in [0.1,0.2,0.3,0.4,0.5]:
                                        stf = ""
                                        stp = ""
                                        for k in range(0,8):
                                            stf = 1,0,0,0,0 + ";"
                                            stp = f"{cost},0.25,{frequency},{radius},{thershold},0,5,5,0" + ";"
                                        econ_copy["strategyFlags"] = stf
                                        econ_copy["strategyParameters"] = stp
                                        appendMasterTriples(econ_copy, configList, f"rogueVariations", masterConfigList, numSets=1)
                else if configtype == 4:
                    # # # # #SPRAY AND ROGUE
                    for x in range(0,1):
                        for e in [600,700,800,900]:
                            for frequency in [45,70,105,135]: #1,6,12,18,24
                                for radius in [0]: #1,8,10,20,40
                                    for cost in [0]: #3,6,9
                                        for thershold in [0.1,0.2,0.3,0.4,0.5]:
                                            stf = ""
                                            stp = ""
                                            for k in range(0,8):
                                                stf = 1,1,0,0,0 + ";"
                                                stp = f"{cost},0.25,{frequency},{radius},{thershold},{e/1000},0.246,5,0" + ";"
                                            econ_copy["strategyFlags"] = stf
                                            econ_copy["strategyParameters"] = stp
                                            appendMasterTriples(econ_copy, configList, f"rogueSprayVariations", masterConfigList, numSets=1)

def runInstance(config):
    biofileName = config[1]["fileName"][:-4]
    econFileName = config[0]["outputFilename"][:-4]
    os.makedirs(f'Agentbasedmodel/configs/{cluster_name}/{config[2]}', exist_ok=True)
    os.makedirs(f'Agentbasedmodel/output/{cluster_name}/{config[2]}', exist_ok=True)
    config[0]["outputFilename"] = f"Agentbasedmodel/output/{cluster_name}/{config[2]}/{econFileName}.csv"
    config[1]["fileName"] = f"Agentbasedmodel/output/{cluster_name}/{config[2]}/{biofileName}.csv"
    config[0]["experimentID"] = 1
    with open(f"Agentbasedmodel/configs/{cluster_name}/{config[2]}/{biofileName}.json", "w") as write_file:
        json.dump(config[1], write_file)
    with open(f"Agentbasedmodel/configs/{cluster_name}/{config[2]}/{econFileName}.json", "w") as write_file:
        json.dump(config[0], write_file) 
    os.system(f"Agentbasedmodel/serial.out \"Agentbasedmodel/configs/{cluster_name}/{config[2]}/{econFileName}.json\" \"Agentbasedmodel/configs/{cluster_name}/{config[2]}/{biofileName}.json\"")

createeconconfig(2):
revisedMaster = []
for config in masterConfigList:
    #expid = getExperimentID(cnx_a, setid_dict[config[2]], config[0], config[1], global_bioid)
    newConfig = copy.deepcopy(config)
    #newConfig.append(expid)
    revisedMaster.append(newConfig)


pool = mp.Pool(mp.cpu_count())

for _ in tqdm.tqdm(pool.imap_unordered(runInstance, revisedMaster), total=len(revisedMaster)):
    pass 
pool.close()



