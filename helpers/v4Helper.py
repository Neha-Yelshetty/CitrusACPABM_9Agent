import json
import os
from pickle import FALSE
import psycopg2
import multiprocessing as mp
import tqdm
import random
import string
import math
import datetime
import copy 

econConfig = None 
with open("configs/econConfig.json", "r") as read_file:
    econConfig = json.load(read_file)

bioconfig_file = None 
with open("configs/bioConfig.json", "r") as read_file:
    bioconfig_file = json.load(read_file)

#List to hold all of our triples
masterConfigList = []
setid_dict = {}
cluster_name = "TESTER"


def runInstance(config):
    biofileName = config[1]["fileName"][:-4]
    econFileName = config[0]["outputFilename"][:-4]
    print("biofileName :{}".format(biofileName))
    print("econFileName :{}".format(econFileName))
    os.makedirs(f'configs/{cluster_name}/{config[2]}', exist_ok=True)
    os.makedirs(f'output/{cluster_name}/{config[2]}', exist_ok=True)
    config[0]["outputFilename"] = f"output/{cluster_name}/{config[2]}/{econFileName}.csv"
    config[1]["fileName"] = f"output/{cluster_name}/{config[2]}/{biofileName}.csv"
    config[0]["experimentID"] = 5
    with open(f"configs/{cluster_name}/{config[2]}/{biofileName}.json", "w") as write_file:
        json.dump(config[1], write_file)
    with open(f"configs/{cluster_name}/{config[2]}/{econFileName}.json", "w") as write_file:
        json.dump(config[0], write_file) 
    os.system(f"./serial.out \"configs/{cluster_name}/{config[2]}/{econFileName}.json\" \"configs/{cluster_name}/{config[2]}/{biofileName}.json\"")




def appendMasterTriples(econ, biocons, folder, master, numSets=0, nameSuffix="",noofiteration = -1):
    if numSets==0:
        #rndBio = random.randint(0, len(biocons)-1)
        bio_copy = copy.deepcopy(biocons)
        econ_copy = copy.deepcopy(econ) 
        econ_copy["outputFilename"] = f"{noofiteration}_econ.csv"
        bio_copy["fileName"] = f"{noofiteration}_bio.csv"
        master.append([econ_copy, bio_copy, folder])
    


for x in range(0,100):
    econ_copy = copy.deepcopy(econConfig)
    #econ_copy["strategyFlags"] = "0,0,0,0"
    #econ_copy["strategyParameters"] = f"0,0,0,0,0,0,0,0,0"
    #econ_copy["noactionread_file_no"] = x

    bio_copy = copy.deepcopy(bioconfig_file)
    invasionDays_lst = "81,446".split(",")  # Split into a list of separate values
    bio_copy["invasionDays"] = ','.join(invasionDays_lst)  # Store it properly

    # Pick a random value and replicate it for each invasion day
    random_value = random.randint(1, 7)
    invasionModality_lst = [random_value] * len(invasionDays_lst)

    # Convert list to string
    invasionModalities = ','.join(map(str, invasionModality_lst))

    invasionGrove = random.randint(1,2)
    print(invasionGrove)
    print(invasionModality_lst)  
    print(invasionModalities)  
    bio_copy["invasionModalities"] = invasionModalities
    appendMasterTriples(econ_copy, bio_copy, f"26grove_varation",  masterConfigList, numSets=0,noofiteration = x)

'''for e in [700,800,900]: #600,
    for x in range(0,10):
        econ_copy = copy.deepcopy(econConfig)
        econ_copy["strategyFlags"] = "0,1,0,0"
        econ_copy["strategyParameters"] = f"5,5,0,0,0.1,{e/1000},0.246,5,0"
        econ_copy["noactionread_file_no"] = x
        no_ofiter = str(x)+"__"+str(e)

        bio_copy = copy.deepcopy(bioconfig_file)
        invasionDays_lst = ["81,446"]
        bio_copy["invasionDays"] = ','.join([x for x in invasionDays_lst])
        invasionModality_lst = random.choices(range(1,7),k=len(invasionDays_lst))
        invasionModalities = ','.join([str(x) for x in invasionModality_lst])
        bio_copy["invasionModalities"] =  invasionModalities
        appendMasterTriples(econ_copy, bio_copy, f"sprayVariations",  masterConfigList, numSets=0,noofiteration = no_ofiter)'''


basePath = f"/home/instr1/repo/test/output/{cluster_name}"
print("Starting data loading...")

pool = mp.Pool(mp.cpu_count())

for _ in tqdm.tqdm(pool.imap_unordered(runInstance, masterConfigList), total=len(masterConfigList)):
    pass 
pool.close()
