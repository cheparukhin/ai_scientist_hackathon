#Python version: 3.6 and above


#Customized example : query server for all model data, based on a smiles-string, and output to out.csv.
#PLEASE NOTE : As seen below, do add quotation marks if you include a SMILES string, since otherwise the = will be mis-interpreted. Likewise,use quotes around the whole query if split drugnames (two words) occur within
#python simple_api.py -t smiles -m ALL_MODELS -o out.csv "CCC(=C(C1=CC=CC=C1)C2=CC=C(C=C2)OCCN(C)C)C3=CC=CC=C3"
#python simple_api.py aspirin,sertraline

import requests;    #Server interaction
import time;        #Timers for waiting time estimation
import argparse;    #Command line switch handling
import json;        #Structured data formatting and transfer
import sys;
import pandas as pd;
import ssl;

ssl._create_default_https_context = ssl._create_unverified_context

#List of all computationally intensive models. Either manually pick from this list, or specify to use ALL_MODELS directly
ALL_MODELS = "dili neuro nephro respi cardio carcino immuno mutagen cyto bbb eco clinical nutri nr_ahr nr_ar nr_ar_lbd nr_aromatase nr_er nr_er_lbd nr_ppar_gamma sr_are sr_hse sr_mmp sr_p53 sr_atad5 mie_thr_alpha mie_thr_beta mie_ttr mie_ryr mie_gabar mie_nmdar mie_ampar mie_kar mie_ache mie_car mie_pxr mie_nadhox mie_vgsc mie_nis CYP1A2 CYP2C19 CYP2C9 CYP2D6 CYP3A4 CYP2E1"

#PROGRAM LOGIC
#---------------------------------
parser=argparse.ArgumentParser(description="Query the VirtualTaste API")
parser.add_argument("-t","--mtype", help="Default: name. specify 'name' (pubchem search using the compound name) or 'smiles' (canonical smiles) as your input data type",type=str,choices=['name','smiles'],default="name")
parser.add_argument("-m","--models",help="Default : 'acute_tox tox_targets'. specify models and data points to compute. You can use all additional models using ALL_MODELS, but be mindful it may incur high calculation times",type=str,default="acute_tox tox_targets")
parser.add_argument("searchterms", help="The actual strings to search for (comma-separated). May be pubchem name (default) or SMILES (use -t smiles). ",type=str)
parser.add_argument("-o","--outfile",help="Default : results.csv. specify where to output your retrieved data (csv format)",type=str,default="results.csv")
parser.add_argument("-r","--retry",help="Default : 5. Retry limit to attempt to retrieve data. Increase this in case of large, unpredictable queries",type=int,default=5)
parser.add_argument("-q","--quiet",help="Suppress all non-error output from the script",action="store_true")
args=parser.parse_args()

input_type=args.mtype;
models=args.models.split(',')
if "ALL_MODELS" in models[0]:
    models[0] = models[0].replace("ALL_MODELS", "")
    models[0]=models[0]+ALL_MODELS
searchterms=args.searchterms.split(',')
outfile=args.outfile;
retry_limit=args.retry;
quiet=args.quiet;

task_id_list=[]

def log(msg):
  if (not quiet):
    print (msg)

def request_data(inputs):
    log("Enqueing request "+inputs+", with models :")

    log(models)
    r=requests.post("https://tox.charite.de/protox3/src/api_enqueue.php",data={'input_type':input_type,'input':inputs,'requested_data':json.dumps(models)}) #encode array, the rest are single strings
    if (r.status_code==200): #Data response

        #Data query response. Add response to our task_id_list
        log("Recieved qualified response with id")
        log(r.text)
        task_id_list.append(r.text) 

        #Set up wait time before next query
        if 'Retry-After' in r.headers:
          wait_time=int(r.headers['Retry-After'])+1 #Wait depending on how long the server thinks it needs
        else:
          wait_time=5 #Something went wrong with transmission, just wait 10s by default

        log("Waiting for "+str(wait_time)+" s till next request")

        time.sleep(wait_time)

    elif (r.status_code==429): #Too many requests. Slow down/wait
            log("Server responds : Too many requests. Slowing down query slightly")
            if 'Retry-After' in r.headers:
                wait_time=r.headers['Retry-After']+1 #Wait depending on how long the server thinks it needs
            else:
                wait_time=5 #Something went wrong with transmission, just wait 10s by default

            log("Waiting for "+wait_time+" s till next request")
            time.sleep(wait_time)

    elif (r.status_code==403): #Daily quota exceeded. Aborting operation for now.
            print ("Daily Quota Exceeded. Terminating process.")
            exit(0);

    else: #Server gone away or different issues, outputting here
        print ("ERROR : Server issue")
        print (r.status_code,r.reason)

for inputs in searchterms:
    time.sleep(2)
    request_data(inputs)

#Once done, retrieve data from the server. Repeat queries every 30s until all data has been retrieved. Append to an outfile if given
log ("All queries have been enqueued. Starting result retrieval...")
time.sleep(10);

while task_id_list: #As long as we still have elements to get :
  response_list=[]
  first_response=1
  i=0
  data=pd.DataFrame()
  for task_id in task_id_list:
    log("Asking for "+task_id)
    tox_data=pd.DataFrame()
    response=pd.DataFrame()
    target_data=pd.DataFrame()
    #print(task_id)
    r=requests.post("https://tox.charite.de/protox3/src/api_retrieve.php",data={'id':task_id})
    if (r.status_code==200): #Data response
      if (r.text==""):
        print ("Warning : Empty response")
      else:
        response_list.append(task_id)
        if "acute_tox" in models[0]:
            tox_data=pd.read_csv("https://tox.charite.de/protox3/csv/"+task_id+"_tox_class.csv", sep='\t')
            tox_data.insert(loc=0, column='type', value="acute toxicity")
        if len(set(ALL_MODELS.split()) & set(models[0].split()))>0:
            response=pd.read_csv("https://tox.charite.de/protox3/csv/"+task_id+"_result.csv", sep='\t')
            response.drop(response.columns[0], axis=1, inplace=True)
            response.insert(loc=0, column='type', value="toxicity model")
        if "tox_targets" in models[0]:
            target_data=pd.read_csv("https://tox.charite.de/protox3/csv/"+task_id+"_tox_targets.csv", sep='\t')
            target_data.insert(loc=0, column='type', value="toxicity target")
        response=pd.concat([tox_data, response, target_data], ignore_index=True)
        response.insert(loc=0, column='input', value=searchterms[i])
        data=pd.concat([data, response], ignore_index=True)
        i=i+1
    elif (r.status_code==404): #Not found, not computed or finished yet. Do nothing
      if (first_response):
        print ("No response yet. Likely cause: computation unfinished (retrying...)")
        first_response=0
    else: #Other codes are not permitted
      print ("Unexpected return from server")
      print (r.status_code,r.reason)
      sys.exit();

  task_id_list=[item for item in task_id_list if item not in response_list] #Remove all found id's
  if (task_id_list):
    print ("Some queries still pending. Retrying in 30s")
    sys.exit();
    time.sleep(30); #Wait 30s before another run if there's still work to do

data.to_csv(outfile, sep='\t', encoding='utf-8')

print ("Completed all operations. Your results are in "+outfile)
