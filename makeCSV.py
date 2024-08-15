#imports
from astropy.table import Table as tbl
from concurrent.futures import ThreadPoolExecutor
from queue import Queue
import pandas as pd
import threading
import os

#<SETUP>
#paths
path = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/images/simple_model/'
saveDirectory = '/work/jyc29/stamps/largerPoolTest2'
scaLocTable = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/Roman_TDS_obseq_11_6_23_radec.fits'
truthTablePath = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/truth/'

numCores = 10

#</SETUP>

#useful methods
def get_pointing(row):
    return (row//18)
def get_sca(row):
    return ((row)%18+1)
def get_filter(row):
    return t['filter'][row//18]

def addCSV(file): #queues new line to add to csv file
    ls = file.split(',')

    x = float(ls[0])
    y = float(ls[1])
    k = int(ls[2][:-4])

    filt = get_filter(k)
    pointing = get_pointing(k)
    sca = get_sca(k)

    #open and check +-5 x,y pixel in truth file

    truth = pd.read_csv(f'{truthTablePath}{filt}/{pointing}/Roman_TDS_index_{filt}_{pointing}_{sca}.txt', sep=' ', skipinitialspace=True, skiprows=1, names=['object_id', 'ra', 'dec', 'x', 'y', 'realized_flux', 'flux', 'mag', 'obj_type'])
    narrow = truth.loc[((truth['x']<x+15) & (truth['x']>x-15)) & ((truth['y']<y+15) & (truth['y']>y-15))]

    bogus = 1 if (len(narrow.loc[narrow['obj_type'] == 'transient'])>0) else 0 #SNe injected or not
    SNeID = int(narrow.loc[narrow['obj_type'] == 'transient']['object_id'].values[0]) if bogus else -99 #if injected, what is the ID

    #take magnitude from truth tables
    if bogus:
        mag = float(narrow.loc[narrow['obj_type'] == 'transient']['mag'].values[0])
    else:
        mag = -99
    #turn into a string and write into file
    q.put(f'{SNeID},{filt},{bogus},{mag},"{file}"')

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

#create csv file (overwrites)
pd.DataFrame(columns=['SNeID','filter','bogus','magnitude','name']).to_csv(saveDirectory+'/!stampID.csv', index=False) #name can be anything, !stampID.csv was abitrary

#table data
t = tbl.read(scaLocTable) #read sca-radec table as an astropy table

csv = open(saveDirectory+'/!stampID.csv','r') #opens with read-only access
csv2 = open(saveDirectory+'/!stampID.csv','a') #opens with append-only access

q = Queue() #creates a thread-safe queue

def writeCSV(q): #writes new line to csv file whenever
    while True:
        csvline = q.get() #reads from queue
        csv2.write(csvline+'\n') #appends new line to csv
        csv2.flush() #refreshes csv (saves new line)
        q.task_done()

writeThread = threading.Thread(target=writeCSV,args=(q,)) #creates writeCSV thread
writeThread.daemon = True
writeThread.start() #starts writeCSV thread

with ThreadPoolExecutor(max_workers=numCores-2) as executor: #creates a threadpool for generating lines to be added
    executor.map(addCSV, os.listdir(saveDirectory)) #maps addCSV function to every stamp in directory
q.join() #waits for queue to finish

df = pd.read_csv(saveDirectory+'/!stampID.csv') #read CSV into dataframe
df['ID'] = [i for i in range(len(df))] #adds an ID value for each row/stamp of CSV
df.to_csv(saveDirectory+'/!stampID.csv', index=False) #turns it back into a CSV

print('done')
