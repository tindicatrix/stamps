#imports
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
from matplotlib import pyplot as plt
from astropy.table import Table as tbl
from collections import defaultdict
from multiprocess import Pool 
from concurrent.futures import ThreadPoolExecutor
from queue import Queue
import pandas as pd
import numpy as np
import threading, itertools, os

#<SETUP>———————————————————————————————————————————————— edit the following before running
numTotalPointings = 57365                                                                       #total number of pointings w/ 18 scas
numCores = 10                                                                                   #number of cpu cores allotted (~10 seemed optimal)

#paths
path = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/images/simple_model/'                         #where are the sims
saveDirectory = '/work/jyc29/stamps/largerPoolTest2'                                            #directory where you want to save the stamps, does not need to be preexisting
scaLocTable = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/Roman_TDS_obseq_11_6_23_radec.fits'    #where is the big sca table
truthTablePath = '/cwork/mat90/RomanDESC_sims_2024/RomanTDS/truth/'                             #where are the truth tables
outputCat = 'output.cat'                                                                        #where/what is the SExtractor catalog file
CSVfile = saveDirectory+'/!stampID.csv'                                                         #where/what is the CSV file

#also edit line 87 where the files are opened (if needed)

#</SETUP>————————————————————————————————————————————————

#useful methods
def dist(radecPairs, ra, dec): #returns an array with the distance between the given radec and all sca centers
        return np.sqrt((radecPairs[:,0]-ra)**2+(radecPairs[:,1]-dec)**2) 

def bound(scaAr): #calculate distance bound (max between distance from center of scas2&12 or 3&11)
    bound = np.fmax(np.sqrt((scaAr[0,:,1]-scaAr[0,:,11])**2+(scaAr[1,:,1]-scaAr[1,:,11])**2)/2,np.sqrt((scaAr[0,:,2]-scaAr[0,:,10])**2+(scaAr[1,:,2]-scaAr[1,:,10])**2)/2) 
    bound = np.atleast_2d(bound).transpose() #make it a 2d array and vertical
    return np.tile(bound,18).reshape(-1) #tile it so that there's 18 copies in each row and reshape so that its all 1 column

def get_pointing(row:int)->int:
    return (row//18)
def get_sca(row:int)->int:
    return ((row)%18+1)
def get_filter(row:int)->str:
    return t['filter'][row//18]

def getCorners(check, scaAr, radecPairs): #returns [[ra1,dec1],[ra2,dec2],[ra3,dec3]] w/ array shape (for example): (3,2,2560)
    check1 = get_pointing(check)
    checkshort = scaAr[:,check1]

    raSca2 = checkshort[0,:,1]
    raSca3 = checkshort[0,:,2]
    raSca12 = checkshort[0,:,11]

    decSca2 = checkshort[1,:,1]
    decSca3 = checkshort[1,:,2]
    decSca12 = checkshort[1,:,11]

    ra12 =  (raSca2-raSca3)/2
    dec12 = (decSca2-decSca3)/2
    ra13 = (raSca12-raSca3)/2
    dec13 = (decSca12-decSca3)/2

    return np.array([[radecPairs[check][:,0] - ra12 - ra13, radecPairs[check][:,1] - dec12 -dec13],[radecPairs[check][:,0] - ra12 + ra13, radecPairs[check][:,1] - dec12 + dec13],[radecPairs[check][:,0] + ra12 - ra13, radecPairs[check][:,1] + dec12 - dec13]])

def slopes(corners): #return shape = (2,2560)
    return np.array([(corners[1,1,:]-corners[0,1,:])/(corners[1,0,:]-corners[0,0,:]),(corners[2,1,:]-corners[0,1,:])/(corners[2,0,:]-corners[0,0,:])])

def withinParallel(corners, x, y): #returns check indexes of all the ones that pass
    slope = slopes(corners)
    check1in = (y < slope[0,:]*(x-corners[0,0,:])+corners[0,1,:]) & (y > slope[0,:]*(x-corners[2,0,:])+corners[2,1,:])
    check2in = (y > slope[1,:]*(x-corners[0,0,:])+corners[0,1,:]) & (y < slope[1,:]*(x-corners[1,0,:])+corners[1,1,:])
    check1out = (y > slope[0,:]*(x-corners[0,0,:])+corners[0,1,:]) & (y < slope[0,:]*(x-corners[2,0,:])+corners[2,1,:])
    check2out = (y < slope[1,:]*(x-corners[0,0,:])+corners[0,1,:]) & (y > slope[1,:]*(x-corners[1,0,:])+corners[1,1,:])

    return np.nonzero((check1out & check2out) | (check1in & check2in) | (check1out & check2in) | (check2out & check1in))[0]

def openFits(group): #opens fitsfiles
    for k in group: #k stands for key (specific sca)
        if k == None:  #zip function fills with None, this prevents errors
            continue

        filt = get_filter(k)
        pointing = get_pointing(k)
        sca = get_sca(k)
        f = fits.open(f'{path}{filt}/{pointing}/Roman_TDS_simple_model_{filt}_{pointing}_{sca}.fits.gz') #———————————————————————————-edit if needed
        fitsfile = f[1]
        fitsfileData = fitsfile.data
        w = WCS(fitsfile.header)

        for radec in checkfinal[k]:
            makeStamp(radec, fitsfileData, w, k)

def makeStamp(radec, fitsfileData, w, k): #makes the stamp w/ matplotlib
    ra_dec_skycoords = SkyCoord(radec[0], radec[1], unit = "deg")
    x,y = w.world_to_pixel(ra_dec_skycoords)
    #make cutout
    try:
        cutout = Cutout2D(fitsfileData, (x,y), 31, copy=True)
    except:
        pass
    else:
        plt.imshow(cutout.data, origin = 'lower')
        plt.savefig(saveDirectory+'/'+f'{x},{y},{k}.png') #stamp saved with name x,y,k where x = x-coord in k and y = y-coord in k (needed to make CSV file)
        plt.close()

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

#SCRIPT

#run SExtractor on Search/Template image
#os.system('sex Roman_TDS_simple_model_H158_1394_12.fits -c config.sex -PARAMETERS_NAME default.param -FILTER_NAME config.sex.conv -CATALOG_NAME output.cat> /dev/null 2>&1') 

#read catalog for RADEC values
catalog = tbl.read(outputCat, hdu=2)
checkfinal = defaultdict(list) #dictionary where specific scas (key) are linked to a list of radecs (values) of each stamp
leng = len(catalog) #for debug/display purposes

#table data
t = tbl.read(scaLocTable) #read sca-radec table as an astropy table
scaAr = np.array([t['ra'],t['dec']]) #convert to numpy array, shape = (2, 57365, 18) array of centers
centers = np.reshape(scaAr, (2,numTotalPointings*18)).transpose() #pairings of radec for all scas, shape = (numPointings*18,2) each row corresponds to one sca

#get scas w/ stamps
for i in range(len(catalog)):
    RA, DEC = catalog[i]['ALPHA_J2000'], catalog[i]['DELTA_J2000']

    #first filter (< dist)
    dists = dist(centers, RA, DEC) #get distances between every sca's center's radec, shape = (1032570,)
    bounds = bound(scaAr) #find the comparison bound for the first check, shape = (1032570,)

    check = np.nonzero(dists<bounds)[0] #compare distances btwn scas to bounds, adds all that are within's index to check

    #second filter (parallel)
    check2 = check[withinParallel(getCorners(check,scaAr,centers),RA,DEC)]  #indices of the centers that match from the radecPairs/centers list

    for key in check2: #append new astro object for given sca into dictionary
        checkfinal[int(key)].append((float(RA),float(DEC)))
    print(f'{i}/{leng}') #for debug/display purposes

#create save directory
if not os.path.exists(saveDirectory): #make folder if doesn't exist
    os.makedirs(saveDirectory)

#make stamps w/ multiprocessing
numBins = -(len(checkfinal) // -numCores) #splits checkfinal dictionary into even groups to be processed in parallel
args = [iter(checkfinal)] * numBins
groups = itertools.zip_longest(*args, fillvalue=None)

if __name__ == '__main__':
    with Pool(numCores) as p:
        p.map(openFits,groups)

#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
#create csv file (overwrites)
pd.DataFrame(columns=['SNeID','filter','bogus','magnitude','name']).to_csv(CSVfile, index=False)

#table data
t = tbl.read(scaLocTable) #read sca-radec table as an astropy table

csv = open(CSVfile,'r') #opens with read-only access
csv2 = open(CSVfile,'a') #opens with append-only access

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

df = pd.read_csv(CSVfile) #read CSV into dataframe
df['ID'] = [i for i in range(len(df))] #adds an ID value for each row/stamp of CSV
df.to_csv(CSVfile, index=False) #turns it back into a CSV

print('\nfinished')