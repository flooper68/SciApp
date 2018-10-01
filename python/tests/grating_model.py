import os

import matplotlib.pyplot as plt
import pandas

#%%

pathToData = 'file://localhost' + os.getcwd() + '/python/tests/data/Grating_840nm.csv'
n_rows = 520

#angle 6 deg,
data1 = pandas.read_csv(pathToData, header=1,usecols=[522, 523], nrows=n_rows)

data2 = pandas.read_csv(pathToData, header=1,usecols=[528, 529], nrows=n_rows)

plt.plot(data1.iloc[:,0], data1.iloc[:,1])
plt.plot(data2.iloc[:,0], data2.iloc[:,1])
plt.show()


#%%

col3 = 196
col4 = 200
data3 = pandas.read_csv(pathToData, header=1,usecols=[col3, col3+1], nrows=n_rows)
data4 = pandas.read_csv(pathToData, header=1,usecols=[col4, col4+1], nrows=n_rows)

print(data3.head())
plt.plot(data3.iloc[:,0], data3.iloc[:,1])
plt.plot(data4.iloc[:,0], data4.iloc[:,1])
plt.show()

#%%

import pymongo

client = pymongo.MongoClient()

db = client.MyDb

collection = db.Data

