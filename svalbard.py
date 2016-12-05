'''
Created on 28. nov. 2016

@author: ELP
'''
import pdb
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import numpy.ma as ma
import numpy as np
from scipy import interpolate
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import pandas as pd

ym = [0,2,5,10,15,20,30,40,50,100,200,300,400,500,1000,1500,
        2000,2500,3000,3500]#range(0,2500,50)#
xm = range(0,2000,50)#[0,100,200,300,400,500,700,750,800,900,1000]
X,Y = np.meshgrid(xm,ym)
X_f = X.flatten()
Y_f = Y.flatten()

def read_netcdf(boundary):

    ncfile = boundary 
    fh = Dataset(ncfile, mode='r')
    date_time = fh.variables['date_time'][:]
    lat =fh.variables['latitude'][:]
    long =fh.variables['longitude'][:]

    distancemasked = fh.variables['var16'][:][:]
    distance =  distancemasked.filled(fill_value= np.nan)
    #same array, but non masked, nan are just the values
    depthmasked = fh.variables['var1'][:][:]
    depth = depthmasked.filled(fill_value= np.nan)
    
    tempmasked = fh.variables['var2'][:][:]
    temp = tempmasked.filled(fill_value= np.nan)

    salmasked = fh.variables['var3'][:][:]
    sal = salmasked.filled(fill_value= np.nan)

    o2masked = fh.variables['var4'][:][:]
    o2 = o2masked.filled(fill_value= np.nan)

    po4masked = fh.variables['var5'][:][:]
    po4 = po4masked.filled(fill_value= np.nan)

    simasked = fh.variables['var6'][:][:]
    si = simasked.filled(fill_value= np.nan)

    no3masked = fh.variables['var7'][:][:]
    no3 = no3masked.filled(fill_value= np.nan)

    no2masked = fh.variables['var8'][:][:]
    no2 = no2masked.filled(fill_value= np.nan)

    pHmasked = fh.variables['var9'][:][:] #almost no data
    pH = pHmasked.filled(fill_value= np.nan)

    chlmasked = fh.variables['var10'][:][:] #almost no data
    chl = chlmasked.filled(fill_value= np.nan)
    
    alkmasked = fh.variables['var12'][:][:] #almost no data
    alk = alkmasked.filled(fill_value= np.nan)   
    fh.close()


    dates = num2date(date_time[:],units='days since 1990-01-01',
                 calendar='proleptic_gregorian' )
    return distance,depth,temp,sal,o2,po4,si,no3,no2,pH,chl,dates,alk
def to_date_tuple(file):
    i = 0  
    time = [] 
    for n in file:
        i=i+1
        t = n.timetuple()
        time.append(t) #create tuple to select dat by month
    return time
     
def select_by_month(month,file,num, name,var,distance_var,depth_var,boundary):
    j=0 
    depths =  []
    distances = []
    var_ = []     

    time = to_date_tuple(file) 
          
    for n in time:
        if n[1] == month: #n[1] - place of month in a time tuple                
            depths.append(depth_var[j])
            distances.append(distance_var[j]) 
            var_.append(var[j]) 
#            years.append(n[0])  
#            months.append(n[1])   
#            days.append(n[2])    
            j =j +1
        else: 
            j =j +1  
              
    return depths,distances,var_

def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]

def y_interp(month,file,num, name,var,distance_var,depth_var,boundary):
    
    var_y = select_by_month(month,file,num, name,var,distance_var,depth_var,boundary)

    depth_3 = np.array(var_y[0])#.flatten() 
    dist_3 = np.array(var_y[1])#.flatten()  
    var_3 = np.array(var_y[2])#.flatten() 

    arr_y = []
    arr_x = []
    arr_z = []

    for s in range(0,len(depth_3)):

        y = depth_3[s]
        z = var_3[s]
        m = dist_3[s]  
    
        if np.isnan(m).all()== True or np.isnan(z).all()== True: 
        #check if the column contains only nan valuers
#            pdb.set_trace()
            continue
        else:       
            nans, x = nan_helper(y)
            y[nans]= np.interp(x(nans), x(~nans), y[~nans])

            nans, x = nan_helper(z)
            z[nans]= np.interp(x(nans), x(~nans), z[~nans])
            zz = np.interp(ym,y,z)
            arr_z.append(zz)
  
            nans, x = nan_helper(m)
            m[nans]= np.interp(x(nans), x(~nans), m[~nans])
            mm = np.interp(ym,y,m)
            arr_x.append(mm)

    var_int = np.array(arr_z)
    var_int_t = var_int.T
    dist_int = np.array(arr_x)
    dist_int_t = dist_int.T
    var_x_y = []    


    for n in range(0,len(dist_int_t)):
        m = interpolate.interp1d(dist_int_t[n], var_int_t[n],
            assume_sorted=False,bounds_error=False, fill_value=np.nan)
        var_x_y_new = m(xm)
        var_x_y.append(var_x_y_new)    

#    np.savetxt('{}_month{}_{}.txt'.format(boundary,month,name), (var_x_y), delimiter=" ")     
               
    plt.subplot(7,7,num)
    plt.ylim(4000,0)
    plt.xlim(0,2000)
    origin = 'lower'

    CS = plt.contourf(X,Y,  var_x_y, 10,
                      #[-1, -0.1, 0, 0.1],
                  #alpha=0.5,
                      cmap=plt.cm.bone,
                      origin=origin)
    plt.colorbar()

    #CS = plt.contour(X_f,Y_f,zi,15,linewidths=0.5,colors='k')
    #plt.scatter(arr_x[0], ym,  s=80, c=arr_z[0], marker=">")
    plt.scatter(X, Y,  s=80, c=var_x_y, marker=">")
    plt.title('{}month{}_{}.txt'.format(boundary,month,name))

def call_boundary(boundary):
    w = read_netcdf(boundary)
    distance = w[0]
    depth = w[1]
    temp = w[2] 
    sal = w[3]
    o2 = w[4]
    po4 = w[5]
    si = w[6]
    no3 = w[7]
    no2 = w[8]
    pH = w[9]
    chl = w[10]
    dates = w[11]
    alk = w[12]
    i = 1
    if boundary == 'border1.nc':
        si_list = [2,3,5,6,7,8,9,11]#[3,6,7,8,9,11] 
        po4_list = [2,3,5,6,7,8,9,11]#[2,3,7,8,11]
        no3_list = [2,3,5,6,7,8,9,11]
        o2_list = [2,3,5,6,7,8,9,11]
        chl_list = [5,6,7,8]  #[6,7,8]
        alk_list = [5,6,7,8] 
        pH_list = [2,5,8]       
        
    elif boundary == 'border2.nc':
        si_list = [2,3,5,6,7,8,9,11] 
        po4_list = [2,3,5,6,7,8,9,11]
        no3_list = [2,3,5,6,7,8,9,11]
        o2_list = [4,6,7,9,11] #,6,7,9,112,]
        chl_list = [5,6,7]   
        alk_list = [6,7,9,10,12] 
        pH_list = [9,12]                     
    elif boundary == 'border3.nc':
#        print dates
        si_list = [6] #,8,9
        po4_list = [6]
        no3_list = [6]
        o2_list = [6,9]
        chl_list = [6]
        alk_list = [9] 
        pH_list = [9]
        
    elif boundary == 'border4_1.nc':
        si_list = [8] #,8,9
        po4_list = [8]
        no3_list = [8]
        o2_list = [8]
        chl_list = [8]  
        alk_list = [5,8] 
        pH_list = []         
                      
    for n in si_list:
        y_interp(n,dates,i,'Si',si,distance,depth,boundary)
        i = i + 1 
    for n in po4_list:
        y_interp(n,dates,i,'po4',po4,distance,depth,boundary)   
        i = i + 1     
    for n in no3_list:
        y_interp(n,dates,i,'no3',no3,distance,depth,boundary)   
        i = i + 1      
    for n in o2_list:
        y_interp(n,dates,i,'O2',o2,distance,depth,boundary)  
        i = i+1 
    for n in chl_list:
        y_interp(n,dates,i,'Chlorophyll',chl,distance,
                 depth,boundary)  
        i = i+1            
    for n in alk_list:
        y_interp(n,dates,i,'Alk',alk,distance,
                 depth,boundary)              
        i = i+1  
        
    for n in pH_list:
        y_interp(n,dates,i,'pH',pH,distance,
                 depth,boundary)              
        i = i+1         
        
        
        
#call_boundary('border1.nc') 
#plt.show() 
#call_boundary('border2.nc')
#plt.show()  
#call_boundary('border4_1.nc')   
#plt.show() 
call_boundary('border3.nc') 
plt.show()    
         

    