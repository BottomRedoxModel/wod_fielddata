# wod_fielddata
Python 2.7 script, preparing field data for using as boundary conditions or for field data relaxation.  
Currently the script produces monthly averaged data along the needed profile. 
Data is interpolated for standard levels and standard distances from the beginning point, in script we use linear interpolation. 
Obtained values are written into text files. There are 2d arrays, x axis is distances, y axis is depths. Each file is for 1 month and 1 variable.
### The main algorithm of preparation data is:
1.	To download data from region of interest from World Ocean Database  https://www.nodc.noaa.gov/OC5/WOD/datageo.html 
2.	Create the Ocean Data View collection, choose more precise region of interest, and produce netCDF file with these data. 
3.	Run the script using Python 
4.	Read the produced text file in BROM

Later in BROM as relaxation data may be used as single columns with data or the whole profile. 
