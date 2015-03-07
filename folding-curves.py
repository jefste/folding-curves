import numpy as np
import scipy

import matplotlib.pyplot as plt
import pylab as pl #use this for resizing figures
from pylab import rcParams #use for resizing figures

import os,string,sys,csv #may need this to load csv files 


from scipy.optimize import fsolve
from scipy.optimize import curve_fit

#features to add:
## grab data from a csv file
    ## change test data to data from Nall dataset?
## clean up some formatting?
## add comments to some code

data_all={}

test_data_x=np.array([0,0.5,0.94,1.33,1.51,1.68,1.85,2,2.15,2.29,2.42,2.55,2.67,2.78,2.89,3,3.1,3.29,3.56,3.86,4.36])
        
test_data_y=np.array([-1573.86424,-373.92125,98.63756,1930.82677,2774.18611,5095.27712,12323.59871,19510.18524,31674.4006,53568.05441,80023.81096,115925.6316,150256.8982,177651.9872,197236.1776,208099.1024,215700.6272,240268.0321,250882.5655,262338.8597,290747.8803])


'''
definition of functions
'''
#reads file that is specified as the argument after folding-curves, otherwise uses the data found above

def getCSVfile():
    x=[]
    y=[]
    if len(sys.argv)>1:
        fileName=sys.argv[1]
        if os.path.isfile(fileName):
            with open(fileName,'rU') as csvFile:
                reader=csv.reader(csvFile,delimiter=',')
                for row in reader:
                    x.append(float(row[0]))
                    y.append(float(row[1]))
        else:
            print('File not found, script exiting')
            sys.exit()
    else:
        print('using data already loaded in routine')
        x=test_data_x
        y=test_data_y
    
    return [np.array(x),np.array(y)] 

#convert to get units of cal/kmol
R=1.98/1000
# temperature in K for 21 C
T=294.15
p_list=['C_f','m_f','C_u','m_u','m_g','D_g']

def fit_folded(x, c_f,m_f,c_u,m_u,m_g,d_g):
    #should x be changed to np.array? 
    unfold_line = c_u+m_u*x
    fold_line = c_f+m_f*x
    exp_term = np.exp(-(d_g-m_g*x)/(R*T))
    function_fit = (fold_line +unfold_line*exp_term)/(1+exp_term)
    return function_fit

def fold_line(x, c_f,m_f,c_u,m_u,m_g,d_g):
    return c_f+m_f*x

def unfold_line(x, c_f,m_f,c_u,m_u,m_g,d_g):
    return c_u+m_u*x

def fold_unfold_fraction_func(x, c_m, c_f,m_f,c_u,m_u,m_g,d_g):
    exp_term=np.exp(-m_g*(c_m-x)/(R*T))
    return exp_term/(1+exp_term)

def fold_unfold_fraction_data(denaturant, observed_y, c_f,m_f,c_u,m_u,m_g,d_g):
    denaturant=np.array(denaturant)
    observed_y=np.array(observed_y)
    k_r_t=(observed_y-(c_u+m_u*denaturant))
    k_r_b=((c_f+m_f*denaturant)-observed_y)
    k_r=k_r_t/k_r_b
    return 1-k_r/(1+k_r)

# creates initial guesses for a given data set (x,y) for the slopes and intercepts
def initial_parameters(x, y):
    # use x and y to estimate the intercept and slope for the folded line
    c_f=y[0]
    m_f=(y[2]-y[0])/(x[2]-x[0])
    # use the last 2 coordinates to estimate the slope of the unfolded line
    m_u=(y[-1]-y[-3])/(x[-1]-x[-3])
    c_u=x[-1]*m_u-y[-1]
    return [c_f,m_f,c_u,m_u]

'''
Reads data from specified csv file
'''
data_all['x_raw'],data_all['y_raw']=getCSVfile()   

'''
fitting of RAW data to a folded/unfolded plot
'''

#make an array of 1000 points
data_all['xfit']=np.linspace(data_all['x_raw'][0],data_all['x_raw'][-1],1000)


#fits data to function fit_folded, returns the parameters in popt_fu
## generate the initial guess from the initial_parameters function
print('initial paramters')
print(initial_parameters(test_data_x,test_data_y)) 
popt_fu,pcov_fu = curve_fit(fit_folded, data_all['x_raw'], data_all['y_raw'],p0=initial_parameters(data_all['x_raw'],data_all['y_raw'])+[3,5])

#generates points for y for the fit parameters
data_all['y_fit_raw']=fit_folded(data_all['xfit'],*popt_fu)

'''
Converted data
'''
#converts the raw data to percentage data
## currently having an issue on this
data_all['y_data_percent_unfold']=fold_unfold_fraction_data(data_all['x_raw'],data_all['y_raw'],*popt_fu)

#fits data to fold_unfold_fraction_func
popt_fup,pcov_fup = curve_fit(fold_unfold_fraction_func, data_all['x_raw'], data_all['y_data_percent_unfold'])

#generates points for y from the fit parameters for the percent unfolded curve
data_all['y_fit_percent_unfold']=fold_unfold_fraction_func(data_all['xfit'],*popt_fup)



'''
first plot
'''
'''
COMMENTING OUT SECOND PLOT, ONLY USING 3rd plot
#add gridlines to plot
plt.grid(b=True, which='both',color='0.5',linestyle='-',alpha=.3)

#plot data
plt.plot(test_data_x,test_data_y,'bo',label='data')
#plot fit
plt.plot(xfit,fit_y,'r-',label='fit')


#adds a green dashed line for the unfolded baseline
plt.plot(xfit,unfold_line(xfit,*popt_fu),'g--',label='unfolded baseline')

#adds a magenta dashed line for the folded baseline
plt.plot(xfit,fold_line(xfit,*popt_fu),'m--',label='folded baseline')

#adds legend to the location upper left
plt.legend(loc='upper left')


#add gridlines to plot
plt.grid(b=True, which='both',color='0.5',linestyle='-',alpha=.3)



plt.show()
'''


'''
write out parameters to be printed to table
'''
#adds C_m to the table
cell_text=[['C_m',np.round_(popt_fup[0],2)]
        ]
#adds all parameters returned from fit of data
for i in range(len(p_list)):
    cell_text.append([p_list[i], np.round_(popt_fu[i],2) ])


'''
second plot
'''

'''
COMMENTING OUT SECOND PLOT, ONLY USING 3rd plot
#plot data
plt.plot(test_data_x,test_data_y_percent_fu,'bo',label='data')
#plot fit
plt.plot(xfit,fit_y_percent_fu,'r-',label='fit')

#adds legend to the location upper left
plt.legend(loc='upper left')

#move plot to make room for table
plt.subplots_adjust(right=0.5)

#add labels
plt.ylabel('Percent Unfolded')
plt.xlabel('Denaturant (GdnHCl in M)')
plt.title('Percent Unfolded vs Denaturant')

#add table to plot on the right hand side

table_plot=plt.table(cellText=cell_text,loc='right')
#increase the font size of the table
table_plot.set_fontsize(24)
#increase the size of the plot (1 for horizontal, 4 for vertical)
table_plot.scale(1,4)


#add gridlines to plot
plt.grid(b=True, which='both',color='0.5',linestyle='-',alpha=.3)

plt.show()
'''



'''
Third plot, One plot containing all info
'''


#change the size of the figure, and the size of the font on the figure.
rcParams['figure.figsize'] = 20, 7
rcParams['font.size'] = 14
#change font size of legend, as it was clipping the fit data
rcParams['legend.fontsize'] = 12


fig,((ax1,ax2,ax3))= plt.subplots(nrows=1,ncols=3)

#plot data
ax1.plot(data_all['x_raw'],data_all['y_raw'],'bo',label='data')
#plot fit
ax1.plot(data_all['xfit'],data_all['y_fit_raw'],'r-',label='fit')


#adds a green dashed line for the unfolded baseline
ax1.plot(data_all['xfit'],unfold_line(data_all['xfit'],*popt_fu),'g--',label='unfolded baseline')
#adds a magenta dashed line for the folded baseline
ax1.plot(data_all['xfit'],fold_line(data_all['xfit'],*popt_fu),'m--',label='folded baseline')

#adds legend to the location upper left
ax1.legend(loc='upper left')

ax1.axes.set_ylabel('Observed Parameter')
ax1.axes.set_xlabel('Denaturant (GdnHCl in M)')
ax1.axes.set_title('Observed Parameter vs Denaturant')

#add gridlines to plot
ax1.grid(b=True, which='both',color='0.5',linestyle='-',alpha=.3)

#second plot with lables and legend
ax2.plot(data_all['x_raw'],data_all['y_data_percent_unfold'],'bo',label='data')
ax2.plot(data_all['xfit'],data_all['y_fit_percent_unfold'],'r-',label='fit')
ax2.legend(loc='upper left')
ax2.axes.set_ylabel('Percent Unfolded')
ax2.axes.set_xlabel('Denaturant (GdnHCl in M)')
ax2.axes.set_title('Percent Unfolded vs Denaturant')

#add gridlines to plot
ax2.grid(b=True, which='both',color='0.5',linestyle='-',alpha=.3)

#third 'plot', is just the table, but hides the x and y labels since they aren't needed
ax3.table(cellText=cell_text,loc='center').scale(1,4)
ax3.axes.get_yaxis().set_visible(False)
ax3.axes.get_xaxis().set_visible(False)
ax3.axes.set_title('Fitting parameters')


# gives too much padding, maybe some way to tweak this?
#plt.tight_layout()
plt.show()





'''
output table to console, use this for debugging
'''

for row in cell_text:
    print row


'''
save file as csv

def writeToCSV():
    raw_input('Save fit and parameters to csv files? (y/n)')
    
    #checks to see if the path and the file name exists. 
    #Is this redunant? Maybe only need to check for file exisiting?
    #Note that when checking exists or isfile, the function does not appear to be case dependant
        if os.path.exists(pathnamePDB) and os.path.isfile(pathnamePDB):



saves data to 3 separate files: paramters, fit curves, and data converted to perecnet unfolded

'''
def save_to_CSV(sample_name):
    answer=raw_input('Save fit and parameters to csv files (y/n)? ')
    if answer=='y':
        with open('fit_CSVs/'+sample_name+'_parameters.csv', 'wb') as csvfile:
            datawriter = csv.writer(csvfile, delimiter=',')
            for row in cell_text:
                datawriter.writerow(row)

        with open('fit_CSVs/'+sample_name+'_fits.csv', 'wb') as csvfile:
            datawriter = csv.writer(csvfile, delimiter=',')
            for i in range(len(xfit)):
                datawriter.writerow([xfit[i],fit_y[i],fit_y_percent_fu[i]])

        with open('fit_CSVs/'+sample_name+'_data_percent_unfolded.csv', 'wb') as csvfile:
            datawriter = csv.writer(csvfile, delimiter=',')
            for i in range(len(test_data_x)):
                datawriter.writerow([test_data_x[i],test_data_y_percent_fu[i],test_data_y[i]])


if len(sys.argv)>1:
    #split ext then calling for [0] index gives the name of the specified file without the file extension
    save_to_CSV(os.path.splitext(sys.argv[1])[0])
else:
    save_to_CSV('test_data')
