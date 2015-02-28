import numpy as np
import scipy

import matplotlib.pyplot as plt

## import os,string,sys #may need this to load csv files 


from scipy.optimize import fsolve
from scipy.optimize import curve_fit

#features to add:
## grab data from a csv file
    ## change test data to data from Nall dataset?
## clean up some formatting?
## add comments to some code

#changed from a tuple to list for the test data, also change to numpy array
test_data_x=np.array([0,0.5,0.94,1.33,1.51,1.68,1.85,2,2.15,2.29,2.42,2.55,2.67,2.78,2.89,3,3.1,3.29,3.56,3.86,4.36])
        
test_data_y=np.array([-1573.86424,-373.92125,98.63756,1930.82677,2774.18611,5095.27712,12323.59871,19510.18524,31674.4006,53568.05441,80023.81096,115925.6316,150256.8982,177651.9872,197236.1776,208099.1024,215700.6272,240268.0321,250882.5655,262338.8597,290747.8803])


#normalize the data to try to have better fitting?
#test_data_y=[item/max(test_data_y) for item in test_data_y]


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
    function_fit = (unfold_line +fold_line*exp_term)/(1+exp_term)
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
    k_r_t=(observed_y-(c_f+m_f*denaturant))
    k_r_b=((c_u+m_u*denaturant)-observed_y)
    k_r=k_r_t/k_r_b
    return k_r/(1+k_r)

# creates initial guesses for a given data set (x,y) for the slopes and intercepts
def initial_parameters(x, y):
    # use x and y to estimate the intercept and slope for the folded line
    c_f=y[0]
    m_f=(y[1]-y[0])/(x[1]-x[0])
    # use the last 2 coordinates to estimate the slope of the unfolded line
    m_u=(y[-1]-y[-2])/(x[-1]-x[-2])
    c_u=x[-1]*m_u-y[-1]
    return [c_f,m_f,c_u,m_u]


#need to find method to fit and return all parameters
## use those parameters to generate a fit plot, fold baseline and unfold baseline
'''
fitting of RAW data to a folded/unfolded plot

'''

#make an array of 1000 points
xfit=np.linspace(test_data_x[0],test_data_x[-1],1000)


#fits data to function fit_folded, returns the parameters in popt_fu
## generate the initial guess from the initial_parameters function 
popt_fu,pcov_fu = curve_fit(fit_folded, test_data_x, test_data_y,p0=initial_parameters(test_data_x,test_data_y)+[1,1])

#generates points for y for the fit parameters
fit_y=fit_folded(xfit,*popt_fu)

'''
Converted data
'''
#converts the raw data to percentage data
## currently having an issue on this
test_data_y_percent_fu=fold_unfold_fraction_data(test_data_x,test_data_y,*popt_fu)

#fits data to fold_unfold_fraction_func
popt_fup,pcov_fup = curve_fit(fold_unfold_fraction_func, test_data_x, test_data_y_percent_fu)

#generates points for y from the fit parameters for the percent unfolded curve
fit_y_percent_fu=fold_unfold_fraction_func(xfit,*popt_fup)



'''
first plot
'''
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
plt.show()



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

plt.table(cellText=cell_text,loc='right')

plt.show()




