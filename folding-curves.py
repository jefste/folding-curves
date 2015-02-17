import numpy as np
import scipy

import matplotlib.pyplot as plt



from scipy.optimize import fsolve
from scipy.optimize import curve_fit

def fit_folded(x, c_f,m_f,c_u,m_u,m_g,d_g):
#convert to get units of cal/kmol
    R=1.98/1000
# temperature in K for 21 C
    T=294.15
#    unfold_line = c_u+m_u*x
#    fold_line = c_f+m_f*x
#    exp_term = exp(-(d_g-m_g*x)/(R*T))
#    function_fit = (unfold_line +fold_line*exp_term)/(1+exp_term)
    function_fit = (c_u+m_u*x+(c_f+m_f*x)*np.exp(-(d_g-m_g*x)/(R*T)))/(1+np.exp(-(d_g-m_g*x)/(R*T)))
    return function_fit

def fold_line(x, c_f,m_f,c_u,m_u,m_g,d_g):
    return c_f+m_f*x

def unfold_line(x, c_f,m_f,c_u,m_u,m_g,d_g):
    return c_u+m_u*x




test_data_x=(0,0.5,0.94,1.33,1.51,1.68,1.85,2,2.15,2.29,2.42,2.55,2.67,2.78,2.89,3,3.1,3.29,3.56,3.86,4.36)
        
test_data_y=(-1573.86424,-373.92125,98.63756,1930.82677,2774.18611,5095.27712,12323.59871,19510.18524,31674.4006,53568.05441,80023.81096,115925.6316,150256.8982,177651.9872,197236.1776,208099.1024,215700.6272,240268.0321,250882.5655,262338.8597,290747.8803)


#need to find method to fit and return all parameters
## use those parameters to generate a fit plot, fold baseline and unfold baseline


#make an array of 100 points
xfit=np.linspace(test_data_x[0],test_data_x[-1],1000)




popt_fu,pcov_fu = curve_fit(fit_folded, test_data_x, test_data_y)


fit_y=fit_folded(xfit,*popt_fu)
#fit_y=test_f(xfit,*popt)

#plot data
plt.plot(test_data_x,test_data_y,'bo',label='data')
#plot fit
plt.plot(xfit,fit_y,'r-',label='fit')



plt.plot(xfit,unfold_line(xfit,*popt_fu),'g--',label='unfolded baseline')
plt.plot(xfit,fold_line(xfit,*popt_fu),'m--',label='folded baseline')
plt.legend(loc='upper left')
plt.show()



#need to put fold baseline on plot
#need to put unfold baseline on plot
#need to put fit on plot


