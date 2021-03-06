import pandas as pd
import molecular_weight as mw
import math

#Bolzmann constant
kB = 1.38064852*(10**-23)

#Calculate radius of protein based on molecular weight, based on analysis in jupyter notebook
def RadiusMW(molecular_weight):
    radius = 1.49 * (molecular_weight**(1/3)) + 0.44
    return radius

#Calculate range of possible diffusion coefficients
def DiffusionCoefficient(my_email, acc_numbers, download=True, file=None):
    
    if download == True:
        data = mw.MolWeights(my_email, acc_numbers, download=True, file=None)
    
    elif download == False:
        data = mw.MolWeights(my_email, acc_numbers, download=False, file=file)
    
    data.insert(2, 'R, nm', RadiusMW(data['kDa']))
    
    D_max = pd.Series(((kB*310.15)/(6*math.pi*0.002*data['R, nm']*(10**-9)))*10**12)
    D_min = pd.Series(((kB*310.15)/(6*math.pi*0.05*data['R, nm']*(10**-9)))*10**12)
    
    data.insert(3,'D max, um^2/s', D_max)
    data.insert(4, 'D min, um^2/s', D_min)

    return data

