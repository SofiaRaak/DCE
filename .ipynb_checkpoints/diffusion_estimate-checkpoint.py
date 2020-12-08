import pandas as pd
import molecular_weight as mw

#Calculate radius of protein based on molecular weight
def RadiusMW(molecular_weight):
    radius = 1.5 * (molecular_weight**(1/3)) + 0.5
    return radius