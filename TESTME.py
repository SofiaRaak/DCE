import diffusion_estimate as de
import sys

my_email = sys.argv[1]

#Dictionary with protein name, diffusion coefficient (data from Tyn and Gusek, 1990)
experimental = {'Bovine serum albumin':5.98, 'Chicken lysozyme':11.4,
                'Yeast phosphoglycerate kinase':6.38}

#lists of example proteins for estimating diffusion coefficients
trial_list = ['CAA76847.1', 'ACL81750.1',  'CAA42329.2']

test = de.DiffusionCoefficient(my_email, trial_list)

BSA = test.iloc[0]
CL = test.iloc[1]
YPGK = test.iloc[2]

print()
print(f"The estimated range of diffusion coefficients for bovine serum albumin is {round(BSA['D min, um^2/s'],2)}-{round(BSA['D max, um^2/s'],2)} um^2/s")
print(f"The experimentally determined diffusion coefficient for BSA is {experimental['Bovine serum albumin']} um^2/s")
print()
print(f"The estimated range of diffusion coefficients for chicken lysozyme is {round(CL['D min, um^2/s'],2)}-{round(CL['D max, um^2/s'],2)} um^2/s")
print(f"The experimentally determined diffusion coefficient for chicken lysozyme is {experimental['Chicken lysozyme']} um^2/s")
print()
print(f"The estimated range of diffusion coefficients for yeast phosphoglycerate kinase is {round(YPGK['D min, um^2/s'],2)}-{round(YPGK['D max, um^2/s'],2)} um^2/s")
print(f"The experimentally determined diffusion coefficient for YPGK is {experimental['Yeast phosphoglycerate kinase']} um^2/s")
print()

