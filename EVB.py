#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 17:15:44 2019

@author: x_chxia
"""

import glob 
import re
import numpy as np
import matplotlib.pyplot as plty
import math
from scipy.signal import savgol_filter
# =============================================================================
# These values has to be altered mannually to calibrate EVB
# =============================================================================
n=237#  no.bins   
alpha=-693#
H12=9.5#
print('alpha=', alpha ,'H12=', H12)

# =============================================================================
# input values read from respected files
# =============================================================================
l=[]
with open("lambda.xvg","r") as file:
    for line in file:
           l=map(float,line.split())
           
fep_energy=[]
with open("fep.xvg","r") as file:
    for line in file:
            row=line.split()
            fep_energy.append(row[-1])
fep_energy=list(map(float,fep_energy))

solv=[]          
with open("solv.xvg", "r") as file:
    for line in file:
        row=line.split()
        solv.append(row[-1])
solv=list(map(float,solv))

# =============================================================================
# Combine alpha(gas shift) and solvation energy to energy     
# =============================================================================
fep_i=[]
fep_i.append(fep_energy[0])
for i in range(1,len(fep_energy)):
    fep_i.append(fep_energy[i]+fep_i[i-1]+(l[i]-l[i-1])*alpha)
    
for i in range(0,len(fep_i)):
    fep_i[i]=fep_i[i]+solv[i]-solv[0]

# =============================================================================
#decide the x axis(E1-E2) range
#in order to ensure all the values are covered, E1 reads from state1 and E2 state2
# =============================================================================
E2_lambda0=[]
E1_lambda1=[]
with open ("0.txt","r") as file:
       for line in file:
               row=line.split()
               E2_lambda0.append(row[-1])
with open ("67.txt","r") as file:
       for line in file:
               row=line.split()               
               E1_lambda1.append(row[-2])

E2_lambda0=list(map(float,E2_lambda0))
E1_lambda1=list(map(float,E1_lambda1))

lower_limit=0-max(E2_lambda0)
upper_limit=max(E1_lambda1)+abs(alpha)
bins = np.arange (lower_limit,upper_limit,(upper_limit-lower_limit)/float(n))
bins=np.append(bins,(upper_limit+2000))

# =============================================================================
# set global variables
# =============================================================================
test=[]
l_total=[]  #array for lambda values, in order of how files are read 
E1=[]
E2=[]
e1=[]
e2=[]
diff=[] #difference between E1 and E2, used to decide where to bin the data point
dG=[]
#create a 2D list to put binned data
hist1=[[] for i in range(n+1)]
#lambda_array=[[] for i in range(len(l))] 
hist_final1=[[] for i in range(n+1)]
#FEP=[[] for i in range(36)]
float_E1=[]
float_E2=[]
diff1=[]

# =============================================================================
# loop through the lambda.txt files and calculate free energy 
# =============================================================================
for fil in glob.glob("*.txt"):
    l_value=re.search(r'\d+',fil).group(0)
    l_value=int(l_value)
    l_total.append(l_value)
    
    with open (fil,"r") as file:
            for line in file:
                row=line.split()
                E1.append(row[-2])
                E2.append(row[-1])
            E1=list(map(float,E1))
            E2=list(map(float,E2))

#rearrange the values to exclude extreme points:
            for i in range (len(E1)):

#exclude the extreme points which are so large that can mess up the whole distribution
#Need to check the $lambda.xvg files to decide how to set these limit
#working on a more universal solution at the moment
               if E1[i] <= 0 and E2[i]>= 1 :
                     diff.append(E1[i]-E2[i]-alpha)
                     float_E1.append(E1[i])
                     float_E2.append(E2[i]+alpha)           
            #distribute the points in bins according to diff(E1-E2)
            bin_number=np.digitize(diff,bins,right=False)

            for i in range(len(bin_number)):
                bn=bin_number[i] #a local variable  
                Eg=0.5*(float_E2[i]+float_E1[i])-0.5*math.sqrt((float_E2[i]-float_E1[i])**2+4*H12*H12)
                dG.append(Eg)               
                vi=(1-l[l_value])*(float_E1[i])+(float_E2[i])*(l[l_value])

                hist1[bn-1].append(np.exp(-(Eg-vi)/2.479)) #2.479 being RT in kj/mol

            for i in range(n):
                if len(hist1[i])>0:
                    average=sum(hist1[i])/len(hist1[i])
                    for j in range(len(hist1[i])):
                        if l[l_value]==0:
                            hist_final1[i].append(0-2.479*np.log(average))
                        else:
                            hist_final1[i].append(fep_i[l_value-1]-2.479*np.log(average))
                       
                        hist1[i][:]=[]

#Empty the variables for the next iteration
            diff[:]=[]
            diff1[:]=[]
            E1[:]=[]
            E2[:]=[]
            float_E1[:]=[]
            float_E2[:]=[]
            dG[:]=[]

# =============================================================================
# Do averages for the sampling and hence calculate the final profile
# =============================================================================
free_energy=[]
x=[]
for i in range(n):
    if len(hist_final1[i])>5:
        aa=0.239*sum(hist_final1[i])/len(hist_final1[i])
        b=0.239*bins[i]
        free_energy.append(aa) 
        x.append(b)        

# =============================================================================
# Print out free energy and activation energy        
# =============================================================================
free_max=max(free_energy)
free_max_index=free_energy.index(free_max)
free_min=min(free_energy)
reac_free=[]
for i in range(free_max_index):
    reac_free.append(free_energy[i])  
reac_min=min(reac_free)
print("Activation free energy", "{0:.2f}".format(free_max-free_min))
print("Reaction free energy","{0:.2f}".format(reac_min-free_min))    



# =============================================================================
# Smooth the data
# =============================================================================
interv=int(len(x)/3)
if (interv %2 ) == 0:
    interv=interv+1
sav = savgol_filter(free_energy,interv,5)
    
  
# =============================================================================
# plot
# =============================================================================
plty.plot(x,sav,'.-k',label='Co(TPP) in water')
plty.xlabel("Reaction Coordinate / e1-e2")
plty.ylabel("Free Energy (kcal/mol)")
plty.legend(loc='upper left')
ax= plty.gca()
ax.invert_xaxis()
# =============================================================================
# Uncomment the section to add a lebel, values have to be typed in mannually
# =============================================================================
#plty.text(124,1.5,' $\Delta$ G  = 4.2 kcal/mol \n $\Delta$ $\mathregular{G^\ddag}$ = 9.8 kcal/mol ',
#          bbox=dict(boxstyle="round", fc="none", ec="gray"))

plty.show
plty.savefig('fig.png',dpi=1000)
