# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, sys
lib_path = os.path.abspath('PyMods/')
sys.path.append(lib_path)

from INSTRUMENTALS import *
from COSMO import *
from UTILITIES import *
from VARIABLES import *

# "Black boxes" #
from cosmocalc import cosmocalc
from math import *
from decimal import *
import numpy as npy
import re

from matplotlib import rc
from matplotlib import rcParams
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot

plt.rcParams.update({'text.usetex' : 'true',
                     'backend': 'pdf',
                     'font.family':'sans-serif',
                     'font.size': 25,
                     'font.weight':'normal',
                     'legend.fontsize': legFontsz,
                     'legend.frameon': True,
                     'legend.loc': 'best',
                     'lines.markersize': 3,
                     'lines.linewidth': 2.5,
                     'axes.linewidth': .5,
                     'axes.labelweight': 'normal',
                     'axes.labelsize':25,
                     'axes.color_cycle': ['#E41A1C', '#377EB8', '#4DAF4A', '#FF7F00']})
########################
# ### User variables ###
########################

# #### Note: Make sure you have files containing the smoothing scale vs skew data in the working directory called:
# 
# 21CMFASTlist 
# GlobMHRlist 
# MHRlist 
# InvMHRlist 
# 
# #### Important there is no white space lines at end of file ####

#---------------- Do you want to normalise the moments? -------------
norm = 1; 
# 0 = no normalisation; 
# 1 = normalise the nth moment with sigma^n;  
# 2 = normalise with sigma^2; 

xRange=[0.003,1.0]
if (norm==1):
    sRange=[-2.5,11.95] #Range for the skew plot
if (norm==2):
    sRange=[-5,20.5] #Range for the skew plot
vRange=[0.1,85] # Range for the variance plot. (8=10, 17 R=30

Revo=10 #smoothing scale you want to plot the evolution of
z2plot=8.75 # redshift to plot for the Rsmooth vs. skew/variance plots

############## Instrumental parameters ###############
inst=3; #0=MWA, 1=Lofar, 2=SKA, 3=HERA (Rsmooth vs. skew/variance plot)

# --- 21cmfast et al. Box details ---
L=300; #Size of one side of box being analysed in Mpc
Res=600; # in pixels per side
Zeta=16; #value of zeta used for simulation

#model legend and corresponding folder 
key=['global-in-out','local-in-out', 'global-out-in', 'local-out-in']
model=['21CMFAST','InvMHR','GlobMHR','MHR']
inst_str=['MWA','LOFAR','HERA', 'SKA']
inst_err=[1,1,1,1]

# ## IN/OUT OPTIONALS ##

#location you require your plots to be saved:
saveDIR = 'Plots/' 

# ##Initialisation##


fig3 = plt.figure(num=None, figsize=figSize, dpi=80, 
                  facecolor='w', edgecolor='k');
ax3 = host_subplot(111, axisbg=BGcol); 

fig4 = plt.figure(num=None, figsize=figSize, dpi=80, 
                  facecolor='w', edgecolor='k');
ax4 = host_subplot(111, axisbg=BGcol); 

fig5 = plt.figure(num=None, figsize=figSize, dpi=80, 
                  facecolor='w', edgecolor='k');
ax5 = host_subplot(111, axisbg=BGcol); 

fig6 = plt.figure(num=None, figsize=figSize, dpi=80, 
                  facecolor='w', edgecolor='k');
ax6 = host_subplot(111, axisbg=BGcol); 

Tint = instParams(inst,'Tint')
Atot = instParams(inst,'Nstat')*instParams(inst,'Aeff')

## THE GUTS OF THE CODE:

i=0
while i<len(model):
    #arrays to store the values at each redshift for the choosen smoothing scale
    ionfstore=[];
    sig2store=[];
    skewstore=[];
    S3errstore=[];
    S2errstore=[];
    
    File=model[i]+'list'
    
    with open(File) as a_file:  
        
        for a_line in a_file:
            noise=[]
            S2_err=[]
            S3_err=[]
            filestr=a_line.rstrip()
            #extract the variables we need from the filename
            try:
                found = re.search('_z(.+?)_nf', filestr).group(1)
            except AttributeError:
                # _z, _use not found in the original string
                found = '' # apply your error handling
            z1=float(found)
            try:
                found = re.search('_nf(.+)', filestr).group(1)
            except AttributeError:
                # _nf, _z not found in the original string
                found = '' # apply your error handling
            nf1=float(found)
            
            (Rsmth, kurt, skew, sig2) = npy.loadtxt(filestr, delimiter='\t', usecols=(0,3,5,7), unpack=True);
   
            x = Cosmology(z1, H_0, Om_m, Om_L)
            FoV=comovSize_Mpc(x,inst) 
            
            j=0
            while (j<len(Rsmth)):
                pixel_size=2*(Rsmth[j]+0.00001)
                dnu= H_0*nu_0*sqrt(Om_m)/(c*sqrt(1.0+z1))*pixel_size;
                NumPixels = int( (FoV/pixel_size)**2.0*B/dnu )+0.0000001;
                ang_res = pixel_size/(pi/10800.0)/x.comovDist_Mpc()
                
                noise.append(getNoise(z1,Atot,ang_res,dnu,Tint));
                
                S2_err.append(  sqrt( getS2var(noise[j],sig2[j], NumPixels) )  );
                if norm==0:
                    S3_err.append(  sqrt( getS3var(noise[j],sig2[j], kurt[j], NumPixels,norm) )  );
                else:
                    S3_err.append(  sqrt( fabs(getS3normVar(noise[j], sig2[j], skew[j], kurt[j], NumPixels,norm)) )  );
                            
                if norm==1:
                    skew[j]/=(sig2[j]**(1.5)); 
                    #print skew[j], "\t", sig2[j];

                if norm==2:
                    skew[j]/=sig2[j];
                
                if (Revo-0.01<= Rsmth[j]<= Revo+0.01):
                    ionfstore.append(1.0-nf1);
                    sig2store.append(sig2[j]);
                    skewstore.append(skew[j]);
                
                j+=1
                               
            if (z1==z2plot):
                ax5.plot(Rsmth, sig2, linestyle=ln[i],marker=mk[i],
                         color=col[i], label=key[i]);
                ax6.plot(Rsmth, skew, linestyle=ln[i],marker=mk[i],
                         color=col[i], label=key[i]);
                if (i==0):
                    errkey=inst_str[inst]+' error'
                    ax5.fill_between(Rsmth, sig2, sig2+S2_err, color=errShade[1])
                    ax5.fill_between(Rsmth, sig2-S2_err, sig2, color=errShade[1])
                    ax5.axvspan(xmin=-10, xmax=-9,ymin=-10, ymax=-9, color=errShade[1], label=errkey)
                    
                    ax6.fill_between(Rsmth, skew, skew+S3_err, color=errShade[1])
                    ax6.fill_between(Rsmth, skew-S3_err, skew, color=errShade[1])   
                    ax6.axvspan(xmin=-10, xmax=-9,ymin=-10, ymax=-9, color=errShade[1], label=errkey)
        
                
    ax3.plot(ionfstore, sig2store, linestyle=ln[i],marker=mk[i],
             color=col[i], label= key[i]);
    ax4.plot(ionfstore, skewstore, linestyle=ln[i],marker=mk[i],
             color=col[i], label= key[i]);     
        
    i+=1
            
# ##loop to make the error shading in the variance/skew vs. neutral fraction plots. Have to reimport the 21cm data##

File=model[0]+'list'

l=0
while (l<len(inst_err)):
    Tint=instParams(l,'Tint')
    Atot=instParams(l,'Nstat')*instParams(l,'Aeff')
    ionfstore=[];
    sig2store=[];
    skewstore=[];
    S3errstore=[];
    S2errstore=[];
    with open(File) as a_file:         
        for a_line in a_file:
            filestr=a_line.rstrip()
            #extract the variables we need from the filename
            try:
                found = re.search('_z(.+?)_nf', filestr).group(1)
            except AttributeError:
                # _z, _use not found in the original string
                found = '' # apply your error handling
            z1=float(found)
            try:
                found = re.search('_nf(.+)', filestr).group(1)
            except AttributeError:
                # _nf, _z not found in the original string
                found = '' # apply your error handling
            nf1=float(found)
            
            (Rsmth, kurt, skew, sig2) = npy.loadtxt(filestr, delimiter='\t', usecols=(0,3,5,7), unpack=True);
           
            if (inst_err[l]==1):
                noise=[];
                S2_err=[];
                S3_err=[];
                
                x = Cosmology(z1, H_0, Om_m, Om_L)
                FoV=comovSize_Mpc(x,l) 

                j=0
                while (j<len(Rsmth)):
                    pixel_size=2*(Rsmth[j]+0.00001)
                    dnu= H_0*nu_0*sqrt(Om_m)/(c*sqrt(1.0+z1))*pixel_size;
                    NumPixels = int( (FoV/pixel_size)**2.0*B/dnu )+0.0000001;
                    ang_res = pixel_size/(pi/10800.0)/x.comovDist_Mpc()
                    noise.append(getNoise(z1,Atot,ang_res,dnu,Tint));
                        
                    S2_err.append(  sqrt( getS2var(noise[j],sig2[j], NumPixels) )  );
                    if norm==0:
                        S3_err.append(  sqrt( getS3var(noise[j],sig2[j], kurt[j], NumPixels,norm) )  );
                    else:
                        S3_err.append(  sqrt( fabs(getS3normVar(noise[j], sig2[j], skew[j], kurt[j], NumPixels,norm)) )  );
                             
                    if norm==1:
                        skew[j]/=(sig2[j]**(1.5));            
                    if norm==2:
                        skew[j]/=sig2[j];
        
                    if (Revo-0.1<= Rsmth[j]<= Revo+0.1):
                        ionfstore.append(1.0-nf1);
                        sig2store.append(sig2[j]);
                        skewstore.append(skew[j]);
                        S3errstore.append(S3_err[j]);
                        S2errstore.append(S2_err[j]);
                    j+=1
                    
        if (inst_err[l]==1):
            errkey=inst_str[l]+' error'
            ax3.fill_between(npy.array(ionfstore), npy.array(sig2store), npy.array(sig2store)+npy.array(S2errstore), color=errShade[l])
            ax3.fill_between(npy.array(ionfstore), npy.array(sig2store)-npy.array(S2errstore), npy.array(sig2store), color=errShade[l])
            ax3.axvspan(xmin=-10, xmax=-9,ymin=-10, ymax=-9, color=errShade[l], label=errkey)
        
            ax4.fill_between(npy.array(ionfstore), npy.array(skewstore), npy.array(skewstore)+npy.array(S3errstore), color=errShade[l])
            ax4.fill_between(npy.array(ionfstore), npy.array(skewstore)-npy.array(S3errstore), skewstore, color=errShade[l])   
            ax4.axvspan(xmin=-10, xmax=-9,ymin=-10, ymax=-9, color=errShade[l], label=errkey) 
    l+=1

# ## Make the plots perty and finish up ##

ax5.set_xlim([0,60])
ax5.set_ylim(vRange)     
ax5.legend(loc='best', fancybox=True,prop={'size':legFontsz});
ax5.grid(addGrid, which='major');

ax5.set_xlabel(r'\boldmath $R_{\mathrm{smooth}}$');
ax5.set_ylabel(r'brightness temperature \boldmath $\sigma^2 \,(\mathrm{mK}^2$) ');   

ax6.set_xlim([0,60])
ax6.set_ylim(sRange)  

ax6.legend(loc='best', fancybox=True,prop={'size':legFontsz});
ax6.grid(addGrid, which='major');

ax6.set_xlabel(r'\boldmath $R_{\mathrm{smooth}}$');

if norm==0:
    ax6.set_ylabel(r'Brightness temperature \boldmath $S_3 \,(\mathrm{mK}^3)$'); 
if norm==1:
    ax6.set_ylabel(r'Brightness temperature \boldmath $S_3/\sigma^3$'); 
if norm==2:
    ax6.set_ylabel(r'Brightness temperature \boldmath $S_3/\sigma^2 \,(\mathrm{mK})$');
    #Save & show plot then clean-up!

PlotFile=saveDIR+'dTvariance_smoothed_z'+str(Decimal(z2plot).quantize(TWOPLACES))+'_'+str(Res)+'_'+str(L)+'Mpc.'+saveFormat
#fig5.savefig(PlotFile, format=saveFormat);
PlotFile=saveDIR+'dTskew_smoothed_z'+str(Decimal(z2plot).quantize(TWOPLACES))+'_norm'+str(norm)+'_'+str(Res)+'_'+str(L)+'Mpc.'+saveFormat
#fig6.savefig(PlotFile, format=saveFormat);

ax3.set_xlim(xRange)
ax3.set_ylim(vRange)     
ax3.legend(loc='best', fancybox=True,prop={'size':legFontsz});
ax3.grid(addGrid, which='major');

ax5 = ax3.twin() # ax2 is responsible for "top" axis and "right" axis
ax5.set_xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
ax5.set_xticklabels(['12.25','','9.75','','8.75','','8.00','','7.25',''])
ax5.set_xlabel(r'Redshift')
ax5.axis["right"].major_ticklabels.set_visible(False) 

ax3.set_xlabel(r'Average ionised fraction');
ax3.set_ylabel(r'brightness temperature \boldmath $\sigma^2 \,(\mathrm{mK}^2$) ');   

ax4.set_xlim(xRange)
ax4.set_ylim(sRange) 
    
ax4.legend(loc='best', fancybox=True,prop={'size':legFontsz});
ax4.grid(addGrid, which='major');

ax6 = ax4.twin() # ax2 is responsible for "top" axis and "right" axis
ax6.set_xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
ax6.set_xticklabels(['12.25','','9.75','','8.75','','8.00','','7.25',''])
ax6.set_xlabel(r'Redshift')
ax6.axis["right"].major_ticklabels.set_visible(False) 

ax4.set_xlabel(r'Average ionised fraction');
    
if norm==0:
    ax4.set_ylabel(r'Brightness temperature \boldmath $S_3 \,(\mathrm{mK}^3)$'); 
if norm==1:
    ax4.set_ylabel(r'Brightness temperature \boldmath $S_3/\sigma^3$'); 
if norm==2:
    ax4.set_ylabel(r'Brightness temperature \boldmath $S_3/\sigma^2 \, (\mathrm{mK})$');
#Save & show plot then clean-up!
    
PlotFile=saveDIR+'dTvariance_smoothed_R'+str(Revo)+'_norm'+str(norm)+'_'+str(Res)+'_'+str(L)+'Mpc.'+saveFormat
fig3.savefig(PlotFile, format=saveFormat);

PlotFile=saveDIR+'dTskew_smoothed_R'+str(Revo)+'_norm'+str(norm)+'_'+str(Res)+'_'+str(L)+'Mpc.'+saveFormat
fig4.savefig(PlotFile, format=saveFormat);

plt.cla();
plt.clf(); # make sure any background plots are wiped



print "------ That's all folks!!! ------"

