# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1

This file should contain all functions to read data from the output of MITgcm to smaller files that are useful for the analysis.

@author: joren
"""

#%%IMPORTING
import xarray as xr
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pandas as pd
import datetime
import matplotlib.colors as colors


def performPCA(data, var, method='all', interest=10, filepath='./figures/', filename=''):
    import sklearn.decomposition as pca
    from scipy import signal
    from numpy import ma
    
    data=data.drop_isel(time=np.where(np.isnan(data.values[0,:,150,150]))[0])
    
    if method=='ens_mean':
        data=data.mean(dim='ens')
        data['time']=data.indexes['time'].to_datetimeindex()
        tc=data.time
    elif method=='all':
        data['time']=data.indexes['time'].to_datetimeindex()
        tc=data.time
        dims=[len(data.ens), len(data.time)]
        data=data.stack(ens_time=['ens', 'time'])
        data=data.transpose("ens_time", "YC", "XC")
        
        
    sla=data.values
    if 'YC' in data.dims:
        lat=data.YC
    else:
        lat=data.YG
        
    if 'XC' in data.dims:
        lon=data.XC
    else:
        lon=data.XG
        
    print(np.shape(sla))
    
    #---------------------------------------------------------------------------------------------------------------------------
    #PERFORM PCA

    #Reshape the grid to one column with position and row corresponding to time
    X = np.reshape(sla, (sla.shape[0], len(lat) * len(lon)), order='F')
    #Create Boolean array corresponding to ocean indices and mask the array where there is land
    X = ma.masked_array(X, np.isnan(X))
    land = X.sum(0).mask
    X[:,land]=0
    total_modes = 10
    EOF_recons = np.ones((total_modes, len(lat) * len(lon))) * -999.

    #Perform PCA and take the EOFs and corersponding Pcs
    skpca = pca.PCA(total_modes)
    skpca.fit(X)
    PCs = skpca.transform(X)
    EOFs = skpca.components_
    #Reconstruct the EOF fields into gridded data EOF_recons =
    np.ones((total_modes, len(lat) * len(lon))) * -999.
    for i in range(total_modes):
        EOF_recons[i,:] = EOFs[i,:]
    EOF_recons = ma.masked_values(np.reshape(EOF_recons, (total_modes, len(lat), len(lon)), order='F'), -999.)
    print(EOF_recons.shape)
    
    
    #---------------------------------------------------------------------------------------------------------------------------
    #VISUALIZE
    
    ### plot loadings and scores
    from matplotlib.colors import Normalize
    import matplotlib.colors as colors

    

    
    for j in range(0,interest):
        #norm = MidpointNormalize(midpoint=0)
        #Plot principal components
        fig, ax = plt.subplots(2,1,figsize=(12,12))
        EOF_recons[j,:,:][EOF_recons[j,:,:]==0]=np.nan
        #im = ax[0].imshow(np.flip(EOF_recons[j,:,:]*max(abs(PCs[j,:])), axis=0), norm=divnorm,
        #cmap=plt.cm.RdBu_r, interpolation='none',extent=[lon[0],lon[-
        #1],lat[0],lat[-1]])
        divnorm = colors.TwoSlopeNorm(vmin=-np.nanquantile(abs(EOF_recons[j,:,:]*max(abs(PCs[j,:]))), 0.95), vmax=np.nanquantile(abs(EOF_recons[j,:,:]*max(abs(PCs[j,:]))), 0.95), vcenter=0)
        im = ax[0].pcolor(lon, lat, EOF_recons[j,:,:]*max(abs(PCs[j,:])), norm=divnorm,
        cmap=plt.cm.RdBu_r)

        ax[0].set_facecolor('grey')
        cbar=plt.colorbar(im, fraction=0.016, pad=0.04, ax=ax[0])
        cbar.set_label('(m)')
        ax[0].set_xlabel(r'Longitude ($^\circ$E)')
        ax[0].set_ylabel(r'Latitude ($^\circ$N)')

        eof=j+1
        ax[0].set_title(r'Empirical Orthogonal Function %i' % eof)
        PCnor=PCs[:,j]/max(abs(PCs[:,j]))
        
        if method=='all':
            PCnor=PCnor.reshape(dims[0], dims[1])
            for k in range(dims[0]):
                if k==0:
                    ax[1].plot(tc,PCnor[k,:], c='grey', linewidth=0.2, label='Individual member PC amplitude')
                else:
                     ax[1].plot(tc,PCnor[k,:], c='grey', linewidth=0.2)
            ax[1].plot(tc,np.mean(PCnor, axis=0), c='black', linewidth=1.2, label='Ensemble averaged PC amplitude')
        else:
            ax[1].plot(tc, PCnor, label='PC')
        #ax[1].plot(tc, (mean-mean.mean())/max(abs(mean-mean.mean())), label='Shelf Mean: corr:{}'.format(str(np.corrcoef(mean, PCnor)[0,1])))
        ax[1].grid()
        ax[1].set_title('Principal Component %i' % eof)
        ax[1].set_xlabel('Time (year)')
        ax[1].set_ylabel('Magnitude')
        #ax[1].set_xlim([tc[0],tc[-1]])
        ax[1].set_ylim([-1,1])
        ax[1].legend()
        
        from datetime import date
        today = date.today()
        today=today.strftime("%Y%m%d")
        
        plt.savefig(filepath+today+'_'+var+"_PCA_component_"+str(eof)+'_'+method+'_'+filename+".png", bbox_inches='tight')
        plt.show()
    
    variance=skpca.explained_variance_ratio_
    variancecumsum=np.cumsum(variance)
    eofnumber=np.arange(1,len(variance)+1)

    f, ax = plt.subplots(figsize=(6,4))
    ax.plot(eofnumber,variancecumsum)
    ax.plot(eofnumber[0:5],variancecumsum[0:5],'ro',markersize=4)
    ax.set_xlabel('EOF modes')
    ax.set_ylabel('Cumulative fraction')
    ax.set_title("Cumulative fraction of explained variance",
    fontsize=14)
    ax.grid()
    #ax.set_xlim([-1,313])
    #ax.set_xticks([1,50,100,150,200,250,300])
    plt.savefig(filepath+today+'_'+var+"cumulativefraction"+method+'_'+filename+".png", bbox_inches='tight')
    plt.show()
    print("Variance contained in first 10 EOFs: % .2f" %
    variancecumsum[9])
    
    np.save('../02_data/principal_components/'+today+'_'+var+"_PCA_component"+method+'_'+filename+".npy", PCs)
    EOF_recons.dump('../02_data/principal_components/'+today+'_'+var+"_EOF_component"+method+'_'+filename+".npy")
    return EOF_recons, PCs, tc


def makePCAAnimation(data_dt, EOF_recons, PCs, tc, var, pcnumber=10, frames=24, start_year=1939, filename='test', filepath='../03_output/pca/'):
    
    #From my code to Kaitlin's code
    #from mitgcm_python_master.plot_latlon import latlon_plot
    #from mitgcm_python_master.grid import ERA5Grid, PACEGrid, Grid, dA_from_latlon, pierre_obs_grid
    #from mitgcm_python_master.utils import real_dir, daily_to_monthly, fix_lon_range, split_longitude, mask_land_ice, moving_average, index_year_start, index_year_end, index_period, mask_2d_to_3d, days_per_month, add_time_dim, z_to_xyz, select_bottom, convert_ismr, mask_except_ice, xy_to_xyz, apply_mask, var_min_max, mask_3d, average_12_months, depth_of_isoline, mask_land, axis_edges
    #from mitgcm_python_master.plot_utils.latlon import shade_mask, overlay_vectors
    #from mitgcm_python_master.plot_utils.windows import finished_plot, set_panels
    
    ### plot loadings and scores
    from matplotlib.colors import Normalize
    import matplotlib.colors as colors
    from matplotlib.animation import FuncAnimation
    import matplotlib.animation as animation
    
    if 'YC' in data_dt.dims:
        lat=data_dt.YC
    else:
        lat=data_dt.YG
        
    if 'XC' in data_dt:
        lon=data_dt.XC
    else:
        lon=data_dt.XG
    
    fig, ax=plt.subplots(figsize=(14,5))
    start_month=int((start_year-1920)*12)
    k=start_month
    
    v=np.max([abs(float(data_dt.quantile(0.98))),abs(float(data_dt.quantile(0.02)))])
    
    divnorm = colors.TwoSlopeNorm(vmin=-v, vmax=v, vcenter=0)

    fig, ax = plt.subplots(1,1,figsize=(18,12))
    
    #Create Initial Frame
    full_recon=np.zeros(np.shape(EOF_recons[0,:,:]))
    for j in range(0, pcnumber):
        #Plot principal components
        EOF_recons[j,:,:][EOF_recons[j,:,:]==0]=np.nan
        full_recon+=EOF_recons[j,:,:]*PCs[k,j]
    
    im = ax.pcolor(lon, lat, full_recon, norm=divnorm,
    cmap=plt.cm.RdBu_r)
    ax.set_facecolor('grey')
    cbar=plt.colorbar(im, fraction=0.016, pad=0.04, ax=ax)
    cbar.set_label('(m)')
    ax.set_xlabel(r'Longitude ($^\circ$E)')
    ax.set_ylabel(r'Latitude ($^\circ$N)')
    ax.set_title(var+': Animation of first {} PCs'.format(str(pcnumber))+'; Time: '+str(tc[start_month].values)[:10])

    # animation function
    def animate(i):
        k=start_month+i
        print(str(tc[k].values)[:10])
        #Create Initial Frame
        full_recon=np.zeros(np.shape(EOF_recons[0,:,:]))
        for j in range(0, pcnumber):
            #Plot principal components
            EOF_recons[j,:,:][EOF_recons[j,:,:]==0]=np.nan
            full_recon+=EOF_recons[j,:,:]*PCs[k,j]
            
        im = ax.pcolor(lon, lat, full_recon, norm=divnorm,
        cmap=plt.cm.RdBu_r)
        ax.set_facecolor('grey')
        #cbar=plt.colorbar(im, fraction=0.016, pad=0.04, ax=ax)
        #cbar.set_label('(m)')
        ax.set_xlabel(r'Longitude ($^\circ$E)')
        ax.set_ylabel(r'Latitude ($^\circ$N)')
        ax.set_title(var+': Animation of first {} PCs'.format(str(pcnumber))+'; Time: '+str(tc[k].values)[:10])

        return 

    
    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    
    anim = animation.FuncAnimation(fig, animate, frames=frames, repeat=False, interval=5000)
    writergif = animation.PillowWriter(fps=2)
    anim.save(filepath+'animation_pca_'+filename+'_'+today+'.gif', writer=writergif)
    return


def plotPCACorrs(PCs, total, real_var, data_dt, var_list=['amundsen_shelf_break_uwind_avg', 'amundsen_shelf_temp_btw_200_700m', 'dotson_to_cosgrove_massloss'], interest=10, filename='', filepath='./figures/'):
    
    #print(cors)
    plt.figure(figsize=(15,10))
    
    for i, var in enumerate(var_list):
        total[var] = total[var].reindex(sorted(total[var].columns), axis=1)
        ts=total[var]['1920':'2013']
        ts=ts.rolling(24,center=True).mean()-ts.rolling(25*12, min_periods=5*12, center=True).mean()
        ts=ts.iloc[12:-11]
        mean=ts.to_numpy().T.flatten()
        
        
        #print(len(mean))
        cors=np.zeros(interest)
        for j in range(0,interest):
            PCnor=PCs[:,j]/max(abs(PCs[:,j]))
            cors[j]=np.corrcoef(mean, PCnor)[0,1]
        
        plt.subplot(len(var_list),1,i+1)
        plt.plot(range(1,interest+1), cors, 'o-', label='Pearson r')
        plt.plot(range(1,interest+1), abs(cors), 'o-', label='Absolute values')
        #plt.yticks([-1,0,1])
        plt.grid()
        plt.legend()
        plt.ylabel('Pearson R')
        plt.xlabel('Principal Component #')
        plt.title('Correlation between '+var+' and PCs from '+real_var)
    
    from datetime import date
    today = date.today()
    today=today.strftime("%Y%m%d")
    
    plt.subplots_adjust(hspace=0.5)
    plt.savefig(filepath+today+"_PCAcorrs_"+real_var+'_'+filename+".png", bbox_inches='tight')
    
    return