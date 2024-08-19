# -*- coding: utf-8 -*-
"""
Created on Fri May  1 10:44:32 2020
Script for performing regression analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import itertools
from scipy import stats
import scipy.optimize
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error, cohen_kappa_score, accuracy_score
from sklearn.model_selection import KFold
import plotly.express as px


def linear(x,a,b):
    return a*x+b


def quadratic(x,a,b,c):
    return a*x**2 + b*x + c


def linear2var(x,a,b,c):
    x0 = x.iloc[:,0]
    x1 = x.iloc[:,1]
    return a*x0+b*x1+c


def linear3var(x,a,b,c,d):
    x0 = x.iloc[:,0]
    x1 = x.iloc[:,1]
    x2 = x.iloc[:,2]
    return a*x0+b*x1+c*x1+d


def kfold_regression(data, x, y, survey, crop, fn, splits=5, bins=None):
    print(f'Running: {", ".join(x)} for {crop} using a {fn.__name__} regression')
    data = data.loc[
                    (data['surveyType'] == survey) & 
                    (data['cdlBack1'] == crop),
                    x+y
                    ].dropna()
    r2s = []
    rmses = []
    kappas = []
    accs = []
    kf = KFold(n_splits=splits, random_state=10051993, shuffle=True)
    for train, test in kf.split(data):
        # Split data
        if data.shape[1] <= 2:
            x_train = data.iloc[train,0]
            x_test = data.iloc[test,0]
            y_train = data.iloc[train,1]
            y_test = data.iloc[test,1]
        else:
            x_train = data.iloc[train,0:-1]
            x_test = data.iloc[test,0:-1]
            y_train = data.iloc[train,-1]
            y_test = data.iloc[test,-1]

        # Perform Regression
        fit_data, covariance = scipy.optimize.curve_fit(fn, xdata=x_train, ydata=y_train, maxfev=1000000)

        # Measure Performance
        df_test = pd.DataFrame({'Test':y_test,'Predict':fn(x_test, *fit_data)})
        if bins != None:
            # df_test['Test'] = pd.cut(df_test['Test'], bins['bins'], labels=bins['categories'])
            # df_test['Predict'] = pd.cut(df_test['Predict'], bins['bins'], labels=bins['categories'])
            # kappa = cohen_kappa_score(df_test['Test'].astype(str), df_test['Predict'].astype(str))
            cat_test = pd.cut(df_test['Test'], bins['bins'], labels=bins['categories']).astype(str)
            cat_predict = pd.cut(df_test['Predict'], bins['bins'], labels=bins['categories']).astype(str)
            kappa = cohen_kappa_score(cat_test, cat_predict)
            kappas.append(kappa)
            acc = accuracy_score(cat_test, cat_predict, normalize=True)
            accs.append(acc)
            
        r2 = r2_score(df_test['Test'], df_test['Predict'])
        r2s.append(r2)
        rmse = np.sqrt(mean_squared_error(df_test['Test'], df_test['Predict']))
        rmses.append(rmse)

    line = pd.Series({'Survey':survey, 'Covariate':', '.join(x), 'Crop':crop, 'Model':fn.__name__,
        'R2':np.mean(r2s), 'RMSE':np.mean(rmses), 'N_OBS':data.shape[0]})

    if bins != None:
        line['Kappa'] = np.mean(kappas)
        line['Accuracy'] = np.mean(accs)

    return line


def drop_low_img(data,xs,n=1):
    if isinstance(xs, str):
        xs = [xs]

    for x in xs:
        if 'SMAP' in x:
            data = data[(data['SMAPcount'] > n)]
        else:
            data = data[(data[f'{x[0:2]}count'] > n)]

    return data


def main():
    # Load initial options and data
    data_path = os.path.join('tabular_data','survey_and_index_data.csv')
    data = pd.read_csv(data_path)
    crops = ['corn', 'soybeans']
    surveys = ['Windshield Tillage', 'Windshield Binned Residue', 'Photo Analysis', 'Photo Analysis Binned']
    y = ['PctResidue']
    tillage = {'bins':[-np.inf, 15, 30, 50, np.inf],'categories':[7.5, 22.5, 40, 75]}

    # For plotting
    crop_color ={'corn':'#F1BE48', 'soybeans':'#C8102E'}
    sns.set('poster', font='serif', style='white')

    # 3.1 Google Earth Engine Performance~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Number of valid images
    s1 = data[(data['S1VVmin'].notna()) & (data['S1count'] > 1)]
    s2 = data[(data['S2NDTImin'].notna()) & (data['S2countfull'] > 1)]
    l7 = data[(data['L7NDTImin'].notna()) & (data['L7countfull'] > 1)]
    l8 = data[(data['L8NDTImin'].notna()) & (data['L8countfull'] > 1)]
    smap = data[(data['SMAPmin'].notna()) & (data['SMAPcount'] > 1)]
    
    valid = s1.append([s2,l7,l8,smap])
    valid = valid[['state', 'surveyType','cdlBack1','S1count','S2count','L7count','L8count',
                    'SMAPcount','S2countfull','L7countfull','L8countfull']]
    
    median_img = valid.groupby(['state', 'surveyType','cdlBack1']).median()
    median_img.to_csv('median_img.csv')

    # 3.2 Model Selection~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # setup
    ndtis = ['L7NDTImin', 'L7NDTImedian', 'L8NDTImin', 'L8NDTImedian','LCNDTImin', 'LCNDTImedian', 'S2NDTImin', 'S2NDTImedian']
    ndvis = ['L7NDVImin', 'L7NDVImedian', 'L8NDVImin', 'L8NDVImedian', 'LCNDVImin', 'LCNDVImedian', 'S2NDVImin', 'S2NDVImedian']
    swcs = ['L7SWCmin', 'L7SWCmedian', 'L8SWCmin', 'L8SWCmedian', 'LCSWCmin', 'LCSWCmedian', 'S2SWCmin', 'S2SWCmedian']
    others = ['SMAPmin','SMAPmedian','S1VVmin', 'S1VVmedian', 'S1VHmin', 'S1VHmedian','S1NDImin', 'S1NDImedian']
    swcs = swcs + others

    # NDTI~~~~~~~~~~~~~~~~~~~~~~#
    ndti_df = pd.DataFrame()
    for x in ndtis:
        data_highcount = drop_low_img(data,x)
        for crop in crops:
            line = kfold_regression(data_highcount, [x], y, 'Photo Analysis', crop, linear, bins=tillage)
            ndti_df = ndti_df.append(line,ignore_index=True)

    print(ndti_df)
    ndti_df = ndti_df.sort_values(['Crop','Accuracy'], ascending=[True,False])
    ndti_means = ndti_df.groupby(['Covariate','Model']).mean().reset_index()

    # Save to excel
    with pd.ExcelWriter('ndtis.xlsx') as writer:
        ndti_df.to_excel(writer, sheet_name='ndti', index=False)
        ndti_means.to_excel(writer, sheet_name='means', index=False)

    # Graph
    # Accuracy bar graph
    rename_x = {'L7NDTImin':'L7 NDTI Min','L8NDTImin':'L8 NDTI Min','LCNDTImin':'L7/8 NDTI Min','S2NDTImin':'S2 NDTI Min'}
    rename_hue = {x:x.capitalize() for x in crop_color}
    crop_color_cap = {x.capitalize():crop_color[x] for x in crop_color}

    graph_df = ndti_df.loc[(ndti_df['Covariate'].str.contains('min'))&(ndti_df['Model']=='linear')].copy()
    graph_df['Accuracy (%)'] = graph_df['Accuracy'].multiply(100).astype(int)
    graph_df['Covariate'] = graph_df['Covariate'].map(rename_x)
    graph_df['Crop'] = graph_df['Crop'].map(rename_hue)

    f, ax = plt.subplots(figsize=(15, 10))
    sns.barplot(data = graph_df, x='Covariate', y='Accuracy (%)', hue='Crop', palette=crop_color_cap, ax=ax)
    f.savefig('ndti_by_sat_acc.svg')

    # Graph all survey types
    y_value = 'PctResidue'
    x_values = {'Landsat 7':'L7NDTImin','Landsat 8':'L8NDTImin','Landsat 7/8':'LCNDTImin','Sentinel2':'S2NDTImin'}
    
    f, axes = plt.subplots(2, 2, figsize=(20, 20))
    subs = [(0,0),(0,1),(1,0),(1,1)]
    for i, x_name in zip(subs,x_values):
        x_value = x_values[x_name]
        ax = axes[i]
        for crop in crops:
            data_highcount = drop_low_img(data,x_value,n=1)
            data_plot = data_highcount.loc[(data_highcount['surveyType'] == 'Photo Analysis')&
                                        (data_highcount['cdlBack1'] == crop),
                                        [x_value,y_value]
                                        ].dropna()
            sns.regplot(x=x_value, y=y_value, data=data_plot, ax=ax, color=crop_color[crop], ci=None, label = crop.capitalize())
        
        ax.set_xlabel('Minimum NDTI')
        ax.set_ylabel('Percent Residue Cover')
        ax.set(ylim=(0, 105), xlim=(0, .17), xticks=[0.05,0.10], title=x_name)
        if i == (0,0):
            ax.legend(loc='upper left')
    f.savefig('ndti_by_sat_reg.svg')

    # Linear vs quadratic regression~~~~~~~~~~~~~~~~~~~~~~#
    x_values = {'corn':'LCNDTImin','soybeans':'L7NDTImin'}
    y_value = 'PctResidue'

    function_df = pd.DataFrame()
    models = [linear, quadratic]
    for model in models:
        for crop in crops:
            x_value = x_values[crop]
            data_highcount = drop_low_img(data,x_value,n=1)
            line = kfold_regression(data_highcount, [x_value], [y_value], 'Photo Analysis', crop, model, bins=tillage)
            function_df = function_df.append(line,ignore_index=True)

    function_means = function_df.groupby('Model').mean().reset_index()
    print(function_df)

    # Save to excel
    with pd.ExcelWriter('functions.xlsx') as writer:
        function_df.to_excel(writer, sheet_name='functions', index=False)
        function_means.to_excel(writer, sheet_name='means', index=False)
    
    # Graph
    models = {'Linear':1, 'Quadratic':2}
    
    f, axes = plt.subplots(1, 2, figsize=(20, 10))
    subs = [0,1]
    for i, model in zip(subs,models):
        ax = axes[i]
        for crop in crops:
            data_plot = data_highcount.loc[(data_highcount['surveyType'] == 'Photo Analysis') & (data_highcount['cdlBack1'] == crop)]
            sns.regplot(x=x_value, y=y_value, data=data_plot, ax=ax, order=models[model], color=crop_color[crop], ci=None, label = crop.capitalize())
        ax.set_xlabel('Minimum Landsat 7/8 NDTI')
        ax.set_ylabel('Percent Residue Cover')
        ax.set(ylim=(0, 105), xlim=(0, .17), xticks=[0.05,0.10], title=model)
        if i == 0:
            ax.legend(loc='upper left')
    f.savefig('functions.svg')


    # NDTI with SWC or NDVI separately~~~~~~~~~~~~~~~~~~~~~~#
    twovar_df = pd.DataFrame()
    xs = [list(a) for a in itertools.product(['L7NDTImin','L8NDTImin','LCNDTImin','S2NDTImin'], swcs+ndvis)]
    for x in xs:
        data_highcount = drop_low_img(data,x)
        for crop in crops:
            line = kfold_regression(data_highcount, x, y, 'Photo Analysis', crop, linear2var, bins=tillage)
            twovar_df = twovar_df.append(line,ignore_index=True)

    twovar_df = twovar_df.sort_values(['Crop','Accuracy'], ascending=[True,False])
    twovar_means = twovar_df.groupby('Covariate').mean().reset_index()
    print(twovar_df)

    # Save to excel
    with pd.ExcelWriter('twovar.xlsx') as writer:
        twovar_df.to_excel(writer, sheet_name='twovar', index=False)
        twovar_means.to_excel(writer, sheet_name='means', index=False)

    # NDTI with SWC and NDVI
    threevar_df = pd.DataFrame()
    xs = [list(a) for a in itertools.product(['L7NDTImin','L8NDTImin','LCNDTImin','S2NDTImin'], swcs, ndvis)]
    for x in xs:
        data_highcount = drop_low_img(data,x)
        for crop in crops:
            line = kfold_regression(data_highcount, x, y, 'Photo Analysis', crop, linear3var, bins=tillage)
            threevar_df = threevar_df.append(line,ignore_index=True)

    threevar_df = threevar_df.sort_values(['Crop','Accuracy'], ascending=[True,False])
    threevar_means = threevar_df.groupby('Covariate').mean().reset_index()
    print(threevar_df)

    # Save to excel
    with pd.ExcelWriter('threevar.xlsx') as writer:
        threevar_df.to_excel(writer, sheet_name='threevar', index=False)
        threevar_means.to_excel(writer, sheet_name='means', index=False)


    # Section 3.3 Utility of survey data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    y_value = 'PctResidue'
    x_soybean_value =['LCNDTImin','L7SWCmedian']
    x_corn_value = ['LCNDTImin','LCNDVImin']
    data = data[~((data['PctResidue'] == 0)&(data['surveyType']=='Windshield Binned Residue'))]

    survey_df = pd.DataFrame()
    for survey in surveys:
        # corn
        data_highcount = drop_low_img(data,x_corn_value,n=1)
        line = kfold_regression(data_highcount, x_corn_value, [y_value], survey, 'corn', linear2var, splits=5, bins=tillage)
        survey_df = survey_df.append(line,ignore_index=True)

        # soybeans
        data_highcount = drop_low_img(data, x_soybean_value,n=1)
        line = kfold_regression(data_highcount, x_soybean_value, [y_value], survey, 'soybeans', linear2var, splits=5, bins=tillage)
        survey_df = survey_df.append(line,ignore_index=True)

    survey_means = survey_df.groupby('Survey').mean().reset_index()

    # Save to excel
    with pd.ExcelWriter('surveys.xlsx') as writer:
        survey_df.to_excel(writer, sheet_name='surveys', index=False)
        survey_means.to_excel(writer, sheet_name='means', index=False)
    print(survey_df)

    # Graph all survey types
    x_value = 'LCNDTImin'
    y_value = 'PctResidue'
    
    f, axes = plt.subplots(2, 2, figsize=(20, 20))
    subs = [(0,0),(0,1),(1,0),(1,1)]
    for i, survey in zip(subs,surveys):
        ax = axes[i]
        for crop in crops:
            data_highcount = drop_low_img(data,x_value,n=1)
            data_plot = data_highcount.loc[(data_highcount['surveyType'] == survey)&
                                        (data_highcount['cdlBack1'] == crop),
                                        [x_value,y_value]
                                        ].dropna()
            sns.regplot(x=x_value, y=y_value, data=data_plot, ax=ax, color=crop_color[crop], ci=None, label = crop.capitalize())
        
        ax.set_xlabel('Minimum Landsat 7/8 NDTI')
        ax.set_ylabel('Percent Residue Cover')
        ax.set(ylim=(0, 105), xlim=(0, .17), xticks=[0.05,0.10], title=survey)
        if i == (0,0):
            ax.legend(loc='upper left')
    f.savefig('surveys.svg')

    # Graph tillage type kde
    x_value = 'LCNDTImin'
    data_highcount = drop_low_img(data,x_value,n=1)
    data_till = data_highcount.loc[(data_highcount['surveyType']=='Windshield Tillage')&
                                    (data_highcount['PctResidue'].isin([3,4,55,72])),
                                    ['cdlBack1','PctResidue','surveyType',x_value]].copy().reset_index()

    data_till = data_till.loc[~((data_till['PctResidue']==55)&(data_till['cdlBack1']=='corn'))]
    data_till['Tillage Group'] = 'No Till'
    data_till.loc[data_till['PctResidue'].isin([3,4]),'Tillage Group'] = 'Conventional Till'

    data_bin = data_highcount.loc[(data_highcount['surveyType']=='Windshield Binned Residue')&
                                    (data_highcount['PctResidue'].isin([5,10,55,65,75,85,95,100])),
                                    ['cdlBack1','PctResidue','surveyType',x_value]].copy().reset_index()
    data_bin['Residue Bin Group'] = r'55% - 100%'
    data_bin.loc[data_bin['PctResidue'].isin([5,10]),'Residue Bin Group'] = r'0% - 10%'

    f, axes = plt.subplots(1, 2, figsize=(20, 10))
    ax0, ax1 = axes

    sns.kdeplot(x=data_till[x_value],hue=data_till['Tillage Group'],common_norm=False,legend=True,ax=ax0)
    ax0.set_xlabel('Minimum Landsat 7/8 NDTI')

    sns.kdeplot(x=data_bin[x_value],hue=data_bin['Residue Bin Group'],common_norm=False,legend=True,ax=ax1)
    ax1.set_xlabel('Minimum Landsat 7/8 NDTI')
    f.savefig('kde.svg')
    
    return None


if __name__ == '__main__':
    main()
