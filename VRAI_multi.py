#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 17:22:53 2023

@author: Ching Ching Lam

extension to VRAI selectivity programme:
    additional processes for treating systems with more than two products that share the same intermediate structure 
    

"""

from VRAI_selectivity_v7 import main1
import pandas as pd 
from datetime import date
import os 
import itertools

def get_binary_combination(stuff): 
    
    ## stuff = a list of items 
    
    binary_ls=[]
    for subset in itertools.combinations(stuff, 2):            
        binary_ls.append(subset)
                
    return binary_ls

def get_file_df(path):
    
    ## Generate an input csv file with input file information – 
    ## has considered the possible binary combinations considering the number of possible products with get_binary_combination()

    
    TS1_gu_files = [f for f in os.listdir(path) if f.endswith('TS1.out')]
    TS2_gu_files = sorted([f for f in os.listdir(path) if f.endswith('TS2.out')])
    int_gu_files = [f for f in os.listdir(path) if f.endswith('int.out')]

    prod_mol_files = sorted([f for f in os.listdir(path) if f.endswith('prod.mol')])
    TS1_mol_files = [f[:-4]+'.mol' for f in TS1_gu_files]
    int_mol_files = [f for f in os.listdir(path) if f.endswith('int.mol')]


    binary_idx_ls=get_binary_combination([i for i in range(0,len(prod_mol_files))])
    binary_ls1=[ls[0] for ls in binary_idx_ls]
    binary_ls2=[ls[1] for ls in binary_idx_ls]

    prod_ls1=[prod_mol_files[idx] for idx in binary_ls1]
    prod_ls2=[prod_mol_files[idx] for idx in binary_ls2]
    
    if len(TS2_gu_files)<1:
        TS2_ls1=['NA']*len(prod_ls1)
        TS2_ls2=['NA']*len(prod_ls1)
    
    else:
        TS2_ls1=[TS2_gu_files[idx] for idx in binary_ls1]
        TS2_ls2=[TS2_gu_files[idx] for idx in binary_ls2]


    int_gls=int_gu_files*len(binary_idx_ls)
    int_mls=int_mol_files*len(binary_idx_ls)
    
    
    file_df_ls=[]
    
    
    for gu,mol in zip(TS1_gu_files,TS1_mol_files):
        
        TS1_gls=[gu]*len(binary_idx_ls)
        TS1_mls=[mol]*len(binary_idx_ls)

        sub_file_df=pd.DataFrame({'TS1':TS1_mls,'int':int_mls,'p1':prod_ls1,'p2':prod_ls2,'TS1_freq':TS1_gls,'int_freq':int_gls,
              'ts2a_freq':TS2_ls1,'ts2b_freq':TS2_ls2})
        
        
        
        file_df_ls.append(sub_file_df)
        
    
    file_df=pd.concat(file_df_ls)
    
    
    return file_df



def get_ratio_df(test_df,prod):
    
    ## given that there's three product: A, B, C
    ## when prod = A; Vrai selectivity calculation for A&B and A&C will be used 
    
    test_df_sub=test_df[(test_df['major'] == prod) | (test_df['minor'] == prod)]

    ls1 = [[j,i] for i,j in zip(test_df_sub['major_perc'], test_df_sub['major'])]
    ls2 = [[j,i] for i,j in zip(test_df_sub['minor_perc'], test_df_sub['minor'])]


    result_ls =[ [i,j] for i,j in zip(ls1,ls2)]

    result_ls1 = [ i for i in result_ls[0]]

    ratio1 = [i for i in result_ls[0] if i[0] == prod][0][1]

    for idx in range(1,len(result_ls)):
        a=result_ls[idx][0]
        b=result_ls[idx][1]
    
        if a[0] == prod:
            try:
                ratio_op= (ratio1/a[1])*b[1]
                result_ls1.append([b[0], ratio_op])
            except ZeroDivisionError:
                print('ZeroDivisionError')
            
    
        elif b[0] == prod:
            try:
                ratio_op= (ratio1/b[1])*a[1]
                result_ls1.append([a[0], ratio_op])
            
            except ZeroDivisionError:
                print('ZeroDivisionError')

    sum_ratio=sum([i[1] for i in result_ls1])

    product_ls=[i[0] for i in result_ls1]
    ratio_ls=[(i[1]/sum_ratio)*100 for i in result_ls1]

    ratio_df=pd.DataFrame({'product':product_ls, 'ratio': ratio_ls})
    ratio_df.sort_values('ratio', ascending=False)
    
    return ratio_df

def all_possible_ratio_df(test_df):
    
    ## repeat the get_ratio_df(test_df,prod) on all possible product
    
    prod_ls=test_df['major'].unique().tolist()+test_df['minor'].unique().tolist()
    prod_ls1=[]
    for i in prod_ls:
        
        if i not in prod_ls1: 
            prod_ls1.append(i)


    ratio_df_ls=[]
    for p in prod_ls1:
        df=get_ratio_df(test_df,p)
        df['based_on'] = p
        df1=df.sort_values('ratio', ascending=False)
    
        ratio_df_ls.append(df1)
    
    ratio_df=pd.concat(ratio_df_ls)

    return ratio_df




class VRAI_multi:
    
    def __init__(self, file_info):
        
        self.file_info=file_info
        
        if file_info[-3:] == 'csv':
            
            self.file_df=pd.read_csv(file_info)
            
        else:
            
            self.file_df=get_file_df(file_info)
            
    
    def perform_vrai(self, path='', intermediate_option=True, weight_option=True, TST_option=True, with_spe=False, temperature = 298.0, energy_cut_off=5):
        
        ## perform vrai selectivity calculation 
        
        columns_ls=self.file_df.columns.tolist()
        
        if self.file_info[-3:] != 'csv':
            path=self.file_info


        all_file_ls=[[path+i for i in self.file_df[n].tolist()] for n in columns_ls]


        result_df_ls=[]

        for idx in range(0,len(all_file_ls[0])):
            
            if intermediate_option==True:
                
                file7=all_file_ls[6][idx]
                file8=all_file_ls[7][idx]
                
                if file7.split('/')[-1] == 'NA' or file8.split('/')[-1] == 'NA' :
                    file7=None
                    file8=None

                df=main1(all_file_ls[0][idx], all_file_ls[1][idx], all_file_ls[2][idx], all_file_ls[3][idx],
                      all_file_ls[4][idx], all_file_ls[5][idx], file7, file8,
                      intermediate_option, weight_option,TST_option, temperature = temperature, with_spe=with_spe, energy_cut_off=energy_cut_off)
                
            elif intermediate_option==False:
                
                TST_option=False
                
                df=main1(all_file_ls[0][idx], all_file_ls[1][idx], all_file_ls[2][idx], all_file_ls[3][idx],
                      all_file_ls[4][idx], all_file_ls[5][idx], None, None,
                      intermediate_option, weight_option,TST_option, temperature = temperature, with_spe=with_spe, energy_cut_off=energy_cut_off)
                
            
            result_df_ls.append(df)

            
        self.raw_result_df=pd.concat(result_df_ls)
        
    
    def clean_result_df(self):
        
        ## clean the raw result df
        
        
        result_df1=self.raw_result_df.dropna()


        major_ls=result_df1['major_per'].tolist()
        major_ls1=[]
        for i in major_ls:
            if i> 100:
                major_ls1.append(100)
            elif i<0:
                major_ls1.append(0)
            else:
                major_ls1.append(i)
            
        minor_ls=result_df1['minor_per'].tolist()

        minor_ls1=[]
        for i in minor_ls:
            if i> 100:
                minor_ls1.append(100)
            elif i<0:
                minor_ls1.append(0)
            else:
                minor_ls1.append(i)
        
        result_df2=result_df1.copy()
        result_df2['major_perc']=major_ls1
        result_df2['minor_perc']=minor_ls1

        result_df3=result_df2[['TS1_name','int_name','major','minor','major_perc','minor_perc','TST_valid']]
        
        
        
        drop_df=result_df3[(result_df3['major_perc'] == 0) & (result_df3['minor_perc'] == 0)]

        df3 = result_df3.merge(drop_df, how='outer', indicator=True)
        df3 = df3.loc[df3['_merge'] == 'left_only']
        df3 = df3.drop(columns='_merge')    
            
        self.result_df=df3
        
        
        
        
    def cal_ratio(self):
        
        ## perform ratio calculations 
        
        self.clean_result_df()
        
        TS1_ls=self.result_df['TS1_name'].unique().tolist()
        
        int_ls=self.result_df['int_name'].unique().tolist()
        
        comb_ls=[]
        for ts1 in TS1_ls:
            for it in int_ls:
                comb_ls.append([ts1,it])
                
        
        processed_result_df_ls=[]
        for ls in comb_ls:
            #try:
            test_df=self.result_df[(self.result_df['TS1_name'] == ls[0])&(self.result_df['int_name'] == ls[1]) ]
            df = all_possible_ratio_df(test_df)
            df['TS1_name'] = ls[0]
            df['int_name'] = ls[1]
            processed_result_df_ls.append(df)
            #except Exception as e:
                #print(e)
    
        self.processed_result_df=pd.concat(processed_result_df_ls)
        
        
    def autopipe(self, path='', intermediate_option=True, weight_option=True, TST_option=True, getcsv=False, with_spe=False, temperature = 298.0,
                 energy_cut_off=5):
        
        self.perform_vrai(path=path, intermediate_option=intermediate_option, weight_option=weight_option, TST_option=TST_option, with_spe=with_spe,
                          temperature = temperature,energy_cut_off=energy_cut_off)
        
        self.cal_ratio()
        
        if getcsv==True: 
            
            today = date.today()
            date_str=today.strftime("%d%m%Y")
            
            csv_name=self.file_info.split('/')[-2]
            
            self.processed_result_df.to_csv(csv_name+'_'+date_str+'.csv')
        


            
        
        





















