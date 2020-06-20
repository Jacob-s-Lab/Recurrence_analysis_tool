#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:07:24 2020

@author: shihyu
"""


import pandas as pd
CLL_WES_list=["CGMH_CLL_014_WES","CGMH_CLL_026_WES","CGMH_CLL_030_WES","CGMH_CLL_031_WES","CGMH_CLL_038_WES","CGMH_CLL_040_WES","CGMH_CLL_042_WES","CGMH_CLL_213_WES","CGMH_CLL_226_WES","CGMH_CLL_244_WES","CGMH_CLL_02A_WES","CGMH_CLL_05A_WES","CGMH_CLL_01A_WES","CGMH_CLL_09A_WES"]

#read file and get header
FILE="ANNOVAR_CGMH_CLL_014_WES.popAF_annotation"
Annotation_file=pd.read_csv(f"{FILE}.txt",sep="\t")
anno_header=Annotation_file.columns.tolist()


#Generate info tags for N/T
N_info_tags=['GT', 'SQ', 'AD', 'AF', 'F1R2', 'F2R1', 'DP', 'SB', 'MB', 'PS']
T_info_tags=['GT', 'SQ', 'AD', 'AF', 'F1R2', 'F2R1', 'DP', 'SB', 'MB', 'PS']
for i in range(10):
    N_info_tags[i]+="_N"
    T_info_tags[i]+="_T"

#processing columns
CLL_WES_list=["CGMH_CLL_014_WES","CGMH_CLL_026_WES","CGMH_CLL_030_WES","CGMH_CLL_031_WES","CGMH_CLL_038_WES","CGMH_CLL_040_WES","CGMH_CLL_042_WES","CGMH_CLL_213_WES","CGMH_CLL_226_WES","CGMH_CLL_244_WES","CGMH_CLL_02A_WES","CGMH_CLL_05A_WES","CGMH_CLL_01A_WES","CGMH_CLL_09A_WES"]

for i in CLL_WES_list:
    file=pd.read_csv(f"{i}.popAF_annotation.txt",sep="\t")
    normal_info=file["Otherinfo13"]
    tumor_info=file["Otherinfo14"]
    n_list=[]
    for j in range(len(normal_info)):
        n=normal_info[j].split(":")
        n_list.append(n)   
     
    t_list=[]
    for k in range(len(tumor_info)):
        t=tumor_info[k].split(":")
        t_list.append(t)
        
    df1=pd.DataFrame(n_list, columns=N_info_tags)
    df2=pd.DataFrame(t_list, columns=T_info_tags)
    final_Annotation_file=pd.concat([file,df1,df2],axis=1)
    final_Annotation_file.to_excel(f"{i}.popAF_annotation_vcfinfo.xlsx", index=False)

#Change the format of union file
pd.read_csv("onlypopAF_union_of_CLL_WES_variants.txt",sep=" ", names=["Chr","Start","End","Ref","Alt","Function","Gene"]).to_excel("nofilter_union_of_CLL_WES_variants.xlsx", index=False)

CLL_WES_list=["CGMH_CLL_014_WES"]

#Calculate VAF
def calculate_AF_T_N(excel_file):
    val_list=[]
    df=pd.read_excel(excel_file)
    for i in range(len(df)):
        val=df.loc[i]["AF_T"]-df.loc[i]["AF_N"]
        val_list.append(val)
    valdf=pd.DataFrame(val_list,columns=["AF_T-AF_N"])
    finaldf=pd.concat([df,valdf],axis=1)
    return finaldf
    #finaldf.to_excel(f"AFfiltered_{excel_file}", index=False)
    
    
#create index  
def create_idx(df):
    idx_list=[]
    L=df
    L["Start"]=L["Start"].apply(str)
    L["End"]=L["End"].apply(str)
    
    for i in range(len(L)):
        idx=L.loc[i]["Chr"]+"/"+L.loc[i]["Start"]+"/"+L.loc[i]["End"]+"/"+L.loc[i]["Ref"]+"/"+L.loc[i]["Alt"]
        idx_list.append(idx)
    idx_df=pd.DataFrame(idx_list,columns=["Idx"])
    finaldf=pd.concat([L,idx_df], axis=1)
    return finaldf
    #finaldf.to_excel(f"Indexed_{excel_file}", index=False)


#generate datasets
for i in CLL_WES_list:
    create_idx(calculate_AF_T_N(f"{i}.popAF_annotation_vcfinfo.xlsx")).to_excel(f"TEST_Indexed_{i}.xlsx", index=False)
    
#process union dataset   
create_idx("nofilter_union_of_CLL_WES_variants.xlsx").to_excel("indexed_nofilter_union_of_CLL_WES_variants.xlsx", index=False)






#Generate array by Idx
U=pd.read_excel("Indexed_onlypopAF_union_of_CLL_WES_variants.xlsx",index_col="Idx")
for i in CLL_WES_list:
    S=pd.read_excel(f"AFfiltered_Indexed_{i}.popAF_annotation_vcfinfo.xlsx",index_col="Idx")
    affiltered_S=S[(S.AF_T  >= 0.05)&(S.AF_T-S.AF_N >= 0.05)]
    ary=pd.concat([U,affiltered_S[["Otherinfo13","Otherinfo14"]]],axis=1,join='outer')
    ary=ary.rename(columns={'Otherinfo14':f'{i}_T','Otherinfo13':f'{i}_N'})
    U=ary

#Count recurrence of samples
ary["Counts_new"]=ary.iloc[:,7:35:2].count(axis=1)
ary.to_excel("nofunctionfilter_Union_variants_array_VAF005_TNinfo.xlsx", index=False)



#Make list of unique Genes
ary=pd.read_excel("VariantBasedArray_popAF_VAF005.xlsx", index_col="Gene")
#ary.drop('Gene.refGene',axis=0,inplace=True) #remove title
idx_list=ary.index.unique().to_list()


#Make Gene-based recurrence array
df=pd.concat([pd.DataFrame(idx_list, columns=['Gene']),pd.DataFrame(columns=CLL_WES_list)],axis=1)
df.set_index("Gene", inplace=True)
for i in idx_list:
    for j in CLL_WES_list:
        n=0
        if isinstance(ary.loc[i][f"{j}_T"],float): #no data in ary[i][j]
            df.at[i,j]=0
        elif isinstance(ary.loc[i][f"{j}_T"],str): #only one position
            df.at[i,j]=1
        else:
            for k,v in ary.loc[i][f"{j}_T"].items(): #more than 2 positions
                if isinstance(v,str):
                    n+=1
            df.at[i,j]=n
    #calculate the occurrence of position for each gene
    if isinstance(ary.loc[i], pd.DataFrame):
        df.at[i,'POS_SUM']=len(ary.loc[i])
    else:
        df.at[i,'POS_SUM']=1
        
#calculate the sum of variants occurrence among samples for each gene
df["Occurrence(variants)"]=df.loc[:,CLL_WES_list].sum(axis=1)
#count the recurrence of samples
df["Counts"]=(df.loc[:,CLL_WES_list] != 0).astype(int).sum(axis=1)
#Output
df.to_excel("GeneBasedArray_CLLWES_popAF_VAF005.xlsx")


    