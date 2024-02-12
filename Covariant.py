import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
import sys
sys.path.insert(1, 'C:/Users/Hasan/HKUST/Haibin SU - group - covid19/SingleSiteAnalysis')
import ImportantFunc as Imp
import time
from datetime import datetime
start_time = time.process_time()

class Covariant:

    def __init__(self,mutList):
        self.MutationList = mutList
        self.MutPair = []
        self.MutSingle = []
        self.dfFinalMutPair = pd.DataFrame({})

    def Unique_pair_func(self,mutList):
        return [ mutList[i]+"|"+mutList[j] for i in range(len(mutList)) for j in range(len(mutList)) if i < j  ]

    def NodeCount(self,dfCountMut):
        Mut = Imp.count_dups(dfCountMut['x'].tolist() + dfCountMut['y'].tolist())
        dfNode = pd.DataFrame({'Mutation':Mut[0],'MutFreq':Mut[1]})
        dfNode.sort_values(['MutFreq'],inplace = True, ascending=[False])
        MutCount = Imp.count_dups(dfNode['MutFreq'].tolist())
        dfP_k = pd.DataFrame({'k':MutCount[0],'P(k)':MutCount[1]})
        MutConcat = Imp.concat_dups(dfNode['MutFreq'].tolist(),dfNode['Mutation'].tolist())
        dfMutConcat = pd.DataFrame({'k':MutConcat[0],'Mutations':MutConcat[1]})
        dfP_k = pd.merge(dfP_k,dfMutConcat,on='k',how='left')
        dfP_k.sort_values(['k'],inplace = True, ascending=[False])
        return dfP_k

    def CovariantMatrix(self):
        ## Collecting the mutation information (including the pairs) in one month
        for mutList in self.MutationList:
            mutList = mutList.split(";");
            self.MutSingle.extend(mutList)
            self.MutPair.extend( self.Unique_pair_func(mutList) )
        ## Creating DataFrame for Counting Mutation Pairs
        M = Imp.count_dups(self.MutPair);
        dfCountPair = pd.DataFrame({'Mut-Pair': M[0],'CountPair': M[1]})
        ## Creating DataFrame for Single Mutation Count
        n = Imp.count_dups(self.MutSingle);
        dfSingleX = pd.DataFrame({'x': n[0],'CountX': n[1]})
        dfSingleY = pd.DataFrame({'y': n[0],'CountY': n[1]})
        ## Retrieving Mutation Index in 1 Month
        MutEle = [ele for Mut in self.MutationList for ele in Mut.split(";")]
        M = Imp.count_dups(sorted(MutEle))
        df = pd.DataFrame({'Mutation': M[0],'Count': M[1]})
        df['Pos'] = df['Mutation'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
        df['Count'].astype(int), df['Pos'].astype(int)
        df.sort_values(['Pos','Count'],inplace = True, ascending=[True, False])
        df = df.reset_index(drop=True)
        MutIndex = df['Mutation'].tolist()
        ## Computing Covariance Matrix and Nullifying self-correlation mutation
        dfMutPair = pd.DataFrame({'Mut-Pair':self.Unique_pair_func(MutIndex)})
        dfMutPair[['x', 'y']] = dfMutPair['Mut-Pair'].str.split('|', n=1, expand=True)
        dfMutPair = pd.merge(dfMutPair,dfCountPair,how='left',on='Mut-Pair')
        dfMutPair = pd.merge(dfMutPair,dfSingleX,how='left',on='x')
        dfMutPair = pd.merge(dfMutPair,dfSingleY,how='left',on='y').fillna(0)
        dfMutPair['Covariant'] = (dfMutPair['CountPair']/len(self.MutationList)) - (dfMutPair['CountX']/len(self.MutationList))*(dfMutPair['CountY']/len(self.MutationList))
        dfMutPair['PosX'] = dfMutPair['x'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
        dfMutPair['PosY'] = dfMutPair['y'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
        dfMutPair['Marker'] = dfMutPair.apply(lambda x: 1 if x['PosX'] == x['PosY'] else 0, axis = 1)
        ## Nulifying the self-correlation sites
        dfMutPair.loc[dfMutPair['Marker'] == 1, 'Covariant'] = 0
        dfMutPair1 = dfMutPair[['x','y','Covariant']].copy(deep=True)
        dfMutPair1['Mut-Pair'] = dfMutPair['x']+"|"+dfMutPair['y']; dfMutPair1.drop(['x','y'],inplace=True,axis=1)
        dfMutPair2 = dfMutPair[['y','x','Covariant']].copy(deep=True)
        dfMutPair2['Mut-Pair'] = dfMutPair['y']+"|"+dfMutPair['x']; dfMutPair2.drop(['x','y'],inplace=True,axis=1)
        dfMutPair = pd.concat([dfMutPair1,dfMutPair2],axis = 0).reset_index(drop=True)
        dfMutPair = dfMutPair[['Mut-Pair','Covariant']]
        self.dfFinalMutPair = dfMutPair.copy(deep=True)
        self.dfFinalMutPair[['x', 'y']] = self.dfFinalMutPair['Mut-Pair'].str.split('|', n=1, expand=True)
        self.dfFinalMutPair = self.dfFinalMutPair.drop(['Mut-Pair'], axis=1)
        self.dfFinalMutPair = self.dfFinalMutPair.loc[(self.dfFinalMutPair['Covariant'] > 0)].reset_index(drop=True)
        print(self.dfFinalMutPair)

dfInp = pd.read_excel('Covariant_Input.xlsx')
print("Current time =", datetime.now().strftime("%H:%M:%S"))
print('TIME TAKEN: ' + str(time.process_time() - start_time) + 's\n')
cov = Covariant(dfInp['mutation info'].tolist())
cov.CovariantMatrix()
cov.dfFinalMutPair.to_excel('Covariant_Output.xlsx')
print("Current time =", datetime.now().strftime("%H:%M:%S"))