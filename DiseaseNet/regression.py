import time
import pandas as pd
import numpy as np
import statsmodels.discrete.conditional_models as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

def defination(sick:str, history:list, inpatient:list) -> int:
    """
    """
    if float(sick) in history or float(sick) in inpatient:
        return 1
    else:
        return 0

def unconditional_logistic(data:pd.DataFrame, result:pd.DataFrame, n_cpus:int, adjustment:str, coe:float) -> pd.DataFrame:
    """
    """
    from sklearn.linear_model import LogisticRegression as LR

    def create_covariates(data:pd.DataFrame, covariates:list[str]) -> pd.DataFrame:
        """
        """
        temp_variates = []
        for variate in covariates:
            temp_df = pd.get_dummies(data[variate], prefix=variate)
            temp_variates += [x for x in temp_df.columns[1::]]
            data = pd.concat([data, temp_variates], axis=1)
        return data

    def logistic_unconditional(d1d2:list, df_matched_group:pd.DataFrame, illness_list:list) -> list:
        time1 = time.time()
        dataset = df_matched_group
        inpatient_variable = 'inpatient_level1'
        eligible_variable = 'd_eligible'
        temp = []
        
        d1 = float(d1d2.split('-')[0])
        d2 = float(d1d2.split('-')[1])
        
        delete_disease_list = illness_list.copy()

        try:
            delete_disease_list.remove(str(d1))
        except:
            delete_disease_list = delete_disease_list
        try:
            delete_disease_list.remove(str(d2))
        except:
            delete_disease_list = delete_disease_list
            
        dataset['flag'] = dataset[eligible_variable].apply(lambda x: d1 in x and d2 in x)
        dataset_d = dataset.loc[dataset['flag']==True]

        # 删除 样本容量中协变量方差等于0的
        var_co_vars_lst_d = dataset_d[covariates].var()
        co_vars_selected = var_co_vars_lst_d[var_co_vars_lst_d != 0].index
        var_covar_lst_d = dataset_d[delDiseaseList].var()
        delDiseaseList = var_covar_lst_d[var_covar_lst_d != 0].index

        dataset_d['d1'] = dataset_d[inpatient_variable].apply(lambda x: 1 if d1 in x.keys() else 0)    
        dataset_d['d2'] = dataset_d[inpatient_variable].apply(lambda x: 1 if d2 in x.keys() else 0)
        dataset_d['constant'] = 1

        # len_d_other, len_cov = len(covar_lst_all_d), 2 + len(co_vars_selected)
        time3 = time.time()
        print("before sklearn model cost %.4f seconds" % (time3-time1))
        # sklearn model
        model = LR(solver='saga', penalty='l1', random_state=42, C=coe**(-1), n_jobs=-1)
        model.fit(np.asarray(dataset_d[['d2','constant'] + list(co_vars_selected) + list(delDiseaseList)], dtype=int),
                            np.asarray(dataset_d['d1'], dtype=int))
        ii = np.flatnonzero(model.coef_)
        ii_ = ii[np.where(ii>(np.asarray(dataset_d[['d2','constant'] + list(co_vars_selected)]).shape[1]-1))]
        covar_lst_all_d1 = dataset_d[['d2','constant'] + list(co_vars_selected) + list(delDiseaseList)].iloc[:,ii_].columns

        #1 d1 L1正则项
        model1 = sm.Logit(np.asarray(dataset_d['d1'],dtype=int),
                                    np.asarray(dataset_d[['d2','constant'] + list(co_vars_selected) + list(covar_lst_all_d1)], dtype=int))
        try:
            try:
                result = model1.fit(maxiter=300)
                temp += [d1d2, result.params[0], result.pvalues[0],'%.2f (%.2f-%.2f)' % (np.exp(result.params[0]),
                                                                                        np.exp(result.conf_int()[0][0]),
                                                                                        np.exp(result.conf_int()[0][1])), np.NaN, 'newton: %i' % (len(list(covar_lst_all_d1)))]
            except:
                result = model1.fit(method='cg', maxiter=300)
                temp += [d1d2, result.params[0], result.pvalues[0],'%.2f (%.2f-%.2f)' % (np.exp(result.params[0]),
                                                                                        np.exp(result.conf_int()[0][0]),
                                                                                        np.exp(result.conf_int()[0][1])), np.NaN, 'cg: %i' % (len(list(covar_lst_all_d1)))]
        except Exception as e:
            print(e)
            temp += [d1d2,np.NaN,np.NaN,np.NaN,e, len(list(covar_lst_all_d1))]

        time2 = time.time()
        print('Spent %0.2f s time in sklearn model and statesmodel' % (time2 - time3))
        return temp

    def main():
        """
        """
        pass


def conditional_logistic(data:pd.DataFrame, result:pd.DataFrame, n_cpus:int, adjustment:str, coe:float) -> pd.DataFrame:
    """
    """
    
    def match(data:pd.DataFrame, time_variate:float) -> pd.DataFrame:
        """
        """
        variate_1 = 'sex'
        variate_2 = 'birth_date'
        date_death = 'time_end'
        match_number = 2
        print('Start matching')
        one_year = 365.25*24*3600
        case = data.loc[~data[time_variate].isna()]
        result = []
        iter_ = 0
        for j in case.index:
            time_ = case.loc[j, time_variate]
            try:
                t = data.loc[(data[variate_1]==case.loc[j,variate_1]) & 
                                (np.abs(data[variate_2]-case.loc[j,variate_2])<=one_year) & 
                                ~(data[date_death]<=time_) & 
                                ~(data[time_variate]<=time_)].sample(match_number,random_state=0)
            except:
                t = data.loc[(data[variate_1]==case.loc[j,variate_1]) & 
                                (np.abs(data[variate_2]-case.loc[j,variate_2])<=one_year) & 
                                ~(data[date_death]<=time_) & 
                                ~(data[time_variate]<=time_)]
                
            result.append([case.loc[j,'eid'], 1, time_,iter_])
            result += [[x, 0, time_, iter_] for x in t.eid.values]
            iter_ += 1
        result = pd.DataFrame(result,columns=['eid','outcome_conlo','date_start_conlo','match_2_conlo'])
        result = pd.merge(result, data, on='eid', how='left')
        return result

    def logistic_conditional(d1d2:list, df_matched_group:pd.DataFrame, illList:list) -> list:
        """
        """

        dataset = df_matched_group.copy()
        inpatient_variable = 'inpatient_level1'
        eligible_variable = 'd_eligible'
        d1 = float(d1d2.split('-')[0])
        d2 = float(d1d2.split('-')[1])
        
        delDiseaseList = illList.copy()
        delDiseaseList.remove(str(d1)), delDiseaseList.remove(str(d2))

        dataset['flag'] = dataset[eligible_variable].apply(lambda x: d1 in x and d2 in x)
        dataset_d = dataset.loc[dataset['flag']==True]
        
        dataset_d['d2_time'] = dataset_d[inpatient_variable].apply(lambda x: x.get(d2, np.NaN))
        dataset_d['d1_time'] = dataset_d[inpatient_variable].apply(lambda x: x.get(d1, np.NaN))
        dataset_d_matched = match(dataset_d, 'd2_time')
        dataset_d_matched['exposure'] = dataset_d_matched.apply(lambda row: 1 if 
                                                                row['d1_time'] < row['date_start_conlo'] else 0,
                                                                axis=1)
        # 排除在匹配样本中 co_var变量方差为0
        var_co_vars_lst_d = dataset_d_matched[co_vars].var()
        co_vars_selected = var_co_vars_lst_d[var_co_vars_lst_d != 0].index
        var_covar_lst_d = dataset_d_matched[delDiseaseList].var()
        delDiseaseList = var_covar_lst_d[var_covar_lst_d != 0].index
        
        # 第一次L1正则项拟合
        model = sm.ConditionalLogit(np.asarray(dataset_d_matched['outcome_conlo'], dtype=int),
                            np.asarray(dataset_d_matched[['exposure'] + list(co_vars_selected) + list(delDiseaseList)],dtype=int),
                            groups=dataset_d_matched['match_2_conlo'].values)
        len_d_other, len_cov = len(delDiseaseList), 1 + len(co_vars_selected)
        result = model.fit_regularized(method='elastic_net',alpha=[0]*len_cov+[coe]*len_d_other)
        
        # 选出非零系数的疾病协变量
        ii = np.flatnonzero(result.params)
        covar_lst_all = dataset_d_matched[['exposure'] + list(co_vars_selected) + list(delDiseaseList)].iloc[:,ii].columns
        covar_lst_all_d = [x for x in covar_lst_all if x in list(delDiseaseList)]
        
        # 删除 L1正则项后的 vif> =5的疾病协变量
        new_dataset = dataset_d_matched[covar_lst_all_d]
        new_dataset['constant'] = 1
        covar_lst_all_d1 = []
        for i in range(len(new_dataset.columns)-1):
            try:
                vif = variance_inflation_factor(new_dataset.values, i)
                if vif < 5:
                    covar_lst_all_d1.append(new_dataset.columns[i])
            except:
                covar_lst_all_d1.append(new_dataset.columns[i])
            
        model_refit = sm.ConditionalLogit(np.asarray(dataset_d_matched['outcome_conlo'],dtype=int),
                                        np.asarray(dataset_d_matched[['exposure'] + list(co_vars_selected) + list(covar_lst_all_d1)],dtype=int),
                                        groups=dataset_d_matched['match_2_conlo'].values)
        try:
            try:
                result_refit = model_refit.fit(maxiter=300)
                return [d1d2,result_refit.params[0],result_refit.pvalues[0],'%.2f (%.2f-%.2f)' % (np.exp(result_refit.params[0]),
                                                            np.exp(result_refit.conf_int()[0][0]),
                                                            np.exp(result_refit.conf_int()[0][1])),np.NaN, 'BFGS: %i' % (len(list(covar_lst_all_d1)))]
            except:
                result_refit = model_refit.fit(method='cg', maxiter=300)
                return [d1d2,result_refit.params[0],result_refit.pvalues[0],'%.2f (%.2f-%.2f)' % (np.exp(result_refit.params[0]),
                                                            np.exp(result_refit.conf_int()[0][0]),
                                                            np.exp(result_refit.conf_int()[0][1])),np.NaN, 'cg: %i' % (len(list(covar_lst_all_d1)))]
        except Exception as e:
            print(e)
            return [d1d2,np.NaN,np.NaN,np.NaN,e,len(list(covar_lst_all_d1))]

    def main():
        pass
