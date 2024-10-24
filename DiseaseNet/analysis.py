import pandas as pd
import numpy as np
import multiprocessing
import time

def phewas(data:pd.DataFrame, n_cpus:int, adjustment:str) -> pd.DataFrame:
    """
    """
    from lifelines import CoxPHFitter
    from statsmodels.duration.hazard_regression import PHReg
    cph = CoxPHFitter()
    phecode_cate = pd.read_csv(r'C:/Users/Administrator/Desktop/data/phecode_1.2/phecode_info.csv')
    phecode_lst = np.array([x for x in phecode_cate.phecode.values])
    result_final = []

    def range_d(x):
        """
        """
        x_str = str(x)
        if x_str.split('.')[1] == '0':
            return round(x+0.99,2)
        elif len(x_str.split('.')[1]) == 1:
            return round(x+0.09,2)
        else:
            return x
        
    def cox(disease:str, data:pd.DataFrame, threshold=1000, covariates=['age', 'sex', 'BMI']) -> float:
        """
        """
        result = [disease]
        dataset_analysis = data[covariates+'eid'+'date_start'+'date_end'+'inpatient'+'outcome'+'match']
        exl_range = phecode_cate.loc[phecode_cate['phecode'] == disease]['phecode_exclude_range'].values[0]
        if pd.isna(exl_range):
            dataset_analysis['flag_exl'] = 0
        else:
            exl_list = []
            for range_ in exl_range.split(','):
                exl_lower,exl_higher = float(range_.split('-')[0]), float(range_.split('-')[1])
                exl_list_index = np.where(np.all([phecode_lst >= exl_lower, phecode_lst <= exl_higher], axis=0))[0]
                exl_list_temp = phecode_lst[exl_list_index]
                exl_list += [x for x in exl_list_temp]
            dataset_analysis['flag_exl'] = dataset_analysis['history'].apply(lambda x: 
                                                                            1 if np.any([i in x for i in exl_list]) 
                                                                            else 0)
                
        dataset_analysis['flag_exl'] = dataset_analysis.apply(lambda row: 1 if disease in row['sex'] 
                                                            else row['flag_exl'],axis=1)

        dataset_analysis = dataset_analysis.loc[dataset_analysis['flag_exl']==0]
        if len(dataset_analysis) == 0:
            result += ['Sex specific']
            return result
        
        disease_upper = range_d(disease)
        disease_lst = [x for x in phecode_lst if x>=disease and x<=disease_upper]
        dataset_analysis['d_time'] = dataset_analysis['inpatient'].apply(lambda x: 
                                                                                np.nanmin([x.get(j,np.NaN) for j in disease_lst]))
            
        dataset_analysis['outcome_'] = dataset_analysis['d_time'].apply(lambda x: 0 if pd.isna(x) else 1)
        length = len(dataset_analysis.loc[(dataset_analysis['outcome']==1) & 
                                                (dataset_analysis['outcome_']==1)])

        result += [length]
        #
        dataset_analysis['date_final'] = np.nanmin(dataset_analysis[['d_time',
                                                                    'date_end']],axis=1)
        dataset_analysis['time'] = (dataset_analysis['date_final'] - 
                                    dataset_analysis['date_start'])/(365.25*24*3600)
        dataset_analysis['time'] = dataset_analysis['time'].astype(float)
        
        #time at risk
        n_exp = len(dataset_analysis.loc[(dataset_analysis['outcome']==1) & (dataset_analysis['outcome_']==1)])
        n_unexp = len(dataset_analysis.loc[(dataset_analysis['outcome']==0) & (dataset_analysis['outcome_']==1)])
        time_exp = dataset_analysis.groupby(by='outcome').sum()['time'].loc[1]/1000
        time_unexp = dataset_analysis.groupby(by='outcome').sum()['time'].loc[0]/1000
        
        if length < threshold:
            result += ['less than threshold','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                    '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
            return result

        match_id = dataset_analysis[dataset_analysis['outcome_']==1]['match']
        dataset_analysis = dataset_analysis[dataset_analysis['match'].isin(match_id)]
        try:
            model = PHReg(np.asarray(dataset_analysis['time'],dtype=np.float32),
                        np.asarray(dataset_analysis[['outcome']+covariates],dtype=np.float32),
                        status=np.asarray(dataset_analysis['outcome_'],dtype=np.int32), 
                        strata=np.asarray(dataset_analysis['match'],dtype=np.float32))
            model_result = model.fit(method='bfgs',maxiter=300,disp=1)
            if pd.isna(model_result.params[0]) or pd.isna(model_result.bse[0]):
                model = cph.fit(dataset_analysis[['time','outcome_','outcome','match']+covariates],
                                fit_options=dict(step_size=0.2), duration_col='time', event_col='outcome_',strata=['match'])
                result_temp = model.summary.loc['outcome']
                result += ['fitted_lifelines','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                            '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
                result += [x for x in result_temp[['coef','se(coef)','p']]]
            else:
                result += ['fitted','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                                    '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
                result += [model_result.params[0],model_result.bse[0],model_result.pvalues[0]]
        except:
            try:
                model = cph.fit(dataset_analysis[['time','outcome_','outcome','match']+covariates],
                                fit_options=dict(step_size=0.2), duration_col='time', event_col='outcome_',strata=['match'])
                result_temp = model.summary.loc['outcome']
                result += ['fitted_lifelines','%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                            '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
                result += [x for x in result_temp[['coef','se(coef)','p']]]
            except Exception as e:
                print(e)
                result += [e,'%i/%.2f (%.2f)' % (n_exp,time_exp,n_exp/time_exp),
                            '%i/%.2f (%.2f)' % (n_unexp,time_unexp,n_unexp/time_unexp)]
        return result
    
    def main():
        for i in range(phecode_lst):
            result_final.append(cox(phecode_lst[i], data))
    pool = multiprocessing.Pool(n_cpus)
    pool.map(main())
    return pd.DataFrame(result_final, columns=['disease','number','describe','exp','unexp','coef','se','p'])

def comorbidity_analysis(data:pd.DataFrame, n_cpus:int, adjustment:str) -> pd.DataFrame:
    """
    """
    threshold = 1
    d1d2_lst = []
    result = []
    d1_list = list(data['d1'])
    d2_list = list(data['d2'])
    for idx in range(len(d2_list)):
        d1d2_lst.append([d1_list[idx],d2_list[idx]])

    def comorbidity(d1d2_lst:list) -> list:
        for d1,d2 in d1d2_lst:
            data['flag'] = data['d_eligible'].apply(lambda x: d1 in x and d2 in x)
            df_ = data.loc[data['flag']==True]
            n = len(df_)
            c = sum([d1 in x and d2 in x for x in df_['inpatient_level1'].values]) #d1d2
            if c<=threshold:
                result.append([d1,d2,c,n,np.NaN,np.NaN,np.NaN,np.NaN,np.NaN])
                continue
            p1 = sum([d1 in x for x in df_['inpatient_level1'].values])
            p2 = sum([d2 in x for x in df_['inpatient_level1'].values])
            rr = (n*c)/(p1*p2) #
            theta = (1/c + 1/((p1*p2)/n) - 1/n - 1/n)**0.5 #RR
            t_ = abs(np.log(rr)/theta) #RRt
            p = (1-t.cdf(t_,n))*2 #RRP
            phi = (c*n-p1*p2)/(((p1*p2)*(n-p1)*(n-p2))**0.5) #phi
            if phi == 1:
                phi = 0.99
            z_phi = 0.5*np.log((1+phi)/(1-phi))
            z_phi_theta = (1/(n-3))**0.5
            z_phi_t = abs(z_phi/z_phi_theta)
            p_phi = (1-t.cdf(z_phi_t,n))*2
            result.append([d1,d2,c,n,rr,theta,p,phi,p_phi])
    pool = multiprocessing.Pool(n_cpus)
    pool.map(comorbidity())


def trajectory_analysis(data:pd.DataFrame, n_cpus:int, adjustment:str) -> pd.DataFrame:
    """
    """
    pool = multiprocessing.Pool(n_cpus)

    def exc_lst(disease:str, phecode_cate:pd.DataFrame, phecode_lst:list) -> set:
        """
        """
        lst = []
        exl_range = phecode_cate.loc[phecode_cate['phecode']==disease]['phecode_exclude_range'].values[0]
        if pd.isna(exl_range):
            return set(lst)
        else:
            for range_ in exl_range.split(','):
                exl_lower, exl_higher = float(range_.split('-')[0]), float(range_.split('-')[1])
                exl_list_index = np.where(np.all([phecode_lst>=exl_lower, phecode_lst<=exl_higher], axis=0))[0]
                exl_list_temp = phecode_lst[exl_list_index]
                lst += [x for x in exl_list_temp]
            return set(lst)
        
    def d1_d2(data:pd.DataFrame) -> pd.DataFrame:
        """
        """
        exposed = data.copy()
        id_ = 'eid'
        inpatient_variable = 'inpatient_level1'
        date_start_variable = 'dia_date'
        eligible_variable = 'd_eligible'

        array = exposed[[id_, date_start_variable, inpatient_variable, eligible_variable]].values
        d1d2_result = []
        total = len(array)
        time0 = time.time()
        for i in range(len(array)):
            if i%3000 == 0:
                print("Progress: %.1f%%" % (i/total*100))
                print("Time spent: %.1f mins" % ((time.time()-time0)/60))
            d1d2 = []
            dict_temp = array[i][2]
            eligible = array[i][-1]
            d_list = [x for x in dict_temp.keys() if x in eligible]
            length = len(d_list)
            if length <= 1:
                d1d2_result.append([])
                continue
            for j in range(length-1):
                for k in (range(j+1,length)):
                    d1 = d_list[j]
                    d2 = d_list[k]
                    if dict_temp.get(d2) > dict_temp.get(d1):
                        d1d2.append("%s-%s" % (d1,d2))
            d1d2_result.append(d1d2)
        exposed['d1d2'] = d1d2_result
        return exposed

    def main():
        pass