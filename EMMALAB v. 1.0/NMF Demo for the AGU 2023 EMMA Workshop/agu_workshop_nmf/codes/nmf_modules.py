# load Packages
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from concurrent.futures import ThreadPoolExecutor, as_completed

# settings
pd.options.mode.chained_assignment = None

def run_single_model(data_norm_path, data_synthetic_path):
    data_norm = pd.read_csv(data_norm_path, header=0, index_col=None)
    data_synthetic = pd.read_csv(data_synthetic_path, header=0, index_col=None)

    model = NMF(n_components=4, init='random',
                solver = 'cd', random_state=233,
                max_iter = 10000,  tol = 1e-6)
    mod = model.fit(data_synthetic)
    W = mod.transform(data_norm) # W and K are the proportion matrix
    K = mod.transform(data_synthetic)
    H = mod.components_ # H is the chemical composition of endmember
    plt.figure()
    sns.heatmap(H, annot=True, cmap="YlGnBu",xticklabels=list(data_norm.columns.values))
    plt.show()
    print("The number of model iteration:", model.n_iter_)
    print("The dimension of matrix W is", W.shape)
    print("The dimension of matrix W is", H.shape)
    #np.savetxt('../output/datasets/'+'example.csv',H,delimiter=",")

# define a nmf model with four components and with All, Cl, SO4, K, Mg, Ca all normalized by Na
def run_model(data_norm, data_synthetic, start=None):
    try:
        # run model
        model = NMF(n_components=4, init='random',
                    solver='cd', random_state=start + 1,
                    max_iter=100000, tol=1e-4)
        # fit synthetic data
        mod = model.fit(data_synthetic)
        # apply model to real data
        W = mod.transform(data_norm)
        H = mod.components_
        data_col_names = list(data_norm.columns.values)
        result = pd.DataFrame(W, columns=['comp1', 'comp2', 'comp3', 'comp4'])
        result['tot'] = result['comp1'] + result['comp2'] + \
                        result['comp3'] + result['comp4']
        # keep samples that add to 1 +/- 0.10
        result = result[result['tot'] > 0.90]
        result = result[result['tot'] < 1.10]
        result['sample'] = result.index + 1
        seed = start + 1

        # define end-members
        em = pd.DataFrame(H, columns=data_col_names)

        # set up rules for delineating end members
        index_cl = np.argmax(np.array(em['Cl_Na']))
        index_so4 = np.argmax(np.array(em['SO4_Na']))
        index_alk = np.argmax(np.array(em['Alk_Na']))
        index4 = 6 - (index_cl + index_so4 + index_alk)

        # end member assignment
        result['em1'] = result.iloc[:, [index_cl]]
        result['em2'] = result.iloc[:, [index_so4]]
        result['em3'] = result.iloc[:, [index_alk]]
        result['em4'] = result.iloc[:, [index4]]
        cols = ['sample', 'tot', 'em1', 'em2', 'em3', 'em4']
        result = result[cols]

        # calculated modeled data for backend filtering
        result['Ca_mod'] = result['em1'] * em.loc[index_cl, 'Ca_Na'] + \
                           result['em2'] * em.loc[index_so4, 'Ca_Na'] + \
                           result['em3'] * em.loc[index_alk, 'Ca_Na'] + \
                           result['em4'] * em.loc[index4, 'Ca_Na']

        result['Mg_mod'] = result['em1'] * em.loc[index_cl, 'Mg_Na'] + \
                           result['em2'] * em.loc[index_so4, 'Mg_Na'] + \
                           result['em3'] * em.loc[index_alk, 'Mg_Na'] + \
                           result['em4'] * em.loc[index4, 'Mg_Na']

        result['K_mod'] = result['em1'] * em.loc[index_cl, 'K_Na'] + \
                          result['em2'] * em.loc[index_so4, 'K_Na'] + \
                          result['em3'] * em.loc[index_alk, 'K_Na'] + \
                          result['em4'] * em.loc[index4, 'K_Na']

        result['SO4_mod'] = result['em1'] * em.loc[index_cl, 'SO4_Na'] + \
                            result['em2'] * em.loc[index_so4, 'SO4_Na'] + \
                            result['em3'] * em.loc[index_alk, 'SO4_Na'] + \
                            result['em4'] * em.loc[index4, 'SO4_Na']

        result['Cl_mod'] = result['em1'] * em.loc[index_cl, 'Cl_Na'] + \
                           result['em2'] * em.loc[index_so4, 'Cl_Na'] + \
                           result['em3'] * em.loc[index_alk, 'Cl_Na'] + \
                           result['em4'] * em.loc[index4, 'Cl_Na']

        result['Alk_mod'] = result['em1'] * em.loc[index_cl, 'Alk_Na'] + \
                           result['em2'] * em.loc[index_so4, 'Alk_Na'] + \
                           result['em3'] * em.loc[index_alk, 'Alk_Na'] + \
                           result['em4'] * em.loc[index4, 'Alk_Na']


        # save end member composition
        result['Ca_Na_em1'] = em.loc[index_cl, 'Ca_Na']
        result['Mg_Na_em1'] = em.loc[index_cl, 'Mg_Na']
        result['K_Na_em1'] = em.loc[index_cl, 'K_Na']
        result['SO4_Na_em1'] = em.loc[index_cl, 'SO4_Na']
        result['Cl_Na_em1'] = em.loc[index_cl, 'Cl_Na']
        result['Alk_Na_em1'] = em.loc[index_cl, 'Alk_Na']

        result['Ca_Na_em2'] = em.loc[index_so4, 'Ca_Na']
        result['Mg_Na_em2'] = em.loc[index_so4, 'Mg_Na']
        result['K_Na_em2'] = em.loc[index_so4, 'K_Na']
        result['SO4_Na_em2'] = em.loc[index_so4, 'SO4_Na']
        result['Cl_Na_em2'] = em.loc[index_so4, 'Cl_Na']
        result['Alk_Na_em2'] = em.loc[index_so4, 'Alk_Na']

        result['Ca_Na_em3'] = em.loc[index_alk, 'Ca_Na']
        result['Mg_Na_em3'] = em.loc[index_alk, 'Mg_Na']
        result['K_Na_em3'] = em.loc[index_alk, 'K_Na']
        result['SO4_Na_em3'] = em.loc[index_alk, 'SO4_Na']
        result['Cl_Na_em3'] = em.loc[index_alk, 'Cl_Na']
        result['Alk_Na_em3'] = em.loc[index_alk, 'Alk_Na']

        result['Ca_Na_em4'] = em.loc[index4, 'Ca_Na']
        result['Mg_Na_em4'] = em.loc[index4, 'Mg_Na']
        result['K_Na_em4'] = em.loc[index4, 'K_Na']
        result['SO4_Na_em4'] = em.loc[index4, 'SO4_Na']
        result['Cl_Na_em4'] = em.loc[index4, 'Cl_Na']
        result['Alk_Na_em4'] = em.loc[index4, 'Alk_Na']

        # save results
        return result
    except:
        # error if none of the delineation rules are followed
        print("error in seed: {}".format(seed))


def main_run_model(suffix, data_norm_path, data_synthetic_path):
    # import data
    data_norm = pd.read_csv(data_norm_path, header=0, index_col=None)
    data_synthetic = pd.read_csv(data_synthetic_path, header=0, index_col=None)

    # define total model dataframe
    final = pd.DataFrame()
    # set up parallel computing.
    pool = ThreadPoolExecutor(max_workers=60)
    futures = []
    # Original times of iterate model with random initiation in Shaughnessy et al is 10,000
    # Modified to 3000 to save time
    for x in range(3000):
        futures.append(pool.submit(run_model, data_norm, data_synthetic, x))
    pool.shutdown()
    counter = 1
    # collect model results as they complete
    for x in as_completed(futures):
        final = pd.concat([final,x.result()])
        #print('model run: {}'.format(counter))
        counter = counter + 1
    # write out model results
    final.to_csv('../output/datasets/' + 'nmf_out_4em_Na_' + suffix + '.csv', index=False)
    print('Done.')
