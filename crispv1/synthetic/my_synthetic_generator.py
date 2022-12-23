import numpy as np
import pandas as pd
from scipy.stats import randint, bernoulli

mu_frac_of_i=0.1
var_frac_of_i=0.001

#def synthetic_generator(n=100, d_layer=3, n_layer=[50, 100, 200], mu=0, sigma=1, n_causal=10, n_env=2):
def synthetic_generator(n=50, d_layer=3, n_layer=[5, 10, 20], mu=0, sigma=1, n_causal=10, n_env=2):

    assert d_layer == len(n_layer)

    max_weight = 1

    d = sum(n_layer)
    # determining childrens
    n_children = randint.rvs(0, 2 * n_causal, size=sum(n_layer[:-1]))
    children = []
    child_weight = []
    count = 0
    for i in range(d_layer - 1):
        count += n_layer[i]
        for j in range(n_layer[i]):
            temp = randint.rvs(count, d - 1, size=n_children[count + j - n_layer[i]])
            children.append(temp)
            child_weight.append(np.random.uniform(-max_weight, max_weight, n_children[count + j - n_layer[i]]))
    # generating features based on the SEM structure
    X = dict()
    X['value'] = np.zeros((n, d + 1))
    X['color'] = np.zeros((n, d + 1))
    for i in range(n):
        for j in range(n_layer[0]):  # Change order of for loops to fix wieghts across samples
            X['value'][i][j] = np.random.normal(mu, sigma)
            X['color'][i][j] = 0
            #X[i][j] = np.random.negative_binomial(n=100000, p=0.2)
            X['value'][i][children[j]] += 0.5 * X['value'][i][j] * child_weight[j]  # np.random.uniform(-2,2)
            X['color'][i][children[j]] = 1
            # print('After:', X[i][children[j]])
        count = n_layer[0]
        for k in range(1, d_layer - 1):
            for j in range(count, count + n_layer[k]):
                #X[i][j] += np.random.normal(mu, sigma)
                X['value'][i][j] += np.random.normal(mu_frac_of_i * j, mu_frac_of_i * j)
                X['color'][i][j] = 2
                X['value'][i][children[j]] += X['value'][i][j] * child_weight[j]  # *np.random.uniform(-2,2)
                X['color'][i][children[j]] = 3
            count += n_layer[k]
        for j in range(count, count + n_layer[-1]):
            X['value'][i][j] += np.random.normal(mu_frac_of_i * j , mu_frac_of_i * j)
            X['color'][i][j] = 4
            #X[i][j] += np.random.normal(mu, sigma)

    # generate target variables only from the top layer
    i_causal = randint.rvs(0, n_layer[0], size=n_causal)
    w_causal = randint.rvs(-max_weight, max_weight, size=n_causal)
    for i in range(n):
        x = np.dot(X['value'][i][i_causal], w_causal) + np.random.normal(mu, sigma)
        #x = np.dot(X[i][i_causal],w_causal) + np.random.negative_binomial(n=1000, p=0.2)
        y = 1 / (1 + np.exp(-x))
        X['value'][i][d] = 1 if y > 0.5 else 0
        X['color'][i][d] = 5


    # generating a dataframe column names
    ii = 0
    columns = []
    for i in range(d):
        if i in i_causal:
            columns.append("Causal_" + str(ii))
            ii += 1
        else:
            columns.append("Non_causal_" + str(i - ii))
    columns.append("Target")
    df = pd.DataFrame(data=X['value'], columns=columns)

    # adding the subject id and environment splits to the dataframe
    ID = 1
    subj_id = [ID]
    env = bernoulli.rvs(0.6)
    # env = [i for i in range(len(n_env))]
    print('env is ', env)
    env_split = [env]
    #env_split = []
    for i in range(1, n):
        if X['value'][i][d] != X['value'][i - 1][d]:
            ID += 1
            env = bernoulli.rvs(0.6)
            #env = np.random.uniform(0, n_env)
            print('env = ', env)
        subj_id.append(ID)
        env_split.append(env)
        print('i: ', str(i), 'color: ', X['color'][i][d])
    df["Subj_ID"] = subj_id
    df["env_split"] = env_split

    #     print(df.shape)

    # display mean variance graph
    return df, X['color']
