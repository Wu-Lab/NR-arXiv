import numpy as np
from scipy.linalg import null_space
from numpy.linalg import inv

'''
Parameters Explanation:
mat: input noisy matrix
m: Parameters of NR_F, controlling the extent of matrix modification
eps1: Value added to the whole matrix in preprocessing step
eps2: Value added to the diagonal in preprocessing step
'''

def NR_B(mat, m, eps1, eps2):
    
    # preprocessing
    n = mat.shape[0]
    mat = mat + eps1 + eps2 * np.eye(n)
    
    # change into transition matrix
    P1 = mat / mat.sum(1).reshape(-1, 1)
    P2 = m * P1.dot(inv(((m-1)*np.eye(n) + P1)))
  
    # minus the smallest negative value from each row with negative numbers
    process = lambda x:min(x,0) 
    row_fac = np.zeros(n)
    for i in range(n):
        row_fac[i] = process(P2.min(1)[i])
    P2 = P2 - row_fac.reshape(-1,1)  
    P2 = P2 / P2.sum(1).reshape(-1, 1)
    
    # compute the stationary distribution and the output matrix
    stationary_d = null_space((P2 - np.eye(n)).T)
    if stationary_d[0].shape[0]==0:
        print('The stationary distribution does not exist！')
    elif stationary_d[0].shape[0]>1:
        print('The stationary distribution is not unique！')
    else :
        net_new = np.diag(abs(stationary_d.T)[0]).dot(P2)
        net_new = net_new + net_new.T
        output_network = net_new-np.diag(np.diag(net_new))
        output_network = output_network / output_network.max()
        return (output_network)
