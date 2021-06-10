"""
Parameters Explanation:
mat: input noisy matrix
m: Parameters of NR_F, controlling the extent of matrix modification
"""
import numpy as np
from scipy.linalg import null_space
from numpy.linalg import inv


def NR_F(mat, m):
    # change into transition matrix
    n = mat.shape[0]
    P1 = mat / mat.sum(1).reshape(-1, 1)

    # diffusion
    P2 = (m - 1) * P1.dot(inv(m * np.eye(n) - P1))

    # compute the stationary distribution and the output matrix
    stationary_d = null_space((P2 - np.eye(n)).T)
    stationary_d = stationary_d / stationary_d.sum()

    if stationary_d[0].shape[0] == 0:
        print('The stationary distribution does not exist！')
    elif stationary_d[0].shape[0] > 1:
        print('The stationary distribution is not unique！')
    else:
        net_new = np.diag(abs(stationary_d.T)[0]).dot(P2)
        net_new = net_new + net_new.T
        output_network = net_new - np.diag(np.diag(net_new))
        output_network = output_network / output_network.max()
        return output_network
