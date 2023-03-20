import numpy as np

# %%
def qubo2ising(Q):
    """
    Convert a qubo problem matrix into an Ising problem
    """
    S = (Q + Q.T) / 2 # making sure Q is symmetric, not necessary in this file

    J = (S - np.diag(np.diag(S))) / 4 # coupling
    h = np.sum(S, axis=0) / 2 # local field
    offset = (np.sum(S) + np.sum(np.diag(S))) / 4 # offset energy, add this back to get the original QUBO energy

    return J, h, offset


# %% [markdown]
# ## Convert to QUBO
# QUBO item variables: {1, 0}
# Objective: Minimize Q.dot(state).dot(state)

# %%
def ising2qubo(J, h):
    """
    Convert a Ising problem into a QUBO problem matrix
    """
    S = (J + J.T) / 2 # making sure Q is symmetric, not necessary in this file

    Q = 4*J + np.diag(2*h - 4*np.sum(J, axis=0))
    offset = np.sum(J) - np.sum(h) # offset energy, add this back to get the original Ising energy

    return Q, offset