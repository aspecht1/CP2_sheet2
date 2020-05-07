import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def w(x):
    """ weight """

    w = np.exp(-x**2)/np.sqrt(np.pi)
    return w


def next_chain_link(x, y):
    """ checks whether y is accepted as next chain link """

    gamma = np.random.rand()
    alpha = w(y)/w(x)

    return alpha >= gamma


def metro_alg(N):
    """ metropolis algorithm that creates markov chain of lenght N """

    chain = []
    chain_removed = []
    chain.append(0)
    chain_removed.append(0)

    for i in range(N):
        j = 0
        y = (np.random.rand()-0.5)*10
        if next_chain_link(chain[i], y):
            chain.append(y)
        else:
            chain.append(chain[i])

        if next_chain_link(chain_removed[j], y):
            chain_removed.append(y)
            j += 1

    return chain, chain_removed


# N = 100000
# chain, chain_removed = metro_alg(N)
#
# x_values = np.linspace(-3, 3, N) #x values to plot w(x)
# sns.distplot(chain, label="chain")
# sns.distplot(chain_removed, label="chain removed")
# plt.plot(x_values, w(x_values), label="weight")
# plt.legend()
# plt.show()

# a) little bump at the peak probably comes from random.rand which creates random number between 0 and whithout 1?
# b) chain-removed has slightly lower peak but very little


#######################################################################################################################

#2 a)

N = 64
kb = 1 #boltzman constant
index = np.arange(1, N+1) #used to create random indices


# def H(lattice, h):
#     """ calculates the energy H({s_l}) """
#
#     H = 0
#     for i in range(1, N+1):
#         for j in range(1, N+1):
#             H -= lattice[i, j]*(lattice[i, j-1] + lattice[i-1, j]) + h*lattice[i, j]
#             H -= 2*lattice[i, j] * (lattice[i, j - 1] + lattice[i - 1, j] + lattice[i, j + 1] + lattice[i + 1, j]) + 2*h * lattice[i, j]
#
#     return H



# def next_chain_link_ising(x, y, T, h):
#     """ checks whether y is accepted as next chain link """
#
#     gamma = np.random.rand()
#     alpha = np.exp(-(H(y, h) - H(x, h))/(kb * T))
#
#     return alpha >= gamma


def transform_lattice(lattice):
    """ transforms random lattice into lattice of +1/2 and -1/2 and sets periodic bounds """

    for i in range(N+1):
        for j in range(N+1):
            if lattice[i, j] >= 0.5:
                lattice[i, j] = 1/2
            else:
                lattice[i, j] = -1/2

    for i in range(N+1):
        lattice[0, i] = lattice[N, i]
        lattice[N+1, i] = lattice[1, i]
        lattice[i, 0] = lattice[i, N]
        lattice[i, N + 1] = lattice[i, 1]

    lattice[0, 0] = lattice[N, N]
    lattice[0, N+1] = lattice[N, 1]
    lattice[N+1, 0] = lattice[1, N]
    lattice[N+1, N+1] = lattice[1, 1]

    return lattice


def H(lattice, i, j, h, T):
    """ checks wether spin flip is accepted """

    gamma = np.random.rand()
    delta_E = -2*lattice[i, j] * (lattice[i, j - 1] + lattice[i - 1, j] + lattice[i, j + 1] + lattice[i + 1, j]) - 2*h * lattice[i, j]

    return not (delta_E > 0 and np.exp(-(delta_E)/(kb * T)) > gamma)


def metro_ising(L, T, h):
    """ creates markov chain of lenght L and calculates magnetization """

    lattice = transform_lattice(np.random.rand(N + 2, N + 2))  # +2 because of periodic bounds
    ising_chain = [lattice]
    m = 0

    for i in range(L):
        rand_row = np.random.choice(index)
        rand_col = np.random.choice(index)

        if H(ising_chain[i], rand_row, rand_col, h, T):
            new_lattice = ising_chain[i].copy()
            new_lattice[rand_row][rand_col] *= -1
            ising_chain.append(transform_lattice(new_lattice))
        else:
            ising_chain.append(ising_chain[i])

        m += np.sum(ising_chain[i][1:N + 1, 1:N + 1])  # magnetization

    return m



chain_lenght = 100 # 10000 is too big
#T = [0.1, 0.5, 2.0, 10.0]
h_arr = [0.1, 0.5, 1, 5]
T = np.linspace(0.1, 30, 10)

# a)

# chain, _ = metro_ising(chain_lenght, T[0], h[0])
# sns.heatmap(chain[chain_lenght-1][1:N, 1:N], xticklabels=False, yticklabels=False, cbar=False)
# plt.title("T = " + str(T[0]))
# plt.legend()
# plt.show()

# b)
# i = 0
# for h in h_arr:
#     m_val = []
#     for temp in T:
#         m_val.append(metro_ising(chain_lenght, temp, h_arr[i]))
#
#     plt.plot(T, m_val, label="h = " + str(h_arr[i]))
#     plt.ylabel("magnetization m")
#     plt.xlabel("Temperature T")
#     i += 1
m_val = []
for temp in T:
    m = metro_ising(chain_lenght, temp, h_arr[0])
    m_val.append(m/chain_lenght)

plt.plot(T, m_val, label="h = " + str(h_arr[0]))
plt.ylabel("magnetization m")
plt.xlabel("Temperature T")
plt.legend()
plt.show()

# state = 2*np.random.randint(2, size=(10,10))-1
# print(state)

