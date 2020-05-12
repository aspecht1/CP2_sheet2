import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# choose esercise

#ex = "1"
ex = "2a"
#ex = "2b"

# EXERCISE 1
# FUNCTIONS
#--------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------


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

    chain = np.zeros(N) # start with x_0 = 0
    chain_removed = np.array([0])
    j = 0
    for i in range(N-1):

        y = (np.random.rand()-0.5)*10
        if next_chain_link(chain[i], y):
            chain[i + 1] = y
        else:
            chain[i + 1] = chain[i]

        if next_chain_link(chain_removed[j], y):
            chain_removed = np.append(chain_removed, y) # append creates new array, does not change array argument
            j += 1

    return chain, chain_removed


# PLOT
#--------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------


if ex == "1":

    N = 100000
    chain, chain_removed = metro_alg(N)

    x_values = np.linspace(-3, 3, N) #x values to plot w(x)
    sns.distplot(chain, label="chain")
    sns.distplot(chain_removed, label="chain removed")
    plt.plot(x_values, w(x_values), label="weight")
    plt.title("N = " + str(N))
    plt.legend()
    plt.show()


#a) little bump at the peak maybe comes from np.random.rand which creates random number between 0 and whithout 1?
#b) chain_removed has  lower peak


# EXERCISE 2
# FUNCTIONS
#--------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------


def update_bounds(lattice):
    """ sets periodic bounds """

    for i in range(N+1): # edges
        lattice[0, i] = lattice[N, i]
        lattice[N+1, i] = lattice[1, i]
        lattice[i, 0] = lattice[i, N]
        lattice[i, N + 1] = lattice[i, 1]

    # corners
    lattice[0, 0] = lattice[N, N]
    lattice[0, N+1] = lattice[N, 1]
    lattice[N+1, 0] = lattice[1, N]
    lattice[N+1, N+1] = lattice[1, 1]

    return lattice



def H(lattice, i, j, h, T):
    """ checks wether spin flip is accepted,
    energy difference only depends on the neighbours of spin flipped site """

    gamma = np.random.rand()
    delta_E = 2*lattice[i, j] * (lattice[i, j - 1] + lattice[i - 1, j] + lattice[i, j + 1] + lattice[i + 1, j]) + 2*h * lattice[i, j]

    return delta_E < 0 or gamma < np.exp(-(delta_E)/(kb * T))



def metro_ising(lattice, L, T, h):
    """ does metropolis algorythm L times and calculates magnetization m """

    m = 0.0

    for i in range(L):
        rand_row = np.random.choice(np.arange(1, N+1))
        rand_col = np.random.choice(np.arange(1, N+1))

        if H(lattice, rand_row, rand_col, h, T):
            lattice[rand_row, rand_col] *= -1

            # check whether bounds have to be updated
            if np.isin(rand_col, [0, N+1]) or np.isin(rand_row, [0, N+1]):
                lattice = update_bounds(lattice)

        m += np.sum(lattice[1:N + 1, 1:N + 1])

    return m/L, lattice



# PLOT a)
#--------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------


chain_lenght = 10000
N = 64
kb = 1 #boltzman constant


# choose a) or b)
#--------------------------------------------------------------------------------------------------------------------------------------


if ex == "2a":

    # a)

    h = 0
    T = [0.1, 0.5, 2, 10]

    plt_num = 1
    for temp in T:
        lattice = update_bounds((np.random.rand(N + 2, N + 2) < 0.5) - 0.5)  # +2 because of periodic bounds
        m, final_lattice = metro_ising(lattice, chain_lenght, temp, h)
        plt.subplot(2, 2, plt_num)
        sns.heatmap(final_lattice[1:N, 1:N], xticklabels=False, yticklabels=False, cbar=False)
        plt.title("T = " + str(temp))
        plt_num += 1

    plt.show()



# PLOT b)
#--------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------


elif ex == "2b":

    #b)

    h_arr = [0, 0.1, 0.5, 1, 5, 10]
    T = np.linspace(0.1, 100, 10)

    i = 0
    for h in h_arr:

        m_val = []
        lattice = update_bounds((np.random.rand(N + 2, N + 2) < 0.5) - 0.5)  # +2 because of periodic bounds
        for temp in T:
            m, _ = metro_ising(lattice, chain_lenght, temp, h_arr[i])
            m_val.append(m)

        plt.plot(T, m_val, label="h = " + str(h_arr[i]))
        i += 1


    plt.ylabel("magnetization m")
    plt.xlabel("Temperature T")
    plt.xlim(T[0], T[-1])
    plt.title("2 b)")
    plt.grid()
    plt.legend()
    plt.show()
