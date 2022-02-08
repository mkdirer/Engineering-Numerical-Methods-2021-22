import matplotlib.pyplot as plt
import numpy as np
import math


########### rysowanie ###############
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams.update({'legend.fontsize': 'large'})


def chart_generator(x, y, z, title, legend, filename):
    plt.ylabel(r"$\it{u}(t)$, $\it{z}(t)$")
    plt.xlabel(r"$\it{t}(s)$")
    plt.title(title)
    figure = plt.gcf()
    figure.set_size_inches(12, 8)
    plt.plot(x, y)
    plt.plot(x, z)
    plt.xticks(np.arange(0, 101, 10))
    plt.yticks(np.arange(0, 501, 50))
    plt.grid(True)
    plt.gca().legend(legend)
    plt.savefig(filename)
    plt.close()


############# funkcja #############
def f(alpha, beta, x):
    return alpha * x - beta * math.pow(x, 2)


############## metoda trapezów ##############


def metoda_trapezow():
    beta = 0.001
    N = 500
    gamma = 0.1
    t_max = 100
    dt = 0.1
    TOL = 1e-6
    MI = 20
    n = int(t_max // dt) + 1
    alpha = (beta * N) - gamma

    time = np.linspace(0, n * dt, num=n, endpoint=True)

    picard = np.ones(time.shape)
    newton = np.ones(time.shape)
    z = np.zeros(time.shape)
    z[0] = N - 1.0

    # Picard
    for i in range(1, n):
        u_prev = picard[i - 1]
        u = 0
        for mi in range(MI+1):
            if abs(u - u_prev) < TOL:
                break

            u = picard[i - 1] + (dt / 2.0) * (f(alpha, beta, picard[i - 1]) +
                                              f(alpha, beta, u_prev))

        picard[i] = u
        z[i] = N - picard[i]

    chart_generator(time, picard, z, "Niejawna Metoda Trapezow z iteraja Picarda", [r"$\it{u}(t)$", r"$\it{z}(t)=N-\it{u}(t)$"],'picard.png')

    # Newton
    for i in range(1, n):
        u_prev = newton[i - 1]
        u = 0

        for mi in range(MI+1):
            if abs(u - u_prev) < TOL:
                break

            u = u_prev - (u_prev - newton[i - 1] - (dt / 2.0) * (f(alpha, beta, newton[i - 1]) +
                                                                 f(alpha, beta, u_prev))) / \
                (1 - (dt / 2.0) * (alpha - 2 * beta * u_prev))

        newton[i] = u
        z[i] = N - newton[i]

    chart_generator(time, newton, z, "Niejawna Metoda Trapezow z iteracja Newtona", [r"$\it{u}(t)$", r"$\it{z}(t)=N-\it{u}(t)$"], 'newton.png',)

############### metoda RK2 ###############


def rk2():
    beta = 0.001
    N = 500
    gamma = 0.1
    t_max = 100
    delta_t = 0.1
    TOL = 1e-6
    MI = 20
    ITERATIONS = int(t_max / delta_t)
    alpha = (beta * N) - gamma

    a = np.array([[0.25, 0.25 - (math.sqrt(3) / 6.0)], [0.25 + (math.sqrt(3) / 6.0), 0.25]])
    b = np.array([0.5, 0.5])
    c = np.array([0.5 - (math.sqrt(3) / 6.0), 0.5 + (math.sqrt(3) / 6.0)])

    def f1(U1, U2, u):
        return U1 - u - (delta_t * (a[0][0] * f(alpha, beta, U1) + a[0][1] * f(alpha, beta, U2)))

    def f2(U1, U2, u):
        return U2 - u - (delta_t * (a[1][0] * f(alpha, beta, U1) + a[1][1] * f(alpha, beta, U2)))

    def m_value(i, j, u):
        res = (-delta_t) * a[i][j] * (alpha - 2.0 * beta * u)
        if i == j:
            res += 1.0

        return res

    time = np.linspace(0, ITERATIONS * delta_t, num=ITERATIONS + 1, endpoint=True)

    U = np.ones(time.shape)
    z = np.zeros(time.shape)
    z[0] = N - 1.0

    for i in range(1, ITERATIONS + 1):
        U_pre = U1 = U2 = U[i - 1]

        for mi in range(MI+1):

            m = np.array([[m_value(0, 0, U1), m_value(0, 1, U2)],
                          [m_value(1, 0, U1), m_value(1, 1, U2)]])

            F1 = f1(U1, U2, U[i - 1])
            F2 = f2(U1, U2, U[i - 1])

            denominator = (m[0][0] * m[1][1]) - (m[0][1] * m[1][0])
            d_U1 = ((F2 * m[0][1]) - (F1 * m[1][1])) / denominator
            d_U2 = ((F1 * m[1][0]) - (F2 * m[0][0])) / denominator

            U1 += d_U1
            U2 += d_U2

            U_akt = U[i-1] + delta_t * (b[0] * f(alpha, beta, U1) + b[1] * f(alpha, beta, U2))

            if abs(U_akt-U_pre) < TOL:
                break
            U_pre = U_akt
        U[i] = U_akt
        z[i] = N - U[i]

    chart_generator(time, U, z, "Niejawna metoda RK2", [r"$\it{u}(t)$", r"$\it{z}(t)=N-\it{u}(t)$"], 'rk2.png',)


if __name__ == "__main__":
    metoda_trapezow()
    rk2()