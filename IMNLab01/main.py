#
# IMN laboratorium nr. 1 Łukasz_Wajda_gr_5 27.10.2021
#

import math
import matplotlib.pyplot as plt

###### Rozwiązanie analityczne #######
def solution(tmax, tmin, step, plambda):
    n = int((tmax - tmin) // step)

    x, y = [[] for _ in range(2)]

    for i in range(n):
        x_n = i * step
        y_n = math.exp(x_n * plambda)
        x.append(x_n)
        y.append(y_n)

    return x, y

############ Rysowanie ###############
plt.rcParams.update({'legend.fontsize': 'large'})
plt.rcParams['axes.titlesize'] = 15
plt.rcParams['axes.labelsize'] = 17


def graph(xl, yl, t, l, n):
    plt.title(t)  
    plt.xlabel(r"t")
    plt.ylabel(r"$\it{u}(t)$")
    fig = plt.gcf()
    fig.set_size_inches(12, 8)
    for i, (x, y) in enumerate(zip(xl, yl)):
      if i<3:
        plt.plot(x, y, marker='o')
      else:
        plt.plot(x, y)
      
    plt.gca().legend(l)
    plt.savefig(n)
    plt.close()


def graph_delta(xl, yl, t, l, n):

    fig, graph1 = plt.subplots()
    fig.set_size_inches(12, 5)
    graph1.set_title(t)
    graph1.set_xlabel(r"t")
    graph1.set_ylabel(r"$\delta u=u_{num}-u_{anal}$")

    for x, y in zip(xl, yl):
        graph1.plot(x, y, marker='o')

    graph1.legend(l)
    plt.savefig(n)
    plt.close()


def graph_delta_mj_euler(xl, yl, t, l, n):
    plt.title(t)
    plt.xlabel(r"t")
    plt.ylabel(r"$\delta u=u_{num}-u_{anal}$")
    fig = plt.gcf()
    fig.set_size_inches(12, 8)

    for x, y in zip(xl, yl):
        plt.plot(x, y, marker='o')

    plt.legend(l)
    plt.savefig(n)
    plt.close()


def graph_RLC(ti, d, title, l, n):
    plt.xlabel(r"t")
    plt.ylabel(r"$\it{Q}$(t)")
    plt.title(title)
    fig = plt.gcf()
    fig.set_size_inches(12, 5)
    for t, data_row in zip(ti, d):
        plt.plot(t, data_row)
    plt.legend(l)
    plt.savefig(n)
    plt.close()

###########  Metoda jawna Eulera  #############

def mj_euler():
    steps = [0.01, 0.1, 1.0]
    tmin = 0
    tmax = 5
    plambda = -1.0
    
    x = [[] for _ in range(4)]
    y = [[] for _ in range(4)] 
    y_delta = [[] for _ in range(3)] 
    for i, step in enumerate(steps):
        n = int((tmax - tmin) // step)
        y_n = 1
        for j in range(n+1):
            y[i].append(y_n)

            x_n = j * step
            x[i].append(x_n)

            y_delta[i].append(y_n - math.exp(x_n * plambda))

            y_n = y_n + step * plambda * y_n

    x[3], y[3] = solution(tmax, tmin, steps[0], plambda)

    graph(x, y, "Task 1 - Euler explicit method - solution", [r"$dt=0.01$ s", r"$dt=0.1$ s", r"$dt=1$ s", "analytical"], 'euler.png' )


    graph_delta_mj_euler(x, y_delta, "Task 1 - Euler explicit method - global error", [r"$dt=0.01$ s", r"$dt=0.1$ s", r"$dt=1$ s"], 'euler_delta.png')




###########  Metoda jawna RK2 (trapezów)  #############

def RK2():
    steps = [0.01, 0.1, 1.0]
    tmin = 0
    tmax = 5
    plambda = -1.0
    x = [[] for _ in range(4)]
    y = [[] for _ in range(4)] 
    y_delta = [[] for _ in range(3)] 
    for i, step in enumerate(steps):

        n = int((tmax - tmin) // step)
        y_n = 1
        for j in range(n+1):
            y[i].append(y_n)
            k1 = plambda * y_n
            k2 = plambda * (y_n + step * k1)

            x_n = j * step
            x[i].append(x_n)

            y_delta[i].append(y_n - math.exp(x_n * plambda))
            y_n = y_n + (step / 2.0) * (k1 + k2)

    x[3], y[3] = solution(tmax, tmin, steps[0], plambda)

    graph(x, y, "Task 2 - RK2 method - solution", [r"$dt=0.01$ s", r"$dt=0.1$ s", r"$dt=1$ s", "analytical"], 'rk2.png')

    graph_delta(x, y_delta, "Task 2 - RK2 method - global error", [r"$dt=0.01$ s", r"$dt=0.1$ s", r"$dt=1$ s"], 'rk2_delta.png',)

###########  Metoda jawna RK4  #############

def RK4():
    steps = [0.01, 0.1, 1.0]
    tmin = 0
    tmax = 5
    plambda = -1.0
    x = [[] for _ in range(4)]
    y = [[] for _ in range(4)] 
    y_delta = [[] for _ in range(3)] 

    for i, step in enumerate(steps):

        n = int((tmax - tmin) // step)
        y_n = 1
        for j in range(n + 1):
            y[i].append(y_n)

            k1 = plambda * y_n
            k2 = plambda * (y_n + (step / 2.0) * k1)
            k3 = plambda * (y_n + (step / 2.0) * k2)
            k4 = plambda * (y_n + step * k1)

            x_n = j * step
            x[i].append(x_n)

            y_delta[i].append(y_n - math.exp(x_n * plambda))
            y_n = y_n + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    x[3], y[3] = solution(tmax, tmin, steps[0], plambda)

    graph(x, y, "Task 3 - RK4 method - solution", [r"$dt=0.01$ s", r"$dt=0.1$ s", r"$dt=1$ s", "analytical"], 'rk4.png',)

    graph_delta(x, y_delta, "Task 3 - RK4 method - global error", [r"$dt=0.01$ s", r"$dt=0.1$ s", r"$dt=1$ s"], 'rk4_delta.png',)


###########   RRZ 2 rzędu na przykładzie obwodu rlc  #############

def RLC():
    R = 100.0
    L = 0.1
    C = 0.001
    step = 1e-4
    w_0 = 1.0 / math.sqrt(L * C)
    T_0 = (2.0 * math.pi) / w_0
    omegas = [0.5 * w_0, 0.8 * w_0, w_0, 1.2 * w_0]
    tmin = 0.0
    tmax = 4.0 * T_0
   
    Q, I, time = [[[] for _ in range(4)] for _ in range(3)]
    n = int((tmax - tmin) // step)

    for i, omega in enumerate(omegas):

        Q_tmp = 0
        I_tmp = 0

        for j in range(n):

            time_tmp = j * step
            time_half = (j + 0.5) * step
            time_next = (j + 1.0) * step

            v = 10.0 * math.sin(time_tmp * omega) 
            v_half = 10.0 * math.sin(time_half * omega)
            v_next = 10.0 * math.sin(time_next * omega)

            Q[i].append(Q_tmp)
            I[i].append(I_tmp)
            time[i].append(time_tmp)

            k1Q = I_tmp
            k1I = (v / L) - ((1.0 / (L * C)) * Q_tmp) - ((R / L) * I_tmp)

            k2Q = I_tmp + ((step / 2.0) * k1I)
            k2I = (v_half / L) - ((1.0 / (L * C)) * (Q_tmp + (step / 2.0) * k1Q)) - \
                  ((R / L) * (I_tmp + (step / 2.0) * k1I))

            k3Q = I_tmp + (step / 2.0) * k2I
            k3I = (v_half / L) - ((1.0 / (L * C)) * (Q_tmp + (step / 2.0) * k2Q)) - \
                  ((R / L) * (I_tmp + (step / 2.0) * k2I))

            k4Q = I_tmp + (step * k3I)
            k4I = (v_next / L) - (1 / (L * C)) * (Q_tmp + step * k3Q) - \
                  ((R / L) * (I_tmp + (step * k3I)))

            Q_tmp = Q_tmp + (step / 6.0) * (k1Q + 2 * k2Q + 2 * k3Q + k4Q)
            I_tmp = I_tmp + (step / 6.0) * (k1I + 2 * k2I + 2 * k3I + k4I)

    graph_RLC(time, Q, "Task 4 - RK4 method - Q(t) diagram", [r"$\omega_v$=0.5$\omega_0$", r"$\omega_v$=0.8$\omega_0$", r"$\omega_v$=1.0$\omega_0$", r"$\omega_v$=1.2$\omega_0$"], 'Q_t.png')


if __name__ == "__main__":
    mj_euler()
    RK2()
    RK4()
    RLC()
