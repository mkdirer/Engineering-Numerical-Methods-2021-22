import matplotlib.pyplot as plt

############### rysowanie ####################

plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams.update({'legend.fontsize': 'large'})


def graph(x_s, y_s, title, legend, x_lab, y_lab, name):
    plt.title(title)
    plt.ylabel(y_lab)
    plt.xlabel(x_lab)
    figure = plt.gcf()
    figure.set_size_inches(12, 8)
    for x, y in zip(x_s, y_s):
        plt.plot(x, y)
    plt.gca().legend(legend)
    plt.savefig(name)
    plt.close()

################# obliczenia ###################

def loop_time(metoda):
    TOL = [1e-2, 1e-5]
    T_M = 40.0
    Alpha = 5.0
    S = 0.75
    p = 2.0
    delta_t_arr = [[] for _ in range(2)]
    x_arr = [[] for _ in range(2)]
    v_arr = [[] for _ in range(2)]
    t_arr = [[] for _ in range(2)]
    for i, tol in enumerate(TOL):
        x = 0.01
        v = 0.0
        t = 0.0
        dt = 1.0
        while t < T_M:

            st_1 = metoda(x, v, dt, Alpha)
            st_1 = metoda(st_1[0], st_1[1], dt, Alpha)

            st_2 = metoda(x, v, 2.0 * dt, Alpha)

            err_x = (st_2[0] - st_1[0]) / (pow(2, p) - 1.0)
            err_v = (st_2[1] - st_1[1]) / (pow(2, p) - 1.0)

            max_err = max(abs(err_x), abs(err_v))

            if max_err < tol:
                t = t + 2 * dt
                x = st_1[0]
                v = st_1[1]
                t_arr[i].append(t)
                x_arr[i].append(x)
                v_arr[i].append(v)
                delta_t_arr[i].append(dt)

            dt *= pow((S * tol) / max_err, (1.0 / (p + 1.0)))

    return x_arr, v_arr, t_arr, delta_t_arr

def g(al, x, v):
    return al * (1.0 - pow(x, 2)) * v - x

######################### metoda trapezów ######################


def count_trapez(x, v, delta_t, alp):
    delta = 1e-10

    x_next = x
    v_next = v

    while True:
        F = x_next - x - (delta_t / 2.0) * (v + v_next)
        G = v_next - v - (delta_t / 2.0) * (g(alp, x, v) + g(alp, x_next, v_next))
        a = [[1.0, -delta_t / 2.0], [(-delta_t / 2.0) * (-2.0 * alp * x_next * v_next - 1.0), 1.0 - (delta_t / 2.0) *
                                     alp * (1.0 - pow(x_next, 2.0))]]

        denominator = a[0][0] * a[1][1] - a[0][1] * a[1][0]

        delta_x = (-F * a[1][1] - (-G) * a[0][1]) / denominator
        delta_v = (a[0][0] * (-G) - a[1][0] * (-F)) / denominator

        x_next += delta_x
        v_next += delta_v

        if delta_x < delta and delta_v < delta:
            break

    return x_next, v_next


def trapez():
    x_arr, v_arr, t_arr, delta_t_arr = loop_time(count_trapez)
    legend = [r"TOL = $10^{-2}$", r"TOL = $10^{-5}$"]
    title = "Niejawana Metoda Trapezow (OscylatorVanderPola a=5)"

    graph(t_arr, x_arr, title, legend, r"$\it{t}(s)$", r"$\it{x}(t)$", 'metoda_trapezow_x_t.png')
    graph(t_arr, v_arr, title, legend, r"$\it{t}(s)$", r"$\it{v}(t)$", 'metoda_trapezow_v_t.png')
    graph(t_arr, delta_t_arr, title, legend, r"$\it{t}(s)$", r"$\it{dt}(t)$", 'metoda_trapezow_dt_t.png')
    graph(x_arr, v_arr, title, legend, r"$\it{x}$", r"$\it{v}$",'metoda_trapezow_v_x.png')


####################### metoda rk2 #####################

def count_rk2(x, v, delta_t, alp):
    k1_x = v
    k1_v = g(alp, x, v)

    k2_x = v + delta_t * k1_v
    k2_v = g(alp, x + delta_t * k1_x, v + delta_t * k1_v)

    return x + (delta_t / 2.0) * (k1_x + k2_x), v + (delta_t / 2.0) * (k1_v + k2_v)


def rk2():
    x_arr, v_arr, t_arr, delta_t_arr = loop_time(count_rk2)
    legend = [r"TOL = $10^{-2}$", r"TOL = $10^{-5}$"]
    title = "RK2 (OscylatorVanderPola a=5)"

    graph(t_arr, x_arr, title, legend, r"$\it{t}(s)$", r"$\it{x}(t)$", 'rk2_x_t.png')
    graph(t_arr, v_arr, title, legend, r"$\it{t}(s)$", r"$\it{v}(t)$", 'rk2_v_t.png')
    graph(t_arr, delta_t_arr, title, legend, r"$\it{t}(s)$", r"$\it{dt}(t)$", 'rk2_dt_t.png')
    graph(x_arr, v_arr, title, legend, r"$\it{x}$", r"$\it{v}$",'rk2_v_x.png')



if __name__ == "__main__":
    rk2()
    trapez()