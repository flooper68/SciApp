import json

import matplotlib.pyplot as plt
import mpld3

from funkce import *

n_terms = 5
n_select = 1
show_figures = False


def get_Q_matrix(perm1, perm2, period, n_terms, fill_factor):
    def q_cap_ref1_integral(x):
        def integrand_function(y):
            if (y < (1 - fill_factor) * period):
                return perm1[0, 0] * np.exp(-1j * x * 2 * np.pi / period * y)
            elif (y < period):
                return perm2[0, 0] * np.exp(-1j * x * 2 * np.pi / period * y)

        return 1 / period * (complex_quadrature(integrand_function, 0, period, 100))

    q_cap_ref1 = get_C_matrix(period, q_cap_ref1_integral, n_terms)

    def q_ref_grating3_integral(x):
        def integrand_function(y):
            if (y < (1 - fill_factor) * period):
                return perm1[2, 2] * np.exp(-1j * x * 2 * np.pi / period * y)
            elif (y < period):
                return perm2[2, 2] * np.exp(-1j * x * 2 * np.pi / period * y)

        return 1 / period * (complex_quadrature(integrand_function, 0, period, 100))

    q_ref_grating3 = get_C_matrix(period, q_ref_grating3_integral, n_terms)

    def q_cap_grating2_integral(x):
        def integrand_function(y):
            if (y < (1 - fill_factor) * period):
                return 1 / perm1[1, 1] * np.exp(-1j * x * 2 * np.pi / period * y)
            elif (y < period):
                return 1 / perm2[1, 1] * np.exp(-1j * x * 2 * np.pi / period * y)

        return 1 / period * (complex_quadrature(integrand_function, 0, period, 100))

    q_cap_grating2 = get_C_matrix(period, q_cap_grating2_integral, n_terms) ** -1

    def q_0_integral(x):
        def integrand_function(y):
            if (y < (1 - fill_factor) * period):
                return 0 * np.exp(-1j * x * 2 * np.pi / period * y)
            elif (y < period):
                return 0 * np.exp(-1j * x * 2 * np.pi / period * y)

        return 1 / period * (complex_quadrature(integrand_function, 0, period, 100))

    q_0 = get_C_matrix(period, q_0_integral, n_terms)

    return q_cap_ref1, q_0, q_0, q_0, q_cap_grating2, q_0, q_0, q_0, q_ref_grating3


def get_incoh_wv(options):
    print('tick_flooper', flush=True)
    n_air = 1

    n_capping = complex(options['n_capping'])
    n_grating = complex(options['n_grating'])
    n_reflective = complex(options['n_reflective'])

    width_capping = float(options['width_capping'])
    width_1 = float(options['width_reflective'])
    width_2 = float(options['width_grating']) - float(options['width_reflective'])
    width_3 = float(options['width_reflective'])

    angle = float(options['angle'])
    period = float(options['grating_period'])
    fill_factor = float(options['fill_factor'])

    divisions = int(options['divisions'])
    start = float(options['divisions_start'])
    end = float(options['divisions_end'])

    perm_grating = np.array([[n_grating ** 2, 0, 0], [0, n_grating ** 2, 0], [0, 0, n_grating ** 2]])
    perm_reflective = np.array([[n_reflective ** 2, 0, 0], [0, n_reflective ** 2, 0], [0, 0, n_reflective ** 2]])
    perm_capping = np.array([[n_capping ** 2, 0, 0], [0, n_capping ** 2, 0], [0, 0, n_capping ** 2]])

    # calculate localy just once
    print('tick_flooper', flush=True)
    q_capping = get_Q_matrix(perm_capping, perm_capping, period, n_terms, 1)
    print('tick_flooper', flush=True)
    q_cap_ref = get_Q_matrix(perm_capping, perm_reflective, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_cap_grating = get_Q_matrix(perm_capping, perm_grating, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_ref_grating = get_Q_matrix(perm_reflective, perm_grating, period, n_terms, fill_factor)

    points = np.linspace(start, end, divisions)

    # Define the function
    def wvDep_incoh(wv):
        print('tick_flooper', flush=True)
        t0 = get_outer_T(n_terms, n_air, angle, wv, period, n_air)
        t_incoherent = get_outer_T(n_terms, n_capping, angle, wv, period, n_air)

        t1, p1 = get_periodic_layer_T_P_aniso(period, q_cap_ref[0], q_cap_ref[1], q_cap_ref[2], q_cap_ref[3],
                                              q_cap_ref[4], q_cap_ref[5], q_cap_ref[6], q_cap_ref[7], q_cap_ref[8],
                                              n_terms, angle, 0, wv, n_air, width_1)
        t2, p2 = get_periodic_layer_T_P_aniso(period, q_cap_grating[0], q_cap_grating[1], q_cap_grating[2],
                                              q_cap_grating[3], q_cap_grating[4], q_cap_grating[5], q_cap_grating[6],
                                              q_cap_grating[7], q_cap_grating[8], n_terms, angle, 0, wv, n_air, width_2)
        t3, p3 = get_periodic_layer_T_P_aniso(period, q_ref_grating[0], q_ref_grating[1], q_ref_grating[2],
                                              q_ref_grating[3], q_ref_grating[4], q_ref_grating[5], q_ref_grating[6],
                                              q_ref_grating[7], q_ref_grating[8], n_terms, angle, 0, wv, n_air, width_3)
        t4 = get_outer_T(n_terms, n_grating, angle, wv, period, n_air)

        s_0 = get_s_boundary(t0, t_incoherent)
        T_uu, R_ud, R_du, T_dd = split_matrix(s_0)
        r_01, t_01, r_10, t_10 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select,
                                                                                     T_dd), select_modes(n_terms,
                                                                                                         n_select,
                                                                                                         R_du), select_modes(
            n_terms, n_select, T_uu)

        s_1 = get_s(t_incoherent, t1, p1)
        s_2 = get_s(t1, t2, p2)
        s_3 = get_s(t2, t3, p3)
        s_boundary = get_s_boundary(t3, t4)

        s = join_s_matrices_recur(s_1, s_2)
        s = join_s_matrices_recur(s, s_3)
        T_uu, R_ud, R_du, T_dd = join_s_matrices(s, s_boundary)
        r_12, t_12 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select, T_dd)

        a = prop_trans(n_select)

        # Propagační matice pro každý mód
        p_13 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)
        p_24 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)

        # Tvorba koherenčních matic
        cr_01 = np.kron(r_01, np.conj(r_01))
        ct_01 = np.kron(t_01, np.conj(t_01))
        cr_10 = np.kron(r_10, np.conj(r_10))
        ct_10 = np.kron(t_10, np.conj(t_10))
        cr_12 = np.kron(r_12, np.conj(r_12))
        ct_12 = np.kron(t_12, np.conj(t_12))
        cp_13 = np.kron(p_13, np.conj(p_13))
        cp_24 = np.kron(p_24, np.conj(p_24))

        coh_reflektance = cr_01 + ct_10 * cp_24 * cr_12 * cp_13 * (
                np.matrix(np.identity((2 * (2 * n_select + 1)) ** 2)) - cr_10 * cp_24 * cr_12 * cp_13) ** -1 * ct_01

        order = 0
        coh_matrix_0 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air, wv,
                                                                                              period)
        order = 1
        if (n_select > 0):
            coh_matrix_1 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_1 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
        order = 2
        if (n_select > 1):
            coh_matrix_2 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_2 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

        return float(np.real(coh_matrix_0[0, 0])), float(np.real(coh_matrix_0[3, 3])), float(
            np.real(coh_matrix_1[0, 0])), float(np.real(coh_matrix_1[3, 3]))

    values = np.array(list(map(wvDep_incoh, points)))

    toSend = {'data': []}

    for x in range(len(points)):
        tempDict = {'point': points[x], 'Rs0': values[x, 0], 'Rp0': values[x, 1]}
        toSend['data'].append(tempDict)

    fig1 = plt.figure(1, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 0], 'r-', label='$R_s$ 0 order')
    plt.plot(points, values[:, 1], 'b-', label='$R_p$ 0 order')
    plt.xlabel('Wavelength [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot1'] = mpld3.fig_to_dict(fig1)

    fig2 = plt.figure(2, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 2], 'r-', label='$R_s$ -1 order')
    plt.plot(points, values[:, 3], 'b-', label='$R_p$ -1 order')
    plt.xlabel('Wavelength [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot2'] = mpld3.fig_to_dict(fig2)

    print('start sending', flush=True)
    print(json.dumps(toSend), flush=True)
    print('end sending', flush=True)


def get_incoh_angle(options):
    print('tick_flooper', flush=True)
    n_air = 1

    n_capping = complex(options['n_capping'])
    n_grating = complex(options['n_grating'])
    n_reflective = complex(options['n_reflective'])

    width_capping = float(options['width_capping'])
    width_1 = float(options['width_reflective'])
    width_2 = float(options['width_grating']) - float(options['width_reflective'])
    width_3 = float(options['width_reflective'])

    wv = float(options['wavelength'])
    period = float(options['grating_period'])
    fill_factor = float(options['fill_factor'])

    divisions = int(options['divisions'])
    start = float(options['divisions_start'])
    end = float(options['divisions_end'])

    perm_grating = np.array([[n_grating ** 2, 0, 0], [0, n_grating ** 2, 0], [0, 0, n_grating ** 2]])
    perm_reflective = np.array([[n_reflective ** 2, 0, 0], [0, n_reflective ** 2, 0], [0, 0, n_reflective ** 2]])
    perm_capping = np.array([[n_capping ** 2, 0, 0], [0, n_capping ** 2, 0], [0, 0, n_capping ** 2]])

    # calculate localy just once
    print('tick_flooper', flush=True)
    q_capping = get_Q_matrix(perm_capping, perm_capping, period, n_terms, 1)
    print('tick_flooper', flush=True)
    q_cap_ref = get_Q_matrix(perm_capping, perm_reflective, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_cap_grating = get_Q_matrix(perm_capping, perm_grating, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_ref_grating = get_Q_matrix(perm_reflective, perm_grating, period, n_terms, fill_factor)

    points = np.linspace(start, end, divisions)

    # Define the function
    def angleDep_incoh(angle):
        angle -= 0.0001
        print('tick_flooper', flush=True)
        t0 = get_outer_T(n_terms, n_air, angle, wv, period, n_air)
        t_incoherent = get_outer_T(n_terms, n_capping, angle, wv, period, n_air)

        t1, p1 = get_periodic_layer_T_P_aniso(period, q_cap_ref[0], q_cap_ref[1], q_cap_ref[2], q_cap_ref[3],
                                              q_cap_ref[4], q_cap_ref[5], q_cap_ref[6], q_cap_ref[7], q_cap_ref[8],
                                              n_terms, angle, 0, wv, n_air, width_1)
        t2, p2 = get_periodic_layer_T_P_aniso(period, q_cap_grating[0], q_cap_grating[1], q_cap_grating[2],
                                              q_cap_grating[3], q_cap_grating[4], q_cap_grating[5], q_cap_grating[6],
                                              q_cap_grating[7], q_cap_grating[8], n_terms, angle, 0, wv, n_air, width_2)
        t3, p3 = get_periodic_layer_T_P_aniso(period, q_ref_grating[0], q_ref_grating[1], q_ref_grating[2],
                                              q_ref_grating[3], q_ref_grating[4], q_ref_grating[5], q_ref_grating[6],
                                              q_ref_grating[7], q_ref_grating[8], n_terms, angle, 0, wv, n_air, width_3)
        t4 = get_outer_T(n_terms, n_grating, angle, wv, period, n_air)

        s_0 = get_s_boundary(t0, t_incoherent)
        T_uu, R_ud, R_du, T_dd = split_matrix(s_0)
        r_01, t_01, r_10, t_10 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select,
                                                                                     T_dd), select_modes(n_terms,
                                                                                                         n_select,
                                                                                                         R_du), select_modes(
            n_terms, n_select, T_uu)

        s_1 = get_s(t_incoherent, t1, p1)
        s_2 = get_s(t1, t2, p2)
        s_3 = get_s(t2, t3, p3)
        s_boundary = get_s_boundary(t3, t4)

        s = join_s_matrices_recur(s_1, s_2)
        s = join_s_matrices_recur(s, s_3)
        T_uu, R_ud, R_du, T_dd = join_s_matrices(s, s_boundary)
        r_12, t_12 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select, T_dd)

        a = prop_trans(n_select)

        # Propagační matice pro každý mód
        p_13 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)
        p_24 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)

        # Tvorba koherenčních matic
        cr_01 = np.kron(r_01, np.conj(r_01))
        ct_01 = np.kron(t_01, np.conj(t_01))
        cr_10 = np.kron(r_10, np.conj(r_10))
        ct_10 = np.kron(t_10, np.conj(t_10))
        cr_12 = np.kron(r_12, np.conj(r_12))
        ct_12 = np.kron(t_12, np.conj(t_12))
        cp_13 = np.kron(p_13, np.conj(p_13))
        cp_24 = np.kron(p_24, np.conj(p_24))

        coh_reflektance = cr_01 + ct_10 * cp_24 * cr_12 * cp_13 * (
                np.matrix(np.identity((2 * (2 * n_select + 1)) ** 2)) - cr_10 * cp_24 * cr_12 * cp_13) ** -1 * ct_01

        order = 0
        coh_matrix_0 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air, wv,
                                                                                              period)
        order = 1
        if (n_select > 0):
            coh_matrix_1 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_1 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
        order = 2
        if (n_select > 1):
            coh_matrix_2 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_2 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

        return float(np.real(coh_matrix_0[0, 0])), float(np.real(coh_matrix_0[3, 3])), float(
            np.real(coh_matrix_1[0, 0])), float(np.real(coh_matrix_1[3, 3]))

    values = np.array(list(map(angleDep_incoh, points)))

    toSend = {'data': []}

    for x in range(len(points)):
        tempDict = {'point': points[x], 'Rs0': values[x, 0], 'Rp0': values[x, 1]}
        toSend['data'].append(tempDict)

    fig1 = plt.figure(1, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 0], 'r-', label='$R_s$ 0 order')
    plt.plot(points, values[:, 1], 'b-', label='$R_p$ 0 order')
    plt.xlabel('Angle [deg]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot1'] = mpld3.fig_to_dict(fig1)

    fig2 = plt.figure(2, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 2], 'r-', label='$R_s$ -1 order')
    plt.plot(points, values[:, 3], 'b-', label='$R_p$ -1 order')
    plt.xlabel('Angle [deg]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot2'] = mpld3.fig_to_dict(fig2)

    print('start sending', flush=True)
    print(json.dumps(toSend), flush=True)
    print('end sending', flush=True)


def get_incoh_fill(options):
    print('tick_flooper', flush=True)
    n_air = 1

    n_capping = complex(options['n_capping'])
    n_grating = complex(options['n_grating'])
    n_reflective = complex(options['n_reflective'])

    width_capping = float(options['width_capping'])
    width_1 = float(options['width_reflective'])
    width_2 = float(options['width_grating']) - float(options['width_reflective'])
    width_3 = float(options['width_reflective'])

    wv = float(options['wavelength'])
    period = float(options['grating_period'])
    angle = float(options['angle'])

    divisions = int(options['divisions'])
    start = float(options['divisions_start'])
    end = float(options['divisions_end'])

    perm_grating = np.array([[n_grating ** 2, 0, 0], [0, n_grating ** 2, 0], [0, 0, n_grating ** 2]])
    perm_reflective = np.array([[n_reflective ** 2, 0, 0], [0, n_reflective ** 2, 0], [0, 0, n_reflective ** 2]])
    perm_capping = np.array([[n_capping ** 2, 0, 0], [0, n_capping ** 2, 0], [0, 0, n_capping ** 2]])

    points = np.linspace(start, end, divisions)

    # Define the function
    def fillDep_incoh(fill_factor):
        print('tick_flooper', flush=True)
        q_capping = get_Q_matrix(perm_capping, perm_capping, period, n_terms, 1)
        q_cap_ref = get_Q_matrix(perm_capping, perm_reflective, period, n_terms, fill_factor)
        q_cap_grating = get_Q_matrix(perm_capping, perm_grating, period, n_terms, fill_factor)
        q_ref_grating = get_Q_matrix(perm_reflective, perm_grating, period, n_terms, fill_factor)

        t0 = get_outer_T(n_terms, n_air, angle, wv, period, n_air)
        t_incoherent = get_outer_T(n_terms, n_capping, angle, wv, period, n_air)

        t1, p1 = get_periodic_layer_T_P_aniso(period, q_cap_ref[0], q_cap_ref[1], q_cap_ref[2], q_cap_ref[3],
                                              q_cap_ref[4], q_cap_ref[5], q_cap_ref[6], q_cap_ref[7], q_cap_ref[8],
                                              n_terms, angle, 0, wv, n_air, width_1)
        t2, p2 = get_periodic_layer_T_P_aniso(period, q_cap_grating[0], q_cap_grating[1], q_cap_grating[2],
                                              q_cap_grating[3], q_cap_grating[4], q_cap_grating[5], q_cap_grating[6],
                                              q_cap_grating[7], q_cap_grating[8], n_terms, angle, 0, wv, n_air, width_2)
        t3, p3 = get_periodic_layer_T_P_aniso(period, q_ref_grating[0], q_ref_grating[1], q_ref_grating[2],
                                              q_ref_grating[3], q_ref_grating[4], q_ref_grating[5], q_ref_grating[6],
                                              q_ref_grating[7], q_ref_grating[8], n_terms, angle, 0, wv, n_air, width_3)
        t4 = get_outer_T(n_terms, n_grating, angle, wv, period, n_air)

        s_0 = get_s_boundary(t0, t_incoherent)
        T_uu, R_ud, R_du, T_dd = split_matrix(s_0)
        r_01, t_01, r_10, t_10 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select,
                                                                                     T_dd), select_modes(n_terms,
                                                                                                         n_select,
                                                                                                         R_du), select_modes(
            n_terms, n_select, T_uu)

        s_1 = get_s(t_incoherent, t1, p1)
        s_2 = get_s(t1, t2, p2)
        s_3 = get_s(t2, t3, p3)
        s_boundary = get_s_boundary(t3, t4)

        s = join_s_matrices_recur(s_1, s_2)
        s = join_s_matrices_recur(s, s_3)
        T_uu, R_ud, R_du, T_dd = join_s_matrices(s, s_boundary)
        r_12, t_12 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select, T_dd)

        a = prop_trans(n_select)

        # Propagační matice pro každý mód
        p_13 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)
        p_24 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)

        # Tvorba koherenčních matic
        cr_01 = np.kron(r_01, np.conj(r_01))
        ct_01 = np.kron(t_01, np.conj(t_01))
        cr_10 = np.kron(r_10, np.conj(r_10))
        ct_10 = np.kron(t_10, np.conj(t_10))
        cr_12 = np.kron(r_12, np.conj(r_12))
        ct_12 = np.kron(t_12, np.conj(t_12))
        cp_13 = np.kron(p_13, np.conj(p_13))
        cp_24 = np.kron(p_24, np.conj(p_24))

        coh_reflektance = cr_01 + ct_10 * cp_24 * cr_12 * cp_13 * (
                np.matrix(np.identity((2 * (2 * n_select + 1)) ** 2)) - cr_10 * cp_24 * cr_12 * cp_13) ** -1 * ct_01

        order = 0
        coh_matrix_0 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air, wv,
                                                                                              period)
        order = 1
        if (n_select > 0):
            coh_matrix_1 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_1 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
        order = 2
        if (n_select > 1):
            coh_matrix_2 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_2 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

        return float(np.real(coh_matrix_0[0, 0])), float(np.real(coh_matrix_0[3, 3])), float(
            np.real(coh_matrix_1[0, 0])), float(np.real(coh_matrix_1[3, 3]))

    values = np.array(list(map(fillDep_incoh, points)))

    toSend = {'data': []}

    for x in range(len(points)):
        tempDict = {'point': points[x], 'Rs0': values[x, 0], 'Rp0': values[x, 1]}
        toSend['data'].append(tempDict)

    fig1 = plt.figure(1, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 0], 'r-', label='$R_s$ 0 order')
    plt.plot(points, values[:, 1], 'b-', label='$R_p$ 0 order')
    plt.xlabel('Fill factor')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot1'] = mpld3.fig_to_dict(fig1)

    fig2 = plt.figure(2, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 2], 'r-', label='$R_s$ -1 order')
    plt.plot(points, values[:, 3], 'b-', label='$R_p$ -1 order')
    plt.xlabel('Fill factor')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot2'] = mpld3.fig_to_dict(fig2)

    print('start sending', flush=True)
    print(json.dumps(toSend), flush=True)
    print('end sending', flush=True)


def get_incoh_period(options):
    print('tick_flooper', flush=True)
    n_air = 1

    n_capping = complex(options['n_capping'])
    n_grating = complex(options['n_grating'])
    n_reflective = complex(options['n_reflective'])

    width_capping = float(options['width_capping'])
    width_1 = float(options['width_reflective'])
    width_2 = float(options['width_grating']) - float(options['width_reflective'])
    width_3 = float(options['width_reflective'])

    wv = float(options['wavelength'])
    fill_factor = float(options['fill_factor'])
    angle = float(options['angle'])

    divisions = int(options['divisions'])
    start = float(options['divisions_start'])
    end = float(options['divisions_end'])

    perm_grating = np.array([[n_grating ** 2, 0, 0], [0, n_grating ** 2, 0], [0, 0, n_grating ** 2]])
    perm_reflective = np.array([[n_reflective ** 2, 0, 0], [0, n_reflective ** 2, 0], [0, 0, n_reflective ** 2]])
    perm_capping = np.array([[n_capping ** 2, 0, 0], [0, n_capping ** 2, 0], [0, 0, n_capping ** 2]])

    points = np.linspace(start, end, divisions)

    # Define the function
    def periodDep_incoh(period):
        print('tick_flooper', flush=True)
        q_capping = get_Q_matrix(perm_capping, perm_capping, period, n_terms, 1)
        q_cap_ref = get_Q_matrix(perm_capping, perm_reflective, period, n_terms, fill_factor)
        q_cap_grating = get_Q_matrix(perm_capping, perm_grating, period, n_terms, fill_factor)
        q_ref_grating = get_Q_matrix(perm_reflective, perm_grating, period, n_terms, fill_factor)

        t0 = get_outer_T(n_terms, n_air, angle, wv, period, n_air)
        t_incoherent = get_outer_T(n_terms, n_capping, angle, wv, period, n_air)

        t1, p1 = get_periodic_layer_T_P_aniso(period, q_cap_ref[0], q_cap_ref[1], q_cap_ref[2], q_cap_ref[3],
                                              q_cap_ref[4], q_cap_ref[5], q_cap_ref[6], q_cap_ref[7], q_cap_ref[8],
                                              n_terms, angle, 0, wv, n_air, width_1)
        t2, p2 = get_periodic_layer_T_P_aniso(period, q_cap_grating[0], q_cap_grating[1], q_cap_grating[2],
                                              q_cap_grating[3], q_cap_grating[4], q_cap_grating[5], q_cap_grating[6],
                                              q_cap_grating[7], q_cap_grating[8], n_terms, angle, 0, wv, n_air, width_2)
        t3, p3 = get_periodic_layer_T_P_aniso(period, q_ref_grating[0], q_ref_grating[1], q_ref_grating[2],
                                              q_ref_grating[3], q_ref_grating[4], q_ref_grating[5], q_ref_grating[6],
                                              q_ref_grating[7], q_ref_grating[8], n_terms, angle, 0, wv, n_air, width_3)
        t4 = get_outer_T(n_terms, n_grating, angle, wv, period, n_air)

        s_0 = get_s_boundary(t0, t_incoherent)
        T_uu, R_ud, R_du, T_dd = split_matrix(s_0)
        r_01, t_01, r_10, t_10 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select,
                                                                                     T_dd), select_modes(n_terms,
                                                                                                         n_select,
                                                                                                         R_du), select_modes(
            n_terms, n_select, T_uu)

        s_1 = get_s(t_incoherent, t1, p1)
        s_2 = get_s(t1, t2, p2)
        s_3 = get_s(t2, t3, p3)
        s_boundary = get_s_boundary(t3, t4)

        s = join_s_matrices_recur(s_1, s_2)
        s = join_s_matrices_recur(s, s_3)
        T_uu, R_ud, R_du, T_dd = join_s_matrices(s, s_boundary)
        r_12, t_12 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select, T_dd)

        a = prop_trans(n_select)

        # Propagační matice pro každý mód
        p_13 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)
        p_24 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)

        # Tvorba koherenčních matic
        cr_01 = np.kron(r_01, np.conj(r_01))
        ct_01 = np.kron(t_01, np.conj(t_01))
        cr_10 = np.kron(r_10, np.conj(r_10))
        ct_10 = np.kron(t_10, np.conj(t_10))
        cr_12 = np.kron(r_12, np.conj(r_12))
        ct_12 = np.kron(t_12, np.conj(t_12))
        cp_13 = np.kron(p_13, np.conj(p_13))
        cp_24 = np.kron(p_24, np.conj(p_24))

        coh_reflektance = cr_01 + ct_10 * cp_24 * cr_12 * cp_13 * (
                np.matrix(np.identity((2 * (2 * n_select + 1)) ** 2)) - cr_10 * cp_24 * cr_12 * cp_13) ** -1 * ct_01

        order = 0
        coh_matrix_0 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air, wv,
                                                                                              period)
        order = 1
        if (n_select > 0):
            coh_matrix_1 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_1 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
        order = 2
        if (n_select > 1):
            coh_matrix_2 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_2 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

        return float(np.real(coh_matrix_0[0, 0])), float(np.real(coh_matrix_0[3, 3])), float(
            np.real(coh_matrix_1[0, 0])), float(np.real(coh_matrix_1[3, 3]))

    values = np.array(list(map(periodDep_incoh, points)))

    toSend = {'data': []}

    for x in range(len(points)):
        tempDict = {'point': points[x], 'Rs0': values[x, 0], 'Rp0': values[x, 1]}
        toSend['data'].append(tempDict)

    fig1 = plt.figure(1, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 0], 'r-', label='$R_s$ 0 order')
    plt.plot(points, values[:, 1], 'b-', label='$R_p$ 0 order')
    plt.xlabel('Period [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot1'] = mpld3.fig_to_dict(fig1)

    fig2 = plt.figure(2, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 2], 'r-', label='$R_s$ -1 order')
    plt.plot(points, values[:, 3], 'b-', label='$R_p$ -1 order')
    plt.xlabel('Period [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot2'] = mpld3.fig_to_dict(fig2)

    print('start sending', flush=True)
    print(json.dumps(toSend), flush=True)
    print('end sending', flush=True)


def get_incoh_ref_thick(options):
    print('tick_flooper', flush=True)
    n_air = 1

    n_capping = complex(options['n_capping'])
    n_grating = complex(options['n_grating'])
    n_reflective = complex(options['n_reflective'])

    width_capping = float(options['width_capping'])

    angle = float(options['angle'])
    wv = float(options['wavelength'])
    period = float(options['grating_period'])
    fill_factor = float(options['fill_factor'])

    divisions = int(options['divisions'])
    start = float(options['divisions_start'])
    end = float(options['divisions_end'])

    perm_grating = np.array([[n_grating ** 2, 0, 0], [0, n_grating ** 2, 0], [0, 0, n_grating ** 2]])
    perm_reflective = np.array([[n_reflective ** 2, 0, 0], [0, n_reflective ** 2, 0], [0, 0, n_reflective ** 2]])
    perm_capping = np.array([[n_capping ** 2, 0, 0], [0, n_capping ** 2, 0], [0, 0, n_capping ** 2]])

    # calculate localy just once
    print('tick_flooper', flush=True)
    q_capping = get_Q_matrix(perm_capping, perm_capping, period, n_terms, 1)
    print('tick_flooper', flush=True)
    q_cap_ref = get_Q_matrix(perm_capping, perm_reflective, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_cap_grating = get_Q_matrix(perm_capping, perm_grating, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_ref_grating = get_Q_matrix(perm_reflective, perm_grating, period, n_terms, fill_factor)

    points = np.linspace(start, end, divisions)

    # Define the function
    def thickDep_incoh(width_ref):
        width_1 = width_ref
        width_2 = float(options['width_grating']) - width_ref
        width_3 = width_ref
        print('tick_flooper', flush=True)
        t0 = get_outer_T(n_terms, n_air, angle, wv, period, n_air)
        t_incoherent = get_outer_T(n_terms, n_capping, angle, wv, period, n_air)

        t1, p1 = get_periodic_layer_T_P_aniso(period, q_cap_ref[0], q_cap_ref[1], q_cap_ref[2], q_cap_ref[3],
                                              q_cap_ref[4], q_cap_ref[5], q_cap_ref[6], q_cap_ref[7], q_cap_ref[8],
                                              n_terms, angle, 0, wv, n_air, width_1)
        t2, p2 = get_periodic_layer_T_P_aniso(period, q_cap_grating[0], q_cap_grating[1], q_cap_grating[2],
                                              q_cap_grating[3], q_cap_grating[4], q_cap_grating[5], q_cap_grating[6],
                                              q_cap_grating[7], q_cap_grating[8], n_terms, angle, 0, wv, n_air, width_2)
        t3, p3 = get_periodic_layer_T_P_aniso(period, q_ref_grating[0], q_ref_grating[1], q_ref_grating[2],
                                              q_ref_grating[3], q_ref_grating[4], q_ref_grating[5], q_ref_grating[6],
                                              q_ref_grating[7], q_ref_grating[8], n_terms, angle, 0, wv, n_air, width_3)
        t4 = get_outer_T(n_terms, n_grating, angle, wv, period, n_air)

        s_0 = get_s_boundary(t0, t_incoherent)
        T_uu, R_ud, R_du, T_dd = split_matrix(s_0)
        r_01, t_01, r_10, t_10 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select,
                                                                                     T_dd), select_modes(n_terms,
                                                                                                         n_select,
                                                                                                         R_du), select_modes(
            n_terms, n_select, T_uu)

        s_1 = get_s(t_incoherent, t1, p1)
        s_2 = get_s(t1, t2, p2)
        s_3 = get_s(t2, t3, p3)
        s_boundary = get_s_boundary(t3, t4)

        s = join_s_matrices_recur(s_1, s_2)
        s = join_s_matrices_recur(s, s_3)
        T_uu, R_ud, R_du, T_dd = join_s_matrices(s, s_boundary)
        r_12, t_12 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select, T_dd)

        a = prop_trans(n_select)

        # Propagační matice pro každý mód
        p_13 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)
        p_24 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)

        # Tvorba koherenčních matic
        cr_01 = np.kron(r_01, np.conj(r_01))
        ct_01 = np.kron(t_01, np.conj(t_01))
        cr_10 = np.kron(r_10, np.conj(r_10))
        ct_10 = np.kron(t_10, np.conj(t_10))
        cr_12 = np.kron(r_12, np.conj(r_12))
        ct_12 = np.kron(t_12, np.conj(t_12))
        cp_13 = np.kron(p_13, np.conj(p_13))
        cp_24 = np.kron(p_24, np.conj(p_24))

        coh_reflektance = cr_01 + ct_10 * cp_24 * cr_12 * cp_13 * (
                np.matrix(np.identity((2 * (2 * n_select + 1)) ** 2)) - cr_10 * cp_24 * cr_12 * cp_13) ** -1 * ct_01

        order = 0
        coh_matrix_0 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air, wv,
                                                                                              period)
        order = 1
        if (n_select > 0):
            coh_matrix_1 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_1 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
        order = 2
        if (n_select > 1):
            coh_matrix_2 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_2 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

        return float(np.real(coh_matrix_0[0, 0])), float(np.real(coh_matrix_0[3, 3])), float(
            np.real(coh_matrix_1[0, 0])), float(np.real(coh_matrix_1[3, 3]))

    values = np.array(list(map(thickDep_incoh, points)))

    toSend = {'data': []}

    for x in range(len(points)):
        tempDict = {'point': points[x], 'Rs0': values[x, 0], 'Rp0': values[x, 1]}
        toSend['data'].append(tempDict)

    fig1 = plt.figure(1, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 0], 'r-', label='$R_s$ 0 order')
    plt.plot(points, values[:, 1], 'b-', label='$R_p$ 0 order')
    plt.xlabel('Reflective layer thickness [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot1'] = mpld3.fig_to_dict(fig1)

    fig2 = plt.figure(2, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 2], 'r-', label='$R_s$ -1 order')
    plt.plot(points, values[:, 3], 'b-', label='$R_p$ -1 order')
    plt.xlabel('Reflective layer thickness [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot2'] = mpld3.fig_to_dict(fig2)

    print('start sending', flush=True)
    print(json.dumps(toSend), flush=True)
    print('end sending', flush=True)


def get_incoh_cap_thick(options):
    print('tick_flooper', flush=True)
    n_air = 1

    n_capping = complex(options['n_capping'])
    n_grating = complex(options['n_grating'])
    n_reflective = complex(options['n_reflective'])

    width_1 = float(options['width_reflective'])
    width_2 = float(options['width_grating']) - float(options['width_reflective'])
    width_3 = float(options['width_reflective'])

    angle = float(options['angle'])
    wv = float(options['wavelength'])
    period = float(options['grating_period'])
    fill_factor = float(options['fill_factor'])

    divisions = int(options['divisions'])
    start = float(options['divisions_start'])
    end = float(options['divisions_end'])

    perm_grating = np.array([[n_grating ** 2, 0, 0], [0, n_grating ** 2, 0], [0, 0, n_grating ** 2]])
    perm_reflective = np.array([[n_reflective ** 2, 0, 0], [0, n_reflective ** 2, 0], [0, 0, n_reflective ** 2]])
    perm_capping = np.array([[n_capping ** 2, 0, 0], [0, n_capping ** 2, 0], [0, 0, n_capping ** 2]])

    # calculate localy just once
    print('tick_flooper', flush=True)
    q_capping = get_Q_matrix(perm_capping, perm_capping, period, n_terms, 1)
    print('tick_flooper', flush=True)
    q_cap_ref = get_Q_matrix(perm_capping, perm_reflective, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_cap_grating = get_Q_matrix(perm_capping, perm_grating, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_ref_grating = get_Q_matrix(perm_reflective, perm_grating, period, n_terms, fill_factor)

    points = np.linspace(start, end, divisions)

    # Define the function
    def thickDep_incoh(width_capping):
        print('tick_flooper', flush=True)
        t0 = get_outer_T(n_terms, n_air, angle, wv, period, n_air)
        t_incoherent = get_outer_T(n_terms, n_capping, angle, wv, period, n_air)

        t1, p1 = get_periodic_layer_T_P_aniso(period, q_cap_ref[0], q_cap_ref[1], q_cap_ref[2], q_cap_ref[3],
                                              q_cap_ref[4], q_cap_ref[5], q_cap_ref[6], q_cap_ref[7], q_cap_ref[8],
                                              n_terms, angle, 0, wv, n_air, width_1)
        t2, p2 = get_periodic_layer_T_P_aniso(period, q_cap_grating[0], q_cap_grating[1], q_cap_grating[2],
                                              q_cap_grating[3], q_cap_grating[4], q_cap_grating[5], q_cap_grating[6],
                                              q_cap_grating[7], q_cap_grating[8], n_terms, angle, 0, wv, n_air, width_2)
        t3, p3 = get_periodic_layer_T_P_aniso(period, q_ref_grating[0], q_ref_grating[1], q_ref_grating[2],
                                              q_ref_grating[3], q_ref_grating[4], q_ref_grating[5], q_ref_grating[6],
                                              q_ref_grating[7], q_ref_grating[8], n_terms, angle, 0, wv, n_air, width_3)
        t4 = get_outer_T(n_terms, n_grating, angle, wv, period, n_air)

        s_0 = get_s_boundary(t0, t_incoherent)
        T_uu, R_ud, R_du, T_dd = split_matrix(s_0)
        r_01, t_01, r_10, t_10 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select,
                                                                                     T_dd), select_modes(n_terms,
                                                                                                         n_select,
                                                                                                         R_du), select_modes(
            n_terms, n_select, T_uu)

        s_1 = get_s(t_incoherent, t1, p1)
        s_2 = get_s(t1, t2, p2)
        s_3 = get_s(t2, t3, p3)
        s_boundary = get_s_boundary(t3, t4)

        s = join_s_matrices_recur(s_1, s_2)
        s = join_s_matrices_recur(s, s_3)
        T_uu, R_ud, R_du, T_dd = join_s_matrices(s, s_boundary)
        r_12, t_12 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select, T_dd)

        a = prop_trans(n_select)

        # Propagační matice pro každý mód
        p_13 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)
        p_24 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)

        # Tvorba koherenčních matic
        cr_01 = np.kron(r_01, np.conj(r_01))
        ct_01 = np.kron(t_01, np.conj(t_01))
        cr_10 = np.kron(r_10, np.conj(r_10))
        ct_10 = np.kron(t_10, np.conj(t_10))
        cr_12 = np.kron(r_12, np.conj(r_12))
        ct_12 = np.kron(t_12, np.conj(t_12))
        cp_13 = np.kron(p_13, np.conj(p_13))
        cp_24 = np.kron(p_24, np.conj(p_24))

        coh_reflektance = cr_01 + ct_10 * cp_24 * cr_12 * cp_13 * (
                np.matrix(np.identity((2 * (2 * n_select + 1)) ** 2)) - cr_10 * cp_24 * cr_12 * cp_13) ** -1 * ct_01

        order = 0
        coh_matrix_0 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air, wv,
                                                                                              period)
        order = 1
        if (n_select > 0):
            coh_matrix_1 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_1 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
        order = 2
        if (n_select > 1):
            coh_matrix_2 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_2 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

        return float(np.real(coh_matrix_0[0, 0])), float(np.real(coh_matrix_0[3, 3])), float(
            np.real(coh_matrix_1[0, 0])), float(np.real(coh_matrix_1[3, 3]))

    values = np.array(list(map(thickDep_incoh, points)))

    toSend = {'data': []}

    for x in range(len(points)):
        tempDict = {'point': points[x], 'Rs0': values[x, 0], 'Rp0': values[x, 1]}
        toSend['data'].append(tempDict)

    fig1 = plt.figure(1, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 0], 'r-', label='$R_s$ 0 order')
    plt.plot(points, values[:, 1], 'b-', label='$R_p$ 0 order')
    plt.xlabel('Capping thickness [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot1'] = mpld3.fig_to_dict(fig1)

    fig2 = plt.figure(2, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 2], 'r-', label='$R_s$ -1 order')
    plt.plot(points, values[:, 3], 'b-', label='$R_p$ -1 order')
    plt.xlabel('Capping thickness [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot2'] = mpld3.fig_to_dict(fig2)

    print('start sending', flush=True)
    print(json.dumps(toSend), flush=True)
    print('end sending', flush=True)


def get_incoh_grating_thick(options):
    print('tick_flooper', flush=True)
    n_air = 1

    n_capping = complex(options['n_capping'])
    n_grating = complex(options['n_grating'])
    n_reflective = complex(options['n_reflective'])

    width_capping = float(options['width_capping'])

    angle = float(options['angle'])
    wv = float(options['wavelength'])
    period = float(options['grating_period'])
    fill_factor = float(options['fill_factor'])

    divisions = int(options['divisions'])
    start = float(options['divisions_start'])
    end = float(options['divisions_end'])

    perm_grating = np.array([[n_grating ** 2, 0, 0], [0, n_grating ** 2, 0], [0, 0, n_grating ** 2]])
    perm_reflective = np.array([[n_reflective ** 2, 0, 0], [0, n_reflective ** 2, 0], [0, 0, n_reflective ** 2]])
    perm_capping = np.array([[n_capping ** 2, 0, 0], [0, n_capping ** 2, 0], [0, 0, n_capping ** 2]])

    # calculate localy just once
    print('tick_flooper', flush=True)
    q_capping = get_Q_matrix(perm_capping, perm_capping, period, n_terms, 1)
    print('tick_flooper', flush=True)
    q_cap_ref = get_Q_matrix(perm_capping, perm_reflective, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_cap_grating = get_Q_matrix(perm_capping, perm_grating, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_ref_grating = get_Q_matrix(perm_reflective, perm_grating, period, n_terms, fill_factor)

    points = np.linspace(start, end, divisions)

    # Define the function
    def thickDep_incoh(width_grating):
        width_1 = float(options['width_reflective'])
        width_2 = width_grating - float(options['width_reflective'])
        width_3 = float(options['width_reflective'])

        print('tick_flooper', flush=True)
        t0 = get_outer_T(n_terms, n_air, angle, wv, period, n_air)
        t_incoherent = get_outer_T(n_terms, n_capping, angle, wv, period, n_air)

        t1, p1 = get_periodic_layer_T_P_aniso(period, q_cap_ref[0], q_cap_ref[1], q_cap_ref[2], q_cap_ref[3],
                                              q_cap_ref[4], q_cap_ref[5], q_cap_ref[6], q_cap_ref[7], q_cap_ref[8],
                                              n_terms, angle, 0, wv, n_air, width_1)
        t2, p2 = get_periodic_layer_T_P_aniso(period, q_cap_grating[0], q_cap_grating[1], q_cap_grating[2],
                                              q_cap_grating[3], q_cap_grating[4], q_cap_grating[5], q_cap_grating[6],
                                              q_cap_grating[7], q_cap_grating[8], n_terms, angle, 0, wv, n_air, width_2)
        t3, p3 = get_periodic_layer_T_P_aniso(period, q_ref_grating[0], q_ref_grating[1], q_ref_grating[2],
                                              q_ref_grating[3], q_ref_grating[4], q_ref_grating[5], q_ref_grating[6],
                                              q_ref_grating[7], q_ref_grating[8], n_terms, angle, 0, wv, n_air, width_3)
        t4 = get_outer_T(n_terms, n_grating, angle, wv, period, n_air)

        s_0 = get_s_boundary(t0, t_incoherent)
        T_uu, R_ud, R_du, T_dd = split_matrix(s_0)
        r_01, t_01, r_10, t_10 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select,
                                                                                     T_dd), select_modes(n_terms,
                                                                                                         n_select,
                                                                                                         R_du), select_modes(
            n_terms, n_select, T_uu)

        s_1 = get_s(t_incoherent, t1, p1)
        s_2 = get_s(t1, t2, p2)
        s_3 = get_s(t2, t3, p3)
        s_boundary = get_s_boundary(t3, t4)

        s = join_s_matrices_recur(s_1, s_2)
        s = join_s_matrices_recur(s, s_3)
        T_uu, R_ud, R_du, T_dd = join_s_matrices(s, s_boundary)
        r_12, t_12 = select_modes(n_terms, n_select, R_ud), select_modes(n_terms, n_select, T_dd)

        a = prop_trans(n_select)

        # Propagační matice pro každý mód
        p_13 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)
        p_24 = get_prop_matrix(period, angle, n_capping, n_air, n_select, width_capping, wv)

        # Tvorba koherenčních matic
        cr_01 = np.kron(r_01, np.conj(r_01))
        ct_01 = np.kron(t_01, np.conj(t_01))
        cr_10 = np.kron(r_10, np.conj(r_10))
        ct_10 = np.kron(t_10, np.conj(t_10))
        cr_12 = np.kron(r_12, np.conj(r_12))
        ct_12 = np.kron(t_12, np.conj(t_12))
        cp_13 = np.kron(p_13, np.conj(p_13))
        cp_24 = np.kron(p_24, np.conj(p_24))

        coh_reflektance = cr_01 + ct_10 * cp_24 * cr_12 * cp_13 * (
                np.matrix(np.identity((2 * (2 * n_select + 1)) ** 2)) - cr_10 * cp_24 * cr_12 * cp_13) ** -1 * ct_01

        order = 0
        coh_matrix_0 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air, wv,
                                                                                              period)
        order = 1
        if (n_select > 0):
            coh_matrix_1 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_1 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
        order = 2
        if (n_select > 1):
            coh_matrix_2 = get_order_coh_matrix(coh_reflektance, n_select, order) * diff_eff_coef(order, angle, n_air,
                                                                                                  wv, period)
        else:
            coh_matrix_2 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

        return float(np.real(coh_matrix_0[0, 0])), float(np.real(coh_matrix_0[3, 3])), float(
            np.real(coh_matrix_1[0, 0])), float(np.real(coh_matrix_1[3, 3]))

    values = np.array(list(map(thickDep_incoh, points)))

    toSend = {'data': []}

    for x in range(len(points)):
        tempDict = {'point': points[x], 'Rs0': values[x, 0], 'Rp0': values[x, 1]}
        toSend['data'].append(tempDict)

    fig1 = plt.figure(1, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 0], 'r-', label='$R_s$ 0 order')
    plt.plot(points, values[:, 1], 'b-', label='$R_p$ 0 order')
    plt.xlabel('Grating thickness [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot1'] = mpld3.fig_to_dict(fig1)

    fig2 = plt.figure(2, figsize=(4.5, 3.7))
    plt.title("Reflectivity of the zero difrraction order")
    plt.plot(points, values[:, 2], 'r-', label='$R_s$ -1 order')
    plt.plot(points, values[:, 3], 'b-', label='$R_p$ -1 order')
    plt.xlabel('Grating thickness [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend['plot2'] = mpld3.fig_to_dict(fig2)

    print('start sending', flush=True)
    print(json.dumps(toSend), flush=True)
    print('end sending', flush=True)


data = {
    'n_capping': '1.3',
    'n_grating': '1.5',
    'n_reflective': '1-6j',
    'width_capping': '2000',
    'width_reflective': '20',
    'width_grating': '100',
    'wavelength': '633',
    'angle': '0',
    'grating_period': '1000',
    'fill_factor': '0.5',
    'divisions': '5',
    'divisions_start': '100',
    'divisions_end': '1000',
}

if __name__ == '__main__':
    show_figures = True
    get_incoh_grating_thick(data)
