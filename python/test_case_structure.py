import json
import math
import time

import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import mpld3
from funkce import *

n_terms = 5;

data = {
    'n_capping': '1',
    'n_grating': '1.5',
    'n_reflective': '1-6j',
    'width_capping': '0',
    'width_reflective': '20',
    'width_grating': '100',
    'wavelength': '633',
    'angle': '-20',
    'grating_period': '1000',
    'fill_factor': '0.5',
    'divisions': '100',
    'divisions_start': '400',
    'divisions_end': '800',
}


class Message:
    inputs = {}
    name = ''
    graphs = {}
    timestamp = ''
    data = []
    duration = ''

    def send(self):
        print(json.dumps({
            'parameters': self.inputs,
            'name': self.name,
            'graphs': self.graphs,
            'timestamp': self.timestamp,
            'data': self.data,
            'duration': self.duration
        }))


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


def get_coh_wv(options):
    startTime = time.time()
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

    print('tick_flooper', flush=True)
    # calculate localy just once
    q_capping = get_Q_matrix(perm_capping, perm_capping, period, n_terms, 1)
    print('tick_flooper', flush=True)
    q_cap_ref = get_Q_matrix(perm_capping, perm_reflective, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_cap_grating = get_Q_matrix(perm_capping, perm_grating, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)
    q_ref_grating = get_Q_matrix(perm_reflective, perm_grating, period, n_terms, fill_factor)
    print('tick_flooper', flush=True)

    points = np.linspace(start, end, divisions)

    # Define the function
    def wvDep(wv):
        print('tick_flooper', flush=True)
        t0 = get_outer_T(n_terms, n_air, angle, wv, period, n_air)
        t_incoherent, p_incoherent = get_periodic_layer_T_P_aniso(period, q_capping[0], q_capping[1], q_capping[2],
                                                                  q_capping[3], q_capping[4], q_capping[5],
                                                                  q_capping[6], q_capping[7], q_capping[8], n_terms,
                                                                  angle, 0, wv, n_air, width_capping)
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

        s_0 = get_s(t0, t_incoherent, p_incoherent)
        s_1 = get_s(t_incoherent, t1, p1)
        s_2 = get_s(t1, t2, p2)
        s_3 = get_s(t2, t3, p3)
        s_boundary = get_s_boundary(t3, t4)

        s = join_s_matrices_recur(s_0, s_1)
        s = join_s_matrices_recur(s, s_2)
        s = join_s_matrices_recur(s, s_3)
        T_uu, R_ud, R_du, T_dd = join_s_matrices(s, s_boundary)

        Rs_order0 = np.abs(R_ud[n_terms, n_terms]) ** 2 * diff_eff_coef(0, angle, n_air, wv, period)
        Rp_order0 = np.abs(R_ud[3 * n_terms + 1, 3 * n_terms + 1]) ** 2 * diff_eff_coef(0, angle, n_air, wv, period)

        tsOrder0 = np.abs(T_dd[n_terms, n_terms]) ** 2 * diff_eff_coef_trans(0, angle, n_air, n_grating, wv, period)
        tpOrder0 = np.abs(T_dd[3 * n_terms + 1, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(0, angle, n_air, n_grating, wv, period)

        Rs_order_neg1 = np.abs(R_ud[n_terms - 1, n_terms]) ** 2 * diff_eff_coef(-1, angle, n_air, wv, period)
        Rs_order_neg2 = np.abs(R_ud[n_terms - 2, n_terms]) ** 2 * diff_eff_coef(-2, angle, n_air, wv, period)
        Rs_order_neg3 = np.abs(R_ud[n_terms - 3, n_terms]) ** 2 * diff_eff_coef(-3, angle, n_air, wv, period)
        Rs_order_neg4 = np.abs(R_ud[n_terms - 4, n_terms]) ** 2 * diff_eff_coef(-4, angle, n_air, wv, period)
        Rs_order_neg5 = np.abs(R_ud[n_terms - 5, n_terms]) ** 2 * diff_eff_coef(-5, angle, n_air, wv, period)

        tsOrderNeg1 = np.abs(T_dd[n_terms - 1, n_terms]) ** 2 * diff_eff_coef_trans(-1, angle, n_air, n_grating, wv, period)
        tsOrderNeg2 = np.abs(T_dd[n_terms - 2, n_terms]) ** 2 * diff_eff_coef_trans(-2, angle, n_air, n_grating, wv, period)
        tsOrderNeg3 = np.abs(T_dd[n_terms - 3, n_terms]) ** 2 * diff_eff_coef_trans(-3, angle, n_air, n_grating, wv, period)
        tsOrderNeg4 = np.abs(T_dd[n_terms - 4, n_terms]) ** 2 * diff_eff_coef_trans(-4, angle, n_air, n_grating, wv, period)
        tsOrderNeg5 = np.abs(T_dd[n_terms - 5, n_terms]) ** 2 * diff_eff_coef_trans(-5, angle, n_air, n_grating, wv, period)

        Rs_order_pos1 = np.abs(R_ud[n_terms + 1, n_terms]) ** 2 * diff_eff_coef(1, angle, n_air, wv, period)
        Rs_order_pos2 = np.abs(R_ud[n_terms + 2, n_terms]) ** 2 * diff_eff_coef(2, angle, n_air, wv, period)
        Rs_order_pos3 = np.abs(R_ud[n_terms + 3, n_terms]) ** 2 * diff_eff_coef(3, angle, n_air, wv, period)
        Rs_order_pos4 = np.abs(R_ud[n_terms + 4, n_terms]) ** 2 * diff_eff_coef(4, angle, n_air, wv, period)
        Rs_order_pos5 = np.abs(R_ud[n_terms + 5, n_terms]) ** 2 * diff_eff_coef(5, angle, n_air, wv, period)

        tsOrderPos1 = np.abs(T_dd[n_terms + 1, n_terms]) ** 2 * diff_eff_coef_trans(1, angle, n_air, n_grating, wv, period)
        tsOrderPos2 = np.abs(T_dd[n_terms + 2, n_terms]) ** 2 * diff_eff_coef_trans(2, angle, n_air, n_grating, wv, period)
        tsOrderPos3 = np.abs(T_dd[n_terms + 3, n_terms]) ** 2 * diff_eff_coef_trans(3, angle, n_air, n_grating, wv, period)
        tsOrderPos4 = np.abs(T_dd[n_terms + 4, n_terms]) ** 2 * diff_eff_coef_trans(4, angle, n_air, n_grating, wv, period)
        tsOrderPos5 = np.abs(T_dd[n_terms + 5, n_terms]) ** 2 * diff_eff_coef_trans(5, angle, n_air, n_grating, wv, period)

        Rp_order_neg1 = np.abs(R_ud[3 * n_terms, 3 * n_terms + 1]) ** 2 * diff_eff_coef(-1, angle, n_air, wv, period)
        Rp_order_neg2 = np.abs(R_ud[3 * n_terms - 1, 3 * n_terms + 1]) ** 2 * diff_eff_coef(-2, angle, n_air, wv, period)
        Rp_order_neg3 = np.abs(R_ud[3 * n_terms - 2, 3 * n_terms + 1]) ** 2 * diff_eff_coef(-3, angle, n_air, wv, period)
        Rp_order_neg4 = np.abs(R_ud[3 * n_terms - 3, 3 * n_terms + 1]) ** 2 * diff_eff_coef(-4, angle, n_air, wv, period)
        Rp_order_neg5 = np.abs(R_ud[3 * n_terms - 4, 3 * n_terms + 1]) ** 2 * diff_eff_coef(-5, angle, n_air, wv, period)

        tpOrderNeg1 = np.abs(T_dd[3 * n_terms, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(-1, angle, n_air, n_grating, wv, period)
        tpOrderNeg2 = np.abs(T_dd[3 * n_terms - 1, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(-2, angle, n_air, n_grating, wv, period)
        tpOrderNeg3 = np.abs(T_dd[3 * n_terms - 2, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(-3, angle, n_air, n_grating, wv, period)
        tpOrderNeg4 = np.abs(T_dd[3 * n_terms - 3, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(-4, angle, n_air, n_grating, wv, period)
        tpOrderNeg5 = np.abs(T_dd[3 * n_terms - 4, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(-5, angle, n_air, n_grating, wv, period)

        Rp_order_pos1 = np.abs(R_ud[3 * n_terms + 2, 3 * n_terms + 1]) ** 2 * diff_eff_coef(1, angle, n_air, wv, period)
        Rp_order_pos2 = np.abs(R_ud[3 * n_terms + 3, 3 * n_terms + 1]) ** 2 * diff_eff_coef(2, angle, n_air, wv, period)
        Rp_order_pos3 = np.abs(R_ud[3 * n_terms + 4, 3 * n_terms + 1]) ** 2 * diff_eff_coef(3, angle, n_air, wv, period)
        Rp_order_pos4 = np.abs(R_ud[3 * n_terms + 5, 3 * n_terms + 1]) ** 2 * diff_eff_coef(4, angle, n_air, wv, period)
        Rp_order_pos5 = np.abs(R_ud[3 * n_terms + 6, 3 * n_terms + 1]) ** 2 * diff_eff_coef(5, angle, n_air, wv, period)

        tpOrderPos1 = np.abs(T_dd[3 * n_terms + 2, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(1, angle, n_air, n_grating, wv, period)
        tpOrderPos2 = np.abs(T_dd[3 * n_terms + 3, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(2, angle, n_air, n_grating, wv, period)
        tpOrderPos3 = np.abs(T_dd[3 * n_terms + 4, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(3, angle, n_air, n_grating, wv, period)
        tpOrderPos4 = np.abs(T_dd[3 * n_terms + 5, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(4, angle, n_air, n_grating, wv, period)
        tpOrderPos5 = np.abs(T_dd[3 * n_terms + 6, 3 * n_terms + 1]) ** 2 * diff_eff_coef_trans(5, angle, n_air, n_grating, wv, period)

        totalRSpecular = (Rs_order0 + Rp_order0) / 2
        totalTSpecular = (tsOrder0 + tpOrder0) / 2
        totalRDiffraction = (
                                        Rs_order_neg1 + Rp_order_neg1 + Rs_order_neg2 + Rp_order_neg2 + Rs_order_neg3 + Rp_order_neg3 + Rs_order_neg4 + Rp_order_neg4 + Rs_order_neg5 + Rp_order_neg5 + Rs_order_pos1 + Rp_order_pos1 + Rs_order_pos2 + Rp_order_pos2 + Rs_order_pos3 + Rp_order_pos3 + Rs_order_pos4 + Rp_order_pos4 + Rs_order_pos5 + Rp_order_pos5) / 2
        totalTDiffraction = (
                                        tsOrderNeg1 + tpOrderNeg1 + tsOrderNeg2 + tpOrderNeg2 + tsOrderNeg3 + tpOrderNeg3 + tsOrderNeg4 + tpOrderNeg4 + tsOrderNeg5 + tpOrderNeg5 + tsOrderPos1 + tpOrderPos1 + tsOrderPos2 + tpOrderPos2 + tsOrderPos3 + tpOrderPos3 + tsOrderPos4 + tpOrderPos4 + tsOrderPos5 + tpOrderPos5) / 2

        total_diff_eff = totalRSpecular + totalTSpecular + totalRDiffraction + totalTDiffraction

        def get_diff_angle(order):
            temp = np.sin(np.deg2rad(angle)) + order * wv / period
            if ((np.rad2deg(temp) > -57.29 and np.rad2deg(temp) < 57.29) and np.rad2deg(temp) != 0):
                return np.rad2deg(np.arcsin(temp))
            elif (np.rad2deg(temp) < -57.29):
                return angle
            else:
                return angle

        return totalRSpecular, totalRDiffraction, Rs_order_neg1, Rp_order_neg1, totalTSpecular, totalTDiffraction, tsOrderNeg1, tpOrderNeg1, total_diff_eff, get_diff_angle(
            -2), get_diff_angle(-1), get_diff_angle(0), get_diff_angle(1), get_diff_angle(2)

    values = np.array(list(map(wvDep, points)))
    stopTime = time.time()

    toSend = Message()

    for x in range(len(points)):
        tempDict = {'point': points[x], 'totalRSpecular': values[x, 0], 'totalRDiffraction': values[x, 1]}
        toSend.data.append(tempDict)

    fig1 = plt.figure(1, figsize=(4.5, 3.7))
    plt.title("Total diffraction efficiency R")
    plt.plot(points, values[:, 0], 'r-', label='Specular')
    plt.plot(points, values[:, 1], 'b-', label='Diffraction')
    plt.xlabel('Wavelength [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend.graphs['totalR'] = mpld3.fig_to_dict(fig1)

    fig2 = plt.figure(2, figsize=(4.5, 3.7))
    plt.title("Diffraction efficiency - Rs, Rp")
    plt.plot(points, values[:, 2], 'r-', label='Rs')
    plt.plot(points, values[:, 3], 'b-', label='Rp')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('-1st diffraction order')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend.graphs['diffR'] = mpld3.fig_to_dict(fig2)

    fig3 = plt.figure(3, figsize=(4.5, 3.7))
    plt.title("Total diffraction efficiency T")
    plt.plot(points, values[:, 4], 'r-', label='Specular')
    plt.plot(points, values[:, 5], 'b-', label='Diffraction')
    plt.xlabel('Wavelength [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend.graphs['totalT'] = mpld3.fig_to_dict(fig3)

    fig5 = plt.figure(5, figsize=(4.5, 3.7))
    plt.title("Diffraction efficiency - Ts, Tp")
    plt.plot(points, values[:, 6], 'r-', label='Ts')
    plt.plot(points, values[:, 7], 'b-', label='Tp')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('-1st diffraction order')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend.graphs['diffT'] = mpld3.fig_to_dict(fig5)

    fig6 = plt.figure(6, figsize=(4.5, 3.7))
    plt.title("Diffraction angles")
    plt.plot(points, values[:, 9], label='-2 order')
    plt.plot(points, values[:, 10], label='-1 order')
    plt.plot(points, values[:, 11], label='0 order')
    plt.plot(points, values[:, 12], label='1 order')
    plt.plot(points, values[:, 13], label='2 order')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Diffraction angles [deg')
    plt.legend(loc='upper right')
    plt.axis([start, end, -100, 100])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend.graphs['diffAngles'] = mpld3.fig_to_dict(fig6)

    fig4 = plt.figure(4, figsize=(4.5, 3.7))
    plt.title("Total diffraction efficiency ")
    plt.plot(points, values[:, 8], 'r-', label='Total')
    plt.xlabel('Wavelength [nm]')
    plt.legend(loc='upper right')
    plt.axis([start, end, 0, 1])
    plt.grid(True)
    if show_figures:
        plt.show()

    toSend.graphs['total'] = mpld3.fig_to_dict(fig4)

    print('start sending', flush=True)
    toSend.duration = math.floor((stopTime - startTime))
    toSend.inputs = options
    toSend.timestamp = time.time()
    toSend.name = 'Calc' + str(toSend.timestamp)
    toSend.send()
    print('end sending', flush=True)

    get_coh_wv(data)
