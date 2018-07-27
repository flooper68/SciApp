import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import mpld3
import numpy as np
import json


class Message:
    options = {}
    structure = {}
    name = ''
    graphs = {}
    timestamp = ''
    data = []
    duration = ''
    _id = ''

    def send(self):
        print(json.dumps({
            'options': self.options,
            'structure': self.structure,
            'graphs': self.graphs,
            'timestamp': self.timestamp,
            'data': self.data,
            'duration': self.duration,
            'name': self.name,
            '_id': self._id
        }))


def show_specular_diffraction_ref(values, options, toSend={}, show=True):
    values = np.real(values)
    specular_ref = np.real((values[:, 0, 0] + values[:, 1, 0]) / 2)

    diffraction_ref = np.real((values[:, 0].sum(axis=1) + values[:, 1].sum(axis=1) - values[:, 0, 0] - values[:, 1, 0]) / 4)

    points = np.linspace(options['divisions_start'], options['divisions_end'], options['divisions'])

    fig1 = plt.figure(1, figsize=(4.5, 3.7))
    plt.title("Total diffraction efficiency R")
    plt.plot(points, specular_ref, 'r-', label='Specular')
    plt.plot(points, diffraction_ref, 'b-', label='Diffraction')
    plt.xlabel(options['dependence']['label'])
    max = np.real(np.array([specular_ref, diffraction_ref]).max())
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4 * max])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()
    toSend.graphs['specular_diff_ref'] = mpld3.fig_to_dict(fig1)


def show_specular_diffraction_trans(values, options, toSend={}, show=True):
    values = np.real(values)
    specular_trans = np.real((values[:, 2, 0] + values[:, 3, 0]) / 2)

    diffraction_trans = np.real((values[:, 2].sum(axis=1) + values[:, 3].sum(axis=1) - values[:, 2, 0] - values[:, 3, 0]) / 4)

    points = np.linspace(options['divisions_start'], options['divisions_end'], options['divisions'])

    fig2 = plt.figure(2, figsize=(4.5, 3.7))
    plt.title("Total diffraction efficiency T")
    plt.plot(points, specular_trans, 'r-', label='Specular')
    plt.plot(points, diffraction_trans, 'b-', label='Diffraction')
    plt.xlabel(options['dependence']['label'])
    max = np.real(np.array([specular_trans, diffraction_trans]).max())
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4 * max])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()
    toSend.graphs['specular_diff_trans'] = mpld3.fig_to_dict(fig2)


def show_total(values, options, toSend={}, show=True):
    values = np.real(values)
    total_trans = (values[:, 0].sum(axis=1) + values[:, 1].sum(axis=1)) / 2
    total_ref = (values[:, 2].sum(axis=1) + values[:, 3].sum(axis=1)) / 2

    points = np.linspace(options['divisions_start'], options['divisions_end'], options['divisions'])

    fig3 = plt.figure(3, figsize=(4.5, 3.7))
    plt.title("Total transmission and reflection")
    plt.plot(points, total_trans, 'r-', label='Total transmission')
    plt.plot(points, total_ref, 'b-', label='Total reflection ')
    plt.xlabel(options['dependence']['label'])
    max = np.array([total_trans, total_ref]).max()
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4 * max])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()

    fig4 = plt.figure(4, figsize=(4.5, 3.7))
    plt.title("Total transmission and reflection and absorbance")
    plt.plot(points, total_ref + total_trans, 'g-', label='Total transmission + reflection')
    plt.plot(points, 1 - total_ref - total_trans, 'y-', label='Absorbance')
    plt.xlabel(options['dependence']['label'])
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()

    toSend.graphs['total_trans_ref'] = mpld3.fig_to_dict(fig3)
    toSend.graphs['total_abs'] = mpld3.fig_to_dict(fig4)


def show_first_negative_orders_trans_pol(values, options, toSend={}, show=True):
    values = np.real(values)
    points = np.linspace(options['divisions_start'], options['divisions_end'], options['divisions'])

    fig5 = plt.figure(5, figsize=(4.5, 3.7))
    plt.title("Diffraction efficiency Ts, Tp")
    plt.plot(points, values[:, 2, 2], 'r-', label='Ts')
    plt.plot(points, values[:, 3, 2], 'b-', label='Tp')
    plt.xlabel(options['dependence']['label'])
    plt.ylabel('-1st diffraction order')
    max = np.array([values[:, 2, 2], values[:, 3, 2]]).max()
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4 * max])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()
    toSend.graphs['first_negative_trans_pol'] = mpld3.fig_to_dict(fig5)


def show_first_negative_orders_trans_unpol(values, options, toSend={}, show=True):
    values = np.real(values)
    points = np.linspace(options['divisions_start'], options['divisions_end'], options['divisions'])

    fig6 = plt.figure(6, figsize=(4.5, 3.7))
    plt.title("Diffraction efficiency T")
    plt.plot(points, (values[:, 2, 2] + values[:, 3, 2]) / 2, 'r-', label='T')
    plt.xlabel(options['dependence']['label'])
    plt.ylabel('-1st diffraction order')
    max = np.array((values[:, 2, 2] + values[:, 3, 2]) / 2).max()
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4 * max])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()

    toSend.graphs['first_negative_trans_unpol'] = mpld3.fig_to_dict(fig6)


def show_first_negative_orders_ref_unpol(values, options, toSend={}, show=True):
    values = np.real(values)
    points = np.linspace(options['divisions_start'], options['divisions_end'], options['divisions'])

    fig7 = plt.figure(7, figsize=(4.5, 3.7))
    plt.title("Diffraction efficiency R")
    plt.plot(points, (values[:, 0, 2] + values[:, 1, 2]) / 2, 'r-', label='R')
    plt.xlabel(options['dependence']['label'])
    plt.ylabel('-1st diffraction order')
    max = np.array((values[:, 0, 2] + values[:, 1, 2]) / 2).max()
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4 * max])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()

    toSend.graphs['first_negative_ref_unpol'] = mpld3.fig_to_dict(fig7)


def show_first_negative_orders_ref_pol(values, options, toSend={}, show=True):
    values = np.real(values)
    points = np.linspace(options['divisions_start'], options['divisions_end'], options['divisions'])

    fig8 = plt.figure(8, figsize=(4.5, 3.7))
    plt.title("Diffraction efficiency Rs, Rp")
    plt.plot(points, values[:, 0, 2], 'r-', label='Rs')
    plt.plot(points, values[:, 1, 2], 'b-', label='Rp')
    plt.xlabel(options['dependence']['label'])
    plt.ylabel('-1st diffraction order')
    max = np.array([values[:, 0, 2], values[:, 1, 2]]).max()
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4 * max])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()

    toSend.graphs['first_negative_ref_pol'] = mpld3.fig_to_dict(fig8)


def show_diff_angles(values, options, toSend={}, show=True):
    points = np.linspace(options['divisions_start'], options['divisions_end'], options['divisions'])

    values = np.real(values)
    fig9 = plt.figure(9, figsize=(4.5, 3.7))
    plt.title("Diffraction angles")
    plt.plot(points, values[:, 4], label='-2 order')
    plt.plot(points, values[:, 2], label='-1 order')
    plt.plot(points, values[:, 0], label='0 order')
    plt.plot(points, values[:, 1], label='1 order')
    plt.plot(points, values[:, 3], label='2 order')
    plt.xlabel('Angle [deg]')
    plt.ylabel('Diffraction angles [deg')
    plt.legend(loc='upper right')
    plt.axis([options['divisions_start'], options['divisions_end'], -100, 100])
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()
    toSend.graphs['diff_angles'] = mpld3.fig_to_dict(fig9)


def show_comprehensive(values, options, toSend={}, show=True):
    values = np.real(values)
    points = np.linspace(options['divisions_start'], options['divisions_end'], options['divisions'])

    fig1 = plt.figure(10, figsize=(4.5, 3.7))
    plt.title("Zero orders - Reflection")
    plt.plot(points, values[:, 0, 0], 'r-', label='Rs')
    plt.plot(points, values[:, 1, 0], 'b-', label='Rp')
    plt.xlabel(options['dependence']['label'])
    max = np.array([values[:, 0, 0], values[:, 1, 0]]).max()
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4 * max])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()

    fig2 = plt.figure(11, figsize=(4.5, 3.7))
    plt.title("Negative orders - Reflection")
    plt.plot(points, values[:, 0, 2], 'r-', label='Rs -1')
    plt.plot(points, values[:, 0, 4], 'r--', label='Rs -2')
    plt.plot(points, values[:, 0, 6], 'r-.', label='Rs -3')
    plt.plot(points, values[:, 1, 2], 'b-', label='Rp -1')
    plt.plot(points, values[:, 1, 4], 'b--', label='Rp -2')
    plt.plot(points, values[:, 1, 6], 'b-.', label='Rp -3')
    plt.xlabel(options['dependence']['label'])
    max = np.array([values[:, 0, 2], values[:, 0, 4], values[:, 0, 6], values[:, 1, 2], values[:, 1, 4], values[:, 1, 6]]).max()
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4 * max])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()

    fig3 = plt.figure(12, figsize=(4.5, 3.7))
    plt.title("Positive orders - Reflection")
    plt.plot(points, values[:, 0, 1], 'r-', label='Rs 1')
    plt.plot(points, values[:, 0, 3], 'r--', label='Rs 2')
    plt.plot(points, values[:, 0, 5], 'r-.', label='Rs 3')
    plt.plot(points, values[:, 1, 1], 'b-', label='Rp 1')
    plt.plot(points, values[:, 1, 3], 'b--', label='Rp 2')
    plt.plot(points, values[:, 1, 5], 'b-.', label='Rp 3')
    plt.xlabel(options['dependence']['label'])
    max = np.array([values[:, 0, 1], values[:, 0, 3], values[:, 0, 5], values[:, 1, 1], values[:, 1, 3], values[:, 1, 5]]).max()
    plt.axis([options['divisions_start'], options['divisions_end'], 0, 1.4 * max])
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()

    toSend.graphs['comp_zero'] = mpld3.fig_to_dict(fig1)
    toSend.graphs['comp_neg'] = mpld3.fig_to_dict(fig2)
    toSend.graphs['comp_pos'] = mpld3.fig_to_dict(fig3)
