from python.funkce import *
from python.zobrazovaci_funkce import *

perm_ref = (-35 - 12j)
perm_capping = 2.25
param = 0.5

structure = {
    'bound_selectors': [
        # ['data', 0, 'permittivity', 4],
        # ['data', 0, 'permittivity', 8],
    ],
    'superstrate': 1,  # indexes, not perm
    'substrate': 1.5,
    'period': 1000,
    'data': [
        {'width': 2000, 'periodic': False, 'coherent': True,
         'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]},
        {'width': 20, 'periodic': True, 'coherent': True, 'materials': [
            {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': param},
            {'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping], 'start': param, 'stop': 1},
        ]},
        {'width': 80, 'periodic': True, 'coherent': True, 'materials': [
            {'permittivity': [2.25, 0, 0, 0, 2.25, 0, 0, 0, 2.25], 'start': 0, 'stop': param},
            {'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping], 'start': param, 'stop': 1},
        ]},
        {'width': 20, 'periodic': True, 'coherent': True, 'materials': [
            {'permittivity': [2.25, 0, 0, 0, 2.25, 0, 0, 0, 2.25], 'start': 0, 'stop': param},
            {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': param, 'stop': 1},
        ]}
    ]}

options = {
    'bound_selectors': [
        ['wavelength'],  # Select the proper attribute to be dependant
    ],
    'n_terms': 10,
    'n_select': 5,
    'wavelength': 833,
    'angle': 0,
    'divisions': 40,
    'divisions_start': 400,
    'divisions_end': 800,
    'dependence': {
        'name': 'wavelength',  ##angle // wavelength // structure
        'label': 'Angle [deg]'
    }
}

values, angles = general_func(structure_defining_func(structure), options_defining_func(options))

options = options_defining_func(options)(None)
show_specular_diffraction_ref(values, options)
