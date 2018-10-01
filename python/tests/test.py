from concurrent.futures import ThreadPoolExecutor

from python.funkce import *
from python.zobrazovaci_funkce import *

pool = ThreadPoolExecutor(4)


def general_func(structure_defining_func=None, options_defining_func=None):
    options = options_defining_func(0)
    structure = structure_defining_func(0)

    def options_func(dep_options):
        print('timer: tick', flush=True)
        for item in coherent_subsystems:
            # Ke kazde podstrukture pridam jeji vypocitanou s-matici
            item['s_matrix'] = calculate_s_matrix(item, dep_options, item['toeplitz_matrix'])

        # Pomoci partial waves algorithm vypocitam vysledne koherencni matice
        coh_matrices = calculate_incoherent_partial_waves(coherent_subsystems, dep_options)

        reflectance_s, reflectance_p, transmittance_s, transmittance_p, diff_angles = get_observables(coh_matrices,
                                                                                                      dep_options,
                                                                                                      structure)
        # Z vypocitanych matic ziskam potrebne veliciny a vratim je
        return {'coefficients': np.array([reflectance_s, reflectance_p, transmittance_s, transmittance_p, ]),
                'diff_angles': diff_angles}

    def structure_func(structure):
        coherent_subsystems = get_coherent_subsets(structure)
        print('timer: tick', flush=True)
        for item in coherent_subsystems:
            # Kdyz neni toeplizova matice predvypocitana, tak ji vypocitam tady
            if not item.get('toeplitz_matrix'):
                item['toeplitz_matrix'] = calculate_toeplitz_matrices(item, options)
            # Ke kazde podstrukture pridam jeji vypocitanou s-matici
            item['s_matrix'] = calculate_s_matrix(item, options, item['toeplitz_matrix'])

        # Pomoci partial waves algorithm vypocitam vysledne koherencni matice
        coh_matrices = calculate_incoherent_partial_waves(coherent_subsystems, options)

        reflectance_s, reflectance_p, transmittance_s, transmittance_p, diff_angles = get_observables(coh_matrices,
                                                                                                      options,
                                                                                                      structure)
        # Z vypocitanych matic ziskam potrebne veliciny a vratim je
        return {'coefficients': np.array([reflectance_s, reflectance_p, transmittance_s, transmittance_p, ]),
                'diff_angles': diff_angles}

    # Ziskam options objekty pro danou zavislost
    points = np.linspace(options['divisions_start'], options['divisions_end'], options['divisions'])
    dependence_points = []
    for x in points:
        if options['dependence']['name'] == 'wavelength' or options['dependence']['name'] == 'angle':
            dependence_points.append(options_defining_func(x))
        else:
            dependence_points.append(structure_defining_func(x))

    # Do not calculate here, if the dependence is not structural parameter...
    if options['dependence']['name'] == 'wavelength' or options['dependence']['name'] == 'angle':
        coherent_subsystems = get_coherent_subsets(structure)
        for item in coherent_subsystems:
            item['toeplitz_matrix'] = calculate_toeplitz_matrices(item, options)
        return get_coeffs(np.array(list(map(options_func, dependence_points))))
    else:
        return get_coeffs(np.array(list(map(structure_func, dependence_points))))


# %%

perm_capping = 1.7**2
perm_grating = 1.3 ** 2
fill_grating = 0.5

structure = {
    'bound_selectors': [
        # ['data', 0, 'permittivity', 4],

        # ['data', 0, 'permittivity', 8],
    ],
    'superstrate': 1,  # indexes, not perm
    'substrate': 1,
    'period': 840,
    'data': [
        {'width': 3000, 'periodic': False, 'coherent': False,
         'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]},
        {'width': 740, 'periodic': True, 'coherent': True, 'materials': [
            {'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating], 'start': 0,
             'stop': fill_grating},
            {'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping], 'start': fill_grating, 'stop': 1},
        ]},
        {'width': 50000, 'periodic': False, 'coherent': False,
         'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]},
    ]}

options = {
    'bound_selectors': [
        ['wavelength'],  # Select the proper attribute to be dependant
    ],
    'n_terms': 20,
    'n_select': 4,
    'wavelength': 833,
    'angle': 6,
    'divisions': 200,
    'divisions_start': 650,
    'divisions_end': 1200,
    'dependence': {
        'name': 'wavelength',  ##angle // wavelength // structure
        'label': 'Wv [nm]'
    }
}

values, angles = general_func(structure_defining_func(structure), options_defining_func(options))

# %%

options = options_defining_func(options)(None)

show_specular_ref_pol(values, options)


