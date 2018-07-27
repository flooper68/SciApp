import math
import sys
import time
import uuid

from funkce import *
from zobrazovaci_funkce import *


# Read data from stdin
def read_in():
    lines = sys.stdin.readlines()
    # Since our input would only be having one line, parse our JSON data from that
    return json.loads(lines[0])


def calc(params):
    startTime = time.time()

    structure = copy.deepcopy(params['structure'])

    options = copy.deepcopy(params['options'])

    structure['superstrate'] = float(structure['superstrate'])
    structure['substrate'] = float(structure['substrate'])
    structure['period'] = float(structure['period'])

    options['angle'] = float(options['angle'])
    options['wavelength'] = float(options['wavelength'])
    options['n_terms'] = int(options['n_terms'])
    options['n_select'] = int(options['n_select'])
    options['divisions_start'] = float(options['divisions_start'])
    options['divisions_end'] = float(options['divisions_end'])
    options['divisions'] = int(options['divisions'])

    for x in structure['data']:
        x['width'] = float(x['width'])
        if x.get('permittivity'):
            for index, item in enumerate(x['permittivity']):
                x['permittivity'][index] = complex(x['permittivity'][index]) ** 2
        else:
            for y in x['materials']:
                for index, item in enumerate(y['permittivity']):
                    y['permittivity'][index] = complex(y['permittivity'][index]) ** 2
                y['start'] = float(y['start'])
                y['stop'] = float(y['stop'])

    if options['dependence']['name'] == 'grating_thick':
        values, angles = general_func(grating_structure(structure), options_defining_func(options))
    elif options['dependence']['name'] == 'ref_thick':
        values, angles = general_func(ref_structure(structure), options_defining_func(options))
    else:
        values, angles = general_func(structure_defining_func(structure), options_defining_func(options))

    stopTime = time.time()
    toSend = Message()

    show_specular_diffraction_ref(values, options, show=False, toSend=toSend)
    show_specular_diffraction_trans(values, options, show=False, toSend=toSend)
    show_total(values, options, show=False, toSend=toSend)
    show_diff_angles(angles, options, show=False, toSend=toSend)
    show_first_negative_orders_trans_pol(values, options, show=False, toSend=toSend)
    show_first_negative_orders_ref_pol(values, options, show=False, toSend=toSend)
    show_first_negative_orders_ref_unpol(values, options, show=False, toSend=toSend)
    show_first_negative_orders_trans_unpol(values, options, show=False, toSend=toSend)
    show_comprehensive(values, options, show=False, toSend=toSend)

    # for x in range(len(points)):
    #     tempDict = {'point': points[x], 'totalRSpecular': values[x, 0], 'totalRDiffraction': values[x, 1]}
    #     toSend.data.append(tempDict)

    print('start sending', flush=True)
    toSend.data = {}
    toSend.options = params['options']
    toSend.duration = math.floor((stopTime - startTime))
    toSend.structure = params['structure']
    toSend.timestamp = time.time()
    toSend.name = 'Calc' + str(toSend.timestamp)
    toSend._id = str(uuid.uuid4())

    toSend.send()
    print('end sending', flush=True)


def main():
    func = sys.argv[1]
    params = read_in()

    calc(params)


def ref_structure(structure):
    def fun(param):
        temp = copy.deepcopy(structure)
        temp['data'][1]['width'] = param
        temp['data'][2]['width'] = temp['data'][2]['width'] - param
        temp['data'][3]['width'] = param
        return temp

    return fun


def grating_structure(structure):
    def fun(param):
        temp = copy.deepcopy(structure)
        temp['data'][2]['width'] = param - temp['data'][1]['width']
        return temp

    return fun

# start process
if __name__ == '__main__':
    main()
