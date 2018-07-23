import sys

from coherent_case import *
from incoherent_case import *


class Message:
    inputs = {}
    name = ''
    graphs = {}
    timestamp = ''
    data = []
    duration = ''

    def send(self):
        print(json.dumps({
            'options': self.options,
            'structure': self.structure,
            'graphs': self.graphs,
            'timestamp': self.timestamp,
            'data': self.data,
            'duration': self.duration,
            'name': self.name
        }))


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

    for x in structure['data']:
        x['width'] = float(x['width'])
        if x.get('permittivity'):
            for index, item in enumerate(x['permittivity']):
                x['permittivity'][index] = complex(x['permittivity'][index])
        else:
            for y in x['materials']:
                for index, item in enumerate(y['permittivity']):
                    y['permittivity'][index] = complex(y['permittivity'][index])
                y['start'] = float(y['start'])
                y['stop'] = float(y['stop'])

    values, angles = general_func(structure_defining_func(structure), options_defining_func(options))

    show_specular_diffraction_ref(values, options, show=False)

    stopTime = time.time()

    toSend = Message()

    # for x in range(len(points)):
    #     tempDict = {'point': points[x], 'totalRSpecular': values[x, 0], 'totalRDiffraction': values[x, 1]}
    #     toSend.data.append(tempDict)

    print('start sending', flush=True)
    toSend.data = {}
    toSend.graphs = []
    toSend.options = params['options']
    toSend.duration = math.floor((stopTime - startTime))
    toSend.structure = params['structure']
    toSend.timestamp = time.time()
    toSend.name = 'Calc' + str(toSend.timestamp)
    toSend.send()
    print('end sending', flush=True)


def main():
    func = sys.argv[1]
    params = read_in()

    calc(params)


# start process
if __name__ == '__main__':
    main()
