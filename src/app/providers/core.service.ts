import {Injectable} from '@angular/core';
import {ElectronService} from './electron.service';
import {Observable} from 'rxjs/internal/Observable';
import {AppState} from '../store/app.reducer';
import {Store} from '@ngrx/store';
import {UpdateGratingsProgress} from '../store/gratings.action';
import {ChildProcess} from 'child_process';
import {Calculation} from '../models/calculation.model';

@Injectable()
export class CoreService {

    isPythonDist = true;
    pyProc: ChildProcess;

    constructor(private electron: ElectronService, public store: Store<AppState>) {

    }


    calculateStructure(structure: any): Observable<Calculation> {

        const settings = setUpStructure(structure);
        console.log(settings)
        return Observable.create((observer) => {
            let corePath: string;
            if (this.electron.os.platform() == 'win32') {
                corePath = __dirname + '/main/main.exe';
            }
            else {
                corePath = __dirname.replace('/view', '') + '/main/main';
            }


            if (!this.isPythonDist) {
                console.log('running script locally');
                this.pyProc = this.electron.childProcess.spawn('/Users/premyslciompa/miniconda3/envs/snakes/bin/python',
                    ['./python/main.py', name]);
            }
            else {
                console.log('running executable at', corePath);
                this.pyProc = this.electron.childProcess.execFile(corePath, [name]);
            }

            this.pyProc.stdin.write(JSON.stringify(settings));

            this.pyProc.stdin.end();

            let counter = 3;
            let sending = false;
            let finishedSending = false;
            let jsonData = '';
            let receivedData;
            let portion = counter / (settings.options.divisions + 2) * 100;
            this.store.dispatch(new UpdateGratingsProgress(portion));

            this.pyProc.stdout.on('data', (data: any) => {
                let msg: string;
                if (this.isPythonDist) {
                    msg = data
                }
                else {
                     msg = String.fromCharCode.apply(null, data);
                }
                if (msg.search('timer: tick') !== -1) {
                    portion = counter / (parseInt(settings.options.divisions) + 2) * 100;
                    this.store.dispatch(new UpdateGratingsProgress(portion));
                    counter++;
                } else {
                    if (msg.search('end sending') !== -1) {
                        jsonData = jsonData + '' + msg.replace('end sending', '');
                        sending = false;
                        finishedSending = true;
                    }
                    if (sending) {
                        jsonData = jsonData + '' + msg;
                    }
                    if (msg.search('start sending') !== -1) {
                        sending = true;
                    }
                }
                if (finishedSending) {
                    console.log(jsonData);
                    receivedData = JSON.parse(jsonData);
                    console.log(receivedData);
                    observer.next(receivedData);
                }
            });

            this.pyProc.stdout.on('end', () => {
                console.log('Msg ended');
                this.store.dispatch(new UpdateGratingsProgress(0));
                observer.complete();
            });

            this.pyProc.stderr.on('data', function (data: any) {
                observer.error(data.toString());
            });

        });

    }

    stopCalculation() {
        this.pyProc.kill();
    }

}

function setUpStructure(structure: any) {
    const perm_ref = structure.reflectiveIndex;
    const perm_capping = structure.cappingIndex;
    const perm_grating = structure.gratingIndex;
    const fill_factor = structure.gratingFill;
    const n_terms = 12;
    const n_select = 4;
    const substrateWidth = 50000;

    if (parseFloat(structure.gratingThickness) > parseFloat(structure.reflectiveThickness)) {
        switch (structure.dependence) {

            case('wavelength'): {
                return {
                    structure: {
                        bound_selectors: [],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.gratingThickness) - parseFloat(structure.reflectiveThickness),
                                'periodic': true,
                                'coherent': true,
                                'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            ['wavelength'],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'wavelength',
                            'label': 'Wavelength [nm]'
                        }
                    }
                };
            }

            case ('angle'): {
                return {
                    structure: {
                        bound_selectors: [],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.gratingThickness) - parseFloat(structure.reflectiveThickness),
                                'periodic': true,
                                'coherent': true,
                                'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            ['angle'],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'angle',
                            'label': 'Angle [deg]'
                        }
                    }
                };
            }

            case ('fill_factor'): {
                return {
                    structure: {
                        bound_selectors: [
                            ['data', 1, 'materials', 0, 'stop'],
                            ['data', 1, 'materials', 1, 'start'],
                            ['data', 2, 'materials', 0, 'stop'],
                            ['data', 2, 'materials', 1, 'start'],
                            ['data', 3, 'materials', 0, 'stop'],
                            ['data', 3, 'materials', 1, 'start'],
                        ],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.gratingThickness) - parseFloat(structure.reflectiveThickness),
                                'periodic': true,
                                'coherent': true,
                                'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            [],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'fill_factor',
                            'label': 'Fill factor'
                        }
                    }
                };
            }

            case ('period'): {
                return {
                    structure: {
                        bound_selectors: [['period']],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.gratingThickness) - parseFloat(structure.reflectiveThickness),
                                'periodic': true,
                                'coherent': true,
                                'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            [],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'period',
                            'label': 'Period [nm]'
                        }
                    }
                };
            }

            case ('capping_thick'): {
                return {
                    structure: {
                        bound_selectors: [
                            ['data', 0, 'width']
                        ],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.gratingThickness) - parseFloat(structure.reflectiveThickness),
                                'periodic': true,
                                'coherent': true,
                                'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            [],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'capping_thick',
                            'label': 'Capping layer thickness [nm]'
                        }
                    }
                };
            }
            case ('ref_thick'): {
                return {
                    structure: {
                        bound_selectors: [],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.gratingThickness), 'periodic': true,
                                'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            [],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'ref_thick',
                            'label': 'Reflective layer thickness [nm]'
                        }
                    }
                };
            }

            case ('grating_thick'): {
                return {
                    structure: {
                        bound_selectors: [],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.gratingThickness), 'periodic': true,
                                'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            [],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'grating_thick',
                            'label': 'Grating thickness [nm]'
                        }
                    }
                };
            }
        }
    }
    else {
        switch (structure.dependence) {


        case('wavelength'): {
                return {
                    structure: {
                        bound_selectors: [],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.gratingThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.reflectiveThickness) - parseFloat(structure.gratingThickness),
                                'periodic': true,
                                'coherent': true,
                                'materials': [
                                    {
                                        'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.gratingThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            ['wavelength'],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'wavelength',
                            'label': 'Wavelength [nm]'
                        }
                    }
                };
            }

            case ('angle'): {
                return {
                    structure: {
                        bound_selectors: [],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.gratingThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.reflectiveThickness) - parseFloat(structure.gratingThickness),
                                'periodic': true,
                                'coherent': true,
                                'materials': [
                                    {
                                        'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.gratingThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            ['angle'],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'angle',
                            'label': 'Angle [deg]'
                        }
                    }
                };
            }

            case ('fill_factor'): {
                return {
                    structure: {
                        bound_selectors: [
                            ['data', 1, 'materials', 0, 'stop'],
                            ['data', 1, 'materials', 1, 'start'],
                            ['data', 2, 'materials', 0, 'stop'],
                            ['data', 2, 'materials', 1, 'start'],
                            ['data', 3, 'materials', 0, 'stop'],
                            ['data', 3, 'materials', 1, 'start'],
                        ],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.gratingThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.reflectiveThickness) - parseFloat(structure.gratingThickness),
                                'periodic': true,
                                'coherent': true,
                                'materials': [
                                    {
                                        'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.gratingThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            [],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'fill_factor',
                            'label': 'Fill factor'
                        }
                    }
                };
            }

            case ('period'): {
                return {
                    structure: {
                        bound_selectors: [['period']],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.gratingThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.reflectiveThickness) - parseFloat(structure.gratingThickness),
                                'periodic': true,
                                'coherent': true,
                                'materials': [
                                    {
                                        'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.gratingThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            [],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'period',
                            'label': 'Period [nm]'
                        }
                    }
                };
            }

            case ('capping_thick'): {
                return {
                    structure: {
                        bound_selectors: [
                            ['data', 0, 'width']
                        ],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.gratingThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.reflectiveThickness) - parseFloat(structure.gratingThickness),
                                'periodic': true,
                                'coherent': true,
                                'materials': [
                                    {
                                        'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.gratingThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            [],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'capping_thick',
                            'label': 'Capping layer thickness [nm]'
                        }
                    }
                };
            }
            case ('ref_thick'): {
                return {
                    structure: {
                        bound_selectors: [],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.gratingThickness), 'periodic': true,
                                'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            [],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'ref_thick',
                            'label': 'Reflective layer thickness [nm]'
                        }
                    }
                };
            }

            case ('grating_thick'): {
                return {
                    structure: {
                        bound_selectors: [],
                        'superstrate': 1,
                        'substrate': 1,
                        'period': structure.gratingPeriod,
                        'data': [
                            {
                                'width': structure.cappingThickness, 'periodic': false, 'coherent': !structure.cappingIncoherent,
                                'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': 0, 'stop': fill_factor},
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': parseFloat(structure.gratingThickness), 'periodic': true,
                                'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {
                                        'permittivity': [perm_capping, 0, 0, 0, perm_capping, 0, 0, 0, perm_capping],
                                        'start': fill_factor,
                                        'stop': 1
                                    },
                                ]
                            },
                            {
                                'width': structure.reflectiveThickness, 'periodic': true, 'coherent': true, 'materials': [
                                    {
                                        'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating],
                                        'start': 0,
                                        'stop': fill_factor
                                    },
                                    {'permittivity': [perm_ref, 0, 0, 0, perm_ref, 0, 0, 0, perm_ref], 'start': fill_factor, 'stop': 1},
                                ]
                            },
                            {
                                'width': substrateWidth, 'periodic': false, 'coherent': false,
                                'permittivity': [perm_grating, 0, 0, 0, perm_grating, 0, 0, 0, perm_grating]
                            }
                        ]
                    },
                    options: {
                        'bound_selectors': [
                            [],
                        ],
                        'n_terms': n_terms,
                        'n_select': n_select,
                        'wavelength': structure.wavelength,
                        'angle': structure.angle,
                        'divisions': structure.divisions,
                        'divisions_start': structure.from,
                        'divisions_end': structure.to,
                        'dependence': {
                            'name': 'grating_thick',
                            'label': 'Grating thickness [nm]'
                        }
                    }
                };
            }
        }
    }


}
