///<reference path="gratings.action.ts"/>
import {
    ADD_GRATING_DATA,
    DELETE_PROJECT_FROM_DB,
    GRATING_DATA_LOADED,
    GratingsAction,
    RESET_ACTIVE_CALCULATION,
    SAVE_TO_DB_GRATING_PROJECT,
    SET_CALCULATION_ACTIVE,
    UPDATE_CALCULATION_NAME,
    UPDATE_GRATINGS_FORM_DATA,
    UPDATE_GRATINGS_PLOTS,
    UPDATE_GRATINGS_PROGRESS
} from './gratings.action';
import {CalculationVM} from '../models/calculation.model';


export interface GratingsState {
    calculations: CalculationVM[];
    viewModel: {
        activeCalculationId: string,
        progress: number,
        plotSettings: {
            showPolarized: boolean,
            showTrans: boolean,
            showDiffAngles: boolean,
        },
        settingsForm: {
            cappingIndex: string,
            cappingThickness: string,
            cappingIncoherent: boolean,
            reflectiveIndex: string,
            reflectiveThickness: string,
            gratingIndex: string,
            gratingThickness: string,
            gratingPeriod: string,
            gratingFill: string,
            wavelength: string,
            angle: string,
            dependence: string,
            from: string,
            to: string,
            divisions: string,
        },
        projectForm: {
            name: string,
        }
    };
}

const initialState: GratingsState = {
    calculations: [],
    viewModel: {
        activeCalculationId: undefined,
        progress: 0,
        plotSettings: {
            showPolarized: true,
            showTrans: false,
            showDiffAngles: true
        },
        settingsForm: {
            cappingIndex: '1.5',
            cappingThickness: '0',
            cappingIncoherent: false,
            reflectiveIndex: '1-6j',
            reflectiveThickness: '20',
            gratingIndex: '1.5',
            gratingThickness: '100',
            gratingPeriod: '1000',
            gratingFill: '0.5',
            wavelength: '833',
            angle: '0',
            dependence: 'wavelength',
            from: '400',
            to: '800',
            divisions: '200',
        },
        projectForm: {
            name: '',
        },
    }
};


export function gratingsReducer(state = initialState, action: GratingsAction) {

    switch (action.type) {

        case ADD_GRATING_DATA:
            return {
                ...state,
                calculations: [...state.calculations, action.payload]
            };


        case SET_CALCULATION_ACTIVE: {
            const activeProject = state.calculations.find(item => item._id === action.payload);
            return {
                ...state,
                viewModel: {
                    ...state.viewModel,
                    activeCalculationId: action.payload,
                    settingsForm: {
                        cappingIndex: activeProject.structure.data[0].permittivity[0],
                        cappingThickness: activeProject.structure.data[0].width,
                        cappingIncoherent: !activeProject.structure.data[0].coherent,
                        reflectiveIndex: activeProject.structure.data[1].materials[0].permittivity[0],
                        reflectiveThickness: activeProject.structure.data[1].width,
                        gratingIndex: activeProject.structure.data[1].materials[1].permittivity[0],
                        gratingThickness: parseFloat(activeProject.structure.data[2].width) + parseFloat(activeProject.structure.data[1].width) + '',
                        gratingPeriod: activeProject.structure.period,
                        gratingFill: activeProject.structure.data[1].materials[0].stop,
                        wavelength: activeProject.options.wavelength,
                        angle: activeProject.options.angle,
                        dependence: activeProject.options.dependence.name,
                        from: activeProject.options.divisions_start,
                        to: activeProject.options.divisions_end,
                        divisions: activeProject.options.divisions,
                    },
                    projectForm: {
                        name: activeProject.name,
                    }
                }
            };
        }

        case UPDATE_GRATINGS_PROGRESS:
            return {
                ...state,
                viewModel: {
                    ...state.viewModel,
                    progress: action.payload
                }
            };


        case UPDATE_GRATINGS_FORM_DATA:
            return {
                ...state,
                viewModel: {
                    ...state.viewModel,
                    settingsForm: action.payload
                }
            };

        case UPDATE_GRATINGS_PLOTS:
            return {
                ...state,
                viewModel: {
                    ...state.viewModel,
                    plotSettings: action.payload
                }
            };


        case GRATING_DATA_LOADED:
            return {
                ...state,
                calculations: [...state.calculations, ...action.payload]
            };


        case DELETE_PROJECT_FROM_DB: {
            const arr = [...state.calculations];
            const index = arr.findIndex((item: CalculationVM) => {
                return item._id === state.viewModel.activeCalculationId;
            });
            arr.splice(index, 1);
            return {
                ...state,
                calculations: arr,

            };
        }

        case RESET_ACTIVE_CALCULATION:
            return {
                ...state,
                viewModel: {
                    ...state.viewModel,
                    activeCalculationId: undefined,
                    projectForm: {
                        name: ''
                    }
                }
            };

        case UPDATE_CALCULATION_NAME:
            return {
                ...state,
                viewModel: {
                    ...state.viewModel,
                    projectForm:
                        {name: action.payload}
                }
            };

        case SAVE_TO_DB_GRATING_PROJECT: {
            if (state.viewModel.activeCalculationId) {
                const index = state.calculations.findIndex(item => item._id === state.viewModel.activeCalculationId);
                const calcs = [...state.calculations];
                calcs[index].inDatabase = true;
                calcs[index].name = state.viewModel.projectForm.name;
                return {
                    ...state,
                    calculations: calcs
                };
            }
            return state;
        }

        default:
            return state;

    }
}

