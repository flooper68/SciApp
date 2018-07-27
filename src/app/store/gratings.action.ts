import {Action} from '@ngrx/store';
import {CalculationVM} from '../models/calculation.model';


export const DO_NEW_CALCULATION = 'DO_NEW_CALCULATION';
export const UPDATE_SHOWN_DATA = 'UPDATE_SHOWN_DATA';
export const UPDATE_GRATINGS_FORM_DATA = 'UPDATE_GRATINGS_FORM_DATA';
export const UPDATE_GRATINGS_PROGRESS = 'UPDATE_GRATINGS_PROGRESS';
export const UPDATE_GRATINGS_PLOTS = 'UPDATE_GRATINGS_PLOTS';
export const STOP_GRATING_CALCULATION = 'STOP_GRATING_CALCULATION';
export const SAVE_GRATING_PROJECT = 'SAVE_GRATING_PROJECT';
export const SAVE_TO_DB_GRATING_PROJECT = 'SAVE_TO_DB_GRATING_PROJECT';
export const DELETE_PROJECT_FROM_DB = 'DELETE_PROJECT_FROM_DB';
export const LOAD_GRATING_PROJECT = 'LOAD_GRATING_PROJECT';
export const LOAD_GRATING_DATA = 'LOAD_GRATING_DATA';
export const GRATING_DATA_LOADED = 'GRATING_DATA_LOADED';
export const UPDATE_GRATING_DATA = 'UPDATE_GRATING_DATA';
export const ADD_GRATING_DATA = 'ADD_GRATING_DATA';
export const SET_CALCULATION_ACTIVE = 'SET_CALCULATION_ACTIVE';
export const UPDATE_CALCULATION_NAME = 'UPDATE_CALCULATION_NAME';
export const RESET_ACTIVE_CALCULATION = 'RESET_ACTIVE_CALCULATION';

export class UpdateGratingData implements Action {

    readonly type = UPDATE_GRATING_DATA;

    constructor(public payload) {

    }
}

export class LoadGratingData implements Action {
    readonly type = LOAD_GRATING_DATA;
}

export class SaveGratingProjectToDb implements Action {
    readonly type = SAVE_TO_DB_GRATING_PROJECT;

    constructor() {
    }
}

export class SaveGratingProjectAction implements Action {
    readonly type = SAVE_GRATING_PROJECT;

    constructor() {
    }
}

export class LoadGratingProjectAction implements Action {
    readonly type = LOAD_GRATING_PROJECT;
}

export class GratingDataLoadedAction implements Action {
    readonly type = GRATING_DATA_LOADED;

    constructor(public payload) {
    }
}

export class AddGratingProjectAction implements Action {
    readonly type = ADD_GRATING_DATA;

    constructor(public payload: CalculationVM) {
    }
}

export class SetCalculationActive implements Action {
    readonly type = SET_CALCULATION_ACTIVE;

    constructor(public payload: string) {
    }
}

export class DoNewCalculationAction implements Action {

    readonly type = DO_NEW_CALCULATION;

    constructor() {

    }

}


export class UpdateGratingsFormData implements Action {

    readonly type = UPDATE_GRATINGS_FORM_DATA;

    constructor(public payload) {

    }
}

export class UpdateGratingsProgress implements Action {

    readonly type = UPDATE_GRATINGS_PROGRESS;

    constructor(public payload) {

    }
}


export class UpdateCalculationName implements Action {

    readonly type = UPDATE_CALCULATION_NAME;

    constructor(public payload: string) {

    }
}

export class UpdateGratingsPlots implements Action {

    readonly type = UPDATE_GRATINGS_PLOTS;

    constructor(public payload) {

    }
}

export class StopGratingCalculation implements Action {

    readonly type = STOP_GRATING_CALCULATION;

}

export class DeleteProjectFromDbAction implements Action {
    readonly type = DELETE_PROJECT_FROM_DB;
}

export class ResetActiveCalculation implements Action {
    readonly type = RESET_ACTIVE_CALCULATION;
}

export type GratingsAction = DoNewCalculationAction | UpdateGratingsFormData
    | UpdateGratingsProgress | StopGratingCalculation | UpdateGratingsPlots | UpdateGratingData
    | LoadGratingData | AddGratingProjectAction | SetCalculationActive | GratingDataLoadedAction
    | DeleteProjectFromDbAction | UpdateCalculationName | SaveGratingProjectToDb | ResetActiveCalculation;


