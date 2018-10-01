import {Injectable} from '@angular/core';
import {Actions, Effect} from '@ngrx/effects';
import {catchError, filter, map, switchMap, tap, withLatestFrom} from 'rxjs/operators';
import {CoreService} from '../providers/core.service';
import {Store} from '@ngrx/store';
import {AppState} from './app.reducer';
import {ElectronService} from '../providers/electron.service';
import {
    ADD_GRATING_DATA,
    DELETE_PROJECT_FROM_DB,
    DO_NEW_CALCULATION,
    GRATING_DATA_LOADED,
    LOAD_GRATING_DATA,
    LOAD_GRATING_PROJECT,
    RESET_ACTIVE_CALCULATION,
    SAVE_GRATING_PROJECT,
    SAVE_TO_DB_GRATING_PROJECT,
    SET_CALCULATION_ACTIVE,
    STOP_GRATING_CALCULATION
} from './gratings.action';
import {Calculation} from '../models/calculation.model';

@Injectable()
export class GratingsEffects {

    @Effect()
    doNewCalculation = this.actions$.ofType(DO_NEW_CALCULATION).pipe(
        withLatestFrom(this.store$.select('gratings')),
        map(([action, state]) => {
            return state.viewModel.settingsForm;
        }),
        switchMap((formState) => {
            return this.core.calculateStructure(formState).pipe(
                catchError((error) => {
                    console.log('Error caught: ', error);
                    return [];
                }),
            );
        }),
        switchMap(data => {
            return [{type: ADD_GRATING_DATA, payload: {...data, inDatabase: false}}, {type: SET_CALCULATION_ACTIVE, payload: data._id}];
        })
    );

    @Effect({dispatch: false})
    stopCalculation = this.actions$.ofType(STOP_GRATING_CALCULATION).pipe(
        tap(() => {
            this.core.stopCalculation();
        })
    );

    @Effect({dispatch: false})
    saveProject = this.actions$.ofType(SAVE_GRATING_PROJECT).pipe(
        withLatestFrom(this.store$.select('gratings')),
        filter(([action, state]) => !!state.viewModel.activeCalculationId),
        map(([action, state]) => {
            const obj = {...state.calculations.find(val => val._id === state.viewModel.activeCalculationId)};
            obj.name = state.viewModel.projectForm.name;
            return obj;
        }),
        tap((data) => {
            this.electron.saveCurrentProject(data);
        })
    );

    @Effect()
    loadProject = this.actions$.ofType(LOAD_GRATING_PROJECT).pipe(
        switchMap(() => {
            return this.electron.loadProject();
        }),
        withLatestFrom(this.store$.select('gratings')),
        filter(([calculation, state]) => {
            return -1 == state.calculations.findIndex(item => item._id === calculation._id);
        }),
        map(([calculation, state]) => calculation),
        switchMap((data: Calculation) => {
            return [{type: ADD_GRATING_DATA, payload: {...data, inDatabase: false}}, {type: SET_CALCULATION_ACTIVE, payload: data._id}];
        }),
    );

    @Effect({dispatch: false})
    saveToDb = this.actions$.ofType(SAVE_TO_DB_GRATING_PROJECT).pipe(
        withLatestFrom(this.store$.select('gratings')),
        map(([action, state]) => state),
        filter(state => !!state.viewModel.activeCalculationId),
        map((state) => {
            const obj = {...state.calculations.find(val => val._id === state.viewModel.activeCalculationId)};
            obj.name = state.viewModel.projectForm.name;
            obj.inDatabase = true;
            return obj;
        }),
        tap((data) => {
            this.electron.saveProjectToDb(data);
        })
    );

    @Effect()
    loadFromDb = this.actions$.ofType(LOAD_GRATING_DATA).pipe(
        switchMap(() => {
            return this.electron.loadAllProjects();
        }),
        switchMap((data) => [{type: GRATING_DATA_LOADED, payload: data}])
    );

    @Effect()
    deleteFromDb = this.actions$.ofType(DELETE_PROJECT_FROM_DB).pipe(
        withLatestFrom(this.store$.select('gratings')),
        map(([action, state]) => state.viewModel.activeCalculationId),
        tap(id => this.electron.deleteProjectFromDb(id)),
        switchMap(() => [{type: RESET_ACTIVE_CALCULATION}])
    );

    constructor(private actions$: Actions, private core: CoreService, private store$: Store<AppState>,
                private electron: ElectronService) {

    }

}
