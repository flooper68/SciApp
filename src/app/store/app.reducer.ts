import {ActionReducerMap} from '@ngrx/store';
import {gratingsReducer, GratingsState} from './gratings.reducer';


export interface AppState {
    gratings: GratingsState;
}

export const appReducer: ActionReducerMap<AppState> = {
    gratings: gratingsReducer,
};

