import {Options} from './options.model';
import {Structure} from './structure.model';
import {ShownGraphs} from './shownGraphs.model';

export interface Calculation {
    _id: string;
    name: string;
    timestamp: number;
    graphs: ShownGraphs;
    data: Object;
    options: Options;
    structure: Structure;
    duration: Number;
}

export interface CalculationVM {
    _id: string;
    name: string;
    timestamp: number;
    graphs: ShownGraphs;
    data: Object;
    options: Options;
    structure: Structure;
    duration: Number;
    inDatabase: boolean;
}
