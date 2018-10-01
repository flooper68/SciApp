import {ChangeDetectorRef, Component, OnDestroy, OnInit, ViewChild} from '@angular/core';
import {ElectronService} from '../../providers/electron.service';
import {Store} from '@ngrx/store';
import {AppState} from '../../store/app.reducer';
import {
    DeleteProjectFromDbAction,
    DoNewCalculationAction,
    LoadGratingData,
    LoadGratingProjectAction,
    SaveGratingProjectAction,
    SaveGratingProjectToDb,
    SetCalculationActive,
    StopGratingCalculation,
    UpdateCalculationName,
    UpdateGratingsFormData,
    UpdateGratingsPlots,
} from '../../store/gratings.action';
import {FormControl, FormGroup} from '@angular/forms';
import {Subscription} from 'rxjs/internal/Subscription';
import {debounceTime, distinctUntilChanged, first, map, tap, withLatestFrom} from 'rxjs/operators';
import {Observable} from 'rxjs/internal/Observable';
import {ShownGraphs} from '../../models/shownGraphs.model';

@Component({
    selector: 'app-gratings',
    templateUrl: './gratings.component.html',
    styleUrls: ['./gratings.component.css']
})
export class GratingsComponent implements OnInit, OnDestroy {

    @ViewChild('matTab') matTab;
    @ViewChild('tabPlots') tabPlots;

    public calcOptions = [
        {value: 'wavelength', text: 'Wavelength'},
        {value: 'angle', text: 'Angle'},
        {value: 'period', text: 'Period'},
        {value: 'fill_factor', text: 'Fill factor'},
        {value: 'capping_thick', text: 'Capping layer Thickness'},
        {value: 'ref_thick', text: 'Reflective layer Thickness'},
        {value: 'grating_thick', text: 'Grating Thickness'},
    ];

    gratingsForm: FormGroup;
    plotsSettingsForm: FormGroup;
    calculationInfoForm: FormGroup;

    progress: Observable<number>;
    gratingData: Observable<any>;
    calcDuration: Observable<any>;
    calcInDb: Observable<boolean>;

    subscriptions: Subscription [] = [];

    graphs: ShownGraphs;
    plotSettings: any;
    plotsShown = false;

    constructor(public store: Store<AppState>, private change: ChangeDetectorRef, private electron: ElectronService) {
        this.gratingsForm = new FormGroup({
            cappingIndex: new FormControl(),
            cappingThickness: new FormControl(),
            cappingIncoherent: new FormControl(),
            reflectiveIndex: new FormControl(),
            reflectiveThickness: new FormControl(),
            gratingIndex: new FormControl(),
            gratingThickness: new FormControl(),
            gratingPeriod: new FormControl(),
            gratingFill: new FormControl(),
            wavelength: new FormControl(),
            angle: new FormControl(),
            dependence: new FormControl(),
            from: new FormControl(),
            to: new FormControl(),
            divisions: new FormControl(),
            substrateWidth: new FormControl(),
        });
        this.plotsSettingsForm = new FormGroup({
            showPolarized: new FormControl(),
            showTrans: new FormControl(),
            showDiffAngles: new FormControl()
        });

        this.calculationInfoForm = new FormGroup({
            name: new FormControl()
        });

        this.subscriptions.push(this.calculationInfoForm.valueChanges.subscribe((data) => {
            this.store.dispatch(new UpdateCalculationName(data.name));
        }));

        this.subscriptions.push(this.store.select('gratings').pipe(
            map(state => state.viewModel.activeCalculationId),
            distinctUntilChanged(),
            withLatestFrom(this.store.select('gratings').pipe(map(state => state.viewModel.projectForm)))
        ).subscribe((state) => {
            if (this.calculationInfoForm.value !== state[1]) {
                this.calculationInfoForm.patchValue(state[1]);
            }
        }));

        this.subscriptions.push(this.plotsSettingsForm.valueChanges.subscribe((data) => {
            this.store.dispatch(new UpdateGratingsPlots(data));
        }));
        this.subscriptions.push(this.store.select('gratings').pipe(
            map(state => state.viewModel.plotSettings),
            first()
        ).subscribe((state) => {
            this.plotsSettingsForm.patchValue(state);
        }));

        this.subscriptions.push(this.gratingsForm.valueChanges.pipe(
            debounceTime(500)
        ).subscribe((data) => {
            this.store.dispatch(new UpdateGratingsFormData(data));
        }));
        this.subscriptions.push(this.store.select('gratings').pipe(
            map(state => state.viewModel.settingsForm)
        ).subscribe((state) => {
            if (this.gratingsForm.value !== state) {
                this.gratingsForm.patchValue(state);
            }
        }))
        ;
        // update progress, needs change detection push
        this.progress = this.store.select('gratings').pipe(
            map(state => state.viewModel.progress),
            tap(() => {
                setTimeout(() => {
                    this.change.detectChanges();
                }, 100);
            })
        );

        // List of calculations
        this.gratingData = this.store.select('gratings').pipe(
            map(state => state.calculations),
        );
        this.calcDuration = this.store.select('gratings').pipe(
            map(state => {
                const index = state.calculations.findIndex(item => item._id === state.viewModel.activeCalculationId);
                return state.calculations[index] ? state.calculations[index].duration : 0;
            })
        );
        this.calcInDb = this.store.select('gratings').pipe(
            map(state => {
                const index = state.calculations.findIndex(item => item._id === state.viewModel.activeCalculationId);
                return state.calculations[index] ? state.calculations[index].inDatabase : false;
            })
        );
        // Show graphs subscription on new shown graphs
        this.subscriptions.push(this.store.select('gratings').pipe(
            map(state => [state.viewModel.activeCalculationId, state.viewModel.plotSettings]),
            distinctUntilChanged((a, b) => {
                return a[1] == b[1] && a[0] == b[0];
            }),
            withLatestFrom(this.store.select('gratings').pipe(map(state => {
                const activeCalc = state.calculations.find(val => val._id === state.viewModel.activeCalculationId);
                return [activeCalc ? activeCalc.graphs : undefined];
            }))),
        ).subscribe(data => {
            this.graphs = data[1][0];
            this.plotSettings = data[0][1];
        }));

        this.store.dispatch(new LoadGratingData());
    }

    ngOnInit() {
    }

    ngOnDestroy() {
        this.subscriptions.forEach(item => item.unsubscribe());
    }

    calculate() {
        this.store.dispatch(new DoNewCalculationAction());
    }

    stopCalculation() {
        this.store.dispatch(new StopGratingCalculation());
    }


    saveProjectToFile() {
        this.store.dispatch(new SaveGratingProjectAction());
    }

    loadProject() {
        this.store.dispatch(new LoadGratingProjectAction());
    }

    saveProjectToDb() {
        this.store.dispatch(new SaveGratingProjectToDb());
    }

    deleteProjectFromDb() {
        this.store.dispatch(new DeleteProjectFromDbAction());
    }

    projectSelected(project) {
        this.store.dispatch(new SetCalculationActive(project._id));
    }
}
