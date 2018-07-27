import {BrowserModule} from '@angular/platform-browser';
import {NgModule} from '@angular/core';

import {AppComponent} from './app.component';
import {MainComponent} from './components/main/main.component';
import {GratingsComponent} from './components/gratings/gratings.component';
import {ElectronService} from './providers/electron.service';
import {BrowserAnimationsModule} from '@angular/platform-browser/animations';
import {StoreModule} from '@ngrx/store';
import {appReducer} from './store/app.reducer';
import {EffectsModule} from '@ngrx/effects';
import {
    MatButtonModule,
    MatCardModule,
    MatCheckboxModule,
    MatFormFieldModule,
    MatInputModule,
    MatListModule,
    MatProgressBarModule,
    MatSelectModule,
    MatSidenavModule,
    MatSlideToggleModule,
    MatTabsModule,
    MatToolbarModule
} from '@angular/material';
import {NavDrawerComponent} from './components/nav-drawer/nav-drawer.component';
import {RouterModule} from '@angular/router';
import {routes} from './app.routing';
import {CoreService} from './providers/core.service';
import {GratingsEffects} from './store/gratings.effects';
import {ReactiveFormsModule} from '@angular/forms';
import {StoreDevtoolsModule} from '@ngrx/store-devtools';
import {PlotsComponent} from './components/gratings/plots/plots.component';

@NgModule({
    declarations: [
        AppComponent,
        MainComponent,
        GratingsComponent,
        PlotsComponent,
        NavDrawerComponent
    ],
    imports: [
        BrowserModule,
        BrowserAnimationsModule,
        ReactiveFormsModule,

        RouterModule.forRoot(routes, {useHash: true}),
        StoreModule.forRoot(appReducer),
        EffectsModule.forRoot([GratingsEffects]),
        StoreDevtoolsModule.instrument(),

        MatSidenavModule,
        MatToolbarModule,
        MatListModule,
        MatFormFieldModule,
        MatInputModule,
        MatCardModule,
        MatButtonModule,
        MatSelectModule,
        MatTabsModule,
        MatProgressBarModule,
        MatCheckboxModule,
        MatSlideToggleModule
    ],
    providers: [ElectronService, CoreService],
    bootstrap: [AppComponent]
})
export class AppModule {
}
