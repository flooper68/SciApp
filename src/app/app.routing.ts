import {Routes} from '@angular/router';
import {GratingsComponent} from './components/gratings/gratings.component';

export const routes: Routes = [
    {path: 'gratings', component: GratingsComponent},
    {path: '**', redirectTo: 'gratings'}
];

