import {Component, Input, OnChanges, OnInit, ViewChild} from '@angular/core';
import {ShownGraphs} from '../../../models/shownGraphs.model';

@Component({
    selector: 'app-plots',
    templateUrl: './plots.component.html',
    styleUrls: ['./plots.component.css']
})
export class PlotsComponent implements OnInit, OnChanges {

    @Input() graphs: ShownGraphs;
    @Input() plotSettings: any;

    @ViewChild('plots') plots;


    constructor() {
    }

    ngOnInit() {
        setTimeout(() => {

            this.showingPlots(this.graphs, this.plotSettings);
        }, 1);
    }

    ngOnChanges() {
        this.showingPlots(this.graphs, this.plotSettings);
    }

    showingPlots(data: ShownGraphs, plotSettings) {
        const plotContainer = document.getElementById('plots');
        if (plotContainer) {
            while (plotContainer.lastChild) {
                plotContainer.removeChild(plotContainer.lastChild);
            }
            if (data) {
                window.mpld3.draw_figure('plots', data.specular_diff_ref);

                if (plotSettings.showTrans) {
                    window.mpld3.draw_figure('plots', data.specular_diff_trans);
                }
                if (plotSettings.showPolarized) {
                    window.mpld3.draw_figure('plots', data.first_negative_ref_pol);
                    if (plotSettings.showTrans) {
                        window.mpld3.draw_figure('plots', data.first_negative_trans_pol);
                    }
                } else {
                    window.mpld3.draw_figure('plots', data.first_negative_ref_unpol);
                    if (plotSettings.showTrans) {
                        window.mpld3.draw_figure('plots', data.first_negative_trans_unpol);
                    }
                }
                if (plotSettings.showDiffAngles) {
                    window.mpld3.draw_figure('plots', data.diff_angles);
                }
                // window.mpld3.draw_figure('plots', data.total_trans_ref);
                // window.mpld3.draw_figure('plots', data.total_abs);
            }
        }

    }
}
