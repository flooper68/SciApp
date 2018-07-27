import {async, ComponentFixture, TestBed} from '@angular/core/testing';

import {GratingsComponent} from './gratings.component';

describe('GratingsComponent', () => {
    let component: GratingsComponent;
    let fixture: ComponentFixture<GratingsComponent>;

    beforeEach(async(() => {
        TestBed.configureTestingModule({
            declarations: [GratingsComponent]
        })
            .compileComponents();
    }));

    beforeEach(() => {
        fixture = TestBed.createComponent(GratingsComponent);
        component = fixture.componentInstance;
        fixture.detectChanges();
    });

    it('should create', () => {
        expect(component).toBeTruthy();
    });
});
