export interface Options {
    bound_selectors: string[];
    n_terms: number;
    n_select: number;
    wavelength: string;
    angle: string;
    divisions: string;
    divisions_start: string;
    divisions_end: string;
    dependence: {
        name: string;
        label: string
    };
    substrateWidth: string;
}
