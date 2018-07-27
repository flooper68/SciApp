declare var nodeModule: NodeModule;

interface NodeModule {
    id: string;
}

declare var window: Window;

interface Window {
    process: any;
    require: any;
    mpld3: any;
    d3: any;
}
