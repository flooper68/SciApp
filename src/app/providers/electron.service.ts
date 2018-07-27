import {Injectable} from '@angular/core';
import {ipcRenderer, remote, webFrame} from 'electron';
import * as childProcess from 'child_process';
import * as fs from 'fs';
import {Observable} from 'rxjs/internal/Observable';
import * as Nedb from 'nedb';
import {Observer} from 'rxjs/internal/types';
import {CalculationVM} from '../models/calculation.model';


@Injectable()
export class ElectronService {

    ipcRenderer: typeof ipcRenderer;
    webFrame: typeof webFrame;
    remote: typeof remote;
    childProcess: typeof childProcess;
    fs: typeof fs;
    gratingsDb: Nedb;

    constructor() {
        // Conditional imports

        this.ipcRenderer = window.require('electron').ipcRenderer;
        this.webFrame = window.require('electron').webFrame;
        this.remote = window.require('electron').remote;

        this.childProcess = window.require('child_process');
        this.fs = window.require('fs');
        const Datastore = window.require('nedb');
        this.gratingsDb = new Datastore({filename: 'dist/db/gratings.db', autoload: true});

    }

    saveCurrentProject(content) {
        this.remote.dialog.showSaveDialog({}, (filename, bookmark) => {
            if (filename) {
                if (filename.search('.json') == -1) {
                    filename = filename + '.json';
                }
                this.fs.writeFile(filename, JSON.stringify(content), function (err) {
                    if (err) {
                        throw err;
                    }
                });
            }
        });

    }

    loadProject(): Observable<CalculationVM> {
        return Observable.create((observer) => {
            this.remote.dialog.showOpenDialog({}, (filename, bookmark) => {
                this.fs.readFile(filename[0], function (err, data) {
                    if (err) {
                        return console.error(err);
                    }
                    let object: CalculationVM = JSON.parse(data.toString());
                    object.inDatabase = false;
                    console.log(object);
                    observer.next(JSON.parse(data.toString()));
                });
            });
        });
    }

    saveProjectToDb(data) {
        this.gratingsDb.remove({_id: data._id});
        this.gratingsDb.insert(data, function (err, newDoc) {   // Callback is optional
            // newDoc is the newly inserted document, including its _id
            // newDoc has no key called notToBeSaved since its value was undefined
        });
    }

    loadAllProjects() {
        return Observable.create((observer: Observer<any>) => {
            this.gratingsDb.find({}, function (err, docs) {
                observer.next(docs);
                observer.complete();
            });
        });

    }

    deleteProjectFromDb(id: string) {
        this.gratingsDb.remove({_id: id});
    }
}
