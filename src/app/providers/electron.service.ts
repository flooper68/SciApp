import {Injectable} from '@angular/core';
import {ipcRenderer, remote, webFrame} from 'electron';
import * as childProcess from 'child_process';
import * as fs from 'fs';
import {Observable} from 'rxjs/internal/Observable';
import {Observer} from 'rxjs/internal/types';
import {CalculationVM} from '../models/calculation.model';
import * as Os from 'os';
import {Collection} from 'mongodb';


@Injectable()
export class ElectronService {

    ipcRenderer: typeof ipcRenderer;
    webFrame: typeof webFrame;
    remote: typeof remote;
    childProcess: typeof childProcess;
    fs: typeof fs;
    gratingsDb: Collection ;
    os: typeof Os;

    constructor() {
        // Conditional imports

        this.ipcRenderer = window.require('electron').ipcRenderer;
        this.webFrame = window.require('electron').webFrame;
        this.remote = window.require('electron').remote;

        this.childProcess = window.require('child_process');
        this.fs = window.require('fs');
        this.os = window.require('os');
        const Datastore = window.require('nedb');
        const dbPath = __dirname.replace('/view', '') + '/db/gratings.db';
        this.gratingsDb = new Datastore({filename: dbPath, autoload: true});

        const dbPathTingo = __dirname.replace('/view', '') + '/db';
        const Db = window.require('tingodb')().Db;
        const db = new Db(dbPathTingo, {});

        this.gratingsDb = db.collection('gratings');

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
        this.gratingsDb.remove({_id: data._id}, () => {
            this.gratingsDb.insert(data, function (err, newDoc) {   // Callback is optional
            // newDoc is the newly inserted document, including its _id
            // newDoc has no key called notToBeSaved since its value was undefined
            console.log('Calc saved', newDoc);
        });
        });
    }

    loadAllProjects() {
        return Observable.create((observer: Observer<any>) => {
            this.gratingsDb.find({}).toArray( function (err, docs) {
                const sorted = docs.sort((a,b) => {
                    return parseFloat(b.timestamp) - parseFloat(a.timestamp)
                });
                observer.next(sorted);
                observer.complete();
            });
        });

    }

    deleteProjectFromDb(id: string) {
        this.gratingsDb.remove({_id: id},() => {
            console.log('calc deleted', id)
        });
    }
}
