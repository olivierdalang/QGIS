
## Scnenarios

- file server (everybody access through network)
  - WAL disabled as per https://github.com/qgis/QGIS/pull/36741)
- personal share (one user accesses locally, the other one through network)
  - WAL will only be disabled for one user. Maybe we can assert there is no WAL file before trying to access DB in non WAL mode ?

## To try...

### Journal mode, locking mode...

Use following pragma :

`OGR_SQLITE_PRAGMA=journal_mode=MEMORY,locking_mode=EXCLUSIVE`

#### Notes
By default, `PRAGMA journal_mode;` returns `wal` (probably specific to DBManager).

By default, `PRAGMA locking_mode;` returns `normal` (even on network ! using DBManager).

#### See

- https://github.com/r-spatial/sf/issues/628

### Synchronous

Use synchronous mode `OGR_SQLITE_SYNCHRONOUS=ON`, which supposedly reduces risk of corrptuion at some performance price

#### Notes
By default, `PRAGMA synchronous;` returns `2` (=full). Try `3` (extra) ?

#### See

- https://trac.osgeo.org/gdal/wiki/ConfigOptions#OGR_SQLITE_SYNCHRONOUS

