import os
import shutil
import ssl
import urllib.request
from zipfile import ZipFile

import xarray as xr

from harmonica import config


class ResourceManager(object):
    """Harmonica resource manager to retrieve and access tide models"""

    # Dictionary of model information
    RESOURCES = {
        'tpxo9': {
            'resource_atts': {
                'url': "ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9_netcdf.tar.gz",
                'archive': 'gz',
            },
            'dataset_atts': {
                'units_multiplier': 1.,  # meters
            },
            'consts': [{  # grouped by dimensionally compatible files
                '2N2': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'K1': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'K2': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'M2': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'M4': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'MF': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'MM': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'MN4': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'MS4': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'N2': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'O1': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'P1': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'Q1': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'S1': 'tpxo9_netcdf/h_tpxo9.v1.nc',
                'S2': 'tpxo9_netcdf/h_tpxo9.v1.nc',
            }, ],
        },
        'tpxo8': {
            'resource_atts': {
                'url': "ftp://ftp.oce.orst.edu/dist/tides/TPXO8_atlas_30_v1_nc/",
                'archive': None,
            },
            'dataset_atts': {
                'units_multiplier': 0.001,  # mm to meter
            },
            'consts': [  # grouped by dimensionally compatible files
                {  # 1/30 degree
                    'K1': 'hf.k1_tpxo8_atlas_30c_v1.nc',
                    'K2': 'hf.k2_tpxo8_atlas_30c_v1.nc',
                    'M2': 'hf.m2_tpxo8_atlas_30c_v1.nc',
                    'M4': 'hf.m4_tpxo8_atlas_30c_v1.nc',
                    'N2': 'hf.n2_tpxo8_atlas_30c_v1.nc',
                    'O1': 'hf.o1_tpxo8_atlas_30c_v1.nc',
                    'P1': 'hf.p1_tpxo8_atlas_30c_v1.nc',
                    'Q1': 'hf.q1_tpxo8_atlas_30c_v1.nc',
                    'S2': 'hf.s2_tpxo8_atlas_30c_v1.nc',
                },
                {  # 1/6 degree
                    'MF': 'hf.mf_tpxo8_atlas_6.nc',
                    'MM': 'hf.mm_tpxo8_atlas_6.nc',
                    'MN4': 'hf.mn4_tpxo8_atlas_6.nc',
                    'MS4': 'hf.ms4_tpxo8_atlas_6.nc',
                },
            ],
        },
        'tpxo7': {
            'resource_atts': {
                'url': "ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo7.2_netcdf.tar.Z",
                'archive': 'gz',  # gzip compression
            },
            'dataset_atts': {
                'units_multiplier': 1.,  # meter
            },
            'consts': [{  # grouped by dimensionally compatible files
                'K1': 'DATA/h_tpxo7.2.nc',
                'K2': 'DATA/h_tpxo7.2.nc',
                'M2': 'DATA/h_tpxo7.2.nc',
                'M4': 'DATA/h_tpxo7.2.nc',
                'MF': 'DATA/h_tpxo7.2.nc',
                'MM': 'DATA/h_tpxo7.2.nc',
                'MN4': 'DATA/h_tpxo7.2.nc',
                'MS4': 'DATA/h_tpxo7.2.nc',
                'N2': 'DATA/h_tpxo7.2.nc',
                'O1': 'DATA/h_tpxo7.2.nc',
                'P1': 'DATA/h_tpxo7.2.nc',
                'Q1': 'DATA/h_tpxo7.2.nc',
                'S2': 'DATA/h_tpxo7.2.nc',
            }, ],
        },
        'leprovost': {
            'resource_atts': {
                'url': 'https://sms.aquaveo.com.s3.amazonaws.com/leprovost_tidal_db.zip',
                'archive': 'zip',  # zip compression
            },
            'dataset_atts': {
                'units_multiplier': 1.,  # meter
                'num_lats': 361,
                'num_lons': 720,
                'min_lon': -180.0
            },
            'consts': [{  # grouped by dimensionally compatible files
                'K1': 'leprovost_tidal_db.nc',
                'K2': 'leprovost_tidal_db.nc',
                'M2': 'leprovost_tidal_db.nc',
                'N2': 'leprovost_tidal_db.nc',
                'O1': 'leprovost_tidal_db.nc',
                'P1': 'leprovost_tidal_db.nc',
                'Q1': 'leprovost_tidal_db.nc',
                'S2': 'leprovost_tidal_db.nc',
                'NU2': 'leprovost_tidal_db.nc',
                'MU2': 'leprovost_tidal_db.nc',
                '2N2': 'leprovost_tidal_db.nc',
                'T2': 'leprovost_tidal_db.nc',
                'L2': 'leprovost_tidal_db.nc',
            }, ],
        },
        'fes2014': {  # Resources must already exist. Licensing restrictions prevent
            'resource_atts': {
                'url': None,
                'archive': None,
            },
            'dataset_atts': {
                'units_multiplier': 1.,  # meter
                'num_lats': 2881,
                'num_lons': 5760,
                'min_lon': 0.0
            },
            'consts': [{  # grouped by dimensionally compatible files
                '2N2': '2n2.nc',
                'EPS2': 'eps2.nc',
                'J1': 'j1.nc',
                'K1': 'k1.nc',
                'K2': 'k2.nc',
                'L2': 'l2.nc',
                'LA2': 'la2.nc',
                'M2': 'm2.nc',
                'M3': 'm3.nc',
                'M4': 'm4.nc',
                'M6': 'm6.nc',
                'M8': 'm8.nc',
                'MF': 'mf.nc',
                'MKS2': 'mks2.nc',
                'MM': 'mm.nc',
                'MN4': 'mn4.nc',
                'MS4': 'ms4.nc',
                'MSF': 'msf.nc',
                'MSQM': 'msqm.nc',
                'MTM': 'mtm.nc',
                'MU2': 'mu2.nc',
                'N2': 'n2.nc',
                'N4': 'n4.nc',
                'NU2': 'nu2.nc',
                'O1': 'o1.nc',
                'P1': 'p1.nc',
                'Q1': 'q1.nc',
                'R2': 'r2.nc',
                'S1': 's1.nc',
                'S2': 's2.nc',
                'S4': 's4.nc',
                'SA': 'sa.nc',
                'SSA': 'ssa.nc',
                'T2': 't2.nc',
            }, ],
        },
        'adcircnwat': {
            'resource_atts': {
                'url': 'http://sms.aquaveo.com/adcircnwattides.zip',
                'archive': 'zip',  # zip compression
            },
            'dataset_atts': {
                'units_multiplier': 1.,  # mete
            },
            'consts': [{  # grouped by dimensionally compatible files
                'M2': 'adcircnwattides.exe',
                'S2': 'adcircnwattides.exe',
                'N2': 'adcircnwattides.exe',
                'K1': 'adcircnwattides.exe',
                'M4': 'adcircnwattides.exe',
                'O1': 'adcircnwattides.exe',
                'M6': 'adcircnwattides.exe',
                'Q1': 'adcircnwattides.exe',
                'K2': 'adcircnwattides.exe',
            }, ],
        },
        'adcircnepac': {
            'resource_atts': {
                'url': 'http://sms.aquaveo.com/adcircnepactides.zip',
                'archive': 'zip',  # zip compression
            },
            'dataset_atts': {
                'units_multiplier': 1.,  # meter
            },
            'consts': [{  # grouped by dimensionally compatible files
                'M2': 'adcircnepactides.exe',
                'S2': 'adcircnepactides.exe',
                'N2': 'adcircnepactides.exe',
                'K1': 'adcircnepactides.exe',
                'M4': 'adcircnepactides.exe',
                'O1': 'adcircnepactides.exe',
                'M6': 'adcircnepactides.exe',
                'Q1': 'adcircnepactides.exe',
                'K2': 'adcircnepactides.exe',
            }, ],
        },
    }
    DEFAULT_RESOURCE = 'tpxo9'

    def __init__(self, model=DEFAULT_RESOURCE):
        if model not in self.RESOURCES:
            raise ValueError('Model not recognized.')
        self.model = model
        self.model_atts = self.RESOURCES[self.model]
        self.datasets = []

    def __del__(self):
        for d in self.datasets:
            d.close()

    def available_constituents(self):
        # get keys from const groups as list of lists and flatten
        return [c for sl in [grp.keys() for grp in self.model_atts['consts']] for c in sl]

    def get_units_multiplier(self):
        return self.model_atts['dataset_atts']['units_multiplier']

    def download(self, resource, destination_dir):
        """Download a specified model resource."""
        if not os.path.isdir(destination_dir):
            os.makedirs(destination_dir)

        rsrc_atts = self.model_atts['resource_atts']
        url = rsrc_atts['url']
        # Check if we can download resources for this model.
        if url is None:
            raise ValueError("Automatic fetching of resources is not available for the {} model.".format(self.model))

        if rsrc_atts['archive'] is None:
            url = "".join((url, resource))

        print('Downloading resource: {}'.format(url))

        path = os.path.join(destination_dir, resource)
        # TODO: We probably don't want to disable SSL certificate verification. Don't pass in a SSL context once
        # TODO:     the link on the Aquaveo website works without bypassing certificate verification.
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        with urllib.request.urlopen(url, context=ctx) as response:
            if rsrc_atts['archive'] is not None:
                if rsrc_atts['archive'] == 'gz':
                    import tarfile
                    try:
                        tar = tarfile.open(mode='r:{}'.format(rsrc_atts['archive']), fileobj=response)
                    except IOError as e:
                        print(str(e))
                    else:
                        rsrcs = set(c for sl in [x.values() for x in self.model_atts['consts']] for c in sl)
                        tar.extractall(path=destination_dir, members=[m for m in tar.getmembers() if m.name in rsrcs])
                        tar.close()
                elif rsrc_atts['archive'] == 'zip':  # Unzip .zip files
                    zip_file = os.path.join(destination_dir, os.path.basename(resource) + '.zip')
                    with open(zip_file, 'wb') as out_file:
                        shutil.copyfileobj(response, out_file)
                    # Unzip the files
                    print("Unzipping files to: {}".format(destination_dir))
                    with ZipFile(zip_file, 'r') as unzipper:
                        # Extract all the files in the archive
                        unzipper.extractall(path=destination_dir)
                    print("Deleting zip file: {}".format(zip_file))
                    os.remove(zip_file)  # delete the zip file
            else:
                with open(path, 'wb') as f:
                    f.write(response.read())

        return path

    def download_model(self, resource_dir=None):
        """Download all of the model's resources for later use."""
        resources = set(r for sl in [grp.values() for grp in self.model_atts['consts']] for r in sl)
        if not resource_dir:
            resource_dir = os.path.join(config['data_dir'], self.model)
        for r in resources:
            path = os.path.join(resource_dir, r)
            if not os.path.exists(path):
                self.download(r, resource_dir)
        return resource_dir

    def remove_model(self):
        """Remove all of the model's resources."""
        resource_dir = os.path.join(config['data_dir'], self.model)
        if os.path.exists(resource_dir):
            import shutil

            shutil.rmtree(resource_dir, ignore_errors=True)

    def get_datasets(self, constituents):
        """Returns a list of xarray datasets."""
        available = self.available_constituents()
        if any(const not in available for const in constituents):
            raise ValueError('Constituent not recognized.')
        # handle compatible files together
        self.datasets = []
        for const_group in self.model_atts['consts']:
            rsrcs = set(const_group[const] for const in set(constituents) & set(const_group))

            paths = set()
            if config['pre_existing_data_dir']:
                missing = set()
                for r in rsrcs:
                    path = os.path.join(config['pre_existing_data_dir'], self.model, r)
                    paths.add(path) if os.path.exists(path) else missing.add(r)
                rsrcs = missing
                if not rsrcs and paths:
                    self.datasets.append(xr.open_mfdataset(paths, engine='netcdf4', concat_dim='nc'))
                    continue

            resource_dir = os.path.join(config['data_dir'], self.model)
            for r in rsrcs:
                path = os.path.join(resource_dir, r)
                if not os.path.exists(path):
                    self.download(r, resource_dir)
                paths.add(path)

            if paths:
                self.datasets.append(xr.open_mfdataset(paths, engine='netcdf4', concat_dim='nc'))

        return self.datasets
