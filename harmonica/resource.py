import os
import shutil
import urllib.request
from zipfile import ZipFile

import xarray as xr

from harmonica import config


class ResourceManager(object):
    """Harmonica resource manager to retrieve and access tide models"""

    # Dictionay of model information
    RESOURCES = {
        'tpxo9': {
            'resource_atts': {
                'url': "ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9_netcdf.tar.gz",
                'archive': 'gz',
            },
            'dataset_atts': {
                'units_multiplier': 1., # meters
            },
            'consts': [{  # grouped by dimensionally compatiable files
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
                'units_multiplier': 0.001, # mm to meter
            },
            'consts': [ # grouped by dimensionally compatiable files 
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
                'archive': 'gz', # gzip compression
            },
            'dataset_atts': {
                'units_multiplier': 1., # meter
            },
            'consts': [{  # grouped by dimensionally compatiable files
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
                'url': "http://sms.aquaveo.com/ADCIRC_Essentials.zip",
                'archive': 'zip',  # gzip compression
            },
            'dataset_atts': {
                'units_multiplier': 1.,  # meter
            },
            'consts': [{  # grouped by dimensionally compatible files
                'K1': 'LeProvost/K1.legi',
                'K2': 'LeProvost/K2.legi',
                'M2': 'LeProvost/M2.legi',
                'N2': 'LeProvost/N2.legi',
                'O1': 'LeProvost/O1.legi',
                'P1': 'LeProvost/P1.legi',
                'Q1': 'LeProvost/Q1.legi',
                'S2': 'LeProvost/S2.legi',
                'NU2': 'LeProvost/NU2.legi',
                'MU2': 'LeProvost/MU2.legi',
                '2N2': 'LeProvost/2N2.legi',
                'T2': 'LeProvost/T2.legi',
                'L2': 'LeProvost/L2.legi',
            }, ],
        },
        'adcircnwat': {
            'resource_atts': {
                'url': 'http://sms.aquaveo.com/adcircnwattides.zip',
                'archive': 'zip',  # gzip compression
            },
            'dataset_atts': {
                'units_multiplier': 1.,  # meter
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
                'archive': 'zip',  # gzip compression
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
        if rsrc_atts['archive'] is None:
            url = "".join((url, resource))

        print('Downloading resource: {}'.format(url))

        path = os.path.join(destination_dir, resource)
        with urllib.request.urlopen(url) as response:
            if rsrc_atts['archive'] is not None:
                if rsrc_atts['archive'] == 'gz':
                    import tarfile

                    try:
                        tar = tarfile.open(mode='r:{}'.format(rsrc_atts['archive']), fileobj=response)
                    except IOError as e:
                        print(str(e))
                    else:
                        rsrcs = set(c for sl in [x.values() for x in self.model_atts['consts']] for c in sl)
                        tar.extractall(path = destination_dir, members = [m for m in tar.getmembers() if m.name in rsrcs])
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
        # handle compatiable files together
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
