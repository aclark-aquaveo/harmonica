from abc import ABCMeta, abstractmethod
import os
import shutil
import urllib.request
from zipfile import ZipFile

import xarray as xr

from harmonica import config


class Resources(object):
    """Abstract base class for model resources

    """
    def __init__(self):
        """Base constructor

        """
        pass

    __metaclass__ = ABCMeta

    @abstractmethod
    def resource_attributes(self):
        """Get the resource attributes of a model (e.g. web url, compression type)

        Returns:
            dict: Dictionary of model resource attributes

        """
        return {}

    @abstractmethod
    def dataset_attributes(self):
        """Get the dataset attributes of a model (e.g. unit multiplier, grid dimensions)

        Returns:
            dict: Dictionary of model dataset attributes

        """
        return {}

    @abstractmethod
    def available_constituents(self):
        """Get all the available constituents of a model

        Returns:
            list: List of all the available constituents

        """
        return []

    @abstractmethod
    def constituent_groups(self):
        """Get all the available constituents of a model grouped by compatible file types

        Returns:
            list(list): 2-D list of available constituents, where the first dimension groups compatible files

        """
        return []

    @abstractmethod
    def constituent_resource(self, con):
        """Get the resource name of a constituent

        Returns:
            str: Name of the constituent's resource

        """
        return None


class Tpxo7Resources(Resources):
    """TPXO7 resources

    """
    TPXO7_CONS = {'K1', 'K2', 'M2', 'M4', 'MF', 'MM', 'MN4', 'MS4', 'N2', 'O1', 'P1', 'Q1', 'S2'}
    DEFAULT_RESOURCE_FILE = 'DATA/h_tpxo7.2.nc'

    def __init__(self):
        super().__init__()

    def resource_attributes(self):
        return {
            'url': 'ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo7.2_netcdf.tar.Z',
            'archive': 'gz',  # gzip compression
        }

    def dataset_attributes(self):
        return {
            'units_multiplier': 1.0,  # meter
        }

    def available_constituents(self):
        return self.TPXO7_CONS

    def constituent_groups(self):
        return [self.available_constituents()]

    def constituent_resource(self, con):
        if con.upper() in self.TPXO7_CONS:
            return self.DEFAULT_RESOURCE_FILE
        else:
            return None


class Tpxo8Resources(Resources):
    """TPXO8 resources

    """
    TPXO8_CONS = [
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
    ]

    def __init__(self):
        super().__init__()

    def resource_attributes(self):
        return {
            'url': "ftp://ftp.oce.orst.edu/dist/tides/TPXO8_atlas_30_v1_nc/",
            'archive': None,
        }

    def dataset_attributes(self):
        return {
            'units_multiplier': 0.001,  # mm to meter
        }

    def available_constituents(self):
        # get keys from const groups as list of lists and flatten
        return [c for sl in [grp.keys() for grp in self.TPXO8_CONS] for c in sl]

    def constituent_groups(self):
        return [self.TPXO8_CONS[0], self.TPXO8_CONS[1]]

    def constituent_resource(self, con):
        con = con.upper()
        for group in self.TPXO8_CONS:
            if con in group:
                return group[con]
        return None


class Tpxo9Resources(Resources):
    """TPXO9 resources

    """
    TPXO9_CONS = {'2N2', 'K1', 'K2', 'M2', 'M4', 'MF', 'MM', 'MN4', 'MS4', 'N2', 'O1', 'P1', 'Q1', 'S1', 'S2'}
    DEFAULT_RESOURCE_FILE = 'tpxo9_netcdf/h_tpxo9.v1.nc'

    def __init__(self):
        super().__init__()

    def resource_attributes(self):
        return {
            'url': "ftp://ftp.oce.orst.edu/dist/tides/Global/tpxo9_netcdf.tar.gz",
            'archive': 'gz',
        }

    def dataset_attributes(self):
        return {
            'units_multiplier': 1.0,  # meter
        }

    def available_constituents(self):
        return self.TPXO9_CONS

    def constituent_groups(self):
        return [self.available_constituents()]

    def constituent_resource(self, con):
        if con.upper() in self.TPXO9_CONS:
            return self.DEFAULT_RESOURCE_FILE
        else:
            return None


class LeProvostResources(Resources):
    """LeProvost resources

    """
    LEPROVOST_CONS = {'K1', 'K2', 'M2', 'N2', 'O1', 'P1', 'Q1', 'S2', 'NU2', 'MU2', '2N2', 'T2', 'L2'}
    DEFAULT_RESOURCE_FILE = 'leprovost_tidal_db.nc'

    def __init__(self):
        super().__init__()

    def resource_attributes(self):
        return {
            'url': 'http://sms.aquaveo.com/leprovost_tidal_db.zip',
            'archive': 'zip',  # zip compression
        }

    def dataset_attributes(self):
        return {
            'units_multiplier': 1.0,  # meter
            'num_lats': 361,
            'num_lons': 720,
            'min_lon': -180.0,
        }

    def available_constituents(self):
        return self.LEPROVOST_CONS

    def constituent_groups(self):
        return [self.available_constituents()]

    def constituent_resource(self, con):
        if con.upper() in self.LEPROVOST_CONS:
            return self.DEFAULT_RESOURCE_FILE
        else:
            return None


class FES2014Resources(Resources):
    """FES2014 resources

    """
    FES2014_CONS = {
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
    }

    def __init__(self):
        super().__init__()

    def resource_attributes(self):
        return {
            'url': None,  # Resources must already exist. Licensing restrictions prevent hosting files.
            'archive': None,
        }

    def dataset_attributes(self):
        return {
            'units_multiplier': 1.0,  # meter
            'num_lats': 2881,
            'num_lons': 5760,
            'min_lon': 0.0,
        }

    def available_constituents(self):
        return self.FES2014_CONS.keys()

    def constituent_groups(self):
        return [self.available_constituents()]

    def constituent_resource(self, con):
        con = con.upper()
        if con in self.FES2014_CONS:
            return self.FES2014_CONS[con]
        else:
            return None


class Adcirc2015Resources(Resources):
    """ADCIRC (v2015) resources

    """
    ADCIRC_CONS = {
        'M2', 'S2', 'N2', 'K1', 'M4', 'O1', 'M6', 'Q1', 'K2', 'L2', '2N2', 'R2', 'T2', 'LAMBDA2', 'MU2',
        'NU2', 'J1', 'M1', 'OO1', 'P1', '2Q1', 'RHO1', 'M8', 'S4', 'S6', 'M3', 'S1', 'MK3', '2MK3', 'MN4',
        'MS4', '2SM2', 'MF', 'MSF', 'MM', 'SA', 'SSA'
    }
    DEFAULT_RESOURCE_FILE = 'all_adcirc.nc'

    def __init__(self):
        super().__init__()

    def resource_attributes(self):
        return {
            'url': 'http://sms.aquaveo.com/',
            'archive': None,  # Uncompressed NetCDF file
        }

    def dataset_attributes(self):
        return {
            'units_multiplier': 1.0,  # meter
        }

    def available_constituents(self):
        return self.ADCIRC_CONS

    def constituent_groups(self):
        return [self.available_constituents()]

    def constituent_resource(self, con):
        if con.upper() in self.ADCIRC_CONS:
            return self.DEFAULT_RESOURCE_FILE
        else:
            return None


class ResourceManager(object):
    """Harmonica resource manager to retrieve and access tide models"""

    RESOURCES = {
        'tpxo7': Tpxo7Resources(),
        'tpxo8': Tpxo8Resources(),
        'tpxo9': Tpxo9Resources(),
        'leprovost': LeProvostResources(),
        'fes2014': FES2014Resources(),
        'adcirc2015': Adcirc2015Resources(),
    }
    TPXO_MODELS = {'tpxo7', 'tpxo8', 'tpxo9'}
    LEPROVOST_MODELS = {'fes2014', 'leprovost'}
    ADCIRC_MODELS = {'adcirc2015'}
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
        return self.model_atts.available_constituents()

    def get_units_multiplier(self):
        return self.model_atts.dataset_attributes()['units_multiplier']

    def download(self, resource, destination_dir):
        """Download a specified model resource."""
        if not os.path.isdir(destination_dir):
            os.makedirs(destination_dir)

        rsrc_atts = self.model_atts.resource_attributes()
        url = rsrc_atts['url']
        # Check if we can download resources for this model.
        if url is None:
            raise ValueError("Automatic fetching of resources is not available for the {} model.".format(self.model))

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
                        rsrcs = set(
                            self.model_atts.constituent_resource(con) for con in
                            self.model_atts.available_constituents()
                        )
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
        resources = set(
            self.model_atts.constituent_resource(con) for con in
            self.model_atts.available_constituents()
        )
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
        for const_group in self.model_atts.constituent_groups():
            rsrcs = set(self.model_atts.constituent_resource(const) for const in set(constituents) & set(const_group))

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
