#! python3

from enum import Enum
from abc import ABCMeta, abstractmethod
import math

import pandas as pd

from .resource import ResourceManager

NCNST = 37

# Dictionary of NOAA constituent speed constants (deg/hr)
# Source: https://tidesandcurrents.noaa.gov
# The speed is the rate change in the phase of a constituent, and is equal to 360 degrees divided by the
# constituent period expressed in hours
NOAA_SPEEDS = {
    'OO1': 16.139101,
    '2Q1': 12.854286,
    '2MK3': 42.92714,
    '2N2': 27.895355,
    '2SM2': 31.015896,
    'K1': 15.041069,
    'K2': 30.082138,
    'J1': 15.5854435,
    'L2': 29.528479,
    'LAM2': 29.455626,
    'M1': 14.496694,
    'M2': 28.984104,
    'M3': 43.47616,
    'M4': 57.96821,
    'M6': 86.95232,
    'M8': 115.93642,
    'MF': 1.0980331,
    'MK3': 44.025173,
    'MM': 0.5443747,
    'MN4': 57.423832,
    'MS4': 58.984104,
    'MSF': 1.0158958,
    'MU2': 27.968208,
    'N2': 28.43973,
    'NU2': 28.512583,
    'O1': 13.943035,
    'P1': 14.958931,
    'Q1': 13.398661,
    'R2': 30.041067,
    'RHO': 13.471515,
    'S1': 15.0,
    'S2': 30.0,
    'S4': 60.0,
    'S6': 90.0,
    'SA': 0.0410686,
    'SSA': 0.0821373,
    'T2': 29.958933,
}


def convert_coords(coords):
    """Convert latitude coordinates to [-180, 180].

    Args:
        coords (:obj:`list` of :obj:`tuple` of :obj:`float`): latitude [-90, 90] and longitude [-180 180] or [0 360]
            of the requested point.

    Returns:
        :obj:`list` of :obj:`tuple` of :obj:`float`: The list of converted coordinates. None if a coordinate out of
            range was encountered

    """
    # Make sure point locations are valid lat/lon
    for idx, pt in enumerate(coords):
        y_lat = pt[0]
        x_lon = pt[1]
        if x_lon < 0.0:
            x_lon = x_lon + 360.0
        if x_lon > 180.0:
            x_lon = x_lon - 360.0
        if x_lon > 180.0 or x_lon < -180.0 or y_lat > 90.0 or y_lat < -90.0:
            # ERROR: Not in latitude/longitude
            return None
        else:
            coords[idx] = (y_lat, x_lon)
    return coords


class TidalDBEnum(Enum):
    TIDAL_DB_LEPROVOST = 0
    TIDAL_DB_ADCIRC = 1


class OrbitVariables(object):
    def __init__(self):
        self.dh = 0.0
        self.di = 0.0
        self.dn = 0.0
        self.dnu = 0.0
        self.dnup = 0.0
        self.dnup2 = 0.0
        self.dp = 0.0
        self.dp1 = 0.0
        self.dpc = 0.0
        self.ds = 0.0
        self.dxi = 0.0
        self.grterm = NCNST*[0.0]
        self.hour = 0.0
        self.nodfac = NCNST*[0.0]
        self.day = 0
        self.month = 0
        self.year = 0


class TidalDB(object):
    """The base class for extracting tidal data from a database.

    Attributes:
        orbit (:obj:`OrbitVariables`): The orbit variables.
        data (:obj:`list` of :obj:`pandas.DataFrame`): List of the constituent component DataFrames with one
            per point location requested from get_components(). Intended return value of get_components().
        resources (:obj:`harmonica.resource.ResourceManager`): Manages fetching of tidal data

    """
    """(:obj:`list` of :obj:`float`): The starting days of the months (non-leap year)."""
    day_t = [0.0, 31.0, 59.0, 90.0, 120.0, 151.0, 181.0, 212.0, 243.0, 273.0, 304.0, 334.0]
    """(float): PI divided by 180.0"""
    pi180 = math.pi/180.0

    def __init__(self, model):
        """Base class constructor for the tidal extractors

        Args:
            model (str): The name of the model. One of 'tpxo9', 'tpxo8', 'tpxo7', 'leprovost, 'adcircnwat', or
                'adcircnepac'

        """
        self.orbit = OrbitVariables()
        self.data = []
        self.resources = None
        self._model = None
        self.change_model(model)

    __metaclass__ = ABCMeta

    @property
    def model(self):
        """str: The name of the model. One of 'tpxo9', 'tpxo8', 'tpxo7', 'leprovost', 'adcircnwat', or 'adcircnepac'

        When setting the model to a different one than the current, required resources are downloaded.

        """
        return self._model

    @model.setter
    def model(self, value):
        self.change_model(value)

    def change_model(self, model):
        """Change the extractor model. If different than the current, required resources are downloaded.

        Args:
            model (str): The name of the model. One of: 'tpxo9', 'tpxo8', 'tpxo7', 'leprovost, 'adcircnwat', or
                'adcircnepac'

        """
        model = model.lower()
        if model == 'tpxo7_2':
            model = 'tpxo7'
        if model != self._model:
            self._model = model
            self.resources = ResourceManager(self._model)

    @abstractmethod
    def get_components(self, locs, cons, positive_ph):
        """Abstract method to get amplitude, phase, and speed of specified constituents at specified point locations.

        Args:
            locs (:obj:`list` of :obj:`tuple` of :obj:`float`): latitude [-90, 90] and longitude [-180 180] or [0 360]
                of the requested points.
            cons (:obj:`list` of :obj:`str`, optional): List of the constituent names to get amplitude and phase for. If
                not supplied, all valid constituents will be extracted.
            positive_ph (bool, optional): Indicate if the returned phase should be all positive [0 360] (True) or
                [-180 180] (False, the default).

        Returns:
           :obj:`list` of :obj:`pandas.DataFrame`: Implementations should return a list of dataframes of constituent
                information including amplitude (meters), phase (degrees) and speed (degrees/hour, UTC/GMT). The list is
                parallel with locs, where each element in the return list is the constituent data for the corresponding
                element in locs. Empty list on error. Note that function uses fluent interface pattern.

        """
        pass

    def have_constituent(self, name):
        """Determine if a constituent is valid for this tidal extractor.

        Args:
            name (str): The name of the constituent.

        Returns:
            bool: True if the constituent is valid, False otherwise

        """
        return name.upper() in self.resources.available_constituents()

    def get_nodal_factor(self, a_names, a_hour, a_day, a_month, a_year):
        """Get the nodal factor for specified constituents at a specified time.

        Args:
            a_names (:obj:`list` of :obj:`str`): Names of the constituents to get nodal factors for
            a_hour (float): The hour of the specified time. Can be fractional
            a_day (int): The day of the specified time.
            a_month (int): The month of the specified time.
            a_year (int): The year of the specified time.

        Returns:
            :obj:`pandas.DataFrame`: Constituent data frames. Each row contains frequency, earth tidal reduction factor,
                amplitude, nodal factor, and equilibrium argument for one of the specified constituents. Rows labeled by
                constituent name.

        """
        con_names = ['M2', 'S2', 'N2', 'K1', 'M4', 'O1', 'M6', 'MK3', 'S4', 'MN4', 'NU2', 'S6', 'MU2', '2N2', 'OO1',
                     'LAMBDA2', 'S1', 'M1', 'J1', 'MM', 'SSA', 'SA', 'MSF', 'MF', 'RHO1', 'Q1', 'T2', 'R2', '2Q1',
                     'P1', '2SM2', 'M3', 'L2', '2MK3', 'K2', 'M8', 'MS4']
        con_freqs = [0.000140518902509, 0.000145444104333, 0.000137879699487, 0.000072921158358, 0.000281037805017,
                     0.000067597744151, 0.000421556708011, 0.000213440061351, 0.000290888208666, 0.000278398601995,
                     0.000138232903707, 0.000436332312999, 0.000135593700684, 0.000135240496464, 0.000078244573050,
                     0.000142804901311, 0.000072722052166, 0.000070281955336, 0.000075560361380, 0.000002639203022,
                     0.000000398212868, 0.000000199106191, 0.000004925201824, 0.000005323414692, 0.000065311745349,
                     0.000064958541129, 0.000145245007353, 0.000145643201313, 0.000062319338107, 0.000072522945975,
                     0.000150369306157, 0.000210778353763, 0.000143158105531, 0.000208116646659, 0.000145842317201,
                     0.000562075610519, 0.000285963006842]
        con_etrf = [0.0690 for _ in range(NCNST)]
        con_etrf[3] = 0.736  # clK1
        con_etrf[5] = 0.695  # O1
        con_etrf[29] = 0.706  # P1
        con_etrf[25] = 0.695  # Q1
        con_etrf[0] = 0.693  # M2
        con_etrf[2] = 0.693  # N2
        con_etrf[1] = 0.693  # S2
        con_etrf[34] = 0.693  # K2
        con_amp = [0.0 for _ in range(NCNST)]
        con_amp[3] = 0.141565  # K1
        con_amp[5] = 0.100514  # O1
        con_amp[29] = 0.046834  # P1
        con_amp[25] = 0.019256  # Q1
        con_amp[0] = 0.242334  # M2
        con_amp[2] = 0.046398  # N2
        con_amp[1] = 0.112841  # S2
        con_amp[34] = 0.030704  # K2

        con_data = pd.DataFrame(columns=["amplitude", "frequency", "earth_tide_reduction_factor",
                                         "equilibrium_argument", "nodal_factor"])
        self.get_eq_args(a_hour, a_day, a_month, a_year)
        for idx, name in enumerate(a_names):
            name = name.upper()
            try:
                name_idx = con_names.index(name)
            except ValueError:
                continue
            equilibrium_arg = 0.0
            nodal_factor = 0.0
            if self.orbit.nodfac[name_idx] != 0.0:
                nodal_factor = self.orbit.nodfac[name_idx]
                equilibrium_arg = self.orbit.grterm[name_idx]
            con_data.loc[name] = [con_amp[name_idx], con_freqs[name_idx], con_etrf[name_idx], equilibrium_arg,
                                  nodal_factor]
        return con_data

    def get_eq_args(self, a_hour, a_day, a_month, a_year):
        """Get equilibrium arguments at a starting time.

        Args:
            a_hour (float): The starting hour.
            a_day (int): The starting day.
            a_month (int): The starting month.
            a_year (int): The starting year.

        """
        day_julian = self.get_day_julian(a_day, a_month, a_year)
        hrm = a_hour / 2.0
        self.nfacs(a_year, day_julian, hrm)
        self.gterms(a_year, day_julian, a_hour, hrm)

    @staticmethod
    def angle(a_number):
        """Converts the angle to be within 0-360.

        Args:
            a_number (float): The angle to convert.

        Returns:
            The angle converted to be within 0-360.

        """
        ret_val = a_number
        while ret_val < 0.0:
            ret_val += 360.0
        while ret_val > 360.0:
            ret_val -= 360.0
        return ret_val

    def get_day_julian(self, a_day, a_month, a_year):
        """Get a float representing the Julian date.

        Args:
            a_day (float): The day.
            a_month (float): The month.
            a_year (float): The year.

        Returns:
            The Julian date with epoch of January 1, 1900.

        """
        days = 12*[0.0]
        days[1] = 31.0
        year_offset = int(a_year-1900.0)
        yrlp = float(abs(year_offset) % 4)
        if year_offset < 0.0:
            yrlp *= -1.0
        d_inc = 0.0
        if yrlp == 0.0:
            d_inc = 1.0
        for i in range(2, 12):
            days[i] = self.day_t[i] + d_inc
        return days[int(a_month)-1]+a_day

    def set_orbit(self, a_year, a_day_julian, a_hour):
        """Determination of primary and secondary orbital functions.

        Args:
            a_year (float): The year.
            a_day_julian (float): The day in julian format.
            a_hour (float): The hour..

        """
        x = int((a_year-1901.0)/4.0)
        dyr = a_year-1900.0
        dday = a_day_julian+x-1.0

        # dn IS THE MOON'S NODE (CAPITAL N, TABLE 1, SCHUREMAN)
        self.orbit.dn = 259.1560564-19.328185764*dyr-0.0529539336*dday-0.0022064139*a_hour
        self.orbit.dn = self.angle(self.orbit.dn)
        n = self.orbit.dn*self.pi180

        # dp IS THE LUNAR PERIGEE (SMALL P, TABLE 1)
        self.orbit.dp = 334.3837214+40.66246584*dyr+0.111404016*dday+0.004641834*a_hour
        self.orbit.dp = self.angle(self.orbit.dp)
        # p = self.orbit.dp*self.pi180
        fi = math.acos(0.9136949-0.0356926*math.cos(n))
        self.orbit.di = self.angle(fi/self.pi180)

        nu = math.asin(0.0897056*math.sin(n)/math.sin(fi))
        self.orbit.dnu = nu/self.pi180
        xi = n-2.0*math.atan(0.64412*math.tan(n/2.0))-nu
        self.orbit.dxi = xi/self.pi180
        self.orbit.dpc = self.angle(self.orbit.dp-self.orbit.dxi)

        # dh IS THE MEAN LONGITUDE OF THE SUN (SMALL H, TABLE 1)
        self.orbit.dh = 280.1895014-0.238724988*dyr+0.9856473288*dday+0.0410686387*a_hour
        self.orbit.dh = self.angle(self.orbit.dh)

        # dp1 IS THE SOLAR PERIGEE (SMALL P1, TABLE 1)
        self.orbit.dp1 = 281.2208569+0.01717836*dyr+0.000047064*dday+0.000001961*a_hour
        self.orbit.dp1 = self.angle(self.orbit.dp1)

        # ds IS THE MEAN LONGITUDE OF THE MOON (SMALL S, TABLE 1)
        self.orbit.ds = 277.0256206+129.38482032*dyr+13.176396768*dday+0.549016532*a_hour
        self.orbit.ds = self.angle(self.orbit.ds)
        nup = math.atan(math.sin(nu)/(math.cos(nu)+0.334766/math.sin(2.0*fi)))
        self.orbit.dnup = nup/self.pi180
        nup2 = math.atan(math.sin(2.0*nu)/(math.cos(2.0*nu)+0.0726184/pow(math.sin(fi), 2.0)))/2.0
        self.orbit.dnup2 = nup2/self.pi180

    def nfacs(self, a_year, a_day_julian, a_hour):
        """Calculates node factors for constituent tidal signal.

        Args:
            a_year (float): The year.
            a_day_julian (float): The day in julian format.
            a_hour (float): The hour.

        Returns:
            The same values as found in table 14 of schureman.

        """
        self.set_orbit(a_year, a_day_julian, a_hour)
        # n = self.orbit.dn*self.pi180
        fi = self.orbit.di*self.pi180
        nu = self.orbit.dnu*self.pi180
        # xi = self.orbit.dxi*self.pi180
        # p = self.orbit.dp*self.pi180
        # pc = self.orbit.dpc*self.pi180
        sini = math.sin(fi)
        sini2 = math.sin(fi/2.0)
        sin2i = math.sin(2.0*fi)
        cosi2 = math.cos(fi/2.0)
        # tani2 = math.tan(fi/2.0)
        # EQUATION 197, SCHUREMAN
        # qainv = math.sqrt(2.310+1.435*math.cos(2.0*pc))
        # EQUATION 213, SCHUREMAN
        # rainv = math.sqrt(1.0-12.0*math.pow(tani2, 2)*math.cos(2.0*pc)+36.0*math.pow(tani2, 4))
        # VARIABLE NAMES REFER TO EQUATION NUMBERS IN SCHUREMAN
        eq73 = (2.0/3.0-math.pow(sini, 2))/0.5021
        eq74 = math.pow(sini, 2.0)/0.1578
        eq75 = sini*math.pow(cosi2, 2.0)/0.37988
        eq76 = math.sin(2*fi)/0.7214
        eq77 = sini*math.pow(sini2, 2.0)/0.0164
        eq78 = math.pow(cosi2, 4.0)/0.91544
        eq149 = math.pow(cosi2, 6.0)/0.8758
        # eq207 = eq75*qainv
        # eq215 = eq78*rainv
        eq227 = math.sqrt(0.8965*math.pow(sin2i, 2.0)+0.6001*sin2i*math.cos(nu)+0.1006)
        eq235 = 0.001+math.sqrt(19.0444*math.pow(sini, 4.0)+2.7702*pow(sini, 2.0)*math.cos(2.0*nu)+0.0981)
        # NODE FACTORS FOR 37 CONSTITUENTS:
        self.orbit.nodfac[0] = eq78
        self.orbit.nodfac[1] = 1.0
        self.orbit.nodfac[2] = eq78
        self.orbit.nodfac[3] = eq227
        self.orbit.nodfac[4] = math.pow(self.orbit.nodfac[0], 2.0)
        self.orbit.nodfac[5] = eq75
        self.orbit.nodfac[6] = math.pow(self.orbit.nodfac[0], 3.0)
        self.orbit.nodfac[7] = self.orbit.nodfac[0]*self.orbit.nodfac[3]
        self.orbit.nodfac[8] = 1.0
        self.orbit.nodfac[9] = math.pow(self.orbit.nodfac[0], 2.0)
        self.orbit.nodfac[10] = eq78
        self.orbit.nodfac[11] = 1.0
        self.orbit.nodfac[12] = eq78
        self.orbit.nodfac[13] = eq78
        self.orbit.nodfac[14] = eq77
        self.orbit.nodfac[15] = eq78
        self.orbit.nodfac[16] = 1.0
        # EQUATION 207 NOT PRODUCING CORRECT ANSWER FOR M1
        # SET NODE FACTOR FOR M1 = 0 UNTIL CAN FURTHER RESEARCH
        self.orbit.nodfac[17] = 0.0
        self.orbit.nodfac[18] = eq76
        self.orbit.nodfac[19] = eq73
        self.orbit.nodfac[20] = 1.0
        self.orbit.nodfac[21] = 1.0
        self.orbit.nodfac[22] = eq78
        self.orbit.nodfac[23] = eq74
        self.orbit.nodfac[24] = eq75
        self.orbit.nodfac[25] = eq75
        self.orbit.nodfac[26] = 1.0
        self.orbit.nodfac[27] = 1.0
        self.orbit.nodfac[28] = eq75
        self.orbit.nodfac[29] = 1.0
        self.orbit.nodfac[30] = eq78
        self.orbit.nodfac[31] = eq149
        # EQUATION 215 NOT PRODUCING CORRECT ANSWER FOR L2
        # SET NODE FACTOR FOR L2 = 0 UNTIL CAN FURTHER RESEARCH
        self.orbit.nodfac[32] = 0.0
        self.orbit.nodfac[33] = math.pow(self.orbit.nodfac[0], 2.0)*self.orbit.nodfac[3]
        self.orbit.nodfac[34] = eq235
        self.orbit.nodfac[35] = math.pow(self.orbit.nodfac[0], 4.0)
        self.orbit.nodfac[36] = eq78

    def gterms(self, a_year, a_day_julian, a_hour, a_hrm):
        """Determines the Greenwich equilibrium terms.

        Args:
            a_year (float): The year.
            a_day_julian (float): The day in julian format.
            a_hour (float): The hour for V0.
            a_hrm (float): The hour for U.

        Returns:
            The same values as found in table 15 of schureman.

        """
        # OBTAINING ORBITAL VALUES AT BEGINNING OF SERIES FOR V0
        self.set_orbit(a_year, a_day_julian, a_hour)
        s = self.orbit.ds
        p = self.orbit.dp
        h = self.orbit.dh
        p1 = self.orbit.dp1
        t = self.angle(180.0+a_hour*(360.0/24.0))

        # OBTAINING ORBITAL VALUES AT MIDDLE OF SERIES FOR U
        self.set_orbit(a_year, a_day_julian, a_hrm)
        nu = self.orbit.dnu
        xi = self.orbit.dxi
        nup = self.orbit.dnup
        nup2 = self.orbit.dnup2

        # SUMMING TERMS TO OBTAIN EQUILIBRIUM ARGUMENTS
        self.orbit.grterm[0] = 2.0*(t-s+h)+2.0*(xi-nu)
        self.orbit.grterm[1] = 2.0*t
        self.orbit.grterm[2] = 2.0*(t+h)-3.*s+p+2.0*(xi-nu)
        self.orbit.grterm[3] = t+h-90.0-nup
        self.orbit.grterm[4] = 4.0*(t-s+h)+4.0*(xi-nu)
        self.orbit.grterm[5] = t-2.0*s+h+90.0+2.0*xi-nu
        self.orbit.grterm[6] = 6.0*(t-s+h)+6.0*(xi-nu)
        self.orbit.grterm[7] = 3.0*(t+h)-2.0*s-90.0+2.0*(xi-nu)-nup
        self.orbit.grterm[8] = 4.0*t
        self.orbit.grterm[9] = 4.0*(t+h)-5.0*s+p+4.0*(xi-nu)
        self.orbit.grterm[10] = 2.0*t-3.0*s+4.0*h-p+2.0*(xi-nu)
        self.orbit.grterm[11] = 6.0*t
        self.orbit.grterm[12] = 2.0*(t+2.0*(h-s))+2.0*(xi-nu)
        self.orbit.grterm[13] = 2.0*(t-2.0*s+h+p)+2.0*(xi-nu)
        self.orbit.grterm[14] = t+2.0*s+h-90.0-2.0*xi-nu
        self.orbit.grterm[15] = 2.0*t-s+p+180.0+2.0*(xi-nu)
        self.orbit.grterm[16] = t
        fi = self.orbit.di*self.pi180
        pc = self.orbit.dpc*self.pi180
        top = (5.0*math.cos(fi)-1.0)*math.sin(pc)
        bottom = (7.0*math.cos(fi)+1.0)*math.cos(pc)
        q = math.degrees(math.atan2(top, bottom))
        self.orbit.grterm[17] = t-s+h-90.0+xi-nu+q
        self.orbit.grterm[18] = t+s+h-p-90.0-nu
        self.orbit.grterm[19] = s-p
        self.orbit.grterm[20] = 2.0*h
        self.orbit.grterm[21] = h
        self.orbit.grterm[22] = 2.0*(s-h)
        self.orbit.grterm[23] = 2.0*s-2.0*xi
        self.orbit.grterm[24] = t+3.0*(h-s)-p+90.0+2.0*xi-nu
        self.orbit.grterm[25] = t-3.0*s+h+p+90.0+2.0*xi-nu
        self.orbit.grterm[26] = 2.0*t-h+p1
        self.orbit.grterm[27] = 2.0*t+h-p1+180.0
        self.orbit.grterm[28] = t-4.0*s+h+2.0*p+90.0+2.0*xi-nu
        self.orbit.grterm[29] = t-h+90.0
        self.orbit.grterm[30] = 2.0*(t+s-h)+2.0*(nu-xi)
        self.orbit.grterm[31] = 3.0*(t-s+h)+3.0*(xi-nu)
        r = math.sin(2.0*pc)/((1.0/6.0)*math.pow((1.0/math.tan(0.5*fi)), 2)-math.cos(2.0*pc))
        r = math.atan(r)/self.pi180
        self.orbit.grterm[32] = 2.0*(t+h)-s-p+180.0+2.0*(xi-nu)-r
        self.orbit.grterm[33] = 3.0*(t+h)-4.0*s+90.0+4.0*(xi-nu)+nup
        self.orbit.grterm[34] = 2.0*(t+h)-2.0*nup2
        self.orbit.grterm[35] = 8.0*(t-s+h)+8.0*(xi-nu)
        self.orbit.grterm[36] = 2.0*(2.0*t-s+h)+2.0*(xi-nu)
        for ih in range(0, 37):
            self.orbit.grterm[ih] = self.angle(self.orbit.grterm[ih])
