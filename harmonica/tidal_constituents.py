from .adcirc_database import AdcircDB
from .leprovost_database import LeProvostDB
from .tpxo_database import TpxoDB


class Constituents:
    tpxo_models = ['tpxo7', 'tpxo8', 'tpxo9']
    leprovost_models = ['fes2014', 'leprovost']

    def __init__(self, model):
        self.current_model = None
        self.switch_model(model.lower())

    def switch_model(self, new_model):
        if self.current_model and self.current_model.model == new_model:
            return  # Already have the correct impl for this model, nothing to do.

        if new_model in self.tpxo_models:  # Switch to a TPXO model
            # If we already have a TPXO impl, change its version if necessary.
            if self.current_model and self.current_model.model in self.tpxo_models:
                self.current_model.change_model(new_model)
            else:  # Construct a new TPXO impl.
                self.current_model = TpxoDB(new_model)
        elif new_model in self.leprovost_models:
            # If we already have a LeProvost impl, change its version if necessary.
            if self.current_model and self.current_model.model in self.leprovost_models:
                self.current_model.change_model(new_model)
            else:  # Construct a new LeProvost impl.
                self.current_model = LeProvostDB(new_model)
        elif new_model == 'adcirc':
            self.current_model = AdcircDB()
        else:
            raise ValueError("Model not supported - {}".format(new_model))

    def get_components(self, locs, cons=None, positive_ph=False, model=None):
        """Abstract method to get amplitude, phase, and speed of specified constituents at specified point locations.

        Args:
            locs (:obj:`list` of :obj:`tuple` of :obj:`float`): latitude [-90, 90] and longitude [-180 180] or [0 360]
                of the requested points.
            cons (:obj:`list` of :obj:`str`, optional): List of the constituent names to get amplitude and phase for. If
                not supplied, all valid constituents will be extracted.
            positive_ph (bool, optional): Indicate if the returned phase should be all positive [0 360] (True) or
                [-180 180] (False, the default).
            model (:obj:`str`, optional): Name of the tidal model to use to query for the data. If not provided, current
                model will be used. If a model other than the current is provided, current model is switched.

        Returns:
           :obj:`list` of :obj:`pandas.DataFrame`: Implementations should return a list of dataframes of constituent
                information including amplitude (meters), phase (degrees) and speed (degrees/hour, UTC/GMT). The list is
                parallel with locs, where each element in the return list is the constituent data for the corresponding
                element in locs. Empty list on error. Note that function uses fluent interface pattern.

        """
        if model and model.lower() != self.current_model.model:
            self.switch_model(model.lower())
        return self.current_model.get_components(locs, cons, positive_ph)

    def get_nodal_factor(self, names, hour, day, month, year):
        """Get the nodal factor for specified constituents at a specified time.

        Args:
            names (:obj:`list` of :obj:`str`): Names of the constituents to get nodal factors for
            hour (float): The hour of the specified time. Can be fractional
            day (int): The day of the specified time.
            month (int): The month of the specified time.
            year (int): The year of the specified time.

        Returns:
            :obj:`pandas.DataFrame`: Constituent data frames. Each row contains frequency, earth tidal reduction factor,
                amplitude, nodal factor, and equilibrium argument for one of the specified constituents. Rows labeled by
                constituent name.

        """
        return self.current_model.get_nodal_factor(names, hour, day, month, year)