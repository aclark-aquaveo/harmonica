from .adcirc_database import AdcircDB
from .leprovost_database import LeProvostDB
from .resource import ResourceManager
from .tpxo_database import TpxoDB


class Constituents:
    """Class for extracting tidal constituent data

    Attributes:
        _current_model (:obj:`tidal_database.TidalDB`): The tidal model currently being used for extraction

    """
    def __init__(self, model=ResourceManager.DEFAULT_RESOURCE):
        """Constructor tidal constituent extractor interface.

        Use this class as opposed to the lower level implementations

        Args:
            model (:obj:`str`, optional): Name of the tidal model. See resource.py for supported models.

        """
        self._current_model = None
        self.change_model(model)

    @property
    def data(self):
        """:obj:`Pandas.DataFrame` Access the underlying Pandas data frame of this constituent object."""
        if self._current_model:
            return self._current_model.data
        return []

    @data.setter
    def data(self, value):
        if self._current_model:
            self._current_model.data = value

    def change_model(self, new_model):
        new_model = new_model.lower()
        if self._current_model and self._current_model.model == new_model:
            return  # Already have the correct impl and resources for this model, nothing to do.

        if new_model in ResourceManager.TPXO_MODELS:  # Switch to a TPXO model
            self._current_model = TpxoDB(new_model)
        elif new_model in ResourceManager.LEPROVOST_MODELS:
            self._current_model = LeProvostDB(new_model)
        elif new_model in ResourceManager.ADCIRC_MODELS:
            self._current_model = AdcircDB()
        else:
            supported_models = (
                ", ".join(ResourceManager.TPXO_MODELS) + ", " +
                ", ".join(ResourceManager.LEPROVOST_MODELS) + ", " +
                ", ".join(ResourceManager.ADCIRC_MODELS)
            )
            raise ValueError("Model not supported: \'{}\'. Must be one of: {}.".format(
                new_model, supported_models.strip()
            ))

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
        if model and model.lower() != self._current_model.model:
            self.change_model(model.lower())
        return self._current_model.get_components(locs, cons, positive_ph)

    def get_nodal_factor(self, names, timestamp, timestamp_middle=None):
        """Get the nodal factor for specified constituents at a specified time.

        Args:
            names (:obj:`list` of :obj:`str`): Names of the constituents to get nodal factors for
            timestamp (:obj:`datetime.datetime`): Stat date and time to extract constituent arguments at
            timestamp_middle (:obj:`datetime.datetime`, optional): Date and time to consider as the middle of the
                series. By default, just uses the start day with half the hours.

        Returns:
            :obj:`pandas.DataFrame`: Constituent data frames. Each row contains frequency, earth tidal reduction factor,
                amplitude, nodal factor, and equilibrium argument for one of the specified constituents. Rows labeled by
                constituent name.

        """
        return self._current_model.get_nodal_factor(names, timestamp, timestamp_middle)
