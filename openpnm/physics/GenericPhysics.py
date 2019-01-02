from openpnm.core import Subdomain, ModelsMixin
from openpnm.utils import Workspace, logging
logger = logging.getLogger(__name__)
ws = Workspace()


class GenericPhysics(Subdomain, ModelsMixin):
    r"""
    This generic class is meant as a starter for custom Physics objects

    It produces a blank object with no pore-scale models attached.  Users can
    add models from the ``models`` module (or create their own).

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this Physics should be attached

    phase : OpenPNM Phase object
        The Phase object to which this Physics applies

    geometry : OpenPNM Geometry object
        The Geometry object that defines the pores/throats where this Physics
        should be applied.

    name : str, optional
        A unique string name to identify the Physics object, typically same as
        instance name but can be anything.  If left blank, and name will be
        generated that include the class name and a random string.

    """

    def __init__(self, project=None, network=None, phase=None,
                 geometry=None, settings={}, **kwargs):

        # Define some default settings
        self.settings.update({'prefix': 'phys'})
        # Overwrite with user supplied settings, if any
        self.settings.update(settings)

        # Deal with network or project arguments
        if network is not None:
            if project is not None:
                assert network is project.network
            else:
                project = network.project

        super().__init__(project=project, **kwargs)

        network = self.project.network
        if network:
            if phase is None:
                raise Exception('All Physics objects must be associated ' +
                                'with a Phase')
            else:
                phase['pore.'+self.name] = False
                phase['throat.'+self.name] = False
            if geometry is None:
                geoms = self.project.geometries().values()
                if len(geoms) == 0:
                    pass
                elif len(geoms) == 1:
                    logger.info('No Geometry given, assigning ' + self.name +
                                ' to all pores and throats')
                else:
                    raise Exception('Multiple Geometry objects are defined; ' +
                                    'must specify which one ' + self.name +
                                    ' should be associated with')
                Ps = network.Ps
                Ts = network.Ts
            else:
                Ps = network.pores(geometry.name)
                Ts = network.throats(geometry.name)
            self._add_locations(pores=Ps, throats=Ts)
