from cdci_data_analysis.analysis.instrument import Instrument

from . import conf_file
from .dataserver_dispatcher import GWDispatcher
from .queries import (GWInstrumentQuery, GWSkymapQuery, 
                      GWSpectrogramQuery, GWStrainQuery)
from cdci_data_analysis.analysis.queries import SourceQuery


def gw_factory():
    src_query = SourceQuery('src_query')
    instr_query = GWInstrumentQuery('instr_query')
    gw_spec_query = GWSpectrogramQuery('gw_spectrogram_query')
    gw_strain_query = GWStrainQuery('gw_strain_query')
    gw_skymap_query = GWSkymapQuery('gw_skymap_query')

    query_dictionary = {}
    query_dictionary['gw_spectrogram'] = 'gw_spectrogram_query'
    query_dictionary['gw_strain'] = 'gw_strain_query'
    query_dictionary['gw_skymap_image'] = 'gw_skymap_query'

    return Instrument('gw', 
                      src_query=src_query,
                      instrumet_query=instr_query,
                      data_serve_conf_file=conf_file,
                      asynch=True, 
                      product_queries_list=[gw_spec_query, gw_strain_query, gw_skymap_query],
                      data_server_query_class=GWDispatcher,
                      query_dictionary=query_dictionary)
