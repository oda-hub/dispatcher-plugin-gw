from .dataserver_dispatcher import GWDispatcher
from cdci_data_analysis.analysis.instrument import Instrument
from .queries import GWSourceQuery, GWInstrumentQuery, GWSpectrogramQuery, GWStrainQuery, GWImageQuery
from . import conf_file

def gw_factory():
    src_query = GWSourceQuery('src_query')
    instr_query = GWInstrumentQuery('instr_query')
    gw_spec_query = GWSpectrogramQuery('gw_spectrogram_query')
    gw_strain_query = GWStrainQuery('gw_strain_query')
    gw_image_query = GWImageQuery('gw_image_query')

    query_dictionary = {}
    query_dictionary['gw_spectrogram'] = 'gw_spectrogram_query'
    query_dictionary['gw_strain'] = 'gw_strain_query'
    query_dictionary['gw_image'] = 'gw_image_query'

    return Instrument('gw', 
                      src_query=src_query,
                      instrumet_query=instr_query,
                      data_serve_conf_file=conf_file,
                      asynch=True, 
                      product_queries_list=[gw_spec_query, gw_strain_query, gw_image_query],
                      data_server_query_class=GWDispatcher,
                      query_dictionary=query_dictionary)