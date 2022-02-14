from cdci_data_analysis.analysis.parameters import (Angle, Integer, Name,
                                                    Parameter, ParameterRange,
                                                    Time)
from cdci_data_analysis.analysis.products import QueryOutput
from cdci_data_analysis.analysis.queries import BaseQuery, ProductQuery
from gwpy.spectrogram import Spectrogram
from gwpy.timeseries.timeseries import TimeSeries

from .products import SkymapProduct, SpectrogramProduct, StrainProduct
from astropy.time import Time as astropyTime

def check_time_interval(T1, T2, maxinterval=60):
    t1 = astropyTime(T1, format='isot')
    t2 = astropyTime(T2, format='isot')
    delta = t2 - t1
    if delta.sec > maxinterval:
        raise ValueError(f'Too long time interval. Current limit is {maxinterval}s')


class Boolean(Parameter):
    def __init__(self, value=None, name=None):

        super().__init__(value=value,
                         name=name
                         )
        self._set_val(value)

    @property
    def value(self):
        return str(self._value).lower() #because passed in json

    @value.setter
    def value(self, v):
        self._set_val(v)

    def _set_val(self, value):
        if value == 'false' or value == "False" or value == "no" \
                            or str(value) == "0" or value is False:
            self._value = False
        elif value == 'true' or value == "True" or value == "yes" \
                            or str(value) == "1" or value is True:
            self._value = True
        else:
            raise ValueError(f'Wrong value for parameter {self.name}')
# class GWSourceQuery(BaseQuery):
#     def __init__(self, name):
#         t1 = Time(value='2017-08-17T12:40:54', name='T1', Time_format_name='T_format')
#         t2 = Time(value='2017-08-17T12:41:10', name='T2', Time_format_name='T_format')

#         t_range = ParameterRange(t1, t2, 'time')

#         token = Name(name_format='str', name='token',value=None)

#         param_list = [t_range, token]

#         super().__init__(name, param_list)

class GWInstrumentQuery(BaseQuery):
    def __init__(self, 
                 name):
        detector = Name(value='H1', name='detector')
        param_list = [detector]
        self.input_prod_list_name = None
        super().__init__(name, param_list)

class GWSpectrogramQuery(ProductQuery):
    def __init__(self, name):
                whiten = Boolean(value=True, name='whiten')
                qmin = Integer(value=4, name='qmin')
                qmax = Integer(value=64, name='qmax')
                parameters_list = [whiten, qmin, qmax]
                super().__init__(name, parameters_list)

    def get_data_server_query(self, instrument, config, **kwargs):
        param_dict = dict(t1 = instrument.get_par_by_name('T1').get_value_in_default_format(),
                          t2 = instrument.get_par_by_name('T2').get_value_in_default_format(),
                          detector = instrument.get_par_by_name('detector').value,
                          whiten = instrument.get_par_by_name('whiten').value,
                          qmin = instrument.get_par_by_name('qmin').value,
                          qmax = instrument.get_par_by_name('qmax').value,
                         ) 
        check_time_interval(param_dict['t1'], param_dict['t2'])
        return instrument.data_server_query_class(instrument=instrument,
                                                  config=config,
                                                  param_dict=param_dict,
                                                  task='/api/v1.0/get/spectrogram')

    def build_product_list(self, instrument, res, out_dir, api=False):
        prod_list = []
        if out_dir is None:
            out_dir = './'
        if 'output' in res.json().keys(): # in synchronous mode
            _o_dict = res.json() 
        else:
            _o_dict = res.json()['data']
        sgram = Spectrogram(_o_dict['output']['value'],
                            unit = 's',
                            t0 = _o_dict['output']['x0'],
                            dt = _o_dict['output']['dx'],
                            frequencies = _o_dict['output']['yindex'])
        dp = SpectrogramProduct(sgram, out_dir=out_dir)
        prod_list.append(dp)       
        return prod_list

    def process_product_method(self, instrument, prod_list, api=False):
        if api is True:
            prod  = prod_list.prod_list[0]
            query_out = QueryOutput()
            query_out.prod_dictionary['gw_spectrogram_product'] = prod.serialize()
        else:
            prod  = prod_list.prod_list[0]
            prod.write()
                        
            plot_dict = {}
            plot_dict['image'] = prod.get_spectrogram_plot()
            plot_dict['header_text'] = ''
            plot_dict['table_text'] = ''
            plot_dict['footer_text'] = ''

            query_out = QueryOutput()
            query_out.prod_dictionary['name'] = 'spectrogram'
            query_out.prod_dictionary['file_name'] = 'spectrogram.h5'
            query_out.prod_dictionary['image'] = plot_dict
            query_out.prod_dictionary['download_file_name'] = 'gw_spectrogram.tar.gz'
            query_out.prod_dictionary['prod_process_message'] = ''
        return query_out


      
        
class GWStrainQuery(ProductQuery):
    def __init__(self, name):
                whiten = Boolean(value=True, name='whiten')
                fmin = Integer(value=30, name='fmin')
                fmax = Integer(value=400, name='fmax')
                parameters_list = [whiten, fmin, fmax]
                super().__init__(name, parameters_list)

    def get_data_server_query(self, instrument, config, **kwargs):
        param_dict = dict(t1 = instrument.get_par_by_name('T1').get_value_in_default_format(),
                          t2 = instrument.get_par_by_name('T2').get_value_in_default_format(),
                          detector = instrument.get_par_by_name('detector').value,
                          whiten = instrument.get_par_by_name('whiten').value,
                          fmin = instrument.get_par_by_name('fmin').value,
                          fmax = instrument.get_par_by_name('fmax').value)
        check_time_interval(param_dict['t1'], param_dict['t2'])
        return instrument.data_server_query_class(instrument=instrument,
                                                  config=config,
                                                  param_dict=param_dict,
                                                  task='/api/v1.0/get/strain')

    def build_product_list(self, instrument, res, out_dir, api=False):
        prod_list = []
        if out_dir is None:
            out_dir = './'
        if 'output' in res.json().keys(): # in synchronous mode
            _o_dict = res.json() 
        else:
            _o_dict = res.json()['data']
        ori_strain = TimeSeries(_o_dict['output']['ori_strain'],
                                t0 = _o_dict['output']['t0'],
                                dt = _o_dict['output']['dt'])
        filt_strain = TimeSeries(_o_dict['output']['bp_strain'],
                                t0 = _o_dict['output']['t0'],
                                dt = _o_dict['output']['dt'])
        dp = StrainProduct(ori_strain=ori_strain,
                            filt_strain=filt_strain,
                            out_dir=out_dir)
        prod_list.append(dp)       
        return prod_list

    def process_product_method(self, instrument, prod_list, api=False):
        if api is True:
            prod  = prod_list.prod_list[0]
            query_out = QueryOutput()
            query_out.prod_dictionary['gw_strain_product_list'] = prod.serialize()
        else:
            prod  = prod_list.prod_list[0]
            prod.write()
            script, div = prod.get_strains_plot()
            html_dict = {}
            html_dict['script'] = script
            html_dict['div'] = div
            plot_dict = {}
            plot_dict['image'] = html_dict
            plot_dict['header_text'] = ''
            plot_dict['table_text'] = ''
            plot_dict['footer_text'] = ''

            query_out = QueryOutput()
            query_out.prod_dictionary['name'] = 'strain'
            query_out.prod_dictionary['file_name'] = ['strain.h5', 'strain_bandpassed.h5']
            query_out.prod_dictionary['image'] = plot_dict
            query_out.prod_dictionary['download_file_name'] = 'gw_strain.tar.gz'
            query_out.prod_dictionary['prod_process_message'] = ''
        return query_out


    
    
class GWSkymapQuery(ProductQuery):
    def __init__(self, name):
        do_cone_search = Boolean(True, name='do_cone_search')
        radius = Angle(value = 0., units='deg', name='radius')
        level_threshold = Integer(10, name='level_threshold')
        contour_levels = Name('50,90', name='contour_levels')
        parameter_list = [do_cone_search, radius, level_threshold, contour_levels]
        super().__init__(name, parameter_list)

    def get_data_server_query(self, instrument, config, **kwargs):
        param_dict = dict(t1 = instrument.get_par_by_name('T1').get_value_in_default_format(),
                          t2 = instrument.get_par_by_name('T2').get_value_in_default_format(),
                          do_cone_search = instrument.get_par_by_name('do_cone_search').value,
                          ra = instrument.get_par_by_name('RA').value, 
                          dec = instrument.get_par_by_name('DEC').value, 
                          radius = instrument.get_par_by_name('radius').value, 
                          level_threshold = instrument.get_par_by_name('level_threshold').value, 
                          contour_levels = instrument.get_par_by_name('contour_levels').value, 
                          )
        return instrument.data_server_query_class(instrument=instrument,
                                                  config=config,
                                                  param_dict=param_dict,
                                                  task='/api/v1.0/get/conesearch')

    def build_product_list(self, instrument, res, out_dir, api=False):
        if out_dir is None:
            out_dir = './'
        if 'output' in res.json().keys(): # in synchronous mode
            _o_dict = res.json() 
        else:
            _o_dict = res.json()['data']
        asciicat = _o_dict['output']['asciicat']
        imagedata = _o_dict['output']['image']
        fits_data = _o_dict['output']['skymap_files']
        contours_data = _o_dict['output']['contours']
        
        prod_list = [SkymapProduct(asciicat, 
                                   imagedata, 
                                   fits_data, 
                                   contours_data, 
                                   out_dir)]
       
        return prod_list

    def process_product_method(self, instrument, prod_list, api=False):
        if api is True:
            skymap = prod_list.prod_list[0]
            query_out = QueryOutput()
            query_out.prod_dictionary['gw_skymap_product'] = skymap.serialize()
            query_out.prod_dictionary['catalog'] = skymap.get_catalog_dict(api=True) 
            
        else:
            skymap = prod_list.prod_list[0]
            
            script, div = skymap.get_plot()
            
            html_dict = {}
            html_dict['script'] = script
            html_dict['div'] = div
            plot_dict = {}
            plot_dict['image'] = html_dict
            plot_dict['header_text'] = ''
            plot_dict['table_text'] = ''
            plot_dict['footer_text'] = ''

            query_out = QueryOutput()
            query_out.prod_dictionary['name'] = 'image'
            query_out.prod_dictionary['file_name'] = skymap.write()
            query_out.prod_dictionary['image'] = plot_dict
            query_out.prod_dictionary['download_file_name'] = 'gw_skymap.tar.gz' 
            query_out.prod_dictionary['prod_process_message'] = ''

            query_out.prod_dictionary['catalog'] = skymap.get_catalog_dict(api=False)

        return query_out
