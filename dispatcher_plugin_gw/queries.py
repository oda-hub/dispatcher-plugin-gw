from cdci_data_analysis.analysis.queries import BaseQuery, ProductQuery
from cdci_data_analysis.analysis.parameters import Time, ParameterRange, Name
from cdci_data_analysis.analysis.products import QueryOutput
from gwpy.spectrogram import Spectrogram
import h5py
from bokeh.plotting import figure
from bokeh.embed import components
import io
import base64

class GWSourceQuery(BaseQuery):
    def __init__(self, name):
        t1 = Time(value='2017-08-17T12:40:54', name='T1', Time_format_name='T_format')
        t2 = Time(value='2017-08-17T12:41:10', name='T2', Time_format_name='T_format')

        t_range = ParameterRange(t1, t2, 'time')

        token = Name(name_format='str', name='token',value=None)

        param_list = [t_range, token]

        super().__init__(name, param_list)

class GWInstrumentQuery(BaseQuery):
    def __init__(self, 
                 name):
        detector = Name(value='H1', name='detector')
        param_list = [detector]
        self.input_prod_list_name = None
        super().__init__(name, param_list)

class GWSpectrogramQuery(ProductQuery):
    def __init__(self,
                 name,
                 parameters_list = [],
                 **kwargs):
                 super().__init__(name, parameters_list, **kwargs)

    def get_data_server_query(self, instrument, config, **kwargs):
        param_dict = dict(t1 = instrument.get_par_by_name('T1').value,
                                                 t2 = instrument.get_par_by_name('T2').value,
                                                 detector = instrument.get_par_by_name('detector').value)
        return instrument.data_server_query_class(instrument=instrument,config=config,param_dict=param_dict,task='/api/v1.0/get/spectrogram')

    def build_product_list(self, instrument, res, out_dir, api=False):
        prod_list = []
        if out_dir is None:
            out_dir = './'
        _o_dict = res.json()
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
            raise NotImplementedError
        else:
            prod  = prod_list.prod_list[0]
            prod.write()
            script, div = prod.get_spectrogram_plot()
            html_dict = {}
            html_dict['script'] = script
            html_dict['div'] = div
            plot_dict = {}
            plot_dict['image'] = html_dict
            plot_dict['header_text'] = ''
            plot_dict['table_text'] = ''
            plot_dict['footer_text'] = ''

            query_out = QueryOutput()
            query_out.prod_dictionary['name'] = 'spectrogram'
            query_out.prod_dictionary['file_name'] = 'spectrogram.h5'
            query_out.prod_dictionary['image'] = plot_dict
            query_out.prod_dictionary['download_file_name'] = 'gw_spectrogram.h5'
            query_out.prod_dictionary['prod_process_message'] = ''

        return query_out


class SpectrogramProduct:
    def __init__(self, sgram, out_dir=None):
        self.sgram = sgram
        self.out_dir = out_dir

    def write(self):
        with h5py.File(self.out_dir + '/spectrogram.h5', 'w') as fd:
            dset = fd.create_dataset('sgram', data=self.sgram.value)
            dset.attrs['x0'] = self.sgram.x0
            dset.attrs['dx'] = self.sgram.dx
            dset.attrs['yindex'] = self.sgram.yindex

    def get_spectrogram_plot(self):
        plot = self.sgram.plot(dpi=60)
        ax = plot.gca()
        plot.colorbar(label="Normalised energy")
        ax.grid(False)
        ax.set_yscale('log')
        bio = io.BytesIO()
        plot.savefig(bio, format='jpg')
        bio.seek(0)
        jpgdata = base64.b64encode(bio.read())
        script = ''
        div = f'<img src="data:image/jpeg;base64,{jpgdata.decode()}">'
        return script, div
