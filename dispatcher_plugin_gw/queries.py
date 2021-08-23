from cdci_data_analysis.analysis.queries import BaseQuery, ProductQuery
from cdci_data_analysis.analysis.parameters import Float, Integer, Time, ParameterRange, Name, Parameter
from cdci_data_analysis.analysis.products import QueryOutput
from gwpy.spectrogram import Spectrogram
import h5py
from bokeh.plotting import show, figure, output_notebook
from bokeh.models import ColorBar, LinearColorMapper
from bokeh.embed import components
from gwpy.time import from_gps

class Boolean(Parameter):
    def __init__(self, value=None, name=None):

        super().__init__(value=value,
                         name=name
                         )
        self._set_val(value)

    @property
    def value(self):
        return self._value

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
    def __init__(self, name):
                whiten = Boolean(value=True, name='whiten')
                qmin = Integer(value=4, name='qmin')
                qmax = Integer(value=64, name='qmax')
                parameters_list = [whiten, qmin, qmax]
                super().__init__(name, parameters_list)

    def get_data_server_query(self, instrument, config, **kwargs):
        param_dict = dict(t1 = instrument.get_par_by_name('T1').value,
                          t2 = instrument.get_par_by_name('T2').value,
                          detector = instrument.get_par_by_name('detector').value,
                          whiten = instrument.get_par_by_name('whiten').value,
                          qmin = instrument.get_par_by_name('qmin').value,
                          qmax = instrument.get_par_by_name('qmax').value)
        return instrument.data_server_query_class(instrument=instrument,
                                                  config=config,
                                                  param_dict=param_dict,
                                                  task='/api/v1.0/get/spectrogram')

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
        
        # simple mpl plot, using spectrogram with logf=False in backend 
        # import io
        # import base64
        # plot = self.sgram.plot(dpi=60)
        # ax = plot.gca()
        # plot.colorbar(label="Normalised energy")
        # ax.grid(False)
        # ax.set_yscale('log')
        # bio = io.BytesIO()
        # plot.savefig(bio, format='jpg')
        # bio.seek(0)
        # jpgdata = base64.b64encode(bio.read())
        # script = ''
        # div = f'<img src="data:image/jpeg;base64,{jpgdata.decode()}">'
        
        evt = from_gps(self.sgram.x0.value)
        
        p = figure(#tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")], 
                y_axis_type="log", 
                y_range=(self.sgram.yindex[0].value, self.sgram.yindex[-1].value), 
                x_range=(0, self.sgram.dx.value * len(self.sgram.xindex)),
                plot_height=300,
                plot_width=600,
                x_axis_label = 'Time [seconds] from %s (%.1f)' % (evt.strftime("%Y-%m-%d %T UTC"), 
                                                                  self.sgram.x0.value),
                y_axis_label = 'Frequency [Hz]')

        #p.x_range.range_padding = 0

        color_map = LinearColorMapper(palette='Plasma256', 
                                      low=self.sgram.min().value,
                                      high=self.sgram.max().value)

        p.image(image=[self.sgram.T], 
                x=0, 
                y=self.sgram.yindex[0].value, 
                dw=self.sgram.dx.value * len(self.sgram.xindex), 
                dh=self.sgram.yindex[-1].value-self.sgram.yindex[0].value, 
                color_mapper=color_map)
        
        # FIXME: doesn't work with bokeh 0.12.16 in cdci_data_analysis requirements
        # does it really need to be such restricted?
        #color_bar = ColorBar(color_mapper=color_map, 
                            #  label_standoff=5, 
                            #  width=10, 
                            #  title="Normalised energy")
        
        #p.add_layout(color_bar, 'right')
        
        script, div = components(p)
        return script, div
