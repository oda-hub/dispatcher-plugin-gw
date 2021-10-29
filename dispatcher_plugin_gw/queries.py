from h5py._hl import base
import numpy as np
from cdci_data_analysis.analysis.queries import BaseQuery, ProductQuery
from cdci_data_analysis.analysis.parameters import Float, Integer, Time, ParameterRange, Name, Parameter
from cdci_data_analysis.analysis.products import QueryOutput, CatalogProduct
from cdci_data_analysis.analysis.catalog import BasicCatalog
from gwpy.spectrogram import Spectrogram
from gwpy.timeseries.timeseries import TimeSeries
import h5py
from bokeh.plotting import figure
from bokeh.models import ColorBar, LinearColorMapper, HoverTool, CustomJS, Slider 
from bokeh.embed import components
from bokeh.layouts import row, widgetbox, column
from gwpy.time import from_gps
import base64
from astropy.io.ascii import read as aread
 
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
        
        #TODO: almost copy of cdci_data_analysis.analysis.plot_tools.Image
        #      would be better to implement needed functionality (axes titles etc.) there
        
        fig = figure(tools=['pan,box_zoom,box_select,wheel_zoom,reset,save,crosshair'], 
                   y_axis_type="log", 
                   y_range=(self.sgram.yindex[0].value, self.sgram.yindex[-1].value), 
                   x_range=(0, self.sgram.dx.value * len(self.sgram.xindex)),
                   plot_height=350,
                   plot_width=650,
                   x_axis_label = 'Time [seconds] from %s (%.1f)' % (evt.strftime("%Y-%m-%d %T UTC"), 
                                                                    self.sgram.x0.value),
                   y_axis_label = 'Frequency [Hz]')

        min_s = self.sgram.min().value
        max_s = self.sgram.max().value
        
        color_mapper = LinearColorMapper(palette='Plasma256', 
                                      low=min_s,
                                      high=max_s)

        fig_im = fig.image(image=[self.sgram.T], 
                           x=0, 
                           y=self.sgram.yindex[0].value, 
                           dw=self.sgram.dx.value * len(self.sgram.xindex), 
                           dh=self.sgram.yindex[-1].value-self.sgram.yindex[0].value, 
                           color_mapper=color_mapper)
        
        hover = HoverTool(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                          renderers=[fig_im])
        
        fig.add_tools(hover)
        
        color_bar = ColorBar(color_mapper=color_mapper, 
                             label_standoff=5, 
                             location=(0,0),
                             width=10)
        
        JS_code_slider = """
                   var vmin = low_slider.value;
                   var vmax = high_slider.value;
                   fig_im.glyph.color_mapper.high = vmax;
                   fig_im.glyph.color_mapper.low = vmin;
               """

        callback = CustomJS(args=dict(fig_im=fig_im), code=JS_code_slider)

        self.graph_min_slider = Slider(title="Norm. En. Min", 
                                       start=min_s, 
                                       end=max_s, 
                                       step=1, 
                                       value=min_s, 
                                       callback=callback,
                                       width=150)
        self.graph_max_slider = Slider(title="Norm. En. Max", 
                                       start=min_s, 
                                       end=max_s, 
                                       step=1, 
                                       value=0.8 * max_s, 
                                       callback=callback,
                                       width=150)

        self.graph_min_slider.on_change('value', self.change_image_contrast)
        self.graph_max_slider.on_change('value', self.change_image_contrast)

        callback.args["low_slider"] = self.graph_min_slider
        callback.args["high_slider"] = self.graph_max_slider

        fig.add_layout(color_bar, 'right')

        layout = row(
            fig, widgetbox(self.graph_min_slider, 
                           self.graph_max_slider, 
                          ),
        )

        script, div = components(layout)

        return script, div

    def change_image_contrast(self, attr, old, new):
        self.fig_im.glyph.color_mapper.update(low=self.graph_min_slider.value, high=self.graph_max_slider.value)
        
        
class GWStrainQuery(ProductQuery):
    def __init__(self, name):
                whiten = Boolean(value=True, name='whiten')
                fmin = Integer(value=30, name='fmin')
                fmax = Integer(value=400, name='fmax')
                parameters_list = [whiten, fmin, fmax]
                super().__init__(name, parameters_list)

    def get_data_server_query(self, instrument, config, **kwargs):
        param_dict = dict(t1 = instrument.get_par_by_name('T1').value,
                          t2 = instrument.get_par_by_name('T2').value,
                          detector = instrument.get_par_by_name('detector').value,
                          whiten = instrument.get_par_by_name('whiten').value,
                          fmin = instrument.get_par_by_name('fmin').value,
                          fmax = instrument.get_par_by_name('fmax').value)
        return instrument.data_server_query_class(instrument=instrument,
                                                  config=config,
                                                  param_dict=param_dict,
                                                  task='/api/v1.0/get/strain')

    def build_product_list(self, instrument, res, out_dir, api=False):
        prod_list = []
        if out_dir is None:
            out_dir = './'
        _o_dict = res.json()
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
            raise NotImplementedError
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
            query_out.prod_dictionary['file_name'] = 'strain.h5'
            query_out.prod_dictionary['image'] = plot_dict
            query_out.prod_dictionary['download_file_name'] = 'gw_strain.h5'
            query_out.prod_dictionary['prod_process_message'] = ''

        return query_out


class StrainProduct:
    def __init__(self, ori_strain, filt_strain=None, out_dir=None):
        self.ori_strain = ori_strain
        self.filt_strain = filt_strain
        self.out_dir = out_dir

    def write(self):
        with h5py.File(self.out_dir + '/strain.h5', 'w') as fd:
            dset = fd.create_dataset('ori_strain', data=self.ori_strain.value)
            dset.attrs['t0'] = self.ori_strain.t0
            dset.attrs['dx'] = self.ori_strain.dt
            if self.filt_strain is not None:
                dset1 = fd.create_dataset('filt_strain', data=self.filt_strain.value)
                dset1.attrs['t0'] = self.filt_strain.t0
                dset1.attrs['dx'] = self.filt_strain.dt
    
    def get_strains_plot(self):
        
        evt = from_gps(self.ori_strain.t0.value)
        times = np.arange(0, 
                          self.ori_strain.dt.value * len(self.ori_strain) + 1,
                          self.ori_strain.dt.value )
       
        hover = HoverTool(tooltips=[("x", "$x"), ("y", "$y")])
                          
        fig = figure(tools=[hover, 'pan,box_zoom,box_select,wheel_zoom,reset,save,crosshair'], 
                   plot_height=250,
                   plot_width=700,
                   x_axis_label = 'Time [seconds] from %s (%.1f)' % (evt.strftime("%Y-%m-%d %T UTC"), 
                                                                    self.ori_strain.t0.value),
                   )
        ln = fig.line(times, self.ori_strain.value, line_color='blue')

        if self.filt_strain is not None:
            hover1 = HoverTool(tooltips=[("x", "$x"), ("y", "$y")])
                          
            fig1 = figure(tools=[hover, 'pan,box_zoom,box_select,wheel_zoom,reset,save,crosshair'], 
                    plot_height=250,
                    plot_width=700,
                    x_axis_label = 'Time [seconds] from %s (%.1f)' % (evt.strftime("%Y-%m-%d %T UTC"), 
                                                                        self.ori_strain.t0.value),
                    )
            ln1 = fig1.line(times, self.filt_strain.value, line_color='blue')

            layout = column(fig, fig1)
        else:
            layout = fig

        script, div = components(layout)

        return script, div
    
    
class GWImageQuery(ProductQuery):
    def __init__(self, name):
        super().__init__(name, [])

    def get_data_server_query(self, instrument, config, **kwargs):
        param_dict = dict(t1 = instrument.get_par_by_name('T1').value,
                          t2 = instrument.get_par_by_name('T2').value)
        return instrument.data_server_query_class(instrument=instrument,
                                                  config=config,
                                                  param_dict=param_dict,
                                                  task='/api/v1.0/get/events')

    def build_product_list(self, instrument, res, out_dir, api=False):
        if out_dir is None:
            out_dir = './'
        _o_dict = res.json()
        asciicat = _o_dict['output']['asciicat']
        jpgdata = _o_dict['output']['jpgdata']
        prod_list = [asciicat, jpgdata]
        return prod_list

    def process_product_method(self, instrument, prod_list, api=False):
        if api is True:
            raise NotImplementedError
        else:
            image  = prod_list.prod_list[1]
            catalog = prod_list.prod_list[0]
            # prod.write() TODO: writing fits file(s)
            with open('image.jpeg', 'wb') as ofd:
                ofd.write(base64.b64decode(image))
            
            with open('catalog.ecsv', 'w') as ofd:
                ofd.write(catalog)
            
            script = ''
            div = f'<img src="data:image/jpeg;base64,{image}">'
            
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
            query_out.prod_dictionary['file_name'] = ['catalog.ecsv', 'image.jpeg'] # TODO: fits instead of png
            query_out.prod_dictionary['image'] = plot_dict
            query_out.prod_dictionary['download_file_name'] = 'gw_image.tar.gz' # FIXME: download doesn't work properly
            query_out.prod_dictionary['prod_process_message'] = ''
            
            catalog_table = aread(catalog)
            catalog_table.add_column(np.arange(len(catalog_table)), name='index', index=0)
            
            column_lists=[catalog_table[name].tolist() for name in catalog_table.colnames]
            for ID,_col in enumerate(column_lists):
                column_lists[ID] = [x if str(x)!='nan' else None for x in _col]

            catalog_dict = dict(cat_column_list=column_lists,
                    cat_column_names=catalog_table.colnames,
                    cat_column_descr=catalog_table.dtype.descr,
                    )
            
            query_out.prod_dictionary['catalog'] = catalog_dict

        return query_out