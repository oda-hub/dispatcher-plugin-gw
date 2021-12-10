import h5py
import numpy as np
from bokeh.embed import components
from bokeh.layouts import column, row, widgetbox
from bokeh.models import (ColorBar, CustomJS, HoverTool, LinearColorMapper,
                          Slider)
from bokeh.plotting import figure
from gwpy.time import from_gps

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
