from functools import lru_cache
import json
import os

import h5py
import numpy as np
from astropy.io import fits
from astropy.io.ascii import read as aread
from bokeh.embed import components
from bokeh.layouts import column, row
from bokeh.models import (ColorBar, CustomJS, HoverTool, LinearColorMapper,
                          Slider)
from bokeh.plotting import figure
from gwpy.time import from_gps
from oda_api.data_products import NumpyDataUnit, NumpyDataProduct
from cdci_data_analysis.analysis.plot_tools import Image 

class SpectrogramProduct:
    def __init__(self, sgram, out_dir=None):
        self.sgram = sgram
        self.out_dir = out_dir

    def serialize(self):
        try:
            dy = self.sgram.dy.value
        except AttributeError:
            dy = 'none'
            
        return dict(
            value = self.sgram.value.tolist(),
            epoch = self.sgram.epoch.value,
            x0 = self.sgram.x0.value,
            dx = self.sgram.dx.value,
            xindex = self.sgram.xindex.value.tolist(),
            y0 = self.sgram.y0.value,
            dy = dy,
            yindex = self.sgram.yindex.value.tolist()
        )
        
    def write(self):
        with h5py.File(self.out_dir + '/spectrogram.h5', 'w') as fd:
            dset = fd.create_dataset('sgram', data=self.sgram.value)
            dset.attrs['x0'] = self.sgram.x0
            dset.attrs['dx'] = self.sgram.dx
            dset.attrs['yindex'] = self.sgram.yindex
    
    def get_spectrogram_plot(self):
        evt = from_gps(self.sgram.x0.value)
        
        spim = Image(self.sgram.T.value, {})
        draw = spim.get_html_draw(w = 650,
                                  h = 350,
                                  y_scale = 'log',
                                  x_range = (0, self.sgram.dx.value * len(self.sgram.xindex)),
                                  y_range = (self.sgram.yindex[0].value, self.sgram.yindex[-1].value),
                                  dw = [ self.sgram.dx.value * len(self.sgram.xindex) ],
                                  dh = [ self.sgram.yindex[-1].value-self.sgram.yindex[0].value ],
                                  y0 = self.sgram.yindex[0].value,
                                  x_label = 'Time [seconds] from %s (%.1f)' % (evt.strftime("%Y-%m-%d %T UTC"), 
                                                                               self.sgram.x0.value),
                                  y_label = 'Frequency [Hz]',
                                  enable_log_cmap = False,
                                  )
        
        return draw
        
        
class StrainProduct:
    def __init__(self, ori_strain, filt_strain=None, out_dir=None):
        self.ori_strain = ori_strain
        self.filt_strain = filt_strain
        self.out_dir = out_dir

    def serialize(self):
        out_list = []
        out_list.append({
            'value': self.ori_strain.value,
            't0': self.ori_strain.t0.value,
            'dt': self.ori_strain.dt.value,
            'name': 'Strain'
            })
        out_list.append({
            'value': self.filt_strain.value,
            't0': self.filt_strain.t0.value,
            'dt': self.filt_strain.dt.value,
            'name': 'Strain_bp'
            })
        return out_list
    
    def write(self):
        self.ori_strain.name = 'Strain'
        self.filt_strain.name = 'Strain_bp'
        self.ori_strain.write(os.path.join(self.out_dir, 'strain.h5'), format='hdf5')  
        if self.filt_strain is not None:
            self.filt_strain.write(os.path.join(self.out_dir, 'strain_bandpassed.h5'), format='hdf5')  
    
    def get_strains_plot(self):
        
        evt = from_gps(self.ori_strain.t0.value)
        times = np.arange(0, 
                          self.ori_strain.dt.value * len(self.ori_strain) + 1,
                          self.ori_strain.dt.value )
       
        hover = HoverTool(tooltips=[("x", "$x"), ("y", "$y")])
                          
        fig = figure(tools=[hover, 'wheel_zoom,box_zoom,pan,reset,save,crosshair'], 
                   plot_height=280,
                   plot_width=700,
                   x_axis_label = 'Time [seconds] from %s (%.1f)' % (evt.strftime("%Y-%m-%d %T UTC"), 
                                                                    self.ori_strain.t0.value),
                   title='Original strain timeseries'
                   )
        ln = fig.line(times, self.ori_strain.value, line_color='blue')

        if self.filt_strain is not None:
            hover1 = HoverTool(tooltips=[("x", "$x"), ("y", "$y")])
                          
            fig1 = figure(tools=[hover, 'wheel_zoom,box_zoom,pan,reset,save,crosshair'], 
                    plot_height=280,
                    plot_width=700,
                    x_axis_label = 'Time [seconds] from %s (%.1f)' % (evt.strftime("%Y-%m-%d %T UTC"), 
                                                                        self.ori_strain.t0.value),
                    x_range = fig.x_range,
                    title='Bandpassed strain timeseries'
                    )
            ln1 = fig1.line(times, self.filt_strain.value, line_color='blue')

            layout = column(fig, fig1)
        else:
            layout = fig

        script, div = components(layout)

        return script, div

class SkymapProduct:
    def __init__(self, asciicat, imagedata, fitsdata, contours_dict=None, out_dir=None):
        self.asciicat = asciicat
        self.imagedata = imagedata
        self.fits_data = fitsdata
        self.contours_dict = contours_dict
        self.out_dir = out_dir if out_dir is not None else '.' 
    
    @staticmethod
    def build_fits_hdu(data_str, head_str):
        head = fits.Header.fromstring(head_str)
        data = np.array(json.loads(data_str))
            
        columns = []
        for i in range(head['TFIELDS']):
            columns.append(fits.Column(array = data[:,i], 
                                        name=head[f'TTYPE{i+1}'], 
                                        format=head[f'TFORM{i+1}'], 
                                        unit=head.get(f'TUNIT{i+1}', None)
                                        )
                            )
        return fits.BinTableHDU.from_columns(columns, head)
        
    def write(self):
        skymaps = {}
        for event in self.fits_data.keys():
            hdu = self.build_fits_hdu(self.fits_data[event]['data'], 
                                 self.fits_data[event]['header'])
            
            skymaps[event] = hdu
            hdu.writeto(os.path.join(self.out_dir, f'{event}_skymap.fits'))
        
        with open(os.path.join(self.out_dir, 'catalog.ecsv'), 'w') as ofd:
            ofd.write(self.asciicat)
        
        return ['catalog.ecsv'] + [f'{x}_skymap.fits' for x in skymaps.keys()]

    def serialize(self):
        _out = {}
        skymaps = {}
        for event in self.fits_data.keys():
            hdu = self.build_fits_hdu(self.fits_data[event]['data'], 
                                      self.fits_data[event]['header'])
            dp = NumpyDataProduct(NumpyDataUnit.from_fits_hdu(hdu), name='skymap_'+event)
            skymaps[event] = dp.encode()
        _out['skymaps'] = skymaps
        if self.contours_dict:
            _out['contours'] = self.contours_dict
        
        return _out
        
    
    def get_plot(self):
        script = ''
        div = f'<img src="data:image/svg+xml;base64,{self.imagedata}">'
        return script, div
    
    def get_catalog_dict(self, api=False):
        catalog_table = aread(self.asciicat)
        
        if not api:   
            with_err  = [x[:-6] for x in catalog_table.columns if 'lower' in x]
            for col in with_err:
                val = catalog_table[col].astype('str')
                low = catalog_table[col+'_lower'].astype('str')
                upp = catalog_table[col+'_upper'].astype('str')
                catalog_table.remove_columns([col, col+'_lower', col+'_upper'])    
                catalog_table.add_column([f'{val[i]}+{upp[i]}{low[i]}' for i in range(len(val))], name=col)
        
        catalog_table.add_column(np.arange(len(catalog_table)), name='index', index=0)
        
        column_lists=[catalog_table[name].tolist() for name in catalog_table.colnames]
        for ID,_col in enumerate(column_lists):
            column_lists[ID] = [x if str(x)!='nan' else None for x in _col]

        catalog_dict = dict(cat_column_list=column_lists,
                cat_column_names=catalog_table.colnames,
                cat_column_descr=catalog_table.dtype.descr,
                cat_meta = catalog_table.meta
                )
        return catalog_dict
    