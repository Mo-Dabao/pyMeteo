# -*- coding: utf-8 -*-
"""
从FY-4A标称数据提取指定范围指定通道

@Time    : 2018/11/14 12:46:47
@Author  : modabao
"""

import xarray as xr
import numpy as np

from projection import latlon2linecolumn


# 各分辨率文件包含的通道号
CONTENTS = {'0500M': ('Channel02',),
            '1000M': ('Channel01', 'Channel02', 'Channel03'),
            '2000M': tuple([f'Channel{x:02d}' for x in range(1, 8)]),
            '4000M': tuple([f'Channel{x:02d}' for x in range(1, 15)])}
# 各分辨率全圆盘数据的行列数
SIZES = {'0500M': 21984,
         '1000M': 10992,
         '2000M': 5496,
         '4000M': 2748}


class AGRI_L1(object):
    """
    FY4A AGRI LEVEL1数据按
    """
    def __init__(self, file_path, geo_desc=None):
        """
        获得L1数据hdf5文件对象
        """
        self.dataset = xr.open_dataset(file_path)
        self.resolution = file_path[-15:-10]
        self.line_begin = self.dataset.attrs['Begin Line Number']
        self.line_end = self.dataset.attrs['End Line Number']
        self.column_begin = self.dataset.attrs['Begin Pixel Number']
        self.column_end = self.dataset.attrs['End Pixel Number']
        self.set_geo_desc(geo_desc)

    def __del__(self):
        """
        确保关闭L1数据hdf5文件
        """
        self.dataset.close()
        
    def set_geo_desc(self, geo_desc):
        if geo_desc is None:
            self.line = self.column = self.geo_desc = None
            return
        # 先乘1000取整是为了减少浮点数的精度误差累积问题
        lat_S, lat_N, lon_W, lon_E, step = [1000 * x for x in geo_desc]
        lat = np.arange(lat_N, lat_S-1, -step) / 1000
        lon = np.arange(lon_W, lon_E+1, step) / 1000
        lon_mesh, lat_mesh = np.meshgrid(lon, lat)
        # 求geo_desc对应的标称全圆盘行列号
        line, column = latlon2linecolumn(lat_mesh, lon_mesh, self.resolution)
        self.line = xr.DataArray(line, coords=(('lat', lat), ('lon', lon)), name='line')
        self.column = xr.DataArray(column, coords=(('lat', lat), ('lon', lon)), name='column')
        self.geo_desc = geo_desc

    def extract(self, channel_name, calibration='reflectance',
                geo_desc=None, interp_method='nearest'):
        """
        按通道名和定标方式提取geo_desc对应的数据
        channel_name：要提取的通道名（如'Channel01'）
        
        calibration: {'dn', 'reflectance', 'radiance', 'brightness_temperature'}
        """
        if geo_desc and geo_desc != self.geo_desc:
            self.set_geo_desc(geo_desc)
        dn_values = self.dataset[f'NOM{channel_name}']
        dn_values = dn_values.rename({dn_values.dims[0]: 'line', dn_values.dims[1]: 'column'})
        dn_values = dn_values.assign_coords(line=range(self.line_begin, self.line_end+1),
                                            column=range(self.column_begin, self.column_end+1))
        if self.geo_desc:
            # 若geo_desc已指定，则插值到对应网格
            dn_values = dn_values.interp(line=self.line, column=self.column, method=interp_method)
            del dn_values.coords['line'], dn_values.coords['column']
        else:
            # 若geo_desc为None，则保留原始NOM网格
            pass
        return self.calibrate(channel_name, calibration, dn_values)

    def calibrate(self, channel_name, calibration, dn_values):
        """
        前面6个通道，用查找表和系数算出来都是反射率，后面用查找表是亮温，用系数是辐亮度。
        """
        if calibration == 'dn':
            dn_values.attrs = {'units': 'DN'}
            return dn_values
        channel_num = int(channel_name[-2:])
        dn_values = dn_values.fillna(dn_values.FillValue)  # 保留缺省值
        if ((calibration == 'reflectance' and channel_num <= 6) or
            (calibration == 'radiance' and channel_num > 6)):
            k, b = self.dataset['CALIBRATION_COEF(SCALE+OFFSET)'].values[channel_num-1]
            data = k * dn_values.where(dn_values != dn_values.FillValue) + b
            data.attrs['units'] = '100%' if calibration == 'reflectance' else 'mW/ (m2 cm-1 sr)'
        elif calibration == 'brightness_temperature' and channel_num > 6:
            cal_table = self.dataset[f'CAL{channel_name}']
            cal_table = cal_table.swap_dims({cal_table.dims[0]: 'dn'})
            data = cal_table.interp(dn=dn_values)
            del data.coords['dn']
            data.attrs = {'units': 'K'}
        else:
            raise ValueError(f'{channel_name}没有{calibration}的定标方式')
        data.name = f'{channel_name}_{calibration}'
        return data


