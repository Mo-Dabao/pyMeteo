# -*- coding: utf-8 -*-
"""
GrADS相关

TODO: “dtype”和站点数据的dtype重名

@Time    : 2020/9/20 9:52
@Author  : modabao
"""

import sys
from datetime import datetime

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import xarray as xr


class GrADS(object):
    """
    根据ctl文件解析binary数据文件
    """
    def __init__(self, ctl_path=None, bin_path=None):
        self.ctl_path = ctl_path
        self.bin_path = bin_path
        if ctl_path:
            self.__ctl = {}
            self.__parse()

    def set_ctl_path(self, ctl_path):
        self.ctl_path = ctl_path
        if ctl_path:
            self.__ctl = {}
            self.__parse()

    def set_bin_path(self, bin_path):
        self.bin_path = bin_path

    def __parse(self):
        """
        """
        with open(self.ctl_path) as f:
            ctl_lines = f.readlines()
        self.ctl_lines = ctl_lines
        checked = -1
        for n, line in enumerate(ctl_lines):
            # if this line (has been parsed) or (is empty) or (starts with '*'), skip
            if (n <= checked) or (not (line := line.strip())) or (line[0] == '*'):
                continue
            first_word, *rest_words = line.lower().split(maxsplit=1)
            if first_word == 'options':
                self.__options(*rest_words)
            elif first_word == 'undef':
                self.__undef(*rest_words)
            elif first_word == 'pdef':
                self.__pdef(*rest_words)
            elif first_word in ('xdef', 'ydef'):
                self.__xydef(*rest_words, first_word)
            elif first_word == 'zdef':
                checked = self.__zdef(*rest_words, n)
            elif first_word == 'tdef':
                self.__tdef(*rest_words)
            elif first_word == 'vars':
                checked = self.__vars(*rest_words, n)

    def __options(self, rest_words):
        if 'big_endian' in rest_words:
            self.__ctl['dtype'] = '>f4'
        elif 'little_endian' in rest_words:
            self.__ctl['dtype'] = '<f4'
        elif 'byteswapped' in rest_words:
            self.__ctl['dtype'] = '>f4' if sys.byteorder == 'big' else '<f4'
        if 'template' in rest_words:
            self.__ctl['template'] = True

    def __undef(self, rest_words):
        self.__ctl['undef'] = float(rest_words)

    def __pdef(self, rest_words):
        """
        TODO: 暂时认为wrf的lcc投影都是基于6370000m的正球体
        TODO：`TransformPoint`这个方法在公开地理坐标系下如`WGS84`传参或返回的地理坐标都是先纬度再经度，
              在自定义地理坐标系下传参或返回的地理坐标却是先经度再维度
        """
        isize, jsize, prj, latref, lonref, iref, jref, slat, nlat, standard_lon, dx, dy = \
            rest_words.split()
        assert prj == 'lcc'
        isize, jsize = int(isize), int(jsize)
        latref, lonref, iref, jref, slat, nlat, standard_lon, dx, dy = (
            float(x) for x in [latref, lonref, iref, jref, slat, nlat, standard_lon, dx, dy]
        )
        for line in self.ctl_lines:
            if 'MOAD_CEN_LAT' in line:
                MOAD_CEN_LAT = float(line.split('=')[-1])
            if 'STAND_LON' in line:
                STAND_LON = float(line.split('=')[-1])
        globe = ccrs.Globe(ellipse='sphere', semimajor_axis=6370000, semiminor_axis=6370000)
        lcc = ccrs.LambertConformal(central_longitude=STAND_LON, central_latitude=MOAD_CEN_LAT,
                                    standard_parallels=(slat, nlat), globe=globe)
        geo = lcc.as_geodetic()
        x_ref, y_ref = lcc.transform_point(lonref, latref, geo)
        false_east = (iref - 1) * dx - x_ref
        false_north = (jref - 1) * dy - y_ref
        lcc = ccrs.LambertConformal(central_longitude=STAND_LON, central_latitude=MOAD_CEN_LAT,
                                    false_easting=false_east, false_northing=false_north,
                                    standard_parallels=(slat, nlat), globe=globe)

        def lcc2geo(x, y):
            lonlat = geo.transform_points(lcc, x, y)
            lon, lat = lonlat[..., 0], lonlat[..., 1]
            return lon, lat

        def geo2lcc(lon, lat):
            xyz = lcc.transform_points(geo, lon, lat)
            x, y = xyz[..., 0], xyz[..., 1]
            return x, y

        self.__ctl['pdef'] = {'lcc2geo': lcc2geo, 'geo2lcc': geo2lcc}
        self.__ctl['x_size'], self.__ctl['y_size'] = isize, jsize
        self.__ctl['x'], self.__ctl['y'] = np.arange(isize) * dx, np.arange(jsize) * dy

    def __xydef(self, rest_words, first_word):
        n, linear, start, step = rest_words.split()
        assert linear == 'linear'
        n = int(n)
        start = float(start)
        step = float(step)
        key = 'lat' if first_word == 'ydef' else 'lon'
        self.__ctl[key] = np.linspace(start, start + (n - 1) * step, n)
        if 'pdef' in self.__ctl:
            return
        else:
            self.__ctl[f'{first_word[0]}_size'] = n

    def __zdef(self, rest_words, n):
        znum, levels, *value_list = rest_words.split()
        assert levels == 'levels'
        znum = int(znum)
        if znum == len(value_list):
            value_list = np.fromiter(value_list, dtype=np.float)
            checked = n
        else:
            checked = n + znum
            value_list = np.fromiter(self.ctl_lines[n + 1: checked + 1], np.float)
        self.__ctl['zdef'] = value_list
        return checked

    def __tdef(self, rest_words):
        """
        TODO: not all circumstances have been handled
        Args:
            rest_words:

        Returns:

        """
        ns_t, linear, t_start, dt = rest_words.split()
        assert linear == 'linear'
        ns_t = int(ns_t)
        t_start = datetime.strptime(t_start, '%HZ%d%b%Y')
        unit_map = {'mn': 'min', 'hr': 'H', 'dy': 'D', 'mo': 'M', 'yr': 'Y'}
        dt = f'{dt[:-2]}{unit_map[dt[-2:]]}'
        t = pd.date_range(start=t_start, periods=ns_t, freq=dt)
        self.__ctl['tdef'] = t

    def __vars(self, rest_words, n):
        n_var = int(rest_words)
        variables = [{'name': y[0], 'layers': int(y[1]), 'describe': y[3]}
                     for x in self.ctl_lines[n + 1: n + n_var + 1]
                     if (y := x.split(maxsplit=3))]
        variables = pd.DataFrame(variables).set_index('name')
        self.__ctl['one_time_levels'] = variables['layers'].sum()
        variables['layer_start'] = variables['layers'].cumsum() - variables['layers']
        variables['count'] = variables['layers'] * self.__ctl['x_size'] * self.__ctl['y_size']
        self.__ctl['vars'] = variables
        return n + n_var

    def reset_time(self, time=None, start_time=None):
        """

        Args:
            time: pd.DatetimeIndex
            start_time: datetime.datetime

        Returns:
            None

        """
        assert 'tdef' in self.__ctl
        if isinstance(time, pd.DatetimeIndex):
            self.__ctl['tdef'] = time
        elif start_time:
            self.__ctl['tdef'] += start_time - self.__ctl['tdef'][0]

    def reset_dtype(self, dtype):
        if dtype in ('>f4', '<f4'):
            print(f'dtype changes from "{self.__ctl["dtype"]}" to "{dtype}"')
            self.__ctl['dtype'] = dtype
        else:
            raise ValueError('dtype must in {">f4", "<f4"}')

    @property
    def ctl(self):
        return self.__ctl.copy()

    def get(self, name, time=None):
        """
        TODO: lazy loading using dask
        Args:
            name:
            time:

        Returns:

        """
        if isinstance(time, int):
            return self.get_element_by_time_id(name, time)
        if time is None:
            time = range(len(self.__ctl['tdef']))
        if hasattr(time, '__iter__'):
            element = xr.concat((self.get_element_by_time_id(name, time_id) for time_id in time),
                                dim='time')
            return element

    def get_element_by_time_id(self, name, time_id):
        """
        读取一个要素一个时次

        Args:
            name: 要素名称
            time_id: 时次

        Returns:

        """
        level_start = self.__ctl['vars'].loc[name, 'layer_start']
        levels = self.__ctl['vars'].loc[name, 'layers']
        count = self.__ctl['vars'].loc[name, 'count']
        offset = ((time_id * self.__ctl['one_time_levels'] + level_start) *
                  self.__ctl['x_size'] * self.__ctl['y_size'] * 4)
        shape = 1, levels, self.__ctl['y_size'], self.__ctl['x_size']
        element = np.fromfile(self.bin_path, dtype=self.__ctl['dtype'],
                              count=count, offset=offset).reshape(shape)
        y, x = ('y', 'x') if 'pdef' in self.__ctl else ('lat', 'lon')
        level = [0] if self.__ctl['vars'].loc[name, 'layers'] == 1 else self.__ctl['zdef']
        coords = (('time', self.__ctl['tdef'][[time_id]]), ('level', level),
                  (y, self.__ctl[y]), (x, self.__ctl[x]))
        element = xr.DataArray(element, coords=coords, name=name)
        return element.where(element != self.__ctl['undef'])
