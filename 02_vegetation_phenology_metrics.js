/**** 成都市耕地 SOM 项目：第4类变量 多时相植被/物候变量（正式版） ****/
/**** 时段：2019-2021 ****/
/**** 输出分辨率：250 m ****/
/**** 输出坐标：EPSG:32648 ****/

// =====================================================
// 0. 参数
// =====================================================
var PROJECT_PATH = 'projects/fit-territory-472114-b0/assets/';
var EXPORT_SCALE = 250;
var EXPORT_CRS = 'EPSG:32648';

// 给 reduceResolution 之前的影像统一一个 10 m 默认投影
var BASE_SCALE = 10;
var BASE_PROJ = ee.Projection(EXPORT_CRS).atScale(BASE_SCALE);

var START_DATE = '2019-01-01';
var END_DATE   = '2021-12-31';
var YEARS = [2019, 2020, 2021];

// 成都先用 4–10 月作为生长季
var GROW_START_MONTH = 4;
var GROW_END_MONTH   = 10;

var MAX_CLOUDY_PIXEL_PERCENTAGE = 80;
var CLOUD_PROB_THRESHOLD = 40;

// =====================================================
// 1. 研究区
// =====================================================
var chengdu = ee.FeatureCollection(PROJECT_PATH + 'chengdushi');
var aoi = chengdu.geometry();

Map.centerObject(chengdu, 8);

// =====================================================
// 2. Sentinel-2 数据
// =====================================================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', MAX_CLOUDY_PIXEL_PERCENTAGE));

var s2CloudProb = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
  .filterBounds(aoi)
  .filterDate(START_DATE, END_DATE);

// 按 system:index 连接云概率
var joined = ee.Join.saveFirst('cloud_prob_img').apply({
  primary: s2,
  secondary: s2CloudProb,
  condition: ee.Filter.equals({
    leftField: 'system:index',
    rightField: 'system:index'
  })
});

var joinedCol = ee.ImageCollection(joined)
  .filter(ee.Filter.notNull(['cloud_prob_img']));

print('S2 matched image count:', joinedCol.size());

// =====================================================
// 3. 工具函数
// =====================================================
function forceProj(img) {
  return ee.Image(img).toFloat().setDefaultProjection(BASE_PROJ);
}

function addPropsAndProj(img, srcImg) {
  var out = ee.Image(img);
  out = ee.Image(out.copyProperties(srcImg, ['system:time_start', 'system:index']));
  out = out.setDefaultProjection(BASE_PROJ);
  return out;
}

// =====================================================
// 4. 云掩膜 + 植被指数 + 时间波段
// =====================================================
function maskAndAddVegIndices(img) {
  var cloudProb = ee.Image(img.get('cloud_prob_img')).select('probability');
  var scl = img.select('SCL');

  var goodMask = cloudProb.lt(CLOUD_PROB_THRESHOLD)
    .and(scl.neq(0))   // No Data
    .and(scl.neq(1))   // Saturated / defective
    .and(scl.neq(3))   // Cloud shadow
    .and(scl.neq(8))   // Cloud medium probability
    .and(scl.neq(9))   // Cloud high probability
    .and(scl.neq(10))  // Cirrus
    .and(scl.neq(11)); // Snow / ice

  var refl = img.select(['B2', 'B4', 'B8'])
    .toFloat()
    .divide(10000)
    .updateMask(goodMask);

  var ndvi = refl.normalizedDifference(['B8', 'B4']).rename('NDVI').toFloat();

  var evi = refl.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
    {
      NIR: refl.select('B8'),
      RED: refl.select('B4'),
      BLUE: refl.select('B2')
    }
  ).rename('EVI').toFloat();

  var doyValue = ee.Number(
    ee.Date(img.get('system:time_start')).getRelative('day', 'year')
  ).add(1);

  var monthValue = ee.Number(
    ee.Date(img.get('system:time_start')).get('month')
  );

  var doy = ndvi.multiply(0).add(doyValue).rename('DOY').toFloat();
  var month = ndvi.multiply(0).add(monthValue).rename('MONTH').toFloat();

  var out = ee.Image.cat([ndvi, evi, doy, month]).toFloat();
  out = addPropsAndProj(out, img);

  return out;
}

// 全年集合：用于振幅、峰值时间
var vegColAll = joinedCol.map(maskAndAddVegIndices);

// 生长季集合：用于生长季统计
var vegColGrow = vegColAll.filter(
  ee.Filter.calendarRange(GROW_START_MONTH, GROW_END_MONTH, 'month')
);

print('All-season image count:', vegColAll.size());
print('Growing-season image count:', vegColGrow.size());

// =====================================================
// 5. 2019-2021 整体生长季统计
// =====================================================
var ndviGrowMean = forceProj(
  vegColGrow.select('NDVI').mean().rename('ndvi_grow_mean')
);

var ndviGrowMedian = forceProj(
  vegColGrow.select('NDVI').median().rename('ndvi_grow_median')
);

var ndviGrowP90 = forceProj(
  vegColGrow.select('NDVI').reduce(ee.Reducer.percentile([90])).rename('ndvi_grow_p90')
);

var ndviGrowStd = forceProj(
  vegColGrow.select('NDVI').reduce(ee.Reducer.stdDev()).rename('ndvi_grow_std')
);

var eviGrowMean = forceProj(
  vegColGrow.select('EVI').mean().rename('evi_grow_mean')
);

var eviGrowMedian = forceProj(
  vegColGrow.select('EVI').median().rename('evi_grow_median')
);

var eviGrowP90 = forceProj(
  vegColGrow.select('EVI').reduce(ee.Reducer.percentile([90])).rename('evi_grow_p90')
);

var eviGrowStd = forceProj(
  vegColGrow.select('EVI').reduce(ee.Reducer.stdDev()).rename('evi_grow_std')
);

// =====================================================
// 6. 各年份年度统计
// =====================================================
function yearlyVegStats(year) {
  var y = Number(year);

  var yearlyAll = vegColAll.filter(ee.Filter.calendarRange(y, y, 'year'));
  var yearlyGrow = vegColGrow.filter(ee.Filter.calendarRange(y, y, 'year'));

  // 年度生长季中位数：供趋势用
  var ndviMedianY = forceProj(
    yearlyGrow.select('NDVI').median().rename('NDVI_median')
  );

  var eviMedianY = forceProj(
    yearlyGrow.select('EVI').median().rename('EVI_median')
  );

  // 年振幅：基于全年
  var ndviMaxY = forceProj(
    yearlyAll.select('NDVI').max().rename('NDVI_max')
  );
  var ndviMinY = forceProj(
    yearlyAll.select('NDVI').min().rename('NDVI_min')
  );
  var ndviAmplitudeY = forceProj(
    ndviMaxY.subtract(ndviMinY).rename('NDVI_amplitude')
  );

  var eviMaxY = forceProj(
    yearlyAll.select('EVI').max().rename('EVI_max')
  );
  var eviMinY = forceProj(
    yearlyAll.select('EVI').min().rename('EVI_min')
  );
  var eviAmplitudeY = forceProj(
    eviMaxY.subtract(eviMinY).rename('EVI_amplitude')
  );

  // 峰值时间：基于全年
  var ndviPeakComp = yearlyAll.qualityMosaic('NDVI');
  var eviPeakComp = yearlyAll.qualityMosaic('EVI');

  var ndviPeakDoyY = forceProj(
    ndviPeakComp.select('DOY').rename('NDVI_peak_doy')
  );
  var eviPeakDoyY = forceProj(
    eviPeakComp.select('DOY').rename('EVI_peak_doy')
  );

  var ndviPeakMonthY = forceProj(
    ndviPeakComp.select('MONTH').rename('NDVI_peak_month')
  );
  var eviPeakMonthY = forceProj(
    eviPeakComp.select('MONTH').rename('EVI_peak_month')
  );

  var yearBand = forceProj(
    ndviMedianY.multiply(0).add(y).rename('year')
  );

  var yearlyOut = ee.Image.cat([
    yearBand,
    ndviMedianY,
    eviMedianY,
    ndviAmplitudeY,
    eviAmplitudeY,
    ndviPeakDoyY,
    eviPeakDoyY,
    ndviPeakMonthY,
    eviPeakMonthY
  ]).toFloat();

  yearlyOut = yearlyOut.setDefaultProjection(BASE_PROJ);
  yearlyOut = ee.Image(yearlyOut.set('year', y));

  return yearlyOut;
}

var yearlyStatsCol = ee.ImageCollection.fromImages(
  YEARS.map(yearlyVegStats)
);

print('Yearly stats image count:', yearlyStatsCol.size());

// =====================================================
// 7. 趋势
// =====================================================
var ndviTrend = forceProj(
  yearlyStatsCol
    .select(['year', 'NDVI_median'])
    .reduce(ee.Reducer.linearFit())
    .select('scale')
    .rename('ndvi_trend')
);

var eviTrend = forceProj(
  yearlyStatsCol
    .select(['year', 'EVI_median'])
    .reduce(ee.Reducer.linearFit())
    .select('scale')
    .rename('evi_trend')
);

// =====================================================
// 8. 多年统计：年度中位数均值/波动
// =====================================================
var ndviAnnualMedianMean = forceProj(
  yearlyStatsCol.select('NDVI_median').mean().rename('ndvi_annual_median_mean')
);

var ndviAnnualMedianStd = forceProj(
  yearlyStatsCol.select('NDVI_median').reduce(ee.Reducer.stdDev()).rename('ndvi_annual_median_std')
);

var eviAnnualMedianMean = forceProj(
  yearlyStatsCol.select('EVI_median').mean().rename('evi_annual_median_mean')
);

var eviAnnualMedianStd = forceProj(
  yearlyStatsCol.select('EVI_median').reduce(ee.Reducer.stdDev()).rename('evi_annual_median_std')
);

// =====================================================
// 9. 多年统计：振幅
// =====================================================
var ndviAmplitudeMean = forceProj(
  yearlyStatsCol.select('NDVI_amplitude').mean().rename('ndvi_amplitude_mean')
);

var ndviAmplitudeStd = forceProj(
  yearlyStatsCol.select('NDVI_amplitude').reduce(ee.Reducer.stdDev()).rename('ndvi_amplitude_std')
);

var eviAmplitudeMean = forceProj(
  yearlyStatsCol.select('EVI_amplitude').mean().rename('evi_amplitude_mean')
);

var eviAmplitudeStd = forceProj(
  yearlyStatsCol.select('EVI_amplitude').reduce(ee.Reducer.stdDev()).rename('evi_amplitude_std')
);

// =====================================================
// 10. 多年统计：峰值 DOY
// =====================================================
var ndviPeakDoyMean = forceProj(
  yearlyStatsCol.select('NDVI_peak_doy').mean().rename('ndvi_peak_doy_mean')
);

var ndviPeakDoyStd = forceProj(
  yearlyStatsCol.select('NDVI_peak_doy').reduce(ee.Reducer.stdDev()).rename('ndvi_peak_doy_std')
);

var eviPeakDoyMean = forceProj(
  yearlyStatsCol.select('EVI_peak_doy').mean().rename('evi_peak_doy_mean')
);

var eviPeakDoyStd = forceProj(
  yearlyStatsCol.select('EVI_peak_doy').reduce(ee.Reducer.stdDev()).rename('evi_peak_doy_std')
);

// =====================================================
// 11. 多年统计：峰值月份
// =====================================================
var ndviPeakMonthMean = forceProj(
  yearlyStatsCol.select('NDVI_peak_month').mean().rename('ndvi_peak_month_mean')
);

var ndviPeakMonthStd = forceProj(
  yearlyStatsCol.select('NDVI_peak_month').reduce(ee.Reducer.stdDev()).rename('ndvi_peak_month_std')
);

var eviPeakMonthMean = forceProj(
  yearlyStatsCol.select('EVI_peak_month').mean().rename('evi_peak_month_mean')
);

var eviPeakMonthStd = forceProj(
  yearlyStatsCol.select('EVI_peak_month').reduce(ee.Reducer.stdDev()).rename('evi_peak_month_std')
);

// =====================================================
// 12. 合并所有输出波段
// =====================================================
var output = ee.Image.cat([
  ndviGrowMean,
  ndviGrowMedian,
  ndviGrowP90,
  ndviGrowStd,
  eviGrowMean,
  eviGrowMedian,
  eviGrowP90,
  eviGrowStd,

  ndviTrend,
  eviTrend,

  ndviAnnualMedianMean,
  ndviAnnualMedianStd,
  eviAnnualMedianMean,
  eviAnnualMedianStd,

  ndviAmplitudeMean,
  ndviAmplitudeStd,
  eviAmplitudeMean,
  eviAmplitudeStd,

  ndviPeakDoyMean,
  ndviPeakDoyStd,
  eviPeakDoyMean,
  eviPeakDoyStd,

  ndviPeakMonthMean,
  ndviPeakMonthStd,
  eviPeakMonthMean,
  eviPeakMonthStd
]).toFloat();

output = output.clip(aoi);
output = output.setDefaultProjection(BASE_PROJ);

// =====================================================
// 13. 聚合到 250 m
// =====================================================
var output250 = output
  .reduceResolution({
    reducer: ee.Reducer.mean(),
    maxPixels: 4096
  })
  .reproject({
    crs: EXPORT_CRS,
    scale: EXPORT_SCALE
  })
  .clip(aoi)
  .toFloat();

output250 = output250.setDefaultProjection(
  ee.Projection(EXPORT_CRS).atScale(EXPORT_SCALE)
);

// =====================================================
// 14. 检查
// =====================================================
print('Final band names:', output250.bandNames());
print('Final band count:', output250.bandNames().size());

Map.addLayer(
  output250.select('ndvi_grow_median'),
  {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green']},
  'ndvi_grow_median_250m'
);

Map.addLayer(
  output250.select('evi_grow_median'),
  {min: 0, max: 0.6, palette: ['brown', 'yellow', 'darkgreen']},
  'evi_grow_median_250m'
);

Map.addLayer(
  output250.select('ndvi_amplitude_mean'),
  {min: 0, max: 0.8, palette: ['white', 'orange', 'red']},
  'ndvi_amplitude_mean_250m'
);

Map.addLayer(
  output250.select('ndvi_peak_doy_mean'),
  {min: 100, max: 320, palette: ['blue', 'cyan', 'yellow', 'red']},
  'ndvi_peak_doy_mean_250m'
);

// =====================================================
// 15. 导出
// =====================================================
Export.image.toDrive({
  image: output250,
  description: 'Chengdu_VegPhenology_Vars_2019_2021_250m',
  folder: 'GEE_Exports',
  fileNamePrefix: 'Chengdu_VegPhenology_Vars_2019_2021_250m',
  region: aoi,
  scale: EXPORT_SCALE,
  crs: EXPORT_CRS,
  maxPixels: 1e13
});