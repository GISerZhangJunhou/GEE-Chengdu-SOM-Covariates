/**** 成都市管理代理变量（250m导出版） ****/
/**** 输出：
  1) 多年灌溉强度均值 / 标准差 / 变异系数 / 趋势
  2) 稻作频率（来自 ChinaCP 2015–2021）
  3) 生长季 LSWI 统计量与趋势
****/

// ==========================
// 0. 参数
// ==========================
var PROJECT_PATH = 'projects/fit-territory-472114-b0/assets/';
var EXPORT_SCALE = 250;
var EXPORT_CRS = 'EPSG:32648';

// Sentinel-2 时段
var START_DATE = '2019-01-01';
var END_DATE   = '2021-12-31';
var GROW_START_MONTH = 5;
var GROW_END_MONTH   = 10;
var CLOUD_PROB_THRESHOLD = 40;

// ==========================
// 1. 研究区
// ==========================
var chengdu = ee.FeatureCollection(PROJECT_PATH + 'chengdushi');
var aoi = chengdu.geometry();

// ==========================
// 2. ChinaCP：2015–2021
// ==========================
function getCP(year) {
  return ee.Image(PROJECT_PATH + 'CDlunzuo' + year).rename('cp');
}

var cpYears = [2015, 2016, 2017, 2018, 2019, 2020, 2021];

var cpCol = ee.ImageCollection.fromImages(
  cpYears.map(function(year) {
    return getCP(year).set('year', year).clip(aoi).toInt16();
  })
);

// ChinaCP 含水稻编码
// 0=休耕地
// 3=三作_水稻+玉米+小麦
// 14=单作_玉米
// 15=单作_水稻
// 16=单作_小麦
// 17=单作_其他作物
// 27=双作_其他作物
// 245=双作_水稻+玉米
// 246=双作_小麦+玉米
// 255=双作_水稻+水稻
// 256=双作_小麦+水稻
var riceCodes = [3, 15, 245, 255, 256];

function toBinaryByCodes(img, codes, bandName) {
  var out = ee.Image(0);
  codes.forEach(function(code) {
    out = out.or(img.eq(code));
  });
  return out.rename(bandName).toFloat();
}

// 稻作频率
var freqRice = cpCol.map(function(img) {
  return toBinaryByCodes(img.select('cp'), riceCodes, 'rice');
}).mean().rename('freq_rice');

// ==========================
// 3. 灌溉：2000–2020
// ==========================
// 如果你的资产名不是 CIrrMap250_2000 这种格式，改这里的前缀即可
function getIrr(year) {
  return ee.Image(PROJECT_PATH + 'CIrrMap250_' + year).rename('irr');
}

var irrYears = [
  2000, 2001, 2002, 2003, 2004, 2005, 2006,
  2007, 2008, 2009, 2010, 2011, 2012, 2013,
  2014, 2015, 2016, 2017, 2018, 2019, 2020
];

var irrCol = ee.ImageCollection.fromImages(
  irrYears.map(function(year) {
    return getIrr(year).set('year', year).clip(aoi).toFloat();
  })
);

function addYearBand(img) {
  var year = ee.Number(img.get('year'));
  var yearBand = ee.Image.constant(year).rename('year').toFloat();
  return yearBand.addBands(img.toFloat());
}

// 多年均值 / 标准差 / 变异系数
var irrMean = irrCol.mean().rename('irr_mean');
var irrStd = irrCol.reduce(ee.Reducer.stdDev()).rename('irr_std');
var irrCV = irrStd.divide(irrMean.add(1e-6)).rename('irr_cv');

// 趋势
var irrTrend = irrCol.map(addYearBand)
  .select(['year', 'irr'])
  .reduce(ee.Reducer.linearFit())
  .select('scale')
  .rename('irr_trend');

// ==========================
// 4. LSWI：Sentinel-2 生长季统计
// ==========================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 80));

var s2CloudProb = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
  .filterBounds(aoi)
  .filterDate(START_DATE, END_DATE);

var joined = ee.Join.saveFirst('cloud_prob_img').apply({
  primary: s2,
  secondary: s2CloudProb,
  condition: ee.Filter.equals({
    leftField: 'system:index',
    rightField: 'system:index'
  })
});

var joinedCol = ee.ImageCollection(joined);

function maskAndAddLSWI(img) {
  var cloudProb = ee.Image(img.get('cloud_prob_img')).select('probability');

  var scaled = img.select(
    ['B8', 'B11', 'SCL'],
    ['B8', 'B11', 'SCL']
  ).toFloat();

  var refl = scaled.select(['B8', 'B11']).divide(10000);
  var scl = scaled.select('SCL');

  var notCloud = cloudProb.lt(CLOUD_PROB_THRESHOLD);
  var notBadSCL = scl.neq(3)
    .and(scl.neq(8))
    .and(scl.neq(9))
    .and(scl.neq(10))
    .and(scl.neq(11));

  var clean = refl.updateMask(notCloud).updateMask(notBadSCL);

  var lswi = clean.normalizedDifference(['B8', 'B11']).rename('LSWI');

  return clean
    .addBands(lswi)
    .copyProperties(img, img.propertyNames());
}

var lswiCol = joinedCol
  .filter(ee.Filter.calendarRange(GROW_START_MONTH, GROW_END_MONTH, 'month'))
  .map(maskAndAddLSWI);

// 生长季整体统计
var lswiMedian = lswiCol.select('LSWI').median().rename('lswi_grow_median');
var lswiP90 = lswiCol.select('LSWI')
  .reduce(ee.Reducer.percentile([90]))
  .rename('lswi_grow_p90');
var lswiStd = lswiCol.select('LSWI')
  .reduce(ee.Reducer.stdDev())
  .rename('lswi_grow_std');

// 每年生长季中位数，再做趋势
function yearlyGrowMedian(year) {
  year = ee.Number(year);

  var yearly = lswiCol
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .select('LSWI')
    .median()
    .rename('LSWI')
    .set('year', year);

  var yearBand = ee.Image.constant(year).rename('year').toFloat();
  return yearBand.addBands(yearly.toFloat());
}

var lswiYears = ee.List.sequence(2019, 2021);

var yearlyLswiCol = ee.ImageCollection.fromImages(
  lswiYears.map(yearlyGrowMedian)
);

var lswiTrend = yearlyLswiCol
  .select(['year', 'LSWI'])
  .reduce(ee.Reducer.linearFit())
  .select('scale')
  .rename('lswi_trend');

// ==========================
// 5. 合并输出
// ==========================
var output = ee.Image.cat([
  freqRice,
  irrMean,
  irrStd,
  irrCV,
  irrTrend,
  lswiMedian,
  lswiP90,
  lswiStd,
  lswiTrend
]).toFloat().clip(aoi);

// ==========================
// 6. 统一到 250m 导出
// ==========================
// 这些输出都是连续变量，所以这里用 bilinear
var output250 = output
  .resample('bilinear')
  .reproject({
    crs: EXPORT_CRS,
    scale: EXPORT_SCALE
  })
  .clip(aoi);

// ==========================
// 7. 可视化
// ==========================
Map.centerObject(chengdu, 8);
Map.addLayer(freqRice, {min: 0, max: 1, palette: ['white', 'blue']}, 'freq_rice');
Map.addLayer(irrMean, {min: 0, max: 100, palette: ['white', 'cyan', 'blue']}, 'irr_mean');
Map.addLayer(lswiMedian, {min: -0.5, max: 0.5, palette: ['brown', 'white', 'green']}, 'lswi_grow_median');

// ==========================
// 8. 导出
// ==========================
Export.image.toDrive({
  image: output250,
  description: 'Chengdu_Management_Proxy_Vars_250m',
  folder: 'GEE_Exports',
  fileNamePrefix: 'Chengdu_Management_Proxy_Vars_250m',
  region: aoi,
  scale: EXPORT_SCALE,
  crs: EXPORT_CRS,
  maxPixels: 1e13
});