/**** 成都市多时相裸土光谱复合（Sentinel-2 SR Harmonized） ****/
/**** 研究区：projects/fit-territory-472114-b0/assets/chengdushi ****/

// ==========================
// 0. 参数区
// ==========================
var START_DATE = '2019-01-01';
var END_DATE   = '2021-12-31';

var CLOUD_PROB_THRESHOLD = 40;

// 裸土经验阈值，后续可调
var NDVI_MAX = 0.25;
var NDWI_MAX = 0.00;
var NBR2_MAX = 0.10;
var BSI_MIN  = 0.00;

var EXPORT_SCALE = 10;
var EXPORT_CRS = 'EPSG:32648';   // 成都项目统一投影

// ==========================
// 1. 成都市边界
// ==========================
var chengdu = ee.FeatureCollection('projects/fit-territory-472114-b0/assets/chengdushi');
var aoi = chengdu.geometry();

// ==========================
// 2. 可选耕地掩膜
// ==========================
// 目前先默认整个成都市范围。
// 后续如果你有耕地掩膜，再把这一行替换掉。
var croplandMask = ee.Image.constant(1).clip(aoi).rename('cropland_mask');

// ==========================
// 3. Sentinel-2 数据
// ==========================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 80));

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

var joinedCol = ee.ImageCollection(joined);

// ==========================
// 4. 云/阴影/SCL 掩膜 + 裸土指数
// ==========================
function maskAndAddIndices(img) {
  var cloudProb = ee.Image(img.get('cloud_prob_img')).select('probability');

  var scaled = img.select(
    ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'SCL'],
    ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'SCL']
  ).toFloat();

  var refl = scaled.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']).divide(10000);
  var scl = scaled.select('SCL');

  // 云/阴影/坏像元掩膜
  var notCloud = cloudProb.lt(CLOUD_PROB_THRESHOLD);
  var notBadSCL = scl.neq(3)   // cloud shadow
    .and(scl.neq(8))           // cloud medium probability
    .and(scl.neq(9))           // cloud high probability
    .and(scl.neq(10))          // cirrus
    .and(scl.neq(11));         // snow/ice

  // 去掉水体
  var notWater = scl.neq(6);

  var clean = refl
    .updateMask(notCloud)
    .updateMask(notBadSCL)
    .updateMask(notWater);

  // 指数
  var ndvi = clean.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndwi = clean.normalizedDifference(['B3', 'B8']).rename('NDWI');
  var nbr2 = clean.normalizedDifference(['B11', 'B12']).rename('NBR2');

  var bsi = clean.expression(
    '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))', {
      SWIR1: clean.select('B11'),
      RED: clean.select('B4'),
      NIR: clean.select('B8'),
      BLUE: clean.select('B2')
    }).rename('BSI');

  // 裸土候选掩膜
  var bareMask = ndvi.lt(NDVI_MAX)
    .and(ndwi.lt(NDWI_MAX))
    .and(nbr2.lt(NBR2_MAX))
    .and(bsi.gt(BSI_MIN))
    .and(croplandMask.gt(0));

  return clean
    .addBands(ndvi)
    .addBands(ndwi)
    .addBands(nbr2)
    .addBands(bsi)
    .updateMask(bareMask)
    .copyProperties(img, img.propertyNames());
}

var bareSoilCol = joinedCol.map(maskAndAddIndices);

// ==========================
// 5. 多时相复合
// ==========================
var bareCompositeMedian = bareSoilCol.median().clip(aoi);
var bareCount = bareSoilCol.select('B2').count().rename('bare_obs_count').clip(aoi);

// 主输出：中位数复合 + 有效观测次数
var output = bareCompositeMedian
  .select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'NDVI', 'NDWI', 'NBR2', 'BSI'])
  .addBands(bareCount)
  .clip(aoi);

// ==========================
// 6. 可视化
// ==========================
Map.centerObject(chengdu, 8);

Map.addLayer(
  output.select(['B4', 'B3', 'B2']),
  {min: 0.02, max: 0.30},
  'Chengdu Bare Soil RGB'
);

Map.addLayer(
  output.select('BSI'),
  {min: -0.3, max: 0.4, palette: ['blue', 'white', 'brown']},
  'BSI'
);

Map.addLayer(
  output.select('bare_obs_count'),
  {min: 0, max: 20, palette: ['black', 'purple', 'cyan', 'yellow']},
  'Bare Obs Count'
);

// ==========================
// 7. 导出
// ==========================
Export.image.toDrive({
  image: output,
  description: 'Chengdu_BareSoil_Composite_S2_2019_2021',
  folder: 'GEE_Exports',
  fileNamePrefix: 'Chengdu_BareSoil_Composite_S2_2019_2021',
  region: aoi,
  scale: EXPORT_SCALE,
  crs: EXPORT_CRS,
  maxPixels: 1e13
});