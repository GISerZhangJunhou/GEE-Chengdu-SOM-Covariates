/**** 成都市 2015–2021 ChinaCP 多年作物制度变量 ****/

// ==========================
// 1. 研究区
// ==========================
var chengdu = ee.FeatureCollection('projects/fit-territory-472114-b0/assets/chengdushi');
var aoi = chengdu.geometry();

// ==========================
// 2. 年度 ChinaCP 资产
// ==========================
var cp2015 = ee.Image('projects/fit-territory-472114-b0/assets/CDlunzuo2015').rename('cp');
var cp2016 = ee.Image('projects/fit-territory-472114-b0/assets/CDlunzuo2016').rename('cp');
var cp2017 = ee.Image('projects/fit-territory-472114-b0/assets/CDlunzuo2017').rename('cp');
var cp2018 = ee.Image('projects/fit-territory-472114-b0/assets/CDlunzuo2018').rename('cp');
var cp2019 = ee.Image('projects/fit-territory-472114-b0/assets/CDlunzuo2019').rename('cp');
var cp2020 = ee.Image('projects/fit-territory-472114-b0/assets/CDlunzuo2020').rename('cp');
var cp2021 = ee.Image('projects/fit-territory-472114-b0/assets/CDlunzuo2021').rename('cp');

var cpCol = ee.ImageCollection.fromImages([
  cp2015.set('year', 2015),
  cp2016.set('year', 2016),
  cp2017.set('year', 2017),
  cp2018.set('year', 2018),
  cp2019.set('year', 2019),
  cp2020.set('year', 2020),
  cp2021.set('year', 2021)
]).map(function(img) {
  return img.clip(aoi).toInt16();
});

// ==========================
// 3. 编码定义
// ==========================
/*
0   = 休耕地
3   = 三作_水稻+玉米+小麦
14  = 单作_玉米
15  = 单作_水稻
16  = 单作_小麦
17  = 单作_其他作物
27  = 双作_其他作物
245 = 双作_水稻+玉米
246 = 双作_小麦+玉米
255 = 双作_水稻+水稻
256 = 双作_小麦+水稻
*/

var singleCodes = [14, 15, 16, 17];
var doubleCodes = [27, 245, 246, 255, 256];
var tripleCodes = [3];

var riceCodes  = [3, 15, 245, 255, 256];
var wheatCodes = [3, 16, 246, 256];
var maizeCodes = [3, 14, 245, 246];

// ==========================
// 4. 工具函数
// ==========================
function toBinaryByCodes(img, codes, bandName) {
  var out = ee.Image(0);
  codes.forEach(function(code) {
    out = out.or(img.eq(code));
  });
  return out.rename(bandName).toFloat();
}

// 复种指数：休耕=0，单作=1，双作=2，三作=3
function toCroppingIntensity(img) {
  var singleMask = toBinaryByCodes(img, singleCodes, 'single');
  var doubleMask = toBinaryByCodes(img, doubleCodes, 'double');
  var tripleMask = toBinaryByCodes(img, tripleCodes, 'triple');

  var intensity = ee.Image(0)
    .where(singleMask.eq(1), 1)
    .where(doubleMask.eq(1), 2)
    .where(tripleMask.eq(1), 3)
    .rename('cropping_intensity')
    .toFloat();

  return intensity;
}

// 轮作稳定性：7年中最常见类别所占比例
function calcRotationStability(ic) {
  var modeImg = ic.mode().rename('cp_mode');
  var agreeCol = ic.map(function(img) {
    return img.eq(modeImg).rename('agree').toFloat();
  });
  return agreeCol.mean().rename('rotation_stability');
}

// ==========================
// 5. 多年频率变量
// ==========================
var freqSingle = cpCol.map(function(img) {
  return toBinaryByCodes(img.select('cp'), singleCodes, 'single');
}).mean().rename('freq_single');

var freqDouble = cpCol.map(function(img) {
  return toBinaryByCodes(img.select('cp'), doubleCodes, 'double');
}).mean().rename('freq_double');

var freqTriple = cpCol.map(function(img) {
  return toBinaryByCodes(img.select('cp'), tripleCodes, 'triple');
}).mean().rename('freq_triple');

var freqRice = cpCol.map(function(img) {
  return toBinaryByCodes(img.select('cp'), riceCodes, 'rice');
}).mean().rename('freq_rice');

var freqWheat = cpCol.map(function(img) {
  return toBinaryByCodes(img.select('cp'), wheatCodes, 'wheat');
}).mean().rename('freq_wheat');

var freqMaize = cpCol.map(function(img) {
  return toBinaryByCodes(img.select('cp'), maizeCodes, 'maize');
}).mean().rename('freq_maize');

// ==========================
// 6. 复种指数
// ==========================
var intensityCol = cpCol.map(function(img) {
  return toCroppingIntensity(img.select('cp'));
});

var meanCroppingIntensity = intensityCol.mean().rename('mean_cropping_intensity');
var stdCroppingIntensity = intensityCol.reduce(ee.Reducer.stdDev()).rename('std_cropping_intensity');

// ==========================
// 7. 轮作稳定性
// ==========================
var rotationStability = calcRotationStability(cpCol.select('cp'));
var rotationAlwaysSame = rotationStability.eq(1).rename('rotation_always_same').toFloat();

// ==========================
// 8. 合并输出
// ==========================
var output = ee.Image.cat([
  freqSingle,
  freqDouble,
  freqTriple,
  freqRice,
  freqWheat,
  freqMaize,
  meanCroppingIntensity,
  stdCroppingIntensity,
  rotationStability,
  rotationAlwaysSame
]).toFloat().clip(aoi);

// ==========================
// 9. 可视化
// ==========================
Map.centerObject(chengdu, 8);

Map.addLayer(freqSingle, {min: 0, max: 1, palette: ['white', 'green']}, 'freq_single');
Map.addLayer(freqDouble, {min: 0, max: 1, palette: ['white', 'orange']}, 'freq_double');
Map.addLayer(freqTriple, {min: 0, max: 1, palette: ['white', 'red']}, 'freq_triple');
Map.addLayer(freqRice, {min: 0, max: 1, palette: ['white', 'blue']}, 'freq_rice');
Map.addLayer(meanCroppingIntensity, {min: 0, max: 3, palette: ['white', 'yellow', 'orange', 'red']}, 'mean_cropping_intensity');
Map.addLayer(rotationStability, {min: 0, max: 1, palette: ['blue', 'white', 'darkgreen']}, 'rotation_stability');

// ==========================
// 10. 导出
// ==========================
Export.image.toDrive({
  image: output,
  description: 'Chengdu_ChinaCP_Multiyear_2015_2021',
  folder: 'GEE_Exports',
  fileNamePrefix: 'Chengdu_ChinaCP_Multiyear_2015_2021',
  region: aoi,
  scale: 500,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});