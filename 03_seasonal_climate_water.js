/**** 成都市耕地 SOM 项目：第5类变量 季节性气候/水热约束变量（正式版） ****/
/**** 时段：2019-2021 ****/
/**** 数据源：TerraClimate（月尺度） ****/
/**** 输出分辨率：250 m（用于与现有主表对齐，不代表原始气候分辨率） ****/
/**** 输出坐标：EPSG:32648 ****/

// =====================================================
// 0. 参数
// =====================================================
var PROJECT_PATH = 'projects/fit-territory-472114-b0/assets/';
var EXPORT_SCALE = 250;
var EXPORT_CRS = 'EPSG:32648';

var START_DATE = '2019-01-01';
var END_DATE   = '2021-12-31';
var YEARS = [2019, 2020, 2021];

// 与前面植被/物候变量保持一致
var GROW_START_MONTH = 4;
var GROW_END_MONTH   = 10;

// =====================================================
// 1. 研究区
// =====================================================
var chengdu = ee.FeatureCollection(PROJECT_PATH + 'chengdushi');
var aoi = chengdu.geometry();

Map.centerObject(chengdu, 8);

// =====================================================
// 2. TerraClimate 数据
// =====================================================
var tcRaw = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
  .filterBounds(aoi)
  .filterDate(START_DATE, END_DATE);

print('TerraClimate monthly image count:', tcRaw.size());

// =====================================================
// 3. 波段预处理
// 单位换算：
// tmmx/tmmn  -> *0.1  (°C)
// pet/aet/def/soil -> *0.1 (mm)
// vpd -> *0.01 (kPa)
// pr 原单位就是 mm
// =====================================================
function prepClimate(img) {
  var tmean = img.expression(
    '((tmax + tmin) / 2) * 0.1',
    {
      tmax: img.select('tmmx'),
      tmin: img.select('tmmn')
    }
  ).rename('tmean').toFloat();

  var pr = img.select('pr').rename('pr').toFloat();

  var pet = img.select('pet')
    .multiply(0.1)
    .rename('pet')
    .toFloat();

  var aet = img.select('aet')
    .multiply(0.1)
    .rename('aet')
    .toFloat();

  var def = img.select('def')
    .multiply(0.1)
    .rename('def')
    .toFloat();

  var soil = img.select('soil')
    .multiply(0.1)
    .rename('soil')
    .toFloat();

  var vpd = img.select('vpd')
    .multiply(0.01)
    .rename('vpd')
    .toFloat();

  var out = ee.Image.cat([
    tmean, pr, pet, aet, def, soil, vpd
  ]).toFloat();

  out = ee.Image(out.copyProperties(img, ['system:time_start', 'system:index']));
  return out;
}

var tc = tcRaw.map(prepClimate);

// =====================================================
// 4. 各年份统计
// 每年计算一组：
// - 生长季平均温度
// - 生长季降水总量
// - 生长季 PET 总量
// - 生长季 AET 总量
// - 生长季水分亏缺 DEF 总量
// - 生长季平均 VPD
// - 生长季平均土壤水分
// - 生长季 P/PET 比
// - 全年月均温度季节性（std）
// - 全年月降水季节性（cv）
// =====================================================
function yearlyClimateStats(year) {
  var y = Number(year);

  var yearlyAll = tc.filter(ee.Filter.calendarRange(y, y, 'year'));
  var yearlyGrow = yearlyAll.filter(
    ee.Filter.calendarRange(GROW_START_MONTH, GROW_END_MONTH, 'month')
  );

  var tmeanGs = yearlyGrow.select('tmean')
    .mean()
    .rename('tmean_gs_mean');

  var prGs = yearlyGrow.select('pr')
    .sum()
    .rename('pr_gs_sum');

  var petGs = yearlyGrow.select('pet')
    .sum()
    .rename('pet_gs_sum');

  var aetGs = yearlyGrow.select('aet')
    .sum()
    .rename('aet_gs_sum');

  var defGs = yearlyGrow.select('def')
    .sum()
    .rename('def_gs_sum');

  var vpdGs = yearlyGrow.select('vpd')
    .mean()
    .rename('vpd_gs_mean');

  var soilGs = yearlyGrow.select('soil')
    .mean()
    .rename('soil_gs_mean');

  var prPetRatioGs = prGs.divide(petGs.add(1e-6))
    .rename('pr_pet_gs_ratio');

  // 全年月均温季节性：月尺度 std
  var tmeanMonthlyStd = yearlyAll.select('tmean')
    .reduce(ee.Reducer.stdDev())
    .rename('tmean_monthly_std');

  // 全年月降水季节性：月尺度变异系数 CV = std / mean
  var prMonthlyMean = yearlyAll.select('pr')
    .mean()
    .rename('pr_monthly_mean_tmp');

  var prMonthlyStd = yearlyAll.select('pr')
    .reduce(ee.Reducer.stdDev())
    .rename('pr_monthly_std_tmp');

  var prMonthlyCv = prMonthlyStd.divide(prMonthlyMean.add(1e-6))
    .rename('pr_monthly_cv');

  var yearlyOut = ee.Image.cat([
    tmeanGs,
    prGs,
    petGs,
    aetGs,
    defGs,
    vpdGs,
    soilGs,
    prPetRatioGs,
    tmeanMonthlyStd,
    prMonthlyCv
  ]).toFloat();

  yearlyOut = ee.Image(yearlyOut.set('year', y));
  return yearlyOut;
}

var yearlyStatsCol = ee.ImageCollection.fromImages(
  YEARS.map(yearlyClimateStats)
);

print('Yearly climate stats image count:', yearlyStatsCol.size());

// =====================================================
// 5. 多年均值/波动（2019-2021）
// =====================================================
function multiYearMeanStd(inBand, outPrefix) {
  var meanImg = yearlyStatsCol.select(inBand)
    .mean()
    .rename(outPrefix + '_mean');

  var stdImg = yearlyStatsCol.select(inBand)
    .reduce(ee.Reducer.stdDev())
    .rename(outPrefix + '_std');

  return ee.Image.cat([meanImg, stdImg]).toFloat();
}

var output = ee.Image.cat([
  multiYearMeanStd('tmean_gs_mean',    'clim_tmean_gs'),
  multiYearMeanStd('pr_gs_sum',        'clim_pr_gs'),
  multiYearMeanStd('pet_gs_sum',       'clim_pet_gs'),
  multiYearMeanStd('aet_gs_sum',       'clim_aet_gs'),
  multiYearMeanStd('def_gs_sum',       'clim_def_gs'),
  multiYearMeanStd('vpd_gs_mean',      'clim_vpd_gs'),
  multiYearMeanStd('soil_gs_mean',     'clim_soil_gs'),
  multiYearMeanStd('pr_pet_gs_ratio',  'clim_prpet_gs_ratio'),
  multiYearMeanStd('tmean_monthly_std','clim_tmean_seasonality'),
  multiYearMeanStd('pr_monthly_cv',    'clim_pr_seasonality')
]).toFloat().clip(aoi);

// =====================================================
// 6. 对齐为 250 m 导出
// 说明：TerraClimate 原始较粗，这里用 bilinear 只是做网格对齐
// =====================================================
var output250 = output
  .resample('bilinear')
  .toFloat();

print('Final band names:', output250.bandNames());
print('Final band count:', output250.bandNames().size());
print('Output nominal scale (source):', output250.projection().nominalScale());

// =====================================================
// 7. 导出
// =====================================================
Export.image.toDrive({
  image: output250,
  description: 'Chengdu_SeasonalClimateWater_Vars_2019_2021_250m',
  folder: 'GEE_Exports',
  fileNamePrefix: 'Chengdu_SeasonalClimateWater_Vars_2019_2021_250m',
  region: aoi,
  scale: EXPORT_SCALE,
  crs: EXPORT_CRS,
  maxPixels: 1e13
});