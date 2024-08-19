var oliCol = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
    etmCol = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2"),
    bounds = 
    /* color: #ffc82d */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-104.18130547118324, 49.418064432479945],
          [-104.18130547118324, 36.92347853342345],
          [-86.78994804930824, 36.92347853342345],
          [-86.78994804930824, 49.418064432479945]]], null, false),
    test_bounds = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-93.17898273137061, 41.462743236916616],
          [-93.17898273137061, 41.44910588831326],
          [-93.14482211735694, 41.44910588831326],
          [-93.14482211735694, 41.462743236916616]]], null, false);

/*
 *Tillage Analysis
 */

// Setup ---------------------------------------------------------------#

// User defined -----------#
var ndviT = ee.Number(0.129);
var ndwiT = ee.Number(-0.081);
var min_images = ee.Number(2);
//-------------------------#

var years = ee.List.sequence(2016, 2021, 1);
var startDate = ee.Date.fromYMD(2017, 1, 1);
var endDate = ee.Date.fromYMD(2021, 12, 31);

// Harmonize ---------------------------------------------------------------#
var coefficients = {
  itcps: ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172])
             .multiply(10000),
  slopes: ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071])
};

// Function to get and rename bands of interest from OLI.
function renameOli(img) {
  return img.select(
      ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL'],
      ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa']);
}

// Function to get and rename bands of interest from ETM+.
function renameEtm(img) {
  return img.select(
      ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
      ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa']);
}

// From https://github.com/google/earthengine-community/blob/master/tutorials/landsat-etm-to-oli-harmonization/index.md
function etmToOli(img) {
  return img.select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
      .multiply(coefficients.slopes)
      .add(coefficients.itcps)
      .round()
      .toShort()
      .addBands(img.select('pixel_qa'));
}

// Cloud Mask ------------------------------------------------------------------#

function fmask(img) {
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var qa = img.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask)
                 .eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return img.updateMask(mask);
}

// Count Mask -----------------------------------------------------------------#
function countMask(img) {
  var mask = img.select('NDTI_MASK').count().gt(2);
  return img.updateMask(mask);
}

// Band Indexes ----------------------------------------------------------------#

function calcIndices(img) {
  var ndvi = img.normalizedDifference(['NIR', 'Red']).rename('NDVI');
  var ndti =  img.normalizedDifference(['SWIR1', 'SWIR2']).rename('NDTI');
  var ndwi = img.normalizedDifference(['NIR', 'SWIR1']).rename('NDWI');
  return ee.Image([ndvi,ndti,ndwi]);
}

function addNDTIinverse(img) {
  var NDTIinv = img.select('NDTI_MASK').divide(-1).rename('NDTI_INV');
  return img.addBands(NDTIinv);
}

// Wrapper Functions ----------------------------------------------------------------#

// Define function to prepare OLI images.
function prepOli(img) {
  var orig = img;
  img = renameOli(img);
  img = fmask(img);
  img = calcIndices(img);
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}

// Define function to prepare ETM+ images.
function prepEtm(img) {
  var orig = img;
  img = renameEtm(img);
  img = fmask(img);
  img = etmToOli(img);
  img = calcIndices(img);
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}

// Add quality Masks ---------------------------------------------------------------#
function addQualMask(image) {
  var ndti = image.select('NDTI');
  var ndvi = image.select('NDVI');
  var ndwi = image.select('NDWI');
  var QualMask = ndvi.lt(ndviT).and(ndwi.lt(ndwiT)).rename('QUALMASK');

  var masked = ee.Image([
                        ndti.mask(QualMask).rename('NDTI_MASK'),
                        ndvi.mask(QualMask).rename('NDVI_MASK'),
                        ndwi.mask(QualMask).rename('NDWI_MASK'),
  ]);
  
  return image.addBands(masked);
}


// Create model by year

function regression_linear(image, covar, slope, intercept){
    var reg = image.expression(
        '(a * covar) + b', {
        'covar': image.select(covar),
        'a': slope,
        'b': intercept,
    });

  var rnd = reg.round().int8();
  var out = rnd.where(rnd.lt(0), 0).where(rnd.gt(100), 100);
  return out;
}

function model(year) {
  var year = ee.Number(year);
  var yearText = year.format('%.0f');
  var year1back = ee.List([ee.Date.fromYMD(year.subtract(1), 1, 1),ee.Date.fromYMD(year.subtract(1), 12, 31)]);
  
  var cdlCornSoy = ee.ImageCollection('USDA/NASS/CDL')
    .filter(ee.Filter.and(ee.Filter.bounds(bounds),
            ee.Filter.date(year1back.get(0), year1back.get(1))))
    .first()
    .remap([1, 12, 13, 225, 226, 228, 237, 241, 5, 26, 239, 240, 254], [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0], NaN, 'cropland')//0:soybean 1:corn
    .toInt8()
    .rename('cdl');

  var colYear = colQual.filter(ee.Filter.date(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year.add(1), 12, 31)));
  var count = colYear.select('NDTI_MASK').count().gte(min_images);
  var minNDTI = colYear.min().updateMask(count);
  
  //Regression for corn fields
  var corn = regression_linear(minNDTI, 'NDTI_MASK', 676, 20.5).rename('corn');
         
  //Regression for soybean fields
  var bean = regression_linear(minNDTI, 'NDTI_MASK', 676, 20.5).rename('soybean');
              
  //Combine Corn and Soybean rasters then clip to geometry
  var combined = ee.Image(-100).where(cdlCornSoy.eq(1), corn).where(cdlCornSoy.eq(0), bean).rename('residue_cover');
  
  var combined = combined.mask(combined.gt(-100)).setMulti({'system:index': year.format('%.0f')});

  return combined;
}

// ACTUAL MODEL RUN -------------------------------------------------------------#
var colFilter = ee.Filter.and(
    ee.Filter.bounds(bounds),
    ee.Filter.calendarRange(60, 166, 'day_of_year'),
    ee.Filter.date(startDate, endDate),
    ee.Filter.lt('CLOUD_COVER', 50), ee.Filter.lt('GEOMETRIC_RMSE_MODEL', 10),
    ee.Filter.or(ee.Filter.eq('IMAGE_QUALITY', 9), ee.Filter.eq('IMAGE_QUALITY_OLI', 9)));

// Filter collections and prepare them for merging.
oliCol = oliCol.filter(colFilter).map(prepOli);
etmCol = etmCol.filter(colFilter).map(prepEtm);

// Merge the collections.
var col = oliCol.merge(etmCol);

// Filter using NDVI and NDWI
var colQual = col.map(addQualMask).select(['NDTI_MASK']);

// Apply model
var residueCoverCol = ee.ImageCollection(years.map(model));
var residueCoverImg = residueCoverCol.toBands();

// View and export results -------------------------------------------------------------#
print(residueCoverImg);


Export.image.toDrive({
  image: residueCoverImg.select('2017_residue_cover').unmask(-100),
  maxPixels: 1e13,
  description: 'rescover_2017',
  scale: 30,
  fileFormat: 'GeoTIFF',
  formatOptions: {noData: -100},
  region: bounds,
  folder:'residue_cover_forpaper'
});

var resViz = {min: 0, max: 100, palette: ['red', 'yellow', 'green']};
Map.addLayer(residueCoverCol.median(), resViz, 'Residue Cover', false);
