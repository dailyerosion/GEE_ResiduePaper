/*
 *Tillage Analysis
 */

var fields_all = ee.FeatureCollection("survey_data_in30");

// Change these values to fit the run
var year = 2016;
var fields = fields_all.filter(ee.Filter.eq('year','2016'));
var outName = 'all'+year+'_30';
var outFolder = 'tillage_analysis';
var scale = 30;

var back1 = year-1;
var back2 = year-2;
var startDate =  year+'-03-01';//changed in Brain meeting -03-01
var endDate = year+'-06-15';//changed in Brain meeting -06-15

/*
Data Validity Periods:
Sentinel2 - July 2015 to present      2016 -
Sentinel1 - October 2014 to present   2015 -
SMAP - April 2015 to present          2015 - 
Landsat8 - 2013 to present            2013 -
Landsat7 - 2000 to present            2000 -
*/

//Mask~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
function addQualMask(image) {
  var ndti = image.select('NDTI');
  var ndvi =image.select('NDVI');
  var swc = image.select('SWC');
  var QualMask = ndvi.lt(0.3).and(swc.lt(1)).rename('QUALMASK');

  var masked = ee.Image([
                        ndti.mask(QualMask).rename('NDTI_MASK'),
                        ndvi.mask(QualMask).rename('NDVI_MASK'),
                        swc.mask(QualMask).rename('SWC_MASK'),
  ]);
  
  return image.addBands(masked);
}


//SMAP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
var smap = ee.ImageCollection("NASA_USDA/HSL/SMAP10KM_soil_moisture")
  .filterDate(startDate, endDate)
  .filterBounds(fields)
  .select('ssm');
  
var smapMin = smap.min();
var smapMedian = smap.median();
var smapCount = smap.count();

print('smap:',smap.size());
//Sentinel 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Load the Sentinel-1 ImageCollection
function S1addD(image) {
  // var d = image.normalizedDifference(['VV','VH']).rename('NDI');
  var ndi = image.select('VV').subtract(image.select('VH')).divide(
    image.select('VV').add(image.select('VH'))).rename('NDI')
  return image.addBands(ndi);
}

var sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .map(S1addD)
  .filterDate(startDate, endDate)
  .filterBounds(fields)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .select(['VV', 'VH', 'NDI']);

var s1Min = sentinel1.min();
var s1Median = sentinel1.median();
var s1Count = sentinel1.count();

print('sentinel1:',sentinel1.size());
//Sentinel 2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Function to mask clouds using the Sentinel-2 QA band.
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));

  // Return the masked and scaled data, without the QA bands.
  return image.updateMask(mask).divide(10000)
      .select("B.*")
      .copyProperties(image, ["system:time_start"]);
}
function s2addNDTI(image) {
  var ndti = image.normalizedDifference(['B11', 'B12']).rename('NDTI');
  return image.addBands(ndti);
}
function s2addNDVI(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  return image.addBands(ndvi);
}
function s2addSWC(image) {
  var nir = image.select('B8');
  var swir2 = image.select('B12');
  var swc = nir.divide(swir2).rename('SWC');
  return image.addBands(swc);
}

var sentinel2 = ee.ImageCollection('COPERNICUS/S2')
                  .filterDate(startDate, endDate)
                  .filterBounds(fields)
                  .map(maskS2clouds)
                  .map(s2addNDVI)
                  .map(s2addNDTI)
                  .map(s2addSWC)
                  .select(['NDVI', 'NDTI', 'SWC']);

var sentinel2Thresh = sentinel2.map(addQualMask).select(['NDVI_MASK', 'NDTI_MASK', 'SWC_MASK']);

var s2Min = sentinel2Thresh.min();
var s2Median = sentinel2Thresh.median();
var s2Count = sentinel2Thresh.count();
var s2CountFull = sentinel2.count();

print('sentinel2:',sentinel2Thresh.size());
//Landsat 8~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}
function l8addNDVI(image) {
  var ndvi = image.normalizedDifference(['B5', 'B4']).rename('NDVI');
  return image.addBands(ndvi);
}
function l8addNDTI(image) {
  var ndti = image.normalizedDifference(['B6', 'B7']).rename('NDTI');
  return image.addBands(ndti);
}
function l8addSWC(image) {
  var b5 = image.select('B5');
  var b7 = image.select('B7');
  var swc = b5.divide(b7).rename('SWC');
  return image.addBands(swc);
}

var landsat8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
                  .filterDate(startDate, endDate)
                  .filterBounds(fields)
                  .map(maskL8sr)
                  .map(l8addNDVI)
                  .map(l8addNDTI)
                  .map(l8addSWC)
                  .select(['NDVI', 'NDTI', 'SWC']);

var landsat8Thresh = landsat8.map(addQualMask).select(['NDVI_MASK', 'NDTI_MASK', 'SWC_MASK']);

var L8Min = landsat8Thresh.min();
var L8Median = landsat8Thresh.median();
var L8Count = landsat8Thresh.count();
var L8CountFull = landsat8.count();

print('landsat8:',landsat8Thresh.size());
//Landsat 7~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
function l7addNDVI(image) {
  var ndvi = image.normalizedDifference(['B4', 'B3']).rename('NDVI');
  return image.addBands(ndvi);
}
function l7addNDTI(image) {
  var ndti = image.normalizedDifference(['B5', 'B7']).rename('NDTI');
  return image.addBands(ndti);
}
function l7addSWC(image) {
  var b4 = image.select('B4');
  var b7 = image.select('B7');
  var swc = b4.divide(b7).rename('SWC');
  return image.addBands(swc);
}
var cloudMaskL457 = function(image) {
  var qa = image.select('pixel_qa');
  // If the cloud bit (5) is set and the cloud confidence (7) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var cloud = qa.bitwiseAnd(1 << 5)
                  .and(qa.bitwiseAnd(1 << 7))
                  .or(qa.bitwiseAnd(1 << 3));
  // Remove edge pixels that don't occur in all bands
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2);
};

var landsat7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
                  .filterDate(startDate, endDate)
                  .filterBounds(fields)
                  .map(cloudMaskL457)
                  .map(l7addNDVI)
                  .map(l7addNDTI)
                  .map(l7addSWC)
                  .select(['NDVI', 'NDTI', 'SWC']);

var landsat7Thresh = landsat7.map(addQualMask).select(['NDVI_MASK', 'NDTI_MASK', 'SWC_MASK']);

var L7Min = landsat7Thresh.min();
var L7Median = landsat7Thresh.median();
var L7Count = landsat7Thresh.count();
var L7CountFull = landsat7.count();

print('landsat7:',landsat7Thresh.size());
//Combine Landsats~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
var landsatc = landsat7.merge(landsat8);
var landsatcThresh = landsat7Thresh.merge(landsat8Thresh);

var LCMin = landsatcThresh.min();
var LCMedian = landsatcThresh.median();
var LCCount = landsatcThresh.count();
var LCCountFull = landsatc.count();

print('landsatC:',landsatcThresh.size());
//Crop Land Data Layer~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
var CDLback1 = ee.ImageCollection('USDA/NASS/CDL')
                  .filter(ee.Filter.date(back1+'-01-01', back1+'-12-31'))
                  .median();
var CDLback2 = ee.ImageCollection('USDA/NASS/CDL')
                  .filter(ee.Filter.date(back2+'-01-01', back2+'-12-31'))
                  .median();

var CDL = ee.Image([
                    CDLback1.select('cropland').rename('cdlBack1'),
                    CDLback2.select('cropland').rename('cdlBack2')
                    ]);

//Combine All~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
var all = ee.Image([
                    s2Min.select('NDVI_MASK').rename('S2NDVImin'),
                    s2Min.select('NDTI_MASK').rename('S2NDTImin'),
                    s2Min.select('SWC_MASK').rename('S2SWCmin'),
                    s2Median.select('NDVI_MASK').rename('S2NDVImedian'),
                    s2Median.select('NDTI_MASK').rename('S2NDTImedian'),
                    s2Median.select('SWC_MASK').rename('S2SWCmedian'),
                    s2Count.select('NDTI_MASK').rename('S2count'),
                    s2CountFull.select('NDTI').rename('S2countfull'),
                    s1Min.select('VV').rename('S1VVmin'),
                    s1Min.select('VH').rename('S1VHmin'),
                    s1Min.select('NDI').rename('S1NDImin'),
                    s1Median.select('VV').rename('S1VVmedian'),
                    s1Median.select('VH').rename('S1VHmedian'),
                    s1Median.select('NDI').rename('S1NDImedian'),
                    s1Count.select('VV').rename('S1count'),
                    smapMin.select('ssm').rename('SMAPmin'),
                    smapMedian.select('ssm').rename('SMAPmedian'),
                    smapCount.select('ssm').rename('SMAPcount'),
                    L8Min.select('NDVI_MASK').rename('L8NDVImin'),
                    L8Min.select('NDTI_MASK').rename('L8NDTImin'),
                    L8Min.select('SWC_MASK').rename('L8SWCmin'),
                    L8Median.select('NDVI_MASK').rename('L8NDVImedian'),
                    L8Median.select('NDTI_MASK').rename('L8NDTImedian'),
                    L8Median.select('SWC_MASK').rename('L8SWCmedian'),
                    L8Count.select('NDTI_MASK').rename('L8count'),
                    L8CountFull.select('NDTI').rename('L8countfull'),
                    L7Min.select('NDVI_MASK').rename('L7NDVImin'),
                    L7Min.select('NDTI_MASK').rename('L7NDTImin'),
                    L7Min.select('SWC_MASK').rename('L7SWCmin'),
                    L7Median.select('NDVI_MASK').rename('L7NDVImedian'),
                    L7Median.select('NDTI_MASK').rename('L7NDTImedian'),
                    L7Median.select('SWC_MASK').rename('L7SWCmedian'),
                    L7Count.select('NDTI_MASK').rename('L7count'),
                    L7CountFull.select('NDTI').rename('L7countfull'),
                    LCMin.select('NDVI_MASK').rename('LCNDVImin'),
                    LCMin.select('NDTI_MASK').rename('LCNDTImin'),
                    LCMin.select('SWC_MASK').rename('LCSWCmin'),
                    LCMedian.select('NDVI_MASK').rename('LCNDVImedian'),
                    LCMedian.select('NDTI_MASK').rename('LCNDTImedian'),
                    LCMedian.select('SWC_MASK').rename('LCSWCmedian'),
                    LCCount.select('NDTI_MASK').rename('LCcount'),
                    LCCountFull.select('NDTI').rename('LCcountfull')
                    ]);

var cdlZonal = CDL.reduceRegions({
    collection: fields, 
    reducer: ee.Reducer.mode(), 
    crs: 'EPSG:26915', //NAD 83 UTM Zone 15 N
    scale: scale
  });
var allZonal = all.reduceRegions({
    collection: cdlZonal, 
    reducer: ee.Reducer.median(), 
    crs: 'EPSG:26915', //NAD 83 UTM Zone 15 N
    scale: scale
  });

Export.table.toDrive({
  collection: allZonal,
  description:outName,
  fileFormat: 'csv',
  folder:outFolder
});
