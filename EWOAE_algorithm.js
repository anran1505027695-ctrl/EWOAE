var table = ee.FeatureCollection("projects/ee-fengfeifei0729/assets/globalreservoir/SixReservoir2025626")
var idList = [2368];
var Reservoir = ee.FeatureCollection("projects/ee-waterresearch/assets/GRanD_reservoirs_v1_3").filter(ee.Filter.inList('GRAND_ID', idList)).first();
var jrc_mth = ee.ImageCollection("JRC/GSW1_4/MonthlyHistory");
var SRTMGL1dataset= ee.Image('USGS/SRTMGL1_003');
var SRTMGL1 = SRTMGL1dataset.select('elevation')
var AW3D30dataset = ee.ImageCollection('JAXA/ALOS/AW3D30/V3_2');
var elevation = AW3D30dataset.select('DSM');
var proj = elevation.first().select(0).projection();
var AW3D30 = elevation.mosaic().setDefaultProjection(proj)
var Mask = AW3D30.updateMask(SRTMGL1.unmask(-9999).eq(-9999))
var dem = SRTMGL1.unmask().blend(Mask);
var prj = jrc_mth.first().projection()
var clearID = table.filter(ee.Filter.eq("ContamClass",1)).filter(ee.Filter.inList('GRAND_ID', idList))
// print("clearID",clearID)
var ContaID = table.filter(ee.Filter.eq("ContamClass",3)).filter(ee.Filter.inList('GRAND_ID', idList))
print("ContaID", ContaID)
// print("ContaID",clearID)
var occr = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select('occurrence').selfMask();


function ImageColl(feature) {
  var year = feature.get("Jrc_Year")
  
  // var lake = lakefeature.geometry() 
  var month = feature.get("Month")
  var wd = feature.get("wd")
  var nodata = feature.get("nodata")//ee.Number(feature.get("nodata")).multiply(100).round()
  var image = jrc_mth.filter(ee.Filter.and(
  ee.Filter.eq('year', year),
  ee.Filter.eq('month', month)
  )).first()//.clip(lake);
  return image.set({"GRAND_ID":feature.get("GRAND_ID"),"wd":wd,"ClearImageYear":year,"ClearImageMonth":month,"CP":nodata})
}
function ImageColl2(feature) {
  var year = feature.get("Jrc_Year")
  
  // var lake = lakefeature.geometry() 
  var month = feature.get("Month")
  var wd = feature.get("wd")
  var nodata = feature.get("nodata")
  var image = jrc_mth.filter(ee.Filter.and(
  ee.Filter.eq('year', year),
  ee.Filter.eq('month', month)
  )).first()//.clip(lake);
  return image.set({"GRAND_ID":feature.get("GRAND_ID"),"ContanImageYear":year,"ContanImageMonth":month,"CP":nodata})
}

var visualization = {
  min: 0.0,
  max: 2.0,
  palette: ['green', 'fffcb8', '0905ff']
};
var ContamImageColl = ContaID.map(ImageColl2).toList(184)
var Conta_Image = ee.Image(ContamImageColl.get(136))


function EWOAE(pollutionImage) {
  
    /// input is an image with single band ['water']
    /// values 0, 1, 2 represents no_data, not_water, water
    /// occr image is the water occurrence image with band ['occurrence']
    /// Sobel operator
  var kx = ee.Kernel.fixed(3, 3, [[-1,0,1],[-2,0,2],[-1,0,1]]);
  var ky = kx.rotate(-1);
  var image_rmp = pollutionImage.remap([0,1,2],[3,0,6]);
  var gdt = image_rmp.convolve(kx).addBands(image_rmp.convolve(ky))
                .expression("sqrt(b(0)*b(0)+b(1)*b(1))");

  var water_edge1 = gdt.gte(14).updateMask(pollutionImage.eq(2)).reproject(prj).rename("Coastline")
  
  var DEM1 = dem.updateMask(water_edge1)
  var hi = DEM1.reduceRegions({
    collection: ee.FeatureCollection([Reservoir]),
		reducer: ee.Reducer.autoHistogram(),
		scale:prj.nominalScale(),
    crs:prj.crs(),
    // tileScale:16
  }).first().get('histogram')
  hi = ee.Array(ee.Algorithms.If(hi, hi, [[-9999, -9999]]));
    // 累加像元个数并计算CDF
    // 提取 occurrence 和 count 列
  var demvalue = hi.slice(1, 0, 1);   // DEM 值
  var cnt = hi.slice(1, 1, 2);   // 对应像元个数
  var cumCnt = cnt.accum(0);                     // 累积和
  var totalCnt = cnt.reduce('sum', [0]);         // 总像元数
  var totalCntScalar = ee.Number(totalCnt.get([0,0])); // 取出标量
  var cdf = cumCnt.divide(totalCntScalar);


    // 找到CDF ≥ 0.9的第一个 occurrence 值
  var mask = cdf.gte(ee.Number(90).divide(100));
  // var mask = cdf.gte(ee.Number(CDFNumber).divide(100));
  var DemThd = demvalue.mask(mask).get([0,0]);
  var maskdem = DEM1.updateMask(DEM1.gte(DemThd)).rename("maskdem").clip(Reservoir)
  var OCC = occr.clip(Reservoir).updateMask(DEM1.gte(DemThd))//.clip(Reservoir)
  var occThd = ee.Number(OCC.reduceRegions({
    collection: ee.FeatureCollection([Reservoir]),
		reducer: ee.Reducer.mean(),
		scale:prj.nominalScale(),
    crs:prj.crs(),
  }).first().get("mean"))
  
  occThd = ee.Number(ee.Algorithms.If(occThd,occThd.round(),-9999))
  var Recover = ee.Image(2).updateMask(pollutionImage.eq(0).and(occr.updateMask(pollutionImage.eq(0)).gte(ee.Number(occThd))))
  var EnhancedImage = pollutionImage.blend(Recover).rename("Enhanced")
  return  EnhancedImage.clip(Reservoir).addBands(pollutionImage.clip(Reservoir))

}
var AllImage = EWOAE(Conta_Image)


Map.centerObject(Reservoir,10)
Map.addLayer(AllImage.select("water"),visualization,"Raw water area", false)
Map.addLayer(AllImage.select("Enhanced"),visualization,"Enhanced water area", false)

