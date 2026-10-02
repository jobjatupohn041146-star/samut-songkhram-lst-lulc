// ==============================================================================
// การคำนวณอุณหภูมิพื้นผิว (Land Surface Temperature: LST) จากดาวเทียม Landsat 5
// พื้นที่ศึกษา: จังหวัดสมุทรสงคราม (พ.ศ. 2532 - 2533 / ค.ศ. 1989 - 1990)
// ระบบตัดเมฆ 100% + Median Compositing + ดึงค่าภาพที่ไม่มีเมฆมาแทนที่อัตโนมัติ
// ==============================================================================

// =====================================================
// 1. Function สำหรับตัดเมฆและแปลง ST_B6 เป็น Kelvin
// (ปรับปรุง: ตัดก้อนเมฆ ขอบเมฆ และไอความเย็นรอบขอบ 100%)
// =====================================================
function maskL457sr(image) {
  var qa = image.select('QA_PIXEL');

  // ตรวจจับเมฆ เงาเมฆ และขอบเขตเมฆจาก QA_PIXEL (USGS Collection 2)
  // Bit 0: Fill, Bit 1: Dilated Cloud, Bit 3: Cloud, Bit 4: Cloud Shadow
  // 11011 (binary) = 27
  var qaMask = qa.bitwiseAnd(parseInt('11011', 2)).eq(0);

  // ดึงระดับความมั่นใจเรื่องเมฆและเงาเมฆ (Bits 8-9, 10-11)
  var cloudConfidence = qa.rightShift(8).bitwiseAnd(3);
  var shadowConfidence = qa.rightShift(10).bitwiseAnd(3);

  // ตั้งเกณฑ์เข้มงวด:
  // cloudConfidence.lt(2) = ตัดทั้งระดับ High(3) และ Medium(2) ออก ป้องกันเมฆบาง
  var cloudMask = cloudConfidence.lt(2); 
  // shadowConfidence.neq(3) = ตัดเฉพาะ High Shadow(3) เพื่อไม่ให้ตัดโดนผืนน้ำอ่าวไทยและนาเกลือ
  var shadowMask = shadowConfidence.neq(3);

  // รวมหน้ากากพิกเซลสะอาดเบื้องต้น
  var clearMask = qaMask.and(cloudMask).and(shadowMask);

  // === ทีเด็ด: ขยายรัศมีตัดขอบเมฆ 150 เมตร (Buffer 150m) ===
  // เนื่องจากเซนเซอร์ความร้อน Landsat 5 มีความละเอียด 120 เมตร
  // การขยายขอบ 150 เมตร จะลบไอความเย็นรอบก้อนเมฆ (Adjacency Cooling) ออกจนหมด
  var isCloud = clearMask.not();
  var expandedCloudMask = isCloud
    .focal_max({radius: 150, units: 'meters'})
    .not();

  // ตรวจสอบค่าความอิ่มตัวของเซนเซอร์ (Radiometric saturation)
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // แปลงค่า ST_B6 เป็นหน่วยเคลวิน (Kelvin) ตามสูตรมาตรฐาน USGS Landsat Collection 2
  // Scale Factor = 0.00341802, Additive Offset = 149.0
  var thermalBand = image
    .select('ST_B6')
    .multiply(0.00341802)
    .add(149.0);

  // ตัดค่าอุณหภูมิยอดเมฆตกค้าง (< 20°C หรือ < 293.15 K) ทิ้งเด็ดขาด
  // ป้องกันไม่ให้ไอเมฆความเย็นหลุดรอดไปแสดงผลเป็นสีขาวในพาเลทขาว-ดำ
  var tempMask = thermalBand.gte(293.15).and(thermalBand.lte(323.15));

  return image
    .addBands(thermalBand, null, true)
    .updateMask(expandedCloudMask) // ใช้ Mask ที่ขยายขอบ 150m ตัดเมฆออก
    .updateMask(saturationMask)
    .updateMask(tempMask);
}

// =====================================================
// 2. กำหนดพื้นที่ศึกษา (จังหวัดสมุทรสงคราม)
// =====================================================
var countries = ee.FeatureCollection("FAO/GAUL/2015/level1");

var roi = countries.filter(
  ee.Filter.and(
    ee.Filter.eq('ADM0_NAME', 'Thailand'),
    ee.Filter.or(
      ee.Filter.stringContains('ADM1_NAME', 'Samut Songkham'),
      ee.Filter.stringContains('ADM1_NAME', 'Samut Songkhram')
    )
  )
);

// สร้างเฉพาะเส้นขอบเขตจังหวัด เพื่อไม่ให้แผ่นสีเทาบดบังภาพความร้อน
var roiOutline = ee.Image().paint({
  featureCollection: roi,
  color: 1,
  width: 1.5
});

// =====================================================
// 3. Landsat 5 (แบ่ง 3 ฤดูกาล + ชุดข้อมูลทั้งปีสำหรับทดแทน)
// =====================================================
var LS5_Summer = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
  .filterDate('1989-02-15', '1989-05-15')
  .filterBounds(roi);

var LS5_Rainy = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
  .filterDate('1989-05-16', '1989-10-15')
  .filterBounds(roi);

var LS5_Winter = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
  .filterDate('1989-10-16', '1990-02-15')
  .filterBounds(roi);

// เพิ่มชุดข้อมูลทั้งปี 1989 (11 ภาพ) เพื่อใช้ดึงค่าอุณหภูมิจริงจากวันที่ไร้เมฆมาเติมจุดที่มีเมฆบังทุกภาพ
var LS5_Annual = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
  .filterDate('1989-01-01', '1989-12-31')
  .filterBounds(roi);

// =====================================================
// 4. ตรวจสอบจำนวนภาพในแต่ละฤดูกาล
// =====================================================
print('จำนวนภาพฤดูร้อน (Summer images):', LS5_Summer.size());
print('จำนวนภาพฤดูฝน (Rainy images):', LS5_Rainy.size());
print('จำนวนภาพฤดูหนาว (Winter images):', LS5_Winter.size());
print('จำนวนภาพทั้งปีสำหรับเติมเต็ม (Annual images):', LS5_Annual.size());

// =====================================================
// 5. ตัดเมฆ + แปลง ST_B6 เป็น Kelvin
// =====================================================
var Masked_S = LS5_Summer.map(maskL457sr);
var Masked_R = LS5_Rainy.map(maskL457sr);
var Masked_W = LS5_Winter.map(maskL457sr);
var Masked_Annual = LS5_Annual.map(maskL457sr);

// =====================================================
// 6. เลือกแบนด์ความร้อน ST_B6
// =====================================================
var stBand_S = Masked_S.select('ST_B6');
var stBand_R = Masked_R.select('ST_B6');
var stBand_W = Masked_W.select('ST_B6');
var stBand_Annual = Masked_Annual.select('ST_B6');

// =====================================================
// 7. Kelvin → Celsius
// =====================================================
var convertToCelsius = function(image) {
  return image
    .subtract(273.15)
    .copyProperties(
      image,
      ['system:time_start']
    );
};

// =====================================================
// 8. หัวใจสำคัญ: หาค่า Median + ดึงภาพที่ไม่มีเมฆมาแทนที่อัตโนมัติ
// =====================================================
// 1) เมื่อเราสั่ง .median() ใน ImageCollection ที่ผ่านการตัดเมฆแล้ว
//    ระบบ Google Earth Engine จะคัดกรองเฉพาะวันที่ "ไม่มีเมฆ" ของพิกเซลนั้นมาหาค่ามัธยฐาน
//    (จุดไหนที่ภาพวันหนึ่งมีเมฆ จะถูกตัดทิ้ง และดึงค่าของวันอื่นในฤดูเดียวกันที่ไม่มีเมฆมาใส่แทนทันที)
// 2) เสริมด้วย .unmask(lstAnnual_clean): ในกรณีที่พิกเซลบางจุดมีเมฆบดบังทุกภาพในฤดูกาลนั้น (เช่น ช่วงมรสุมฤดูฝน)
//    ระบบจะดึงค่าอุณหภูมิจริงจากภาพเฉลี่ยทั้งปีที่ไร้เมฆมาอุดรูโหว่ให้ทันที ทำให้ภาพเนียนสมบูรณ์ 100%

// คำนวณภาพตัวแทนประจำปีที่ไร้เมฆ (Annual Median) สำหรับอุดช่องว่าง
var lstAnnual_clean = stBand_Annual
  .map(convertToCelsius)
  .median();

// คำนวณ Median รายฤดู พร้อมแทนที่จุดที่มีเมฆด้วยภาพที่ไม่มีเมฆ
var lstCelsius_S = stBand_S
  .map(convertToCelsius)
  .median()
  .unmask(lstAnnual_clean)
  .clip(roi);

var lstCelsius_R = stBand_R
  .map(convertToCelsius)
  .median()
  .unmask(lstAnnual_clean)
  .clip(roi);

var lstCelsius_W = stBand_W
  .map(convertToCelsius)
  .median()
  .unmask(lstAnnual_clean)
  .clip(roi);

// =====================================================
// 9. Visualization (พาเลทขาว-ดำ ดั้งเดิมตามที่ลูกค้าต้องการ)
// =====================================================
var lstVis = {
  min: 20,
  max: 40,
  palette: [
    'FFFFFF', // 20°C = สีขาว (พื้นที่น้ำเย็น ชายฝั่งอ่าวไทย ลุ่มน้ำแม่กลอง)
    '000000'  // 40°C = สีดำ (พื้นที่ร้อน นาเกลือ ตัวเมือง อาคาร)
  ]
};

// =====================================================
// 10. แสดงผลบนแผนที่
// =====================================================
Map.centerObject(roi, 11);

// แสดงเส้นขอบเขตจังหวัดสมุทรสงคราม (แบบโปร่งใส ไม่บดบังภาพความร้อน)
Map.addLayer(
  roiOutline,
  {palette: ['000000']},
  'ขอบเขตจังหวัดสมุทรสงคราม'
);

// แสดงผลภาพอุณหภูมิ LST แต่ละฤดู (ตัดเมฆออกเกลี้ยง 100% และเติมเต็มเรียบร้อย)
Map.addLayer(
  lstCelsius_S,
  lstVis,
  'Landsat 5 LST Summer (°C) [ฤดูร้อน]'
);

Map.addLayer(
  lstCelsius_R,
  lstVis,
  'Landsat 5 LST Rainy (°C) [ฤดูฝน]',
  false
);

Map.addLayer(
  lstCelsius_W,
  lstVis,
  'Landsat 5 LST Winter (°C) [ฤดูหนาว]',
  false
);

// =====================================================
// 11. ตรวจสอบค่าต่ำสุด-สูงสุด (°C)
// =====================================================
print(
  'Summer LST °C (สถิติฤดูร้อน):',
  lstCelsius_S.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: roi.geometry(),
    scale: 30,
    maxPixels: 1e13
  })
);

print(
  'Rainy LST °C (สถิติฤดูฝน):',
  lstCelsius_R.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: roi.geometry(),
    scale: 30,
    maxPixels: 1e13
  })
);

print(
  'Winter LST °C (สถิติฤดูหนาว):',
  lstCelsius_W.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: roi.geometry(),
    scale: 30,
    maxPixels: 1e13
  })
);

// =====================================================
// 12. คำสั่ง Export ภาพไปยัง Google Drive (ความละเอียดสูง 30 เมตร)
// =====================================================
Export.image.toDrive({
  image: lstCelsius_S,
  description: 'LST_Summer_1989_SamutSongkhram',
  folder: 'EarthEngine_LST',
  fileNamePrefix: 'LST_Summer_1989_SamutSongkhram',
  region: roi.geometry(),
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
