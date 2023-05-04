! caviar_continuum.F90 - Compute continuum absorption using CAVIAR model

module caviar_continuum

public
  
contains
  
  subroutine calc_caviar_continuum(nlev, nwav, pressure, temperature, mole_fraction, &
       wavenumber_cm1, continuum)

    use parkind1,      only : jprb
    use interpolation, only : interpolate

    ! Avogadro's number
    real(jprb), parameter :: NAVOGADRO = 6.02214076e23_jprb
    
    ! CAVIAR self and foreign continua
    integer,    parameter :: nwav_lut = 2002
    real(jprb), parameter :: temperature_lut = 296.0_jprb
    real(jprb), parameter :: wavenumber_lut_cm1(nwav_lut) = [ &
      &  0_jprb, 10_jprb, 20_jprb, 30_jprb, 40_jprb, 50_jprb, 60_jprb, 70_jprb, &
      &  80_jprb, 90_jprb, 100_jprb, 110_jprb, 120_jprb, 130_jprb, 140_jprb, 150_jprb, &
      &  160_jprb, 170_jprb, 180_jprb, 190_jprb, 200_jprb, 210_jprb, 220_jprb, 230_jprb, &
      &  240_jprb, 250_jprb, 260_jprb, 270_jprb, 280_jprb, 290_jprb, 300_jprb, 310_jprb, &
      &  320_jprb, 330_jprb, 340_jprb, 350_jprb, 360_jprb, 370_jprb, 380_jprb, 390_jprb, &
      &  400_jprb, 410_jprb, 420_jprb, 430_jprb, 440_jprb, 450_jprb, 460_jprb, 470_jprb, &
      &  480_jprb, 490_jprb, 500_jprb, 510_jprb, 520_jprb, 530_jprb, 540_jprb, 550_jprb, &
      &  560_jprb, 570_jprb, 580_jprb, 590_jprb, 600_jprb, 610_jprb, 620_jprb, 630_jprb, &
      &  640_jprb, 650_jprb, 660_jprb, 670_jprb, 680_jprb, 690_jprb, 700_jprb, 710_jprb, &
      &  720_jprb, 730_jprb, 740_jprb, 750_jprb, 760_jprb, 770_jprb, 780_jprb, 790_jprb, &
      &  800_jprb, 810_jprb, 820_jprb, 830_jprb, 840_jprb, 850_jprb, 860_jprb, 870_jprb, &
      &  880_jprb, 890_jprb, 900_jprb, 910_jprb, 920_jprb, 930_jprb, 940_jprb, 950_jprb, &
      &  960_jprb, 970_jprb, 980_jprb, 990_jprb, 1000_jprb, 1010_jprb, 1020_jprb, 1030_jprb, &
      &  1040_jprb, 1050_jprb, 1060_jprb, 1070_jprb, 1080_jprb, 1090_jprb, 1100_jprb, 1110_jprb, &
      &  1120_jprb, 1130_jprb, 1140_jprb, 1150_jprb, 1160_jprb, 1170_jprb, 1180_jprb, 1190_jprb, &
      &  1200_jprb, 1210_jprb, 1220_jprb, 1230_jprb, 1240_jprb, 1250_jprb, 1260_jprb, 1270_jprb, &
      &  1280_jprb, 1290_jprb, 1300_jprb, 1310_jprb, 1320_jprb, 1330_jprb, 1340_jprb, 1350_jprb, &
      &  1360_jprb, 1370_jprb, 1380_jprb, 1390_jprb, 1400_jprb, 1410_jprb, 1420_jprb, 1430_jprb, &
      &  1440_jprb, 1450_jprb, 1460_jprb, 1470_jprb, 1480_jprb, 1490_jprb, 1500_jprb, 1510_jprb, &
      &  1520_jprb, 1530_jprb, 1540_jprb, 1550_jprb, 1560_jprb, 1570_jprb, 1580_jprb, 1590_jprb, &
      &  1600_jprb, 1610_jprb, 1620_jprb, 1630_jprb, 1640_jprb, 1650_jprb, 1660_jprb, 1670_jprb, &
      &  1680_jprb, 1690_jprb, 1700_jprb, 1710_jprb, 1720_jprb, 1730_jprb, 1740_jprb, 1750_jprb, &
      &  1760_jprb, 1770_jprb, 1780_jprb, 1790_jprb, 1800_jprb, 1810_jprb, 1820_jprb, 1830_jprb, &
      &  1840_jprb, 1850_jprb, 1860_jprb, 1870_jprb, 1880_jprb, 1890_jprb, 1900_jprb, 1910_jprb, &
      &  1920_jprb, 1930_jprb, 1940_jprb, 1950_jprb, 1960_jprb, 1970_jprb, 1980_jprb, 1990_jprb, &
      &  2000_jprb, 2010_jprb, 2020_jprb, 2030_jprb, 2040_jprb, 2050_jprb, 2060_jprb, 2070_jprb, &
      &  2080_jprb, 2090_jprb, 2100_jprb, 2110_jprb, 2120_jprb, 2130_jprb, 2140_jprb, 2150_jprb, &
      &  2160_jprb, 2170_jprb, 2180_jprb, 2190_jprb, 2200_jprb, 2210_jprb, 2220_jprb, 2230_jprb, &
      &  2240_jprb, 2250_jprb, 2260_jprb, 2270_jprb, 2280_jprb, 2290_jprb, 2300_jprb, 2310_jprb, &
      &  2320_jprb, 2330_jprb, 2340_jprb, 2350_jprb, 2360_jprb, 2370_jprb, 2380_jprb, 2390_jprb, &
      &  2400_jprb, 2410_jprb, 2420_jprb, 2430_jprb, 2440_jprb, 2450_jprb, 2460_jprb, 2470_jprb, &
      &  2480_jprb, 2490_jprb, 2500_jprb, 2510_jprb, 2520_jprb, 2530_jprb, 2540_jprb, 2550_jprb, &
      &  2560_jprb, 2570_jprb, 2580_jprb, 2590_jprb, 2600_jprb, 2610_jprb, 2620_jprb, 2630_jprb, &
      &  2640_jprb, 2650_jprb, 2660_jprb, 2670_jprb, 2680_jprb, 2690_jprb, 2700_jprb, 2710_jprb, &
      &  2720_jprb, 2730_jprb, 2740_jprb, 2750_jprb, 2760_jprb, 2770_jprb, 2780_jprb, 2790_jprb, &
      &  2800_jprb, 2810_jprb, 2820_jprb, 2830_jprb, 2840_jprb, 2850_jprb, 2860_jprb, 2870_jprb, &
      &  2880_jprb, 2890_jprb, 2900_jprb, 2910_jprb, 2920_jprb, 2930_jprb, 2940_jprb, 2950_jprb, &
      &  2960_jprb, 2970_jprb, 2980_jprb, 2990_jprb, 3000_jprb, 3010_jprb, 3020_jprb, 3030_jprb, &
      &  3040_jprb, 3050_jprb, 3060_jprb, 3070_jprb, 3080_jprb, 3090_jprb, 3100_jprb, 3110_jprb, &
      &  3120_jprb, 3130_jprb, 3140_jprb, 3150_jprb, 3160_jprb, 3170_jprb, 3180_jprb, 3190_jprb, &
      &  3200_jprb, 3210_jprb, 3220_jprb, 3230_jprb, 3240_jprb, 3250_jprb, 3260_jprb, 3270_jprb, &
      &  3280_jprb, 3290_jprb, 3300_jprb, 3310_jprb, 3320_jprb, 3330_jprb, 3340_jprb, 3350_jprb, &
      &  3360_jprb, 3370_jprb, 3380_jprb, 3390_jprb, 3400_jprb, 3410_jprb, 3420_jprb, 3430_jprb, &
      &  3440_jprb, 3450_jprb, 3460_jprb, 3470_jprb, 3480_jprb, 3490_jprb, 3500_jprb, 3510_jprb, &
      &  3520_jprb, 3530_jprb, 3540_jprb, 3550_jprb, 3560_jprb, 3570_jprb, 3580_jprb, 3590_jprb, &
      &  3600_jprb, 3610_jprb, 3620_jprb, 3630_jprb, 3640_jprb, 3650_jprb, 3660_jprb, 3670_jprb, &
      &  3680_jprb, 3690_jprb, 3700_jprb, 3710_jprb, 3720_jprb, 3730_jprb, 3740_jprb, 3750_jprb, &
      &  3760_jprb, 3770_jprb, 3780_jprb, 3790_jprb, 3800_jprb, 3810_jprb, 3820_jprb, 3830_jprb, &
      &  3840_jprb, 3850_jprb, 3860_jprb, 3870_jprb, 3880_jprb, 3890_jprb, 3900_jprb, 3910_jprb, &
      &  3920_jprb, 3930_jprb, 3940_jprb, 3950_jprb, 3960_jprb, 3970_jprb, 3980_jprb, 3990_jprb, &
      &  4000_jprb, 4010_jprb, 4020_jprb, 4030_jprb, 4040_jprb, 4050_jprb, 4060_jprb, 4070_jprb, &
      &  4080_jprb, 4090_jprb, 4100_jprb, 4110_jprb, 4120_jprb, 4130_jprb, 4140_jprb, 4150_jprb, &
      &  4160_jprb, 4170_jprb, 4180_jprb, 4190_jprb, 4200_jprb, 4210_jprb, 4220_jprb, 4230_jprb, &
      &  4240_jprb, 4250_jprb, 4260_jprb, 4270_jprb, 4280_jprb, 4290_jprb, 4300_jprb, 4310_jprb, &
      &  4320_jprb, 4330_jprb, 4340_jprb, 4350_jprb, 4360_jprb, 4370_jprb, 4380_jprb, 4390_jprb, &
      &  4400_jprb, 4410_jprb, 4420_jprb, 4430_jprb, 4440_jprb, 4450_jprb, 4460_jprb, 4470_jprb, &
      &  4480_jprb, 4490_jprb, 4500_jprb, 4510_jprb, 4520_jprb, 4530_jprb, 4540_jprb, 4550_jprb, &
      &  4560_jprb, 4570_jprb, 4580_jprb, 4590_jprb, 4600_jprb, 4610_jprb, 4620_jprb, 4630_jprb, &
      &  4640_jprb, 4650_jprb, 4660_jprb, 4670_jprb, 4680_jprb, 4690_jprb, 4700_jprb, 4710_jprb, &
      &  4720_jprb, 4730_jprb, 4740_jprb, 4750_jprb, 4760_jprb, 4770_jprb, 4780_jprb, 4790_jprb, &
      &  4800_jprb, 4810_jprb, 4820_jprb, 4830_jprb, 4840_jprb, 4850_jprb, 4860_jprb, 4870_jprb, &
      &  4880_jprb, 4890_jprb, 4900_jprb, 4910_jprb, 4920_jprb, 4930_jprb, 4940_jprb, 4950_jprb, &
      &  4960_jprb, 4970_jprb, 4980_jprb, 4990_jprb, 5000_jprb, 5010_jprb, 5020_jprb, 5030_jprb, &
      &  5040_jprb, 5050_jprb, 5060_jprb, 5070_jprb, 5080_jprb, 5090_jprb, 5100_jprb, 5110_jprb, &
      &  5120_jprb, 5130_jprb, 5140_jprb, 5150_jprb, 5160_jprb, 5170_jprb, 5180_jprb, 5190_jprb, &
      &  5200_jprb, 5210_jprb, 5220_jprb, 5230_jprb, 5240_jprb, 5250_jprb, 5260_jprb, 5270_jprb, &
      &  5280_jprb, 5290_jprb, 5300_jprb, 5310_jprb, 5320_jprb, 5330_jprb, 5340_jprb, 5350_jprb, &
      &  5360_jprb, 5370_jprb, 5380_jprb, 5390_jprb, 5400_jprb, 5410_jprb, 5420_jprb, 5430_jprb, &
      &  5440_jprb, 5450_jprb, 5460_jprb, 5470_jprb, 5480_jprb, 5490_jprb, 5500_jprb, 5510_jprb, &
      &  5520_jprb, 5530_jprb, 5540_jprb, 5550_jprb, 5560_jprb, 5570_jprb, 5580_jprb, 5590_jprb, &
      &  5600_jprb, 5610_jprb, 5620_jprb, 5630_jprb, 5640_jprb, 5650_jprb, 5660_jprb, 5670_jprb, &
      &  5680_jprb, 5690_jprb, 5700_jprb, 5710_jprb, 5720_jprb, 5730_jprb, 5740_jprb, 5750_jprb, &
      &  5760_jprb, 5770_jprb, 5780_jprb, 5790_jprb, 5800_jprb, 5810_jprb, 5820_jprb, 5830_jprb, &
      &  5840_jprb, 5850_jprb, 5860_jprb, 5870_jprb, 5880_jprb, 5890_jprb, 5900_jprb, 5910_jprb, &
      &  5920_jprb, 5930_jprb, 5940_jprb, 5950_jprb, 5960_jprb, 5970_jprb, 5980_jprb, 5990_jprb, &
      &  6000_jprb, 6010_jprb, 6020_jprb, 6030_jprb, 6040_jprb, 6050_jprb, 6060_jprb, 6070_jprb, &
      &  6080_jprb, 6090_jprb, 6100_jprb, 6110_jprb, 6120_jprb, 6130_jprb, 6140_jprb, 6150_jprb, &
      &  6160_jprb, 6170_jprb, 6180_jprb, 6190_jprb, 6200_jprb, 6210_jprb, 6220_jprb, 6230_jprb, &
      &  6240_jprb, 6250_jprb, 6260_jprb, 6270_jprb, 6280_jprb, 6290_jprb, 6300_jprb, 6310_jprb, &
      &  6320_jprb, 6330_jprb, 6340_jprb, 6350_jprb, 6360_jprb, 6370_jprb, 6380_jprb, 6390_jprb, &
      &  6400_jprb, 6410_jprb, 6420_jprb, 6430_jprb, 6440_jprb, 6450_jprb, 6460_jprb, 6470_jprb, &
      &  6480_jprb, 6490_jprb, 6500_jprb, 6510_jprb, 6520_jprb, 6530_jprb, 6540_jprb, 6550_jprb, &
      &  6560_jprb, 6570_jprb, 6580_jprb, 6590_jprb, 6600_jprb, 6610_jprb, 6620_jprb, 6630_jprb, &
      &  6640_jprb, 6650_jprb, 6660_jprb, 6670_jprb, 6680_jprb, 6690_jprb, 6700_jprb, 6710_jprb, &
      &  6720_jprb, 6730_jprb, 6740_jprb, 6750_jprb, 6760_jprb, 6770_jprb, 6780_jprb, 6790_jprb, &
      &  6800_jprb, 6810_jprb, 6820_jprb, 6830_jprb, 6840_jprb, 6850_jprb, 6860_jprb, 6870_jprb, &
      &  6880_jprb, 6890_jprb, 6900_jprb, 6910_jprb, 6920_jprb, 6930_jprb, 6940_jprb, 6950_jprb, &
      &  6960_jprb, 6970_jprb, 6980_jprb, 6990_jprb, 7000_jprb, 7010_jprb, 7020_jprb, 7030_jprb, &
      &  7040_jprb, 7050_jprb, 7060_jprb, 7070_jprb, 7080_jprb, 7090_jprb, 7100_jprb, 7110_jprb, &
      &  7120_jprb, 7130_jprb, 7140_jprb, 7150_jprb, 7160_jprb, 7170_jprb, 7180_jprb, 7190_jprb, &
      &  7200_jprb, 7210_jprb, 7220_jprb, 7230_jprb, 7240_jprb, 7250_jprb, 7260_jprb, 7270_jprb, &
      &  7280_jprb, 7290_jprb, 7300_jprb, 7310_jprb, 7320_jprb, 7330_jprb, 7340_jprb, 7350_jprb, &
      &  7360_jprb, 7370_jprb, 7380_jprb, 7390_jprb, 7400_jprb, 7410_jprb, 7420_jprb, 7430_jprb, &
      &  7440_jprb, 7450_jprb, 7460_jprb, 7470_jprb, 7480_jprb, 7490_jprb, 7500_jprb, 7510_jprb, &
      &  7520_jprb, 7530_jprb, 7540_jprb, 7550_jprb, 7560_jprb, 7570_jprb, 7580_jprb, 7590_jprb, &
      &  7600_jprb, 7610_jprb, 7620_jprb, 7630_jprb, 7640_jprb, 7650_jprb, 7660_jprb, 7670_jprb, &
      &  7680_jprb, 7690_jprb, 7700_jprb, 7710_jprb, 7720_jprb, 7730_jprb, 7740_jprb, 7750_jprb, &
      &  7760_jprb, 7770_jprb, 7780_jprb, 7790_jprb, 7800_jprb, 7810_jprb, 7820_jprb, 7830_jprb, &
      &  7840_jprb, 7850_jprb, 7860_jprb, 7870_jprb, 7880_jprb, 7890_jprb, 7900_jprb, 7910_jprb, &
      &  7920_jprb, 7930_jprb, 7940_jprb, 7950_jprb, 7960_jprb, 7970_jprb, 7980_jprb, 7990_jprb, &
      &  8000_jprb, 8010_jprb, 8020_jprb, 8030_jprb, 8040_jprb, 8050_jprb, 8060_jprb, 8070_jprb, &
      &  8080_jprb, 8090_jprb, 8100_jprb, 8110_jprb, 8120_jprb, 8130_jprb, 8140_jprb, 8150_jprb, &
      &  8160_jprb, 8170_jprb, 8180_jprb, 8190_jprb, 8200_jprb, 8210_jprb, 8220_jprb, 8230_jprb, &
      &  8240_jprb, 8250_jprb, 8260_jprb, 8270_jprb, 8280_jprb, 8290_jprb, 8300_jprb, 8310_jprb, &
      &  8320_jprb, 8330_jprb, 8340_jprb, 8350_jprb, 8360_jprb, 8370_jprb, 8380_jprb, 8390_jprb, &
      &  8400_jprb, 8410_jprb, 8420_jprb, 8430_jprb, 8440_jprb, 8450_jprb, 8460_jprb, 8470_jprb, &
      &  8480_jprb, 8490_jprb, 8500_jprb, 8510_jprb, 8520_jprb, 8530_jprb, 8540_jprb, 8550_jprb, &
      &  8560_jprb, 8570_jprb, 8580_jprb, 8590_jprb, 8600_jprb, 8610_jprb, 8620_jprb, 8630_jprb, &
      &  8640_jprb, 8650_jprb, 8660_jprb, 8670_jprb, 8680_jprb, 8690_jprb, 8700_jprb, 8710_jprb, &
      &  8720_jprb, 8730_jprb, 8740_jprb, 8750_jprb, 8760_jprb, 8770_jprb, 8780_jprb, 8790_jprb, &
      &  8800_jprb, 8810_jprb, 8820_jprb, 8830_jprb, 8840_jprb, 8850_jprb, 8860_jprb, 8870_jprb, &
      &  8880_jprb, 8890_jprb, 8900_jprb, 8910_jprb, 8920_jprb, 8930_jprb, 8940_jprb, 8950_jprb, &
      &  8960_jprb, 8970_jprb, 8980_jprb, 8990_jprb, 9000_jprb, 9010_jprb, 9020_jprb, 9030_jprb, &
      &  9040_jprb, 9050_jprb, 9060_jprb, 9070_jprb, 9080_jprb, 9090_jprb, 9100_jprb, 9110_jprb, &
      &  9120_jprb, 9130_jprb, 9140_jprb, 9150_jprb, 9160_jprb, 9170_jprb, 9180_jprb, 9190_jprb, &
      &  9200_jprb, 9210_jprb, 9220_jprb, 9230_jprb, 9240_jprb, 9250_jprb, 9260_jprb, 9270_jprb, &
      &  9280_jprb, 9290_jprb, 9300_jprb, 9310_jprb, 9320_jprb, 9330_jprb, 9340_jprb, 9350_jprb, &
      &  9360_jprb, 9370_jprb, 9380_jprb, 9390_jprb, 9400_jprb, 9410_jprb, 9420_jprb, 9430_jprb, &
      &  9440_jprb, 9450_jprb, 9460_jprb, 9470_jprb, 9480_jprb, 9490_jprb, 9500_jprb, 9510_jprb, &
      &  9520_jprb, 9530_jprb, 9540_jprb, 9550_jprb, 9560_jprb, 9570_jprb, 9580_jprb, 9590_jprb, &
      &  9600_jprb, 9610_jprb, 9620_jprb, 9630_jprb, 9640_jprb, 9650_jprb, 9660_jprb, 9670_jprb, &
      &  9680_jprb, 9690_jprb, 9700_jprb, 9710_jprb, 9720_jprb, 9730_jprb, 9740_jprb, 9750_jprb, &
      &  9760_jprb, 9770_jprb, 9780_jprb, 9790_jprb, 9800_jprb, 9810_jprb, 9820_jprb, 9830_jprb, &
      &  9840_jprb, 9850_jprb, 9860_jprb, 9870_jprb, 9880_jprb, 9890_jprb, 9900_jprb, 9910_jprb, &
      &  9920_jprb, 9930_jprb, 9940_jprb, 9950_jprb, 9960_jprb, 9970_jprb, 9980_jprb, 9990_jprb, &
      &  10000_jprb, 10010_jprb, 10020_jprb, 10030_jprb, 10040_jprb, 10050_jprb, 10060_jprb, 10070_jprb, &
      &  10080_jprb, 10090_jprb, 10100_jprb, 10110_jprb, 10120_jprb, 10130_jprb, 10140_jprb, 10150_jprb, &
      &  10160_jprb, 10170_jprb, 10180_jprb, 10190_jprb, 10200_jprb, 10210_jprb, 10220_jprb, 10230_jprb, &
      &  10240_jprb, 10250_jprb, 10260_jprb, 10270_jprb, 10280_jprb, 10290_jprb, 10300_jprb, 10310_jprb, &
      &  10320_jprb, 10330_jprb, 10340_jprb, 10350_jprb, 10360_jprb, 10370_jprb, 10380_jprb, 10390_jprb, &
      &  10400_jprb, 10410_jprb, 10420_jprb, 10430_jprb, 10440_jprb, 10450_jprb, 10460_jprb, 10470_jprb, &
      &  10480_jprb, 10490_jprb, 10500_jprb, 10510_jprb, 10520_jprb, 10530_jprb, 10540_jprb, 10550_jprb, &
      &  10560_jprb, 10570_jprb, 10580_jprb, 10590_jprb, 10600_jprb, 10610_jprb, 10620_jprb, 10630_jprb, &
      &  10640_jprb, 10650_jprb, 10660_jprb, 10670_jprb, 10680_jprb, 10690_jprb, 10700_jprb, 10710_jprb, &
      &  10720_jprb, 10730_jprb, 10740_jprb, 10750_jprb, 10760_jprb, 10770_jprb, 10780_jprb, 10790_jprb, &
      &  10800_jprb, 10810_jprb, 10820_jprb, 10830_jprb, 10840_jprb, 10850_jprb, 10860_jprb, 10870_jprb, &
      &  10880_jprb, 10890_jprb, 10900_jprb, 10910_jprb, 10920_jprb, 10930_jprb, 10940_jprb, 10950_jprb, &
      &  10960_jprb, 10970_jprb, 10980_jprb, 10990_jprb, 11000_jprb, 11010_jprb, 11020_jprb, 11030_jprb, &
      &  11040_jprb, 11050_jprb, 11060_jprb, 11070_jprb, 11080_jprb, 11090_jprb, 11100_jprb, 11110_jprb, &
      &  11120_jprb, 11130_jprb, 11140_jprb, 11150_jprb, 11160_jprb, 11170_jprb, 11180_jprb, 11190_jprb, &
      &  11200_jprb, 11210_jprb, 11220_jprb, 11230_jprb, 11240_jprb, 11250_jprb, 11260_jprb, 11270_jprb, &
      &  11280_jprb, 11290_jprb, 11300_jprb, 11310_jprb, 11320_jprb, 11330_jprb, 11340_jprb, 11350_jprb, &
      &  11360_jprb, 11370_jprb, 11380_jprb, 11390_jprb, 11400_jprb, 11410_jprb, 11420_jprb, 11430_jprb, &
      &  11440_jprb, 11450_jprb, 11460_jprb, 11470_jprb, 11480_jprb, 11490_jprb, 11500_jprb, 11510_jprb, &
      &  11520_jprb, 11530_jprb, 11540_jprb, 11550_jprb, 11560_jprb, 11570_jprb, 11580_jprb, 11590_jprb, &
      &  11600_jprb, 11610_jprb, 11620_jprb, 11630_jprb, 11640_jprb, 11650_jprb, 11660_jprb, 11670_jprb, &
      &  11680_jprb, 11690_jprb, 11700_jprb, 11710_jprb, 11720_jprb, 11730_jprb, 11740_jprb, 11750_jprb, &
      &  11760_jprb, 11770_jprb, 11780_jprb, 11790_jprb, 11800_jprb, 11810_jprb, 11820_jprb, 11830_jprb, &
      &  11840_jprb, 11850_jprb, 11860_jprb, 11870_jprb, 11880_jprb, 11890_jprb, 11900_jprb, 11910_jprb, &
      &  11920_jprb, 11930_jprb, 11940_jprb, 11950_jprb, 11960_jprb, 11970_jprb, 11980_jprb, 11990_jprb, &
      &  12000_jprb, 12010_jprb, 12020_jprb, 12030_jprb, 12040_jprb, 12050_jprb, 12060_jprb, 12070_jprb, &
      &  12080_jprb, 12090_jprb, 12100_jprb, 12110_jprb, 12120_jprb, 12130_jprb, 12140_jprb, 12150_jprb, &
      &  12160_jprb, 12170_jprb, 12180_jprb, 12190_jprb, 12200_jprb, 12210_jprb, 12220_jprb, 12230_jprb, &
      &  12240_jprb, 12250_jprb, 12260_jprb, 12270_jprb, 12280_jprb, 12290_jprb, 12300_jprb, 12310_jprb, &
      &  12320_jprb, 12330_jprb, 12340_jprb, 12350_jprb, 12360_jprb, 12370_jprb, 12380_jprb, 12390_jprb, &
      &  12400_jprb, 12410_jprb, 12420_jprb, 12430_jprb, 12440_jprb, 12450_jprb, 12460_jprb, 12470_jprb, &
      &  12480_jprb, 12490_jprb, 12500_jprb, 12510_jprb, 12520_jprb, 12530_jprb, 12540_jprb, 12550_jprb, &
      &  12560_jprb, 12570_jprb, 12580_jprb, 12590_jprb, 12600_jprb, 12610_jprb, 12620_jprb, 12630_jprb, &
      &  12640_jprb, 12650_jprb, 12660_jprb, 12670_jprb, 12680_jprb, 12690_jprb, 12700_jprb, 12710_jprb, &
      &  12720_jprb, 12730_jprb, 12740_jprb, 12750_jprb, 12760_jprb, 12770_jprb, 12780_jprb, 12790_jprb, &
      &  12800_jprb, 12810_jprb, 12820_jprb, 12830_jprb, 12840_jprb, 12850_jprb, 12860_jprb, 12870_jprb, &
      &  12880_jprb, 12890_jprb, 12900_jprb, 12910_jprb, 12920_jprb, 12930_jprb, 12940_jprb, 12950_jprb, &
      &  12960_jprb, 12970_jprb, 12980_jprb, 12990_jprb, 13000_jprb, 13010_jprb, 13020_jprb, 13030_jprb, &
      &  13040_jprb, 13050_jprb, 13060_jprb, 13070_jprb, 13080_jprb, 13090_jprb, 13100_jprb, 13110_jprb, &
      &  13120_jprb, 13130_jprb, 13140_jprb, 13150_jprb, 13160_jprb, 13170_jprb, 13180_jprb, 13190_jprb, &
      &  13200_jprb, 13210_jprb, 13220_jprb, 13230_jprb, 13240_jprb, 13250_jprb, 13260_jprb, 13270_jprb, &
      &  13280_jprb, 13290_jprb, 13300_jprb, 13310_jprb, 13320_jprb, 13330_jprb, 13340_jprb, 13350_jprb, &
      &  13360_jprb, 13370_jprb, 13380_jprb, 13390_jprb, 13400_jprb, 13410_jprb, 13420_jprb, 13430_jprb, &
      &  13440_jprb, 13450_jprb, 13460_jprb, 13470_jprb, 13480_jprb, 13490_jprb, 13500_jprb, 13510_jprb, &
      &  13520_jprb, 13530_jprb, 13540_jprb, 13550_jprb, 13560_jprb, 13570_jprb, 13580_jprb, 13590_jprb, &
      &  13600_jprb, 13610_jprb, 13620_jprb, 13630_jprb, 13640_jprb, 13650_jprb, 13660_jprb, 13670_jprb, &
      &  13680_jprb, 13690_jprb, 13700_jprb, 13710_jprb, 13720_jprb, 13730_jprb, 13740_jprb, 13750_jprb, &
      &  13760_jprb, 13770_jprb, 13780_jprb, 13790_jprb, 13800_jprb, 13810_jprb, 13820_jprb, 13830_jprb, &
      &  13840_jprb, 13850_jprb, 13860_jprb, 13870_jprb, 13880_jprb, 13890_jprb, 13900_jprb, 13910_jprb, &
      &  13920_jprb, 13930_jprb, 13940_jprb, 13950_jprb, 13960_jprb, 13970_jprb, 13980_jprb, 13990_jprb, &
      &  14000_jprb, 14010_jprb, 14020_jprb, 14030_jprb, 14040_jprb, 14050_jprb, 14060_jprb, 14070_jprb, &
      &  14080_jprb, 14090_jprb, 14100_jprb, 14110_jprb, 14120_jprb, 14130_jprb, 14140_jprb, 14150_jprb, &
      &  14160_jprb, 14170_jprb, 14180_jprb, 14190_jprb, 14200_jprb, 14210_jprb, 14220_jprb, 14230_jprb, &
      &  14240_jprb, 14250_jprb, 14260_jprb, 14270_jprb, 14280_jprb, 14290_jprb, 14300_jprb, 14310_jprb, &
      &  14320_jprb, 14330_jprb, 14340_jprb, 14350_jprb, 14360_jprb, 14370_jprb, 14380_jprb, 14390_jprb, &
      &  14400_jprb, 14410_jprb, 14420_jprb, 14430_jprb, 14440_jprb, 14450_jprb, 14460_jprb, 14470_jprb, &
      &  14480_jprb, 14490_jprb, 14500_jprb, 14510_jprb, 14520_jprb, 14530_jprb, 14540_jprb, 14550_jprb, &
      &  14560_jprb, 14570_jprb, 14580_jprb, 14590_jprb, 14600_jprb, 14610_jprb, 14620_jprb, 14630_jprb, &
      &  14640_jprb, 14650_jprb, 14660_jprb, 14670_jprb, 14680_jprb, 14690_jprb, 14700_jprb, 14710_jprb, &
      &  14720_jprb, 14730_jprb, 14740_jprb, 14750_jprb, 14760_jprb, 14770_jprb, 14780_jprb, 14790_jprb, &
      &  14800_jprb, 14810_jprb, 14820_jprb, 14830_jprb, 14840_jprb, 14850_jprb, 14860_jprb, 14870_jprb, &
      &  14880_jprb, 14890_jprb, 14900_jprb, 14910_jprb, 14920_jprb, 14930_jprb, 14940_jprb, 14950_jprb, &
      &  14960_jprb, 14970_jprb, 14980_jprb, 14990_jprb, 15000_jprb, 15010_jprb, 15020_jprb, 15030_jprb, &
      &  15040_jprb, 15050_jprb, 15060_jprb, 15070_jprb, 15080_jprb, 15090_jprb, 15100_jprb, 15110_jprb, &
      &  15120_jprb, 15130_jprb, 15140_jprb, 15150_jprb, 15160_jprb, 15170_jprb, 15180_jprb, 15190_jprb, &
      &  15200_jprb, 15210_jprb, 15220_jprb, 15230_jprb, 15240_jprb, 15250_jprb, 15260_jprb, 15270_jprb, &
      &  15280_jprb, 15290_jprb, 15300_jprb, 15310_jprb, 15320_jprb, 15330_jprb, 15340_jprb, 15350_jprb, &
      &  15360_jprb, 15370_jprb, 15380_jprb, 15390_jprb, 15400_jprb, 15410_jprb, 15420_jprb, 15430_jprb, &
      &  15440_jprb, 15450_jprb, 15460_jprb, 15470_jprb, 15480_jprb, 15490_jprb, 15500_jprb, 15510_jprb, &
      &  15520_jprb, 15530_jprb, 15540_jprb, 15550_jprb, 15560_jprb, 15570_jprb, 15580_jprb, 15590_jprb, &
      &  15600_jprb, 15610_jprb, 15620_jprb, 15630_jprb, 15640_jprb, 15650_jprb, 15660_jprb, 15670_jprb, &
      &  15680_jprb, 15690_jprb, 15700_jprb, 15710_jprb, 15720_jprb, 15730_jprb, 15740_jprb, 15750_jprb, &
      &  15760_jprb, 15770_jprb, 15780_jprb, 15790_jprb, 15800_jprb, 15810_jprb, 15820_jprb, 15830_jprb, &
      &  15840_jprb, 15850_jprb, 15860_jprb, 15870_jprb, 15880_jprb, 15890_jprb, 15900_jprb, 15910_jprb, &
      &  15920_jprb, 15930_jprb, 15940_jprb, 15950_jprb, 15960_jprb, 15970_jprb, 15980_jprb, 15990_jprb, &
      &  16000_jprb, 16010_jprb, 16020_jprb, 16030_jprb, 16040_jprb, 16050_jprb, 16060_jprb, 16070_jprb, &
      &  16080_jprb, 16090_jprb, 16100_jprb, 16110_jprb, 16120_jprb, 16130_jprb, 16140_jprb, 16150_jprb, &
      &  16160_jprb, 16170_jprb, 16180_jprb, 16190_jprb, 16200_jprb, 16210_jprb, 16220_jprb, 16230_jprb, &
      &  16240_jprb, 16250_jprb, 16260_jprb, 16270_jprb, 16280_jprb, 16290_jprb, 16300_jprb, 16310_jprb, &
      &  16320_jprb, 16330_jprb, 16340_jprb, 16350_jprb, 16360_jprb, 16370_jprb, 16380_jprb, 16390_jprb, &
      &  16400_jprb, 16410_jprb, 16420_jprb, 16430_jprb, 16440_jprb, 16450_jprb, 16460_jprb, 16470_jprb, &
      &  16480_jprb, 16490_jprb, 16500_jprb, 16510_jprb, 16520_jprb, 16530_jprb, 16540_jprb, 16550_jprb, &
      &  16560_jprb, 16570_jprb, 16580_jprb, 16590_jprb, 16600_jprb, 16610_jprb, 16620_jprb, 16630_jprb, &
      &  16640_jprb, 16650_jprb, 16660_jprb, 16670_jprb, 16680_jprb, 16690_jprb, 16700_jprb, 16710_jprb, &
      &  16720_jprb, 16730_jprb, 16740_jprb, 16750_jprb, 16760_jprb, 16770_jprb, 16780_jprb, 16790_jprb, &
      &  16800_jprb, 16810_jprb, 16820_jprb, 16830_jprb, 16840_jprb, 16850_jprb, 16860_jprb, 16870_jprb, &
      &  16880_jprb, 16890_jprb, 16900_jprb, 16910_jprb, 16920_jprb, 16930_jprb, 16940_jprb, 16950_jprb, &
      &  16960_jprb, 16970_jprb, 16980_jprb, 16990_jprb, 17000_jprb, 17010_jprb, 17020_jprb, 17030_jprb, &
      &  17040_jprb, 17050_jprb, 17060_jprb, 17070_jprb, 17080_jprb, 17090_jprb, 17100_jprb, 17110_jprb, &
      &  17120_jprb, 17130_jprb, 17140_jprb, 17150_jprb, 17160_jprb, 17170_jprb, 17180_jprb, 17190_jprb, &
      &  17200_jprb, 17210_jprb, 17220_jprb, 17230_jprb, 17240_jprb, 17250_jprb, 17260_jprb, 17270_jprb, &
      &  17280_jprb, 17290_jprb, 17300_jprb, 17310_jprb, 17320_jprb, 17330_jprb, 17340_jprb, 17350_jprb, &
      &  17360_jprb, 17370_jprb, 17380_jprb, 17390_jprb, 17400_jprb, 17410_jprb, 17420_jprb, 17430_jprb, &
      &  17440_jprb, 17450_jprb, 17460_jprb, 17470_jprb, 17480_jprb, 17490_jprb, 17500_jprb, 17510_jprb, &
      &  17520_jprb, 17530_jprb, 17540_jprb, 17550_jprb, 17560_jprb, 17570_jprb, 17580_jprb, 17590_jprb, &
      &  17600_jprb, 17610_jprb, 17620_jprb, 17630_jprb, 17640_jprb, 17650_jprb, 17660_jprb, 17670_jprb, &
      &  17680_jprb, 17690_jprb, 17700_jprb, 17710_jprb, 17720_jprb, 17730_jprb, 17740_jprb, 17750_jprb, &
      &  17760_jprb, 17770_jprb, 17780_jprb, 17790_jprb, 17800_jprb, 17810_jprb, 17820_jprb, 17830_jprb, &
      &  17840_jprb, 17850_jprb, 17860_jprb, 17870_jprb, 17880_jprb, 17890_jprb, 17900_jprb, 17910_jprb, &
      &  17920_jprb, 17930_jprb, 17940_jprb, 17950_jprb, 17960_jprb, 17970_jprb, 17980_jprb, 17990_jprb, &
      &  18000_jprb, 18010_jprb, 18020_jprb, 18030_jprb, 18040_jprb, 18050_jprb, 18060_jprb, 18070_jprb, &
      &  18080_jprb, 18090_jprb, 18100_jprb, 18110_jprb, 18120_jprb, 18130_jprb, 18140_jprb, 18150_jprb, &
      &  18160_jprb, 18170_jprb, 18180_jprb, 18190_jprb, 18200_jprb, 18210_jprb, 18220_jprb, 18230_jprb, &
      &  18240_jprb, 18250_jprb, 18260_jprb, 18270_jprb, 18280_jprb, 18290_jprb, 18300_jprb, 18310_jprb, &
      &  18320_jprb, 18330_jprb, 18340_jprb, 18350_jprb, 18360_jprb, 18370_jprb, 18380_jprb, 18390_jprb, &
      &  18400_jprb, 18410_jprb, 18420_jprb, 18430_jprb, 18440_jprb, 18450_jprb, 18460_jprb, 18470_jprb, &
      &  18480_jprb, 18490_jprb, 18500_jprb, 18510_jprb, 18520_jprb, 18530_jprb, 18540_jprb, 18550_jprb, &
      &  18560_jprb, 18570_jprb, 18580_jprb, 18590_jprb, 18600_jprb, 18610_jprb, 18620_jprb, 18630_jprb, &
      &  18640_jprb, 18650_jprb, 18660_jprb, 18670_jprb, 18680_jprb, 18690_jprb, 18700_jprb, 18710_jprb, &
      &  18720_jprb, 18730_jprb, 18740_jprb, 18750_jprb, 18760_jprb, 18770_jprb, 18780_jprb, 18790_jprb, &
      &  18800_jprb, 18810_jprb, 18820_jprb, 18830_jprb, 18840_jprb, 18850_jprb, 18860_jprb, 18870_jprb, &
      &  18880_jprb, 18890_jprb, 18900_jprb, 18910_jprb, 18920_jprb, 18930_jprb, 18940_jprb, 18950_jprb, &
      &  18960_jprb, 18970_jprb, 18980_jprb, 18990_jprb, 19000_jprb, 19010_jprb, 19020_jprb, 19030_jprb, &
      &  19040_jprb, 19050_jprb, 19060_jprb, 19070_jprb, 19080_jprb, 19090_jprb, 19100_jprb, 19110_jprb, &
      &  19120_jprb, 19130_jprb, 19140_jprb, 19150_jprb, 19160_jprb, 19170_jprb, 19180_jprb, 19190_jprb, &
      &  19200_jprb, 19210_jprb, 19220_jprb, 19230_jprb, 19240_jprb, 19250_jprb, 19260_jprb, 19270_jprb, &
      &  19280_jprb, 19290_jprb, 19300_jprb, 19310_jprb, 19320_jprb, 19330_jprb, 19340_jprb, 19350_jprb, &
      &  19360_jprb, 19370_jprb, 19380_jprb, 19390_jprb, 19400_jprb, 19410_jprb, 19420_jprb, 19430_jprb, &
      &  19440_jprb, 19450_jprb, 19460_jprb, 19470_jprb, 19480_jprb, 19490_jprb, 19500_jprb, 19510_jprb, &
      &  19520_jprb, 19530_jprb, 19540_jprb, 19550_jprb, 19560_jprb, 19570_jprb, 19580_jprb, 19590_jprb, &
      &  19600_jprb, 19610_jprb, 19620_jprb, 19630_jprb, 19640_jprb, 19650_jprb, 19660_jprb, 19670_jprb, &
      &  19680_jprb, 19690_jprb, 19700_jprb, 19710_jprb, 19720_jprb, 19730_jprb, 19740_jprb, 19750_jprb, &
      &  19760_jprb, 19770_jprb, 19780_jprb, 19790_jprb, 19800_jprb, 19810_jprb, 19820_jprb, 19830_jprb, &
      &  19840_jprb, 19850_jprb, 19860_jprb, 19870_jprb, 19880_jprb, 19890_jprb, 19900_jprb, 19910_jprb, &
      &  19920_jprb, 19930_jprb, 19940_jprb, 19950_jprb, 19960_jprb, 19970_jprb, 19980_jprb, 19990_jprb, &
      &  20000_jprb, 20010_jprb ]
    real(jprb), parameter :: self_lut(nwav_lut) = [ &
      &  0.0_jprb, 5.5601e-22_jprb, 2.253e-21_jprb, 4.8967e-21_jprb, 8.5013e-21_jprb, 1.2478e-20_jprb, &
      &  1.7553e-20_jprb, 2.2294e-20_jprb, 2.7818e-20_jprb, 3.2573e-20_jprb, 3.6179e-20_jprb, 4.0683e-20_jprb, &
      &  4.2921e-20_jprb, 4.5786e-20_jprb, 4.7545e-20_jprb, 4.5708e-20_jprb, 4.8159e-20_jprb, 4.5341e-20_jprb, &
      &  4.6647e-20_jprb, 4.5043e-20_jprb, 4.2331e-20_jprb, 4.2622e-20_jprb, 3.979e-20_jprb, 3.9074e-20_jprb, &
      &  3.5933e-20_jprb, 3.415e-20_jprb, 3.2099e-20_jprb, 2.9862e-20_jprb, 2.7895e-20_jprb, 2.6198e-20_jprb, &
      &  2.378e-20_jprb, 2.2268e-20_jprb, 2.0331e-20_jprb, 1.8928e-20_jprb, 1.7176e-20_jprb, 1.6008e-20_jprb, &
      &  1.4504e-20_jprb, 1.3196e-20_jprb, 1.2074e-20_jprb, 1.0978e-20_jprb, 1.0061e-20_jprb, 9.1346e-21_jprb, &
      &  8.4084e-21_jprb, 7.5948e-21_jprb, 6.913e-21_jprb, 6.2929e-21_jprb, 5.7394e-21_jprb, 5.2098e-21_jprb, &
      &  4.7495e-21_jprb, 4.3204e-21_jprb, 3.9268e-21_jprb, 3.5772e-21_jprb, 3.264e-21_jprb, 2.9682e-21_jprb, &
      &  2.7019e-21_jprb, 2.4628e-21_jprb, 2.2456e-21_jprb, 2.0516e-21_jprb, 1.8722e-21_jprb, 1.7132e-21_jprb, &
      &  1.5701e-21_jprb, 1.4411e-21_jprb, 1.3236e-21_jprb, 1.218e-21_jprb, 1.1236e-21_jprb, 1.0379e-21_jprb, &
      &  9.6086e-22_jprb, 8.9145e-22_jprb, 8.2835e-22_jprb, 7.7208e-22_jprb, 7.2118e-22_jprb, 6.7473e-22_jprb, &
      &  6.3286e-22_jprb, 5.9474e-22_jprb, 5.5993e-22_jprb, 5.286e-22_jprb, 4.9984e-22_jprb, 4.7345e-22_jprb, &
      &  4.4922e-22_jprb, 4.2693e-22_jprb, 4.0639e-22_jprb, 3.8725e-22_jprb, 3.7055e-22_jprb, 3.5601e-22_jprb, &
      &  3.4238e-22_jprb, 3.3024e-22_jprb, 3.1809e-22_jprb, 3.0601e-22_jprb, 2.9461e-22_jprb, 2.8322e-22_jprb, &
      &  2.7204e-22_jprb, 2.6274e-22_jprb, 2.4995e-22_jprb, 2.381e-22_jprb, 2.251e-22_jprb, 2.147e-22_jprb, &
      &  2.0607e-22_jprb, 1.9819e-22_jprb, 1.9076e-22_jprb, 1.836e-22_jprb, 1.7682e-22_jprb, 1.7035e-22_jprb, &
      &  1.6431e-22_jprb, 1.585e-22_jprb, 1.5305e-22_jprb, 1.4797e-22_jprb, 1.4328e-22_jprb, 1.3889e-22_jprb, &
      &  1.3502e-22_jprb, 1.3148e-22_jprb, 1.2838e-22_jprb, 1.2576e-22_jprb, 1.235e-22_jprb, 1.2196e-22_jprb, &
      &  1.2105e-22_jprb, 1.2067e-22_jprb, 1.2096e-22_jprb, 1.2216e-22_jprb, 1.2418e-22_jprb, 1.2728e-22_jprb, &
      &  1.3125e-22_jprb, 1.3635e-22_jprb, 1.4322e-22_jprb, 1.514e-22_jprb, 1.6143e-22_jprb, 1.7335e-22_jprb, &
      &  2.044e-22_jprb, 2.1547e-22_jprb, 2.3326e-22_jprb, 2.5656e-22_jprb, 2.8449e-22_jprb, 3.1663e-22_jprb, &
      &  3.5336e-22_jprb, 3.971e-22_jprb, 4.5312e-22_jprb, 5.2761e-22_jprb, 6.2334e-22_jprb, 7.3745e-22_jprb, &
      &  8.6516e-22_jprb, 1.0062e-21_jprb, 1.168e-21_jprb, 1.3627e-21_jprb, 1.6019e-21_jprb, 1.8941e-21_jprb, &
      &  2.2422e-21_jprb, 2.6394e-21_jprb, 3.0675e-21_jprb, 3.5086e-21_jprb, 3.9574e-21_jprb, 4.4142e-21_jprb, &
      &  4.8677e-21_jprb, 5.2904e-21_jprb, 5.6475e-21_jprb, 5.9009e-21_jprb, 6.0173e-21_jprb, 5.9974e-21_jprb, &
      &  5.905e-21_jprb, 5.8444e-21_jprb, 5.8906e-21_jprb, 6.0341e-21_jprb, 6.1926e-21_jprb, 6.2745e-21_jprb, &
      &  6.2389e-21_jprb, 6.11e-21_jprb, 5.9466e-21_jprb, 5.8031e-21_jprb, 5.7119e-21_jprb, 5.6776e-21_jprb, &
      &  5.6704e-21_jprb, 5.6334e-21_jprb, 5.5153e-21_jprb, 5.3008e-21_jprb, 5.0091e-21_jprb, 4.6695e-21_jprb, &
      &  4.3038e-21_jprb, 3.9294e-21_jprb, 3.5626e-21_jprb, 3.2098e-21_jprb, 2.862e-21_jprb, 2.5048e-21_jprb, &
      &  2.1378e-21_jprb, 1.783e-21_jprb, 1.4707e-21_jprb, 1.2177e-21_jprb, 1.0188e-21_jprb, 8.5763e-22_jprb, &
      &  7.2193e-22_jprb, 6.0744e-22_jprb, 5.1391e-22_jprb, 4.4097e-22_jprb, 3.8657e-22_jprb, 3.4695e-22_jprb, &
      &  3.1763e-22_jprb, 2.9483e-22_jprb, 2.7578e-22_jprb, 2.5792e-22_jprb, 2.3902e-22_jprb, 2.1792e-22_jprb, &
      &  1.9475e-22_jprb, 1.7005e-22_jprb, 1.4451e-22_jprb, 1.1943e-22_jprb, 9.6806e-23_jprb, 7.8393e-23_jprb, &
      &  6.4697e-23_jprb, 5.4949e-23_jprb, 4.7872e-23_jprb, 4.2391e-23_jprb, 3.7912e-23_jprb, 3.4198e-23_jprb, &
      &  3.111e-23_jprb, 2.8428e-23_jprb, 2.5901e-23_jprb, 2.3406e-23_jprb, 2.1025e-23_jprb, 1.8943e-23_jprb, &
      &  1.7248e-23_jprb, 1.586e-23_jprb, 1.4629e-23_jprb, 1.3485e-23_jprb, 1.2458e-23_jprb, 1.1606e-23_jprb, &
      &  1.0943e-23_jprb, 1.0429e-23_jprb, 1.0009e-23_jprb, 9.6413e-24_jprb, 9.306e-24_jprb, 9.0006e-24_jprb, &
      &  8.7331e-24_jprb, 8.5105e-24_jprb, 8.3277e-24_jprb, 8.1707e-24_jprb, 8.0243e-24_jprb, 7.8742e-24_jprb, &
      &  7.7039e-24_jprb, 7.5007e-24_jprb, 7.2674e-24_jprb, 7.024e-24_jprb, 6.7952e-24_jprb, 6.5954e-24_jprb, &
      &  6.4246e-24_jprb, 6.2745e-24_jprb, 6.1379e-24_jprb, 6.012e-24_jprb, 5.8965e-24_jprb, 5.7893e-24_jprb, &
      &  5.6846e-24_jprb, 5.5762e-24_jprb, 5.4626e-24_jprb, 5.3503e-24_jprb, 5.2493e-24_jprb, 5.1645e-24_jprb, &
      &  5.0935e-24_jprb, 5.0325e-24_jprb, 4.9831e-24_jprb, 4.9502e-24_jprb, 4.9326e-24_jprb, 4.9196e-24_jprb, &
      &  4.8985e-24_jprb, 4.8649e-24_jprb, 4.8235e-24_jprb, 4.7815e-24_jprb, 4.7443e-24_jprb, 4.7188e-24_jprb, &
      &  4.7184e-24_jprb, 4.7588e-24_jprb, 4.8476e-24_jprb, 4.9732e-24_jprb, 5.1057e-24_jprb, 5.2126e-24_jprb, &
      &  5.2794e-24_jprb, 5.3176e-24_jprb, 5.352e-24_jprb, 5.3985e-24_jprb, 5.4552e-24_jprb, 5.512e-24_jprb, &
      &  5.5662e-24_jprb, 5.626e-24_jprb, 5.7043e-24_jprb, 5.8132e-24_jprb, 5.9634e-24_jprb, 6.1644e-24_jprb, &
      &  6.423e-24_jprb, 6.7393e-24_jprb, 7.1025e-24_jprb, 7.4901e-24_jprb, 7.8717e-24_jprb, 8.2203e-24_jprb, &
      &  8.5216e-24_jprb, 8.7774e-24_jprb, 9.0045e-24_jprb, 9.2358e-24_jprb, 9.5182e-24_jprb, 9.8952e-24_jprb, &
      &  1.0391e-23_jprb, 1.1028e-23_jprb, 1.1877e-23_jprb, 1.3085e-23_jprb, 1.4868e-23_jprb, 1.7426e-23_jprb, &
      &  2.0818e-23_jprb, 2.4845e-23_jprb, 2.9048e-23_jprb, 3.2903e-23_jprb, 3.6103e-23_jprb, 3.871e-23_jprb, &
      &  4.1068e-23_jprb, 4.3612e-23_jprb, 4.6634e-23_jprb, 5.0014e-23_jprb, 5.3208e-23_jprb, 5.5834e-23_jprb, &
      &  5.8393e-23_jprb, 6.2322e-23_jprb, 6.926e-23_jprb, 8.0144e-23_jprb, 9.4568e-23_jprb, 1.1048e-22_jprb, &
      &  1.2453e-22_jprb, 1.3327e-22_jprb, 1.3484e-22_jprb, 1.2996e-22_jprb, 1.213e-22_jprb, 1.1176e-22_jprb, &
      &  1.0312e-22_jprb, 9.5889e-23_jprb, 8.9985e-23_jprb, 8.5239e-23_jprb, 8.1363e-23_jprb, 7.8021e-23_jprb, &
      &  7.5221e-23_jprb, 7.361e-23_jprb, 7.4255e-23_jprb, 7.8068e-23_jprb, 8.5399e-23_jprb, 9.5936e-23_jprb, &
      &  1.0867e-22_jprb, 1.217e-22_jprb, 1.3237e-22_jprb, 1.3862e-22_jprb, 1.4106e-22_jprb, 1.433e-22_jprb, &
      &  1.5026e-22_jprb, 1.6658e-22_jprb, 1.9713e-22_jprb, 2.4742e-22_jprb, 3.1995e-22_jprb, 4.0804e-22_jprb, &
      &  4.9535e-22_jprb, 5.6626e-22_jprb, 6.1883e-22_jprb, 6.6746e-22_jprb, 7.3407e-22_jprb, 8.4108e-22_jprb, &
      &  1.0189e-21_jprb, 1.3242e-21_jprb, 1.8459e-21_jprb, 2.6722e-21_jprb, 3.8123e-21_jprb, 5.1268e-21_jprb, &
      &  6.3343e-21_jprb, 7.1195e-21_jprb, 7.2867e-21_jprb, 6.8681e-21_jprb, 6.1129e-21_jprb, 5.3543e-21_jprb, &
      &  4.8525e-21_jprb, 4.7192e-21_jprb, 4.9458e-21_jprb, 5.4583e-21_jprb, 6.1302e-21_jprb, 6.7719e-21_jprb, &
      &  7.1673e-21_jprb, 7.1709e-21_jprb, 6.7871e-21_jprb, 6.1566e-21_jprb, 5.4688e-21_jprb, 4.8745e-21_jprb, &
      &  4.445e-21_jprb, 4.1768e-21_jprb, 4.019e-21_jprb, 3.9046e-21_jprb, 3.7776e-21_jprb, 3.6045e-21_jprb, &
      &  3.3712e-21_jprb, 3.0759e-21_jprb, 2.7291e-21_jprb, 2.3517e-21_jprb, 1.968e-21_jprb, 1.6023e-21_jprb, &
      &  1.2803e-21_jprb, 1.0289e-21_jprb, 8.6711e-22_jprb, 7.9363e-22_jprb, 7.8116e-22_jprb, 7.8693e-22_jprb, &
      &  7.7527e-22_jprb, 7.3373e-22_jprb, 6.703e-22_jprb, 5.9934e-22_jprb, 5.3098e-22_jprb, 4.6931e-22_jprb, &
      &  4.1445e-22_jprb, 3.6464e-22_jprb, 3.1707e-22_jprb, 2.6865e-22_jprb, 2.1813e-22_jprb, 1.6812e-22_jprb, &
      &  1.2442e-22_jprb, 9.2116e-23_jprb, 7.2181e-23_jprb, 6.1583e-23_jprb, 5.6019e-23_jprb, 5.2306e-23_jprb, &
      &  4.909e-23_jprb, 4.628e-23_jprb, 4.4098e-23_jprb, 4.2463e-23_jprb, 4.1028e-23_jprb, 3.9488e-23_jprb, &
      &  3.7792e-23_jprb, 3.6089e-23_jprb, 3.4516e-23_jprb, 3.3053e-23_jprb, 3.1567e-23_jprb, 2.9965e-23_jprb, &
      &  2.8286e-23_jprb, 2.6653e-23_jprb, 2.5179e-23_jprb, 2.3899e-23_jprb, 2.2775e-23_jprb, 2.1744e-23_jprb, &
      &  2.0773e-23_jprb, 1.9876e-23_jprb, 1.9066e-23_jprb, 1.8312e-23_jprb, 1.7569e-23_jprb, 1.6823e-23_jprb, &
      &  1.6104e-23_jprb, 1.544e-23_jprb, 1.4824e-23_jprb, 1.4234e-23_jprb, 1.3655e-23_jprb, 1.3096e-23_jprb, &
      &  1.2574e-23_jprb, 1.2097e-23_jprb, 1.1659e-23_jprb, 1.1251e-23_jprb, 1.0864e-23_jprb, 1.0497e-23_jprb, &
      &  1.015e-23_jprb, 9.8244e-24_jprb, 9.5229e-24_jprb, 9.2467e-24_jprb, 8.996e-24_jprb, 8.7673e-24_jprb, &
      &  8.5579e-24_jprb, 8.3712e-24_jprb, 8.2139e-24_jprb, 8.0841e-24_jprb, 7.967e-24_jprb, 7.8438e-24_jprb, &
      &  7.7061e-24_jprb, 7.5646e-24_jprb, 7.4417e-24_jprb, 7.3555e-24_jprb, 7.3099e-24_jprb, 7.2989e-24_jprb, &
      &  7.3188e-24_jprb, 7.3674e-24_jprb, 7.432e-24_jprb, 7.4806e-24_jprb, 7.4719e-24_jprb, 7.3755e-24_jprb, &
      &  7.1878e-24_jprb, 6.9332e-24_jprb, 6.6538e-24_jprb, 6.3908e-24_jprb, 6.1659e-24_jprb, 5.978e-24_jprb, &
      &  5.8156e-24_jprb, 5.6712e-24_jprb, 5.5455e-24_jprb, 5.4421e-24_jprb, 5.3644e-24_jprb, 5.3155e-24_jprb, &
      &  5.3e-24_jprb, 5.3236e-24_jprb, 5.387e-24_jprb, 5.479e-24_jprb, 5.5793e-24_jprb, 5.6766e-24_jprb, &
      &  5.7872e-24_jprb, 5.9573e-24_jprb, 6.2443e-24_jprb, 6.6884e-24_jprb, 7.2827e-24_jprb, 7.9588e-24_jprb, &
      &  8.6087e-24_jprb, 9.1645e-24_jprb, 9.7027e-24_jprb, 1.0508e-23_jprb, 1.202e-23_jprb, 1.4624e-23_jprb, &
      &  1.8377e-23_jprb, 2.2935e-23_jprb, 2.7829e-23_jprb, 3.2813e-23_jprb, 3.7973e-23_jprb, 4.3569e-23_jprb, &
      &  4.9865e-23_jprb, 5.7059e-23_jprb, 6.5321e-23_jprb, 7.4831e-23_jprb, 8.5818e-23_jprb, 9.8555e-23_jprb, &
      &  1.1335e-22_jprb, 1.305e-22_jprb, 1.5022e-22_jprb, 1.7253e-22_jprb, 1.9724e-22_jprb, 2.2401e-22_jprb, &
      &  2.5255e-22_jprb, 2.8264e-22_jprb, 3.1396e-22_jprb, 3.4594e-22_jprb, 3.778e-22_jprb, 4.088e-22_jprb, &
      &  4.3853e-22_jprb, 4.6682e-22_jprb, 4.9341e-22_jprb, 5.176e-22_jprb, 5.3839e-22_jprb, 5.5465e-22_jprb, &
      &  5.6523e-22_jprb, 5.6912e-22_jprb, 5.6578e-22_jprb, 5.5532e-22_jprb, 5.3835e-22_jprb, 5.1576e-22_jprb, &
      &  4.8857e-22_jprb, 4.578e-22_jprb, 4.2446e-22_jprb, 3.8955e-22_jprb, 3.541e-22_jprb, 3.1909e-22_jprb, &
      &  2.8537e-22_jprb, 2.536e-22_jprb, 2.242e-22_jprb, 1.9736e-22_jprb, 1.7312e-22_jprb, 1.5143e-22_jprb, &
      &  1.322e-22_jprb, 1.1533e-22_jprb, 1.0066e-22_jprb, 8.7981e-23_jprb, 7.7042e-23_jprb, 6.759e-23_jprb, &
      &  5.9396e-23_jprb, 5.2248e-23_jprb, 4.5915e-23_jprb, 4.0105e-23_jprb, 3.4479e-23_jprb, 2.8788e-23_jprb, &
      &  2.3098e-23_jprb, 1.7861e-23_jprb, 1.3668e-23_jprb, 1.0845e-23_jprb, 9.269e-24_jprb, 8.5449e-24_jprb, &
      &  8.2758e-24_jprb, 8.1894e-24_jprb, 8.1242e-24_jprb, 7.991e-24_jprb, 7.7539e-24_jprb, 7.4036e-24_jprb, &
      &  6.9253e-24_jprb, 6.307e-24_jprb, 5.5898e-24_jprb, 4.8798e-24_jprb, 4.2816e-24_jprb, 3.828e-24_jprb, &
      &  3.4894e-24_jprb, 3.2295e-24_jprb, 3.032e-24_jprb, 2.8892e-24_jprb, 2.7913e-24_jprb, 2.7263e-24_jprb, &
      &  2.6792e-24_jprb, 2.6313e-24_jprb, 2.5682e-24_jprb, 2.4886e-24_jprb, 2.4007e-24_jprb, 2.3106e-24_jprb, &
      &  2.22e-24_jprb, 2.1351e-24_jprb, 2.0711e-24_jprb, 2.0448e-24_jprb, 2.0615e-24_jprb, 2.1091e-24_jprb, &
      &  2.1616e-24_jprb, 2.192e-24_jprb, 2.1861e-24_jprb, 2.1476e-24_jprb, 2.0915e-24_jprb, 2.0328e-24_jprb, &
      &  1.9807e-24_jprb, 1.9374e-24_jprb, 1.9014e-24_jprb, 1.8698e-24_jprb, 1.8406e-24_jprb, 1.8126e-24_jprb, &
      &  1.7853e-24_jprb, 1.7582e-24_jprb, 1.7305e-24_jprb, 1.7012e-24_jprb, 1.6698e-24_jprb, 1.6369e-24_jprb, &
      &  1.6041e-24_jprb, 1.5729e-24_jprb, 1.544e-24_jprb, 1.517e-24_jprb, 1.4909e-24_jprb, 1.4654e-24_jprb, &
      &  1.4417e-24_jprb, 1.4218e-24_jprb, 1.4075e-24_jprb, 1.3992e-24_jprb, 1.3955e-24_jprb, 1.3933e-24_jprb, &
      &  1.3895e-24_jprb, 1.3819e-24_jprb, 1.3703e-24_jprb, 1.3559e-24_jprb, 1.3401e-24_jprb, 1.3241e-24_jprb, &
      &  1.3087e-24_jprb, 1.2944e-24_jprb, 1.2815e-24_jprb, 1.2706e-24_jprb, 1.2626e-24_jprb, 1.2592e-24_jprb, &
      &  1.2621e-24_jprb, 1.2738e-24_jprb, 1.2964e-24_jprb, 1.3323e-24_jprb, 1.3823e-24_jprb, 1.4463e-24_jprb, &
      &  1.5234e-24_jprb, 1.6126e-24_jprb, 1.7135e-24_jprb, 1.826e-24_jprb, 1.9511e-24_jprb, 2.0894e-24_jprb, &
      &  2.2413e-24_jprb, 2.4068e-24_jprb, 2.5858e-24_jprb, 2.7786e-24_jprb, 2.986e-24_jprb, 3.2089e-24_jprb, &
      &  3.4489e-24_jprb, 3.7091e-24_jprb, 3.9972e-24_jprb, 4.3322e-24_jprb, 4.7542e-24_jprb, 5.3303e-24_jprb, &
      &  6.151e-24_jprb, 7.3119e-24_jprb, 8.8861e-24_jprb, 1.0895e-23_jprb, 1.329e-23_jprb, 1.596e-23_jprb, &
      &  1.8759e-23_jprb, 2.1541e-23_jprb, 2.4172e-23_jprb, 2.6528e-23_jprb, 2.8491e-23_jprb, 2.995e-23_jprb, &
      &  3.0827e-23_jprb, 3.1105e-23_jprb, 3.087e-23_jprb, 3.0312e-23_jprb, 2.9702e-23_jprb, 2.9343e-23_jprb, &
      &  2.9526e-23_jprb, 3.0466e-23_jprb, 3.2244e-23_jprb, 3.4774e-23_jprb, 3.781e-23_jprb, 4.0996e-23_jprb, &
      &  4.3961e-23_jprb, 4.6423e-23_jprb, 4.8272e-23_jprb, 4.9584e-23_jprb, 5.0569e-23_jprb, 5.1499e-23_jprb, &
      &  5.2645e-23_jprb, 5.4239e-23_jprb, 5.6467e-23_jprb, 5.9489e-23_jprb, 6.3479e-23_jprb, 6.8664e-23_jprb, &
      &  7.5288e-23_jprb, 8.3516e-23_jprb, 9.3269e-23_jprb, 1.0403e-22_jprb, 1.1481e-22_jprb, 1.2444e-22_jprb, &
      &  1.3218e-22_jprb, 1.383e-22_jprb, 1.4412e-22_jprb, 1.5134e-22_jprb, 1.6125e-22_jprb, 1.7415e-22_jprb, &
      &  1.8952e-22_jprb, 2.066e-22_jprb, 2.2496e-22_jprb, 2.4438e-22_jprb, 2.6452e-22_jprb, 2.8476e-22_jprb, &
      &  3.0421e-22_jprb, 3.2198e-22_jprb, 3.3752e-22_jprb, 3.507e-22_jprb, 3.6156e-22_jprb, 3.6989e-22_jprb, &
      &  3.7521e-22_jprb, 3.7714e-22_jprb, 3.7551e-22_jprb, 3.7029e-22_jprb, 3.6137e-22_jprb, 3.4861e-22_jprb, &
      &  3.3207e-22_jprb, 3.1195e-22_jprb, 2.8839e-22_jprb, 2.6145e-22_jprb, 2.3153e-22_jprb, 1.9994e-22_jprb, &
      &  1.6892e-22_jprb, 1.4067e-22_jprb, 1.1627e-22_jprb, 9.5566e-23_jprb, 7.7963e-23_jprb, 6.3124e-23_jprb, &
      &  5.0965e-23_jprb, 4.1384e-23_jprb, 3.4094e-23_jprb, 2.8657e-23_jprb, 2.4597e-23_jprb, 2.1508e-23_jprb, &
      &  1.9094e-23_jprb, 1.7158e-23_jprb, 1.5579e-23_jprb, 1.427e-23_jprb, 1.316e-23_jprb, 1.2168e-23_jprb, &
      &  1.1217e-23_jprb, 1.025e-23_jprb, 9.254e-24_jprb, 8.2653e-24_jprb, 7.3457e-24_jprb, 6.5512e-24_jprb, &
      &  5.9068e-24_jprb, 5.4008e-24_jprb, 4.9982e-24_jprb, 4.6615e-24_jprb, 4.3653e-24_jprb, 4.0984e-24_jprb, &
      &  3.859e-24_jprb, 3.6483e-24_jprb, 3.466e-24_jprb, 3.309e-24_jprb, 3.1716e-24_jprb, 3.0474e-24_jprb, &
      &  2.9316e-24_jprb, 2.8223e-24_jprb, 2.7201e-24_jprb, 2.6264e-24_jprb, 2.5411e-24_jprb, 2.4612e-24_jprb, &
      &  2.3815e-24_jprb, 2.2973e-24_jprb, 2.2069e-24_jprb, 2.1121e-24_jprb, 2.0171e-24_jprb, 1.9258e-24_jprb, &
      &  1.8414e-24_jprb, 1.7654e-24_jprb, 1.6979e-24_jprb, 1.6382e-24_jprb, 1.5852e-24_jprb, 1.5381e-24_jprb, &
      &  1.4962e-24_jprb, 1.4583e-24_jprb, 1.4231e-24_jprb, 1.3898e-24_jprb, 1.3576e-24_jprb, 1.3272e-24_jprb, &
      &  1.2993e-24_jprb, 1.2748e-24_jprb, 1.2539e-24_jprb, 1.2365e-24_jprb, 1.2228e-24_jprb, 1.2134e-24_jprb, &
      &  1.2094e-24_jprb, 1.2115e-24_jprb, 1.2201e-24_jprb, 1.2354e-24_jprb, 1.2567e-24_jprb, 1.2831e-24_jprb, &
      &  1.3133e-24_jprb, 1.3459e-24_jprb, 1.3798e-24_jprb, 1.4141e-24_jprb, 1.4486e-24_jprb, 1.4835e-24_jprb, &
      &  1.5195e-24_jprb, 1.5586e-24_jprb, 1.6028e-24_jprb, 1.6537e-24_jprb, 1.7117e-24_jprb, 1.7753e-24_jprb, &
      &  1.8413e-24_jprb, 1.9059e-24_jprb, 1.9665e-24_jprb, 2.0218e-24_jprb, 2.0726e-24_jprb, 2.1204e-24_jprb, &
      &  2.1668e-24_jprb, 2.2128e-24_jprb, 2.259e-24_jprb, 2.3058e-24_jprb, 2.3534e-24_jprb, 2.4021e-24_jprb, &
      &  2.4523e-24_jprb, 2.5053e-24_jprb, 2.5634e-24_jprb, 2.6306e-24_jprb, 2.712e-24_jprb, 2.813e-24_jprb, &
      &  2.9377e-24_jprb, 3.0877e-24_jprb, 3.2621e-24_jprb, 3.4573e-24_jprb, 3.6685e-24_jprb, 3.8913e-24_jprb, &
      &  4.1229e-24_jprb, 4.3632e-24_jprb, 4.6148e-24_jprb, 4.8819e-24_jprb, 5.1688e-24_jprb, 5.4773e-24_jprb, &
      &  5.8049e-24_jprb, 6.1428e-24_jprb, 6.4778e-24_jprb, 6.7946e-24_jprb, 7.0796e-24_jprb, 7.3237e-24_jprb, &
      &  7.5217e-24_jprb, 7.672e-24_jprb, 7.7766e-24_jprb, 7.8416e-24_jprb, 7.8776e-24_jprb, 7.8996e-24_jprb, &
      &  7.9285e-24_jprb, 7.9987e-24_jprb, 8.1649e-24_jprb, 8.4949e-24_jprb, 9.0419e-24_jprb, 9.8145e-24_jprb, &
      &  1.0776e-23_jprb, 1.1875e-23_jprb, 1.3068e-23_jprb, 1.4314e-23_jprb, 1.5564e-23_jprb, 1.6771e-23_jprb, &
      &  1.7904e-23_jprb, 1.8948e-23_jprb, 1.9884e-23_jprb, 2.0691e-23_jprb, 2.1354e-23_jprb, 2.1862e-23_jprb, &
      &  2.2182e-23_jprb, 2.2268e-23_jprb, 2.2093e-23_jprb, 2.1672e-23_jprb, 2.1045e-23_jprb, 2.0243e-23_jprb, &
      &  1.9278e-23_jprb, 1.8154e-23_jprb, 1.6881e-23_jprb, 1.5486e-23_jprb, 1.4016e-23_jprb, 1.2541e-23_jprb, &
      &  1.1129e-23_jprb, 9.832e-24_jprb, 8.6691e-24_jprb, 7.6346e-24_jprb, 6.7141e-24_jprb, 5.8975e-24_jprb, &
      &  4.968e-24_jprb, 4.3671e-24_jprb, 3.8624e-24_jprb, 3.4278e-24_jprb, 3.0528e-24_jprb, 2.7114e-24_jprb, &
      &  2.4263e-24_jprb, 2.1768e-24_jprb, 1.9377e-24_jprb, 1.7307e-24_jprb, 1.5515e-24_jprb, 1.3875e-24_jprb, &
      &  1.2312e-24_jprb, 1.0993e-24_jprb, 9.7798e-25_jprb, 8.6083e-25_jprb, 7.764e-25_jprb, 6.85e-25_jprb, &
      &  6.1313e-25_jprb, 5.4488e-25_jprb, 4.8456e-25_jprb, 4.3388e-25_jprb, 3.8466e-25_jprb, 3.4714e-25_jprb, &
      &  3.0381e-25_jprb, 2.7093e-25_jprb, 2.4159e-25_jprb, 2.146e-25_jprb, 1.895e-25_jprb, 1.7131e-25_jprb, &
      &  1.5252e-25_jprb, 1.3686e-25_jprb, 1.2209e-25_jprb, 1.0972e-25_jprb, 9.7977e-26_jprb, 8.7675e-26_jprb, &
      &  7.9204e-26_jprb, 7.1362e-26_jprb, 6.3953e-26_jprb, 5.7457e-26_jprb, 5.2245e-26_jprb, 4.848e-26_jprb, &
      &  4.5517e-26_jprb, 4.2869e-26_jprb, 4.028e-26_jprb, 3.7885e-26_jprb, 3.5664e-26_jprb, 3.3581e-26_jprb, &
      &  3.1663e-26_jprb, 2.9884e-26_jprb, 2.8225e-26_jprb, 2.6676e-26_jprb, 2.5247e-26_jprb, 2.3911e-26_jprb, &
      &  2.2677e-26_jprb, 2.1535e-26_jprb, 2.0478e-26_jprb, 1.9504e-26_jprb, 1.8614e-26_jprb, 1.7789e-26_jprb, &
      &  1.705e-26_jprb, 1.6366e-26_jprb, 1.5777e-26_jprb, 1.5273e-26_jprb, 1.4807e-26_jprb, 1.4417e-26_jprb, &
      &  1.4162e-26_jprb, 1.3983e-26_jprb, 1.3939e-26_jprb, 1.3944e-26_jprb, 1.4113e-26_jprb, 1.4186e-26_jprb, &
      &  1.4376e-26_jprb, 1.4731e-26_jprb, 1.498e-26_jprb, 1.5502e-26_jprb, 1.6124e-26_jprb, 1.7098e-26_jprb, &
      &  1.7712e-26_jprb, 1.8405e-26_jprb, 1.8943e-26_jprb, 1.9885e-26_jprb, 2.0406e-26_jprb, 2.1331e-26_jprb, &
      &  2.213e-26_jprb, 2.2882e-26_jprb, 2.3743e-26_jprb, 2.429e-26_jprb, 2.5056e-26_jprb, 2.5526e-26_jprb, &
      &  2.5888e-26_jprb, 2.6965e-26_jprb, 2.7597e-26_jprb, 2.8459e-26_jprb, 2.9204e-26_jprb, 2.9432e-26_jprb, &
      &  2.999e-26_jprb, 3.1455e-26_jprb, 3.2305e-26_jprb, 3.5854e-26_jprb, 3.838e-26_jprb, 4.0871e-26_jprb, &
      &  4.3627e-26_jprb, 4.7592e-26_jprb, 5.1686e-26_jprb, 5.6762e-26_jprb, 6.2312e-26_jprb, 6.7952e-26_jprb, &
      &  7.6265e-26_jprb, 8.4968e-26_jprb, 9.5213e-26_jprb, 1.0535e-25_jprb, 1.179e-25_jprb, 1.3584e-25_jprb, &
      &  1.5393e-25_jprb, 1.7113e-25_jprb, 2.0686e-25_jprb, 2.4022e-25_jprb, 2.8616e-25_jprb, 3.3586e-25_jprb, &
      &  4.1514e-25_jprb, 4.7936e-25_jprb, 6.2873e-25_jprb, 7.7462e-25_jprb, 9.7864e-25_jprb, 1.1634e-24_jprb, &
      &  1.3123e-24_jprb, 1.4008e-24_jprb, 1.4762e-24_jprb, 1.5847e-24_jprb, 1.7222e-24_jprb, 1.8486e-24_jprb, &
      &  1.967e-24_jprb, 2.1032e-24_jprb, 2.2169e-24_jprb, 2.3722e-24_jprb, 2.4564e-24_jprb, 2.5852e-24_jprb, &
      &  2.705e-24_jprb, 2.8115e-24_jprb, 2.9806e-24_jprb, 3.1605e-24_jprb, 3.2677e-24_jprb, 3.5514e-24_jprb, &
      &  3.8263e-24_jprb, 4.0097e-24_jprb, 4.4235e-24_jprb, 4.9712e-24_jprb, 5.5754e-24_jprb, 6.0454e-24_jprb, &
      &  7.0108e-24_jprb, 7.584e-24_jprb, 8.1172e-24_jprb, 9.4328e-24_jprb, 9.9276e-24_jprb, 1.093e-23_jprb, &
      &  1.17e-23_jprb, 1.2568e-23_jprb, 1.3183e-23_jprb, 1.4275e-23_jprb, 1.3886e-23_jprb, 1.444e-23_jprb, &
      &  1.4486e-23_jprb, 1.4489e-23_jprb, 1.429e-23_jprb, 1.377e-23_jprb, 1.3368e-23_jprb, 1.3177e-23_jprb, &
      &  1.2367e-23_jprb, 1.1834e-23_jprb, 1.0658e-23_jprb, 9.6508e-24_jprb, 8.4795e-24_jprb, 7.4209e-24_jprb, &
      &  6.807e-24_jprb, 5.9813e-24_jprb, 5.1939e-24_jprb, 4.5212e-24_jprb, 3.9832e-24_jprb, 3.5348e-24_jprb, &
      &  3.1633e-24_jprb, 2.8463e-24_jprb, 2.5871e-24_jprb, 2.3739e-24_jprb, 2.1908e-24_jprb, 2.0506e-24_jprb, &
      &  1.9244e-24_jprb, 1.812e-24_jprb, 1.7441e-24_jprb, 1.6847e-24_jprb, 1.6753e-24_jprb, 1.6278e-24_jprb, &
      &  1.6544e-24_jprb, 1.6165e-24_jprb, 1.5863e-24_jprb, 1.6567e-24_jprb, 1.6385e-24_jprb, 1.6773e-24_jprb, &
      &  1.7184e-24_jprb, 1.7188e-24_jprb, 1.7083e-24_jprb, 1.7275e-24_jprb, 1.7456e-24_jprb, 1.7317e-24_jprb, &
      &  1.698e-24_jprb, 1.6918e-24_jprb, 1.4887e-24_jprb, 1.3273e-24_jprb, 1.1568e-24_jprb, 9.8612e-25_jprb, &
      &  7.8544e-25_jprb, 6.5371e-25_jprb, 5.6034e-25_jprb, 4.5956e-25_jprb, 4.0806e-25_jprb, 3.3372e-25_jprb, &
      &  2.9262e-25_jprb, 2.4931e-25_jprb, 2.1063e-25_jprb, 1.7557e-25_jprb, 1.5098e-25_jprb, 1.2824e-25_jprb, &
      &  1.1321e-25_jprb, 9.9195e-26_jprb, 8.8875e-26_jprb, 7.9639e-26_jprb, 7.3348e-26_jprb, 6.5625e-26_jprb, &
      &  5.9919e-26_jprb, 5.5129e-26_jprb, 5.0375e-26_jprb, 4.6303e-26_jprb, 4.2269e-26_jprb, 3.9179e-26_jprb, &
      &  3.6345e-26_jprb, 3.3619e-26_jprb, 3.1047e-26_jprb, 2.8789e-26_jprb, 2.6754e-26_jprb, 2.491e-26_jprb, &
      &  2.3188e-26_jprb, 2.1553e-26_jprb, 2.0145e-26_jprb, 1.8837e-26_jprb, 1.7629e-26_jprb, 1.6522e-26_jprb, &
      &  1.5528e-26_jprb, 1.4613e-26_jprb, 1.3753e-26_jprb, 1.3018e-26_jprb, 1.2316e-26_jprb, 1.1683e-26_jprb, &
      &  1.1114e-26_jprb, 1.0602e-26_jprb, 1.0145e-26_jprb, 9.7424e-27_jprb, 9.3925e-27_jprb, 9.0929e-27_jprb, &
      &  8.8251e-27_jprb, 8.6705e-27_jprb, 8.5272e-27_jprb, 8.3917e-27_jprb, 8.4768e-27_jprb, 8.5353e-27_jprb, &
      &  8.8324e-27_jprb, 8.9717e-27_jprb, 9.1951e-27_jprb, 9.9323e-27_jprb, 1.0583e-26_jprb, 1.1603e-26_jprb, &
      &  1.2671e-26_jprb, 1.4532e-26_jprb, 1.7217e-26_jprb, 2.0363e-26_jprb, 2.3081e-26_jprb, 2.5051e-26_jprb, &
      &  2.6813e-26_jprb, 2.9107e-26_jprb, 3.2018e-26_jprb, 3.445e-26_jprb, 3.691e-26_jprb, 3.9693e-26_jprb, &
      &  4.3037e-26_jprb, 4.6255e-26_jprb, 4.9681e-26_jprb, 5.3337e-26_jprb, 5.8494e-26_jprb, 6.2412e-26_jprb, &
      &  6.8369e-26_jprb, 7.4348e-26_jprb, 7.9302e-26_jprb, 8.7467e-26_jprb, 9.5181e-26_jprb, 1.059e-25_jprb, &
      &  1.2e-25_jprb, 1.336e-25_jprb, 1.4962e-25_jprb, 1.7093e-25_jprb, 1.9288e-25_jprb, 2.2637e-25_jprb, &
      &  2.4924e-25_jprb, 2.7995e-25_jprb, 3.3067e-25_jprb, 3.6824e-25_jprb, 4.2405e-25_jprb, 4.585e-25_jprb, &
      &  5.0423e-25_jprb, 5.522e-25_jprb, 6.2417e-25_jprb, 6.8115e-25_jprb, 7.2685e-25_jprb, 7.8727e-25_jprb, &
      &  8.2767e-25_jprb, 8.7469e-25_jprb, 8.8161e-25_jprb, 9.3506e-25_jprb, 9.5237e-25_jprb, 9.4914e-25_jprb, &
      &  9.3189e-25_jprb, 8.9036e-25_jprb, 8.723e-25_jprb, 8.5189e-25_jprb, 8.0689e-25_jprb, 7.5936e-25_jprb, &
      &  6.9646e-25_jprb, 6.6542e-25_jprb, 5.8407e-25_jprb, 5.1399e-25_jprb, 4.4699e-25_jprb, 3.9365e-25_jprb, &
      &  3.503e-25_jprb, 3.0898e-25_jprb, 2.703e-25_jprb, 2.3994e-25_jprb, 2.1595e-25_jprb, 1.9352e-25_jprb, &
      &  1.7625e-25_jprb, 1.6106e-25_jprb, 1.4819e-25_jprb, 1.3741e-25_jprb, 1.2834e-25_jprb, 1.2061e-25_jprb, &
      &  1.1435e-25_jprb, 1.1134e-25_jprb, 1.0671e-25_jprb, 1.0286e-25_jprb, 1.0114e-25_jprb, 9.9598e-26_jprb, &
      &  9.828e-26_jprb, 9.6685e-26_jprb, 9.7025e-26_jprb, 9.4288e-26_jprb, 9.2948e-26_jprb, 9.3649e-26_jprb, &
      &  9.1467e-26_jprb, 9.0561e-26_jprb, 9.0846e-26_jprb, 8.975e-26_jprb, 8.9129e-26_jprb, 8.1004e-26_jprb, &
      &  7.0371e-26_jprb, 6.3239e-26_jprb, 4.9925e-26_jprb, 4.2146e-26_jprb, 3.2169e-26_jprb, 2.7084e-26_jprb, &
      &  2.2294e-26_jprb, 2.1349e-26_jprb, 1.9109e-26_jprb, 1.5748e-26_jprb, 1.3729e-26_jprb, 1.2297e-26_jprb, &
      &  1.013e-26_jprb, 8.8219e-27_jprb, 7.3765e-27_jprb, 6.463e-27_jprb, 5.6808e-27_jprb, 5.0697e-27_jprb, &
      &  4.5382e-27_jprb, 4.1286e-27_jprb, 3.7619e-27_jprb, 3.4356e-27_jprb, 3.1678e-27_jprb, 2.9227e-27_jprb, &
      &  2.7041e-27_jprb, 2.5058e-27_jprb, 2.3289e-27_jprb, 2.1647e-27_jprb, 2.0131e-27_jprb, 1.8794e-27_jprb, &
      &  1.7583e-27_jprb, 1.6499e-27_jprb, 1.5504e-27_jprb, 1.4598e-27_jprb, 1.3794e-27_jprb, 1.3041e-27_jprb, &
      &  1.237e-27_jprb, 1.177e-27_jprb, 1.1235e-27_jprb, 1.0762e-27_jprb, 1.0345e-27_jprb, 9.9839e-28_jprb, &
      &  9.6752e-28_jprb, 9.4168e-28_jprb, 9.2075e-28_jprb, 9.0489e-28_jprb, 8.937e-28_jprb, 8.8745e-28_jprb, &
      &  8.8617e-28_jprb, 8.8986e-28_jprb, 8.9853e-28_jprb, 9.1246e-28_jprb, 9.3178e-28_jprb, 9.5678e-28_jprb, &
      &  9.8773e-28_jprb, 1.025e-27_jprb, 1.069e-27_jprb, 1.12e-27_jprb, 1.179e-27_jprb, 1.2465e-27_jprb, &
      &  1.324e-27_jprb, 1.4135e-27_jprb, 1.5177e-27_jprb, 1.6299e-27_jprb, 1.7569e-27_jprb, 1.9093e-27_jprb, &
      &  2.0354e-27_jprb, 2.1962e-27_jprb, 2.409e-27_jprb, 2.6846e-27_jprb, 3.0683e-27_jprb, 3.4899e-27_jprb, &
      &  3.9454e-27_jprb, 4.4909e-27_jprb, 5.0945e-27_jprb, 5.7792e-27_jprb, 6.5892e-27_jprb, 7.4725e-27_jprb, &
      &  8.515e-27_jprb, 9.7533e-27_jprb, 1.1036e-26_jprb, 1.2668e-26_jprb, 1.4467e-26_jprb, 1.6505e-26_jprb, &
      &  1.9044e-26_jprb, 2.1708e-26_jprb, 2.4807e-26_jprb, 2.8381e-26_jprb, 3.2433e-26_jprb, 3.6801e-26_jprb, &
      &  4.3133e-26_jprb, 4.8568e-26_jprb, 5.6149e-26_jprb, 6.4078e-26_jprb, 7.3251e-26_jprb, 8.2492e-26_jprb, &
      &  9.1788e-26_jprb, 1.0437e-25_jprb, 1.1604e-25_jprb, 1.3075e-25_jprb, 1.4239e-25_jprb, 1.5801e-25_jprb, &
      &  1.7652e-25_jprb, 1.9409e-25_jprb, 2.1497e-25_jprb, 2.3737e-25_jprb, 2.5858e-25_jprb, 2.8434e-25_jprb, &
      &  3.1956e-25_jprb, 3.5649e-25_jprb, 3.7812e-25_jprb, 4.1144e-25_jprb, 4.6168e-25_jprb, 4.8741e-25_jprb, &
      &  5.3517e-25_jprb, 5.7516e-25_jprb, 6.3489e-25_jprb, 6.6729e-25_jprb, 6.9837e-25_jprb, 7.4838e-25_jprb, &
      &  7.6466e-25_jprb, 8.0982e-25_jprb, 8.2879e-25_jprb, 8.2939e-25_jprb, 8.3884e-25_jprb, 8.2324e-25_jprb, &
      &  8.115e-25_jprb, 7.8518e-25_jprb, 7.6659e-25_jprb, 7.4645e-25_jprb, 6.8666e-25_jprb, 6.521e-25_jprb, &
      &  5.8046e-25_jprb, 5.4508e-25_jprb, 4.7717e-25_jprb, 4.1808e-25_jprb, 3.6268e-25_jprb, 3.1474e-25_jprb, &
      &  2.7205e-25_jprb, 2.3923e-25_jprb, 2.1098e-25_jprb, 1.8717e-25_jprb, 1.6684e-25_jprb, 1.4984e-25_jprb, &
      &  1.3511e-25_jprb, 1.2325e-25_jprb, 1.1258e-25_jprb, 1.044e-25_jprb, 9.6927e-26_jprb, 9.1345e-26_jprb, &
      &  8.6884e-26_jprb, 8.3559e-26_jprb, 8.0315e-26_jprb, 7.9891e-26_jprb, 7.9707e-26_jprb, 8.0669e-26_jprb, &
      &  8.1859e-26_jprb, 8.3235e-26_jprb, 8.888e-26_jprb, 9.0206e-26_jprb, 9.2655e-26_jprb, 1.0142e-25_jprb, &
      &  1.0581e-25_jprb, 1.1418e-25_jprb, 1.2083e-25_jprb, 1.3029e-25_jprb, 1.3598e-25_jprb, 1.3866e-25_jprb, &
      &  1.4537e-25_jprb, 1.489e-25_jprb, 1.5072e-25_jprb, 1.5627e-25_jprb, 1.5652e-25_jprb, 1.5362e-25_jprb, &
      &  1.5e-25_jprb, 1.3321e-25_jprb, 1.1564e-25_jprb, 9.5992e-26_jprb, 8.0959e-26_jprb, 7.0986e-26_jprb, &
      &  5.0227e-26_jprb, 4.7899e-26_jprb, 4.0751e-26_jprb, 3.505e-26_jprb, 2.9616e-26_jprb, 2.4262e-26_jprb, &
      &  2.0403e-26_jprb, 1.7943e-26_jprb, 1.5899e-26_jprb, 1.4009e-26_jprb, 1.2367e-26_jprb, 1.1387e-26_jprb, &
      &  1.0306e-26_jprb, 8.888e-27_jprb, 7.6277e-27_jprb, 6.5781e-27_jprb, 5.9099e-27_jprb, 5.5118e-27_jprb, &
      &  5.1059e-27_jprb, 4.7724e-27_jprb, 4.3158e-27_jprb, 3.8512e-27_jprb, 3.6024e-27_jprb, 3.204e-27_jprb, &
      &  2.9075e-27_jprb, 2.6238e-27_jprb, 2.4013e-27_jprb, 2.1873e-27_jprb, 2.0302e-27_jprb, 1.8671e-27_jprb, &
      &  1.714e-27_jprb, 1.5916e-27_jprb, 1.4749e-27_jprb, 1.3724e-27_jprb, 1.2819e-27_jprb, 1.2032e-27_jprb, &
      &  1.1333e-27_jprb, 1.0711e-27_jprb, 1.016e-27_jprb, 9.6712e-28_jprb, 9.2396e-28_jprb, 8.8638e-28_jprb, &
      &  8.5408e-28_jprb, 8.2677e-28_jprb, 8.0448e-28_jprb, 7.8705e-28_jprb, 7.745e-28_jprb, 7.6685e-28_jprb, &
      &  7.6453e-28_jprb, 7.6773e-28_jprb, 7.7689e-28_jprb, 7.9217e-28_jprb, 8.1418e-28_jprb, 8.4892e-28_jprb, &
      &  8.8938e-28_jprb, 9.3946e-28_jprb, 1.0257e-27_jprb, 1.1004e-27_jprb, 1.1993e-27_jprb, 1.3427e-27_jprb, &
      &  1.5465e-27_jprb, 1.7562e-27_jprb, 2.0202e-27_jprb, 2.4394e-27_jprb, 3.0065e-27_jprb, 3.773e-27_jprb, &
      &  4.524e-27_jprb, 5.0906e-27_jprb, 5.5615e-27_jprb, 6.2141e-27_jprb, 6.8992e-27_jprb, 7.5928e-27_jprb, &
      &  8.5579e-27_jprb, 9.4275e-27_jprb, 1.0418e-26_jprb, 1.2114e-26_jprb, 1.2925e-26_jprb, 1.4706e-26_jprb, &
      &  1.6789e-26_jprb, 1.8562e-26_jprb, 2.0596e-26_jprb, 2.3195e-26_jprb, 2.7016e-26_jprb, 2.9303e-26_jprb, &
      &  3.3238e-26_jprb, 3.628e-26_jprb, 4.0622e-26_jprb, 4.439e-26_jprb, 4.8193e-26_jprb, 5.1466e-26_jprb, &
      &  5.4759e-26_jprb, 5.9066e-26_jprb, 6.1985e-26_jprb, 6.4325e-26_jprb, 6.3906e-26_jprb, 6.8415e-26_jprb, &
      &  6.8598e-26_jprb, 6.7689e-26_jprb, 6.6703e-26_jprb, 6.3868e-26_jprb, 6.2878e-26_jprb, 6.0253e-26_jprb, &
      &  5.7409e-26_jprb, 5.4221e-26_jprb, 5.0998e-26_jprb, 4.3986e-26_jprb, 3.9361e-26_jprb, 3.4498e-26_jprb, &
      &  3.0867e-26_jprb, 2.6844e-26_jprb, 2.3064e-26_jprb, 2.0116e-26_jprb, 1.7755e-26_jprb, 1.5592e-26_jprb, &
      &  1.3584e-26_jprb, 1.2328e-26_jprb, 1.0945e-26_jprb, 9.8449e-27_jprb, 8.9647e-27_jprb, 8.1848e-27_jprb, &
      &  7.566e-27_jprb, 7.0276e-27_jprb, 6.6947e-27_jprb, 6.2958e-27_jprb, 6.098e-27_jprb, 5.9846e-27_jprb, &
      &  5.8897e-27_jprb, 5.9609e-27_jprb, 6.0196e-27_jprb, 6.0328e-27_jprb, 6.6631e-27_jprb, 6.8841e-27_jprb, &
      &  7.2029e-27_jprb, 7.6243e-27_jprb, 8.199e-27_jprb, 8.4798e-27_jprb, 9.2779e-27_jprb, 1.0047e-26_jprb, &
      &  1.0631e-26_jprb, 1.1318e-26_jprb, 1.1716e-26_jprb, 1.2374e-26_jprb, 1.2715e-26_jprb, 1.2854e-26_jprb, &
      &  1.2916e-26_jprb, 1.265e-26_jprb, 1.2344e-26_jprb, 1.095e-26_jprb, 9.5756e-27_jprb, 7.7861e-27_jprb, &
      &  6.166e-27_jprb, 4.8812e-27_jprb, 4.0325e-27_jprb, 3.5922e-27_jprb, 3.2773e-27_jprb, 2.5026e-27_jprb, &
      &  2.2312e-27_jprb, 1.835e-27_jprb, 1.478e-27_jprb, 1.2322e-27_jprb, 1.0298e-27_jprb, 8.8279e-28_jprb, &
      &  7.6608e-28_jprb, 6.7358e-28_jprb, 5.9685e-28_jprb, 5.3543e-28_jprb, 4.8421e-28_jprb, 4.4096e-28_jprb, &
      &  4.0393e-28_jprb, 3.72e-28_jprb, 3.4422e-28_jprb, 3.2011e-28_jprb, 2.9903e-28_jprb, 2.8082e-28_jprb, &
      &  2.6502e-28_jprb, 2.5146e-28_jprb, 2.3965e-28_jprb, 2.2994e-28_jprb, 2.2167e-28_jprb, 2.1517e-28_jprb, &
      &  2.1011e-28_jprb, 2.0652e-28_jprb, 2.0421e-28_jprb, 2.0336e-28_jprb, 2.0381e-28_jprb, 2.0572e-28_jprb, &
      &  2.0878e-28_jprb, 2.1346e-28_jprb, 2.1945e-28_jprb, 2.2692e-28_jprb, 2.3602e-28_jprb, 2.4661e-28_jprb, &
      &  2.59e-28_jprb, 2.7304e-28_jprb, 2.8922e-28_jprb, 3.0738e-28_jprb, 3.2556e-28_jprb, 3.4541e-28_jprb, &
      &  3.7248e-28_jprb, 4.1237e-28_jprb, 4.6642e-28_jprb, 5.3021e-28_jprb, 6.0114e-28_jprb, 6.8234e-28_jprb, &
      &  7.7547e-28_jprb, 8.8172e-28_jprb, 1.0031e-27_jprb, 1.1414e-27_jprb, 1.299e-27_jprb, 1.4785e-27_jprb, &
      &  1.683e-27_jprb, 1.9185e-27_jprb, 2.1872e-27_jprb, 2.496e-27_jprb, 2.8449e-27_jprb, 3.2521e-27_jprb, &
      &  3.7045e-27_jprb, 4.2419e-27_jprb, 4.853e-27_jprb, 5.5925e-27_jprb, 6.3777e-27_jprb, 7.3566e-27_jprb, &
      &  8.541e-27_jprb, 9.6421e-27_jprb, 1.1144e-26_jprb, 1.2697e-26_jprb, 1.4406e-26_jprb, 1.7037e-26_jprb, &
      &  1.8948e-26_jprb, 2.1697e-26_jprb, 2.4733e-26_jprb, 2.7939e-26_jprb, 3.1283e-26_jprb, 3.5568e-26_jprb, &
      &  3.9456e-26_jprb, 4.3265e-26_jprb, 5.023e-26_jprb, 5.5157e-26_jprb, 6.212e-26_jprb, 6.7361e-26_jprb, &
      &  7.2559e-26_jprb, 7.6099e-26_jprb, 8.4117e-26_jprb, 9.169e-26_jprb, 9.4489e-26_jprb, 9.7056e-26_jprb, &
      &  1.0057e-25_jprb, 1.0255e-25_jprb, 1.0477e-25_jprb, 1.06e-25_jprb, 1.0515e-25_jprb, 1.0378e-25_jprb, &
      &  1.0152e-25_jprb, 9.8431e-26_jprb, 9.1357e-26_jprb, 9.0682e-26_jprb, 8.3901e-26_jprb, 7.8554e-26_jprb, &
      &  6.9907e-26_jprb, 6.3067e-26_jprb, 5.7324e-26_jprb, 4.972e-26_jprb, 4.3401e-26_jprb, 3.8011e-26_jprb, &
      &  3.287e-26_jprb, 2.8627e-26_jprb, 2.4908e-26_jprb, 2.2003e-26_jprb, 1.93e-26_jprb, 1.701e-26_jprb, &
      &  1.4985e-26_jprb, 1.3277e-26_jprb, 1.1813e-26_jprb, 1.0552e-26_jprb, 9.4698e-27_jprb, 8.503e-27_jprb, &
      &  7.7031e-27_jprb, 6.9916e-27_jprb, 6.4047e-27_jprb, 5.9099e-27_jprb, 5.5126e-27_jprb, 5.1733e-27_jprb, &
      &  4.9043e-27_jprb, 4.7658e-27_jprb, 4.6341e-27_jprb, 4.5954e-27_jprb, 4.6084e-27_jprb, 4.656e-27_jprb, &
      &  4.8142e-27_jprb, 4.9467e-27_jprb, 5.2056e-27_jprb, 5.5011e-27_jprb, 5.8697e-27_jprb, 6.2163e-27_jprb, &
      &  6.6221e-27_jprb, 7.0337e-27_jprb, 7.6679e-27_jprb, 8.1691e-27_jprb, 8.9281e-27_jprb, 9.4949e-27_jprb, &
      &  1.0536e-26_jprb, 1.1026e-26_jprb, 1.1799e-26_jprb, 1.2583e-26_jprb, 1.3458e-26_jprb, 1.4189e-26_jprb, &
      &  1.4825e-26_jprb, 1.4979e-26_jprb, 1.5749e-26_jprb, 1.5771e-26_jprb, 1.5748e-26_jprb, 1.5687e-26_jprb, &
      &  1.5484e-26_jprb, 1.5039e-26_jprb, 1.4316e-26_jprb, 1.4212e-26_jprb, 1.3121e-26_jprb, 1.2489e-26_jprb, &
      &  1.1531e-26_jprb, 1.0543e-26_jprb, 9.1027e-27_jprb, 7.7907e-27_jprb, 6.9617e-27_jprb, 5.8884e-27_jprb, &
      &  5.1191e-27_jprb, 4.4849e-27_jprb, 3.9117e-27_jprb, 3.4368e-27_jprb, 3.0339e-27_jprb, 2.6588e-27_jprb, &
      &  2.37e-27_jprb, 2.0986e-27_jprb, 1.8606e-27_jprb, 1.6553e-27_jprb, 1.4737e-27_jprb, 1.2979e-27_jprb, &
      &  1.1633e-27_jprb, 1.0371e-27_jprb, 9.2065e-28_jprb, 8.1834e-28_jprb, 7.2731e-28_jprb, 6.4811e-28_jprb, &
      &  5.7862e-28_jprb, 5.1778e-28_jprb, 4.6455e-28_jprb, 4.1787e-28_jprb, 3.772e-28_jprb, 3.4203e-28_jprb, &
      &  3.1183e-28_jprb, 2.8588e-28_jprb, 2.642e-28_jprb, 2.4608e-28_jprb, 2.3171e-28_jprb, 2.2054e-28_jprb, &
      &  2.1259e-28_jprb, 2.0786e-28_jprb, 2.06e-28_jprb, 2.0737e-28_jprb, 2.1162e-28_jprb, 2.193e-28_jprb, &
      &  2.304e-28_jprb, 2.4494e-28_jprb, 2.6363e-28_jprb, 2.8686e-28_jprb, 3.157e-28_jprb, 3.5017e-28_jprb, &
      &  3.8937e-28_jprb, 4.391e-28_jprb, 4.9684e-28_jprb, 5.5211e-28_jprb, 6.3947e-28_jprb, 7.1842e-28_jprb, &
      &  8.1015e-28_jprb, 9.1393e-28_jprb, 1.063e-27_jprb, 1.2095e-27_jprb, 1.366e-27_jprb, 1.5157e-27_jprb, &
      &  1.7111e-27_jprb, 1.95e-27_jprb, 2.2077e-27_jprb, 2.4456e-27_jprb, 2.7257e-27_jprb, 3.0335e-27_jprb, &
      &  3.3233e-27_jprb, 3.7029e-27_jprb, 4.0775e-27_jprb, 4.4889e-27_jprb, 4.9521e-27_jprb, 5.3553e-27_jprb, &
      &  5.7004e-27_jprb, 6.1448e-27_jprb, 6.5805e-27_jprb, 6.9984e-27_jprb, 7.2186e-27_jprb, 7.4519e-27_jprb, &
      &  7.6176e-27_jprb, 7.8973e-27_jprb, 7.7619e-27_jprb, 7.8176e-27_jprb, 7.8881e-27_jprb, 7.6899e-27_jprb, &
      &  7.5411e-27_jprb, 7.2651e-27_jprb, 6.9556e-27_jprb, 6.6678e-27_jprb, 6.2653e-27_jprb, 5.951e-27_jprb, &
      &  5.1578e-27_jprb, 4.9276e-27_jprb, 4.2717e-27_jprb, 3.7612e-27_jprb, 3.2225e-27_jprb, 2.8444e-27_jprb, &
      &  2.4528e-27_jprb, 2.1351e-27_jprb, 1.8194e-27_jprb, 1.6063e-27_jprb, 1.4009e-27_jprb, 1.2305e-27_jprb, &
      &  1.0848e-27_jprb, 9.5153e-28_jprb, 8.3846e-28_jprb, 7.4054e-28_jprb, 6.5613e-28_jprb, 5.8039e-28_jprb, &
      &  5.2192e-28_jprb, 4.6675e-28_jprb, 4.1843e-28_jprb, 3.7623e-28_jprb, 3.4015e-28_jprb, 3.0928e-28_jprb, &
      &  2.8323e-28_jprb, 2.6128e-28_jprb, 2.4343e-28_jprb, 2.2969e-28_jprb, 2.1912e-28_jprb, 2.1097e-28_jprb, &
      &  2.1559e-28_jprb, 2.0801e-28_jprb, 2.1206e-28_jprb, 2.1989e-28_jprb, 2.3205e-28_jprb, 2.4611e-28_jprb, &
      &  2.5604e-28_jprb, 2.9086e-28_jprb, 3.1798e-28_jprb, 3.3853e-28_jprb, 3.8289e-28_jprb, 4.1747e-28_jprb, &
      &  4.6135e-28_jprb, 5.0168e-28_jprb, 5.3506e-28_jprb, 5.855e-28_jprb, 6.2123e-28_jprb, 6.7291e-28_jprb, &
      &  6.9716e-28_jprb, 7.2143e-28_jprb, 7.2959e-28_jprb, 7.7251e-28_jprb, 7.7767e-28_jprb, 7.6078e-28_jprb, &
      &  7.4254e-28_jprb, 6.8965e-28_jprb, 6.9172e-28_jprb, 6.8999e-28_jprb, 6.5624e-28_jprb, 6.0719e-28_jprb, &
      &  5.5981e-28_jprb, 4.7267e-28_jprb, 4.2249e-28_jprb, 3.5678e-28_jprb, 3.0477e-28_jprb, 2.6208e-28_jprb, &
      &  2.2815e-28_jprb, 2.0012e-28_jprb, 1.7648e-28_jprb, 1.5662e-28_jprb, 1.3982e-28_jprb, 1.2569e-28_jprb, &
      &  1.1384e-28_jprb, 1.0402e-28_jprb, 9.6004e-29_jprb, 8.9631e-29_jprb, 8.4771e-29_jprb, 8.1331e-29_jprb, &
      &  7.9274e-29_jprb, 7.8525e-29_jprb, 7.9087e-29_jprb, 8.0979e-29_jprb, 8.4244e-29_jprb, 8.8942e-29_jprb, &
      &  9.5151e-29_jprb, 1.0295e-28_jprb, 1.1252e-28_jprb, 1.2398e-28_jprb, 1.3753e-28_jprb, 1.5341e-28_jprb, &
      &  1.7186e-28_jprb, 1.932e-28_jprb, 2.1786e-28_jprb, 2.4612e-28_jprb, 2.7848e-28_jprb, 3.1593e-28_jprb, &
      &  3.5925e-28_jprb, 4.0981e-28_jprb, 4.7132e-28_jprb, 5.3971e-28_jprb, 6.1265e-28_jprb, 7.2191e-28_jprb, &
      &  8.3265e-28_jprb, 9.2107e-28_jprb, 1.0465e-27_jprb, 1.188e-27_jprb, 1.3864e-27_jprb, 1.5577e-27_jprb, &
      &  1.805e-27_jprb, 2.0509e-27_jprb, 2.3241e-27_jprb, 2.6172e-27_jprb, 3.0909e-27_jprb, 3.4749e-27_jprb, &
      &  3.9574e-27_jprb, 4.4167e-27_jprb, 5.1653e-27_jprb, 5.7201e-27_jprb, 6.5979e-27_jprb, 7.5139e-27_jprb, &
      &  8.2125e-27_jprb, 9.4729e-27_jprb, 1.035e-26_jprb, 1.1509e-26_jprb, 1.2201e-26_jprb, 1.3426e-26_jprb, &
      &  1.4051e-26_jprb, 1.4554e-26_jprb, 1.5144e-26_jprb, 1.5409e-26_jprb, 1.5719e-26_jprb, 1.5903e-26_jprb, &
      &  1.5749e-26_jprb, 1.5549e-26_jprb, 1.4901e-26_jprb, 1.4797e-26_jprb, 1.3741e-26_jprb, 1.2964e-26_jprb, &
      &  1.2307e-26_jprb, 1.0772e-26_jprb, 9.4967e-27_jprb, 8.1768e-27_jprb, 7.3212e-27_jprb, 6.3035e-27_jprb, &
      &  5.5497e-27_jprb, 4.8549e-27_jprb, 4.2751e-27_jprb, 3.7087e-27_jprb, 3.2495e-27_jprb, 2.8238e-27_jprb, &
      &  2.4715e-27_jprb, 2.1849e-27_jprb, 1.9116e-27_jprb, 0.0_jprb ]
    real(jprb), parameter :: foreign_lut(nwav_lut) = [ &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 1.12e-23_jprb, 1.58e-23_jprb, 2.26e-23_jprb, 3.34e-23_jprb, &
      &  4.42e-23_jprb, 5.965e-23_jprb, 7.51e-23_jprb, 1.0055e-22_jprb, 1.26e-22_jprb, 1.52e-22_jprb, &
      &  1.74e-22_jprb, 1.79e-22_jprb, 2.48e-22_jprb, 2.78e-22_jprb, 2.78e-22_jprb, 3.575e-22_jprb, &
      &  4.37e-22_jprb, 4.6e-22_jprb, 4.89e-22_jprb, 5.17e-22_jprb, 5.82e-22_jprb, 5.57e-22_jprb, &
      &  4.29e-22_jprb, 3.79e-22_jprb, 3.34e-22_jprb, 3.57e-22_jprb, 3.58e-22_jprb, 3.36e-22_jprb, &
      &  3.46e-22_jprb, 3.86e-22_jprb, 4.3e-22_jprb, 4.42e-22_jprb, 4.38e-22_jprb, 4.75e-22_jprb, &
      &  5.19e-22_jprb, 5.38e-22_jprb, 4.33e-22_jprb, 3.71e-22_jprb, 3.43e-22_jprb, 3.05e-22_jprb, &
      &  2.81e-22_jprb, 2.73e-22_jprb, 2.17e-22_jprb, 1.88e-22_jprb, 1.54e-22_jprb, 1.23e-22_jprb, &
      &  1.09054e-22_jprb, 9.17568e-23_jprb, 7.62162e-23_jprb, 6.62162e-23_jprb, 5.91892e-23_jprb, 4.83784e-23_jprb, &
      &  4.24324e-23_jprb, 3.48649e-23_jprb, 2.93243e-23_jprb, 2.59459e-23_jprb, 2.24324e-23_jprb, 2.01351e-23_jprb, &
      &  1.71622e-23_jprb, 1.40541e-23_jprb, 1.27027e-23_jprb, 1.15e-23_jprb, 1.02027e-23_jprb, 9.51351e-24_jprb, &
      &  8.56757e-24_jprb, 7.5e-24_jprb, 6.64865e-24_jprb, 6e-24_jprb, 5.02703e-24_jprb, 4.17568e-24_jprb, &
      &  3.63514e-24_jprb, 3.24324e-24_jprb, 2.64865e-24_jprb, 2.32432e-24_jprb, 1.90541e-24_jprb, 1.59459e-24_jprb, &
      &  1.40541e-24_jprb, 1.22027e-24_jprb, 1.1e-24_jprb, 1.02027e-24_jprb, 9.21622e-25_jprb, 8.32432e-25_jprb, &
      &  7.5e-25_jprb, 6.94595e-25_jprb, 6.17568e-25_jprb, 5.90541e-25_jprb, 5.7027e-25_jprb, 5.64865e-25_jprb, &
      &  5.45946e-25_jprb, 5.22973e-25_jprb, 5.31081e-25_jprb, 5.06757e-25_jprb, 5.09459e-25_jprb, 5.09459e-25_jprb, &
      &  4.95946e-25_jprb, 4.82432e-25_jprb, 4.94595e-25_jprb, 4.81081e-25_jprb, 4.67568e-25_jprb, 4.55405e-25_jprb, &
      &  4.41892e-25_jprb, 4.2973e-25_jprb, 4.17568e-25_jprb, 4.06757e-25_jprb, 3.95946e-25_jprb, 3.85135e-25_jprb, &
      &  3.78378e-25_jprb, 3.7973e-25_jprb, 3.78378e-25_jprb, 3.7027e-25_jprb, 3.58108e-25_jprb, 3.44595e-25_jprb, &
      &  3.36486e-25_jprb, 3.24324e-25_jprb, 3.06757e-25_jprb, 2.78378e-25_jprb, 2.66216e-25_jprb, 2.44595e-25_jprb, &
      &  2.25676e-25_jprb, 2.10811e-25_jprb, 1.91892e-25_jprb, 1.86486e-25_jprb, 1.72973e-25_jprb, 1.63514e-25_jprb, &
      &  1.60811e-25_jprb, 1.55405e-25_jprb, 1.52703e-25_jprb, 1.54054e-25_jprb, 1.44595e-25_jprb, 1.44595e-25_jprb, &
      &  1.41892e-25_jprb, 1.37838e-25_jprb, 1.40541e-25_jprb, 1.37838e-25_jprb, 1.36486e-25_jprb, 1.3e-25_jprb, &
      &  1.3e-25_jprb, 1.37838e-25_jprb, 1.37838e-25_jprb, 1.37838e-25_jprb, 1.30946e-25_jprb, 1.39189e-25_jprb, &
      &  1.43243e-25_jprb, 1.39189e-25_jprb, 1.44595e-25_jprb, 1.44595e-25_jprb, 1.51351e-25_jprb, 1.48649e-25_jprb, &
      &  1.55405e-25_jprb, 1.51351e-25_jprb, 1.51351e-25_jprb, 1.52703e-25_jprb, 1.67568e-25_jprb, 1.78378e-25_jprb, &
      &  1.82432e-25_jprb, 2.24324e-25_jprb, 2.71622e-25_jprb, 2.95946e-25_jprb, 3.7027e-25_jprb, 4.33784e-25_jprb, &
      &  4.82432e-25_jprb, 5.81081e-25_jprb, 7.17568e-25_jprb, 8.59459e-25_jprb, 1.02973e-24_jprb, 1.17027e-24_jprb, &
      &  1.32973e-24_jprb, 1.66216e-24_jprb, 2.06757e-24_jprb, 2.36486e-24_jprb, 2.67568e-24_jprb, 3.04054e-24_jprb, &
      &  3.21622e-24_jprb, 3.32432e-24_jprb, 3.56757e-24_jprb, 3.82432e-24_jprb, 3.86486e-24_jprb, 3.91892e-24_jprb, &
      &  3.63514e-24_jprb, 3.22973e-24_jprb, 3.36486e-24_jprb, 2.86486e-24_jprb, 2.81081e-24_jprb, 2.95946e-24_jprb, &
      &  3.21622e-24_jprb, 3.48649e-24_jprb, 3.85135e-24_jprb, 4.04054e-24_jprb, 4.22973e-24_jprb, 4.44595e-24_jprb, &
      &  4.66216e-24_jprb, 4.64865e-24_jprb, 4.60811e-24_jprb, 4.51351e-24_jprb, 4.41892e-24_jprb, 4.32432e-24_jprb, &
      &  4.22973e-24_jprb, 4.2027e-24_jprb, 4.55405e-24_jprb, 4.71622e-24_jprb, 5.09459e-24_jprb, 5.72973e-24_jprb, &
      &  6.52703e-24_jprb, 6.72973e-24_jprb, 7.33784e-24_jprb, 8.25676e-24_jprb, 8.40541e-24_jprb, 9.32432e-24_jprb, &
      &  1.05946e-23_jprb, 1.22973e-23_jprb, 1.37838e-23_jprb, 1.52703e-23_jprb, 1.77027e-23_jprb, 1.93243e-23_jprb, &
      &  2.36486e-23_jprb, 2.66216e-23_jprb, 2.90541e-23_jprb, 3.18919e-23_jprb, 3.52703e-23_jprb, 3.90541e-23_jprb, &
      &  4.24324e-23_jprb, 4.71622e-23_jprb, 5.25676e-23_jprb, 5.90541e-23_jprb, 6.86486e-23_jprb, 6.75901e-23_jprb, &
      &  6.65315e-23_jprb, 6.5473e-23_jprb, 6.44144e-23_jprb, 6.33559e-23_jprb, 6.22973e-23_jprb, 6.12387e-23_jprb, &
      &  6.01802e-23_jprb, 5.91216e-23_jprb, 5.80631e-23_jprb, 5.70045e-23_jprb, 5.59459e-23_jprb, 5.48874e-23_jprb, &
      &  5.38288e-23_jprb, 5.27703e-23_jprb, 5.17117e-23_jprb, 5.06532e-23_jprb, 4.95946e-23_jprb, 4.8536e-23_jprb, &
      &  4.74775e-23_jprb, 4.64189e-23_jprb, 4.53604e-23_jprb, 4.43018e-23_jprb, 4.32432e-23_jprb, 4.21847e-23_jprb, &
      &  4.11261e-23_jprb, 4.00676e-23_jprb, 3.9009e-23_jprb, 3.79505e-23_jprb, 3.68919e-23_jprb, 3.58333e-23_jprb, &
      &  3.47748e-23_jprb, 3.37162e-23_jprb, 3.26577e-23_jprb, 3.15991e-23_jprb, 3.05405e-23_jprb, 2.59459e-23_jprb, &
      &  2.05405e-23_jprb, 1.63514e-23_jprb, 1.3e-23_jprb, 1.1e-23_jprb, 9.67568e-24_jprb, 8.44595e-24_jprb, &
      &  7.21622e-24_jprb, 6.13514e-24_jprb, 5.21622e-24_jprb, 4.59459e-24_jprb, 4.06757e-24_jprb, 3.5e-24_jprb, &
      &  3.17568e-24_jprb, 3e-24_jprb, 2.67568e-24_jprb, 2.5e-24_jprb, 2.40541e-24_jprb, 2.17568e-24_jprb, &
      &  1.78378e-24_jprb, 1.7027e-24_jprb, 1.51351e-24_jprb, 1.39189e-24_jprb, 1.29054e-24_jprb, 1.22027e-24_jprb, &
      &  1.12027e-24_jprb, 1.04054e-24_jprb, 9.54054e-25_jprb, 9.31081e-25_jprb, 8.82432e-25_jprb, 8.27027e-25_jprb, &
      &  7.77027e-25_jprb, 6.68919e-25_jprb, 5.94595e-25_jprb, 5.66216e-25_jprb, 5.2973e-25_jprb, 4.91892e-25_jprb, &
      &  4.68919e-25_jprb, 4.2973e-25_jprb, 4.22973e-25_jprb, 4.36486e-25_jprb, 3.93243e-25_jprb, 3.48649e-25_jprb, &
      &  3.36486e-25_jprb, 3.25676e-25_jprb, 3.21622e-25_jprb, 3.13514e-25_jprb, 2.85135e-25_jprb, 2.66216e-25_jprb, &
      &  2.59459e-25_jprb, 2.48649e-25_jprb, 2.35135e-25_jprb, 2.24324e-25_jprb, 2.13514e-25_jprb, 2.01351e-25_jprb, &
      &  2e-25_jprb, 1.94595e-25_jprb, 1.85135e-25_jprb, 1.77027e-25_jprb, 1.62162e-25_jprb, 1.58108e-25_jprb, &
      &  1.60811e-25_jprb, 1.45946e-25_jprb, 1.35e-25_jprb, 1.44595e-25_jprb, 1.44595e-25_jprb, 1.36486e-25_jprb, &
      &  1.27973e-25_jprb, 1.22973e-25_jprb, 1.2e-25_jprb, 1.17027e-25_jprb, 1.22027e-25_jprb, 1.20946e-25_jprb, &
      &  1.17973e-25_jprb, 1.12027e-25_jprb, 1.07973e-25_jprb, 1.04054e-25_jprb, 9.77027e-26_jprb, 9.22973e-26_jprb, &
      &  9.77027e-26_jprb, 9.27027e-26_jprb, 7.94595e-26_jprb, 7.66216e-26_jprb, 8.36486e-26_jprb, 8.40541e-26_jprb, &
      &  6.97297e-26_jprb, 7.59459e-26_jprb, 8.48649e-26_jprb, 7.55405e-26_jprb, 6.55405e-26_jprb, 6.59459e-26_jprb, &
      &  6.67568e-26_jprb, 6.77027e-26_jprb, 6.83784e-26_jprb, 8.14865e-26_jprb, 8.98649e-26_jprb, 9.67568e-26_jprb, &
      &  8.86486e-26_jprb, 8.36486e-26_jprb, 9.2973e-26_jprb, 1.09054e-25_jprb, 1.15e-25_jprb, 1.48649e-25_jprb, &
      &  1.90541e-25_jprb, 2.05405e-25_jprb, 2.28378e-25_jprb, 2.59459e-25_jprb, 2.78378e-25_jprb, 3.48649e-25_jprb, &
      &  4.36486e-25_jprb, 5.41892e-25_jprb, 7.05405e-25_jprb, 8.74324e-25_jprb, 1.05946e-24_jprb, 1.32973e-24_jprb, &
      &  1.75676e-24_jprb, 2.21622e-24_jprb, 2.74324e-24_jprb, 3.51351e-24_jprb, 4.28378e-24_jprb, 5.51351e-24_jprb, &
      &  7.21622e-24_jprb, 8.90541e-24_jprb, 1.04054e-23_jprb, 1.34054e-23_jprb, 1.75676e-23_jprb, 2.16216e-23_jprb, &
      &  2.13793e-23_jprb, 2.1137e-23_jprb, 2.08947e-23_jprb, 2.06524e-23_jprb, 2.04101e-23_jprb, 2.01678e-23_jprb, &
      &  1.99254e-23_jprb, 1.96831e-23_jprb, 1.94408e-23_jprb, 1.91985e-23_jprb, 1.89562e-23_jprb, 1.87139e-23_jprb, &
      &  1.84716e-23_jprb, 1.82293e-23_jprb, 1.7987e-23_jprb, 1.77446e-23_jprb, 1.75023e-23_jprb, 1.726e-23_jprb, &
      &  1.70177e-23_jprb, 1.67754e-23_jprb, 1.65331e-23_jprb, 1.62908e-23_jprb, 1.60485e-23_jprb, 1.58062e-23_jprb, &
      &  1.55638e-23_jprb, 1.53215e-23_jprb, 1.50792e-23_jprb, 1.48369e-23_jprb, 1.45946e-23_jprb, 1.32973e-23_jprb, &
      &  9.67568e-24_jprb, 7.71622e-24_jprb, 5.74324e-24_jprb, 4.27027e-24_jprb, 3.17568e-24_jprb, 2.93243e-24_jprb, &
      &  2.32432e-24_jprb, 1.78378e-24_jprb, 1.37838e-24_jprb, 1.05946e-24_jprb, 8.44595e-25_jprb, 7.31081e-25_jprb, &
      &  5.89189e-25_jprb, 4.74324e-25_jprb, 3.83784e-25_jprb, 3.77027e-25_jprb, 3.17568e-25_jprb, 2.7027e-25_jprb, &
      &  2.2973e-25_jprb, 2e-25_jprb, 1.90541e-25_jprb, 1.67568e-25_jprb, 1.5e-25_jprb, 1.35e-25_jprb, &
      &  1.17027e-25_jprb, 1.10946e-25_jprb, 1.04054e-25_jprb, 9.43243e-26_jprb, 8.90541e-26_jprb, 8.41892e-26_jprb, &
      &  7.95946e-26_jprb, 7.54054e-26_jprb, 7.13514e-26_jprb, 6.74324e-26_jprb, 6.40541e-26_jprb, 6.09459e-26_jprb, &
      &  5.85135e-26_jprb, 5.62162e-26_jprb, 5.40541e-26_jprb, 5.21622e-26_jprb, 5.01351e-26_jprb, 4.40541e-26_jprb, &
      &  4.31081e-26_jprb, 4.21622e-26_jprb, 4.12162e-26_jprb, 4.02703e-26_jprb, 3.94595e-26_jprb, 3.86486e-26_jprb, &
      &  3.78378e-26_jprb, 3.7027e-26_jprb, 3.62162e-26_jprb, 3.55405e-26_jprb, 3.48649e-26_jprb, 3.40541e-26_jprb, &
      &  3.33784e-26_jprb, 3.28378e-26_jprb, 3.21622e-26_jprb, 3.17568e-26_jprb, 3.12162e-26_jprb, 3.06757e-26_jprb, &
      &  3.07389e-26_jprb, 3.08022e-26_jprb, 3.08654e-26_jprb, 3.09287e-26_jprb, 3.09919e-26_jprb, 3.10552e-26_jprb, &
      &  3.11185e-26_jprb, 3.11817e-26_jprb, 3.1245e-26_jprb, 3.13082e-26_jprb, 3.13715e-26_jprb, 3.14347e-26_jprb, &
      &  3.1498e-26_jprb, 3.15612e-26_jprb, 3.16245e-26_jprb, 3.16878e-26_jprb, 3.1751e-26_jprb, 3.18143e-26_jprb, &
      &  3.18775e-26_jprb, 3.19408e-26_jprb, 3.2004e-26_jprb, 3.20673e-26_jprb, 3.21305e-26_jprb, 3.21938e-26_jprb, &
      &  3.2257e-26_jprb, 3.23203e-26_jprb, 3.23836e-26_jprb, 3.24468e-26_jprb, 3.25101e-26_jprb, 3.25733e-26_jprb, &
      &  3.26366e-26_jprb, 3.26998e-26_jprb, 3.27631e-26_jprb, 3.28263e-26_jprb, 3.28896e-26_jprb, 3.29528e-26_jprb, &
      &  3.30161e-26_jprb, 3.30794e-26_jprb, 3.31426e-26_jprb, 3.32059e-26_jprb, 3.32691e-26_jprb, 3.33324e-26_jprb, &
      &  3.33956e-26_jprb, 3.34589e-26_jprb, 3.35221e-26_jprb, 3.35854e-26_jprb, 3.36486e-26_jprb, 3.78378e-26_jprb, &
      &  4.28378e-26_jprb, 5.18919e-26_jprb, 6.39189e-26_jprb, 8e-26_jprb, 1.00946e-25_jprb, 1.3e-25_jprb, &
      &  1.55405e-25_jprb, 1.86486e-25_jprb, 2.24324e-25_jprb, 2.87838e-25_jprb, 3.71622e-25_jprb, 4.82432e-25_jprb, &
      &  6.25676e-25_jprb, 6.95946e-25_jprb, 8.77027e-25_jprb, 1.24054e-24_jprb, 1.67568e-24_jprb, 1.48649e-24_jprb, &
      &  1.48649e-24_jprb, 1.94595e-24_jprb, 1.78378e-24_jprb, 1.66216e-24_jprb, 1.40541e-24_jprb, 1.10946e-24_jprb, &
      &  1.04054e-24_jprb, 1.14054e-24_jprb, 9.12162e-25_jprb, 1.1e-24_jprb, 1.67568e-24_jprb, 1.89189e-24_jprb, &
      &  2.10811e-24_jprb, 2.45946e-24_jprb, 2.71622e-24_jprb, 2.51351e-24_jprb, 2.08108e-24_jprb, 2.32432e-24_jprb, &
      &  2.68919e-24_jprb, 3.7973e-24_jprb, 3.95946e-24_jprb, 4.13514e-24_jprb, 4.21622e-24_jprb, 4.78378e-24_jprb, &
      &  5.33784e-24_jprb, 6e-24_jprb, 6.63514e-24_jprb, 7.45946e-24_jprb, 8.28378e-24_jprb, 9.24324e-24_jprb, &
      &  1.02973e-23_jprb, 1.15946e-23_jprb, 1.29054e-23_jprb, 1.41892e-23_jprb, 1.38345e-23_jprb, 1.34797e-23_jprb, &
      &  1.3125e-23_jprb, 1.27703e-23_jprb, 1.24155e-23_jprb, 1.20608e-23_jprb, 1.17061e-23_jprb, 1.13514e-23_jprb, &
      &  1.09966e-23_jprb, 1.06419e-23_jprb, 1.02872e-23_jprb, 9.93243e-24_jprb, 9.5777e-24_jprb, 9.22297e-24_jprb, &
      &  8.86824e-24_jprb, 8.51351e-24_jprb, 8.15878e-24_jprb, 7.80405e-24_jprb, 7.44932e-24_jprb, 7.09459e-24_jprb, &
      &  6.73986e-24_jprb, 6.38514e-24_jprb, 6.03041e-24_jprb, 5.67568e-24_jprb, 5.32095e-24_jprb, 4.96622e-24_jprb, &
      &  4.61149e-24_jprb, 4.25676e-24_jprb, 3.90203e-24_jprb, 3.5473e-24_jprb, 3.19257e-24_jprb, 2.83784e-24_jprb, &
      &  1.98649e-24_jprb, 1.58108e-24_jprb, 1.25946e-24_jprb, 1.00946e-24_jprb, 8.39189e-25_jprb, 7.5e-25_jprb, &
      &  5.51351e-25_jprb, 4.89189e-25_jprb, 4.33784e-25_jprb, 3.85135e-25_jprb, 3.40541e-25_jprb, 3.02703e-25_jprb, &
      &  2.67568e-25_jprb, 2.37838e-25_jprb, 2.05405e-25_jprb, 1.72973e-25_jprb, 1.44595e-25_jprb, 1.27027e-25_jprb, &
      &  1.12027e-25_jprb, 9.78378e-26_jprb, 8.59459e-26_jprb, 7.55405e-26_jprb, 6.63514e-26_jprb, 5.86486e-26_jprb, &
      &  5.17568e-26_jprb, 4.58108e-26_jprb, 4.12162e-26_jprb, 3.81081e-26_jprb, 3.52703e-26_jprb, 3.27027e-26_jprb, &
      &  3.04054e-26_jprb, 2.82432e-26_jprb, 2.63514e-26_jprb, 2.62462e-26_jprb, 2.61411e-26_jprb, 2.6036e-26_jprb, &
      &  2.59309e-26_jprb, 2.58258e-26_jprb, 2.57207e-26_jprb, 2.56156e-26_jprb, 2.55105e-26_jprb, 2.54054e-26_jprb, &
      &  2.53003e-26_jprb, 2.51952e-26_jprb, 2.50901e-26_jprb, 2.4985e-26_jprb, 2.48799e-26_jprb, 2.47748e-26_jprb, &
      &  2.46697e-26_jprb, 2.45646e-26_jprb, 2.44595e-26_jprb, 2.43544e-26_jprb, 2.42492e-26_jprb, 2.41441e-26_jprb, &
      &  2.4039e-26_jprb, 2.39339e-26_jprb, 2.38288e-26_jprb, 2.37237e-26_jprb, 2.36186e-26_jprb, 2.35135e-26_jprb, &
      &  2.34084e-26_jprb, 2.33033e-26_jprb, 2.31982e-26_jprb, 2.30931e-26_jprb, 2.2988e-26_jprb, 2.28829e-26_jprb, &
      &  2.27778e-26_jprb, 2.26727e-26_jprb, 2.25676e-26_jprb, 2.24625e-26_jprb, 2.23574e-26_jprb, 2.22523e-26_jprb, &
      &  2.21471e-26_jprb, 2.2042e-26_jprb, 2.19369e-26_jprb, 2.18318e-26_jprb, 2.17267e-26_jprb, 2.16216e-26_jprb, &
      &  2.35135e-26_jprb, 2.56757e-26_jprb, 2.82432e-26_jprb, 3.09459e-26_jprb, 3.40541e-26_jprb, 3.74324e-26_jprb, &
      &  4.14865e-26_jprb, 4.59459e-26_jprb, 5.08108e-26_jprb, 5.67568e-26_jprb, 5.90333e-26_jprb, 6.13098e-26_jprb, &
      &  6.35863e-26_jprb, 6.58628e-26_jprb, 6.81393e-26_jprb, 7.04158e-26_jprb, 7.26923e-26_jprb, 7.49688e-26_jprb, &
      &  7.72453e-26_jprb, 7.95218e-26_jprb, 8.17983e-26_jprb, 8.40748e-26_jprb, 8.63514e-26_jprb, 8.56757e-26_jprb, &
      &  8.40541e-26_jprb, 8.04054e-26_jprb, 7.7973e-26_jprb, 8.13514e-26_jprb, 8.47297e-26_jprb, 9.40541e-26_jprb, &
      &  1.04054e-25_jprb, 1.15946e-25_jprb, 1.29054e-25_jprb, 1.44595e-25_jprb, 1.60811e-25_jprb, 1.98649e-25_jprb, &
      &  2.37838e-25_jprb, 3.21622e-25_jprb, 4.04054e-25_jprb, 4.39189e-25_jprb, 5.04054e-25_jprb, 6.27027e-25_jprb, &
      &  8.33784e-25_jprb, 1.02027e-24_jprb, 1.12973e-24_jprb, 1.2e-24_jprb, 1.25946e-24_jprb, 1.40541e-24_jprb, &
      &  1.59459e-24_jprb, 1.67568e-24_jprb, 1.63514e-24_jprb, 1.59459e-24_jprb, 1.59459e-24_jprb, 1.5e-24_jprb, &
      &  1.27973e-24_jprb, 1.20946e-24_jprb, 1.27973e-24_jprb, 1.29054e-24_jprb, 1.39189e-24_jprb, 1.44595e-24_jprb, &
      &  1.44595e-24_jprb, 1.48649e-24_jprb, 1.59459e-24_jprb, 1.71622e-24_jprb, 1.87838e-24_jprb, 2.01351e-24_jprb, &
      &  2.08108e-24_jprb, 2.05405e-24_jprb, 1.94595e-24_jprb, 1.83784e-24_jprb, 1.72973e-24_jprb, 1.5e-24_jprb, &
      &  1.35e-24_jprb, 1.10946e-24_jprb, 1.00946e-24_jprb, 9e-25_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb, &
      &  0.0_jprb, 0.0_jprb, 0.0_jprb, 0.0_jprb ]


    ! Number of levels and wavenumbers
    integer, intent(in) :: nlev, nwav
    
    ! Pressure (Pa), temperature (K) and water vapour mole fraction
    ! (mol/mol)
    real(jprb), intent(in) :: pressure(nlev), temperature(nlev), mole_fraction(nlev)
    
    ! Wavenumber (cm-1)
    real(jprb), intent(in) :: wavenumber_cm1(nwav)

    ! Continuum molar absorption (m2 mol-1)
    real(jprb), intent(out) :: continuum(nwav,nlev)

    ! Continuum molar absorption (m2 mol-1) at spectral resolution of look-up table
    real(jprb) :: continuum_coarse(nwav_lut)
    
    ! Loop index for level
    integer :: jlev

    do jlev = 1,nlev
      ! MODIFY THIS using temperature(jlev), pressure(jlev), mole_fraction(jlev)
      continuum_coarse = self_lut + foreign_lut

      ! Convert from cm2 per molecule to m2 per mole
      continuum_coarse = continuum_coarse * (0.0001_jprb * NAVOGADRO)

      ! Interpolate to the input wavenumber
      call interpolate(wavenumber_lut_cm1, continuum_coarse, wavenumber_cm1, continuum(:,jlev))
    end do
    
  end subroutine calc_caviar_continuum
  
end module caviar_continuum
