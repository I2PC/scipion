#ifndef _UNITCELL_HH
#define _UNITCELL_HH
/*
#include "matrix1d.h"
#include "matrix2d.h"
#include "xmipp_funcs.h"
#include "args.h"
#include "grids.h"
*/
#include "symmetries.h"
#include <vector>
#include <iterator>
#include "xmipp_image_generic.h"
#include "data/xmipp_image.h"

#define VALIDNUMBERS 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 30, 32, 33, 34, 35, 36, 38, 39, 40, 42, 44, 45, 48, 49, 50, 51, 52, 54, 55, 56, 57, 60, 63, 64, 65, 66, 68, 70, 72, 75, 76, 77, 78, 80, 81, 84, 85, 88, 90, 91, 95, 96, 98, 99, 100, 102, 104, 105, 108, 110, 112, 114, 117, 119, 120, 121, 125, 126, 128, 130, 132, 133, 135, 136, 140, 143, 144, 147, 150, 152, 153, 154, 156, 160, 162, 165, 168, 169, 170, 171, 175, 176, 180, 182, 187, 189, 190, 192, 195, 196, 198, 200, 204, 208, 209, 210, 216, 220, 221, 224, 225, 228, 231, 234, 238, 240, 242, 243, 245, 247, 250, 252, 255, 256, 260, 264, 266, 270, 272, 273, 275, 280, 285, 286, 288, 289, 294, 297, 300, 304, 306, 308, 312, 315, 320, 323, 324, 325, 330, 336, 338, 340, 342, 343, 350, 351, 352, 357, 360, 361, 363, 364, 374, 375, 378, 380, 384, 385, 390, 392, 396, 399, 400, 405, 408, 416, 418, 420, 425, 429, 432, 440, 441, 442, 448, 450, 455, 456, 459, 462, 468, 475, 476, 480, 484, 486, 490, 494, 495, 500, 504, 507, 510, 512, 513, 520, 525, 528, 532, 539, 540, 544, 546, 550, 560, 561, 567, 570, 572, 576, 578, 585, 588, 594, 595, 600, 605, 608, 612, 616, 624, 625, 627, 630, 637, 640, 646, 648, 650, 660, 663, 665, 672, 675, 676, 680, 684, 686, 693, 700, 702, 704, 714, 715, 720, 722, 726, 728, 729, 735, 741, 748, 750, 756, 760, 765, 768, 770, 780, 784, 792, 798, 800, 810, 816, 819, 825, 832, 833, 836, 840, 845, 847, 850, 855, 858, 864, 867, 875, 880, 882, 884, 891, 896, 900, 910, 912, 918, 924, 931, 935, 936, 945, 950, 952, 960, 968, 969, 972, 975, 980, 988, 990, 1000, 1001, 1008, 1014, 1020, 1024, 1026, 1029, 1040, 1045, 1050, 1053, 1056, 1064, 1071, 1078, 1080, 1083, 1088, 1089, 1092, 1100, 1105, 1120, 1122, 1125, 1134, 1140, 1144, 1152, 1155, 1156, 1170, 1176, 1183, 1188, 1190, 1197, 1200, 1210, 1215, 1216, 1224, 1225, 1232, 1235, 1248, 1250, 1254, 1260, 1274, 1275, 1280, 1287, 1292, 1296, 1300, 1309, 1320, 1323, 1326, 1330, 1331, 1344, 1350, 1352, 1360, 1365, 1368, 1372, 1375, 1377, 1386, 1400, 1404, 1408, 1425, 1428, 1430, 1440, 1444, 1445, 1452, 1456, 1458, 1463, 1470, 1482, 1485, 1496, 1500, 1512, 1520, 1521, 1530, 1536, 1539, 1540, 1547, 1560, 1568, 1573, 1575, 1584, 1596, 1600, 1615, 1617, 1620, 1625, 1632, 1638, 1650, 1664, 1666, 1672, 1680, 1683, 1690, 1694, 1700, 1701, 1710, 1715, 1716, 1728, 1729, 1734, 1750, 1755, 1760, 1764, 1768, 1782, 1785, 1792, 1800, 1805, 1815, 1820, 1824, 1836, 1848, 1859, 1862, 1870, 1872, 1875, 1881, 1890, 1900, 1904, 1911, 1920, 1925, 1936, 1938, 1944, 1950, 1960, 1976, 1980, 1989, 1995, 2000, 2002, 2016, 2023, 2025, 2028, 2040, 2048, 2052, 2057, 2058, 2079, 2080, 2090, 2100, 2106, 2112, 2125, 2128, 2142, 2145, 2156, 2160, 2166, 2176, 2178, 2184, 2187, 2197, 2200, 2205, 2210, 2223, 2240, 2244, 2250, 2261, 2268, 2275, 2280, 2288, 2295, 2299, 2304, 2310, 2312, 2340, 2352, 2366, 2375, 2376, 2380, 2394, 2400, 2401, 2420, 2430, 2431, 2432, 2448, 2450, 2457, 2464, 2470, 2475, 2496, 2499, 2500, 2508, 2520, 2527, 2535, 2541, 2548, 2550, 2560, 2565, 2574, 2584, 2592, 2600, 2601, 2618, 2625, 2640, 2646, 2652, 2660, 2662, 2673, 2688, 2695, 2700, 2704, 2717, 2720, 2730, 2736, 2744, 2750, 2754, 2772, 2793, 2800, 2805, 2808, 2816, 2835, 2850, 2856, 2860, 2873, 2880, 2888, 2890, 2904, 2907, 2912, 2916, 2925, 2926, 2940, 2964, 2970, 2975, 2992, 3000, 3003, 3024, 3025, 3040, 3042, 3060, 3072, 3078, 3080, 3087, 3094, 3120, 3125, 3135, 3136, 3146, 3150, 3159, 3168, 3179, 3185, 3192, 3200, 3211, 3213, 3230, 3234, 3240, 3249, 3250, 3264, 3267, 3276, 3300, 3315, 3325, 3328, 3332, 3344, 3360, 3366, 3375, 3380, 3388, 3400, 3402, 3420, 3430, 3432, 3456, 3458, 3465, 3468, 3500, 3510, 3520, 3528, 3536, 3549, 3553, 3564, 3570, 3575, 3584, 3591, 3600, 3610, 3630, 3640, 3645, 3648, 3672, 3675, 3696, 3705, 3718, 3724, 3740, 3744, 3750, 3757, 3762, 3773, 3780, 3800, 3808, 3822, 3825, 3840, 3850, 3861, 3872, 3876, 3888, 3900, 3920, 3927, 3952, 3960, 3969, 3971, 3978, 3990, 3993, 4000, 4004, 4032, 4046, 4050, 4056, 4080, 4095, 4096, 4104, 4114, 4116, 4125, 4131, 4158, 4160, 4165, 4180, 4199, 4200, 4212, 4224, 4225, 4235, 4250, 4256, 4275, 4284, 4290, 4312, 4320, 4332, 4335, 4352, 4356, 4368, 4374, 4375, 4389, 4394, 4400, 4410, 4420, 4446, 4455, 4459, 4480, 4488, 4500, 4522, 4536, 4550, 4560, 4563, 4576, 4590, 4598, 4608, 4617, 4620, 4624, 4641, 4655, 4675, 4680, 4693, 4704, 4719, 4725, 4732, 4750, 4752, 4760, 4788, 4800, 4802, 4840, 4845, 4851, 4860, 4862, 4864, 4875, 4896, 4900, 4913, 4914, 4928, 4940, 4950, 4992, 4998, 5000, 5005, 5016, 5040, 5049, 5054, 5070, 5082, 5096, 5100, 5103, 5120, 5130, 5145, 5148, 5168, 5184, 5187, 5200, 5202, 5225, 5236, 5250, 5265, 5280, 5292, 5304, 5320, 5324, 5346, 5355, 5376, 5390, 5400, 5408, 5415, 5434, 5440, 5445, 5460, 5472, 5488, 5491, 5500, 5508, 5525, 5544, 5577, 5586, 5600, 5610, 5616, 5625, 5632, 5643, 5670, 5700, 5712, 5720, 5733, 5746, 5760, 5775, 5776, 5780, 5808, 5814, 5824, 5831, 5832, 5850, 5852, 5880, 5915, 5928, 5929, 5940, 5950, 5967, 5984, 5985, 6000, 6006, 6048, 6050, 6069, 6075, 6080, 6084, 6120, 6125, 6137, 6144, 6156, 6160, 6171, 6174, 6175, 6188, 6237, 6240, 6250, 6270, 6272, 6292, 6300, 6318, 6336, 6358, 6370, 6375, 6384, 6400, 6422, 6426, 6435, 6460, 6468, 6480, 6498, 6500, 6517, 6528, 6534, 6545, 6552, 6561, 6591, 6600, 6615, 6630, 6650, 6655, 6656, 6664, 6669, 6688, 6720, 6732, 6750, 6760, 6776, 6783, 6800, 6804, 6825, 6840, 6859, 6860, 6864, 6875, 6885, 6897, 6912, 6916, 6930, 6936, 7000, 7007, 7020, 7040, 7056, 7072, 7098, 7106, 7125, 7128, 7140, 7150, 7168, 7182, 7200, 7203, 7220, 7225, 7260, 7280, 7290, 7293, 7296, 7315, 7344, 7350, 7371, 7392, 7410, 7425, 7436, 7448, 7480, 7488, 7497, 7500, 7514, 7524, 7546, 7560, 7581, 7600, 7605, 7616, 7623, 7644, 7650, 7680, 7695, 7700, 7722, 7735, 7744, 7752, 7776, 7800, 7803, 7840, 7854, 7865, 7875, 7904, 7920, 7938, 7942, 7956, 7980, 7986, 8000, 8008, 8019, 8064, 8075, 8085, 8092, 8100, 8112, 8125, 8151, 8160, 8190, 8192, 8208, 8228, 8232, 8250, 8262, 8281, 8316, 8320, 8330, 8360, 8379, 8398, 8400, 8415, 8424, 8448, 8450, 8470, 8500, 8505, 8512, 8550, 8568, 8575, 8580, 8619, 8624, 8640, 8645, 8664, 8670, 8704, 8712, 8721, 8736, 8748, 8750, 8775, 8778, 8788, 8800, 8820, 8840, 8892, 8910, 8918, 8925, 8960, 8976, 9000, 9009, 9025, 9044, 9072, 9075, 9100, 9120, 9126, 9152, 9163, 9180, 9196, 9216, 9234, 9240, 9248, 9261, 9282, 9295, 9310, 9317, 9350, 9360, 9375, 9386, 9405, 9408, 9438, 9450, 9464, 9477, 9500, 9504, 9520, 9537, 9555, 9576, 9600, 9604, 9625, 9633, 9639, 9680, 9690, 9702, 9720, 9724, 9728, 9747, 9750, 9792, 9800, 9801, 9826, 9828, 9856, 9880, 9900, 9945, 9975, 9984, 9996

class UnitCell{

/** Typename to contain a list of plane_vector */

private:
	//point group symmetry
	int symmetry;
	//pointgroup symmetry order
	int sym_order;
	//symmetry object
	SymList SL;
	//min max radii for masking crown
	double rmin, rmax;
	//offset for unitcell
	double offset;
	//sampling rate
	double sampling;
	//vectors defining expansion direction
	double expanded;
	//origin x coordinate introduced by the user with the input volume
	double x_origin;
	//origin y coordinate introduced by the user with the input volume
	double y_origin;
	//origin z coordinate introduced by the user with the input volume
	double z_origin;
	//new origin after expansion
	Matrix1D<double>  newOriginAfterExpansion;
	//
	std::vector<Matrix1D<double> > vectExpansion;
	// normal vector to the polyhedron that define a unit cell.
	std::vector<Matrix1D<double> > planeVectors;
	//vectors with corner of expanded unit cell
	std::vector<Matrix1D<double> > expandedUnitCell;
	//above vector multiplied by rmin
	std::vector<Matrix1D<double> > expandedUnitCellMin;
	//above vector multiplied by rmax
	std::vector<Matrix1D<double> > expandedUnitCellMax;
	//reference point always inside the unitcell
	Matrix1D<double> unitCellPoint;

	//handle cyclic symmetry (Cn)
	//input 3 vertices that form a triangle
	//the function fills out the variable planeVectors (normal vectors to the polyhedron that define a unit cell)
	void cyclicSymmetry(const Matrix1D<double> & _centre,
			            const Matrix1D<double> & _2f,
			            const Matrix1D<double> & _2fp,
			            const int order,
			            const double expanded,
			            const double offset);
	//handle dihedral symmetry (Dn)
	//input 3 vertices that form the triangle
	//the function fills out the variable planeVectors (normal vectors to the polyhedron that define a unit cell)
	void dihedralSymmetry(const Matrix1D<double> & _centre,
			              const Matrix1D<double> & _2f,
			              const Matrix1D<double> & _2fp,
			              const int order,
			              const double expanded,
			              const double offset);
	//handle tetrahedral symmetry (T)
	//input centroid and 3 vertices that form the unit cell
	//the function fills out the variable planeVectors (normal vectors to the polyhedron that define a unit cell)
	void tetrahedralSymmetry(const Matrix1D<double> & _centroid,
				             const Matrix1D<double> & _3f,
				             const Matrix1D<double> & _3fp,
				             const Matrix1D<double> & _3fpp,
				             const double expanded);
	//handle tetrahedral symmetry (T)
	//input centre and 3 vertices that form the unit cell
	//the function fills out the variable planeVectors (normal vectors to the polyhedron that define a unit cell)
	void octahedralSymmetry(const Matrix1D<double> & _centre,
					        const Matrix1D<double> & _4f,
					        const Matrix1D<double> & _4fp,
					        const Matrix1D<double> & _4fpp,
					        const double expanded);
	//handle icosahedral case.
	//input centre plus 3 five fold vertices that form a triangle
	//the function fills out the variable planeVectors (normal vectors to the polyhedron that define a unit cell)
	void icoSymmetry(const Matrix1D<double> & _centre,
					 const Matrix1D<double> & _5f,
		             const Matrix1D<double> & _5fp,
		             const Matrix1D<double> & _5fpp,
		             const double expanded);
private:
        //maxmimum/minimm zvalue
        double _minZ, _maxZ;
public:
	/** set to zero everything that is not the unit cell */
	void maskUnitCell(ImageGeneric & in3DDmap,
			          ImageGeneric & out3DDmap);
	/** close debug file */
	//~UnitCell();
	/** unit auxiliary vectors */
	UnitCell(String sym, double rmin, double rmax, double expanded,
			double offset, double sampling, double x_origin, double y_origin, double z_origin);
	/** process CN symmetry */
	void doCyclic(int order=1);
	/** process DN symmetry */
	void doDihedral(int order=1);
	/** process T symmetry */
	void doTetrahedral(int symmetry=pg_T);
	/** process O symmetry */
	void doOctahedral(int symmetry=pg_O);
	/** process IN symmetry */
	void doIcosahedral(int symmetry=pg_I1);

};//end unitcell class

#endif
