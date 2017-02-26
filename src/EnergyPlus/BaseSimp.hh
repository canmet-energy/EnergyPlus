#ifndef BaseSimp_hh_INCLUDED
#define BaseSimp_hh_INCLUDED

// C++ Headers
#include <memory>

// ObjexxFCL Headers
#include <ObjexxFCL/Array1D.hh>

// EnergyPlus Headers
#include <EnergyPlus.hh>
#include <DataGlobals.hh>

namespace EnergyPlus {

namespace BaseSimp {

struct BSFoundationSpecs
{
	// Public Members
	std::string Name; // user identifier
	std::string BSFoundationName;			// Name of BASESIMP object
	std::string BSZoneName;					// Name of corresponding foundation zone in Energylus model
	int BSZoneNumber;						// Number of of corresponding foundation zone in Energylus model
	int BSFndType;							// Foundation Type (1=Basement or 2=Slab-in-grade)
	int BSModelType;						// Model Type (1=Foundation description or 2=BASECALC Coefficients)
	int BSFndMat;							// Foundation Material (1=Concrete, 2=Wood, 3=Wood/Concrete)
	int BSFndConfig;						// Foundation Configuration
	Real64 BSExposedFraction;				// Fraction of the foundation exposed to ambient/soil
	Real64 BSHeight;						// Foundation height (m)
	Real64 BSDepth;							// Foundation depth (m)
	Real64 BSLength;						// Foundation length (m)
	Real64 BSWidth;							// Foundation width (m)
	Real64 BSOverlap;						// Insulation overlap (m)
	Real64 BSRSI;							// Insulation resistance in RSI (m2K/W)
	Real64 BSSag;							// BASECALC above grade heat losses coefficient (W/K)
	Real64 BSSbgavg;						// BASECALC below grade average heat losses coefficient (W/K)
	Real64 BSSbgvar;						// BASECALC below grade variable heat losses coefficient (W/K)
	Real64 BSPhase;							// BASECALC thermal response of the foundation/soil system (rad)
	Real64 BSSoilK;							// Soil thermal conductivity (W/mK)
	Real64 BSWTD;							// Water Table Depth (m)
	Real64 BSTGavg;							// Anually averaged soil temperature (deg-C)
	Real64 BSTGamp;							// Amplitude of ground temperature's annual sine wave (deg-C)
	Real64 BSTGps;							// Phase lag of ground temperature's annual sine wave (rad)
	Real64 BSFoundationLosses;				// Total heat losses for this foundation (W)
	Real64 BSIntGainsToZone;				// Internal gains to pass to the zone (W)
	Real64 BSPhaseAdj;						// Phase Lag adjusted by Pi/2
	Real64 BSZoneAvgTemp;					// Average zone temperature during the last 7 days
	Array1D< Real64 > BSZoneInstantTemps;	// Array containing the zone temperature for the last 7 days

	// Default Constructor
	BSFoundationSpecs() :
		BSExposedFraction(1.0),
		BSHeight(2.0),
		BSDepth(1.8),
		BSLength(10.0),
		BSWidth(5.0),
		BSOverlap(0.0),
		BSRSI(0.0),
		BSSoilK(2.0),
		BSWTD(10.0),
		BSSag(0.0),
		BSSbgavg(0.0),	
		BSSbgvar(0.0),	
		BSPhase(0.0),	
		BSTGavg(0.0),
		BSTGamp(0.0),
		BSTGps(0.0),
		BSFoundationLosses(0.0),
		BSIntGainsToZone(0.0),
		BSPhaseAdj(0.0),
		BSZoneAvgTemp(0.0)
		{}
	};


// Object Data
extern Array1D< BSFoundationSpecs > BaseSimp; // BASESIMP data 

// ********************************  Functions
void
clear_state();

// Controling routine: sends to initialisation on first call and sends to calculations for all iterations
void
CalcBasesimpGains();

// Routine calculating the gains (negative) to each zone
void
CalcBSGainsToZone();

// Routine calculating the heat losses from one specific foundation
static void
CalcBSFoundationHeatLosses(int const BSFoundationNum);

// Routine calculating the avarage zone temperature and updating the array of zone temperatures
static void
UpdateZonesTemps(int const BSFoundationNum);

// Controling routine for the initialisation 
static void
InitAllBSFoundations();

// Read and check all input parameters for a specific foundation
static void
GetBSFoundationInput(int const BSFoundationNum);

// Evaluates the heat loss factors for a foundation
// The main difference with the next routine (CalcBSFactors) is that this one does the exponential interpollation
// So this routine uses CalcBSFactors() to calculate coefficients and this one interpollates the final coefficients
static void
BSFactCtrl(int const BSFoundationNum);

// Evaluates the heat loss factors for a foundation with a specific configuration
static void
CalcBSFactors(int const BSFoundationNum, int const BSFoundationCfg);

// Evaluates the coefficients used to calculate the heat loss factors
// The value of these coefficient depends on the configuration
static void
GetBSCoeff(int const BSFoundationNum, int const BSFoundationCfg);

// Set the link between the heat losses calculated here and the internal gains in the zone energy balance
static void
SetBasesimpGainLink();

// Initialize the array containing the corner correction factors
static void
InitBSCornerCoeff();


// *********************   Static data 

// Intermediate coefficients
static Real64 BSa1, BSb1, BSc1, BSd1, BSe1, BSf1, BSg1, BSh1, BSi1, BSj1, BSq2,
BSr2, BSu2, BSv2, BSw2, BSx2, BSs2, BSt2, BSy2, BSa2, BSb2, BSc2, BSd2, BSe2, BSf2, BSg2, BSh2,
BSa3, BSb3, BSc3, BSe3, BSf3, BSg3, BSh3, BSi3, BSa4, BSb4, BSc4; 
static int ICol;

// Table containing corner correction factors 
// -> Initiated by one element more than required in each direction
// This is to start accessing from 1 instead of 0 to avoid potential errors in conversion from Fortran
static Real64	BSCornerCoeff[17][20];  // Corner coefficients

// Omega constant (period)
static Real64	BSOmega;

// Indicate if BASESIMP initialization is required
static bool InitBasesimpFlag(true);

// Indicates if BASESIMP environment has yet not been called
static bool BSEnvrnFlag(true);

// Number of time steps in a week
static int BSNumTimeStepsInWeek(0);		

// Zones initial temperature
static const Real64	BSTzoneInit(20.0);  // Initial zone temperature to determine heat exchange potential with soil & ambient

// Time (Hour of the year)
static Real64 HourOfYear(0.0);		// Hour of the year

// Number of BASESIMP objects in model
// Make public to access from gTester
extern int NumBSFoundation;

} // BaseSimp

} // EnergyPlus

#endif