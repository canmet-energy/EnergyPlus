// EnergyPlus, Copyright (c) 1996-2016, The Board of Trustees of the University of Illinois and
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights
// reserved.
//
// If you have questions about your rights to use or distribute this software, please contact
// Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
//
// NOTICE: This Software was developed under funding from the U.S. Department of Energy and the
// U.S. Government consequently retains certain rights. As such, the U.S. Government has been
// granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
// worldwide license in the Software to reproduce, distribute copies to the public, prepare
// derivative works, and perform publicly and display publicly, and to permit others to do so.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice, this list of
//     conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory,
//     the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be
//     used to endorse or promote products derived from this software without specific prior
//     written permission.
//
// (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form
//     without changes from the version obtained under this License, or (ii) Licensee makes a
//     reference solely to the software portion of its product, Licensee must refer to the
//     software as "EnergyPlus version X" software, where "X" is the version number Licensee
//     obtained under this License and may not use a different name for the software. Except as
//     specifically required in this Section (4), Licensee shall not use in a company name, a
//     product name, in advertising, publicity, or other promotional activities any name, trade
//     name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly
//     similar designation, without Lawrence Berkeley National Laboratory's prior written consent.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the
// features, functionality or performance of the source code ("Enhancements") to anyone; however,
// if you choose to make your Enhancements available either publicly, or directly to Lawrence
// Berkeley National Laboratory, without imposing a separate written license agreement for such
// Enhancements, then you hereby grant the following license: a non-exclusive, royalty-free
// perpetual license to install, use, modify, prepare derivative works, incorporate into other
// computer software, distribute, and sublicense such enhancements or derivative works thereof,
// in binary and source code form.

// C++ Headers
#include <cassert>
#include <cmath>

// ObjexxFCL Headers
#include <ObjexxFCL/Array.functions.hh>
#include <ObjexxFCL/Fmath.hh>

// EnergyPlus Headers
#include <BaseSimp.hh>
#include <InputProcessor.hh>
#include <OutputProcessor.hh>
#include <DataIPShortCuts.hh>
#include <General.hh>
#include <GlobalNames.hh>
#include <DataEnvironment.hh>
#include <DataGlobals.hh>
#include <DataHeatBalFanSys.hh>
#include <DataHeatBalance.hh>
#include <HeatBalanceInternalHeatGains.hh>


namespace EnergyPlus {

	namespace BaseSimp {

		// BaseSimp Foundation Model

		// MODULE INFORMATION:
		//       AUTHOR     Based on ESP-r model developped by Ian Beausoleil Morrison    
		//					Converted to CPP and integrated in EnergyPlus by Patrice Pinel
		//       DATE WRITTEN   2016
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS MODULE:
		// Model residential Foundations using the BASESIMP Canadian model

		// METHODOLOGY EMPLOYED:
		// The BASESIMP model evaluates regressions coefficients for Foundation heat losses as a function of indoor and outdoor conditions
		// It then evaluates Foundation heat losses at each time step by calculating the regression in question
		// The heat losses are passed to an EnergyPlus zone in the form of internal gains
		// For Basements:
		// The EnergyPlus basement zone does not need to incorporate walls and floors, only a ceiling connecting the basement to the zones above
		// The correct basement zone volume must be described in the "Zone" object representing a basement
		// For slab-on/in-grade Foundations:
		// The floor/slab should not be described in the zone located above the slab: if it is, it should be modelled as adiabatic

		// REFERENCES: Pinel P, Integration of the BASESIMP model in EnergyPlus, Report for PO # , Natural Resources Canada

		// OTHER NOTES: See report for references describing the BASESIMP method

		// USE STATEMENTS:
		// Use statements for data only modules
		// Using/Aliasing
		
		//USE FunctionFluidProperties
		// Use statements for access to subroutines in other modules

		// Data
		// MODULE PARAMETER DEFINITIONS

		// DERIVED TYPE DEFINITIONS

		// Notes:
		//	1- There are "AnyBasementInModel" and "AnySlabsInModel" flags in the class DataGlobals.
		//	   Their purpose seems to be to be able to identify foundations in the model and prevent the use of other Site:GroundDomain objects
		//	   At the moment, I don't see any possible negative consequence of runing BASESIMP and the other Site:GroundDomain objects
		//     So the flag is not initiated
		
		// MODULE VARIABLE DECLARATIONS:
		
		
		// SUBROUTINE SPECIFICATIONS FOR MODULE BASESIMP

		// Object Data
		Array1D< BSFoundationSpecs > BaseSimp; // BaseSimp data

		// MODULE SUBROUTINES:

		// Beginning of BaseSimp Module Driver Subroutines
		//*************************************************************************

		// Functions

		void
		clear_state()
		{
			BaseSimp.deallocate();
		}

		void
			CalcBasesimpGains()
		{
			// SUBROUTINE INFORMATION:
			//       AUTHOR         Patrice Pinel
			//       DATE WRITTEN   March 2016
			//       MODIFIED       na
			//       RE-ENGINEERED  na

			// PURPOSE OF THIS SUBROUTINE:
			// Calculates the zone internal gains injected to reproduce heat losses through a BASESIMP foundation
			
			// METHODOLOGY EMPLOYED:
			// This routine and its method are heavily influenced by the Water Heater model
			// See: CalcWaterThermalTankZoneGains() routine by Peter Graham Ellis in the WaterThermalTanks.cc file

			// This is the controling routing of the BASESIMP implementation
			// It is called from the InternalHeatGains.cc file whenever an energy balance is performed on the zone
			// If it is the first call to BASESIMP, this subroutine calls the initialisation routines
			// For all calls, this subroutine calls the calculation routines to assess the foundation heat losses 
			// and convert them into equivalent zone heat gains
			
			// Using/Aliasing
			using DataGlobals::BeginEnvrnFlag;
						
			// FIRST Call initialization
			if (InitBasesimpFlag) {

				// Initiate all BASESIMP foundations: Allocate memory, read parameters, and calculate coefficients
				InitAllBSFoundations();

				// Ensure the initiation is not repeated by setting the flag to false
				InitBasesimpFlag = false;
				
			}

			// Do calculation of equivalent zone internal heat gains to account for heat losses from foundations 
			// Only if there are foundations
			if (NumBSFoundation > 0) CalcBSGainsToZone();

		}  // End routine

		
		void
			CalcBSGainsToZone()
		{
			// SUBROUTINE INFORMATION:
			//       AUTHOR:          Patrice Pinel
			//       DATE WRITTEN:    March 2016
			//       MODIFIED:        
			//       RE-ENGINEERED:   

			// PURPOSE OF THIS SUBROUTINE:
			// Evaluate internal gains to send to zones to account for BASESIMP foundations 

			// METHODOLOGY EMPLOYED:
			// 1- Evaluate the time (hour of the year)
			// 2- Update the average temperature for all zones including a BASESIMP foundation
			// 3- Assess the heat losses from all foundations
			// 4- Convert these heat losses into equivalent zone internal gains to pass to the zone energy balance routines

			// REFERENCES: ESP-r source code

			//Using/Aliasing
			using DataEnvironment::DayOfYear;
			using DataGlobals::HourOfDay;
			
			// Local variables
			int BSFoundationNum;				// Foundation counter for loops

			// Evaluate the time
			HourOfYear = (DayOfYear - 1.0)*24. + HourOfDay;

			// Calculate basement heat loss for this call
			for (BSFoundationNum = 1; BSFoundationNum <= NumBSFoundation; ++BSFoundationNum) {

				// Evaluate zones temperature (average of the last two weeks)
				UpdateZonesTemps(BSFoundationNum);

				// Assess the heat losses from the foundation
				CalcBSFoundationHeatLosses(BSFoundationNum);

				// Internal gains to zone are negative for heat losses
				BaseSimp(BSFoundationNum).BSIntGainsToZone = -BaseSimp(BSFoundationNum).BSFoundationLosses;

			}

		} //  End routine


		void
			CalcBSFoundationHeatLosses(int const BSFoundationNum)
		{
			// SUBROUTINE INFORMATION:
			//       AUTHOR:          Patrice Pinel
			//       DATE WRITTEN:    March 2016
			//       MODIFIED:        
			//       RE-ENGINEERED:   

			// PURPOSE OF THIS SUBROUTINE:
			// Calculate the heat losses (negative for gain) from a foundation modelled with BASESIMP 

			// METHODOLOGY EMPLOYED:
			// Translates the BSHEAT routine from file "basefimp.F" of the ESP-r source code

			// REFERENCES: ESP-r source code

			//Using/Aliasing
			using DataEnvironment::OutDryBulbTemp;

			// BASESIMP individual heatlosses
			Real64 BSQag;			// Above grade
			Real64 BSQbgAvg;		// Below grade average
			Real64 BSQbgVar;		// Below grade variable

			int ZoneNum;			// Zone Number
			
			// Above - grade heat loss(W).The outdoor air temperature is set in Class DataEnvironment.
			BSQag = BaseSimp(BSFoundationNum).BSSag*(BaseSimp(BSFoundationNum).BSZoneAvgTemp - OutDryBulbTemp);

			// Steady component of below - grade heat loss(W).
			// Will now vary since the zone temp is allowed to varry
			BSQbgAvg = BaseSimp(BSFoundationNum).BSSbgavg*(BaseSimp(BSFoundationNum).BSZoneAvgTemp - BaseSimp(BSFoundationNum).BSTGavg);

			// Varying component of below - grade heat loss(W). `HourOfYear' measures the time in hours from the beginning of the calendar year.
			// BSOmega was calculated once at the begining of the simulation
			// BSPhaseAdj, the phase lag adjusted by Pi/2, was calculated for each foundatiuon at the begining of the simulation
			BSQbgVar = BaseSimp(BSFoundationNum).BSSbgvar * BaseSimp(BSFoundationNum).BSTGamp * sin(BSOmega*HourOfYear - BaseSimp(BSFoundationNum).BSPhaseAdj);

			// Total heat losses
			BaseSimp(BSFoundationNum).BSFoundationLosses = (BSQag + BSQbgAvg + BSQbgVar)*BaseSimp(BSFoundationNum).BSExposedFraction;
		
		} //  End routine



		void
			UpdateZonesTemps(int const BSFoundationNum)
		{
			// SUBROUTINE INFORMATION:
			//       AUTHOR:          Patrice Pinel
			//       DATE WRITTEN:    March 2016
			//       MODIFIED:        
			//       RE-ENGINEERED:   

			// PURPOSE OF THIS SUBROUTINE:
			// Calculate the average temperature of the zones where the foundations are located
			// Maintain an array of temperatures, containing a week of data, for each of these zones

			// METHODOLOGY EMPLOYED:
			// Update array when the first iteration of time step
			// Update the average at every iteration
			// Average updated by adding new terms and removing exiting terms to it: not by adding all terms and dividing by number of terms
			// This requires much less calculations

			// REFERENCES: ESP-r source code

			//Using/Aliasing
			using DataHeatBalFanSys::MAT;
			using DataGlobals::BeginTimeStepFlag;
			
			int BSTimeStepNum;		// Time step counter for loops
			
			if (BeginTimeStepFlag) {  // new time step -> update array and average

				// Update average temperature
				// Add term for current temperature and remove term for last temperature in array
				BaseSimp(BSFoundationNum).BSZoneAvgTemp = BaseSimp(BSFoundationNum).BSZoneAvgTemp + (MAT(BaseSimp(BSFoundationNum).BSZoneNumber) - BaseSimp(BSFoundationNum).BSZoneInstantTemps(BSNumTimeStepsInWeek)) / static_cast<float>(BSNumTimeStepsInWeek);

				// Update array
				for (BSTimeStepNum = 1; BSTimeStepNum <= BSNumTimeStepsInWeek - 1; ++BSTimeStepNum) {

					// Move all items in the array to the right by one position 
					BaseSimp(BSFoundationNum).BSZoneInstantTemps(BSTimeStepNum + 1) = BaseSimp(BSFoundationNum).BSZoneInstantTemps(BSTimeStepNum);
				}

				// First item is the the present temp
				BaseSimp(BSFoundationNum).BSZoneInstantTemps(1) = MAT(BaseSimp(BSFoundationNum).BSZoneNumber);
				
			}
			else {  // Only update average and first item of array
				
				// Average requires removing first term of the array and replacing it with the current temperature
				BaseSimp(BSFoundationNum).BSZoneAvgTemp = BaseSimp(BSFoundationNum).BSZoneAvgTemp + (MAT(BaseSimp(BSFoundationNum).BSZoneNumber) - BaseSimp(BSFoundationNum).BSZoneInstantTemps(1)) / static_cast<float>(BSNumTimeStepsInWeek);

				// Change first item of array for actual temp
				BaseSimp(BSFoundationNum).BSZoneInstantTemps(1) = MAT(BaseSimp(BSFoundationNum).BSZoneNumber);

			} // End else loop for first iter of the time step

		}	// end routine



		void
			InitAllBSFoundations()
		{
			// SUBROUTINE INFORMATION:
			//       AUTHOR:          Patrice Pinel
			//       DATE WRITTEN:    March 2016
			//       MODIFIED:        
			//       RE-ENGINEERED:   

			// PURPOSE OF THIS SUBROUTINE:
			// Handles all initiations for BASESIMP foundations at the begining of the simulation including:
			// - Evaluation of parameters required for memory allocation
			// - Reading of input parameter (from idf file) for each foundation
			// - Evaluation of BASESIMP heat loss factors for each foundation

			// METHODOLOGY EMPLOYED:

			// REFERENCES: na

			// Note: This routine should be called only once at the begining of the simulation

			//Using/Aliasing
			using InputProcessor::GetNumObjectsFound;
			using namespace DataIPShortCuts;				// Data for field names, blank numerics
			using DataGlobals::NumOfTimeStepInHour;
			using DataGlobals::Pi;
			using DataHeatBalance::IntGainTypeOf_Basesimp;	// Connection to zone internal gains
								
			//LOCAL VARIABLES
			int BSFoundationNum;				// Foundation counter for loops
			int BSTimeStepNum;					// Time step counter for loops
			Array1D_bool CheckEquipName;

			//Identify the EnergyPlus object being read
			cCurrentModuleObject = "Site:GroundDomain:BASESIMP";
		
			//GET NUMBER OF BASESIMP Models
			NumBSFoundation = GetNumObjectsFound(cCurrentModuleObject);
			CheckEquipName.dimension(NumBSFoundation, true);

			// If there are no foundations, exit
			if (NumBSFoundation==0) return;
			// If BASESIMP has already been initiated, exit
			if (allocated(BaseSimp)) return;
			
			// Otherwise, allocate the memory and read the variables
			BaseSimp.allocate(NumBSFoundation);
			//BaseSimpReport.allocate(NumBSFoundation);
			CheckEquipName.dimension(NumBSFoundation, true);
			
			// Evaluate number of time steps in a week
			BSNumTimeStepsInWeek = 7 * 24 * NumOfTimeStepInHour;
						
			// Initiate the omega constant
			// Pi comes from the DataGlobals class
			BSOmega = 2.*Pi / (365.*24.);

			// Initiate corner correction factors table
			InitBSCornerCoeff();
			
			//Evaluate coefficients for every BASESIMP foundation -> Loop
			for (BSFoundationNum = 1; BSFoundationNum <= NumBSFoundation; ++BSFoundationNum) {
				// Allocate memory for zone temperatures array and initialize it
				if (!allocated(BaseSimp(BSFoundationNum).BSZoneInstantTemps))		BaseSimp(BSFoundationNum).BSZoneInstantTemps.allocate(BSNumTimeStepsInWeek);
				for (BSTimeStepNum = 1; BSTimeStepNum <= BSNumTimeStepsInWeek; ++BSTimeStepNum) {
					BaseSimp(BSFoundationNum).BSZoneInstantTemps(BSTimeStepNum) = BSTzoneInit;
				}

				// Initialize zone average temperature
				BaseSimp(BSFoundationNum).BSZoneAvgTemp = BSTzoneInit;


				// Read the data from the input file (idf)
				GetBSFoundationInput(BSFoundationNum);

				// Evaluate the BASESIMP factors (Sag, Sbgavg, Sbgvar and phase lag)
				BSFactCtrl(BSFoundationNum);

				// Adjust the phase lag by Pi/2 and the phase to avoid having to do this every time the heat loss calculation is performed
				// Note that phase lag must be evaluated first so this operation must be performed after call to BSFactCtrl()
				BaseSimp(BSFoundationNum).BSPhaseAdj = BaseSimp(BSFoundationNum).BSTGps + BaseSimp(BSFoundationNum).BSPhase  + Pi / 2.;

				
				// Establish link between foundation heat losses and corresponding zone internal gains
				SetupZoneInternalGain(BaseSimp(BSFoundationNum).BSZoneNumber, cCurrentModuleObject, BaseSimp(BSFoundationNum).Name, IntGainTypeOf_Basesimp, BaseSimp(BSFoundationNum).BSIntGainsToZone);

				// Establish output variable for BASESIMP foundation heat losses
				SetupOutputVariable("Total heat losses from foundation [W]", BaseSimp(BSFoundationNum).BSFoundationLosses, "Zone", "Average", BaseSimp(BSFoundationNum).BSFoundationName);
							
			} // End loop through foundations
		} // End routine

		
		void
			GetBSFoundationInput(int const BSFoundationNum)
		{
			// SUBROUTINE INFORMATION:
			//       AUTHOR:          Patrice Pinel
			//       DATE WRITTEN:    March 2016
			//       MODIFIED:        
			//       RE-ENGINEERED:   

			// PURPOSE OF THIS SUBROUTINE:
			// get data for one BASESIMP Foundation from input file

			// METHODOLOGY EMPLOYED:
			// standard EnergyPlus input retrieval using input Processor

			// REFERENCES: na

			//Using/Aliasing
			// using namespace InputProcessor;
			using InputProcessor::GetObjectItem;
			using InputProcessor::SameString;
			using InputProcessor::VerifyName;
			using InputProcessor::FindItemInList;
			using namespace DataIPShortCuts; // Data for field names, blank numerics
			using DataHeatBalance::Zone;
								
			//LOCAL VARIABLES
			static std::string const RoutineName("GetBSFoundationInput: ");
			int NumAlphas;						// Number of elements in the alpha array
			int NumNums;						// Number of elements in the numeric array
			int IOStat;							// IO Status when calling get input subroutine
			static bool ErrorsFound(false);		// Flag to show errors were found during GetInput
			bool IsBlank;						// Flag for blank name
			bool errFlag;						// Flag to show errors were found during function call

			//Identify the EnergyPlus object being read
			cCurrentModuleObject = "Site:GroundDomain:BASESIMP";
						
			//LOAD ARRAY WITH BaseSimp DATA
			GetObjectItem(cCurrentModuleObject, BSFoundationNum, cAlphaArgs, NumAlphas, rNumericArgs, NumNums, IOStat, lNumericFieldBlanks, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames);

			IsBlank = false;
			// Called this was copied for used IsNotOK boolean and then made ErrorsFound true if IsNotOK
			// I bypassed this here and sent ErrorsFound directly to VerifyName
			VerifyName(cAlphaArgs(1), BaseSimp, BSFoundationNum - 1, ErrorsFound, IsBlank, cCurrentModuleObject + " Name");
			
		
			// Read the name of the foundation
			BaseSimp(BSFoundationNum).BSFoundationName = cAlphaArgs(1);

			//Read the name of the zone affected by the foundation
			BaseSimp(BSFoundationNum).BSZoneName = cAlphaArgs(2);
			BaseSimp(BSFoundationNum).BSZoneNumber = FindItemInList(cAlphaArgs(2), Zone);

			// Read numeric value used for both modelling method (BASECALC coef or foundation description)
			BaseSimp(BSFoundationNum).BSExposedFraction = rNumericArgs(1);
			BaseSimp(BSFoundationNum).BSTGps = rNumericArgs(15);

			// Read data related to soil temperature
			BaseSimp(BSFoundationNum).BSTGavg = rNumericArgs(14);
			BaseSimp(BSFoundationNum).BSTGamp = rNumericArgs(15);
			BaseSimp(BSFoundationNum).BSTGps = rNumericArgs(16);

			// Check the model type
			if (SameString(cAlphaArgs(3), "BASECALC coefficients")) {  // BASECALC Coefficients
			
				// Set model type
				BaseSimp(BSFoundationNum).BSModelType = 2;

				// Read Sag, Sbgavg, and Sbgvar
				BaseSimp(BSFoundationNum).BSSag = rNumericArgs(10);
				BaseSimp(BSFoundationNum).BSSbgavg = rNumericArgs(11);
				BaseSimp(BSFoundationNum).BSSbgvar = rNumericArgs(12);
				BaseSimp(BSFoundationNum).BSPhase = rNumericArgs(13);
							

			}  // end model that works with BASECALC coefficients
			else if (SameString(cAlphaArgs(3), "Foundation description"))	{	// Description of the foundation
			
				// Set model type
				BaseSimp(BSFoundationNum).BSModelType = 1;

				// Read the BASESIMP Foundation type, name, and configuration
				// Check the type of foundation
				if (SameString(cAlphaArgs(4), "Basement"))	{

					BaseSimp(BSFoundationNum).BSFndType = 1;			// Basement

					// Check Material
					if (SameString(cAlphaArgs(5), "Concrete"))	{

						BaseSimp(BSFoundationNum).BSFndMat = 1;			// Concrete
						
						// Check configuration
						if (SameString(cAlphaArgs(6), "BCIN_1"))		BaseSimp(BSFoundationNum).BSFndConfig = 1;		// conc Wl&Fl; intrn ins full - height; *any constr
						else if (SameString(cAlphaArgs(6), "BCIN_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 2;		// conc Wl&Fl; intrn ins to 0.2 above floor; *any constr
						else if (SameString(cAlphaArgs(6), "BCIN_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 3;		// conc Wl&Fl; intrn ins to 0.6 below grade; *brick veneer
						else if (SameString(cAlphaArgs(6), "BCIN_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 4;		// conc Wl&Fl; intrn ins to 0.6 below grade; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCEN_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 5;		// conc Wl&Fl; extrn ins full - height; *brick veneer
						else if (SameString(cAlphaArgs(6), "BCEN_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 6;		// conc Wl&Fl; extrn ins full - height; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCEN_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 7;		// conc Wl&Fl; extrn ins below grade; *brick veneer
						else if (SameString(cAlphaArgs(6), "BCEN_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 8;		// conc Wl&Fl; extrn ins below grade; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCNN_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 9;		// conc Wl&Fl; no ins; *brick veneer
						else if (SameString(cAlphaArgs(6), "BCNN_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 10;		// conc Wl&Fl; no ins; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCCN_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 11;		// conc Wl&Fl; extrn below grade with overlap; *brick on wall
						else if (SameString(cAlphaArgs(6), "BCCN_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 12;		// conc Wl&Fl; extrn below grade with overlap; *non - brick
						else if (SameString(cAlphaArgs(6), "BCIB_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 19;		// conc Wl&Fl; intrn ins full; ins sub - surface; *any constr
						else if (SameString(cAlphaArgs(6), "BCIB_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 20;		// conc Wl&Fl; intrn ins full; 0.6m perm floor ins; *any constr
						else if (SameString(cAlphaArgs(6), "BCIB_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 21;		// conc Wl&Fl; intrn ins full; 1.0m perm floor ins; *any constr
						else if (SameString(cAlphaArgs(6), "BCIB_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 22;		// conc Wl&Fl; intrn ins full; full floor ins; *any constr
						else if (SameString(cAlphaArgs(6), "BCIB_5"))	BaseSimp(BSFoundationNum).BSFndConfig = 23;		// conc Wl&Fl; intrn ins full; 0.6m perm floor ins; TB; *any
						else if (SameString(cAlphaArgs(6), "BCIB_6"))	BaseSimp(BSFoundationNum).BSFndConfig = 24;		// conc Wl&Fl; intrn ins full; 1.0m perm floor ins; TB; *any
						else if (SameString(cAlphaArgs(6), "BCEB_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 25;		// conc Wl&Fl; extrn ins full; full floor ins; *any constr
						else if (SameString(cAlphaArgs(6), "BCEB_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 26;		// conc Wl&Fl; extrn ins full; 0.6m perm floor ins; *any constr
						else if (SameString(cAlphaArgs(6), "BCCN_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 68;		// conc Wl&Fl; no slab ins; TB; full ins both sides; *any
						else if (SameString(cAlphaArgs(6), "BCCN_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 69;		// conc Wl&Fl; full slab; TB; full ins both sides; *any constr
						else if (SameString(cAlphaArgs(6), "BCEA_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 70;		// conc Wl&Fl; full top slab; extrn ins full; *brick veneer
						else if (SameString(cAlphaArgs(6), "BCEA_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 71;		// conc Wl&Fl; full top slab; extrn ins full; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCIA_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 72;		// conc Wl&Fl; full top slab; intern ins full; *brick veneer
						else if (SameString(cAlphaArgs(6), "BCIA_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 73;		// conc Wl&Fl; full top slab; intern ins full; *non - brick
						else if (SameString(cAlphaArgs(6), "BCEA_7"))	BaseSimp(BSFoundationNum).BSFndConfig = 74;		// conc Wl&Fl; full top slab; extrn ins blw_grd; *brick on slb
						else if (SameString(cAlphaArgs(6), "BCEA_8"))	BaseSimp(BSFoundationNum).BSFndConfig = 75;		// conc Wl&Fl; full top slab; extrn blw_grd; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCEB_8"))	BaseSimp(BSFoundationNum).BSFndConfig = 76;		// conc Wl&Fl; full slab ins; extrn ins blw_grd; *brick veneer
						else if (SameString(cAlphaArgs(6), "BCEB_9"))	BaseSimp(BSFoundationNum).BSFndConfig = 77;		// conc Wl&Fl; full slab ins; TB; extrn 0.6 blw_grd; *non - brick
						else if (SameString(cAlphaArgs(6), "BCCB_8"))	BaseSimp(BSFoundationNum).BSFndConfig = 92;		// conc Wl&Fl; 0.6m perm slab; intrn&extrn ins full; *any con
						else if (SameString(cAlphaArgs(6), "BCCA_7"))	BaseSimp(BSFoundationNum).BSFndConfig = 93;		// conc Wl&Fl; top slab; intrn full; ext 0.6 blw_grd; *any con
						else if (SameString(cAlphaArgs(6), "BCCA_8"))	BaseSimp(BSFoundationNum).BSFndConfig = 94;		// conc Wl&Fl; top slab; intrn above 0.2; extrn full; *any con
						else if (SameString(cAlphaArgs(6), "BCCN_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 95;		// conc Wl&Fl; intrn full; extrn 0.6 blw_grd; *any constr
						else if (SameString(cAlphaArgs(6), "BCCN_5"))	BaseSimp(BSFoundationNum).BSFndConfig = 96;		// conc Wl&Fl; intrn top to 0.2; extrn ins full; *any constr
						else if (SameString(cAlphaArgs(6), "BCEA_5"))	BaseSimp(BSFoundationNum).BSFndConfig = 97;		// conc Wl&Fl; top slab; extrn 0.6 blw_grd; *brick veneer
						else if (SameString(cAlphaArgs(6), "BCEA_6"))	BaseSimp(BSFoundationNum).BSFndConfig = 98;		// conc Wl&Fl; top slab; extrn 0.6 blw_grd; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCEB_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 99;		// conc Wl&Fl; full slab ins; extrn full; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCEB_5"))	BaseSimp(BSFoundationNum).BSFndConfig = 100;	// conc Wl&Fl; 0.6m perm slab; extrn full; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCEB_6"))	BaseSimp(BSFoundationNum).BSFndConfig = 101;	// conc Wl&Fl; 1.0m perm slab; extrn full; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCEN_5"))	BaseSimp(BSFoundationNum).BSFndConfig = 109;	// conc Wl&Fl; extrn ins to 0.6 blw_grd; *brick on wall
						else if (SameString(cAlphaArgs(6), "BCEN_6"))	BaseSimp(BSFoundationNum).BSFndConfig = 110;	// conc Wl&Fl; extrn ins to 0.6 blw_grd; *non - brick veneer
						else if (SameString(cAlphaArgs(6), "BCCB_9"))	BaseSimp(BSFoundationNum).BSFndConfig = 114;	// conc Wl&Fl; full slab; intrn full; extrn 0.6 blw_grd; *any
						else if (SameString(cAlphaArgs(6), "BCCB_10"))	BaseSimp(BSFoundationNum).BSFndConfig = 115;	// conc Wl&Fl; 0.6m perm; intrn full; extrn 0.6 blw_grd; *any
						else if (SameString(cAlphaArgs(6), "BCCA_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 116;	// conc Wl&Fl; full top; int overlap; extrn blw_grd; *brick
						else if (SameString(cAlphaArgs(6), "BCCA_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 117;	// conc Wl&Fl; full top; int overlap; extrn blw_grd; *non - brick
						else if (SameString(cAlphaArgs(6), "BCIB_7"))	BaseSimp(BSFoundationNum).BSFndConfig = 118;	// conc Wl&Fl; full slab; TB; intrn 0.6 blw_grd; *brick on wall
						else if (SameString(cAlphaArgs(6), "BCIB_8"))	BaseSimp(BSFoundationNum).BSFndConfig = 119;	// conc Wl&Fl; full slab; TB; intrn 0.6 blw_grd; *non - brick
						else if (SameString(cAlphaArgs(6), "BCIA_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 121;	// conc Wl&Fl; 1.0m perm top; intrn full; *brick on wall
						else if (SameString(cAlphaArgs(6), "BCIA_5"))	BaseSimp(BSFoundationNum).BSFndConfig = 122;	// conc Wl&Fl; 0.6m perm top; intrn 0.6 blw_grd; *brick on wall
						else if (SameString(cAlphaArgs(6), "BCIA_6"))	BaseSimp(BSFoundationNum).BSFndConfig = 123;	// conc Wl&Fl; 0.6m perm top; intrn 0.6 blw_grd; *non - brick
						else if (SameString(cAlphaArgs(6), "BCIB_9"))	BaseSimp(BSFoundationNum).BSFndConfig = 124;	// conc Wl&Fl; 0.6m perm; TB; intrn 0.6 blw_grd; *brick on wall
						else if (SameString(cAlphaArgs(6), "BCIB_10"))	BaseSimp(BSFoundationNum).BSFndConfig = 125;	// conc Wl&Fl; 0.6m perm; TB; intrn 0.6 blw_grd; *non - brick
						else if (SameString(cAlphaArgs(6), "BCEB_10"))	BaseSimp(BSFoundationNum).BSFndConfig = 126;	// conc Wl&Fl; 0.6m perm; TB; extrn 0.6 blw_grd; *brick on w
						else if (SameString(cAlphaArgs(6), "BCEB_11"))	BaseSimp(BSFoundationNum).BSFndConfig = 127;	// conc Wl&Fl; 0.6m perm; TB; extrn 0.6 blw_grd; *non - brick
						else if (SameString(cAlphaArgs(6), "BCEA_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 128;	// conc Wl&Fl; 1.0m perm top; extrn full; *brick on wall
						else if (SameString(cAlphaArgs(6), "BCEA_9"))	BaseSimp(BSFoundationNum).BSFndConfig = 129;	// conc Wl&Fl; 1.0m perm top; extrn full; *non - brick
						else if (SameString(cAlphaArgs(6), "BCEA_10"))	BaseSimp(BSFoundationNum).BSFndConfig = 130;	// conc Wl&Fl; 0.6m perm top; extrn 0.6 blw_grd; *brick on w
						else if (SameString(cAlphaArgs(6), "BCEA_11"))	BaseSimp(BSFoundationNum).BSFndConfig = 131;	// conc Wl&Fl; 0.6m perm top; extrn 0.6 blw_grd; *non - brick
						else { // No configuration -> error
						
							ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
							ShowContinueError("Invalid " + cAlphaFieldNames(6) + '=' + cAlphaArgs(6));
							//Assign temporary value to avoid crash
							BaseSimp(BSFoundationNum).BSFndConfig = 1;
							ErrorsFound = true;
						}  // end error
					}	// end concrete basement
					else if (SameString(cAlphaArgs(5), "Wood"))	{

						BaseSimp(BSFoundationNum).BSFndMat = 2;			// Wood
						
						// Check configuration
						if (SameString(cAlphaArgs(7), "BWNN_1"))		BaseSimp(BSFoundationNum).BSFndConfig = 13;		// wood Wl&Fl; no ins; *any constr
						else if (SameString(cAlphaArgs(7), "BWIN_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 14;		// wood Wl&Fl; intrn ins full - height; *any constr
						else if (SameString(cAlphaArgs(7), "BWIN_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 15;		// wood Wl&Fl; intrn ins to 0.2 above floor; *any constr
						else if (SameString(cAlphaArgs(7), "BWIN_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 16;		// wood Wl&Fl; intrn ins to 0.6 below grade; *any constr
						else if (SameString(cAlphaArgs(7), "BWEN_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 17;		// wood Wl&Fl; extrn ins full - height; *any constr
						else if (SameString(cAlphaArgs(7), "BWEN_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 18;		// wood Wl&Fl; extrn ins below grade; *any constr
						else if (SameString(cAlphaArgs(7), "BWEN_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 78;		// wood Wl&Fl; extrn to 0.6 blw_grd; *any constr
						else if (SameString(cAlphaArgs(7), "BWIA_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 102;		// wood Wl&Fl; 0.6m perm top floor; intrn full; *any constr
						else if (SameString(cAlphaArgs(7), "BWIA_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 103;		// wood Wl&Fl; top floor ins; intrn full; *any constr
						else if (SameString(cAlphaArgs(7), "BWIB_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 132;		// wood Wl&Fl; 1.0m perm floor; intrn full; *brick on wall
						else if (SameString(cAlphaArgs(7), "BWIB_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 133;		// wood Wl&Fl; full floor; intrn full; *brick on wall
						else if (SameString(cAlphaArgs(7), "BWIB_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 134;		// wood Wl&Fl; 0.6m perm floor; intrn full; *brick on wall
						else if (SameString(cAlphaArgs(7), "BWIA_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 135;		// wood Wl&Fl; 0.6m perm top; intrn 0.6 blw_grd; *brick on wall
						else if (SameString(cAlphaArgs(7), "BWEB_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 136;		// wood Wl&Fl; 1.0m perm floor; extrn full; *brick on wall
						else if (SameString(cAlphaArgs(7), "BWEB_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 137;		// wood Wl&Fl; full floor; extrn full; *brick on wall
						else if (SameString(cAlphaArgs(7), "BWEB_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 138;		// wood Wl&Fl; 0.6m perm floor; extrn 0.6 blw_grd; *brick on w
						else if (SameString(cAlphaArgs(7), "BWEB_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 139;		// wood Wl&Fl; 0.6m perm floor; extrn full; *brick on wall
						else {// No configuration -> error
							ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
							ShowContinueError("Invalid " + cAlphaFieldNames(7) + '=' + cAlphaArgs(7));
							//Assign temporary value to avoid crash
							BaseSimp(BSFoundationNum).BSFndConfig = 1;
							ErrorsFound = true;
						}  // end error
					}  // end wood basement
					else if (SameString(cAlphaArgs(5), "Concrete & Wood"))	{

						BaseSimp(BSFoundationNum).BSFndMat = 3;			// Wood & Concrete

						// Check configuration
						if (SameString(cAlphaArgs(8), "BBIB_3"))		BaseSimp(BSFoundationNum).BSFndConfig = 79;		// wood Wl conc Fl; full slab; TB; intrn full; *any constr
						else if (SameString(cAlphaArgs(8), "BBIB_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 80;		// wood Wl conc Fl; full slab; TB; intrn 0.6 blw_grd; *brick
						else if (SameString(cAlphaArgs(8), "BBEB_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 87;		// wood Wl conc Fl; full slab; extrn 0.6 blw_grd; *brick on w
						else if (SameString(cAlphaArgs(8), "BBEN_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 88;		// wood Wl conc Fl; no slab; extrn full; *brick on wall
						else if (SameString(cAlphaArgs(8), "BBEN_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 89;		// wood Wl conc Fl; no slab; extrn 0.6 blw_grd; *brick on wal
						else if (SameString(cAlphaArgs(8), "BBIA_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 90;		// wood Wl conc Fl; top slab; intrn ins full; *any constr
						else if (SameString(cAlphaArgs(8), "BBIN_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 91;		// wood Wl conc Fl; intrn to 0.6 blw_grd; *non - brick veneer
						else if (SameString(cAlphaArgs(8), "BBIN_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 108;		// wood Wl conc Fl; intrn ins full; *non - brick veneer
						else if (SameString(cAlphaArgs(8), "BBIA_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 111;		// wood Wl conc Fl; top slab; intrn full; *any constr
						else if (SameString(cAlphaArgs(8), "BBIB_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 112;		// wood Wl conc Fl; 0.6m perm slab; intrn full; *any constr
						else if (SameString(cAlphaArgs(8), "BBIB_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 113;		// wood Wl conc Fl; full slab; intrn full; *any constr
						else if (SameString(cAlphaArgs(8), "BBEB_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 120;		// wood W conc F; full slab; TB; extrn 0.6 blw_grd; *brick on w
						else if (SameString(cAlphaArgs(8), "BBIB_5"))	BaseSimp(BSFoundationNum).BSFndConfig = 140;		// wood Wl conc Fl; 1.0m perm; intrn ins full; *brick on wall
						else if (SameString(cAlphaArgs(8), "BBIB_6"))	BaseSimp(BSFoundationNum).BSFndConfig = 141;		// wood Wl conc Fl; 1.0m perm; intrn 0.6 blw_grd; *brick on w
						else if (SameString(cAlphaArgs(8), "BBEB_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 142;		// wood Wl conc Fl; 1.0m perm; extrn full; *brick on wall
						else if (SameString(cAlphaArgs(8), "BBEB_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 143;		// wood Wl conc Fl; 0.6m perm; extrn full; *brick on wall
						else if (SameString(cAlphaArgs(8), "BBEB_5"))	BaseSimp(BSFoundationNum).BSFndConfig = 144;		// wood Wl conc Fl; 0.6m perm; extrn 0.6 blw_grd; *brick on w
						else if (SameString(cAlphaArgs(8), "BBEA_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 145;		// wood Wl conc Fl; full top slab; extrn full; *brick on wall
						else { // No configuration -> error
							ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
							ShowContinueError("Invalid " + cAlphaFieldNames(8) + '=' + cAlphaArgs(8));
							//Assign temporary value to avoid crash
							BaseSimp(BSFoundationNum).BSFndConfig = 1;
							ErrorsFound = true;
						}  // end error

					} // end concrete and wood basement
					else { // No material -> error
						ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
						ShowContinueError("Invalid " + cAlphaFieldNames(5) + '=' + cAlphaArgs(5));
						//Assign temporary value to avoid crash
						BaseSimp(BSFoundationNum).BSFndMat = 1;
						ErrorsFound = true;
					}
				}  // end basement
				else if (SameString(cAlphaArgs(4), "Slab-in-grade")) {
					BaseSimp(BSFoundationNum).BSFndType = 2;				// Slab
					// Check Material
					if (SameString(cAlphaArgs(5), "Concrete")) {
						BaseSimp(BSFoundationNum).BSFndMat = 1;			// Concrete

						// Check configuration
						if (SameString(cAlphaArgs(9), "SCN_1"))			BaseSimp(BSFoundationNum).BSFndConfig = 28;		// conc / soil FL; no ins; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCN_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 29;		// conc / soil FL; no ins; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCN_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 30;		// conc / soil FL; no ins; TB; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCN_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 31;		// conc / soil FL; no ins; TB; *brick veneer
						else if (SameString(cAlphaArgs(9), "SCN_7"))	BaseSimp(BSFoundationNum).BSFndConfig = 32;		// conc / soil FL; no ins; vert skirt; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCN_8"))	BaseSimp(BSFoundationNum).BSFndConfig = 33;		// conc / soil FL; no ins; vert skirt; TB; *brick veneer
						else if (SameString(cAlphaArgs(9), "SCB_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 34;		// conc / soil FL; 0.6m perm slab ins; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCB_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 35;		// conc / soil FL; 0.6m perm slab ins; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_3"))	BaseSimp(BSFoundationNum).BSFndConfig = 36;		// conc / soil FL; 0.6m perm slab & footing ins; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCB_4"))	BaseSimp(BSFoundationNum).BSFndConfig = 37;		// conc / soil FL; 0.6m perm slab & footing ins; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_5"))	BaseSimp(BSFoundationNum).BSFndConfig = 38;		// conc / soil FL; no slab or footing ins; TB; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCB_6"))	BaseSimp(BSFoundationNum).BSFndConfig = 39;		// conc / soil FL; no slab or footing ins; TB; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_9"))	BaseSimp(BSFoundationNum).BSFndConfig = 40;		// conc / soil FL; 0.6m perm slab ins; TB; vert skirt; *non - brick
						else if (SameString(cAlphaArgs(9), "SCB_10"))	BaseSimp(BSFoundationNum).BSFndConfig = 41;		// conc / soil FL; 0.6m perm slab; TB; vert skirt;  *brick on slb
						else if (SameString(cAlphaArgs(9), "SCB_11"))	BaseSimp(BSFoundationNum).BSFndConfig = 42;		// conc / soil FL; 0.6m perm slab; horiz skirt; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCB_12"))	BaseSimp(BSFoundationNum).BSFndConfig = 43;		// conc / soil FL; 0.6m perm slab; horiz skirt; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_13"))	BaseSimp(BSFoundationNum).BSFndConfig = 44;		// conc / soil FL; 1.0m perm slab; horiz skirt; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_14"))	BaseSimp(BSFoundationNum).BSFndConfig = 45;		// conc / soil FL; 1.0m perm slab; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_17"))	BaseSimp(BSFoundationNum).BSFndConfig = 46;		// conc / soil FL; 1.0m perm slab; TB; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCB_18"))	BaseSimp(BSFoundationNum).BSFndConfig = 47;		// conc / soil FL; 1.0m perm slab; TB; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_21"))	BaseSimp(BSFoundationNum).BSFndConfig = 48;		// conc / soil FL; 1.0m perm slab; TB; vert skirt; *non - brick
						else if (SameString(cAlphaArgs(9), "SCB_22"))	BaseSimp(BSFoundationNum).BSFndConfig = 49;		// conc / soil FL; 1.0m perm slab; TB; vert skirt; *brick on slb
						else if (SameString(cAlphaArgs(9), "SCB_23"))	BaseSimp(BSFoundationNum).BSFndConfig = 50;		// conc / soil FL; 1.0m perm slab; TB; horiz skirt; *non - brick
						else if (SameString(cAlphaArgs(9), "SCB_24"))	BaseSimp(BSFoundationNum).BSFndConfig = 51;		// conc / soil FL; 1.0m perm slab; TB; horiz skirt; *brick on slb
						else if (SameString(cAlphaArgs(9), "SCB_25"))	BaseSimp(BSFoundationNum).BSFndConfig = 52;		// conc / soil FL; full slab ins; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCB_26"))	BaseSimp(BSFoundationNum).BSFndConfig = 53;		// conc / soil FL; full slab ins; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_29"))	BaseSimp(BSFoundationNum).BSFndConfig = 54;		// conc / soil FL; full slab ins; TB; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCB_30"))	BaseSimp(BSFoundationNum).BSFndConfig = 55;		// conc / soil FL; full slab ins; TB; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_33"))	BaseSimp(BSFoundationNum).BSFndConfig = 56;		// conc / soil FL; full slab ins; TB; vert skirt; *non - brick
						else if (SameString(cAlphaArgs(9), "SCB_34"))	BaseSimp(BSFoundationNum).BSFndConfig = 57;		// conc / soil FL; full slab ins; TB; vert skirt; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_35"))	BaseSimp(BSFoundationNum).BSFndConfig = 58;		// conc / soil FL; full slab ins; horiz skirt; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCB_36"))	BaseSimp(BSFoundationNum).BSFndConfig = 59;		// conc / soil FL; full slab ins; horiz skirt; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCA_17"))	BaseSimp(BSFoundationNum).BSFndConfig = 60;		// conc / soil FL; full top slab; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCA_18"))	BaseSimp(BSFoundationNum).BSFndConfig = 61;		// conc / soil FL; full top slab; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCA_19"))	BaseSimp(BSFoundationNum).BSFndConfig = 62;		// conc / soil FL; full top slab; TB; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCA_20"))	BaseSimp(BSFoundationNum).BSFndConfig = 63;		// conc / soil FL; full top slab; TB; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCA_21"))	BaseSimp(BSFoundationNum).BSFndConfig = 64;		// conc / soil FL; full top slab; TB; vert skirt; *non - brick venr
						else if (SameString(cAlphaArgs(9), "SCA_22"))	BaseSimp(BSFoundationNum).BSFndConfig = 65;		// conc / soil FL; full top slab; TB; vert skirt; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCA_23"))	BaseSimp(BSFoundationNum).BSFndConfig = 66;		// conc / soil FL; full top slab; TB; horiz skirt; *non - brick
						else if (SameString(cAlphaArgs(9), "SCA_24"))	BaseSimp(BSFoundationNum).BSFndConfig = 67;		// conc / soil FL; full top slab; TB; horiz skirt; brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_31"))	BaseSimp(BSFoundationNum).BSFndConfig = 81;		// conc / soil FL; full slab & footings; TB; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCB_32"))	BaseSimp(BSFoundationNum).BSFndConfig = 82;		// conc / soil FL; full slab & footings; TB; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCB_37"))	BaseSimp(BSFoundationNum).BSFndConfig = 83;		// conc / soil FL; 0.35m perm slab & footings; TB; *non - brick
						else if (SameString(cAlphaArgs(9), "SCB_38"))	BaseSimp(BSFoundationNum).BSFndConfig = 84;		// conc / soil FL; 0.35m perm slab & footings; TB; *brick on slb
						else if (SameString(cAlphaArgs(9), "SCB_39"))	BaseSimp(BSFoundationNum).BSFndConfig = 85;		// conc / soil FL; 0.75m perm slab & footings; TB; *non - brick
						else if (SameString(cAlphaArgs(9), "SCB_40"))	BaseSimp(BSFoundationNum).BSFndConfig = 86;		// conc / soil FL; 0.75m perm slab & footings; TB; *brick on slb
						else if (SameString(cAlphaArgs(9), "SCA_1"))	BaseSimp(BSFoundationNum).BSFndConfig = 104;	// conc / soil FL; 0.6m perm top slab ins; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCA_2"))	BaseSimp(BSFoundationNum).BSFndConfig = 105;	// conc / soil FL; 0.6m perm top slab ins; *brick on slab
						else if (SameString(cAlphaArgs(9), "SCA_9"))	BaseSimp(BSFoundationNum).BSFndConfig = 106;	// conc / soil FL; 1.0m perm top slab ins; *non - brick veneer
						else if (SameString(cAlphaArgs(9), "SCA_10"))	BaseSimp(BSFoundationNum).BSFndConfig = 107;	// conc / soil FL; 1.0m perm top slab ins; *brick on slab
						
						else { // No configuration -> error
							ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
							ShowContinueError("Invalid " + cAlphaFieldNames(9) + '=' + cAlphaArgs(9));
							//Assign temporary value to avoid crash
							BaseSimp(BSFoundationNum).BSFndConfig = 1;
							ErrorsFound = true;
						}  // end error

					}  // end concrete slab
					else { // Slabs can only be concrete -> error
						ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
						ShowContinueError("Invalid " + cAlphaFieldNames(5) + '=' + cAlphaArgs(5));
						//Assign temporary value to avoid crash
						BaseSimp(BSFoundationNum).BSFndMat = 1;
						ErrorsFound = true;
					}
				}  // end foundation type (basement/slab)
				else { // No valid foundation type -> Error
					ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
					ShowContinueError("Invalid " + cAlphaFieldNames(4) + '=' + cAlphaArgs(4));
					//Assign temporary value to avoid crash
					BaseSimp(BSFoundationNum).BSFndType = 1;
					ErrorsFound = true;
				}
				
				// Read rest of foundation description
				BaseSimp(BSFoundationNum).BSHeight = rNumericArgs(2);
				BaseSimp(BSFoundationNum).BSDepth = rNumericArgs(3);
				BaseSimp(BSFoundationNum).BSLength = rNumericArgs(4);
				BaseSimp(BSFoundationNum).BSWidth = rNumericArgs(5);
				BaseSimp(BSFoundationNum).BSOverlap = rNumericArgs(6);
				BaseSimp(BSFoundationNum).BSRSI = rNumericArgs(7);
				BaseSimp(BSFoundationNum).BSSoilK = rNumericArgs(8);
				BaseSimp(BSFoundationNum).BSWTD = rNumericArgs(9);
				
				// At the moment I am assuming that the absolute values of the numerical inputs
				// are automatically checked by the input processor from the limits indicated in the dictionary
				// However, a few limits that are relative to the values other parameters can not be specified in the dictionary
				// So they are checked manually here
					
				// Verify that the height is larger than the depth
				if ((BaseSimp(BSFoundationNum).BSHeight - BaseSimp(BSFoundationNum).BSDepth) < 0.1)	{
					ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
					ShowContinueError("Height must be 0.1m greater than depth");
					ErrorsFound = true;
				}
					
				// Verify the depth limit for basements
				if ((BaseSimp(BSFoundationNum).BSDepth < 0.65) & (BaseSimp(BSFoundationNum).BSFndType == 1)) {
				// Basement with depth < 0.65 m
					ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
					ShowContinueError("Basement must have a depth larger than 0.65 m");
					ErrorsFound = true;
				}
				// Verify the depth limit for slabs
				if ((BaseSimp(BSFoundationNum).BSDepth > 0.05) & (BaseSimp(BSFoundationNum).BSFndType == 2)) {
				// Slab with depth > 0.05 m
					ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
					ShowContinueError("Slab-in-grade foundations must have a depth smaller than 0.05 m");
					ErrorsFound = true;
				}
				// Verify that the width is not larger than the length
				if (BaseSimp(BSFoundationNum).BSWidth > BaseSimp(BSFoundationNum).BSLength) 	{
					ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
					ShowContinueError("Width of foundations can not be larger than their length");
					ErrorsFound = true;
				}
				// Check that "combination" overlap is not greater than depth
				if (BaseSimp(BSFoundationNum).BSOverlap > BaseSimp(BSFoundationNum).BSDepth) {
					ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
					ShowContinueError("Overlap of insulation can not be larger than foundation depth");
					ErrorsFound = true;
				}
					
			}	// End of case where model is based on description of foundation
			else { // Error with the model type description
			
				ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs(1) + "\",");
				ShowContinueError("Invalid " + cAlphaFieldNames(3) + '=' + cAlphaArgs(3));
				//Assign temporary value to avoid crash
				BaseSimp(BSFoundationNum).BSModelType = 1;
				ErrorsFound = true;
			}
			

			if (ErrorsFound) 
			{
				ShowFatalError(RoutineName + "Errors found in processing " + cCurrentModuleObject + " input.");
			}
		}		// end of subroutine GetBSFoundationInput()
		

		void
			BSFactCtrl(int const BSFoundationNum)
		{

			// SUBROUTINE INFORMATION:
			//       AUTHOR         Patrice Pinel
			//       DATE WRITTEN   March. 2016
			//       MODIFIED       
			//       RE-ENGINEERED  na

			// PURPOSE OF THIS SUBROUTINE:
			// This subrountine calculates the BASESIMP/BASECALC heat loss factors

			// METHODOLOGY EMPLOYED: 
			// This routine reproduces the BSFCCTL ESP-r routine from the basesimp.F source file

			// REFERENCES: ESP-r source code
			// This subroutine controls the calcualtion of the BASESIMP heat loss
			// factors for the foundation under consideration.It makes appropriate
			// calls to establish the BASESIMP and corner - correction coefficients
			// and to calculate the four heat loss factors.

			// If the BASESIMP configuration is insulated and if the insulation level
			// is outside the range of the BASESIMP correlations(ie.RSI < 1.5), then the
			// "exponential-interpolation" method must be applied to calculate the
			// BASESIMP heat loss factors.This requires two passes at setting coefficients
			// and calculating the heat loss factors.Only a single pass is required when
			// the foundation is uninsulated or when the insulation is within the RSI range.
			// intflg controls the application of the "exponential-interpolation" method:

			//Using/Aliasing
			using InputProcessor::GetNumObjectsFound;
			using namespace DataIPShortCuts; // Data for field names, blank numerics

			//LOCAL VARIABLES
			int FndConfig;						// Configuration number
			
			//Identify the EnergyPlus object being read
			cCurrentModuleObject = "Site:GroundDomain:BASESIMP";

			// Removed code for models with BASESIMP coefficients and old input files
			// This was not pertinent to this model
						
			// Do calculations only if the model type is a Foundation Description
			// Not required if coefficients supplied by user
			if (BaseSimp(BSFoundationNum).BSModelType == 1)	{

				// Initialise the boolean variables defining the foundation type
				// Recognised BASESIMP configurations with no insulation
				bool IsNotInsulated = (BaseSimp(BSFoundationNum).BSFndConfig == 9) || (BaseSimp(BSFoundationNum).BSFndConfig == 10)
					|| (BaseSimp(BSFoundationNum).BSFndConfig == 13) || (BaseSimp(BSFoundationNum).BSFndConfig == 28)
					|| (BaseSimp(BSFoundationNum).BSFndConfig == 29);  
				// Recognised Concrete basement
				// Note: Replaced the ESP-r code that checked all configurations by this more compact and efficient form
				bool IsConcreteBasement = (BaseSimp(BSFoundationNum).BSFndType == 1) && (BaseSimp(BSFoundationNum).BSFndMat == 1);

				// Get foundation configuration
				FndConfig = BaseSimp(BSFoundationNum).BSFndConfig;

				if ((!IsNotInsulated) && (BaseSimp(BSFoundationNum).BSRSI < 1.5)) {
					// Insulated but insulation lower than BASESIMP limit -> Exponential interpolation
				
					// ******************   Step 1: Evaluate factors for RSI=1.5
					Real64 TempRSI = BaseSimp(BSFoundationNum).BSRSI;	// Variable to temporarily hold the RSI value

					BaseSimp(BSFoundationNum).BSRSI = 1.5;				// Change the RSI value for the foundation to 1.5

					// Initiate intermediate coefficients
					GetBSCoeff(BSFoundationNum, FndConfig);

					CalcBSFactors(BSFoundationNum, FndConfig);

					// Transfer the calculated factors to temporary variables
					Real64 BSSag15 = BaseSimp(BSFoundationNum).BSSag;
					Real64 BSSbgavg15 = BaseSimp(BSFoundationNum).BSSbgavg;
					Real64 BSSbgvar15 = BaseSimp(BSFoundationNum).BSSbgvar;
					Real64 BSPhase15 = BaseSimp(BSFoundationNum).BSPhase;

					// Replace the original RSI
					BaseSimp(BSFoundationNum).BSRSI = TempRSI;

					//********************  Step 2: Evaluate factors for no insulation
					int TmpCfg;					// Temporary configuration
					// Concrete basement -> Calculate for BCNN_1
					if (IsConcreteBasement)	TmpCfg = 9;
					// Wood or mixed basement -> Calculate for BWNN_1
					else if (BaseSimp(BSFoundationNum).BSFndType == 1)	TmpCfg = 13;
					// Slab -> Calculate for SCN_1
					else TmpCfg = 28;

					// Evaluate the intermediate coefficients for this uninsulated configuration
					GetBSCoeff(BSFoundationNum, TmpCfg);

					// Evaluate new factors for this uninsulated configuration
					CalcBSFactors(BSFoundationNum, TmpCfg);

					// Transfer the calculated factors to temporary variables
					Real64 BSSag0 = BaseSimp(BSFoundationNum).BSSag;
					Real64 BSSbgavg0 = BaseSimp(BSFoundationNum).BSSbgavg;
					Real64 BSSbgvar0 = BaseSimp(BSFoundationNum).BSSbgvar;
					Real64 BSPhase0 = BaseSimp(BSFoundationNum).BSPhase;

					// ***************** Step 3: Perform interpolation
					
					// Assign Interpolation factor
					Real64 WInt;		
					//Basement
					if (BaseSimp(BSFoundationNum).BSFndType == 1) WInt = 2.29;
					// Slab
					else WInt = 1.77;

					// Interpolate
					BaseSimp(BSFoundationNum).BSSag = BSSag15 + (BSSag0 - BSSag15) / (exp(WInt*TempRSI));
					BaseSimp(BSFoundationNum).BSSbgavg = BSSbgavg15 + (BSSbgavg0 - BSSbgavg15) / (exp(WInt*TempRSI));
					BaseSimp(BSFoundationNum).BSSbgvar = BSSbgvar15 + (BSSbgvar0 - BSSbgvar15) / (exp(WInt*TempRSI));
					BaseSimp(BSFoundationNum).BSPhase = BSPhase15 + (BSPhase0 - BSPhase15) / (exp(WInt*TempRSI));
				
				}	// if insulation < 1.5
				else {
				// Everything normal, just calculate the factors normally
				
					// Assign the intermediate coefficients for this configuration
					GetBSCoeff(BSFoundationNum, FndConfig);

					// Evaluate new factors for this configuration
					CalcBSFactors(BSFoundationNum, FndConfig);
				
				}  // else from insulation < 1.5
			
			} // If model is a description of the foundation

			// Do nothing if model is the BASECALC factors
						
		}	// Routine


		
		void
		CalcBSFactors(int const BSFoundationNum, int const BSFoundationCfg)
		{
			// SUBROUTINE INFORMATION:
			//       AUTHOR         Patrice Pinel
			//       DATE WRITTEN   March. 2016
			//       MODIFIED       
			//       RE-ENGINEERED  na

			// PURPOSE OF THIS SUBROUTINE:
			// This subrountine calculates the BASESIMP/BASECALC heat loss factors

			// METHODOLOGY EMPLOYED: 
			// This routine reproduces the BSFACS ESP-r routine from the basesimp.F source file
			
			// REFERENCES: ESP-r source code

			using InputProcessor::GetNumObjectsFound;
			using namespace DataIPShortCuts; // Data for field names, blank numerics
						
			Real64	RPart1, RPart2, RPart3, RPart4;		// Parts of the heat transfer factors
			Real64	Rr1, Rr2, Rr3, Rr4;					// Parts of the corner correction factors
			Real64	Sumuo, Sumur;
			Real64	HminusD;		// Height minus depth
						
			// Check if user input rsi = 0.  If so, set rsi to a small value to
			// avoid `divide by zero' values in the sumuo, sumur, atten, and phase
			// correlations. A user input of rsi = 0 signifies an uninsulated
			// foundation. The correlation coefficients for uninsulated foundations
			// (eg.BCNN_1, SCN_1) nullify the rsi input(ie.any_value^0. = 1.).
			// In other words the correlations are insensitive to rsi. However,
			// some compilers result in 0.^0. = 0., which causes `divide by zero'
			// problems.This simply avoids this problem without affecting the results.
			if (BaseSimp(BSFoundationNum).BSRSI < 0.01)	{
				BaseSimp(BSFoundationNum).BSRSI = 0.01;
			}

			// Set icol if the BASESIMP system is BCCN_1 or BCCN_2.These two systems
			// do not have fixed icol as they do not correspond to any of the eight
			// modelled for the corner - correction method.icol is 4 if there is less
			// than 0.6m of overlap; it is 5 if there is greater than 0.6m of overlap
			// and the exterior coverage is greater than the interior coverage; and
			// it is 3 if the interior coverage is greater than the exterior coverage.
			
			HminusD = BaseSimp(BSFoundationNum).BSHeight - BaseSimp(BSFoundationNum).BSDepth;
			
			if (ICol == 99)	{
				
				Real64 Wilen = HminusD + BaseSimp(BSFoundationNum).BSOverlap;
				Real64 Welen = 0.1 + BaseSimp(BSFoundationNum).BSDepth;
				
				if ((BaseSimp(BSFoundationNum).BSOverlap / 0.6) > 0.9999)	{
					if ((Welen / Wilen) > 1.0)	ICol = 5;
					else ICol = 3;
				}
				else ICol = 4;
			}

			// Calculate sumuo.
			RPart1 = (BSa1 + BSb1*HminusD + BSc1 / BaseSimp(BSFoundationNum).BSSoilK) / pow(BaseSimp(BSFoundationNum).BSRSI, BSd1);
			RPart2 = 1 / (BSe1 + BSi1*pow(BaseSimp(BSFoundationNum).BSOverlap, BSf1)*pow(BaseSimp(BSFoundationNum).BSRSI, BSg1)*pow(HminusD, BSh1));
			RPart3 = BSj1;
			Sumuo = RPart1*RPart2 + RPart3;

			// Calculate Sag.
			BaseSimp(BSFoundationNum).BSSag = Sumuo*2.*(BaseSimp(BSFoundationNum).BSLength + BaseSimp(BSFoundationNum).BSWidth);

			// Calculate sumur.
			RPart1 = (BSq2 + BSr2*BaseSimp(BSFoundationNum).BSWidth)*(BSu2 + BSv2*BaseSimp(BSFoundationNum).BSSoilK)*(BSw2 + BSx2*BaseSimp(BSFoundationNum).BSDepth);
			RPart2 = pow(BaseSimp(BSFoundationNum).BSWTD, BSs2 + BSt2*BaseSimp(BSFoundationNum).BSWidth + BSy2*BaseSimp(BSFoundationNum).BSDepth);
			RPart3 = BSa2*pow(BaseSimp(BSFoundationNum).BSDepth, BSb2)*pow(BaseSimp(BSFoundationNum).BSSoilK, BSc2);
			RPart4 = pow(BaseSimp(BSFoundationNum).BSWTD, BSd2)*pow(BaseSimp(BSFoundationNum).BSRSI, BSe2 + BSf2*BaseSimp(BSFoundationNum).BSSoilK + BSg2*BaseSimp(BSFoundationNum).BSDepth + BSh2*BaseSimp(BSFoundationNum).BSOverlap);
			Sumur = (RPart1 / RPart2) + (RPart3 / RPart4);

			// Calculate the steady corner factor.
			// Create dummy properties
			Real64 Soil = BaseSimp(BSFoundationNum).BSSoilK;
			Real64 Dept = BaseSimp(BSFoundationNum).BSDepth;
			Real64 Wtabl = BaseSimp(BSFoundationNum).BSWTD;
			Real64 Rs = BaseSimp(BSFoundationNum).BSRSI;
			Real64 Widt = BaseSimp(BSFoundationNum).BSWidth;

			// Set limits on these dummy properties
			if (BaseSimp(BSFoundationNum).BSRSI > 5.)		Rs = 5.;
			if (BaseSimp(BSFoundationNum).BSWidth > 10.)	Widt = 10.;
			if (BaseSimp(BSFoundationNum).BSDepth > 2.)		Dept = 2.;

			Real64 Wby2 = Widt / 2.;
			
			if (ICol == 98)	{
				ICol = 3;
				Rs = 0.;
			}
			
			int IUse = 2 * (ICol - 1) + 1;	// Column to use when looking for corner correction factors in the array
			
			// Evaluate average conrner correction factor
			Rr1 = BSCornerCoeff[IUse][1] + BSCornerCoeff[IUse][2] * Rs + BSCornerCoeff[IUse][3] * Soil
				+ BSCornerCoeff[IUse][4] * Wby2 + BSCornerCoeff[IUse][5] * Dept + BSCornerCoeff[IUse][6] * Wtabl;
			Rr2 = BSCornerCoeff[IUse][7] * pow(Rs, 2) + BSCornerCoeff[IUse][8] * Soil * Rs
				+ BSCornerCoeff[IUse][9] * Wby2 * Rs + BSCornerCoeff[IUse][10] * Wby2 * Soil
				+ BSCornerCoeff[IUse][11] * pow(Wby2, 2);
			Rr3 = BSCornerCoeff[IUse][12] * Dept * Rs + BSCornerCoeff[IUse][13] * Dept * Soil +
				BSCornerCoeff[IUse][14] * Dept * Wby2 + BSCornerCoeff[IUse][15] * pow(Dept, 2);
			Rr4 = BSCornerCoeff[IUse][16] * Wtabl * Rs + BSCornerCoeff[IUse][17] * Wtabl * Soil +
				BSCornerCoeff[IUse][18] * Wtabl * Wby2 + BSCornerCoeff[IUse][19] * Wtabl * Dept;
			Real64 Fcs = Rr1 + Rr2 + Rr3 + Rr4;
			
			// Calculate Sbgavg.
			BaseSimp(BSFoundationNum).BSSbgavg = Sumur*(2.*(BaseSimp(BSFoundationNum).BSLength - BaseSimp(BSFoundationNum).BSWidth) + 4.*Fcs*BaseSimp(BSFoundationNum).BSWidth);

			// Calculate Atten.
			RPart1 = BSa3 + BSb3*BaseSimp(BSFoundationNum).BSSoilK + BSc3*BaseSimp(BSFoundationNum).BSDepth;
			RPart2 = BSe3 + BSf3*BaseSimp(BSFoundationNum).BSSoilK + BSg3*BaseSimp(BSFoundationNum).BSDepth;
			RPart3 = pow(BaseSimp(BSFoundationNum).BSRSI, BSh3 + BSi3*BaseSimp(BSFoundationNum).BSOverlap);
			Real64 Atten = RPart1 + RPart2 / RPart3;

			// Calculate the variable corner factor.
			IUse = 2 * (ICol - 1) + 2;
			Rr1 = BSCornerCoeff[IUse][1] + BSCornerCoeff[IUse][2] * Rs +
				BSCornerCoeff[IUse][3] * Soil + BSCornerCoeff[IUse][4] * Wby2 +
				BSCornerCoeff[IUse][5] * Dept + BSCornerCoeff[IUse][6] * Wtabl;
			Rr2 = BSCornerCoeff[IUse][7] * pow(Rs, 2) + BSCornerCoeff[IUse][8] * Soil * Rs +
				BSCornerCoeff[IUse][9] * Wby2 * Rs + BSCornerCoeff[IUse][10] * Wby2 * Soil +
				BSCornerCoeff[IUse][11] * pow(Wby2, 2);
			Rr3 = BSCornerCoeff[IUse][12] * Dept * Rs + BSCornerCoeff[IUse][13] * Dept * Soil +
				BSCornerCoeff[IUse][14] * Dept * Wby2 + BSCornerCoeff[IUse][15] * pow(Dept, 2);
			Rr4 = BSCornerCoeff[IUse][16] * Wtabl * Rs + BSCornerCoeff[IUse][17] * Wtabl * Soil +
				BSCornerCoeff[IUse][18] * Wtabl * Wby2 + BSCornerCoeff[IUse][19] * Wtabl * Dept;
			Real64 Fcv = Rr1 + Rr2 + Rr3 + Rr4;

			// Calculate Sbgvar.
			BaseSimp(BSFoundationNum).BSSbgvar = Atten*(2.*(BaseSimp(BSFoundationNum).BSLength - BaseSimp(BSFoundationNum).BSWidth) + 4. * BaseSimp(BSFoundationNum).BSWidth * Fcv);

			// Calculate phase.
			BaseSimp(BSFoundationNum).BSPhase = BSa4 + BSb4 / pow(BaseSimp(BSFoundationNum).BSRSI, BSc4);
		
		} // Routine

		
		void
			GetBSCoeff(int const BSFoundationNum, int const BSFoundationCfg)
		{
			// SUBROUTINE INFORMATION:
			//       AUTHOR         Patrice Pinel
			//       DATE WRITTEN   March. 2016
			//       MODIFIED       
			//       RE-ENGINEERED  na

			// PURPOSE OF THIS SUBROUTINE:
			// This subrountine assigns the BASESIMP coefficients depending on configuration

			// METHODOLOGY EMPLOYED: 
			// This routine merges the following ESP-r routines from the bscoeff_extended.F source file:
			// COEFEXT1     BASESIMP coefficients for remaining configs from original set of 67.
			// COEFEXT2     BASESIMP coeffs for configs generated by Debra Haltrect, 1998.
			// COEFEXT3     BASESIMP coeffs for configs generated by Julia Purdy, summer 1999.
			// COEFEXT4     BASESIMP coeffs for configs generated by Kamel Haddad, summer 1999.
			// COEFEXT5		BASESIMP coeffs for configs generated by Julia Purdy, October 1999.
			// As well as the following routine from the bscoeff.F source file:
			// COEFBAS		BASESIMP coeffs for basic configurations.

			// REFERENCES: ESP-r source code

			using namespace DataIPShortCuts; // Data for field names, blank numerics
			
			//LOCAL VARIABLES
			// Variable indicating if configuration has been found
			bool BSFndConfigFound(true);

			
			if (BSFoundationCfg < 50) {
			// Configs cut by blocks of 50 to avoid nested too deeply -> compiler error
				if (BSFoundationCfg == 1)	{
					BSa1 = 0.021400;
					BSb1 = 0.706000;
					BSc1 = 0.102000;
					BSd1 = 0.704000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.291000;
					BSr2 = 0.318000;
					BSu2 = 0.229000;
					BSv2 = 0.620000;
					BSw2 = 0.711000;
					BSx2 = 0.500000;
					BSs2 = -0.055000;
					BSt2 = 0.023500;
					BSy2 = 0.179000;
					BSa2 = 0.749000;
					BSb2 = 0.712000;
					BSc2 = 0.452000;
					BSd2 = 0.263000;
					BSe2 = 3.000000;
					BSf2 = -0.035300;
					BSg2 = -1.011000;
					BSh2 = 0.000000;
					BSa3 = 0.258000;
					BSb3 = 0.317000;
					BSc3 = -0.188000;
					BSe3 = 0.006760;
					BSf3 = 0.110000;
					BSg3 = 0.219000;
					BSh3 = 0.769000;
					BSi3 = 0.000000;
					BSa4 = 2.415000;
					BSb4 = 0.488000;
					BSc4 = 0.259000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 2)	{
					BSa1 = 0.022500;
					BSb1 = 0.698000;
					BSc1 = 0.117000;
					BSd1 = 0.643000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.439000;
					BSr2 = 0.281000;
					BSu2 = 0.266000;
					BSv2 = 0.620000;
					BSw2 = 0.826000;
					BSx2 = 0.500000;
					BSs2 = -0.072000;
					BSt2 = 0.023000;
					BSy2 = 0.178000;
					BSa2 = 0.706000;
					BSb2 = 0.865000;
					BSc2 = 0.533000;
					BSd2 = 0.358000;
					BSe2 = 3.486000;
					BSf2 = -0.064700;
					BSg2 = -1.211000;
					BSh2 = 0.000000;
					BSa3 = 0.439000;
					BSb3 = 0.377000;
					BSc3 = -0.278000;
					BSe3 = -0.045700;
					BSf3 = 0.104000;
					BSg3 = 0.227000;
					BSh3 = 0.756000;
					BSi3 = 0.000000;
					BSa4 = 2.645000;
					BSb4 = 0.284000;
					BSc4 = 0.325000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 3) {
					BSa1 = -0.102000;
					BSb1 = 0.735000;
					BSc1 = 0.133000;
					BSd1 = 0.764000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.095200;
					BSq2 = 1.253000;
					BSr2 = 0.546000;
					BSu2 = 0.256000;
					BSv2 = 0.620000;
					BSw2 = 0.177000;
					BSx2 = 0.500000;
					BSs2 = -0.061700;
					BSt2 = 0.021900;
					BSy2 = 0.187000;
					BSa2 = 0.104000;
					BSb2 = 1.350000;
					BSc2 = 0.765000;
					BSd2 = -0.233000;
					BSe2 = 0.769000;
					BSf2 = -0.013000;
					BSg2 = -0.283000;
					BSh2 = 0.000000;
					BSa3 = -0.218000;
					BSb3 = 0.609000;
					BSc3 = 0.270000;
					BSe3 = 0.180000;
					BSf3 = 0.064600;
					BSg3 = -0.030600;
					BSh3 = 0.812000;
					BSi3 = 0.000000;
					BSa4 = 2.978000;
					BSb4 = -0.008520;
					BSc4 = -0.026600;
					ICol = 2;
				}
				else if (BSFoundationCfg == 4)	{
					BSa1 = -0.107000;
					BSb1 = 0.756000;
					BSc1 = 0.120000;
					BSd1 = 0.727000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.080800;
					BSq2 = 1.217000;
					BSr2 = 0.533000;
					BSu2 = 0.260000;
					BSv2 = 0.620000;
					BSw2 = 0.184000;
					BSx2 = 0.500000;
					BSs2 = -0.072900;
					BSt2 = 0.022100;
					BSy2 = 0.188000;
					BSa2 = 0.111000;
					BSb2 = 1.446000;
					BSc2 = 0.763000;
					BSd2 = -0.180000;
					BSe2 = 0.810000;
					BSf2 = -0.014100;
					BSg2 = -0.297000;
					BSh2 = 0.000000;
					BSa3 = -0.216000;
					BSb3 = 0.608000;
					BSc3 = 0.273000;
					BSe3 = 0.182000;
					BSf3 = 0.064900;
					BSg3 = -0.029500;
					BSh3 = 0.811000;
					BSi3 = 0.000000;
					BSa4 = 3.052000;
					BSb4 = -0.092800;
					BSc4 = -0.135000;
					ICol = 2;
				}
				else if (BSFoundationCfg == 5)	{
					BSa1 = -0.083500;
					BSb1 = 0.749000;
					BSc1 = 0.003340;
					BSd1 = 0.885000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.357000;
					BSq2 = 0.312000;
					BSr2 = 0.349000;
					BSu2 = 0.186000;
					BSv2 = 0.620000;
					BSw2 = 0.678000;
					BSx2 = 0.500000;
					BSs2 = -0.036000;
					BSt2 = 0.023200;
					BSy2 = 0.184000;
					BSa2 = 0.824000;
					BSb2 = 0.722000;
					BSc2 = 0.358000;
					BSd2 = 0.271000;
					BSe2 = 2.917000;
					BSf2 = -0.061300;
					BSg2 = -0.967000;
					BSh2 = 0.000000;
					BSa3 = 0.136000;
					BSb3 = 0.330000;
					BSc3 = -0.181000;
					BSe3 = 0.043400;
					BSf3 = 0.082800;
					BSg3 = 0.246000;
					BSh3 = 0.728000;
					BSi3 = 0.000000;
					BSa4 = 1.574000;
					BSb4 = 1.302000;
					BSc4 = 0.108000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 6)	{
					BSa1 = -0.037600;
					BSb1 = 0.765000;
					BSc1 = 0.001760;
					BSd1 = 0.888000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.120000;
					BSq2 = 0.313000;
					BSr2 = 0.350000;
					BSu2 = 0.187000;
					BSv2 = 0.620000;
					BSw2 = 0.675000;
					BSx2 = 0.500000;
					BSs2 = -0.036400;
					BSt2 = 0.023200;
					BSy2 = 0.184000;
					BSa2 = 0.821000;
					BSb2 = 0.731000;
					BSc2 = 0.354000;
					BSd2 = 0.263000;
					BSe2 = 2.888000;
					BSf2 = -0.060800;
					BSg2 = -0.953000;
					BSh2 = 0.000000;
					BSa3 = 0.137000;
					BSb3 = 0.330000;
					BSc3 = -0.182000;
					BSe3 = 0.039900;
					BSf3 = 0.083300;
					BSg3 = 0.255000;
					BSh3 = 0.727000;
					BSi3 = 0.000000;
					BSa4 = 1.508000;
					BSb4 = 1.371000;
					BSc4 = 0.103000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 7)	{
					BSa1 = 0.202000;
					BSb1 = 2.921000;
					BSc1 = 0.004110;
					BSd1 = -0.002890;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.439000;
					BSr2 = 0.268000;
					BSu2 = 0.353000;
					BSv2 = 0.620000;
					BSw2 = 0.781000;
					BSx2 = 0.500000;
					BSs2 = -0.070300;
					BSt2 = 0.023100;
					BSy2 = 0.170000;
					BSa2 = 0.643000;
					BSb2 = 0.878000;
					BSc2 = 0.391000;
					BSd2 = 0.266000;
					BSe2 = 3.079000;
					BSf2 = -0.081900;
					BSg2 = -0.996000;
					BSh2 = 0.000000;
					BSa3 = 0.728000;
					BSb3 = 0.318000;
					BSc3 = -0.349000;
					BSe3 = -0.078100;
					BSf3 = 0.077200;
					BSg3 = 0.299000;
					BSh3 = 0.654000;
					BSi3 = 0.000000;
					BSa4 = 2.866000;
					BSb4 = 0.116000;
					BSc4 = 0.411000;
					ICol = 4;
				}
				else if (BSFoundationCfg == 8)	{
					BSa1 = 0.128000;
					BSb1 = 2.951000;
					BSc1 = 0.003960;
					BSd1 = -0.002840;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.443000;
					BSr2 = 0.268000;
					BSu2 = 0.355000;
					BSv2 = 0.620000;
					BSw2 = 0.777000;
					BSx2 = 0.500000;
					BSs2 = -0.069800;
					BSt2 = 0.023100;
					BSy2 = 0.169000;
					BSa2 = 0.631000;
					BSb2 = 0.890000;
					BSc2 = 0.386000;
					BSd2 = 0.258000;
					BSe2 = 3.021000;
					BSf2 = -0.071300;
					BSg2 = -0.980000;
					BSh2 = 0.000000;
					BSa3 = 0.726000;
					BSb3 = 0.317000;
					BSc3 = -0.342000;
					BSe3 = -0.079500;
					BSf3 = 0.077500;
					BSg3 = 0.301000;
					BSh3 = 0.653000;
					BSi3 = 0.000000;
					BSa4 = 2.874000;
					BSb4 = 0.109000;
					BSc4 = 0.426000;
					ICol = 4;
				}
				else if (BSFoundationCfg == 9) {
					BSa1 = 0.001340;
					BSb1 = 2.936000;
					BSc1 = 0.095700;
					BSd1 = 0.000000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 1.302000;
					BSr2 = 0.279000;
					BSu2 = 0.344000;
					BSv2 = 0.620000;
					BSw2 = 0.690000;
					BSx2 = 0.500000;
					BSs2 = -0.034700;
					BSt2 = 0.020300;
					BSy2 = 0.086000;
					BSa2 = 0.000000;
					BSb2 = 0.000000;
					BSc2 = 0.000000;
					BSd2 = 0.000000;
					BSe2 = 0.000000;
					BSf2 = 0.000000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.459000;
					BSb3 = 0.829000;
					BSc3 = 0.174000;
					BSe3 = 0.000000;
					BSf3 = 0.000000;
					BSg3 = 0.000000;
					BSh3 = 0.000000;
					BSi3 = 0.000000;
					BSa4 = 3.064000;
					BSb4 = 0.000000;
					BSc4 = 0.000000;
					ICol = 1;
				}
				else if (BSFoundationCfg == 10)	{
					BSa1 = -0.050500;
					BSb1 = 2.959000;
					BSc1 = 0.083000;
					BSd1 = 0.000000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 1.309000;
					BSr2 = 0.279000;
					BSu2 = 0.348000;
					BSv2 = 0.620000;
					BSw2 = 0.686000;
					BSx2 = 0.500000;
					BSs2 = -0.034200;
					BSt2 = 0.020200;
					BSy2 = 0.085000;
					BSa2 = 0.000000;
					BSb2 = 0.000000;
					BSc2 = 0.000000;
					BSd2 = 0.000000;
					BSe2 = 0.000000;
					BSf2 = 0.000000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.453000;
					BSb3 = 0.830000;
					BSc3 = 0.186000;
					BSe3 = 0.000000;
					BSf3 = 0.000000;
					BSg3 = 0.000000;
					BSh3 = 0.000000;
					BSi3 = 0.000000;
					BSa4 = 3.065000;
					BSb4 = 0.000000;
					BSc4 = 0.000000;
					ICol = 1;
				}
				else if (BSFoundationCfg == 11)	{
					BSa1 = 1.038000;
					BSb1 = 1.412000;
					BSc1 = 0.038300;
					BSd1 = 0.279000;
					BSe1 = 2.690000;
					BSf1 = 0.468000;
					BSg1 = 1.210000;
					BSh1 = -0.940000;
					BSi1 = 1.000000;
					BSj1 = 0.163000;
					BSq2 = 0.123000;
					BSr2 = 0.551000;
					BSu2 = 0.204000;
					BSv2 = 0.620000;
					BSw2 = 0.707000;
					BSx2 = 0.500000;
					BSs2 = 0.654000;
					BSt2 = 0.013500;
					BSy2 = 0.025400;
					BSa2 = 0.410000;
					BSb2 = -0.077700;
					BSc2 = 0.590000;
					BSd2 = -0.318000;
					BSe2 = -0.011800;
					BSf2 = -0.035800;
					BSg2 = 0.096100;
					BSh2 = 0.255000;
					BSa3 = 0.119000;
					BSb3 = 0.295000;
					BSc3 = -0.063000;
					BSe3 = 0.460000;
					BSf3 = 0.086400;
					BSg3 = -0.058300;
					BSh3 = 0.237000;
					BSi3 = 1.889000;
					BSa4 = 2.532000;
					BSb4 = 0.352000;
					BSc4 = 0.338000;
					ICol = 99;
				}
				else if (BSFoundationCfg == 12)	{
					BSa1 = 0.900000;
					BSb1 = 1.562000;
					BSc1 = 0.038000;
					BSd1 = 0.282000;
					BSe1 = 2.667000;
					BSf1 = 0.409000;
					BSg1 = 1.114000;
					BSh1 = -0.873000;
					BSi1 = 1.000000;
					BSj1 = 0.141000;
					BSq2 = 0.121000;
					BSr2 = 0.553000;
					BSu2 = 0.206000;
					BSv2 = 0.620000;
					BSw2 = 0.712000;
					BSx2 = 0.500000;
					BSs2 = 0.667000;
					BSt2 = 0.013200;
					BSy2 = 0.023400;
					BSa2 = 0.418000;
					BSb2 = -0.069000;
					BSc2 = 0.590000;
					BSd2 = -0.316000;
					BSe2 = -0.009290;
					BSf2 = -0.035800;
					BSg2 = 0.092800;
					BSh2 = 0.252000;
					BSa3 = 0.116000;
					BSb3 = 0.295000;
					BSc3 = -0.059600;
					BSe3 = 0.458000;
					BSf3 = 0.086400;
					BSg3 = -0.052900;
					BSh3 = 0.234000;
					BSi3 = 1.848000;
					BSa4 = 2.541000;
					BSb4 = 0.346000;
					BSc4 = 0.339000;
					ICol = 99;
				}
				else if (BSFoundationCfg == 13)	{
					BSa1 = -0.018000;
					BSb1 = 1.570000;
					BSc1 = 0.003180;
					BSd1 = 0.000000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.612000;
					BSr2 = 0.125000;
					BSu2 = 0.741000;
					BSv2 = 0.620000;
					BSw2 = 0.587000;
					BSx2 = 0.500000;
					BSs2 = -0.047600;
					BSt2 = 0.016700;
					BSy2 = 0.064500;
					BSa2 = 0.000000;
					BSb2 = 0.000000;
					BSc2 = 0.000000;
					BSd2 = 0.000000;
					BSe2 = 0.000000;
					BSf2 = 0.000000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.159000;
					BSb3 = 0.392000;
					BSc3 = 0.284000;
					BSe3 = 0.000000;
					BSf3 = 0.000000;
					BSg3 = 0.000000;
					BSh3 = 0.000000;
					BSi3 = 0.000000;
					BSa4 = 2.983000;
					BSb4 = 0.000000;
					BSc4 = 0.000000;
					ICol = 1;
				}
				else if (BSFoundationCfg == 14)	{
					BSa1 = -0.008910;
					BSb1 = 0.647000;
					BSc1 = 0.003130;
					BSd1 = 0.798000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.003850;
					BSr2 = 0.102000;
					BSu2 = 0.570000;
					BSv2 = 0.620000;
					BSw2 = 1.068000;
					BSx2 = 0.500000;
					BSs2 = -0.042400;
					BSt2 = 0.018800;
					BSy2 = 0.134000;
					BSa2 = 0.494000;
					BSb2 = 0.797000;
					BSc2 = 0.347000;
					BSd2 = 0.064100;
					BSe2 = 1.105000;
					BSf2 = 0.094900;
					BSg2 = -0.246000;
					BSh2 = 0.000000;
					BSa3 = 0.104000;
					BSb3 = 0.112000;
					BSc3 = -0.087000;
					BSe3 = 0.007320;
					BSf3 = 0.086600;
					BSg3 = 0.244000;
					BSh3 = 0.659000;
					BSi3 = 0.000000;
					BSa4 = -8.379000;
					BSb4 = 11.286000;
					BSc4 = 0.010300;
					ICol = 3;
				}
				else if (BSFoundationCfg == 15)	{
					BSa1 = -0.008370;
					BSb1 = 0.647000;
					BSc1 = 0.003180;
					BSd1 = 0.797000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.094500;
					BSr2 = 0.093700;
					BSu2 = 0.656000;
					BSv2 = 0.620000;
					BSw2 = 1.115000;
					BSx2 = 0.500000;
					BSs2 = -0.063400;
					BSt2 = 0.018600;
					BSy2 = 0.137000;
					BSa2 = 0.366000;
					BSb2 = 0.947000;
					BSc2 = 0.404000;
					BSd2 = 0.049700;
					BSe2 = 1.356000;
					BSf2 = 0.074900;
					BSg2 = -0.330000;
					BSh2 = 0.000000;
					BSa3 = 0.205000;
					BSb3 = 0.141000;
					BSc3 = -0.129000;
					BSe3 = -0.028100;
					BSf3 = 0.084500;
					BSg3 = 0.245000;
					BSh3 = 0.625000;
					BSi3 = 0.000000;
					BSa4 = -2.996000;
					BSb4 = 5.837000;
					BSc4 = 0.023600;
					ICol = 3;
				}
				else if (BSFoundationCfg == 16)	{
					BSa1 = -0.004050;
					BSb1 = 0.635000;
					BSc1 = 0.002900;
					BSd1 = 0.782000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.126000;
					BSr2 = 0.295000;
					BSu2 = 0.783000;
					BSv2 = 0.620000;
					BSw2 = 0.128000;
					BSx2 = 0.500000;
					BSs2 = 0.377000;
					BSt2 = 0.006610;
					BSy2 = 0.163000;
					BSa2 = 0.508000;
					BSb2 = 0.507000;
					BSc2 = 0.619000;
					BSd2 = -0.161000;
					BSe2 = 0.166000;
					BSf2 = -0.008010;
					BSg2 = -0.041300;
					BSh2 = 0.000000;
					BSa3 = -0.276000;
					BSb3 = 0.306000;
					BSc3 = 0.279000;
					BSe3 = 0.214000;
					BSf3 = 0.050600;
					BSg3 = -0.018200;
					BSh3 = 0.656000;
					BSi3 = 0.000000;
					BSa4 = 2.697000;
					BSb4 = 0.196000;
					BSc4 = 0.402000;
					ICol = 2;
				}
				else if (BSFoundationCfg == 17)	{
					BSa1 = -0.030900;
					BSb1 = 0.651000;
					BSc1 = -0.000065;
					BSd1 = 0.809000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075600;
					BSq2 = 0.015100;
					BSr2 = 0.103000;
					BSu2 = 0.561000;
					BSv2 = 0.620000;
					BSw2 = 1.051000;
					BSx2 = 0.500000;
					BSs2 = -0.039000;
					BSt2 = 0.018700;
					BSy2 = 0.136000;
					BSa2 = 0.506000;
					BSb2 = 0.791000;
					BSc2 = 0.340000;
					BSd2 = 0.067000;
					BSe2 = 1.134000;
					BSf2 = 0.084400;
					BSg2 = -0.258000;
					BSh2 = 0.000000;
					BSa3 = 0.092500;
					BSb3 = 0.111000;
					BSc3 = -0.088900;
					BSe3 = 0.014500;
					BSf3 = 0.087000;
					BSg3 = 0.250000;
					BSh3 = 0.633000;
					BSi3 = 0.000000;
					BSa4 = 2.343000;
					BSb4 = 0.645000;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 18)	{
					BSa1 = -0.011000;
					BSb1 = 1.569000;
					BSc1 = -0.000076;
					BSd1 = -0.000364;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.046100;
					BSr2 = 0.091400;
					BSu2 = 0.651000;
					BSv2 = 0.620000;
					BSw2 = 1.135000;
					BSx2 = 0.500000;
					BSs2 = -0.052100;
					BSt2 = 0.018900;
					BSy2 = 0.130000;
					BSa2 = 0.470000;
					BSb2 = 0.838000;
					BSc2 = 0.345000;
					BSd2 = 0.063300;
					BSe2 = 1.093000;
					BSf2 = 0.098500;
					BSg2 = -0.245000;
					BSh2 = 0.000000;
					BSa3 = 0.188000;
					BSb3 = 0.107000;
					BSc3 = -0.121000;
					BSe3 = 0.002400;
					BSf3 = 0.086600;
					BSg3 = 0.268000;
					BSh3 = 0.557000;
					BSi3 = 0.000000;
					BSa4 = 2.607000;
					BSb4 = 0.342000;
					BSc4 = 1.000000;
					ICol = 4;
				}
				else if (BSFoundationCfg == 19)	{
					BSa1 = -0.003740;
					BSb1 = 0.724000;
					BSc1 = 0.116000;
					BSd1 = 0.757000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.996000;
					BSr2 = -0.057600;
					BSu2 = 0.620000;
					BSv2 = 0.620000;
					BSw2 = -2.850000;
					BSx2 = 0.735000;
					BSs2 = 0.746000;
					BSt2 = 0.007420;
					BSy2 = -0.201000;
					BSa2 = 3.429000;
					BSb2 = 0.176000;
					BSc2 = 0.522000;
					BSd2 = 0.269000;
					BSe2 = 0.155000;
					BSf2 = 0.047800;
					BSg2 = 0.039400;
					BSh2 = 0.000000;
					BSa3 = 0.235000;
					BSb3 = 0.050600;
					BSc3 = -0.082100;
					BSe3 = 0.025500;
					BSf3 = 0.244000;
					BSg3 = 0.176000;
					BSh3 = 0.711000;
					BSi3 = 0.000000;
					BSa4 = 2.813000;
					BSb4 = 0.118000;
					BSc4 = 0.729000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 20)	{
					BSa1 = 0.016000;
					BSb1 = 0.709000;
					BSc1 = 0.103000;
					BSd1 = 0.715000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 11.396000;
					BSr2 = -0.640000;
					BSu2 = 0.095400;
					BSv2 = 0.620000;
					BSw2 = -1.869000;
					BSx2 = 0.500000;
					BSs2 = 1.213000;
					BSt2 = -0.004700;
					BSy2 = -0.278000;
					BSa2 = 8.129000;
					BSb2 = 0.157000;
					BSc2 = 0.775000;
					BSd2 = 0.537000;
					BSe2 = 0.047700;
					BSf2 = -0.015100;
					BSg2 = 0.019700;
					BSh2 = 0.000000;
					BSa3 = 0.237000;
					BSb3 = 0.288000;
					BSc3 = -0.171000;
					BSe3 = 0.013900;
					BSf3 = 0.113000;
					BSg3 = 0.219000;
					BSh3 = 0.770000;
					BSi3 = 0.000000;
					BSa4 = 2.344000;
					BSb4 = 0.558000;
					BSc4 = 0.233000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 21)	{
					BSa1 = 0.012800;
					BSb1 = 0.711000;
					BSc1 = 0.105000;
					BSd1 = 0.722000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 10.942000;
					BSr2 = -0.637000;
					BSu2 = 0.086300;
					BSv2 = 0.620000;
					BSw2 = -1.867000;
					BSx2 = 0.500000;
					BSs2 = 1.206000;
					BSt2 = -0.000041;
					BSy2 = -0.284000;
					BSa2 = 7.541000;
					BSb2 = 0.168000;
					BSc2 = 0.761000;
					BSd2 = 0.530000;
					BSe2 = 0.058200;
					BSf2 = -0.016300;
					BSg2 = 0.020200;
					BSh2 = 0.000000;
					BSa3 = 0.220000;
					BSb3 = 0.251000;
					BSc3 = -0.153000;
					BSe3 = 0.027100;
					BSf3 = 0.123000;
					BSg3 = 0.214000;
					BSh3 = 0.768000;
					BSi3 = 0.000000;
					BSa4 = 2.246000;
					BSb4 = 0.658000;
					BSc4 = 0.199000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 22)	{
					BSa1 = -0.002110;
					BSb1 = 0.748000;
					BSc1 = 0.093800;
					BSd1 = 0.845000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -0.140000;
					BSr2 = 0.009540;
					BSu2 = 0.662000;
					BSv2 = 0.620000;
					BSw2 = 8.476000;
					BSx2 = 0.500000;
					BSs2 = 0.453000;
					BSt2 = 0.007740;
					BSy2 = -0.007940;
					BSa2 = 3.002000;
					BSb2 = 0.248000;
					BSc2 = 0.506000;
					BSd2 = 0.256000;
					BSe2 = 0.363000;
					BSf2 = 0.058000;
					BSg2 = 0.013500;
					BSh2 = 0.000000;
					BSa3 = -0.013600;
					BSb3 = -0.026300;
					BSc3 = -0.003810;
					BSe3 = 0.166000;
					BSf3 = 0.267000;
					BSg3 = 0.158000;
					BSh3 = 0.679000;
					BSi3 = 0.000000;
					BSa4 = 2.526000;
					BSb4 = 0.388000;
					BSc4 = 0.255000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 23)	{
					BSa1 = 0.020900;
					BSb1 = 0.723000;
					BSc1 = 0.081600;
					BSd1 = 0.776000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -0.068300;
					BSr2 = 0.375000;
					BSu2 = 0.136000;
					BSv2 = 0.620000;
					BSw2 = 0.630000;
					BSx2 = 0.500000;
					BSs2 = -0.003790;
					BSt2 = 0.021900;
					BSy2 = 0.177000;
					BSa2 = 0.785000;
					BSb2 = 0.567000;
					BSc2 = 0.342000;
					BSd2 = 0.142000;
					BSe2 = 2.234000;
					BSf2 = 0.003520;
					BSg2 = -0.686000;
					BSh2 = 0.000000;
					BSa3 = 0.037400;
					BSb3 = 0.254000;
					BSc3 = -0.117000;
					BSe3 = 0.103000;
					BSf3 = 0.107000;
					BSg3 = 0.222000;
					BSh3 = 0.778000;
					BSi3 = 0.000000;
					BSa4 = 16.468000;
					BSb4 = -13.584000;
					BSc4 = -0.012400;
					ICol = 3;
				}
				else if (BSFoundationCfg == 24) {
					BSa1 = 0.017100;
					BSb1 = 0.727000;
					BSc1 = 0.082500;
					BSd1 = 0.789000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -0.325000;
					BSr2 = 0.374000;
					BSu2 = 0.131000;
					BSv2 = 0.620000;
					BSw2 = 0.637000;
					BSx2 = 0.500000;
					BSs2 = -0.001340;
					BSt2 = 0.022800;
					BSy2 = 0.170000;
					BSa2 = 0.830000;
					BSb2 = 0.513000;
					BSc2 = 0.346000;
					BSd2 = 0.121000;
					BSe2 = 1.849000;
					BSf2 = 0.043600;
					BSg2 = -0.530000;
					BSh2 = 0.000000;
					BSa3 = 0.010700;
					BSb3 = 0.210000;
					BSc3 = -0.095400;
					BSe3 = 0.124000;
					BSf3 = 0.117000;
					BSg3 = 0.218000;
					BSh3 = 0.778000;
					BSi3 = 0.000000;
					BSa4 = 4.396000;
					BSb4 = -1.511000;
					BSc4 = -0.113000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 25)	{
					BSa1 = -0.084100;
					BSb1 = 0.748000;
					BSc1 = 0.004410;
					BSd1 = 0.884000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.356000;
					BSq2 = -0.006140;
					BSr2 = 0.015300;
					BSu2 = 0.592000;
					BSv2 = 0.620000;
					BSw2 = 3.497000;
					BSx2 = 0.500000;
					BSs2 = -0.228000;
					BSt2 = 0.019900;
					BSy2 = 0.088500;
					BSa2 = 2.157000;
					BSb2 = 0.358000;
					BSc2 = 0.556000;
					BSd2 = 0.312000;
					BSe2 = 0.533000;
					BSf2 = 0.150000;
					BSg2 = -0.040800;
					BSh2 = 0.000000;
					BSa3 = 0.140000;
					BSb3 = 0.097800;
					BSc3 = -0.095100;
					BSe3 = 0.047200;
					BSf3 = 0.198000;
					BSg3 = 0.213000;
					BSh3 = 0.682000;
					BSi3 = 0.000000;
					BSa4 = 2.188000;
					BSb4 = 0.711000;
					BSc4 = 0.132000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 26)	{
					BSa1 = -0.083600;
					BSb1 = 0.749000;
					BSc1 = 0.003390;
					BSd1 = 0.885000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.357000;
					BSq2 = 0.188000;
					BSr2 = 0.354000;
					BSu2 = 0.197000;
					BSv2 = 0.620000;
					BSw2 = 0.666000;
					BSx2 = 0.500000;
					BSs2 = -0.025900;
					BSt2 = 0.022300;
					BSy2 = 0.183000;
					BSa2 = 0.786000;
					BSb2 = 0.683000;
					BSc2 = 0.353000;
					BSd2 = 0.229000;
					BSe2 = 2.850000;
					BSf2 = -0.047800;
					BSg2 = -0.939000;
					BSh2 = 0.000000;
					BSa3 = 0.129000;
					BSb3 = 0.304000;
					BSc3 = -0.169000;
					BSe3 = 0.044400;
					BSf3 = 0.085800;
					BSg3 = 0.247000;
					BSh3 = 0.730000;
					BSi3 = 0.000000;
					BSa4 = 1.232000;
					BSb4 = 1.644000;
					BSc4 = 0.088100;
					ICol = 5;
				}
				else if (BSFoundationCfg == 27)	{
					BSa1 = -0.083700;
					BSb1 = 0.749000;
					BSc1 = 0.003520;
					BSd1 = 0.885000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.356000;
					BSq2 = 0.039000;
					BSr2 = 0.345000;
					BSu2 = 0.214000;
					BSv2 = 0.620000;
					BSw2 = 0.668000;
					BSx2 = 0.500000;
					BSs2 = -0.033200;
					BSt2 = 0.022900;
					BSy2 = 0.179000;
					BSa2 = 0.794000;
					BSb2 = 0.617000;
					BSc2 = 0.359000;
					BSd2 = 0.199000;
					BSe2 = 2.673000;
					BSf2 = -0.024200;
					BSg2 = -0.865000;
					BSh2 = 0.000000;
					BSa3 = 0.121000;
					BSb3 = 0.273000;
					BSc3 = -0.155000;
					BSe3 = 0.051300;
					BSf3 = 0.093300;
					BSg3 = 0.245000;
					BSh3 = 0.731000;
					BSi3 = 0.000000;
					BSa4 = 0.408000;
					BSb4 = 2.469000;
					BSc4 = 0.059100;
					ICol = 5;
				}
				else if (BSFoundationCfg == 28)	{
					BSa1 = 0.024900;
					BSb1 = 0.000000;
					BSc1 = 0.004650;
					BSd1 = 0.000000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.398000;
					BSr2 = 0.423000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.592000;
					BSt2 = 0.010300;
					BSy2 = 0.000000;
					BSa2 = 2.193000;
					BSb2 = 0.000000;
					BSc2 = 0.774000;
					BSd2 = 0.176000;
					BSe2 = 0.000000;
					BSf2 = 0.000000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.144000;
					BSb3 = 0.502000;
					BSc3 = 0.000000;
					BSe3 = 0.000000;
					BSf3 = 0.000000;
					BSg3 = 0.000000;
					BSh3 = 0.000000;
					BSi3 = 0.000000;
					BSa4 = 2.878000;
					BSb4 = 0.000000;
					BSc4 = 0.000000;
					ICol = 98;
				}
				else if (BSFoundationCfg == 29) {
					BSa1 = 0.059300;
					BSb1 = 0.000000;
					BSc1 = 0.008790;
					BSd1 = 0.000000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.441000;
					BSr2 = 0.423000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.588000;
					BSt2 = 0.010400;
					BSy2 = 0.000000;
					BSa2 = 2.241000;
					BSb2 = 0.000000;
					BSc2 = 0.747000;
					BSd2 = 0.173000;
					BSe2 = 0.000000;
					BSf2 = 0.000000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.193000;
					BSb3 = 0.496000;
					BSc3 = 0.000000;
					BSe3 = 0.000000;
					BSf3 = 0.000000;
					BSg3 = 0.000000;
					BSh3 = 0.000000;
					BSi3 = 0.000000;
					BSa4 = 2.897000;
					BSb4 = 0.000000;
					BSc4 = 0.000000;
					ICol = 98;
				}
				else if (BSFoundationCfg == 30)	{
					BSa1 = 0.027200;
					BSb1 = 0.000000;
					BSc1 = 0.005820;
					BSd1 = -0.026400;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.424000;
					BSr2 = 0.432000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.597000;
					BSt2 = 0.010800;
					BSy2 = 0.000000;
					BSa2 = 2.102000;
					BSb2 = 0.000000;
					BSc2 = 0.811000;
					BSd2 = 0.180000;
					BSe2 = 0.003440;
					BSf2 = -0.000649;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.018200;
					BSb3 = 0.517000;
					BSc3 = 0.000000;
					BSe3 = 0.049300;
					BSf3 = -0.011400;
					BSg3 = 0.000000;
					BSh3 = 0.495000;
					BSi3 = 0.000000;
					BSa4 = 2.821000;
					BSb4 = 0.020400;
					BSc4 = 0.512000;
					ICol = 98;
				}
				else if (BSFoundationCfg == 31)	{
					BSa1 = 0.070100;
					BSb1 = 0.000000;
					BSc1 = 0.014500;
					BSd1 = -0.032800;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.366000;
					BSr2 = 0.431000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.602000;
					BSt2 = 0.010500;
					BSy2 = 0.000000;
					BSa2 = 2.116000;
					BSb2 = 0.000000;
					BSc2 = 0.794000;
					BSd2 = 0.177000;
					BSe2 = 0.005250;
					BSf2 = -0.001490;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.021700;
					BSb3 = 0.510000;
					BSc3 = 0.000000;
					BSe3 = 0.066500;
					BSf3 = -0.010900;
					BSg3 = 0.000000;
					BSh3 = 0.466000;
					BSi3 = 0.000000;
					BSa4 = 2.826000;
					BSb4 = 0.026400;
					BSc4 = 0.455000;
					ICol = 98;
				}
				else if (BSFoundationCfg == 32)	{
					BSa1 = 0.040300;
					BSb1 = 0.000000;
					BSc1 = 0.001420;
					BSd1 = -0.067100;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.354000;
					BSr2 = 0.431000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.613000;
					BSt2 = 0.009350;
					BSy2 = 0.000000;
					BSa2 = 1.955000;
					BSb2 = 0.000000;
					BSc2 = 0.822000;
					BSd2 = 0.219000;
					BSe2 = 0.058300;
					BSf2 = -0.018500;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.101000;
					BSb3 = 0.357000;
					BSc3 = 0.000000;
					BSe3 = 0.159000;
					BSf3 = 0.003710;
					BSg3 = 0.000000;
					BSh3 = 0.669000;
					BSi3 = 0.000000;
					BSa4 = 2.428000;
					BSb4 = 0.296000;
					BSc4 = 0.362000;
					ICol = 98;
				}
				else if (BSFoundationCfg == 33)	{
					BSa1 = 0.106000;
					BSb1 = 0.000000;
					BSc1 = -0.000446;
					BSd1 = -0.071100;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.350000;
					BSr2 = 0.434000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.623000;
					BSt2 = 0.009410;
					BSy2 = 0.000000;
					BSa2 = 2.016000;
					BSb2 = 0.000000;
					BSc2 = 0.785000;
					BSd2 = 0.212000;
					BSe2 = 0.045800;
					BSf2 = -0.013600;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.072800;
					BSb3 = 0.352000;
					BSc3 = 0.000000;
					BSe3 = 0.160000;
					BSf3 = 0.007080;
					BSg3 = 0.000000;
					BSh3 = 0.627000;
					BSi3 = 0.000000;
					BSa4 = 2.511000;
					BSb4 = 0.244000;
					BSc4 = 0.356000;
					ICol = 98;
				}
				else if (BSFoundationCfg == 34)	{
					BSa1 = 0.024300;
					BSb1 = 0.000000;
					BSc1 = 0.004060;
					BSd1 = 0.002510;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.313000;
					BSr2 = 0.412000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.562000;
					BSt2 = 0.011300;
					BSy2 = 0.000000;
					BSa2 = 2.075000;
					BSb2 = 0.000000;
					BSc2 = 0.763000;
					BSd2 = 0.180000;
					BSe2 = 0.010100;
					BSf2 = -0.001210;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.108000;
					BSb3 = 0.432000;
					BSc3 = 0.000000;
					BSe3 = 0.024800;
					BSf3 = 0.009120;
					BSg3 = 0.000000;
					BSh3 = 0.685000;
					BSi3 = 0.000000;
					BSa4 = 2.826000;
					BSb4 = 0.016200;
					BSc4 = 0.636000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 35)	{
					BSa1 = 0.056500;
					BSb1 = 0.000000;
					BSc1 = 0.007810;
					BSd1 = 0.004430;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.345000;
					BSr2 = 0.411000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.558000;
					BSt2 = 0.011300;
					BSy2 = 0.000000;
					BSa2 = 2.120000;
					BSb2 = 0.000000;
					BSc2 = 0.737000;
					BSd2 = 0.178000;
					BSe2 = 0.011500;
					BSf2 = -0.002150;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.148000;
					BSb3 = 0.427000;
					BSc3 = 0.000000;
					BSe3 = 0.027900;
					BSf3 = 0.008520;
					BSg3 = 0.000000;
					BSh3 = 0.671000;
					BSi3 = 0.000000;
					BSa4 = 2.848000;
					BSb4 = 0.015100;
					BSc4 = 0.615000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 36)	{
					BSa1 = 0.025400;
					BSb1 = 0.000000;
					BSc1 = 0.004240;
					BSd1 = -0.014400;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.269000;
					BSr2 = 0.404000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.546000;
					BSt2 = 0.011500;
					BSy2 = 0.000000;
					BSa2 = 2.051000;
					BSb2 = 0.000000;
					BSc2 = 0.771000;
					BSd2 = 0.183000;
					BSe2 = 0.014800;
					BSf2 = -0.002720;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.082000;
					BSb3 = 0.429000;
					BSc3 = 0.000000;
					BSe3 = 0.032200;
					BSf3 = 0.007170;
					BSg3 = 0.000000;
					BSh3 = 0.730000;
					BSi3 = 0.000000;
					BSa4 = 2.818000;
					BSb4 = 0.014000;
					BSc4 = 0.768000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 37)	{
					BSa1 = 0.057000;
					BSb1 = 0.000000;
					BSc1 = 0.007130;
					BSd1 = -0.011000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.286000;
					BSr2 = 0.402000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.545000;
					BSt2 = 0.011300;
					BSy2 = 0.000000;
					BSa2 = 2.094000;
					BSb2 = 0.000000;
					BSc2 = 0.743000;
					BSd2 = 0.180000;
					BSe2 = 0.014000;
					BSf2 = -0.002790;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.125000;
					BSb3 = 0.423000;
					BSc3 = 0.000000;
					BSe3 = 0.030500;
					BSf3 = 0.007520;
					BSg3 = 0.000000;
					BSh3 = 0.741000;
					BSi3 = 0.000000;
					BSa4 = 2.844000;
					BSb4 = 0.010100;
					BSc4 = 0.823000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 38)	{
					BSa1 = 0.026600;
					BSb1 = 0.000000;
					BSc1 = 0.005120;
					BSd1 = -0.023900;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.329000;
					BSr2 = 0.419000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.566000;
					BSt2 = 0.011600;
					BSy2 = 0.000000;
					BSa2 = 1.998000;
					BSb2 = 0.000000;
					BSc2 = 0.796000;
					BSd2 = 0.185000;
					BSe2 = 0.015200;
					BSf2 = -0.003010;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.004810;
					BSb3 = 0.443000;
					BSc3 = 0.000000;
					BSe3 = 0.060500;
					BSf3 = 0.001520;
					BSg3 = 0.000000;
					BSh3 = 0.619000;
					BSi3 = 0.000000;
					BSa4 = 2.756000;
					BSb4 = 0.044700;
					BSc4 = 0.556000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 39)	{
					BSa1 = 0.067100;
					BSb1 = 0.000000;
					BSc1 = 0.012900;
					BSd1 = -0.028900;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.270000;
					BSr2 = 0.418000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.570000;
					BSt2 = 0.011300;
					BSy2 = 0.000000;
					BSa2 = 2.007000;
					BSb2 = 0.000000;
					BSc2 = 0.781000;
					BSd2 = 0.183000;
					BSe2 = 0.017400;
					BSf2 = -0.003980;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.006410;
					BSb3 = 0.436000;
					BSc3 = 0.000000;
					BSe3 = 0.076000;
					BSf3 = 0.001740;
					BSg3 = 0.000000;
					BSh3 = 0.584000;
					BSi3 = 0.000000;
					BSa4 = 2.761000;
					BSb4 = 0.051600;
					BSc4 = 0.507000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 40)	{
					BSa1 = 0.038600;
					BSb1 = 0.000000;
					BSc1 = 0.001500;
					BSd1 = -0.064300;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.298000;
					BSr2 = 0.422000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.587000;
					BSt2 = 0.010300;
					BSy2 = 0.000000;
					BSa2 = 1.881000;
					BSb2 = 0.000000;
					BSc2 = 0.813000;
					BSd2 = 0.222000;
					BSe2 = 0.061500;
					BSf2 = -0.018800;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.097700;
					BSb3 = 0.323000;
					BSc3 = 0.000000;
					BSe3 = 0.145000;
					BSf3 = 0.007050;
					BSg3 = 0.000000;
					BSh3 = 0.701000;
					BSi3 = 0.000000;
					BSa4 = 2.346000;
					BSb4 = 0.335000;
					BSc4 = 0.375000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 41)	{
					BSa1 = 0.099500;
					BSb1 = 0.000000;
					BSc1 = 0.000022;
					BSd1 = -0.066600;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.283000;
					BSr2 = 0.423000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.594000;
					BSt2 = 0.010400;
					BSy2 = 0.000000;
					BSa2 = 1.931000;
					BSb2 = 0.000000;
					BSc2 = 0.779000;
					BSd2 = 0.215000;
					BSe2 = 0.050500;
					BSf2 = -0.014500;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.074900;
					BSb3 = 0.319000;
					BSc3 = 0.000000;
					BSe3 = 0.148000;
					BSf3 = 0.009450;
					BSg3 = 0.000000;
					BSh3 = 0.660000;
					BSi3 = 0.000000;
					BSa4 = 2.438000;
					BSb4 = 0.277000;
					BSc4 = 0.370000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 42)	{
					BSa1 = 0.026600;
					BSb1 = 0.000000;
					BSc1 = 0.004020;
					BSd1 = -0.003710;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.396000;
					BSr2 = 0.420000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.588000;
					BSt2 = 0.010100;
					BSy2 = 0.000000;
					BSa2 = 2.088000;
					BSb2 = 0.000000;
					BSc2 = 0.765000;
					BSd2 = 0.206000;
					BSe2 = 0.018500;
					BSf2 = -0.004170;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.072500;
					BSb3 = 0.369000;
					BSc3 = 0.000000;
					BSe3 = 0.053100;
					BSf3 = 0.007400;
					BSg3 = 0.000000;
					BSh3 = 0.755000;
					BSi3 = 0.000000;
					BSa4 = 2.766000;
					BSb4 = 0.031800;
					BSc4 = 0.724000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 43)	{
					BSa1 = 0.061100;
					BSb1 = 0.000000;
					BSc1 = 0.008110;
					BSd1 = -0.001550;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.419000;
					BSr2 = 0.416000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.583000;
					BSt2 = 0.009900;
					BSy2 = 0.000000;
					BSa2 = 2.148000;
					BSb2 = 0.000000;
					BSc2 = 0.732000;
					BSd2 = 0.204000;
					BSe2 = 0.019300;
					BSf2 = -0.004850;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.118000;
					BSb3 = 0.362000;
					BSc3 = 0.000000;
					BSe3 = 0.054800;
					BSf3 = 0.007470;
					BSg3 = 0.000000;
					BSh3 = 0.745000;
					BSi3 = 0.000000;
					BSa4 = 2.797000;
					BSb4 = 0.027900;
					BSc4 = 0.715000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 44)	{
					BSa1 = 0.023100;
					BSb1 = 0.000000;
					BSc1 = 0.004580;
					BSd1 = 0.009590;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.203000;
					BSr2 = 0.397000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.521000;
					BSt2 = 0.013500;
					BSy2 = 0.000000;
					BSa2 = 1.998000;
					BSb2 = 0.000000;
					BSc2 = 0.733000;
					BSd2 = 0.180000;
					BSe2 = 0.027100;
					BSf2 = -0.003600;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.100000;
					BSb3 = 0.364000;
					BSc3 = 0.000000;
					BSe3 = 0.048600;
					BSf3 = 0.026500;
					BSg3 = 0.000000;
					BSh3 = 0.677000;
					BSi3 = 0.000000;
					BSa4 = 2.791000;
					BSb4 = 0.034700;
					BSc4 = 0.631000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 45)	{
					BSa1 = 0.053200;
					BSb1 = 0.000000;
					BSc1 = 0.009350;
					BSd1 = 0.012500;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.234000;
					BSr2 = 0.396000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.518000;
					BSt2 = 0.013400;
					BSy2 = 0.000000;
					BSa2 = 2.042000;
					BSb2 = 0.000000;
					BSc2 = 0.707000;
					BSd2 = 0.178000;
					BSe2 = 0.028700;
					BSf2 = -0.004570;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.137000;
					BSb3 = 0.359000;
					BSc3 = 0.000000;
					BSe3 = 0.054000;
					BSf3 = 0.025600;
					BSg3 = 0.000000;
					BSh3 = 0.671000;
					BSi3 = 0.000000;
					BSa4 = 2.816000;
					BSb4 = 0.031700;
					BSc4 = 0.627000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 46)	{
					BSa1 = 0.025200;
					BSb1 = 0.000000;
					BSc1 = 0.005710;
					BSd1 = -0.016900;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.215000;
					BSr2 = 0.403000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.526000;
					BSt2 = 0.013700;
					BSy2 = 0.000000;
					BSa2 = 1.925000;
					BSb2 = 0.000000;
					BSc2 = 0.767000;
					BSd2 = 0.184000;
					BSe2 = 0.032400;
					BSf2 = -0.005470;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.004970;
					BSb3 = 0.377000;
					BSc3 = 0.000000;
					BSe3 = 0.077200;
					BSf3 = 0.019500;
					BSg3 = 0.000000;
					BSh3 = 0.651000;
					BSi3 = 0.000000;
					BSa4 = 2.708000;
					BSb4 = 0.073500;
					BSc4 = 0.546000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 47)	{
					BSa1 = 0.063300;
					BSb1 = 0.000000;
					BSc1 = 0.014600;
					BSd1 = -0.021100;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.156000;
					BSr2 = 0.401000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.527000;
					BSt2 = 0.013500;
					BSy2 = 0.000000;
					BSa2 = 1.931000;
					BSb2 = 0.000000;
					BSc2 = 0.752000;
					BSd2 = 0.182000;
					BSe2 = 0.035200;
					BSf2 = -0.006660;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.006420;
					BSb3 = 0.370000;
					BSc3 = 0.000000;
					BSe3 = 0.092500;
					BSf3 = 0.019200;
					BSg3 = 0.000000;
					BSh3 = 0.629000;
					BSi3 = 0.000000;
					BSa4 = 2.713000;
					BSb4 = 0.080500;
					BSc4 = 0.508000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 48)	{
					BSa1 = 0.036600;
					BSb1 = 0.000000;
					BSc1 = 0.002540;
					BSd1 = -0.056900;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.200000;
					BSr2 = 0.408000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.552000;
					BSt2 = 0.011900;
					BSy2 = 0.000000;
					BSa2 = 1.815000;
					BSb2 = 0.000000;
					BSc2 = 0.790000;
					BSd2 = 0.221000;
					BSe2 = 0.075400;
					BSf2 = -0.021200;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.091600;
					BSb3 = 0.285000;
					BSc3 = 0.000000;
					BSe3 = 0.146000;
					BSf3 = 0.014000;
					BSg3 = 0.000000;
					BSh3 = 0.736000;
					BSi3 = 0.000000;
					BSa4 = 2.241000;
					BSb4 = 0.416000;
					BSc4 = 0.361000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 49) {
					BSa1 = 0.093900;
					BSb1 = 0.000000;
					BSc1 = 0.003020;
					BSd1 = -0.058100;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.176000;
					BSr2 = 0.408000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.555000;
					BSt2 = 0.012200;
					BSy2 = 0.000000;
					BSa2 = 1.857000;
					BSb2 = 0.000000;
					BSc2 = 0.756000;
					BSd2 = 0.213000;
					BSe2 = 0.065700;
					BSf2 = -0.017300;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.071900;
					BSb3 = 0.280000;
					BSc3 = 0.000000;
					BSe3 = 0.153000;
					BSf3 = 0.015500;
					BSg3 = 0.000000;
					BSh3 = 0.696000;
					BSi3 = 0.000000;
					BSa4 = 2.352000;
					BSb4 = 0.340000;
					BSc4 = 0.362000;
					ICol = 3;
				}
				else {
					// Config < 50 but not in list
					BSFndConfigFound = false;
				}
			}
			else if (BSFoundationCfg < 100)	{
				if (BSFoundationCfg == 50)	{
					BSa1 = 0.025300;
					BSb1 = 0.000000;
					BSc1 = 0.004640;
					BSd1 = 0.003570;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.283000;
					BSr2 = 0.406000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.549000;
					BSt2 = 0.012200;
					BSy2 = 0.000000;
					BSa2 = 2.004000;
					BSb2 = 0.000000;
					BSc2 = 0.739000;
					BSd2 = 0.203000;
					BSe2 = 0.033500;
					BSf2 = -0.006230;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.067700;
					BSb3 = 0.315000;
					BSc3 = 0.000000;
					BSe3 = 0.069100;
					BSf3 = 0.020400;
					BSg3 = 0.000000;
					BSh3 = 0.737000;
					BSi3 = 0.000000;
					BSa4 = 2.736000;
					BSb4 = 0.045400;
					BSc4 = 0.729000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 51)	{
					BSa1 = 0.057600;
					BSb1 = 0.000000;
					BSc1 = 0.009710;
					BSd1 = 0.006590;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.306000;
					BSr2 = 0.401000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.544000;
					BSt2 = 0.012000;
					BSy2 = 0.000000;
					BSa2 = 2.062000;
					BSb2 = 0.000000;
					BSc2 = 0.706000;
					BSd2 = 0.201000;
					BSe2 = 0.034400;
					BSf2 = -0.006960;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.110000;
					BSb3 = 0.308000;
					BSc3 = 0.000000;
					BSe3 = 0.073100;
					BSf3 = 0.020100;
					BSg3 = 0.000000;
					BSh3 = 0.730000;
					BSi3 = 0.000000;
					BSa4 = 2.770000;
					BSb4 = 0.039300;
					BSc4 = 0.732000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 52)	{
					BSa1 = 0.017300;
					BSb1 = 0.000000;
					BSc1 = 0.009190;
					BSd1 = 0.083200;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.583000;
					BSr2 = 0.113000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.309000;
					BSt2 = 0.007560;
					BSy2 = 0.000000;
					BSa2 = 2.130000;
					BSb2 = 0.000000;
					BSc2 = 0.423000;
					BSd2 = 0.156000;
					BSe2 = 0.120000;
					BSf2 = 0.094800;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.243000;
					BSb3 = 0.010100;
					BSc3 = 0.000000;
					BSe3 = -0.003700;
					BSf3 = 0.226000;
					BSg3 = 0.000000;
					BSh3 = 0.610000;
					BSi3 = 0.000000;
					BSa4 = 2.644000;
					BSb4 = 0.186000;
					BSc4 = -0.296000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 53)	{
					BSa1 = 0.040200;
					BSb1 = 0.000000;
					BSc1 = 0.019900;
					BSd1 = 0.088900;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.611000;
					BSr2 = 0.113000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.308000;
					BSt2 = 0.007390;
					BSy2 = 0.000000;
					BSa2 = 2.169000;
					BSb2 = 0.000000;
					BSc2 = 0.402000;
					BSd2 = 0.153000;
					BSe2 = 0.118000;
					BSf2 = 0.094200;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.279000;
					BSb3 = 0.001160;
					BSc3 = 0.000000;
					BSe3 = 0.003570;
					BSf3 = 0.226000;
					BSg3 = 0.000000;
					BSh3 = 0.612000;
					BSi3 = 0.000000;
					BSa4 = 2.634000;
					BSb4 = 0.220000;
					BSc4 = -0.263000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 54)	{
					BSa1 = 0.019000;
					BSb1 = 0.000000;
					BSc1 = 0.010600;
					BSd1 = 0.054100;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.536000;
					BSr2 = 0.114000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.311000;
					BSt2 = 0.007740;
					BSy2 = 0.000000;
					BSa2 = 2.033000;
					BSb2 = 0.000000;
					BSc2 = 0.452000;
					BSd2 = 0.158000;
					BSe2 = 0.133000;
					BSf2 = 0.093800;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.149000;
					BSb3 = 0.030600;
					BSc3 = 0.000000;
					BSe3 = 0.022700;
					BSf3 = 0.216000;
					BSg3 = 0.000000;
					BSh3 = 0.599000;
					BSi3 = 0.000000;
					BSa4 = 2.709000;
					BSb4 = 0.071500;
					BSc4 = -0.505000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 55)	{
					BSa1 = 0.048200;
					BSb1 = 0.000000;
					BSc1 = 0.026700;
					BSd1 = 0.051400;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.484000;
					BSr2 = 0.114000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.311000;
					BSt2 = 0.007770;
					BSy2 = 0.000000;
					BSa2 = 2.026000;
					BSb2 = 0.000000;
					BSc2 = 0.441000;
					BSd2 = 0.156000;
					BSe2 = 0.138000;
					BSf2 = 0.094600;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.151000;
					BSb3 = 0.025500;
					BSc3 = 0.000000;
					BSe3 = 0.036700;
					BSf3 = 0.213000;
					BSg3 = 0.000000;
					BSh3 = 0.601000;
					BSi3 = 0.000000;
					BSa4 = 2.720000;
					BSb4 = 0.074200;
					BSc4 = -0.500000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 56)	{
					BSa1 = 0.028000;
					BSb1 = 0.000000;
					BSc1 = 0.009720;
					BSd1 = 0.018900;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.409000;
					BSr2 = 0.115000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.324000;
					BSt2 = 0.006290;
					BSy2 = 0.000000;
					BSa2 = 1.837000;
					BSb2 = 0.000000;
					BSc2 = 0.476000;
					BSd2 = 0.177000;
					BSe2 = 0.190000;
					BSf2 = 0.080200;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.010300;
					BSb3 = 0.030500;
					BSc3 = 0.000000;
					BSe3 = 0.111000;
					BSf3 = 0.148000;
					BSg3 = 0.000000;
					BSh3 = 0.657000;
					BSi3 = 0.000000;
					BSa4 = 2.549000;
					BSb4 = 0.168000;
					BSc4 = 2.465000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 57)	{
					BSa1 = 0.072300;
					BSb1 = 0.000000;
					BSc1 = 0.021300;
					BSd1 = 0.021400;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.385000;
					BSr2 = 0.116000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.328000;
					BSt2 = 0.006910;
					BSy2 = 0.000000;
					BSa2 = 1.870000;
					BSb2 = 0.000000;
					BSc2 = 0.452000;
					BSd2 = 0.174000;
					BSe2 = 0.180000;
					BSf2 = 0.087200;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.037500;
					BSb3 = 0.019900;
					BSc3 = 0.000000;
					BSe3 = 0.110000;
					BSf3 = 0.154000;
					BSg3 = 0.000000;
					BSh3 = 0.654000;
					BSi3 = 0.000000;
					BSa4 = 2.670000;
					BSb4 = -0.018700;
					BSc4 = 1.000000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 58)	{
					BSa1 = 0.018800;
					BSb1 = 0.000000;
					BSc1 = 0.009840;
					BSd1 = 0.081700;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.620000;
					BSr2 = 0.115000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.326000;
					BSt2 = 0.006580;
					BSy2 = 0.000000;
					BSa2 = 2.129000;
					BSb2 = 0.000000;
					BSc2 = 0.428000;
					BSd2 = 0.171000;
					BSe2 = 0.121000;
					BSf2 = 0.091400;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.201000;
					BSb3 = 0.021800;
					BSc3 = 0.000000;
					BSe3 = 0.009230;
					BSf3 = 0.183000;
					BSg3 = 0.000000;
					BSh3 = 0.679000;
					BSi3 = 0.000000;
					BSa4 = 2.977000;
					BSb4 = -0.262000;
					BSc4 = 1.000000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 59) {
					BSa1 = 0.043500;
					BSb1 = 0.000000;
					BSc1 = 0.021300;
					BSd1 = 0.086800;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.659000;
					BSr2 = 0.114000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.323000;
					BSt2 = 0.006260;
					BSy2 = 0.000000;
					BSa2 = 2.188000;
					BSb2 = 0.000000;
					BSc2 = 0.401000;
					BSd2 = 0.169000;
					BSe2 = 0.117000;
					BSf2 = 0.090700;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.242000;
					BSb3 = 0.010700;
					BSc3 = 0.000000;
					BSe3 = 0.015100;
					BSf3 = 0.184000;
					BSg3 = 0.000000;
					BSh3 = 0.681000;
					BSi3 = 0.000000;
					BSa4 = 2.562000;
					BSb4 = 0.256000;
					BSc4 = -0.246000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 60)	{
					BSa1 = 0.010800;
					BSb1 = 0.000000;
					BSc1 = 0.007830;
					BSd1 = 0.490000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.620000;
					BSr2 = 0.098400;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.230000;
					BSt2 = 0.008990;
					BSy2 = 0.000000;
					BSa2 = 2.185000;
					BSb2 = 0.000000;
					BSc2 = 0.405000;
					BSd2 = 0.177000;
					BSe2 = 0.196000;
					BSf2 = 0.111000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.000555;
					BSb3 = -0.040600;
					BSc3 = 0.000000;
					BSe3 = 0.143000;
					BSf3 = 0.241000;
					BSg3 = 0.000000;
					BSh3 = 0.556000;
					BSi3 = 0.000000;
					BSa4 = 2.472000;
					BSb4 = 0.226000;
					BSc4 = 0.838000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 61)	{
					BSa1 = 0.024800;
					BSb1 = 0.000000;
					BSc1 = 0.017200;
					BSd1 = 0.563000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.601000;
					BSr2 = 0.098000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.227000;
					BSt2 = 0.008790;
					BSy2 = 0.000000;
					BSa2 = 2.171000;
					BSb2 = 0.000000;
					BSc2 = 0.400000;
					BSd2 = 0.174000;
					BSe2 = 0.204000;
					BSf2 = 0.109000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.004710;
					BSb3 = -0.039800;
					BSc3 = 0.000000;
					BSe3 = 0.168000;
					BSf3 = 0.236000;
					BSg3 = 0.000000;
					BSh3 = 0.560000;
					BSi3 = 0.000000;
					BSa4 = 2.442000;
					BSb4 = 0.270000;
					BSc4 = 0.623000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 62)	{
					BSa1 = 0.011400;
					BSb1 = 0.000000;
					BSc1 = 0.009890;
					BSd1 = 0.467000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.626000;
					BSr2 = 0.099500;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.237000;
					BSt2 = 0.008750;
					BSy2 = 0.000000;
					BSa2 = 2.162000;
					BSb2 = 0.000000;
					BSc2 = 0.418000;
					BSd2 = 0.181000;
					BSe2 = 0.191000;
					BSf2 = 0.113000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.001230;
					BSb3 = -0.037400;
					BSc3 = 0.000000;
					BSe3 = 0.103000;
					BSf3 = 0.246000;
					BSg3 = 0.000000;
					BSh3 = 0.564000;
					BSi3 = 0.000000;
					BSa4 = 2.385000;
					BSb4 = 0.272000;
					BSc4 = 0.720000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 63)	{
					BSa1 = 0.028200;
					BSb1 = 0.000000;
					BSc1 = 0.024900;
					BSd1 = 0.524000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.594000;
					BSr2 = 0.099400;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.237000;
					BSt2 = 0.008870;
					BSy2 = 0.000000;
					BSa2 = 2.141000;
					BSb2 = 0.000000;
					BSc2 = 0.417000;
					BSd2 = 0.182000;
					BSe2 = 0.202000;
					BSf2 = 0.112000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.002830;
					BSb3 = -0.034000;
					BSc3 = 0.000000;
					BSe3 = 0.110000;
					BSf3 = 0.240000;
					BSg3 = 0.000000;
					BSh3 = 0.579000;
					BSi3 = 0.000000;
					BSa4 = 2.332000;
					BSb4 = 0.329000;
					BSc4 = 0.546000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 64)	{
					BSa1 = 0.018700;
					BSb1 = 0.000000;
					BSc1 = 0.010500;
					BSd1 = 0.400000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.658000;
					BSr2 = 0.104000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.270000;
					BSt2 = 0.007570;
					BSy2 = 0.000000;
					BSa2 = 2.110000;
					BSb2 = 0.000000;
					BSc2 = 0.429000;
					BSd2 = 0.206000;
					BSe2 = 0.187000;
					BSf2 = 0.110000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.008890;
					BSb3 = -0.010100;
					BSc3 = 0.000000;
					BSe3 = 0.075700;
					BSf3 = 0.179000;
					BSg3 = 0.000000;
					BSh3 = 0.673000;
					BSi3 = 0.000000;
					BSa4 = 1.849000;
					BSb4 = 0.674000;
					BSc4 = 0.432000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 65)	{
					BSa1 = 0.049100;
					BSb1 = 0.000000;
					BSc1 = 0.020900;
					BSd1 = 0.457000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.616000;
					BSr2 = 0.103000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.265000;
					BSt2 = 0.008140;
					BSy2 = 0.000000;
					BSa2 = 2.099000;
					BSb2 = 0.000000;
					BSc2 = 0.421000;
					BSd2 = 0.203000;
					BSe2 = 0.203000;
					BSf2 = 0.108000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.009400;
					BSb3 = -0.011100;
					BSc3 = 0.000000;
					BSe3 = 0.088700;
					BSf3 = 0.177000;
					BSg3 = 0.000000;
					BSh3 = 0.675000;
					BSi3 = 0.000000;
					BSa4 = 1.919000;
					BSb4 = 0.635000;
					BSc4 = 0.417000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 66)	{
					BSa1 = 0.012300;
					BSb1 = 0.000000;
					BSc1 = 0.008490;
					BSd1 = 0.477000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.684000;
					BSr2 = 0.101000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.251000;
					BSt2 = 0.008190;
					BSy2 = 0.000000;
					BSa2 = 2.220000;
					BSb2 = 0.000000;
					BSc2 = 0.404000;
					BSd2 = 0.193000;
					BSe2 = 0.188000;
					BSf2 = 0.109000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.001580;
					BSb3 = -0.028100;
					BSc3 = 0.000000;
					BSe3 = 0.128000;
					BSf3 = 0.203000;
					BSg3 = 0.000000;
					BSh3 = 0.588000;
					BSi3 = 0.000000;
					BSa4 = 2.395000;
					BSb4 = 0.253000;
					BSc4 = 0.893000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 67)	{
					BSa1 = 0.028200;
					BSb1 = 0.000000;
					BSc1 = 0.018400;
					BSd1 = 0.543000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.663000;
					BSr2 = 0.100000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.246000;
					BSt2 = 0.007990;
					BSy2 = 0.000000;
					BSa2 = 2.207000;
					BSb2 = 0.000000;
					BSc2 = 0.397000;
					BSd2 = 0.189000;
					BSe2 = 0.196000;
					BSf2 = 0.107000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.005900;
					BSb3 = -0.028000;
					BSc3 = 0.000000;
					BSe3 = 0.155000;
					BSf3 = 0.197000;
					BSg3 = 0.000000;
					BSh3 = 0.589000;
					BSi3 = 0.000000;
					BSa4 = 2.373000;
					BSb4 = 0.292000;
					BSc4 = 0.674000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 68)	{
					BSa1 = 0.066000;
					BSb1 = 0.384000;
					BSc1 = 0.021000;
					BSd1 = 1.072000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.080000;
					BSq2 = 0.134000;
					BSr2 = 0.352000;
					BSu2 = 0.118000;
					BSv2 = 0.620000;
					BSw2 = 0.690000;
					BSx2 = 0.500000;
					BSs2 = -0.012000;
					BSt2 = 0.024000;
					BSy2 = 0.176000;
					BSa2 = 0.948000;
					BSb2 = 0.449000;
					BSc2 = 0.262000;
					BSd2 = 0.372000;
					BSe2 = 3.845000;
					BSf2 = -0.116000;
					BSg2 = -1.323000;
					BSh2 = 0.000000;
					BSa3 = 0.040000;
					BSb3 = 0.275000;
					BSc3 = -0.123000;
					BSe3 = 0.060000;
					BSf3 = 0.052000;
					BSg3 = 0.122000;
					BSh3 = 0.948000;
					BSi3 = 0.000000;
					BSa4 = 1.962000;
					BSb4 = 0.775000;
					BSc4 = 0.349000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 69)	{
					BSa1 = 0.095000;
					BSb1 = 0.372000;
					BSc1 = 0.034000;
					BSd1 = 1.030000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.048000;
					BSq2 = 0.596000;
					BSr2 = -0.037000;
					BSu2 = 0.672000;
					BSv2 = 0.620000;
					BSw2 = -3.480000;
					BSx2 = 0.735000;
					BSs2 = 0.631000;
					BSt2 = 0.006000;
					BSy2 = -0.151000;
					BSa2 = 3.009000;
					BSb2 = 0.159000;
					BSc2 = 0.519000;
					BSd2 = 0.296000;
					BSe2 = 0.285000;
					BSf2 = 0.049000;
					BSg2 = 0.035000;
					BSh2 = 0.000000;
					BSa3 = 0.014000;
					BSb3 = -0.002000;
					BSc3 = -0.025000;
					BSe3 = 0.099000;
					BSf3 = 0.192000;
					BSg3 = 0.074000;
					BSh3 = 0.779000;
					BSi3 = 0.000000;
					BSa4 = -5.952000;
					BSb4 = 8.699000;
					BSc4 = 0.024000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 70)	{
					BSa1 = -0.084000;
					BSb1 = 0.748000;
					BSc1 = 0.005000;
					BSd1 = 0.883000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.355000;
					BSq2 = -0.024000;
					BSr2 = 0.021000;
					BSu2 = 0.619000;
					BSv2 = 0.620000;
					BSw2 = 2.348000;
					BSx2 = 0.500000;
					BSs2 = -0.214000;
					BSt2 = 0.019000;
					BSy2 = 0.095000;
					BSa2 = 2.091000;
					BSb2 = 0.343000;
					BSc2 = 0.571000;
					BSd2 = 0.285000;
					BSe2 = 0.469000;
					BSf2 = 0.174000;
					BSg2 = -0.027000;
					BSh2 = 0.000000;
					BSa3 = 0.117000;
					BSb3 = 0.068000;
					BSc3 = -0.076000;
					BSe3 = 0.048000;
					BSf3 = 0.214000;
					BSg3 = 0.208000;
					BSh3 = 0.683000;
					BSi3 = 0.000000;
					BSa4 = 1.793000;
					BSb4 = 1.096000;
					BSc4 = 0.092000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 71)	{
					BSa1 = -0.037000;
					BSb1 = 0.765000;
					BSc1 = 0.002900;
					BSd1 = 0.888130;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.132000;
					BSq2 = 1.189000;
					BSr2 = -0.072000;
					BSu2 = 0.584000;
					BSv2 = 0.620000;
					BSw2 = -2.634000;
					BSx2 = 0.735000;
					BSs2 = 0.822000;
					BSt2 = 0.006000;
					BSy2 = -0.239000;
					BSa2 = 3.472000;
					BSb2 = 0.199000;
					BSc2 = 0.533000;
					BSd2 = 0.283000;
					BSe2 = 0.161000;
					BSf2 = 0.035000;
					BSg2 = 0.035000;
					BSh2 = 0.000000;
					BSa3 = 0.118000;
					BSb3 = 0.068000;
					BSc3 = -0.078000;
					BSe3 = 0.044000;
					BSf3 = 0.214000;
					BSg3 = 0.218000;
					BSh3 = 0.683000;
					BSi3 = 0.000000;
					BSa4 = 1.695000;
					BSb4 = 1.197000;
					BSc4 = 0.084000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 72)	{
					BSa1 = -0.020400;
					BSb1 = 0.760000;
					BSc1 = 0.100000;
					BSd1 = 0.861000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.008000;
					BSq2 = 0.443000;
					BSr2 = -0.029000;
					BSu2 = 0.690000;
					BSv2 = 0.620000;
					BSw2 = -4.225000;
					BSx2 = 0.735000;
					BSs2 = 0.556000;
					BSt2 = 0.006500;
					BSy2 = -0.117000;
					BSa2 = 3.150000;
					BSb2 = 0.222000;
					BSc2 = 0.507000;
					BSd2 = 0.254000;
					BSe2 = 0.330000;
					BSf2 = 0.054400;
					BSg2 = 0.022400;
					BSh2 = 0.000000;
					BSa3 = -0.011000;
					BSb3 = -0.031000;
					BSc3 = -0.006000;
					BSe3 = 0.193000;
					BSf3 = 0.277000;
					BSg3 = 0.145000;
					BSh3 = 0.664000;
					BSi3 = 0.000000;
					BSa4 = 2.523000;
					BSb4 = 0.385000;
					BSc4 = 0.276000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 73)	{
					BSa1 = -0.071230;
					BSb1 = 0.760347;
					BSc1 = 0.131294;
					BSd1 = 0.785350;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.036898;
					BSq2 = 0.773493;
					BSr2 = -0.027347;
					BSu2 = 0.450354;
					BSv2 = 0.620000;
					BSw2 = -4.899600;
					BSx2 = 0.735000;
					BSs2 = 0.590526;
					BSt2 = 0.007582;
					BSy2 = -0.086467;
					BSa2 = 4.715540;
					BSb2 = 0.113200;
					BSc2 = 0.551280;
					BSd2 = 0.319280;
					BSe2 = 0.126000;
					BSf2 = 0.027226;
					BSg2 = 0.027579;
					BSh2 = 0.000000;
					BSa3 = 0.294070;
					BSb3 = 0.081500;
					BSc3 = -0.106280;
					BSe3 = 0.062620;
					BSf3 = 0.243660;
					BSg3 = 0.144810;
					BSh3 = 0.724765;
					BSi3 = 0.000000;
					BSa4 = 2.876000;
					BSb4 = 0.064696;
					BSc4 = 0.276000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 74)	{
					BSa1 = 0.198990;
					BSb1 = 2.918970;
					BSc1 = 0.006260;
					BSd1 = -0.002425;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.010190;
					BSr2 = 0.020820;
					BSu2 = 0.881412;
					BSv2 = 0.620000;
					BSw2 = 1.954236;
					BSx2 = 0.500000;
					BSs2 = -0.275610;
					BSt2 = 0.020197;
					BSy2 = 0.086238;
					BSa2 = 2.211346;
					BSb2 = 0.331075;
					BSc2 = 0.530859;
					BSd2 = 0.316900;
					BSe2 = 0.228122;
					BSf2 = 0.228748;
					BSg2 = 0.045537;
					BSh2 = 0.000000;
					BSa3 = 0.704441;
					BSb3 = 0.061763;
					BSc3 = -0.222700;
					BSe3 = -0.077969;
					BSf3 = 0.201850;
					BSg3 = 0.246550;
					BSh3 = 0.666425;
					BSi3 = 0.000000;
					BSa4 = 3.032800;
					BSb4 = -0.049830;
					BSc4 = 0.411000;
					ICol = 4;
				}
				else if (BSFoundationCfg == 75)	{
					BSa1 = 0.124817;
					BSb1 = 2.949200;
					BSc1 = 0.006070;
					BSd1 = -0.002375;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = 0.010425;
					BSr2 = 0.020505;
					BSu2 = 0.885840;
					BSv2 = 0.620000;
					BSw2 = 1.985315;
					BSx2 = 0.500000;
					BSs2 = -0.274739;
					BSt2 = 0.020236;
					BSy2 = 0.083996;
					BSa2 = 2.217098;
					BSb2 = 0.335890;
					BSc2 = 0.528369;
					BSd2 = 0.317068;
					BSe2 = 0.231890;
					BSf2 = 0.228317;
					BSg2 = 0.041176;
					BSh2 = 0.000000;
					BSa3 = 0.701859;
					BSb3 = 0.061754;
					BSc3 = -0.215070;
					BSe3 = -0.078714;
					BSf3 = 0.201920;
					BSg3 = 0.278666;
					BSh3 = 0.666900;
					BSi3 = 0.000000;
					BSa4 = 3.036870;
					BSb4 = -0.052290;
					BSc4 = 0.426000;
					ICol = 4;
				}
				else if (BSFoundationCfg == 76)	{
					BSa1 = -0.083310;
					BSb1 = 0.746580;
					BSc1 = 0.006320;
					BSd1 = 0.878840;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.353700;
					BSq2 = -0.069721;
					BSr2 = 0.046983;
					BSu2 = 4.049199;
					BSv2 = 0.620000;
					BSw2 = -0.029300;
					BSx2 = 0.500000;
					BSs2 = -0.023519;
					BSt2 = 0.012950;
					BSy2 = 0.328888;
					BSa2 = 1.520370;
					BSb2 = 0.466830;
					BSc2 = 0.724090;
					BSd2 = 0.110352;
					BSe2 = 0.561002;
					BSf2 = -0.009356;
					BSg2 = -0.165070;
					BSh2 = 0.000000;
					BSa3 = -0.418600;
					BSb3 = 0.277700;
					BSc3 = 0.418100;
					BSe3 = 0.307104;
					BSf3 = 0.135589;
					BSg3 = -0.089250;
					BSh3 = 0.747110;
					BSi3 = 0.000000;
					BSa4 = 2.854400;
					BSb4 = 0.072350;
					BSc4 = 1.299000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 77)	{
					BSa1 = -0.038213;
					BSb1 = 0.763880;
					BSc1 = 0.003235;
					BSd1 = 0.887379;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.119581;
					BSq2 = 0.660969;
					BSr2 = 0.168863;
					BSu2 = 0.112584;
					BSv2 = 0.620000;
					BSw2 = -0.020200;
					BSx2 = 0.500000;
					BSs2 = -0.206723;
					BSt2 = 0.015950;
					BSy2 = 0.095326;
					BSa2 = 2.152850;
					BSb2 = 0.024900;
					BSc2 = 0.455866;
					BSd2 = 0.360834;
					BSe2 = 0.324333;
					BSf2 = 0.359352;
					BSg2 = -0.118260;
					BSh2 = 0.000000;
					BSa3 = -0.418948;
					BSb3 = 0.278115;
					BSc3 = 0.418326;
					BSe3 = 0.304350;
					BSf3 = 0.135980;
					BSg3 = -0.081449;
					BSh3 = 0.746450;
					BSi3 = 0.000000;
					BSa4 = 2.854540;
					BSb4 = 0.074000;
					BSc4 = 1.254000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 78)	{
					BSa1 = -0.030788;
					BSb1 = 0.650499;
					BSc1 = -0.000330;
					BSd1 = 0.810000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075968;
					BSq2 = -0.162143;
					BSr2 = 0.215933;
					BSu2 = 0.694206;
					BSv2 = 0.620000;
					BSw2 = 0.293126;
					BSx2 = 0.500000;
					BSs2 = 0.067465;
					BSt2 = 0.014496;
					BSy2 = 0.206205;
					BSa2 = 0.615051;
					BSb2 = 0.713764;
					BSc2 = 0.576454;
					BSd2 = -0.010456;
					BSe2 = 0.317440;
					BSf2 = -0.019935;
					BSg2 = -0.097595;
					BSh2 = 0.000000;
					BSa3 = -0.267450;
					BSb3 = 0.224320;
					BSc3 = 0.277777;
					BSe3 = 0.207546;
					BSf3 = 0.046110;
					BSg3 = 0.013254;
					BSh3 = 0.677812;
					BSi3 = 0.000000;
					BSa4 = 2.560000;
					BSb4 = 0.347400;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 79)	{
					BSa1 = -0.012492;
					BSb1 = 0.649800;
					BSc1 = 0.003403;
					BSd1 = 0.804034;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.001523;
					BSq2 = 0.393200;
					BSr2 = -0.024024;
					BSu2 = 0.681964;
					BSv2 = 0.620000;
					BSw2 = -4.663550;
					BSx2 = 0.735000;
					BSs2 = 0.561570;
					BSt2 = 0.006560;
					BSy2 = -0.108203;
					BSa2 = 3.118051;
					BSb2 = 0.213268;
					BSc2 = 0.488352;
					BSd2 = 0.263121;
					BSe2 = 0.315049;
					BSf2 = 0.052232;
					BSg2 = 0.018485;
					BSh2 = 0.000000;
					BSa3 = 0.006841;
					BSb3 = -0.040390;
					BSc3 = -0.006733;
					BSe3 = 0.134274;
					BSf3 = 0.253267;
					BSg3 = 0.145490;
					BSh3 = 0.587305;
					BSi3 = 0.000000;
					BSa4 = 2.466150;
					BSb4 = 0.485450;
					BSc4 = 1.000000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 80)	{
					BSa1 = -0.013814;
					BSb1 = 0.650135;
					BSc1 = 0.002998;
					BSd1 = 0.808695;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.003843;
					BSq2 = 0.222760;
					BSr2 = 0.149651;
					BSu2 = 0.388961;
					BSv2 = 0.620000;
					BSw2 = -0.204761;
					BSx2 = 0.500000;
					BSs2 = -0.454127;
					BSt2 = 0.017890;
					BSy2 = 0.205281;
					BSa2 = 2.210061;
					BSb2 = -0.005805;
					BSc2 = 0.520561;
					BSd2 = 0.327860;
					BSe2 = 0.770500;
					BSf2 = 0.225836;
					BSg2 = -0.285055;
					BSh2 = 0.000000;
					BSa3 = -0.470070;
					BSb3 = 0.135941;
					BSc3 = 0.495440;
					BSe3 = 0.397302;
					BSf3 = 0.174450;
					BSg3 = -0.176874;
					BSh3 = 0.640170;
					BSi3 = 0.000000;
					BSa4 = 2.671512;
					BSb4 = 0.236990;
					BSc4 = 1.000000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 81)	{
					BSa1 = 0.037138;
					BSb1 = 0.000000;
					BSc1 = 0.000693;
					BSd1 = -0.100380;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.735360;
					BSr2 = 0.096844;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.243830;
					BSt2 = 0.000000;
					BSy2 = 0.000000;
					BSa2 = 2.378200;
					BSb2 = 0.000000;
					BSc2 = 0.387290;
					BSd2 = 0.176850;
					BSe2 = 0.161380;
					BSf2 = 0.108647;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.006050;
					BSb3 = -0.060180;
					BSc3 = 0.000000;
					BSe3 = 0.137650;
					BSf3 = 0.276020;
					BSg3 = 0.000000;
					BSh3 = 0.483710;
					BSi3 = 0.000000;
					BSa4 = 2.657300;
					BSb4 = 0.156290;
					BSc4 = 0.000000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 82)	{
					BSa1 = 0.085100;
					BSb1 = 0.000000;
					BSc1 = 0.001800;
					BSd1 = -0.090300;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -1.556000;
					BSr2 = 0.099600;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.254100;
					BSt2 = 0.008000;
					BSy2 = 0.000000;
					BSa2 = 2.227000;
					BSb2 = 0.000000;
					BSc2 = 0.401800;
					BSd2 = 0.170700;
					BSe2 = 0.160700;
					BSf2 = 0.112400;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.047000;
					BSb3 = -0.055000;
					BSc3 = 0.000000;
					BSe3 = 0.119000;
					BSf3 = 0.269000;
					BSg3 = 0.000000;
					BSh3 = 0.500700;
					BSi3 = 0.000000;
					BSa4 = 2.737000;
					BSb4 = 0.043200;
					BSc4 = 0.000000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 83)	{
					BSa1 = 0.039517;
					BSb1 = 0.000000;
					BSc1 = 0.000185;
					BSd1 = -0.083530;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.292589;
					BSr2 = 0.417510;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.521346;
					BSt2 = 0.016540;
					BSy2 = 0.000000;
					BSa2 = 1.980200;
					BSb2 = 0.000000;
					BSc2 = 0.826860;
					BSd2 = 0.190353;
					BSe2 = 0.031480;
					BSf2 = -0.008430;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.095980;
					BSb3 = 0.454060;
					BSc3 = 0.000000;
					BSe3 = 0.116340;
					BSf3 = 0.006612;
					BSg3 = 0.000000;
					BSh3 = 0.459470;
					BSi3 = 0.000000;
					BSa4 = 2.615900;
					BSb4 = 0.159950;
					BSc4 = 0.296200;
					ICol = 3;
				}
				else if (BSFoundationCfg == 84)	{
					BSa1 = 0.088380;
					BSb1 = 0.000000;
					BSc1 = 0.002409;
					BSd1 = -0.068600;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.266900;
					BSr2 = 0.420750;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.535850;
					BSt2 = 0.015980;
					BSy2 = 0.000000;
					BSa2 = 2.018850;
					BSb2 = 0.000000;
					BSc2 = 0.801370;
					BSd2 = 0.184610;
					BSe2 = 0.019588;
					BSf2 = -0.005032;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.064300;
					BSb3 = 0.453470;
					BSc3 = 0.000000;
					BSe3 = 0.109440;
					BSf3 = 0.004804;
					BSg3 = 0.000000;
					BSh3 = 0.440260;
					BSi3 = 0.000000;
					BSa4 = 2.664500;
					BSb4 = 0.131514;
					BSc4 = 0.285300;
					ICol = 3;
				}
				else if (BSFoundationCfg == 85)	{
					BSa1 = 0.040500;
					BSb1 = 0.000000;
					BSc1 = -0.000900;
					BSd1 = -0.096100;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.127100;
					BSr2 = 0.386500;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.522600;
					BSt2 = 0.010600;
					BSy2 = 0.000000;
					BSa2 = 1.925000;
					BSb2 = 0.000000;
					BSc2 = 0.798000;
					BSd2 = 0.198600;
					BSe2 = 0.058700;
					BSf2 = -0.014200;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.105000;
					BSb3 = 0.376000;
					BSc3 = 0.000000;
					BSe3 = 0.141500;
					BSf3 = 0.027700;
					BSg3 = 0.000000;
					BSh3 = 0.539800;
					BSi3 = 0.000000;
					BSa4 = 2.401500;
					BSb4 = 0.341100;
					BSc4 = 0.246600;
					ICol = 3;
				}
				else if (BSFoundationCfg == 86)	{
					BSa1 = 0.090763;
					BSb1 = 0.000000;
					BSc1 = -0.002720;
					BSd1 = -0.081400;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.106820;
					BSr2 = 0.390350;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.535410;
					BSt2 = 0.010430;
					BSy2 = 0.000000;
					BSa2 = 1.960170;
					BSb2 = 0.000000;
					BSc2 = 0.773957;
					BSd2 = 0.192200;
					BSe2 = 0.042120;
					BSf2 = -0.009680;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.075860;
					BSb3 = 0.379150;
					BSc3 = 0.000000;
					BSe3 = 0.135920;
					BSf3 = 0.022610;
					BSg3 = 0.000000;
					BSh3 = 0.534300;
					BSi3 = 0.000000;
					BSa4 = 2.512200;
					BSb4 = 0.253350;
					BSc4 = 0.266500;
					ICol = 3;
				}
				else if (BSFoundationCfg == 87)	{
					BSa1 = -0.030700;
					BSb1 = 0.650000;
					BSc1 = -0.000065;
					BSd1 = 0.809000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075400;
					BSq2 = 0.604000;
					BSr2 = -0.038200;
					BSu2 = 0.675000;
					BSv2 = 0.620000;
					BSw2 = -3.429000;
					BSx2 = 0.730000;
					BSs2 = 0.645000;
					BSt2 = 0.007040;
					BSy2 = -0.156000;
					BSa2 = 3.101000;
					BSb2 = 0.211000;
					BSc2 = 0.491000;
					BSd2 = 0.264000;
					BSe2 = 0.274000;
					BSf2 = 0.056700;
					BSg2 = 0.024500;
					BSh2 = 0.000000;
					BSa3 = 0.038500;
					BSb3 = -0.032000;
					BSc3 = -0.020300;
					BSe3 = 0.106000;
					BSf3 = 0.248000;
					BSg3 = 0.160000;
					BSh3 = 0.572000;
					BSi3 = 0.000000;
					BSa4 = 3.109000;
					BSb4 = -0.292000;
					BSc4 = -0.313000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 88)	{
					BSa1 = -0.030500;
					BSb1 = 0.650000;
					BSc1 = -0.000394;
					BSd1 = 0.809000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075700;
					BSq2 = 8.622000;
					BSr2 = -0.503000;
					BSu2 = 0.097300;
					BSv2 = 0.620000;
					BSw2 = -2.535000;
					BSx2 = 0.730000;
					BSs2 = 1.361000;
					BSt2 = -0.003860;
					BSy2 = -0.334000;
					BSa2 = 7.643000;
					BSb2 = 0.188000;
					BSc2 = 0.775000;
					BSd2 = 0.562000;
					BSe2 = 0.071200;
					BSf2 = -0.024400;
					BSg2 = 0.018700;
					BSh2 = 0.000000;
					BSa3 = 0.085000;
					BSb3 = 0.269000;
					BSc3 = -0.142000;
					BSe3 = 0.056200;
					BSf3 = 0.088800;
					BSg3 = 0.217000;
					BSh3 = 0.632000;
					BSi3 = 0.000000;
					BSa4 = 13.055000;
					BSb4 = -10.259000;
					BSc4 = -0.015400;
					ICol = 5;
				}
				else if (BSFoundationCfg == 89)	{
					BSa1 = -0.030400;
					BSb1 = 0.650000;
					BSc1 = -0.000727;
					BSd1 = 0.810000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.076200;
					BSq2 = 7.651000;
					BSr2 = -0.505000;
					BSu2 = 0.098300;
					BSv2 = 0.620000;
					BSw2 = -2.685000;
					BSx2 = 0.730000;
					BSs2 = 1.407000;
					BSt2 = -0.000380;
					BSy2 = -0.338000;
					BSa2 = 6.665000;
					BSb2 = 0.285000;
					BSc2 = 0.781000;
					BSd2 = 0.514000;
					BSe2 = 0.072800;
					BSf2 = -0.013600;
					BSg2 = -0.014500;
					BSh2 = 0.000000;
					BSa3 = -0.236000;
					BSb3 = 0.400000;
					BSc3 = 0.194000;
					BSe3 = 0.221000;
					BSf3 = 0.045400;
					BSg3 = 0.002140;
					BSh3 = 0.662000;
					BSi3 = 0.000000;
					BSa4 = 2.488000;
					BSb4 = 0.334000;
					BSc4 = 0.333000;
					ICol = 4;
				}
				else if (BSFoundationCfg == 90)	{
					BSa1 = -0.012700;
					BSb1 = 0.649800;
					BSc1 = 0.002400;
					BSd1 = 0.803900;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.002700;
					BSq2 = -0.033900;
					BSr2 = 0.347000;
					BSu2 = 0.200000;
					BSv2 = 0.620000;
					BSw2 = 0.690000;
					BSx2 = 0.500000;
					BSs2 = -0.054600;
					BSt2 = 0.025400;
					BSy2 = 0.179400;
					BSa2 = 0.730900;
					BSb2 = 0.685000;
					BSc2 = 0.372000;
					BSd2 = 0.221000;
					BSe2 = 2.535000;
					BSf2 = -0.025100;
					BSg2 = -0.826000;
					BSh2 = 0.000000;
					BSa3 = 0.133000;
					BSb3 = 0.256000;
					BSc3 = -0.151000;
					BSe3 = 0.035700;
					BSf3 = 0.101000;
					BSg3 = 0.210000;
					BSh3 = 0.647000;
					BSi3 = 0.000000;
					BSa4 = -2.296000;
					BSb4 = 5.091000;
					BSc4 = 0.031400;
					ICol = 3;
				}
				else if (BSFoundationCfg == 91)	{
					BSa1 = -0.017500;
					BSb1 = 0.650000;
					BSc1 = 0.001940;
					BSd1 = 0.808200;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.003950;
					BSq2 = 7.542000;
					BSr2 = -0.481000;
					BSu2 = 0.105000;
					BSv2 = 0.620000;
					BSw2 = -2.685000;
					BSx2 = 0.730000;
					BSs2 = 1.359000;
					BSt2 = -0.001270;
					BSy2 = -0.322000;
					BSa2 = 6.883000;
					BSb2 = 0.259000;
					BSc2 = 0.775000;
					BSd2 = 0.501000;
					BSe2 = 0.057600;
					BSf2 = -0.009680;
					BSg2 = -0.012600;
					BSh2 = 0.000000;
					BSa3 = -0.167000;
					BSb3 = 0.429000;
					BSc3 = 0.192000;
					BSe3 = 0.197000;
					BSf3 = 0.045300;
					BSg3 = -0.008610;
					BSh3 = 0.662000;
					BSi3 = 0.000000;
					BSa4 = 2.650000;
					BSb4 = 0.198000;
					BSc4 = 0.426000;
					ICol = 2;
				}
				else if (BSFoundationCfg == 92)	{
					BSa1 = 0.031600;
					BSb1 = 0.393000;
					BSc1 = 0.030400;
					BSd1 = 1.090500;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.103000;
					BSq2 = 0.064300;
					BSr2 = 0.354000;
					BSu2 = 0.161000;
					BSv2 = 0.620000;
					BSw2 = 0.709000;
					BSx2 = 0.500000;
					BSs2 = -0.024700;
					BSt2 = 0.024200;
					BSy2 = 0.175000;
					BSa2 = 0.952000;
					BSb2 = 0.470000;
					BSc2 = 0.299000;
					BSd2 = 0.421000;
					BSe2 = 4.241000;
					BSf2 = -0.110000;
					BSg2 = -1.500000;
					BSh2 = 0.000000;
					BSa3 = 0.094600;
					BSb3 = 0.291000;
					BSc3 = -0.138000;
					BSe3 = 0.055600;
					BSf3 = 0.046500;
					BSg3 = 0.121000;
					BSh3 = 0.914000;
					BSi3 = 0.000000;
					BSa4 = 2.231000;
					BSb4 = 0.526000;
					BSc4 = 0.422000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 93)	{
					BSa1 = 0.092100;
					BSb1 = 0.371000;
					BSc1 = 0.048300;
					BSd1 = 0.962000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.030100;
					BSq2 = 0.528000;
					BSr2 = -0.033800;
					BSu2 = 0.660000;
					BSv2 = 0.620000;
					BSw2 = -3.895000;
					BSx2 = 0.730000;
					BSs2 = 0.591000;
					BSt2 = 0.006700;
					BSy2 = -0.128000;
					BSa2 = 3.058000;
					BSb2 = 0.218000;
					BSc2 = 0.535000;
					BSd2 = 0.282000;
					BSe2 = 0.280000;
					BSf2 = 0.056400;
					BSg2 = 0.033000;
					BSh2 = 0.000000;
					BSa3 = 0.000966;
					BSb3 = -0.018800;
					BSc3 = -0.011700;
					BSe3 = 0.027200;
					BSf3 = 0.246000;
					BSg3 = 0.138000;
					BSh3 = 0.665000;
					BSi3 = 0.000000;
					BSa4 = 2.069000;
					BSb4 = 0.693000;
					BSc4 = 0.272000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 94)	{
					BSa1 = 0.024600;
					BSb1 = 0.397000;
					BSc1 = 0.037800;
					BSd1 = 1.099600;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.107000;
					BSq2 = 0.846000;
					BSr2 = -0.046800;
					BSu2 = 0.568000;
					BSv2 = 0.620000;
					BSw2 = -3.490000;
					BSx2 = 0.730000;
					BSs2 = 0.708000;
					BSt2 = 0.007200;
					BSy2 = -0.160000;
					BSa2 = 3.464000;
					BSb2 = 0.130000;
					BSc2 = 0.546000;
					BSd2 = 0.315000;
					BSe2 = 0.170000;
					BSf2 = 0.037100;
					BSg2 = 0.026500;
					BSh2 = 0.000000;
					BSa3 = 0.091800;
					BSb3 = 0.053600;
					BSc3 = -0.056200;
					BSe3 = 0.097300;
					BSf3 = 0.175000;
					BSg3 = 0.068700;
					BSh3 = 0.775000;
					BSi3 = 0.000000;
					BSa4 = 2.426000;
					BSb4 = 0.363000;
					BSc4 = 0.500000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 95)	{
					BSa1 = 0.042300;
					BSb1 = 0.383000;
					BSc1 = 0.044200;
					BSd1 = 0.975000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.076600;
					BSq2 = 8.141000;
					BSr2 = -0.456000;
					BSu2 = 0.083100;
					BSv2 = 0.620000;
					BSw2 = -2.771000;
					BSx2 = 0.730000;
					BSs2 = 1.303000;
					BSt2 = -0.003200;
					BSy2 = -0.298000;
					BSa2 = 8.005000;
					BSb2 = 0.192000;
					BSc2 = 0.815000;
					BSd2 = 0.575000;
					BSe2 = 0.035200;
					BSf2 = -0.012200;
					BSg2 = 0.015200;
					BSh2 = 0.000000;
					BSa3 = 0.072600;
					BSb3 = 0.320000;
					BSc3 = -0.096700;
					BSe3 = -0.015700;
					BSf3 = 0.072800;
					BSg3 = 0.156000;
					BSh3 = 0.805000;
					BSi3 = 0.000000;
					BSa4 = 2.441000;
					BSb4 = 0.338000;
					BSc4 = 0.492000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 96)	{
					BSa1 = 0.009130;
					BSb1 = 0.406000;
					BSc1 = 0.027100;
					BSd1 = 1.115000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.128000;
					BSq2 = 0.243000;
					BSr2 = 0.341000;
					BSu2 = 0.169000;
					BSv2 = 0.620000;
					BSw2 = 0.723000;
					BSx2 = 0.500000;
					BSs2 = -0.024900;
					BSt2 = 0.023200;
					BSy2 = 0.177000;
					BSa2 = 0.993000;
					BSb2 = 0.460000;
					BSc2 = 0.319000;
					BSd2 = 0.451000;
					BSe2 = 4.253000;
					BSf2 = -0.120000;
					BSg2 = -1.519000;
					BSh2 = 0.000000;
					BSa3 = 0.130000;
					BSb3 = 0.322000;
					BSc3 = -0.166000;
					BSe3 = 0.078600;
					BSf3 = 0.045400;
					BSg3 = 0.115000;
					BSh3 = 0.864000;
					BSi3 = 0.000000;
					BSa4 = 2.300000;
					BSb4 = 0.484000;
					BSc4 = 0.383000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 97)	{
					BSa1 = -0.083500;
					BSb1 = 0.745000;
					BSc1 = 0.007900;
					BSd1 = 0.875000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.352000;
					BSq2 = -0.090000;
					BSr2 = 0.065400;
					BSu2 = 3.141000;
					BSv2 = 0.620000;
					BSw2 = -0.045500;
					BSx2 = 0.500000;
					BSs2 = -0.226000;
					BSt2 = 0.011900;
					BSy2 = 0.350000;
					BSa2 = 1.497000;
					BSb2 = 0.468000;
					BSc2 = 0.741000;
					BSd2 = 0.099100;
					BSe2 = 0.498000;
					BSf2 = 0.000540;
					BSg2 = -0.150000;
					BSh2 = 0.000000;
					BSa3 = -0.438000;
					BSb3 = 0.312000;
					BSc3 = 0.445000;
					BSe3 = 0.309000;
					BSf3 = 0.135000;
					BSg3 = -0.107000;
					BSh3 = 0.746000;
					BSi3 = 0.000000;
					BSa4 = 2.871000;
					BSb4 = 0.060200;
					BSc4 = 1.619000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 98) {
					BSa1 = -0.038600;
					BSb1 = 0.763000;
					BSc1 = 0.004000;
					BSd1 = 0.886000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.118900;
					BSq2 = 0.726000;
					BSr2 = 0.175000;
					BSu2 = 0.125000;
					BSv2 = 0.620000;
					BSw2 = -0.015000;
					BSx2 = 0.500000;
					BSs2 = -0.218000;
					BSt2 = 0.015700;
					BSy2 = 0.104000;
					BSa2 = 2.171000;
					BSb2 = 0.009670;
					BSc2 = 0.479000;
					BSd2 = 0.383000;
					BSe2 = 0.271000;
					BSf2 = 0.375000;
					BSg2 = -0.105000;
					BSh2 = 0.000000;
					BSa3 = -0.439000;
					BSb3 = 0.312000;
					BSc3 = 0.446000;
					BSe3 = 0.306000;
					BSf3 = 0.135000;
					BSg3 = -0.100000;
					BSh3 = 0.746000;
					BSi3 = 0.000000;
					BSa4 = 2.871000;
					BSb4 = 0.061400;
					BSc4 = 1.553000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 99)	{
					BSa1 = -0.038100;
					BSb1 = 0.764000;
					BSc1 = 0.002590;
					BSd1 = 0.888000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.120000;
					BSq2 = -0.005890;
					BSr2 = 0.014900;
					BSu2 = 0.592000;
					BSv2 = 0.620000;
					BSw2 = 3.604000;
					BSx2 = 0.500000;
					BSs2 = -0.226000;
					BSt2 = 0.019900;
					BSy2 = 0.087700;
					BSa2 = 2.162000;
					BSb2 = 0.363000;
					BSc2 = 0.554000;
					BSd2 = 0.310000;
					BSe2 = 0.534000;
					BSf2 = 0.149000;
					BSg2 = -0.040000;
					BSh2 = 0.000000;
					BSa3 = 0.141000;
					BSb3 = 0.097800;
					BSc3 = -0.096100;
					BSe3 = 0.043500;
					BSf3 = 0.198000;
					BSg3 = 0.221000;
					BSh3 = 0.681000;
					BSi3 = 0.000000;
					BSa4 = 2.148000;
					BSb4 = 0.753000;
					BSc4 = 0.123000;
					ICol = 8;
				}
				else // 50 < Config < 100 but not in list
				{
					BSFndConfigFound = false;
				}
			}
			else	// Config >= 100
			{
				if (BSFoundationCfg == 100) {
					BSa1 = -0.037600;
					BSb1 = 0.764000;
					BSc1 = 0.001800;
					BSd1 = 0.888000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.120000;
					BSq2 = 0.189000;
					BSr2 = 0.354000;
					BSu2 = 0.197000;
					BSv2 = 0.620000;
					BSw2 = 0.663000;
					BSx2 = 0.500000;
					BSs2 = -0.026300;
					BSt2 = 0.022300;
					BSy2 = 0.183000;
					BSa2 = 0.782000;
					BSb2 = 0.693000;
					BSc2 = 0.348000;
					BSd2 = 0.222000;
					BSe2 = 2.819000;
					BSf2 = -0.047200;
					BSg2 = -0.924000;
					BSh2 = 0.000000;
					BSa3 = 0.129000;
					BSb3 = 0.304000;
					BSc3 = -0.170000;
					BSe3 = 0.040800;
					BSf3 = 0.086300;
					BSg3 = 0.256000;
					BSh3 = 0.729000;
					BSi3 = 0.000000;
					BSa4 = 1.127000;
					BSb4 = 1.751000;
					BSc4 = 0.082300;
					ICol = 5;
				}
				else if (BSFoundationCfg == 101) {
					BSa1 = -0.037700;
					BSb1 = 0.765000;
					BSc1 = 0.001920;
					BSd1 = 0.888000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.120000;
					BSq2 = 0.038200;
					BSr2 = 0.345000;
					BSu2 = 0.213000;
					BSv2 = 0.620000;
					BSw2 = 0.666000;
					BSx2 = 0.500000;
					BSs2 = -0.034300;
					BSt2 = 0.022900;
					BSy2 = 0.179000;
					BSa2 = 0.792000;
					BSb2 = 0.632000;
					BSc2 = 0.354000;
					BSd2 = 0.195000;
					BSe2 = 2.635000;
					BSf2 = -0.023400;
					BSg2 = -0.848000;
					BSh2 = 0.000000;
					BSa3 = 0.121000;
					BSb3 = 0.273000;
					BSc3 = -0.156000;
					BSe3 = 0.047800;
					BSf3 = 0.093800;
					BSg3 = 0.253000;
					BSh3 = 0.729000;
					BSi3 = 0.000000;
					BSa4 = 0.158000;
					BSb4 = 2.721000;
					BSc4 = 0.053400;
					ICol = 5;
				}
				else if (BSFoundationCfg == 102) {
					BSa1 = -0.013200;
					BSb1 = 0.649000;
					BSc1 = 0.003230;
					BSd1 = 0.804100;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.002020;
					BSq2 = -0.070000;
					BSr2 = 0.092400;
					BSu2 = 0.577000;
					BSv2 = 0.620000;
					BSw2 = 1.153000;
					BSx2 = 0.500000;
					BSs2 = -0.039000;
					BSt2 = 0.018800;
					BSy2 = 0.127000;
					BSa2 = 0.594000;
					BSb2 = 0.703000;
					BSc2 = 0.344000;
					BSd2 = 0.073600;
					BSe2 = 0.878000;
					BSf2 = 0.109000;
					BSg2 = -0.163000;
					BSh2 = 0.000000;
					BSa3 = 0.074900;
					BSb3 = 0.085600;
					BSc3 = -0.075000;
					BSe3 = 0.029500;
					BSf3 = 0.100000;
					BSg3 = 0.240000;
					BSh3 = 0.643000;
					BSi3 = 0.000000;
					BSa4 = 3.642000;
					BSb4 = -0.819000;
					BSc4 = -0.217000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 103) {
					BSa1 = -0.012700;
					BSb1 = 0.649000;
					BSc1 = 0.003450;
					BSd1 = 0.804000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.001500;
					BSq2 = -0.126000;
					BSr2 = 0.007840;
					BSu2 = 0.809000;
					BSv2 = 0.620000;
					BSw2 = 7.878000;
					BSx2 = 0.500000;
					BSs2 = 0.357000;
					BSt2 = 0.006830;
					BSy2 = 0.010100;
					BSa2 = 2.585000;
					BSb2 = 0.243000;
					BSc2 = 0.427000;
					BSd2 = 0.228000;
					BSe2 = 0.292000;
					BSf2 = 0.039000;
					BSg2 = 0.024600;
					BSh2 = 0.000000;
					BSa3 = 0.014300;
					BSb3 = -0.035100;
					BSc3 = -0.026400;
					BSe3 = 0.092400;
					BSf3 = 0.208800;
					BSg3 = 0.192000;
					BSh3 = 0.537000;
					BSi3 = 0.000000;
					BSa4 = 3.106000;
					BSb4 = -0.278000;
					BSc4 = -0.380000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 104) {
					BSa1 = 0.019400;
					BSb1 = 0.000000;
					BSc1 = 0.001960;
					BSd1 = 0.054000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.178000;
					BSr2 = 0.398000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.517000;
					BSt2 = 0.013700;
					BSy2 = 0.000000;
					BSa2 = 1.912000;
					BSb2 = 0.000000;
					BSc2 = 0.795000;
					BSd2 = 0.185000;
					BSe2 = 0.026700;
					BSf2 = -0.004200;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.016500;
					BSb3 = 0.400000;
					BSc3 = 0.000000;
					BSe3 = 0.049400;
					BSf3 = 0.019600;
					BSg3 = 0.000000;
					BSh3 = 0.787000;
					BSi3 = 0.000000;
					BSa4 = 2.706000;
					BSb4 = 0.066300;
					BSc4 = 0.000000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 105) {
					BSa1 = 0.045300;
					BSb1 = 0.000000;
					BSc1 = 0.002370;
					BSd1 = 0.056400;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.198000;
					BSr2 = 0.396000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.513000;
					BSt2 = 0.013800;
					BSy2 = 0.000000;
					BSa2 = 1.936000;
					BSb2 = 0.000000;
					BSc2 = 0.780000;
					BSd2 = 0.184000;
					BSe2 = 0.029900;
					BSf2 = -0.005670;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = 0.031300;
					BSb3 = 0.400000;
					BSc3 = 0.000000;
					BSe3 = 0.059500;
					BSf3 = 0.017600;
					BSg3 = 0.000000;
					BSh3 = 0.791000;
					BSi3 = 0.000000;
					BSa4 = 2.721000;
					BSb4 = 0.068700;
					BSc4 = 0.000000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 106) {
					BSa1 = 0.017800;
					BSb1 = 0.000000;
					BSc1 = 0.001760;
					BSd1 = 0.092500;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.033000;
					BSr2 = 0.375000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.461000;
					BSt2 = 0.015900;
					BSy2 = 0.000000;
					BSa2 = 1.831000;
					BSb2 = 0.000000;
					BSc2 = 0.774000;
					BSd2 = 0.190000;
					BSe2 = 0.056300;
					BSf2 = -0.010000;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.021000;
					BSb3 = 0.331000;
					BSc3 = 0.000000;
					BSe3 = 0.084900;
					BSf3 = 0.038200;
					BSg3 = 0.000000;
					BSh3 = 0.744000;
					BSi3 = 0.000000;
					BSa4 = 2.583000;
					BSb4 = 0.154000;
					BSc4 = 0.000000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 107) {
					BSa1 = 0.041000;
					BSb1 = 0.000000;
					BSc1 = 0.002200;
					BSd1 = 0.101000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.000000;
					BSq2 = -3.047500;
					BSr2 = 0.374000;
					BSu2 = 1.000000;
					BSv2 = 0.000000;
					BSw2 = 1.000000;
					BSx2 = 0.000000;
					BSs2 = 0.458000;
					BSt2 = 0.015900;
					BSy2 = 0.000000;
					BSa2 = 1.852000;
					BSb2 = 0.000000;
					BSc2 = 0.760000;
					BSd2 = 0.190000;
					BSe2 = 0.060100;
					BSf2 = -0.011800;
					BSg2 = 0.000000;
					BSh2 = 0.000000;
					BSa3 = -0.013300;
					BSb3 = 0.432000;
					BSc3 = 0.000000;
					BSe3 = 0.098600;
					BSf3 = 0.035200;
					BSg3 = 0.000000;
					BSh3 = 0.749000;
					BSi3 = 0.000000;
					BSa4 = 2.597000;
					BSb4 = 0.159000;
					BSc4 = 0.000000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 108) {
					BSa1 = -0.012900;
					BSb1 = 0.649900;
					BSc1 = 0.002300;
					BSd1 = 0.803769;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.002992;
					BSq2 = 0.269209;
					BSr2 = 0.314870;
					BSu2 = 0.214797;
					BSv2 = 0.620000;
					BSw2 = 0.727638;
					BSx2 = 0.500000;
					BSs2 = -0.049231;
					BSt2 = 0.023589;
					BSy2 = 0.178418;
					BSa2 = 0.715430;
					BSb2 = 0.796889;
					BSc2 = 0.391086;
					BSd2 = 0.293721;
					BSe2 = 2.930860;
					BSf2 = -0.057227;
					BSg2 = -0.986300;
					BSh2 = 0.000000;
					BSa3 = 0.195559;
					BSb3 = 0.310135;
					BSc3 = -0.180597;
					BSe3 = 0.009348;
					BSf3 = 0.086330;
					BSg3 = 0.213013;
					BSh3 = 0.634290;
					BSi3 = 0.000000;
					BSa4 = 1.781950;
					BSb4 = 1.029925;
					BSc4 = 0.131120;
					ICol = 3;
				}
				else if (BSFoundationCfg == 109) {
					BSa1 = -0.082700;
					BSb1 = 0.747400;
					BSc1 = 0.004940;
					BSd1 = 0.880900;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.354900;
					BSq2 = 0.765400;
					BSr2 = 0.457780;
					BSu2 = 0.193889;
					BSv2 = 0.620000;
					BSw2 = 0.376500;
					BSx2 = 0.500000;
					BSs2 = -0.002930;
					BSt2 = 0.022400;
					BSy2 = 0.143538;
					BSa2 = 0.347228;
					BSb2 = 0.872426;
					BSc2 = 0.624580;
					BSd2 = 0.123450;
					BSe2 = 1.337000;
					BSf2 = -0.098390;
					BSg2 = -0.464230;
					BSh2 = 0.000000;
					BSa3 = -0.399600;
					BSb3 = 0.515790;
					BSc3 = 0.266200;
					BSe3 = 0.253600;
					BSf3 = 0.031260;
					BSg3 = -0.013035;
					BSh3 = 0.797335;
					BSi3 = 0.000000;
					BSa4 = 2.711700;
					BSb4 = 0.191900;
					BSc4 = 0.498570;
					ICol = 5;
				}
				else if (BSFoundationCfg == 110) {
					BSa1 = -0.037530;
					BSb1 = 0.764620;
					BSc1 = 0.002210;
					BSd1 = 0.887900;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.120195;
					BSq2 = 0.770237;
					BSr2 = 0.459580;
					BSu2 = 0.194470;
					BSv2 = 0.620000;
					BSw2 = 0.373630;
					BSx2 = 0.500000;
					BSs2 = -0.003087;
					BSt2 = 0.022389;
					BSy2 = 0.144460;
					BSa2 = 0.347840;
					BSb2 = 0.879500;
					BSc2 = 0.615230;
					BSd2 = 0.113790;
					BSe2 = 1.334200;
					BSf2 = -0.100500;
					BSg2 = -0.458800;
					BSh2 = 0.000000;
					BSa3 = -0.340350;
					BSb3 = 0.516170;
					BSc3 = 0.266680;
					BSe3 = 0.251100;
					BSf3 = 0.031627;
					BSg3 = -0.005400;
					BSh3 = 0.794950;
					BSi3 = 0.000000;
					BSa4 = 2.710660;
					BSb4 = 0.195060;
					BSc4 = 0.493950;
					ICol = 5;
				}
				else if (BSFoundationCfg == 111) {
					BSa1 = -0.012480;
					BSb1 = 0.649810;
					BSc1 = 0.003444;
					BSd1 = 0.804040;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.001512;
					BSq2 = 0.422447;
					BSr2 = -0.024940;
					BSu2 = 0.681265;
					BSv2 = 0.620000;
					BSw2 = -4.687400;
					BSx2 = 0.735000;
					BSs2 = 0.529814;
					BSt2 = 0.006118;
					BSy2 = -0.096650;
					BSa2 = 3.277222;
					BSb2 = 0.192607;
					BSc2 = 0.493670;
					BSd2 = 0.265127;
					BSe2 = 0.293738;
					BSf2 = 0.049900;
					BSg2 = 0.021455;
					BSh2 = 0.000000;
					BSa3 = 0.008850;
					BSb3 = -0.044900;
					BSc3 = -0.008220;
					BSe3 = 0.151468;
					BSf3 = 0.265000;
					BSg3 = 0.134320;
					BSh3 = 0.575730;
					BSi3 = 0.000000;
					BSa4 = 2.458470;
					BSb4 = 0.483090;
					BSc4 = 1.000000;
					ICol = 6;
				}
				else if (BSFoundationCfg == 112) {
					BSa1 = -0.012880;
					BSb1 = 0.649900;
					BSc1 = 0.002346;
					BSd1 = 0.803777;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.002950;
					BSq2 = 0.103510;
					BSr2 = 0.330900;
					BSu2 = 0.214500;
					BSv2 = 0.620000;
					BSw2 = 0.722600;
					BSx2 = 0.500000;
					BSs2 = -0.062340;
					BSt2 = 0.025200;
					BSy2 = 0.177360;
					BSa2 = 0.718550;
					BSb2 = 0.761550;
					BSc2 = 0.382000;
					BSd2 = 0.278350;
					BSe2 = 2.910330;
					BSf2 = -0.051470;
					BSg2 = -0.976010;
					BSh2 = 0.000000;
					BSa3 = 0.185250;
					BSb3 = 0.298400;
					BSc3 = -0.173000;
					BSe3 = 0.012660;
					BSf3 = 0.087635;
					BSg3 = 0.212620;
					BSh3 = 0.635980;
					BSi3 = 0.000000;
					BSa4 = 1.671200;
					BSb4 = 1.140200;
					BSc4 = 0.120000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 113) {
					BSa1 = -0.013000;
					BSb1 = 0.649850;
					BSc1 = 0.003350;
					BSd1 = 0.803740;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.002180;
					BSq2 = 0.928299;
					BSr2 = -0.052500;
					BSu2 = 0.622920;
					BSv2 = 0.620000;
					BSw2 = -2.922000;
					BSx2 = 0.735000;
					BSs2 = 0.711757;
					BSt2 = 0.008103;
					BSy2 = -0.187300;
					BSa2 = 3.347700;
					BSb2 = 0.165700;
					BSc2 = 0.510060;
					BSd2 = 0.270720;
					BSe2 = 0.175190;
					BSf2 = 0.044550;
					BSg2 = 0.037400;
					BSh2 = 0.000000;
					BSa3 = 0.173880;
					BSb3 = 0.034540;
					BSc3 = -0.069000;
					BSe3 = 0.034100;
					BSf3 = 0.226800;
					BSg3 = 0.163500;
					BSh3 = 0.593769;
					BSi3 = 0.000000;
					BSa4 = 2.434300;
					BSb4 = 0.403470;
					BSc4 = 0.201620;
					ICol = 7;
				}
				else if (BSFoundationCfg == 114) {
					BSa1 = 0.092179;
					BSb1 = 0.371259;
					BSc1 = 0.048327;
					BSd1 = 0.962350;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.030159;
					BSq2 = 0.524800;
					BSr2 = -0.033659;
					BSu2 = 0.659850;
					BSv2 = 0.620000;
					BSw2 = -3.919500;
					BSx2 = 0.735000;
					BSs2 = 0.591498;
					BSt2 = 0.006708;
					BSy2 = -0.128440;
					BSa2 = 3.057845;
					BSb2 = 0.218556;
					BSc2 = 0.535304;
					BSd2 = 0.282197;
					BSe2 = 0.280364;
					BSf2 = 0.056400;
					BSg2 = 0.033100;
					BSh2 = 0.000000;
					BSa3 = 0.000966;
					BSb3 = -0.018830;
					BSc3 = -0.011740;
					BSe3 = 0.027200;
					BSf3 = 0.245900;
					BSg3 = 0.138780;
					BSh3 = 0.665040;
					BSi3 = 0.000000;
					BSa4 = 2.420990;
					BSb4 = 0.417256;
					BSc4 = 1.000000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 115) {
					BSa1 = 0.062090;
					BSb1 = 0.378913;
					BSc1 = 0.039857;
					BSd1 = 0.978853;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.060326;
					BSq2 = 7.487250;
					BSr2 = -0.452186;
					BSu2 = 0.066454;
					BSv2 = 0.620000;
					BSw2 = -2.818540;
					BSx2 = 0.735000;
					BSs2 = 1.244550;
					BSt2 = -0.000196;
					BSy2 = -0.281688;
					BSa2 = 7.284925;
					BSb2 = 0.192060;
					BSc2 = 0.810250;
					BSd2 = 0.570720;
					BSe2 = 0.056399;
					BSf2 = -0.017679;
					BSg2 = 0.018369;
					BSh2 = 0.000000;
					BSa3 = 0.041629;
					BSb3 = 0.249299;
					BSc3 = -0.096260;
					BSe3 = -0.019194;
					BSf3 = 0.093779;
					BSg3 = 0.178270;
					BSh3 = 0.794877;
					BSi3 = 0.000000;
					BSa4 = 2.374770;
					BSb4 = 0.454750;
					BSc4 = 1.000000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 116) {
					BSa1 = 1.248000;
					BSb1 = 1.583000;
					BSc1 = 0.058500;
					BSd1 = 0.271000;
					BSe1 = 3.092000;
					BSf1 = 0.292000;
					BSg1 = 1.523000;
					BSh1 = -0.918000;
					BSi1 = 1.000000;
					BSj1 = 0.142000;
					BSq2 = -0.093900;
					BSr2 = 0.080600;
					BSu2 = 1.261000;
					BSv2 = 0.620000;
					BSw2 = -0.034700;
					BSx2 = 0.500000;
					BSs2 = -0.041300;
					BSt2 = 0.014800;
					BSy2 = 0.083200;
					BSa2 = 1.789000;
					BSb2 = 0.196000;
					BSc2 = 0.561000;
					BSd2 = 0.165000;
					BSe2 = 0.030100;
					BSf2 = 0.121000;
					BSg2 = 0.144000;
					BSh2 = 0.676000;
					BSa3 = -0.148000;
					BSb3 = 0.057800;
					BSc3 = 0.127000;
					BSe3 = 0.702000;
					BSf3 = 0.136000;
					BSg3 = -0.151000;
					BSh3 = 0.183000;
					BSi3 = 1.242000;
					BSa4 = 2.720000;
					BSb4 = 0.178000;
					BSc4 = 0.499000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 117) {
					BSa1 = 1.066200;
					BSb1 = 1.743000;
					BSc1 = 0.057200;
					BSd1 = 0.277000;
					BSe1 = 3.036000;
					BSf1 = 0.233000;
					BSg1 = 1.426000;
					BSh1 = -0.873000;
					BSi1 = 1.000000;
					BSj1 = 0.124000;
					BSq2 = -0.090300;
					BSr2 = 0.079600;
					BSu2 = 1.280000;
					BSv2 = 0.620000;
					BSw2 = -0.033500;
					BSx2 = 0.500000;
					BSs2 = -0.039800;
					BSt2 = 0.014900;
					BSy2 = 0.081500;
					BSa2 = 1.792000;
					BSb2 = 0.199000;
					BSc2 = 0.560000;
					BSd2 = 0.165000;
					BSe2 = 0.032800;
					BSf2 = 0.121000;
					BSg2 = 0.141000;
					BSh2 = 0.676000;
					BSa3 = -0.155000;
					BSb3 = 0.057100;
					BSc3 = 0.132000;
					BSe3 = 0.704000;
					BSf3 = 0.136000;
					BSg3 = -0.147000;
					BSh3 = 0.179000;
					BSi3 = 1.219000;
					BSa4 = 2.732000;
					BSb4 = 0.169000;
					BSc4 = 0.520000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 118) {
					BSa1 = -0.094600;
					BSb1 = 0.766000;
					BSc1 = 0.132000;
					BSd1 = 0.992000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.085500;
					BSq2 = 1.185000;
					BSr2 = 0.199000;
					BSu2 = 0.332000;
					BSv2 = 0.620000;
					BSw2 = -0.260000;
					BSx2 = 0.500000;
					BSs2 = -0.427000;
					BSt2 = 0.014200;
					BSy2 = 0.226000;
					BSa2 = 2.210000;
					BSb2 = -0.211000;
					BSc2 = 0.566000;
					BSd2 = 0.360000;
					BSe2 = 0.653000;
					BSf2 = 0.300000;
					BSg2 = -0.248000;
					BSh2 = 0.000000;
					BSa3 = -0.727000;
					BSb3 = 0.248000;
					BSc3 = 0.749000;
					BSe3 = 0.483000;
					BSf3 = 0.177000;
					BSg3 = -0.252000;
					BSh3 = 0.728000;
					BSi3 = 0.000000;
					BSa4 = 2.827000;
					BSb4 = 0.128000;
					BSc4 = 0.411000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 119) {
					BSa1 = -0.102000;
					BSb1 = 0.776000;
					BSc1 = 0.124000;
					BSd1 = 0.964000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.074900;
					BSq2 = 1.181000;
					BSr2 = 0.198000;
					BSu2 = 0.337000;
					BSv2 = 0.620000;
					BSw2 = -0.260000;
					BSx2 = 0.500000;
					BSs2 = -0.427000;
					BSt2 = 0.014200;
					BSy2 = 0.225000;
					BSa2 = 2.215000;
					BSb2 = -0.206000;
					BSc2 = 0.565000;
					BSd2 = 0.359000;
					BSe2 = 0.655000;
					BSf2 = 0.297000;
					BSg2 = -0.248000;
					BSh2 = 0.000000;
					BSa3 = -0.729000;
					BSb3 = 0.247000;
					BSc3 = 0.754000;
					BSe3 = 0.482000;
					BSf3 = 0.177000;
					BSg3 = -0.249000;
					BSh3 = 0.728000;
					BSi3 = 0.000000;
					BSa4 = 2.828000;
					BSb4 = 0.128000;
					BSc4 = 0.411000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 120) {
					BSa1 = -0.030759;
					BSb1 = 0.651000;
					BSc1 = 0.000081;
					BSd1 = 0.809100;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075460;
					BSq2 = 0.504520;
					BSr2 = -0.032000;
					BSu2 = 0.675600;
					BSv2 = 0.620000;
					BSw2 = -3.911500;
					BSx2 = 0.730000;
					BSs2 = 0.617160;
					BSt2 = 0.007200;
					BSy2 = -0.134500;
					BSa2 = 3.076800;
					BSb2 = 0.222400;
					BSc2 = 0.490900;
					BSd2 = 0.265700;
					BSe2 = 0.280100;
					BSf2 = 0.056400;
					BSg2 = 0.022360;
					BSh2 = 0.000000;
					BSa3 = 0.034480;
					BSb3 = -0.032400;
					BSc3 = -0.018500;
					BSe3 = 0.092230;
					BSf3 = 0.244000;
					BSg3 = 0.168100;
					BSh3 = 0.574200;
					BSi3 = 0.000000;
					BSa4 = 2.502000;
					BSb4 = 0.432500;
					BSc4 = 1.000000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 121) {
					BSa1 = -0.035830;
					BSb1 = 0.756150;
					BSc1 = 0.094100;
					BSd1 = 0.831300;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.025030;
					BSq2 = 6.654780;
					BSr2 = -0.423500;
					BSu2 = 0.072140;
					BSv2 = 0.620000;
					BSw2 = -2.802500;
					BSx2 = 0.735000;
					BSs2 = 1.174800;
					BSt2 = 0.004120;
					BSy2 = -0.278130;
					BSa2 = 6.769600;
					BSb2 = 0.193100;
					BSc2 = 0.758400;
					BSd2 = 0.527400;
					BSe2 = 0.100900;
					BSf2 = -0.026750;
					BSg2 = 0.020180;
					BSh2 = 0.000000;
					BSa3 = 0.071520;
					BSb3 = 0.207100;
					BSc3 = -0.111500;
					BSe3 = 0.102980;
					BSf3 = 0.131900;
					BSg3 = 0.212360;
					BSh3 = 0.769790;
					BSi3 = 0.000000;
					BSa4 = 2.471400;
					BSb4 = 0.540800;
					BSc4 = 1.000000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 122) {
					BSa1 = -0.102755;
					BSb1 = 0.758400;
					BSc1 = 0.127850;
					BSd1 = 0.923830;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.095400;
					BSq2 = 0.356082;
					BSr2 = 1.175170;
					BSu2 = 0.193035;
					BSv2 = 0.620000;
					BSw2 = 0.019890;
					BSx2 = 0.500000;
					BSs2 = 0.385205;
					BSt2 = 0.010492;
					BSy2 = 0.231285;
					BSa2 = 0.533470;
					BSb2 = 0.529600;
					BSc2 = 0.678140;
					BSd2 = -0.226610;
					BSe2 = 0.194250;
					BSf2 = -0.008800;
					BSg2 = -0.063400;
					BSh2 = 0.000000;
					BSa3 = -0.505580;
					BSb3 = 0.463320;
					BSc3 = 0.516750;
					BSe3 = 0.318440;
					BSf3 = 0.071850;
					BSg3 = -0.106340;
					BSh3 = 0.820841;
					BSi3 = 0.000000;
					BSa4 = 2.785130;
					BSb4 = 0.177090;
					BSc4 = 1.000000;
					ICol = 2;
				}
				else if (BSFoundationCfg == 123) {
					BSa1 = -0.110115;
					BSb1 = 0.768570;
					BSc1 = 0.119650;
					BSd1 = 0.896200;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.084336;
					BSq2 = 0.322508;
					BSr2 = 1.177253;
					BSu2 = 0.192520;
					BSv2 = 0.620000;
					BSw2 = 0.018815;
					BSx2 = 0.500000;
					BSs2 = 0.383320;
					BSt2 = 0.010505;
					BSy2 = 0.232525;
					BSa2 = 0.546050;
					BSb2 = 0.534980;
					BSc2 = 0.673230;
					BSd2 = -0.220550;
					BSe2 = 0.193790;
					BSf2 = -0.008839;
					BSg2 = -0.062920;
					BSh2 = 0.000000;
					BSa3 = -0.507432;
					BSb3 = 0.462610;
					BSc3 = 0.521787;
					BSe3 = 0.317700;
					BSf3 = 0.072040;
					BSg3 = -0.103580;
					BSh3 = 0.820152;
					BSi3 = 0.000000;
					BSa4 = 2.785986;
					BSb4 = 0.177500;
					BSc4 = 1.000000;
					ICol = 2;
				}
				else if (BSFoundationCfg == 124) {
					BSa1 = -0.095730;
					BSb1 = 0.759500;
					BSc1 = 0.120900;
					BSd1 = 0.944400;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.094136;
					BSq2 = 0.392119;
					BSr2 = 1.244337;
					BSu2 = 0.176210;
					BSv2 = 0.620000;
					BSw2 = -0.017980;
					BSx2 = 0.500000;
					BSs2 = 0.334740;
					BSt2 = 0.010732;
					BSy2 = 0.258040;
					BSa2 = 0.525216;
					BSb2 = 0.588802;
					BSc2 = 0.676886;
					BSd2 = -0.218160;
					BSe2 = 0.210820;
					BSf2 = -0.012470;
					BSg2 = -0.068500;
					BSh2 = 0.000000;
					BSa3 = -0.581040;
					BSb3 = 0.467500;
					BSc3 = 0.554980;
					BSe3 = 0.345205;
					BSf3 = 0.063890;
					BSg3 = -0.116690;
					BSh3 = 0.834930;
					BSi3 = 0.000000;
					BSa4 = 2.756980;
					BSb4 = 0.214150;
					BSc4 = 1.000000;
					ICol = 2;
				}
				else if (BSFoundationCfg == 125) {
					BSa1 = -0.103200;
					BSb1 = 0.769500;
					BSc1 = 0.112700;
					BSd1 = 0.915800;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.083070;
					BSq2 = 0.347370;
					BSr2 = 1.247860;
					BSu2 = 0.175410;
					BSv2 = 0.620000;
					BSw2 = -0.019260;
					BSx2 = 0.500000;
					BSs2 = 0.332680;
					BSt2 = 0.010750;
					BSy2 = 0.259590;
					BSa2 = 0.539770;
					BSb2 = 0.593620;
					BSc2 = 0.672010;
					BSd2 = -0.211000;
					BSe2 = 0.209830;
					BSf2 = -0.012460;
					BSg2 = -0.067828;
					BSh2 = 0.000000;
					BSa3 = -0.582970;
					BSb3 = 0.466810;
					BSc3 = 0.560070;
					BSe3 = 0.344520;
					BSf3 = 0.064070;
					BSg3 = -0.113950;
					BSh3 = 0.834110;
					BSi3 = 0.000000;
					BSa4 = 2.757860;
					BSb4 = 0.214560;
					BSc4 = 1.000000;
					ICol = 2;
				}
				else if (BSFoundationCfg == 126) {
					BSa1 = -0.082866;
					BSb1 = 0.747000;
					BSc1 = 0.005117;
					BSd1 = 0.880190;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.354580;
					BSq2 = 0.473306;
					BSr2 = 0.512520;
					BSu2 = 0.208830;
					BSv2 = 0.620000;
					BSw2 = 0.311590;
					BSx2 = 0.500000;
					BSs2 = 0.016306;
					BSt2 = 0.021280;
					BSy2 = 0.172150;
					BSa2 = 0.346200;
					BSb2 = 0.930900;
					BSc2 = 0.668600;
					BSd2 = -0.013540;
					BSe2 = 0.792270;
					BSf2 = -0.067900;
					BSg2 = -0.262030;
					BSh2 = 0.000000;
					BSa3 = -0.347432;
					BSb3 = 0.471450;
					BSc3 = 0.299225;
					BSe3 = 0.252390;
					BSf3 = 0.037600;
					BSg3 = -0.006287;
					BSh3 = 0.802489;
					BSi3 = 0.000000;
					BSa4 = 2.669500;
					BSb4 = 0.230600;
					BSc4 = 0.500000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 127) {
					BSa1 = -0.037710;
					BSb1 = 0.764220;
					BSc1 = 0.002345;
					BSd1 = 0.887697;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.120000;
					BSq2 = 0.471953;
					BSr2 = 0.512454;
					BSu2 = 0.208830;
					BSv2 = 0.620000;
					BSw2 = 0.311880;
					BSx2 = 0.500000;
					BSs2 = 0.016550;
					BSt2 = 0.021280;
					BSy2 = 0.172080;
					BSa2 = 0.346870;
					BSb2 = 0.929640;
					BSc2 = 0.668740;
					BSd2 = -0.013460;
					BSe2 = 0.790300;
					BSf2 = -0.067800;
					BSg2 = -0.261300;
					BSh2 = 0.000000;
					BSa3 = -0.374325;
					BSb3 = 0.471450;
					BSc3 = 0.299225;
					BSe3 = 0.252390;
					BSf3 = 0.037000;
					BSg3 = -0.006287;
					BSh3 = 0.802500;
					BSi3 = 0.000000;
					BSa4 = 2.669500;
					BSb4 = 0.230570;
					BSc4 = 0.500000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 128) {
					BSa1 = -0.083846;
					BSb1 = 0.748270;
					BSc1 = 0.003780;
					BSd1 = 0.884560;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.356000;
					BSq2 = -0.013900;
					BSr2 = 0.346370;
					BSu2 = 0.224479;
					BSv2 = 0.620000;
					BSw2 = 0.658920;
					BSx2 = 0.500000;
					BSs2 = -0.028296;
					BSt2 = 0.022880;
					BSy2 = 0.179770;
					BSa2 = 0.780258;
					BSb2 = 0.620340;
					BSc2 = 0.384045;
					BSd2 = 0.189720;
					BSe2 = 2.519800;
					BSf2 = 0.003260;
					BSg2 = -0.815290;
					BSh2 = 0.000000;
					BSa3 = 0.119523;
					BSb3 = 0.253770;
					BSc3 = -0.147840;
					BSe3 = 0.035197;
					BSf3 = 0.103024;
					BSg3 = 0.255490;
					BSh3 = 0.737528;
					BSi3 = 0.000000;
					BSa4 = 2.523800;
					BSb4 = 0.449460;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 129) {
					BSa1 = -0.037960;
					BSb1 = 0.764560;
					BSc1 = 0.002095;
					BSd1 = 0.882560;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.120206;
					BSq2 = -0.012624;
					BSr2 = 0.346870;
					BSu2 = 0.225017;
					BSv2 = 0.620000;
					BSw2 = 0.656788;
					BSx2 = 0.500000;
					BSs2 = -0.028665;
					BSt2 = 0.022879;
					BSy2 = 0.179115;
					BSa2 = 0.779820;
					BSb2 = 0.631300;
					BSc2 = 0.380318;
					BSd2 = 0.184787;
					BSe2 = 2.496030;
					BSf2 = 0.003646;
					BSg2 = -0.803790;
					BSh2 = 0.000000;
					BSa3 = 0.119500;
					BSb3 = 0.253770;
					BSc3 = -0.147800;
					BSe3 = 0.035200;
					BSf3 = 0.103020;
					BSg3 = 0.255500;
					BSh3 = 0.737530;
					BSi3 = 0.000000;
					BSa4 = 2.523800;
					BSb4 = 0.449500;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 130) {
					BSa1 = -0.082930;
					BSb1 = 0.747060;
					BSc1 = 0.005237;
					BSd1 = 0.880240;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.354600;
					BSq2 = 6.365400;
					BSr2 = -0.440800;
					BSu2 = 0.062176;
					BSv2 = 0.620000;
					BSw2 = -3.020130;
					BSx2 = 0.730000;
					BSs2 = 1.285500;
					BSt2 = 0.000656;
					BSy2 = -0.275300;
					BSa2 = 6.466700;
					BSb2 = 0.308730;
					BSc2 = 0.795500;
					BSd2 = 0.488500;
					BSe2 = 0.076124;
					BSf2 = -0.013000;
					BSg2 = -0.017600;
					BSh2 = 0.000000;
					BSa3 = -0.358950;
					BSb3 = 0.470600;
					BSc3 = 0.291980;
					BSe3 = 0.251810;
					BSf3 = 0.041790;
					BSg3 = -0.014440;
					BSh3 = 0.803190;
					BSi3 = 0.000000;
					BSa4 = 2.735900;
					BSb4 = 0.186500;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 131) {
					BSa1 = -0.037760;
					BSb1 = 0.764260;
					BSc1 = 0.002430;
					BSd1 = 0.887720;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.120050;
					BSq2 = 6.338190;
					BSr2 = -0.440300;
					BSu2 = 0.062640;
					BSv2 = 0.620000;
					BSw2 = -3.029940;
					BSx2 = 0.735000;
					BSs2 = 1.287370;
					BSt2 = 0.000713;
					BSy2 = -0.276570;
					BSa2 = 6.454200;
					BSb2 = 0.310885;
					BSc2 = 0.793500;
					BSd2 = 0.487150;
					BSe2 = 0.076922;
					BSf2 = -0.013430;
					BSg2 = -0.017130;
					BSh2 = 0.000000;
					BSa3 = -0.359600;
					BSb3 = 0.470970;
					BSc3 = 0.292400;
					BSe3 = 0.249220;
					BSf3 = 0.042160;
					BSg3 = -0.006776;
					BSh3 = 0.800740;
					BSi3 = 0.000000;
					BSa4 = 2.669680;
					BSb4 = 0.228820;
					BSc4 = 0.500000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 132) {
					BSa1 = -0.013317;
					BSb1 = 0.649880;
					BSc1 = 0.003184;
					BSd1 = 0.804315;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.002130;
					BSq2 = -0.097633;
					BSr2 = 0.099027;
					BSu2 = 0.550320;
					BSv2 = 0.620000;
					BSw2 = 1.116420;
					BSx2 = 0.500000;
					BSs2 = -0.048976;
					BSt2 = 0.019566;
					BSy2 = 0.130299;
					BSa2 = 0.599600;
					BSb2 = 0.717490;
					BSc2 = 0.340250;
					BSd2 = 0.081380;
					BSe2 = 0.893277;
					BSf2 = 0.113290;
					BSg2 = -0.169400;
					BSh2 = 0.000000;
					BSa3 = 0.072319;
					BSb3 = 0.084716;
					BSc3 = -0.070685;
					BSe3 = 0.023877;
					BSf3 = 0.097430;
					BSg3 = 0.241229;
					BSh3 = 0.649840;
					BSi3 = 0.000000;
					BSa4 = 2.312320;
					BSb4 = 0.697238;
					BSc4 = 1.000000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 133) {
					BSa1 = -0.012950;
					BSb1 = 0.649780;
					BSc1 = 0.003517;
					BSd1 = 0.804050;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.001567;
					BSq2 = -1.186390;
					BSr2 = 0.114417;
					BSu2 = 0.921700;
					BSv2 = 0.620000;
					BSw2 = 0.000000;
					BSx2 = 0.500000;
					BSs2 = -0.067360;
					BSt2 = 0.005783;
					BSy2 = 0.355200;
					BSa2 = 1.866000;
					BSb2 = 0.309700;
					BSc2 = 0.398500;
					BSd2 = 0.165950;
					BSe2 = 0.333000;
					BSf2 = 0.062570;
					BSg2 = 0.028100;
					BSh2 = 0.000000;
					BSa3 = 0.041710;
					BSb3 = -0.022650;
					BSc3 = -0.032800;
					BSe3 = 0.057337;
					BSf3 = 0.184820;
					BSg3 = 0.207860;
					BSh3 = 0.552620;
					BSi3 = 0.000000;
					BSa4 = 2.500000;
					BSb4 = 0.468100;
					BSc4 = 1.000000;
					ICol = 7;
				}
				else if (BSFoundationCfg == 134) {
					BSa1 = -0.013181;
					BSb1 = 0.649820;
					BSc1 = 0.003152;
					BSd1 = 0.804095;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.002091;
					BSq2 = -0.058239;
					BSr2 = 0.106230;
					BSu2 = 0.544950;
					BSv2 = 0.620000;
					BSw2 = 1.079680;
					BSx2 = 0.500000;
					BSs2 = -0.043880;
					BSt2 = 0.019430;
					BSy2 = 0.131390;
					BSa2 = 0.521400;
					BSb2 = 0.769000;
					BSc2 = 0.335500;
					BSd2 = 0.065850;
					BSe2 = 1.007430;
					BSf2 = 0.110560;
					BSg2 = -0.213200;
					BSh2 = 0.000000;
					BSa3 = 0.089900;
					BSb3 = 0.103910;
					BSc3 = -0.080670;
					BSe3 = 0.011959;
					BSf3 = 0.088945;
					BSg3 = 0.243870;
					BSh3 = 0.658276;
					BSi3 = 0.000000;
					BSa4 = 2.329980;
					BSb4 = 0.659540;
					BSc4 = 1.000000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 135) {
					BSa1 = -0.014170;
					BSb1 = 0.650050;
					BSc1 = 0.002820;
					BSd1 = 0.808200;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.004060;
					BSq2 = -0.892810;
					BSr2 = 0.365590;
					BSu2 = 0.659620;
					BSv2 = 0.620000;
					BSw2 = 0.043142;
					BSx2 = 0.500000;
					BSs2 = 0.074252;
					BSt2 = 0.011677;
					BSy2 = 0.303140;
					BSa2 = 0.870480;
					BSb2 = 0.653690;
					BSc2 = 0.556110;
					BSd2 = 0.019229;
					BSe2 = 0.220000;
					BSf2 = -0.008110;
					BSg2 = -0.067740;
					BSh2 = 0.000000;
					BSa3 = -0.350977;
					BSb3 = 0.226320;
					BSc3 = 0.372970;
					BSe3 = 0.249970;
					BSf3 = 0.055950;
					BSg3 = -0.036000;
					BSh3 = 0.678300;
					BSi3 = 0.000000;
					BSa4 = 1.597820;
					BSb4 = 1.252800;
					BSc4 = 0.105500;
					ICol = 2;
				}
				else if (BSFoundationCfg == 136) {
					BSa1 = -0.030890;
					BSb1 = 0.650530;
					BSc1 = -0.000019;
					BSd1 = 0.809187;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075536;
					BSq2 = -0.100000;
					BSr2 = 0.100700;
					BSu2 = 0.538523;
					BSv2 = 0.620000;
					BSw2 = 1.126370;
					BSx2 = 0.500000;
					BSs2 = -0.039893;
					BSt2 = 0.019840;
					BSy2 = 0.127870;
					BSa2 = 0.615260;
					BSb2 = 0.694830;
					BSc2 = 0.333130;
					BSd2 = 0.079298;
					BSe2 = 0.882230;
					BSf2 = 0.111220;
					BSg2 = -0.167070;
					BSh2 = 0.000000;
					BSa3 = 0.063040;
					BSb3 = 0.081414;
					BSc3 = -0.072290;
					BSe3 = 0.030850;
					BSf3 = 0.099640;
					BSg3 = 0.246660;
					BSh3 = 0.625840;
					BSi3 = 0.000000;
					BSa4 = 2.307880;
					BSb4 = 0.715339;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 137) {
					BSa1 = -0.030900;
					BSb1 = 0.650520;
					BSc1 = 0.000126;
					BSd1 = 0.809110;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075400;
					BSq2 = -0.100666;
					BSr2 = 0.006160;
					BSu2 = 0.811340;
					BSv2 = 0.620000;
					BSw2 = 9.354000;
					BSx2 = 0.500000;
					BSs2 = 0.359300;
					BSt2 = 0.006200;
					BSy2 = 0.004830;
					BSa2 = 2.451100;
					BSb2 = 0.257300;
					BSc2 = 0.401400;
					BSd2 = 0.223140;
					BSe2 = 0.276600;
					BSf2 = 0.034240;
					BSg2 = 0.026860;
					BSh2 = 0.000000;
					BSa3 = 0.032563;
					BSb3 = -0.023900;
					BSc3 = -0.033700;
					BSe3 = 0.068500;
					BSf3 = 0.185210;
					BSg3 = 0.211670;
					BSh3 = 0.541690;
					BSi3 = 0.000000;
					BSa4 = 2.491800;
					BSb4 = 0.493390;
					BSc4 = 1.000000;
					ICol = 8;
				}
				else if (BSFoundationCfg == 138) {
					BSa1 = -0.030800;
					BSb1 = 0.650500;
					BSc1 = -0.000310;
					BSd1 = 0.810000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075960;
					BSq2 = -0.564500;
					BSr2 = 0.250510;
					BSu2 = 0.611998;
					BSv2 = 0.620000;
					BSw2 = 0.277700;
					BSx2 = 0.500000;
					BSs2 = 0.102009;
					BSt2 = 0.013950;
					BSy2 = 0.218030;
					BSa2 = 0.838520;
					BSb2 = 0.643200;
					BSc2 = 0.541440;
					BSd2 = 0.041040;
					BSe2 = 0.261800;
					BSf2 = -0.015200;
					BSg2 = -0.078700;
					BSh2 = 0.000000;
					BSa3 = -0.283600;
					BSb3 = 0.215440;
					BSc3 = 0.288500;
					BSe3 = 0.214720;
					BSf3 = 0.049360;
					BSg3 = 0.009120;
					BSh3 = 0.678600;
					BSi3 = 0.000000;
					BSa4 = 2.548800;
					BSb4 = 0.367100;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 139) {
					BSa1 = -0.030880;
					BSb1 = 0.650540;
					BSc1 = -0.000045;
					BSd1 = 0.809210;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075558;
					BSq2 = -0.057123;
					BSr2 = 0.107770;
					BSu2 = 0.540350;
					BSv2 = 0.620000;
					BSw2 = 1.068860;
					BSx2 = 0.500000;
					BSs2 = -0.040303;
					BSt2 = 0.019384;
					BSy2 = 0.133140;
					BSa2 = 0.539700;
					BSb2 = 0.756590;
					BSc2 = 0.331970;
					BSd2 = 0.069257;
					BSe2 = 1.018200;
					BSf2 = 0.103580;
					BSg2 = -0.220000;
					BSh2 = 0.000000;
					BSa3 = 0.080160;
					BSb3 = 0.100000;
					BSc3 = -0.082570;
					BSe3 = 0.018670;
					BSf3 = 0.091450;
					BSg3 = 0.249800;
					BSh3 = 0.631300;
					BSi3 = 0.000000;
					BSa4 = 2.326200;
					BSb4 = 0.676240;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 140) {
					BSa1 = -0.012590;
					BSb1 = 0.649880;
					BSc1 = 0.002510;
					BSd1 = 0.804100;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.002360;
					BSq2 = -0.259990;
					BSr2 = 0.362950;
					BSu2 = 0.153770;
					BSv2 = 0.620000;
					BSw2 = 0.665620;
					BSx2 = 0.500000;
					BSs2 = -0.032500;
					BSt2 = 0.025330;
					BSy2 = 0.173580;
					BSa2 = 0.771200;
					BSb2 = 0.567300;
					BSc2 = 0.313740;
					BSd2 = 0.154600;
					BSe2 = 2.042900;
					BSf2 = 0.002390;
					BSg2 = -0.615050;
					BSh2 = 0.000000;
					BSa3 = 0.035480;
					BSb3 = 0.213550;
					BSc3 = -0.105500;
					BSe3 = 0.085500;
					BSf3 = 0.099230;
					BSg3 = 0.205840;
					BSh3 = 0.674560;
					BSi3 = 0.000000;
					BSa4 = 2.303130;
					BSb4 = 0.647160;
					BSc4 = 1.000000;
					ICol = 3;
				}
				else if (BSFoundationCfg == 141) {
					BSa1 = -0.013690;
					BSb1 = 0.650050;
					BSc1 = 0.002022;
					BSd1 = 0.807080;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.004660;
					BSq2 = 7.336860;
					BSr2 = -0.471580;
					BSu2 = 0.096200;
					BSv2 = 0.620000;
					BSw2 = -2.740600;
					BSx2 = 0.730000;
					BSs2 = 1.311300;
					BSt2 = 0.006077;
					BSy2 = -0.309480;
					BSa2 = 6.739200;
					BSb2 = 0.261210;
					BSc2 = 0.774910;
					BSd2 = 0.501790;
					BSe2 = 0.060400;
					BSf2 = -0.009997;
					BSg2 = -0.013330;
					BSh2 = 0.000000;
					BSa3 = -0.179200;
					BSb3 = 0.419230;
					BSc3 = 0.201200;
					BSe3 = 0.201400;
					BSf3 = 0.046400;
					BSg3 = -0.010700;
					BSh3 = 0.663900;
					BSi3 = 0.000000;
					BSa4 = 2.712500;
					BSb4 = 0.160500;
					BSc4 = 1.000000;
					ICol = 2;
				}
				else if (BSFoundationCfg == 142) {
					BSa1 = -0.030640;
					BSb1 = 0.650530;
					BSc1 = -0.000307;
					BSd1 = 0.809200;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075730;
					BSq2 = 7.627630;
					BSr2 = -0.467320;
					BSu2 = 0.078184;
					BSv2 = 0.620000;
					BSw2 = -2.630500;
					BSx2 = 0.730000;
					BSs2 = 1.249332;
					BSt2 = 0.007903;
					BSy2 = -0.306510;
					BSa2 = 6.940040;
					BSb2 = 0.191590;
					BSc2 = 0.758668;
					BSd2 = 0.553020;
					BSe2 = 0.091920;
					BSf2 = -0.028970;
					BSg2 = 0.020170;
					BSh2 = 0.000000;
					BSa3 = 0.051800;
					BSb3 = 0.212390;
					BSc3 = -0.111600;
					BSe3 = 0.073330;
					BSf3 = 0.099610;
					BSg3 = 0.214240;
					BSh3 = 0.648060;
					BSi3 = 0.000000;
					BSa4 = 2.347400;
					BSb4 = 0.595840;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 143) {
					BSa1 = -0.030585;
					BSb1 = 0.650520;
					BSc1 = -0.000376;
					BSd1 = 0.809200;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075760;
					BSq2 = 8.625330;
					BSr2 = -0.498760;
					BSu2 = 0.081750;
					BSv2 = 0.620000;
					BSw2 = -2.553000;
					BSx2 = 0.730000;
					BSs2 = 1.308400;
					BSt2 = 0.007158;
					BSy2 = -0.327700;
					BSa2 = 7.506358;
					BSb2 = 0.188600;
					BSc2 = 0.775200;
					BSd2 = 0.569300;
					BSe2 = 0.077900;
					BSf2 = -0.026120;
					BSg2 = 0.019240;
					BSh2 = 0.000000;
					BSa3 = 0.070580;
					BSb3 = 0.253220;
					BSc3 = -0.132400;
					BSe3 = 0.062050;
					BSf3 = 0.092810;
					BSg3 = 0.216200;
					BSh3 = 0.636400;
					BSi3 = 0.000000;
					BSa4 = 2.394550;
					BSb4 = 0.526970;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 144) {
					BSa1 = -0.030440;
					BSb1 = 0.650000;
					BSc1 = -0.000700;
					BSd1 = 0.810000;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.076200;
					BSq2 = 6.980800;
					BSr2 = -0.462000;
					BSu2 = 0.089760;
					BSv2 = 0.620000;
					BSw2 = -2.833000;
					BSx2 = 0.735000;
					BSs2 = 1.316000;
					BSt2 = 0.007800;
					BSy2 = -0.306200;
					BSa2 = 6.558000;
					BSb2 = 0.285200;
					BSc2 = 0.780500;
					BSd2 = 0.514480;
					BSe2 = 0.076770;
					BSf2 = -0.014130;
					BSg2 = -0.015600;
					BSh2 = 0.000000;
					BSa3 = -0.253500;
					BSb3 = 0.389750;
					BSc3 = 0.206770;
					BSe3 = 0.228810;
					BSf3 = 0.046850;
					BSg3 = -0.001900;
					BSh3 = 0.665640;
					BSi3 = 0.000000;
					BSa4 = 2.622900;
					BSb4 = 0.240800;
					BSc4 = 1.000000;
					ICol = 5;
				}
				else if (BSFoundationCfg == 145) {
					BSa1 = -0.030756;
					BSb1 = 0.650500;
					BSc1 = 0.000055;
					BSd1 = 0.809140;
					BSe1 = 1.000000;
					BSf1 = 1.000000;
					BSg1 = 1.000000;
					BSh1 = 1.000000;
					BSi1 = 0.000000;
					BSj1 = 0.075449;
					BSq2 = 0.588859;
					BSr2 = -0.037886;
					BSu2 = 0.644873;
					BSv2 = 0.620000;
					BSw2 = -3.702740;
					BSx2 = 0.730000;
					BSs2 = 0.653646;
					BSt2 = 0.007464;
					BSy2 = -0.145560;
					BSa2 = 3.113500;
					BSb2 = 0.219430;
					BSc2 = 0.507580;
					BSd2 = 0.272440;
					BSe2 = 0.252230;
					BSf2 = 0.052930;
					BSg2 = 0.020660;
					BSh2 = 0.000000;
					BSa3 = 0.035720;
					BSb3 = -0.013200;
					BSc3 = -0.018100;
					BSe3 = 0.087370;
					BSf3 = 0.233800;
					BSg3 = 0.165200;
					BSh3 = 0.592780;
					BSi3 = 0.000000;
					BSa4 = 2.454450;
					BSb4 = 0.482300;
					BSc4 = 1.000000;
					ICol = 7;
				}
				else {
					// Configuration not in list -> Error
					BSFndConfigFound = false;
				}
			}  // else loop for config >= 100

			
			// Declare error if configuration not found
			if (!BSFndConfigFound){
				// Declare variable to hold the routine name
				static std::string const RoutineName("GetBSCoeff: ");

				// Identify the EnergyPlus object being read
				cCurrentModuleObject = "Site:GroundDomain:BASESIMP";

				// Write error message
				ShowSevereError(RoutineName + cCurrentModuleObject + "=\"" + BaseSimp(BSFoundationNum).BSFoundationName + "\",");
				ShowContinueError("Invalid Configuration" + '=' + BSFoundationCfg);
				// Process error
				ShowFatalError(RoutineName + "Errors found in processing " + cCurrentModuleObject + " input.");
			}

		} // routine GetBSCoeff


		void
		InitBSCornerCoeff()
		{
			// SUBROUTINE INFORMATION:
			//       AUTHOR         Patrice Pinel
			//       DATE WRITTEN   March. 2016
			//       MODIFIED       
			//       RE-ENGINEERED  na

			// PURPOSE OF THIS SUBROUTINE:
			// This subrountine assigns the BASESIMP corner coefficients values in their table

			// METHODOLOGY EMPLOYED: 
			// This routine reproduces the BSCORNER ESP-r routine from the bscoeff.F source file:

			// REFERENCES: ESP-r source code
			
			// Note, array contains one element in each indices more than needed
			// Indices are initiated from 1 rather than 0 in order to avoid errors in conversion for Fortran
			BSCornerCoeff[1][1] = 0.78404336;
			BSCornerCoeff[1][2] = 0.00000000;
			BSCornerCoeff[1][3] = -0.06269160;
			BSCornerCoeff[1][4] = 0.01588997;
			BSCornerCoeff[1][5] = 0.18362466;
			BSCornerCoeff[1][6] = -0.00280010;
			BSCornerCoeff[1][7] = 0.00000000;
			BSCornerCoeff[1][8] = 0.00000000;
			BSCornerCoeff[1][9] = 0.00000000;
			BSCornerCoeff[1][10] = -0.00058830;
			BSCornerCoeff[1][11] = -0.00323320;
			BSCornerCoeff[1][12] = 0.00000000;
			BSCornerCoeff[1][13] = 0.01234568;
			BSCornerCoeff[1][14] = -0.01884040;
			BSCornerCoeff[1][15] = -0.04484290;
			BSCornerCoeff[1][16] = 0.00000000;
			BSCornerCoeff[1][17] = 0.00115767;
			BSCornerCoeff[1][18] = 0.00186534;
			BSCornerCoeff[1][19] = 0.00302408;
			BSCornerCoeff[2][1] = 0.75572785;
			BSCornerCoeff[2][2] = 0.00000000;
			BSCornerCoeff[2][3] = -0.04127260;
			BSCornerCoeff[2][4] = 0.06788635;
			BSCornerCoeff[2][5] = 0.17354284;
			BSCornerCoeff[2][6] = -0.00142870;
			BSCornerCoeff[2][7] = 0.00000000;
			BSCornerCoeff[2][8] = 0.00000000;
			BSCornerCoeff[2][9] = 0.00000000;
			BSCornerCoeff[2][10] = 0.00369666;
			BSCornerCoeff[2][11] = -0.00508140;
			BSCornerCoeff[2][12] = 0.00000000;
			BSCornerCoeff[2][13] = 0.01781545;
			BSCornerCoeff[2][14] = -0.01948280;
			BSCornerCoeff[2][15] = -0.03993470;
			BSCornerCoeff[2][16] = 0.00000000;
			BSCornerCoeff[2][17] = -0.00069420;
			BSCornerCoeff[2][18] = -0.00023500;
			BSCornerCoeff[2][19] = 0.00259176;
			BSCornerCoeff[3][1] = 0.89231276;
			BSCornerCoeff[3][2] = -0.00590050;
			BSCornerCoeff[3][3] = -0.07535490;
			BSCornerCoeff[3][4] = -0.02696010;
			BSCornerCoeff[3][5] = 0.06201586;
			BSCornerCoeff[3][6] = -0.00249220;
			BSCornerCoeff[3][7] = 0.00039837;
			BSCornerCoeff[3][8] = -0.00062370;
			BSCornerCoeff[3][9] = -0.00089240;
			BSCornerCoeff[3][10] = 0.00308021;
			BSCornerCoeff[3][11] = 0.00185631;
			BSCornerCoeff[3][12] = 0.00237192;
			BSCornerCoeff[3][13] = 0.01246938;
			BSCornerCoeff[3][14] = -0.02075720;
			BSCornerCoeff[3][15] = 0.00000000;
			BSCornerCoeff[3][16] = 0.00023050;
			BSCornerCoeff[3][17] = 0.00063693;
			BSCornerCoeff[3][18] = 0.00189477;
			BSCornerCoeff[3][19] = 0.00490457;
			BSCornerCoeff[4][1] = 0.85653204;
			BSCornerCoeff[4][2] = 0.00560896;
			BSCornerCoeff[4][3] = -0.04706910;
			BSCornerCoeff[4][4] = 0.03933353;
			BSCornerCoeff[4][5] = 0.07332308;
			BSCornerCoeff[4][6] = -0.00260210;
			BSCornerCoeff[4][7] = -0.00049210;
			BSCornerCoeff[4][8] = -0.00142660;
			BSCornerCoeff[4][9] = -0.00063020;
			BSCornerCoeff[4][10] = 0.00551828;
			BSCornerCoeff[4][11] = -0.00129990;
			BSCornerCoeff[4][12] = 0.00224003;
			BSCornerCoeff[4][13] = 0.01512380;
			BSCornerCoeff[4][14] = -0.01987420;
			BSCornerCoeff[4][15] = 0.00000000;
			BSCornerCoeff[4][16] = -0.00002580;
			BSCornerCoeff[4][17] = -0.00098480;
			BSCornerCoeff[4][18] = -0.00034900;
			BSCornerCoeff[4][19] = 0.00407644;
			BSCornerCoeff[5][1] = 0.71725273;
			BSCornerCoeff[5][2] = -0.02203810;
			BSCornerCoeff[5][3] = -0.06561350;
			BSCornerCoeff[5][4] = 0.03467739;
			BSCornerCoeff[5][5] = 0.11542548;
			BSCornerCoeff[5][6] = -0.00243060;
			BSCornerCoeff[5][7] = 0.00239333;
			BSCornerCoeff[5][8] = 0.00017063;
			BSCornerCoeff[5][9] = 0.00078584;
			BSCornerCoeff[5][10] = 0.00538004;
			BSCornerCoeff[5][11] = -0.00690080;
			BSCornerCoeff[5][12] = -0.00684610;
			BSCornerCoeff[5][13] = -0.00592680;
			BSCornerCoeff[5][14] = -0.01821230;
			BSCornerCoeff[5][15] = -0.01805300;
			BSCornerCoeff[5][16] = 0.00021603;
			BSCornerCoeff[5][17] = 0.00007409;
			BSCornerCoeff[5][18] = 0.00219304;
			BSCornerCoeff[5][19] = 0.00495485;
			BSCornerCoeff[6][1] = 0.66286115;
			BSCornerCoeff[6][2] = -0.01117000;
			BSCornerCoeff[6][3] = -0.04668110;
			BSCornerCoeff[6][4] = 0.10415677;
			BSCornerCoeff[6][5] = 0.15558250;
			BSCornerCoeff[6][6] = -0.00099990;
			BSCornerCoeff[6][7] = 0.00123304;
			BSCornerCoeff[6][8] = -0.00138840;
			BSCornerCoeff[6][9] = 0.00144637;
			BSCornerCoeff[6][10] = 0.00813609;
			BSCornerCoeff[6][11] = -0.01029190;
			BSCornerCoeff[6][12] = -0.00276280;
			BSCornerCoeff[6][13] = 0.00179876;
			BSCornerCoeff[6][14] = -0.01555330;
			BSCornerCoeff[6][15] = -0.02351600;
			BSCornerCoeff[6][16] = -0.00013170;
			BSCornerCoeff[6][17] = -0.00130010;
			BSCornerCoeff[6][18] = -0.00017010;
			BSCornerCoeff[6][19] = 0.00211170;
			BSCornerCoeff[7][1] = 0.79399414;
			BSCornerCoeff[7][2] = -0.00590040;
			BSCornerCoeff[7][3] = -0.06301890;
			BSCornerCoeff[7][4] = 0.01801793;
			BSCornerCoeff[7][5] = 0.16194821;
			BSCornerCoeff[7][6] = -0.00370510;
			BSCornerCoeff[7][7] = 0.00081849;
			BSCornerCoeff[7][8] = -0.00032200;
			BSCornerCoeff[7][9] = 0.00039037;
			BSCornerCoeff[7][10] = -0.00088960;
			BSCornerCoeff[7][11] = -0.00412530;
			BSCornerCoeff[7][12] = -0.00348100;
			BSCornerCoeff[7][13] = 0.00266370;
			BSCornerCoeff[7][14] = -0.01568740;
			BSCornerCoeff[7][15] = -0.04236650;
			BSCornerCoeff[7][16] = 0.00004032;
			BSCornerCoeff[7][17] = 0.00137467;
			BSCornerCoeff[7][18] = 0.00205882;
			BSCornerCoeff[7][19] = 0.00345549;
			BSCornerCoeff[8][1] = 0.75304947;
			BSCornerCoeff[8][2] = -0.00394460;
			BSCornerCoeff[8][3] = -0.03874130;
			BSCornerCoeff[8][4] = 0.07195393;
			BSCornerCoeff[8][5] = 0.15915441;
			BSCornerCoeff[8][6] = -0.00103200;
			BSCornerCoeff[8][7] = 0.00051130;
			BSCornerCoeff[8][8] = -0.00007610;
			BSCornerCoeff[8][9] = 0.00059102;
			BSCornerCoeff[8][10] = 0.00414472;
			BSCornerCoeff[8][11] = -0.00626760;
			BSCornerCoeff[8][12] = -0.00225830;
			BSCornerCoeff[8][13] = 0.01422931;
			BSCornerCoeff[8][14] = -0.01478020;
			BSCornerCoeff[8][15] = -0.03936080;
			BSCornerCoeff[8][16] = -0.00007220;
			BSCornerCoeff[8][17] = -0.00103310;
			BSCornerCoeff[8][18] = -0.00015610;
			BSCornerCoeff[8][19] = 0.00204608;
			BSCornerCoeff[9][1] = 0.86826744;
			BSCornerCoeff[9][2] = -0.01540410;
			BSCornerCoeff[9][3] = -0.06757750;
			BSCornerCoeff[9][4] = -0.03465790;
			BSCornerCoeff[9][5] = 0.10827631;
			BSCornerCoeff[9][6] = -0.00003390;
			BSCornerCoeff[9][7] = 0.00177775;
			BSCornerCoeff[9][8] = 0.00086900;
			BSCornerCoeff[9][9] = 0.00039226;
			BSCornerCoeff[9][10] = 0.00834285;
			BSCornerCoeff[9][11] = 0.00051283;
			BSCornerCoeff[9][12] = -0.00852540;
			BSCornerCoeff[9][13] = -0.00794040;
			BSCornerCoeff[9][14] = -0.01895280;
			BSCornerCoeff[9][15] = -0.01331900;
			BSCornerCoeff[9][16] = 0.00024774;
			BSCornerCoeff[9][17] = -0.00044170;
			BSCornerCoeff[9][18] = 0.00193629;
			BSCornerCoeff[9][19] = 0.00554321;
			BSCornerCoeff[10][1] = 0.85992881;
			BSCornerCoeff[10][2] = -0.00117880;
			BSCornerCoeff[10][3] = -0.05706770;
			BSCornerCoeff[10][4] = 0.03259717;
			BSCornerCoeff[10][5] = 0.16789688;
			BSCornerCoeff[10][6] = -0.00017742;
			BSCornerCoeff[10][7] = 0.00041268;
			BSCornerCoeff[10][8] = -0.00120900;
			BSCornerCoeff[10][9] = 0.00075719;
			BSCornerCoeff[10][10] = 0.01013219;
			BSCornerCoeff[10][11] = -0.00276250;
			BSCornerCoeff[10][12] = -0.00340710;
			BSCornerCoeff[10][13] = 0.00021417;
			BSCornerCoeff[10][14] = -0.01798040;
			BSCornerCoeff[10][15] = -0.02177810;
			BSCornerCoeff[10][16] = -0.00015070;
			BSCornerCoeff[10][17] = -0.00123310;
			BSCornerCoeff[10][18] = -0.00005570;
			BSCornerCoeff[10][19] = 0.00222115;
			BSCornerCoeff[11][1] = 0.59593229;
			BSCornerCoeff[11][2] = -0.04732680;
			BSCornerCoeff[11][3] = -0.03112210;
			BSCornerCoeff[11][4] = 0.03256568;
			BSCornerCoeff[11][5] = 0.21272148;
			BSCornerCoeff[11][6] = 0.00202659;
			BSCornerCoeff[11][7] = 0.00500298;
			BSCornerCoeff[11][8] = 0.00590225;
			BSCornerCoeff[11][9] = -0.00119010;
			BSCornerCoeff[11][10] = -0.00301830;
			BSCornerCoeff[11][11] = -0.00349060;
			BSCornerCoeff[11][12] = -0.00087780;
			BSCornerCoeff[11][13] = -0.00436520;
			BSCornerCoeff[11][14] = -0.01597470;
			BSCornerCoeff[11][15] = -0.03786950;
			BSCornerCoeff[11][16] = -0.00082380;
			BSCornerCoeff[11][17] = -0.00140310;
			BSCornerCoeff[11][18] = 0.00155153;
			BSCornerCoeff[11][19] = 0.00220878;
			BSCornerCoeff[12][1] = 0.55161620;
			BSCornerCoeff[12][2] = -0.04071620;
			BSCornerCoeff[12][3] = -0.01822780;
			BSCornerCoeff[12][4] = 0.10359240;
			BSCornerCoeff[12][5] = 0.24491063;
			BSCornerCoeff[12][6] = 0.00285987;
			BSCornerCoeff[12][7] = 0.00422790;
			BSCornerCoeff[12][8] = 0.00407121;
			BSCornerCoeff[12][9] = -0.00045380;
			BSCornerCoeff[12][10] = -0.00146840;
			BSCornerCoeff[12][11] = -0.00608740;
			BSCornerCoeff[12][12] = 0.00187184;
			BSCornerCoeff[12][13] = 0.00199230;
			BSCornerCoeff[12][14] = -0.01458080;
			BSCornerCoeff[12][15] = -0.04651480;
			BSCornerCoeff[12][16] = -0.00103380;
			BSCornerCoeff[12][17] = -0.00224000;
			BSCornerCoeff[12][18] = -0.00109160;
			BSCornerCoeff[12][19] = 0.00074973;
			BSCornerCoeff[13][1] = 0.67369793;
			BSCornerCoeff[13][2] = -0.00747010;
			BSCornerCoeff[13][3] = -0.04190240;
			BSCornerCoeff[13][4] = 0.02388243;
			BSCornerCoeff[13][5] = 0.13744754;
			BSCornerCoeff[13][6] = 0.00223426;
			BSCornerCoeff[13][7] = 0.00279771;
			BSCornerCoeff[13][8] = 0.00605931;
			BSCornerCoeff[13][9] = -0.00144370;
			BSCornerCoeff[13][10] = -0.00151650;
			BSCornerCoeff[13][11] = -0.00385710;
			BSCornerCoeff[13][12] = -0.00620850;
			BSCornerCoeff[13][13] = 0.00184524;
			BSCornerCoeff[13][14] = -0.01126230;
			BSCornerCoeff[13][15] = -0.02505730;
			BSCornerCoeff[13][16] = -0.00075940;
			BSCornerCoeff[13][17] = -0.00146740;
			BSCornerCoeff[13][18] = 0.00142162;
			BSCornerCoeff[13][19] = 0.00275213;
			BSCornerCoeff[14][1] = 0.65043818;
			BSCornerCoeff[14][2] = -0.00626730;
			BSCornerCoeff[14][3] = -0.02744980;
			BSCornerCoeff[14][4] = 0.08847391;
			BSCornerCoeff[14][5] = 0.16448990;
			BSCornerCoeff[14][6] = 0.00098788;
			BSCornerCoeff[14][7] = 0.00209113;
			BSCornerCoeff[14][8] = 0.00376914;
			BSCornerCoeff[14][9] = -0.00137280;
			BSCornerCoeff[14][10] = 0.00063718;
			BSCornerCoeff[14][11] = -0.00626790;
			BSCornerCoeff[14][12] = -0.00473840;
			BSCornerCoeff[14][13] = 0.00404005;
			BSCornerCoeff[14][14] = -0.00997580;
			BSCornerCoeff[14][15] = -0.03069250;
			BSCornerCoeff[14][16] = -0.00043140;
			BSCornerCoeff[14][17] = -0.00210410;
			BSCornerCoeff[14][18] = -0.00083370;
			BSCornerCoeff[14][19] = 0.00182356;
			BSCornerCoeff[15][1] = 0.77184789;
			BSCornerCoeff[15][2] = -0.03865610;
			BSCornerCoeff[15][3] = -0.04342940;
			BSCornerCoeff[15][4] = -0.03465390;
			BSCornerCoeff[15][5] = 0.21522972;
			BSCornerCoeff[15][6] = 0.00393668;
			BSCornerCoeff[15][7] = 0.00409718;
			BSCornerCoeff[15][8] = 0.00543320;
			BSCornerCoeff[15][9] = -0.00118130;
			BSCornerCoeff[15][10] = 0.00001733;
			BSCornerCoeff[15][11] = 0.00349585;
			BSCornerCoeff[15][12] = 0.00013448;
			BSCornerCoeff[15][13] = -0.00067790;
			BSCornerCoeff[15][14] = -0.01977780;
			BSCornerCoeff[15][15] = -0.03873230;
			BSCornerCoeff[15][16] = -0.00097200;
			BSCornerCoeff[15][17] = -0.00200060;
			BSCornerCoeff[15][18] = 0.00151529;
			BSCornerCoeff[15][19] = 0.00269907;
			BSCornerCoeff[16][1] = 0.77682841;
			BSCornerCoeff[16][2] = -0.02884030;
			BSCornerCoeff[16][3] = -0.03734350;
			BSCornerCoeff[16][4] = 0.03387784;
			BSCornerCoeff[16][5] = 0.26006713;
			BSCornerCoeff[16][6] = 0.00155403;
			BSCornerCoeff[16][7] = 0.00309614;
			BSCornerCoeff[16][8] = 0.00338795;
			BSCornerCoeff[16][9] = -0.00062250;
			BSCornerCoeff[16][10] = 0.00132589;
			BSCornerCoeff[16][11] = 0.00080919;
			BSCornerCoeff[16][12] = 0.00325629;
			BSCornerCoeff[16][13] = 0.00538607;
			BSCornerCoeff[16][14] = -0.01955600;
			BSCornerCoeff[16][15] = -0.05034970;
			BSCornerCoeff[16][16] = -0.00123340;
			BSCornerCoeff[16][17] = -0.00247470;
			BSCornerCoeff[16][18] = -0.00076300;
			BSCornerCoeff[16][19] = 0.00120785;
		}	// InitBSCornerCoeff
		
		
	} // BaseSimp

} // EnergyPlus







