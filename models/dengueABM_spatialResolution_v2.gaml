/*
* Name:        Dengue ABM - Spatial Resolution (v2)
* Author:      Katharina Wöhs, Department of Geoinformatics - Z_GIS, University of Salzburg
* Description: Prototype ABM exploring the impact of grid cell size on dengue fever spread
*              in a given area of interest. NOT a prediction model.
* Based on:    salzburgDengueABM-v1 (Wöhs, 2025, master thesis)
* Tags:        epidemiological modelling, dengue, ABM, spatial simulation, MAUP
* Version:     2.0
* Last edit:   10-04-2026
*
* Documentation:
*   - README.md          — project overview, v2 changelog, how to run
*   - ODD.md             — ODD protocol description (Grimm et al. 2020)
*   - OPEN_QUESTIONS.md  — known limitations and discussion points
*/

// Model
model dengueABM_MAUP
																						// SECTION for INLINE DOCUMENTATION and COMMENTS
// GLOBAL
global {
	// MODEL AREA 																		// initialize field using a tiff raster file, see documentation: https://gama-platform.org/wiki/DataTypes#field
	// -- MAUP
	//string raster_path <- "../includes/masked_map_greenness_id79229-Elsbethen.tiff";  // risk raster from Sen2Cube in 10x10m; swap for different resolutions (10m/50m/100m/1000m/single-cell)
	string raster_path <- "../includes/habitat_50_mask.tif";							// factor 5 aggregation, mean
	//string raster_path <- "../includes/habitat_100_mask.tif";							// factor 10 aggregation, mean
	//string raster_path <- "../includes/habitat_500_mask.tif";							// factor 50 aggregation, mean
	//string raster_path <- "../includes/habitat_1000_mask.tif";						// factor 100 aggregation, mean
	//string raster_path <- "../includes/habitat_3000_mask_singlepixle.tif";				// factor 300 aggregation - single cell, minimum
	
	field field_from_sen2cube <- field(grid_file(raster_path));						// greenness May - Oct 2024
	file elsbethen_shp <- file("../includes/boundariesElsbethen.shp");					// add Elsbethen catastral community boundaries to create a geometry within which the simulation is defined
	file elsbethen_settlement_shp <- file ("../includes/permanentSettlementAreaElsbethen.shp"); // add Elsbethen permanent settlement areas to create human habitat geometry
	geometry shape <- envelope(elsbethen_shp);											// define shape with Elsbethen-shapefile
	geometry elsbethen_geom <- geometry(elsbethen_shp);									// define geometry for movement boundaries
	geometry settlement_geom <- geometry(elsbethen_settlement_shp);						// define geometry for movement boundaries - humans


	// --- MODEL PARAMETERS FOR SIMULATION ---
	float step <- 8 #hour;																// each time step equals 8h


	// --- HUMAN PARAMETERS ---
	// sir
	int nb_humans <- 3252;																// number of inhabitants
	int nb_infected_init <- 1;															// number of infected humans at the beginning of the simulation
	int nb_recovered_init <- 0;															// number of humans recovered at beginning (declared before nb_humans_recovered)
	int nb_humans_infected <- nb_infected_init update: humans count(each.is_infectious); // update the number of infectious humans (I state) after each time step
	int nb_humans_not_infected <- (nb_humans - nb_humans_infected);						// number of healthy people
	int nb_humans_recovered <- nb_recovered_init update: humans count(each.is_recovered); // update the number of recovered humans after each time step
	float infected_rate update: nb_humans_infected/nb_humans;							// rate of people currently infected with DENV
	// transmission probabilities
	float proba_infection_human <- 0.1;													// probability of human to get infected with bite
	float proba_getinfectious_human <- 0.5;												// probability of falling sick with a high viraemic phase and symptomes
	// speed
	float humans_min_speed <- 0.0 #m / #h;												// minimal speed for humans
	float humans_max_speed <- 200.0 #m / #h;											// maximal speed for humans
	
	
	// --- MOSQUITO PARAMETERS ---
	// sir
	int nb_mosquitos_infected;															// variable for infected mosquitoes
	int nb_mosquitos_infectious_init <- 0;												// number of infectious mosquitoes at the beginning
	int nb_mosquitos_infectious <- nb_mosquitos_infectious_init update: mosquitos count(each.is_infectious); // update number of infectious mosquitoes
	// transmission probabilities
	int eip_mosquito <- 8;																// extrinsic incubation period for mosquitoes in days (dengue: 8-12 days at optimal temperature)
	int gonotrophic_cycle <- 3;															// days between blood meals (Aedes, ~3 days at optimal temperature)
	float proba_infection_vector <- 0.9;												// probability of mosquito to get infected when biting a human
	float proba_getinfectious_vector <- 0.37;											// probability of virus to survive in mosquito, making it a vector
	float proba_mosquito_bite <- 0.1;													// probability of a mosquito biting an infectious human per 8h step in a risk area
	// speed 
	float vector_min_speed <- 0.0 #m / #h;												// minimal flying speed for mosquitoes
	float vector_max_speed  <- 50.0 #m / #h;											// maximal flying speed for mosquitoes
	float dispersal_radius <- 250.0 #m;													// home range radius; >250m documented for European Aedes populations (Vavassori, Saddler & Müller, 2019)
	float proba_dispersal <- 0.05;														// probability of breaking out of home range per 8h step; updates birth_location on dispersal


	// --- TRACKING & OUTPUT ---
	int time_of_peak <- 0;																// cycle at which nb_humans_infected was highest
	int peak_infected <- 0;																// highest nb_humans_infected observed during simulation
	bool localInfections <- false;														// flags simulation runs with local infection
	bool noLocalInfections <- false;													// flags simulation runs without local infection
	int nb_localInfections <- 0;														// cumulative count of mosquito-to-human transmission events
	int nb_nolocalInfections <- 0;														// count number of runs without local infection
	list<point> transmission_locations <- [];											// list of all locations where virus transmission to a human occurred
	string results_folder <- "../results/";												// output folder for CSV export (create this folder in your project directory)
	bool sim_over <- false;																// flag to signal clean simulation stop; used in until: condition to avoid do die auto-restart


	// --- INIT: SET UP WORLD AND AGENTS ---
	init {
        create humans number: nb_humans {
			location <- any_location_in(settlement_geom);								// create humans anywhere in settlement area
		}
		ask nb_infected_init among humans {												// set the infected traveller to infectious / change health status
			is_infected <- true;
			is_infectious <- true;
			is_susceptible <- false;
			healthstatus <- "I";
		}
    }
	action save_transmission_output {													// **ACTION SAVE_TRANSMISSION_OUTPUT**: export transmission locations to CSV
		string filepath <- results_folder + "transmission_" + string(seed) + ".csv";
		save "x,y" to: filepath;														// write header; always creates the file, even for runs with no transmission events
		loop pt over: transmission_locations {
			point real_pt <- point(CRS_transform(pt, "EPSG:3035"));						// convert from GAMA internal coordinates to ETRS89 LAEA
			save ("" + real_pt.x + "," + real_pt.y) to: filepath rewrite: false;		// append real-world coordinates
		}
	}

	reflex track_peak when: nb_humans_infected > peak_infected {						// **REFLEX TRACK_PEAK**: record cycle of epidemic peak
		peak_infected <- nb_humans_infected;
		time_of_peak <- cycle;
	}

	reflex end_season when: cycle = 539 and sim_over = false {							// **REFLEX END_SEASON**: fires at cycle 539, one cycle before until:540 terminates the run
		if nb_localInfections = 0 {														// partition runs based on actual transmission count, not on current infectious status
			noLocalInfections <- true;
			nb_nolocalInfections <- nb_nolocalInfections + 1;
			write "End of season. No local infection occurred.";
		} else {
			write "End of season. " + nb_localInfections + " local infection(s) occurred.";
		}
		do save_transmission_output;
		sim_over <- true;
	}
																						// **REFLEX END_SIMULATION_NOINFECTION**: disease died out before end of season
	reflex end_simulation_noinfection when: nb_humans_infected = 0 and nb_mosquitos_infectious = 0 and (mosquitos count (each.is_infected)) = 0 and sim_over = false {
		if nb_localInfections = 0 {														// outbreak never started — no mosquito ever bit an infectious human successfully
			noLocalInfections <- true;
			nb_nolocalInfections <- nb_nolocalInfections + 1;
			write "Outbreak never started. No local infection.";
		} else {																		// outbreak started but has now died out — localInfections flag was set in reflex infect, keep it
			write "Outbreak ended. " + nb_localInfections + " local infection(s) occurred.";
		}
		do save_transmission_output;
		sim_over <- true;
	}
}																						// END GLOBAL


// SPECIES: HUMANS
species humans skills: [moving] {
	rgb color;																			// color in display chart
	float size <- 1.0;																	// size in display chart
	float action_radius <- rnd(humans_min_speed, humans_max_speed);						// define action radius
	bool is_susceptible <- true;														// susceptible
	bool is_infectious <- false;														// infected with symptomes
	bool is_infected <- false;															// infected without symptomes - Neglectable
	bool is_quarantine <- false;														// infected with symptomes, out of system
	bool is_recovered <- false;															// recovered
	string healthstatus;																// create variable for health status
	int infectionDay <- 0;																// day when infection starts
	bool riskArea <- false;																// risk area initially set to false

	reflex move {																		// **REFLEX MOVE**
		riskArea <- false;																// reset risk area each step before re-evaluating current position
		ask cell overlapping self.location {											// ask cell at current position of human
			if self.grid_value <= 5 {													// all values of 5 and below are risk areas
				myself.riskArea <- true;												// set risk area to true
			}
    	//write self.location;															// console output for check-up
		//write self.bands;
		//write self.grid_value;
		//write myself.riskArea;
		}

		if riskArea = true and self.is_infectious = true and flip(proba_mosquito_bite) {	// TRIGGER **ACTION MAKE_MOSQUITO_AGENT**
			do make_mosquito_agent;
		}
		// Soft preference for settlement: humans live there, so when outside they usually head home,
		// but excursions into the wider Elsbethen area are allowed. Health status no longer branches
		// movement (the three previous branches were identical).
		if settlement_geom covers self.location {										// inside settlement — free local wander within Elsbethen
			do wander amplitude: 60.0 bounds: elsbethen_geom speed: action_radius;
		}
		else if flip(0.6) {																// outside settlement — usually head home (they live there)
			do goto target: any_location_in(settlement_geom) speed: action_radius;
		}
		else {																			// occasional longer excursion — keep wandering freely within Elsbethen
			do wander amplitude: 60.0 bounds: elsbethen_geom speed: action_radius;
		}
	}

	action make_mosquito_agent {														// **ACTION MAKE_MOSQUITO_AGENT**
		create mosquitos number: 1 {													// create only 1 mosquito
			set location <- myself.location;											// create mosquito at my location
			nb_mosquitos_infected <- nb_mosquitos_infected + 1;							// increase number infected mosquito for end condition
			birth_location <- location;													// record birth location as home range centre
			is_infected <- flip(proba_infection_vector);								// infect mosquito with given probability
			// is_infectious not set here: mosquito must complete EIP before becoming infectious
		}
	}

	reflex updateSIR {																	// **REFLEX UPDATE SIR for humans**

		if is_infectious = false and is_infected = false and is_recovered = false and (infectionDay/24) <= 0 {	// 'S' - susceptible; is_recovered guard prevents asymptomatic-R humans (infectionDay=0) being reset to S
			is_susceptible <- true;														// make agent susceptible
			healthstatus <- "S";
			color <- #limegreen;														// change colour for display
		}

		else if is_infectious = true and is_infected = true and (infectionDay/24) < 9 { // 'I' - infectious
			healthstatus <- "I";
			infectionDay <- infectionDay + 8;											// increase infection day by 8h per time step
			color <- #red;																// change colour for display
		}

		else if is_infectious = true and is_infected = true and (infectionDay/24) > 8 { // 'R' - recovered and immune
			healthstatus <- "R";														// make agent healthy again
			is_infectious <- false;
			is_infected <- false;
			is_susceptible <- false;
			is_recovered <- true;
			color <- #lightblue;														// change colour for display
		}
	}

	aspect base {																		// **ASPECT HUMANS**
		draw circle(size) color: color;													// health status is displayed via colours defined in if-elseif
	}
}																						// END HUMAN


// SPECIES: MOSQUITOS
species mosquitos skills: [moving] {
	float size <- 0.1;																	// size of mosquito
	float action_radius <- rnd(vector_min_speed, vector_max_speed);						// flying distance of mosquito
	humans my_human;																	// target to feed on
	geometry perception_area;															// geometry for movement within one timestep
	bool is_infected;																	// infected status
	bool is_infectious;																	// infectious status
	int eipDayVector <- 0;																// days spent in extrinsic incubation period (EIP)
	int age <- 0;																		// age of mosquito in hours since creation
	int hours_since_last_feed <- 0;														// hours since last blood meal; 0 at creation (make_mosquito_agent represents a feeding event)
	point birth_location;																// location where mosquito was created; defines home range centre
	bool is_susceptible;
	string healthstatus;																// variable for health status
	rgb color;																			// variable for color

	reflex updateVectorSIR {															// **REFLEX UPDATE SIR for mosquitoes**
		ask self {
			if self.is_infectious = true {
				color <- #red;															// change colour for display
				healthstatus <- "I";													// set health status
			}

			else if self.is_infected = true and self.is_infectious = false {			// EIP phase: infected but not yet infectious
				color <- #orange;														// change colour for display
				healthstatus <- "E";													// E = exposed / in EIP
				eipDayVector <- eipDayVector + 8;										// count hours spent in EIP
				if (eipDayVector / 24) >= eip_mosquito {								// EIP complete: mosquito may become infectious
					is_infectious <- flip(proba_getinfectious_vector);					// virus survives and mosquito becomes infectious with given probability
				}
			}

			else if self.is_infected = false {
				color <- #limegreen;													// change colour for display
				healthstatus <- "S";													// set health status
				is_susceptible <- true;
			}
		}
	}																					// mosquitoes do not recover: virus persists until mosquito dies
	
	reflex lifespan {																	// **REFLEX LIFESPAN**: mosquito dies after 45 days
		age <- age + 8;																	// increment age by 8h per time step
		hours_since_last_feed <- hours_since_last_feed + 8;								// increment gonotrophic cycle counter by 8h per time step
		if (age / 24) >= 45 {
			do die;
		}
	}

	reflex move {																		// **REFLEX MOVE**
		if flip(proba_dispersal) {														// breakout: move freely within Elsbethen and update home range
			do wander amplitude: 60.0 bounds: elsbethen_geom speed: action_radius;
			birth_location <- location;													// permanently relocate home range centre
		} else {																		// normal: move within 250m home range (unclipped circle reduces edge effects)
			geometry home_range <- circle(dispersal_radius) at_location birth_location;
			do wander amplitude: 60.0 bounds: home_range speed: action_radius;
		}
	}

	reflex hunt when: (hours_since_last_feed / 24) >= gonotrophic_cycle {				// **REFLEX HUNT**: only hunt when gonotrophic cycle cooldown is over
		perception_area <- (self + 80 #m) intersection cone (int(heading - 60), int(heading + 60)); // set heading and perception area for sensing
		my_human <- (humans overlapping perception_area) closest_to self;				// target human that is closest
		if (my_human != nil) {															// when human in range
			do goto target: my_human speed:action_radius;								// go to my_human
			//write "I " + self + "Got you: " + my_human;								// console check-up for tracing
		}
		else {																			// if no human in range, move within home range (unclipped circle reduces edge effects)
			geometry home_range <- circle(dispersal_radius) at_location birth_location;
			do wander amplitude: 60.0 bounds: home_range speed: action_radius;
		}
	}

	reflex infect when: my_human != nil {												// **REFLEX INFECT**
		ask my_human {																	// ask if human is susceptible
			if self.healthstatus = "I" or self.healthstatus = "R" {
				do wander;																// if not, nothing happens
			}
			else if myself.is_infectious = true and self.healthstatus = "S" {					// ask if mosquito is infectious and human is susceptible
				self.is_infected <- flip(proba_infection_human);						// probability human gets infected
				if self.is_infected {													// only if virus actually transferred to human
					nb_localInfections <- nb_localInfections + 1;						// increment cumulative local infection count
					localInfections <- true;											// flag this run as having local transmission
					transmission_locations <- transmission_locations + [myself.location]; // record location where virus transmission occurs (covers both symptomatic and asymptomatic)
					if flip(proba_getinfectious_human) {								// symptomatic path: transition to I state
						self.is_infectious <- true;
					} else {															// asymptomatic path: skip I, go straight to R (cleared virus without symptoms)
						self.is_infected <- false;
						self.is_susceptible <- false;
						self.is_recovered <- true;
						self.healthstatus <- "R";
						self.color <- #lightblue;
					}
				}
				//write "Help!" + self + "I got infected by" + myself + ". Now, I need to rest."; // console check-up for tracing
			}
		}
		hours_since_last_feed <- 0;														// reset gonotrophic cycle counter after feeding
		my_human <- nil;																// clear target so mosquito re-hunts after cooldown
	}

	aspect base {																		// **ASPECT MOSQUITOS**
		draw triangle(size) color: color;												// health status as defined in updateVectorSIR reflex
	}
}																						// END MOSQUITOS


// SPECIES: GRID																		// create grid to test if the field location matches the cell location
grid cell file: grid_file(raster_path) {												// grid dimensions and grid_value auto-populated from raster; adapts to any resolution (10m/50m/100m/1000m)
}																						// proven by > "visualisationOfTransmissionProcess" > Inspect cell
																						// END GRID


// EXPERIMENT: BATCH
experiment simulationBatch type: batch repeat: 16 until: cycle = 540 parallel: 16 {	// 100 Monte Carlo runs, parallelized across 16 cores; 540 cycles = 180 days = 6 months; do die handles early exit
	//method stochastic;

	permanent {
	    display Distribution background: #white type: 2d {
			chart "Histogram no local vs. local infections" type: histogram size: {0.5,0.5} position: {0, 0} {
				data "No local infections" value: simulations count (each.noLocalInfections) color: #green;
		     	data "Local infections" value: simulations count (each.localInfections) color: #red;
		    }
		}
  	}																					// close permanent
}																						// END EXPERIMENT "simulationBatch"


// EXPERIMENT: SINGLE RUN
experiment singleSim type: gui until: cycle = 540 {									// hard cap at 540 cycles = 180 days = 6 months; do die handles early exit
	output {
		monitor "Current hour" value: current_date.hour;
		monitor "Infected people rate" value: nb_infected_init;
		monitor "Time of peak" value: time_of_peak;
		monitor "Peak infected" value: peak_infected;
		monitor "Local infections" value: nb_localInfections;

		// map: agents in AoI
		display display_grid type: 2d {													// creates a visual display for the Tracing to follow agents around
			grid cell border: #black;													// creates a grid in the visualisation
			mesh field_from_sen2cube color: #lightgreen scale:0.0;						// displays the raster in different green colour values depending on the pixel value
			species humans aspect:base;													// visualises humans in the display grid
			species mosquitos aspect:base;												// visualises mosquitos in the display grid
			graphics "transmission locations" {											// display all transmission locations as points
				loop pt over: transmission_locations {
					draw circle(30) at: pt color: #pink border: #black;
				}
			}
		}

		// diagram: numbers infected
		display chart refresh: every(1 #cycles) type: 2d {								// create chart with statistics on infection state
			chart "Disease spread" type: series y_range: {0, nb_humans} {
				data "infected humans" value: nb_humans_infected color: #red marker: false; // display number of infected humans in red
			}
		}
		
		//map: locations transmission 
		display transmission_map type: 2d {												// standalone map of transmission locations for results output
			mesh field_from_sen2cube color: #lightgreen scale: 0.0;						// background raster for spatial reference
			graphics "transmission locations" refresh: every(1 #cycles) {				// accumulates transmission points over the full simulation run
				loop pt over: transmission_locations {
					draw circle(30) at: pt color: #pink border: #black;
				}
			}
		}
	}
}																						// END EXPERIMENT "visualisationOfTransmissionProcess"


// EXPERIMENT: OPENMOLE (minimal, headless — called by OpenMOLE via GAMATask)
// No displays: OpenMOLE only reads the experiment section for parameters and monitors.
// Input:  raster_path (exposed as parameter so OpenMOLE can sweep across resolutions)
// Output: Time of peak, Peak infected, Local infections (read by OpenMOLE via monitor LABELS)
experiment openmoleEx type: gui until: sim_over or cycle >= 540 {
	parameter "raster_path" var: raster_path;											// expose raster input for OpenMOLE direct sampling
	output {
		monitor "Time of peak" value: time_of_peak;
		monitor "Peak infected" value: peak_infected;
		monitor "Local infections" value: nb_localInfections;
	
		//map: locations transmission - to check if raster-path swap works
//		display transmission_map type: 2d {												// standalone map of transmission locations for results output
//			mesh field_from_sen2cube color: #lightgreen scale: 0.0;						// background raster for spatial reference
//			graphics "transmission locations" refresh: every(1 #cycles) {				// accumulates transmission points over the full simulation run
//				loop pt over: transmission_locations {
//					draw circle(30) at: pt color: #pink border: #black;
//				}
//			}
//		}
	}
	
	
}																						// END EXPERIMENT "openmole"
																						// END OF SIMULATION model dengueABM_MAUP
