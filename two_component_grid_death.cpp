#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <cstdio>
#include <fstream>

using namespace std;

// global parameters
double tau1; // incubation time
double tau2; // time from infection to recovery; tau2 > tau1
double tau3; // characteristic vehicle infectivity decay time
double personContactRate;
double vehicleDemandRate;
double disinfectionRate;
double lambda; // exponent to compute infection probability
double deathRate; // infected person dies with rate deathRate times
		  // the difference of activity from 1

// nearest neighbor offsets--can add diagonal neighbors or next
// nearest neighbors easily
vector< pair<int,int> > offsets = {{-1,0},{0,1},{1,0},{0,-1}};

// random number generator
default_random_engine generator;
uniform_real_distribution<double> randDist(0.,1.);

// function to compute whether transmission happens
bool transmit(double infectivity) {
  return randDist(generator) > exp(-lambda*infectivity);
}

class Person {
  double severity; // determines how much activity is suppressed
public:
  bool vulnerable;
  bool dead;
  double infectedTime;
  Person() {
    vulnerable = true;
    dead = false;
    // severity = 1.;
    severity = randDist(generator); // random uniform [0:1)
  }
  double infectivity(double time) {
    if (vulnerable or dead) return 0.;
    double deltaT = time - infectedTime; // time since infection
    if (deltaT >= tau2) return 0.; // disease ran its course
    if (deltaT < tau1) {
      return deltaT/tau1;
    } else {
      return (tau2 - deltaT)/(tau2 - tau1);
    }
  }
  double activity(double time) {
    if (dead) return 0.;
    return 1. - this->infectivity(time)*severity;
  }
  void infect(double time) {
    vulnerable = false;
    infectedTime = time;
  }
  bool infectious(double time) {
    return not vulnerable and (time - infectedTime < tau2);
  }
  bool willInfect(double time) {
    // use the global transmission function
    return transmit(this->infectivity(time));
  }
};

class Vehicle {
  double baseInfectivity;
public:
  double lastUpdateTime;
  Vehicle() {
    baseInfectivity = 0.;
    lastUpdateTime = 0.;
  }
  double infectivity(double time) {
    return baseInfectivity*exp((lastUpdateTime - time)/tau3);
  }
  void updateInfectivity(double time,double infectivity) {
    baseInfectivity = this->infectivity(time) + infectivity;
    lastUpdateTime = time;
  }
  void disinfect(double time) {
    baseInfectivity = 0.;
    lastUpdateTime = time;
  }
  bool infectious() {
    return baseInfectivity > 0.;
  }
  bool willInfect(double time) {
    return transmit(this->infectivity(time));
  }
};

class System {
  int numPeopleX;
  int numPeopleY;
  int numPeople;
  int numVehicles;
  vector<Person*> people; // rows first
  vector<Vehicle*> vehicles;
  vector<double> rates;
public:
  double time;
  System(int NPx,int NPy,int NV) {
    numPeopleX = NPx;
    numPeopleY = NPy;
    numPeople = NPx*NPy;
    numVehicles = NV;
    // simulation time
    time = 0.;
    // init people and vehicles
    for (int row = 0; row < numPeopleY; ++row) {
      for (int col = 0; col < numPeopleX; ++col) {
	Person* person = new Person();
	people.push_back(person);
      }
    }
    for (int i = 0; i < numVehicles; ++i) {
      Vehicle* vehicle = new Vehicle();
      vehicles.push_back(vehicle);
    }
    // infect a random person
    this->people[randomPersonIndex()]->infect(time);
  }
  double fracAffected() {
    int numAffected = 0;
    for (auto person : people) {
      if (not person->vulnerable) ++numAffected;
    }
    return double(numAffected)/numPeople;
  }
  double avgVehicleInfectivity() {
    double avg = 0.;
    for (auto vehicle : vehicles) {
      avg += vehicle->infectivity(time);
    }
    return avg/numVehicles;
  }
  double lastInfectionTime() {
    double lastTime = 0.;
    for (auto person : people) {
      if (person->infectedTime > lastTime)
	lastTime = person->infectedTime;
    }
    return lastTime;
  }
  double totalActivity() {
    double activity = 0;
    for (auto person : people) {
      activity += person->activity(time);
    }
    return activity;
  }
  double totalSeverity() {
    double severity = 0;
    for (auto person : people) {
      if (not person->dead) {
	severity += 1. - person->activity(time);
      }
    }
    return severity;
  }
  double fracDead() {
    int n = 0;
    for (auto person : people) {
      if (person->dead) ++n;
    }
    return double(n)/numPeople;
  }
  int randomPersonIndex() {
    vector<double> activities;
    activities.clear();
    double total = 0.;
    for (auto person : people) {
      total += person->activity(time);
      activities.push_back(total);
    }
    double shot = total*randDist(generator);
    int i = 0;
    while (activities[i] < shot) ++i;
    return i;
  }
  int randomPersonToDieIndex() {
    vector<double> severities;
    severities.clear();
    double total = 0.;
    for (auto person : people) {
      if (not person->dead) {
	total += 1. - person->activity(time);
      }
      severities.push_back(total);
    }
    double shot = total*randDist(generator);
    int i = 0;
    while (severities[i] < shot) ++i;
    return i;
  }
  pair<int,int> coords(int idx) { // row,col
    int row = idx/numPeopleX;
    int col = idx - row*numPeopleX;
    return make_pair(row,col);
  }
  int offsetIndex(int centerIdx,pair<int,int> offset) {
    // periodic boundary conditions in both directions (simularion is
    // on a torus)

    // given the index of the center and a pair of x,y offsets, return
    // the index of the location that is offset by offset.first in x
    // and offset.second in y
    // coordinates that correspont to index idx
    auto centerCoords = this->coords(centerIdx);
    int row = (centerCoords.first + offset.first + numPeopleY) % numPeopleY;
    int col = (centerCoords.second + offset.second + numPeopleX) % numPeopleX;
    return numPeopleX*row + col;
  }
  int randomNeighborIndex(int idx) {
    // make two vectors: 1) indices of neighbors and 2) cumulative
    // activity
    vector<int> neighbors;
    vector<double> cumActivity;
    double totAct = 0.;
    for (auto offset : offsets) {
      int neighbor = this->offsetIndex(idx,offset);
      neighbors.push_back(neighbor);
      totAct += people[neighbor]->activity(time);
      cumActivity.push_back(totAct);
    }
    // when all neighbors are dead return -1
    if (totAct == 0.) return -1;
    // pick a random neighbor with probability proportional to
    // activity
    double shot = totAct*randDist(generator);
    int i = 0;
    while (cumActivity[i] < shot) ++i;
    return neighbors[i];
  }    
  bool continueSimulation() {
    // stop if there are no vulnerable people or when nothing and
    // nobody is infectious
    bool existInfectious = false;
    bool existVulnerable = false;
    for (auto person : people) {
      if (person->infectious(time)) existInfectious = true;
      if (person->vulnerable) existVulnerable = true;
      if (existVulnerable and existInfectious) break;
    }
    // check vehicles if all people are not infectious and there are
    // still vulnerable people
    if (not existInfectious and existVulnerable) {
      for (auto vehicle : vehicles) {
	if (vehicle->infectious()) {
	  existInfectious = true;
	  break;
	}
      }
    }
    return existVulnerable and existInfectious;
  }
  void nextEvent() {
    // compute the cumuative rates
    double totActivity = this->totalActivity();
    rates.clear();
    rates.push_back(personContactRate*totActivity);
    rates.push_back(rates.back() + vehicleDemandRate*totActivity);
    rates.push_back(rates.back() + disinfectionRate*numVehicles);
    rates.push_back(rates.back() + deathRate*this->totalSeverity());
    double totRate = rates.back();
    // for (auto r : rates) cerr << r << " "; cerr << endl;
    
    // advance simulation time; make sure random number is not 0;
    double shot;
    do {
      shot = randDist(generator);
    } while (shot == 0.);
    time -= log(shot)/totRate;;

    // decide which event happens
    shot = randDist(generator)*totRate;

    int i = 0;
    while (rates[i] < shot) ++i;

    if (i == 0) {
      // cerr << "person to person\n";
      // person to person contact, select a person and the nearest
      // neighbor with probabilities proportional to their activities
      int i1,i2;
      do {
	i1 = this->randomPersonIndex();
	i2 = this->randomNeighborIndex(i1);
      } while (i2 < 0);
      auto p1 = people[i1];
      auto p2 = people[i2];
      // infection may happen if p1 is infectious and the p2 is
      // vulnerable or vice versa
      if (p1->infectious(time) and p2->vulnerable) {
	if (p1->willInfect(time)) p2->infect(time);
      }
      if (p2->infectious(time) and p1->vulnerable) {
	if (p2->willInfect(time)) p1->infect(time);
      }
      return;
    }
    if (i == 1) {
      // cerr << "person to vehicle\n";
      // person to vehicle contact
      auto person = people[this->randomPersonIndex()];
      auto vehicle = vehicles[int(randDist(generator)*numVehicles)];
      // if person is infectious, update vehicle infectivity
      if (person->infectious(time)) {
	vehicle->updateInfectivity(time,person->infectivity(time));
      }
      // if person is vulnerable maybe infect it
      if (person->vulnerable) {
	if (vehicle->willInfect(time)) person->infect(time);
      }
      return;
    }
    if (i == 2) {
      // cerr << "vehicle disinfection\n";
      // vehicle disinfection; TODO: a growing disinfection rate to
      // approximate regular disinfections; for now, select a random
      // vehicle to disinfect
      auto vehicle = vehicles[int(randDist(generator)*numVehicles)];
      vehicle->disinfect(time);
      return;
    }
    if (i == 3) {
      // cerr << "death\n";
      // pick a random person to kill
      people[this->randomPersonToDieIndex()]->dead = true;
    }
  }
  void outputIntermediate() {
    double avgPersonInfectivity = 0.;
    int vulnerable = 0;
    int dead = 0;
    for (auto person : people) {
      if (person->vulnerable) {
	++vulnerable;
      } else if (person->dead) {
	++dead;
      } else {
	avgPersonInfectivity += person->infectivity(time);
      }
    }
    cout << "Day " << int(time) << ": "
	 << vulnerable << " vulnerable; "
	 << dead << " dead; "
	 << "average person infectivity: " 
	 << avgPersonInfectivity/(numPeople - dead - vulnerable)
	 << "; average vehicle infectivity: " 
	 << this->avgVehicleInfectivity() << endl;
  }
  void outputVehiclesInfectivity() {
    cout << time;;
    for (auto vehicle : vehicles) {
      cout << " " << vehicle->infectivity(time);
    }
    cout << endl;
  }
  void makeFrame(int frameIndex) {
    // make the frame file name
    char fname[15];
    sprintf(fname,"frame_%04d.gp",frameIndex);
    fstream fh;
    fh.open(fname,fstream::out);

    // reuse the variable to make the output name
    sprintf(fname,"frame_%04d.png",frameIndex);

    // set up
    double xbuffer = 0.05;
    fh << "set size nosquare " << endl
       << "set xran [0:" << 1.333 + xbuffer << "]" << endl
       << "set yran [0:1]" << endl
       << "unset border" << endl
       << "unset xtics" << endl
       << "unset ytics" << endl
       << "unset key" << endl
       << "unset colorbox" << endl
       << "load 'moreland.pal'" << endl
       << "set title 'Day " << int(time) + 1 << "'" << endl
       << "set term pngcairo" << endl
       << "set out '" << fname << "'" << endl;

    // set up the grid dimensions and the object scales; devote the
    // upper 80% of real estate to people and the lower 20% to
    // vehicles
    double peopleSpacing = 1./numPeopleY;
    double peopleRadius = peopleSpacing*0.4;
    double rightFrac = 1.333 - xbuffer - numPeopleX*peopleSpacing;
    double vehicleAspectRatio = 2.;
    int vehicleNumRows = int(sqrt(numVehicles/rightFrac)) + 1;
    int vehicleNumCols = numVehicles/vehicleNumRows + 1;
    double vehicleYspacing = 1./vehicleNumRows;
    double vehicleXspacing = rightFrac/vehicleNumCols;
    double vehicleYsize = vehicleYspacing*0.8;
    double vehicleXsize = vehicleYsize/vehicleAspectRatio;

    // write a circle for each person, affect with fs solid 1 vulnerable
    // fs solid 0.5; vehicle as rectangles
    int idx;
    for (idx = 0; idx < numPeople; ++idx) {
      auto coords = this->coords(idx);
      int row = coords.first;
      int col = coords.second;
      auto person = people[idx];
      fh << "set object " << idx + 1 << " circle center ";
      fh << (col + 0.5)*peopleSpacing << "," << (row + 0.5)*peopleSpacing;
      fh << " radius " << peopleRadius;
      fh << " fs solid 1 ";
      if (person->vulnerable) {
	fh << "fc \"gray\"" << endl;
      } else if (person->dead) {
	fh << "fc \"black\"" << endl;
      } else {
	fh << "fc palette frac " << person->infectivity(time) << endl;
      }
    }
    // write rectangles for vehicles
    idx = 0;
    for (auto vehicle : vehicles) {
      int row = idx/vehicleNumCols;
      int col = idx - row*vehicleNumCols;
      ++idx;
      fh << "set object " << numPeople + idx << " rectangle center ";
      fh << 1.333 - rightFrac + xbuffer + (col + 0.5)*vehicleXspacing
	 << "," << (row + 0.5)*vehicleYspacing
	 << " size " << vehicleXsize << "," << vehicleYsize;
      double frac = vehicle->infectivity(time);
      if (frac > 1.) frac = 1.;
      fh << " fs solid fc palette frac " << frac << endl;
    }
    fh << "plot 1.05 lt -1 lw 0.01" << endl;
    fh.close();
  }    
};

int main(int argc, char* argv[]) {
  if (argc != 12) {
    cout << "Usage: " << argv[0]
	 << " numPeopleX numPeopleY numVehicles tau1 tau2 tau3"
	 << " personContactRate vehicleDemandRate"
	 << " disinfectionRate deathRate lambda" << endl;
    return 1;
  }

  // assign parameters
  int numPeopleX = atoi(argv[1]);
  int numPeopleY = atoi(argv[2]);
  int numVehicles = atoi(argv[3]);
  tau1 = atof(argv[4]);
  tau2 = atof(argv[5]);
  tau3 = atof(argv[6]);
  personContactRate = atof(argv[7]);
  vehicleDemandRate = atof(argv[8]);
  disinfectionRate = atof(argv[9]);
  deathRate = atof(argv[10]);
  lambda = atof(argv[11]);

  // swap numPeopleX and numPeopleY if needed to make sure numPeopleX
  // <= numPeopleY;
  if (numPeopleY < numPeopleX) swap(numPeopleX,numPeopleY);

  // output parameters
  cout << "Number of people X " << numPeopleX << endl;
  cout << "Number of people Y " << numPeopleY << endl;
  cout << "Number of vehicles " << numVehicles << endl;
  cout << "tau1 " << tau1 << endl;
  cout << "tau2 " << tau2 << endl;
  cout << "tau3 " << tau3 << endl;
  cout << "personContactRate " << personContactRate << endl;
  cout << "vehicleDemandRate " << vehicleDemandRate << endl;
  cout << "disinfectionRate " << disinfectionRate << endl;
  cout << "deathRate " << deathRate << endl;
  cout << "lambda " << lambda << endl;
  
  // random seed
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  generator.seed(seed);
  
  // init the system
  System system(numPeopleX,numPeopleY,numVehicles);

  // frame time in days
  double frameTime = 0.1;
  double time = 0.;
  int frameIndex = 0;
  system.makeFrame(frameIndex++);
  // simulate
  while (system.continueSimulation()) {
    system.nextEvent();
    if (system.time > time + frameTime) {
      system.makeFrame(frameIndex++);
      // system.outputVehiclesInfectivity();
      // system.outputIntermediate();
      time = system.time;
    }
  }

  // final output
  cout << "Epidemic took " << system.lastInfectionTime() + tau2 << " days; "
       << "Fraction infected: " << system.fracAffected()
       << "; Fraction dead: " << system.fracDead() << endl;
}
