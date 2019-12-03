#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <stack>
using namespace std;

const int FLOAT_OUTPUT_PRECISION = 7;
const int LATTICE_SIZE = 100;                           // Has to be an even integer
const long long LATTICE_POINTS = LATTICE_SIZE * LATTICE_SIZE; // Just for readability
const long long PRINT_STEP = 50 * LATTICE_POINTS;
const long long VISUALIZE_STEP = 10 * LATTICE_POINTS;
const float MF_STEP = 0.01; //0.05

// Equilibrium will be assumed to have arrived, when the difference between the average of the hydrocarbon production in the last EQUILIBRIUM_VERIFICATION_STEPS and
// the average of the latest and the oldest steps among the last EQUILIBRIUM_VERIFICATION_STEPS is less than EQUILIBRIUM_VERIFICATION_THRESHOLD
const int EQUILIBRIUM_VERIFICATION_STEPS = 100 * LATTICE_POINTS;
const int EQUILIBRIUM_VERIFICATION_THRESHOLD = -1; //Completely removing early equilibrium criteria
const int MAX_ITERATIONS_PER_SIMULATION = 5000 * LATTICE_POINTS;

bool visualize = true, info_messages = true, essential_info_messages = true;
default_random_engine generator;
long long hydrocarbon_production = 0;
array<array<string, LATTICE_SIZE>, LATTICE_SIZE> lattice;
array<array<array<pair<int, int>, 6>, LATTICE_SIZE>, LATTICE_SIZE> nearest_neighbours;
array<array<set<pair<int, int>>, LATTICE_SIZE>, LATTICE_SIZE> bonded_H_atoms;
array<array<pair<int, int>, LATTICE_SIZE>, LATTICE_SIZE> is_connected_to;
map<int, int> hydrocarbon_sizes_seen;
vector<string> product_smiles;

// Some useful overloads
ostream &operator<<(ostream &out, array<array<string, LATTICE_SIZE>, LATTICE_SIZE> &lattice)
{
	// for(auto &x: lattice)
	//     for(auto &y:x)
	//         out << y << ' ';
	for(int i=0; i<LATTICE_SIZE; i++)
		for(int j=0; j<LATTICE_SIZE; j++)
			out << lattice[i][j] <<(bonded_H_atoms[i][j].size()>0?("H"+(bonded_H_atoms[i][j].size()>1?to_string(bonded_H_atoms[i][j].size()):"")):"") << ' ';
	return out;
}

class LatticeOperations
{
public:
	// To be called once in the beginning of the program
	static void initialize()
	{
		lattice_operations = {
			initiate_reaction_with_CO_if_possible,
			initiate_reaction_with_CHx_if_possible,
			initiate_desorption_on_CHx_if_possible
		};
		lattice_operations_queues.resize(lattice_operations.size());
	}

	static void process_site(pair<int,int> site)
	{
		lattice_operations_queues[0].push(site);
		process_all_queues();
	};

private:
	// Vector of all lattice operations in their priority order
	static vector<function<void(pair<int, int>)>> lattice_operations;

	// Vector of queues of sites to be processed, one queue per operation
	static vector<queue<pair<int, int>>> lattice_operations_queues;

	static void process_all_queues()
	{
		bool processed = false;
		do
		{
			processed = true;
			for(int i=0; i<lattice_operations.size(); i++)
			{
				if(!lattice_operations_queues[i].empty())
				{
					pair<int, int> site = lattice_operations_queues[i].front();
					lattice_operations_queues[i].pop();
					lattice_operations[i](site);
					if(i+1 != lattice_operations.size())
						lattice_operations_queues[i+1].push(site);
					processed = false;
					break;
				}
			}
		} while (!processed);
	}

	// Initiates reaction on a CO site if possible
	static void initiate_reaction_with_CO_if_possible(pair<int, int> site)
	{
		// If two H atoms are nearest neighbours to a CO site, form and desorb H2O
		if (lattice[site.first][site.second] == "CO")
		{
			vector<pair<int, int>> nearest_sites_with_H;
			for (auto nearest_site : nearest_neighbours[site.first][site.second])
				if (lattice[nearest_site.first][nearest_site.second] == "H") // && is_connected_to[nearest_site.first][nearest_site.second] == make_pair(-1,-1))
					nearest_sites_with_H.push_back(nearest_site);
			if (nearest_sites_with_H.size() >= 2)
			{
				lattice[site.first][site.second] = "C";
				shuffle(nearest_sites_with_H.begin(), nearest_sites_with_H.end(), generator);
				for (int i = 0; i < 2; i++)
				{
					lattice[nearest_sites_with_H[i].first][nearest_sites_with_H[i].second] = "empty";
					if (is_connected_to[nearest_sites_with_H[i].first][nearest_sites_with_H[i].second] != make_pair(-1, -1))
					{
						pair<int, int> connected_C = is_connected_to[nearest_sites_with_H[i].first][nearest_sites_with_H[i].second];
						bonded_H_atoms[connected_C.first][connected_C.second].erase(nearest_sites_with_H[i]);
						// Pushing to queue 0 conserves the interchangability of the priority of functions without affecting the algorithm
						lattice_operations_queues[0].push(connected_C);
					}
					is_connected_to[nearest_sites_with_H[i].first][nearest_sites_with_H[i].second] = make_pair(-1, -1);
				}
			}
		}
	}

	// Initiates reaction on a C, CH, CH2, CH3 site if possible
	static void initiate_reaction_with_CHx_if_possible(pair<int, int> site)
	{
		// If one or more H atoms are nearest neighbours to a CHx site, increase x and possibly desorb CH4
		if (lattice[site.first][site.second] == "C" && bonded_H_atoms[site.first][site.second].size() != 4)
		{
			vector<pair<int, int>> nearest_sites_with_H;
			for (auto neighbour_site : nearest_neighbours[site.first][site.second])
				if (lattice[neighbour_site.first][neighbour_site.second] == "H" && is_connected_to[neighbour_site.first][neighbour_site.second] == make_pair(-1, -1))
					nearest_sites_with_H.push_back(neighbour_site);
			shuffle(nearest_sites_with_H.begin(), nearest_sites_with_H.end(), generator);
			for (auto H_site : nearest_sites_with_H)
			{
				bonded_H_atoms[site.first][site.second].insert(H_site);
				is_connected_to[H_site.first][H_site.second] = site;
				if (bonded_H_atoms[site.first][site.second].size() == 4)
					break;
			}
		}
	}

	// Search branched chains
	// Each element of unsatisfied_valency_sites is a pair of a site and number of its unsatisfied valencies
	static bool search_branched_chains(stack<pair<pair<int, int>, int>> &unsatisfied_valency_sites, set<pair<int, int>> &vis, string &smile)
	{
		if(unsatisfied_valency_sites.empty())
			return true;

		pair<pair<int, int>, int> stack_top = unsatisfied_valency_sites.top();
		pair<int, int> site = stack_top.first;
		int num_unsatisfied_valency = stack_top.second;
		
		if(num_unsatisfied_valency == 0)
		{
			unsatisfied_valency_sites.pop();
			smile += ")";
			if(search_branched_chains(unsatisfied_valency_sites, vis, smile))
				return true;
			unsatisfied_valency_sites.push(stack_top);
			smile.pop_back();
		}
		else
		{
			auto shuffled_neighbours = nearest_neighbours[site.first][site.second];
			shuffle(shuffled_neighbours.begin(), shuffled_neighbours.end(), generator);
			for(auto neighbour_site : shuffled_neighbours)
			{
				if(vis.count(neighbour_site) || lattice[neighbour_site.first][neighbour_site.second] != "C")
					continue;
				for(int outgoing_valency = 1; outgoing_valency != 4 && outgoing_valency <= min(num_unsatisfied_valency, 4 - (int) bonded_H_atoms[neighbour_site.first][neighbour_site.second].size()); outgoing_valency++)
				{
					// Remove it
					// if(outgoing_valency != num_unsatisfied_valency)
					// 	continue;
					unsatisfied_valency_sites.top().second -= outgoing_valency;
					unsatisfied_valency_sites.push(make_pair(neighbour_site, 4 - bonded_H_atoms[neighbour_site.first][neighbour_site.second].size() - outgoing_valency));
					vis.insert(neighbour_site);
					smile += "(";
					if(outgoing_valency == 2)
						smile += "=";
					else if(outgoing_valency == 3)
						smile += "#";
					smile += "C";
					if(search_branched_chains(unsatisfied_valency_sites, vis, smile))
						return true;
					unsatisfied_valency_sites.pop();
					unsatisfied_valency_sites.top().second += outgoing_valency;
					vis.erase(neighbour_site);
					for(int i = 0; i < 2 + (outgoing_valency == 2 || outgoing_valency == 3); i++)
						smile.pop_back();
				}
			}
		}

		return false;
	}

	// Initiates desorption on a C, CH, CH2, CH3 site if possible
	static void initiate_desorption_on_CHx_if_possible(pair<int, int> site)
	{
		if(lattice[site.first][site.second] != "C")
			return;

		set<pair<int, int>> vis_nodes{site};
		stack<pair<pair<int, int>, int>> unsatisfied_valency_sites({make_pair(site, 4 - bonded_H_atoms[site.first][site.second].size())});
		string smile("C");
		
		if(!search_branched_chains(unsatisfied_valency_sites, vis_nodes, smile))
			return;

		smile.pop_back();

		hydrocarbon_production++;
		hydrocarbon_sizes_seen[vis_nodes.size()]++;
		product_smiles.push_back(smile);
		for (auto C_site : vis_nodes)
		{
			for (auto H_site : bonded_H_atoms[C_site.first][C_site.second])
			{
				lattice[H_site.first][H_site.second] = "empty";
				is_connected_to[H_site.first][H_site.second] = make_pair(-1, -1);
			}
			lattice[C_site.first][C_site.second] = "empty";
			bonded_H_atoms[C_site.first][C_site.second].clear();
		}
	}
};

vector<function<void(pair<int, int>)>> LatticeOperations::lattice_operations;

// Initialize the vector of queues
vector<queue<pair<int, int>>> LatticeOperations::lattice_operations_queues(3);

int main(int argc, char *argv[])
{
	auto start_time = chrono::high_resolution_clock::now();

	ios_base::sync_with_stdio(0);

	LatticeOperations::initialize();

	// Process command line options
	for (int i = 1; i < argc; i++)
		if (strcmp(argv[i], "--novis") == 0)
			visualize = false;
		else if (strcmp(argv[i], "--noinfo") == 0)
			info_messages = false;
		else if (strcmp(argv[i], "--noessentialinfo") == 0)
			essential_info_messages = false;
		else
			return cerr << "Error: unknown option '" << argv[i] << "' \n", 1;

	ofstream hydrocarbon_production_outfile("hydrocarbon_production.csv"), hydrocarbon_sizes_outfile("hydrocarbon_sizes_distribution.txt"), product_smiles_outfile("product_smiles.txt");
	hydrocarbon_production_outfile << fixed;
	hydrocarbon_sizes_outfile << fixed;
	hydrocarbon_production_outfile.precision(FLOAT_OUTPUT_PRECISION);
	hydrocarbon_sizes_outfile.precision(FLOAT_OUTPUT_PRECISION);

	hydrocarbon_production_outfile << "Mole Fraction of CO,Hydrocarbon production,Coverage fraction of CO,Coverage fraction of H,Coverage fraction of C,Coverage fraction of empty sites,Average hydrocarbon production\n";

	if (essential_info_messages)
	{
		cerr << fixed;
		cerr.precision(FLOAT_OUTPUT_PRECISION);
	}

	// Initialize the random generator with a time based seed (namely the current time)
	generator = default_random_engine(
		chrono::system_clock::now().time_since_epoch().count());

	// Initialize nearest neighbours
	for (int i = 0; i < LATTICE_SIZE; i++)
		for (int j = 0; j < LATTICE_SIZE; j++)
		{
			// Right
			nearest_neighbours[i][j][0] = (make_pair(i, (j + 1) % LATTICE_SIZE));
			// Left
			nearest_neighbours[i][j][1] = (make_pair(i, (j - 1 + LATTICE_SIZE) % LATTICE_SIZE));
			// Bottom-right
			nearest_neighbours[i][j][2] = (make_pair((i + 1) % LATTICE_SIZE, (j + (i & 1)) % LATTICE_SIZE));
			// Bottom-left
			nearest_neighbours[i][j][3] = (make_pair((i + 1) % LATTICE_SIZE, (j + (i & 1) - 1 + LATTICE_SIZE) % LATTICE_SIZE));
			// Up-right
			nearest_neighbours[i][j][4] = (make_pair((i - 1 + LATTICE_SIZE) % LATTICE_SIZE, (j + (i & 1)) % LATTICE_SIZE));
			// Up-left
			nearest_neighbours[i][j][5] = (make_pair((i - 1 + LATTICE_SIZE) % LATTICE_SIZE, (j + (i & 1) - 1 + LATTICE_SIZE) % LATTICE_SIZE));
		}

	if (visualize)
		// Print Lattice Size (One Time)
		cout << LATTICE_SIZE << '\n';

	enum class Reactant
	{
		CO,
		H2
	};

	enum class Direction
	{
		right,
		bottom_right,
		bottom_left
	};

	for (float mole_fraction_of_CO = MF_STEP; mole_fraction_of_CO <= 1.001; mole_fraction_of_CO += MF_STEP)
	{
		if (essential_info_messages)
			cerr << "Mole fraction of CO: " << mole_fraction_of_CO << '\n';

		// Initialize the lattice properties
		for (auto &row : lattice)
			for (auto &element : row)
				element = "empty";

		for (auto &row : is_connected_to)
			for (auto &element : row)
				element = make_pair(-1, -1);

		for (auto &row : bonded_H_atoms)
			for (auto &element : row)
				element.clear();

		// Variables for Main Loop
		int n_trials = 0, n_succesful_trials = 0;
		uniform_real_distribution<float> reactant_distribution(0, 1);
		uniform_int_distribution<int> index_distribution(0, LATTICE_SIZE - 1);
		uniform_int_distribution<int> pair_direction_distribution(0, 2);
		deque<long long> latest_hydrocarbon_productions;
		long long latest_hydrocarbon_productions_sum = 0;
		long long latest_hydrocarbon_prouctions_average = 0;
		// int base_HC_production;

		hydrocarbon_production = 0;
		hydrocarbon_sizes_seen.clear();
		product_smiles.clear();

		// Main Loop

		while (n_trials < MAX_ITERATIONS_PER_SIMULATION
				&& !(latest_hydrocarbon_productions.size()==EQUILIBRIUM_VERIFICATION_STEPS && abs(latest_hydrocarbon_productions_sum - (latest_hydrocarbon_productions.front() + latest_hydrocarbon_productions.back()) * EQUILIBRIUM_VERIFICATION_STEPS / 2) <= EQUILIBRIUM_VERIFICATION_THRESHOLD))
		{
			// Sleep(10);
			n_trials++;
			Reactant reactant = static_cast<Reactant>(reactant_distribution(generator) > mole_fraction_of_CO);

			if (reactant == Reactant::CO)
			{
				// Choose a site randomly
				int site_row = index_distribution(generator), site_col = index_distribution(generator);

				// Proceed if the site is empty
				if (lattice[site_row][site_col] == "empty")
				{
					n_succesful_trials++;
					lattice[site_row][site_col] = "CO";
					LatticeOperations::process_site(make_pair(site_row, site_col));
				}
				// else
				// {
				//     iter_num--;
				//     continue;
				// }
			}
			else
			{
				array<pair<int, int>, 2> sites;

				// Choose first site randomly
				sites[0] = make_pair(index_distribution(generator), index_distribution(generator));
				// Choose direction of the second site
				Direction direction = static_cast<Direction>(pair_direction_distribution(generator));

				// Choose second site based on the random variable direction
				if (direction == Direction::right)
					sites[1] = make_pair(sites[0].first, (sites[0].second + 1) % LATTICE_SIZE);
				else if (direction == Direction::bottom_right)
					sites[1] = make_pair((sites[0].first + 1) % LATTICE_SIZE, (sites[0].second + (sites[0].first & 1)) % LATTICE_SIZE);
				else if (direction == Direction::bottom_left)
					sites[1] = make_pair((sites[0].first + 1) % LATTICE_SIZE, (sites[0].second + (sites[0].first & 1) - 1 + LATTICE_SIZE) % LATTICE_SIZE);

				bool ok = true;
				for (auto x : sites)
					if (lattice[x.first][x.second] != "empty")
					{
						ok = false;
						break;
					}

				// Proceed if both the sites are empty
				if (ok)
				{
					n_succesful_trials++;
					for (auto x : sites)
						lattice[x.first][x.second] = "H";
					set<pair<int, int>> sites_set;
					sites_set.insert(nearest_neighbours[sites[0].first][sites[0].second].begin(), nearest_neighbours[sites[0].first][sites[0].second].end());
					sites_set.insert(nearest_neighbours[sites[1].first][sites[1].second].begin(), nearest_neighbours[sites[1].first][sites[1].second].end());
					sites_set.erase(sites[0]);
					sites_set.erase(sites[1]);
					vector<pair<int, int>> nearest_sites(sites_set.begin(), sites_set.end());
					shuffle(nearest_sites.begin(), nearest_sites.end(), generator);
					for (pair<int, int> nearest_site : nearest_sites)
						LatticeOperations::process_site(nearest_site);
				}
				// else
				// {
				//     iter_num--;
				//     continue;
				// }
			}

			if(latest_hydrocarbon_productions.size()==EQUILIBRIUM_VERIFICATION_STEPS)
			{
				latest_hydrocarbon_productions_sum -= latest_hydrocarbon_productions.front();
				latest_hydrocarbon_productions.pop_front();
			}
			latest_hydrocarbon_productions_sum += hydrocarbon_production;
			latest_hydrocarbon_productions.push_back(hydrocarbon_production);

			// latest_hydrocarbon_productions_average = double(latest_hydrocarbon_productions_sum) / EQUILIBRIUM_VERIFICATION_STEPS;
			// extremes_average = (latest_hydrocarbon_productions.front() + latest_hydrocarbon_productions.back()) / 2;

			// cerr<<latest_hydrocarbon_productions_average<<' '<<extremes_average<<' '<<latest_hydrocarbon_productions.size()<<' '<<EQUILIBRIUM_VERIFICATION_THRESHOLD<<' '<<abs(latest_hydrocarbon_productions_average - extremes_average)<<'\n';

			if (info_messages && n_trials % PRINT_STEP == 0)
			{
				cerr << "MCS: " << n_trials / LATTICE_POINTS
						<< ", Trials: " << n_trials
						<< ", Successful Trials: " << n_succesful_trials
						<< ", Unsuccessful Trials: " << n_trials - n_succesful_trials;
					//  << '\n';
				cerr<<' '<<latest_hydrocarbon_productions.size()<<' '<<EQUILIBRIUM_VERIFICATION_THRESHOLD<<' '<<latest_hydrocarbon_productions_sum<<' '<<latest_hydrocarbon_productions.front()<<' '<<latest_hydrocarbon_productions.back()<<' '<<(latest_hydrocarbon_productions.front() + latest_hydrocarbon_productions.back()) * EQUILIBRIUM_VERIFICATION_STEPS / 2<<' '<<abs(latest_hydrocarbon_productions_sum - (latest_hydrocarbon_productions.front() + latest_hydrocarbon_productions.back()) * EQUILIBRIUM_VERIFICATION_STEPS / 2)<<'\n';
			}

			if (visualize && n_trials % VISUALIZE_STEP == 0)
			{
				// Visualize current state
				cout << mole_fraction_of_CO<< '\n'
						<< double(n_trials) / LATTICE_POINTS << '\n'
						<< lattice << '\n'
						<< hydrocarbon_production << '\n';
			}

		}

		if(essential_info_messages)
		{
			cerr << "Ran for " << n_trials / LATTICE_POINTS << '.' << n_trials % LATTICE_POINTS << " Monte Carlo Steps\n";
			cerr<<"Hydrocarbon sizes seen: ";
			for(pair<int, int> x:hydrocarbon_sizes_seen)
				cerr<<x.first<<": "<<x.second<<(x.first == hydrocarbon_sizes_seen.rbegin()->first ? "" : ", ");
			cerr<<'\n';
		}

		hydrocarbon_production_outfile << mole_fraction_of_CO << ',' << latest_hydrocarbon_productions.back() - latest_hydrocarbon_productions.front() << ',';
		float f_CO = 0, f_H = 0, f_C = 0, f_E = 0;
		for (auto &x : lattice)
			for (auto &y : x)
				if (y == "CO")
					f_CO++;
				else if (y == "H")
					f_H++;
				else if (y == "C")
					f_C++;
				else if (y == "empty")
					f_E++;
		f_CO /= LATTICE_POINTS;
		f_H /= LATTICE_POINTS;
		f_C /= LATTICE_POINTS;
		f_E /= LATTICE_POINTS;
		hydrocarbon_production_outfile << f_CO << ',' << f_H << ',' << f_C << ',' << f_E << ',' << double(latest_hydrocarbon_productions.back() - latest_hydrocarbon_productions.front()) / EQUILIBRIUM_VERIFICATION_STEPS << ",\n";

		hydrocarbon_sizes_outfile << mole_fraction_of_CO << ' ' << hydrocarbon_sizes_seen.size() << '\n';
		for(auto x : hydrocarbon_sizes_seen)
			hydrocarbon_sizes_outfile << x.first << ' ' << double(x.second) / n_trials << '\n';

		product_smiles_outfile << mole_fraction_of_CO << ' ' << product_smiles.size() << '\n';
		for(string smile : product_smiles)
			product_smiles_outfile << smile << '\n';
	}
	auto end_time = chrono::high_resolution_clock::now();
	chrono::high_resolution_clock::duration time_duration = end_time - start_time;

	auto h = chrono::duration_cast<chrono::hours>(time_duration); time_duration -= h;
	auto m = chrono::duration_cast<chrono::minutes>(time_duration); time_duration -= m;
	auto s = chrono::duration_cast<chrono::seconds>(time_duration); time_duration -= s;
	auto ms = chrono::duration_cast<chrono::milliseconds>(time_duration);

	cerr<<"Program finished.\n";
	cerr<<"Execution Time: "<<h.count()<<" hours "<<m.count()<<" minutes "<<s.count()<<" seconds "<<ms.count()<<" milliseconds";
}