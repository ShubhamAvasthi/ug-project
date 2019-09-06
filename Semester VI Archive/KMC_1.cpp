// Recommended Compilation Command:
// g++ -O3 -march=native

#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <windows.h>
using namespace std;

const int FLOAT_OUTPUT_PRECISION = 5;
const int LATTICE_SIZE = 100;                           // Has to be an even integer
const int LATTICE_POINTS = LATTICE_SIZE * LATTICE_SIZE; // Just for readability
const int PRINT_STEP = LATTICE_POINTS;
const int VISUALIZE_STEP = 200 * LATTICE_POINTS;
const int STEPS_BEFORE_EQUILIBRIUM = 10000 * LATTICE_POINTS;
const int STEPS_AFTER_EQUILIBRIUM = 10000 * LATTICE_POINTS;
const float MF_STEP = 0.05;

bool visualize = true, info_messages = true, essential_info_messages = true;
default_random_engine generator;
int hydrocarbon_production = 0;
array<array<string, LATTICE_SIZE>, LATTICE_SIZE> lattice;
array<array<array<pair<int, int>, 6>, LATTICE_SIZE>, LATTICE_SIZE> nearest_neighbours;
array<array<set<pair<int, int>>, LATTICE_SIZE>, LATTICE_SIZE> bonded_H_atoms;
array<array<pair<int, int>, LATTICE_SIZE>, LATTICE_SIZE> is_connected_to;

// Some useful overloads
ostream &operator<<(ostream &out, array<array<string, LATTICE_SIZE>, LATTICE_SIZE> &lattice)
{
    for (auto &x : lattice)
        for (auto &y : x)
            out << y << ' ';
    return out;
}

// Initiates reaction on a CO site if possible
void initiate_reaction_with_CO_if_possible(pair<int, int> site)
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
                }
                is_connected_to[nearest_sites_with_H[i].first][nearest_sites_with_H[i].second] = make_pair(-1, -1);
            }
        }
    }
}

// Initiates reaction on a C, CH, CH2, CH3 site if possible
void initiate_reaction_with_CHx_if_possible(pair<int, int> site)
{
    // If one or more H atoms are nearest neighbours to a CHx site, increase x and possibly desorb CH4
    if (lattice[site.first][site.second] == "C")
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

        // Remove CH4 if possible
        if (bonded_H_atoms[site.first][site.second].size() == 4)
        {
            hydrocarbon_production++;
            for (auto x : bonded_H_atoms[site.first][site.second])
            {
                lattice[x.first][x.second] = "empty";
                is_connected_to[x.first][x.second] = make_pair(-1, -1);
            }
            lattice[site.first][site.second] = "empty";
            bonded_H_atoms[site.first][site.second].clear();
        }
    }
}

// Initiates desorption on a C, CH, CH2, CH3 site if possible
void initiate_desorption_on_CHx_if_possible(pair<int, int> site)
{
    if (lattice[site.first][site.second] != "C")
        return;

    vector<pair<int, int>> nearest_sites_with_C;
    for (auto neighbour_site : nearest_neighbours[site.first][site.second])
        if (lattice[neighbour_site.first][neighbour_site.second] == "C" && bonded_H_atoms[neighbour_site.first][neighbour_site.second].size() == bonded_H_atoms[site.first][site.second].size())
            nearest_sites_with_C.push_back(neighbour_site);
    shuffle(nearest_sites_with_C.begin(), nearest_sites_with_C.end(), generator);

    if (nearest_sites_with_C.empty())
        return;

    for (auto H_site : bonded_H_atoms[site.first][site.second])
    {
        lattice[H_site.first][H_site.second] = "empty";
        is_connected_to[H_site.first][H_site.second] = make_pair(-1, -1);
    }
    for (auto H_site : bonded_H_atoms[nearest_sites_with_C[0].first][nearest_sites_with_C[0].second])
    {
        lattice[H_site.first][H_site.second] = "empty";
        is_connected_to[H_site.first][H_site.second] = make_pair(-1, -1);
    }

    lattice[site.first][site.second] = "empty";
    lattice[nearest_sites_with_C[0].first][nearest_sites_with_C[0].second] = "empty";
    bonded_H_atoms[site.first][site.second].clear();
    bonded_H_atoms[nearest_sites_with_C[0].first][nearest_sites_with_C[0].second].clear();
    hydrocarbon_production++;
}

int main(int argc, char *argv[])
{
    ios_base::sync_with_stdio(0);

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

    ofstream outfile("out.csv");
    outfile << "Mole Fraction of CO,Hydrocarbon production,Coverage fraction of CO,Coverage fraction of H,Coverage fraction of C,Coverage fraction of empty sites\n";

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

    for (float mole_fraction_of_CO = MF_STEP; mole_fraction_of_CO <= 1.001; mole_fraction_of_CO += MF_STEP)
    {
        if (essential_info_messages)
            cerr << "Mole fraction of CO: " << mole_fraction_of_CO << '\n';

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

        // Initialize the lattice
        for (auto &row : lattice)
            for (auto &element : row)
                element = "empty";

        for (auto &row : is_connected_to)
            for (auto &element : row)
                element = make_pair(-1, -1);

        // Variables for Main Loop
        hydrocarbon_production = 0;
        int n_trials = 0, n_succesful_trials = 0;
        uniform_real_distribution<float> reactant_distribution(0, 1);
        uniform_int_distribution<int> index_distribution(0, LATTICE_SIZE - 1);
        uniform_int_distribution<int> pair_direction_distribution(0, 2);
        int base_HC_production;

        // Main Loop

        for (int equilibrium_approached = 0; equilibrium_approached < 2; equilibrium_approached++)
        {
            int iter_steps = (equilibrium_approached ? STEPS_AFTER_EQUILIBRIUM : STEPS_BEFORE_EQUILIBRIUM);

            // When the iteration with iter_num = 0 ends, equilibrium is supposed to be achieved
            // All required values are calculated in the iteration when iter_num = 1
            for (int iter_num = 0; iter_num < iter_steps; iter_num++)
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
                        initiate_reaction_with_CO_if_possible(make_pair(site_row, site_col));
                        initiate_reaction_with_CHx_if_possible(make_pair(site_row, site_col));
                        initiate_desorption_on_CHx_if_possible(make_pair(site_row, site_col));
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

                    bool ok = 1;
                    for (auto x : sites)
                        if (lattice[x.first][x.second] != "empty")
                        {
                            ok = 0;
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
                            initiate_reaction_with_CO_if_possible(nearest_site);
                        shuffle(nearest_sites.begin(), nearest_sites.end(), generator);
                        for (pair<int, int> nearest_site : nearest_sites)
                        {
                            initiate_reaction_with_CHx_if_possible(nearest_site);
                            initiate_desorption_on_CHx_if_possible(nearest_site);
                        }
                    }
                    // else
                    // {
                    //     iter_num--;
                    //     continue;
                    // }
                }

#define count_trials n_trials
                if (info_messages && count_trials % PRINT_STEP == 0)
                {
                    cerr << "MCS: " << count_trials / LATTICE_POINTS
                         << ", Trials: " << n_trials
                         << ", Successful Trials: " << n_succesful_trials
                         << ", Unsuccessful Trials: " << n_trials - n_succesful_trials
                         << '\n';
                }

                if (visualize && count_trials % VISUALIZE_STEP == 0)
                {
                    // Visualize current state
                    cout << count_trials / LATTICE_POINTS << '\n'
                         << lattice << '\n'
                         << hydrocarbon_production << '\n';
                }
#undef count_trials
            }
            
            if(!equilibrium_approached)
                base_HC_production = hydrocarbon_production;

            if(equilibrium_approached)
            {
                outfile << mole_fraction_of_CO << ',' << hydrocarbon_production - base_HC_production << ',';
                float f_CO=0, f_H=0, f_C=0, f_E=0;
                for(auto &x:lattice)
                    for(auto &y: x)
                        if(y == "CO")
                            f_CO++;
                        else if(y == "H")
                            f_H++;
                        else if(y == "C")
                            f_C++;
                        else if(y == "empty")
                            f_E++;
                f_CO/=LATTICE_POINTS;
                f_H/=LATTICE_POINTS;
                f_C/=LATTICE_POINTS;
                f_E/=LATTICE_POINTS;
                outfile<<f_CO<<','<<f_H<<','<<f_C<<','<<f_E<<",\n";
            }
        }
    }
}