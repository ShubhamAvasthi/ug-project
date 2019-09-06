// Recommended Compilation Command:
// g++ -O3 -march=native <filename>.cpp -o <filename>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <windows.h>
using namespace std;

const int FLOAT_OUTPUT_PRECISION = 5;
const int LATTICE_SIZE = 100;                           // Has to be an even integer
const long long LATTICE_POINTS = LATTICE_SIZE * LATTICE_SIZE; // Just for readability
const long long PRINT_STEP = 50 * LATTICE_POINTS;
const long long VISUALIZE_STEP = 10 * LATTICE_POINTS;
const float MF_STEP = 0.001; //0.05

// Equilibrium will be assumed to have arrived, when the difference between the average of the hydrocarbon production in the last EQUILIBRIUM_VERIFICATION_STEPS and
// the average of the latest and the oldest steps is less than EQUILIBRIUM_VERIFICATION_THRESHOLD
const int EQUILIBRIUM_VERIFICATION_STEPS = 500000;
const int EQUILIBRIUM_VERIFICATION_THRESHOLD = 0;
// const int STEPS_BEFORE_EQUILIBRIUM = 100/*000*/ * LATTICE_POINTS;
// const int STEPS_AFTER_EQUILIBRIUM = 100/*000*/ * LATTICE_POINTS;

bool visualize = true, info_messages = true, essential_info_messages = true;
default_random_engine generator;
long long hydrocarbon_production = 0;
array<array<string, LATTICE_SIZE>, LATTICE_SIZE> lattice;
array<array<array<pair<int, int>, 6>, LATTICE_SIZE>, LATTICE_SIZE> nearest_neighbours;
array<array<set<pair<int, int>>, LATTICE_SIZE>, LATTICE_SIZE> bonded_H_atoms;
array<array<pair<int, int>, LATTICE_SIZE>, LATTICE_SIZE> is_connected_to;

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

// // Search chains
// void search_chains(pair<int, int> site, int incoming_valency, vector<vector<pair<int, int>>> &ans, vector<vector<int>> &vis, vector<pair<int, int>> &chain)
// {
//     if (lattice[site.first][site.second] != "C" || vis[site.first][site.second] || incoming_valency == 4 || incoming_valency + bonded_H_atoms[site.first][site.second].size() > 4)
//         return;

//     // cerr<<site.first<<','<<site.second<<'\n';
//     chain.push_back(site);
//     vis[site.first][site.second] = true;

//     if (incoming_valency + bonded_H_atoms[site.first][site.second].size() == 4)
//         ans.push_back(chain);
//     else
//         for (auto neighbour_site : nearest_neighbours[site.first][site.second])
//         {
//             // cerr<<neighbour_site.first<<','<<neighbour_site.second<<' '<<4 - incoming_valency - bonded_H_atoms[site.first][site.second].size()<<'\n';
//             search_chains(neighbour_site, 4 - incoming_valency - bonded_H_atoms[site.first][site.second].size(), ans, vis, chain);
//         }
//     chain.pop_back();
//     vis[site.first][site.second] = false;
// }

// // Initiates desorption on a C, CH, CH2, CH3 site if possible
// void initiate_desorption_on_CHx_if_possible(pair<int, int> site)
// {
//     vector<pair<int, int>> chain;
//     vector<vector<pair<int, int>>> ans;
//     vector<vector<int>> vis(LATTICE_SIZE, vector<int>(LATTICE_SIZE));
//     search_chains(site, 0, ans, vis, chain);

//     if (ans.empty())
//         return;

//     hydrocarbon_production++;
//     int index = uniform_int_distribution<int>(0, ans.size())(generator);
//     for (auto C_site : ans[index])
//     {
//         for (auto H_site : bonded_H_atoms[C_site.first][C_site.second])
//         {
//             lattice[H_site.first][H_site.second] = "empty";
//             is_connected_to[H_site.first][H_site.second] = make_pair(-1, -1);
//         }
//     }
//     for (auto C_site : ans[index])
//     {
//         lattice[C_site.first][C_site.second] = "empty";
//         bonded_H_atoms[C_site.first][C_site.second].clear();
//     }
// }

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
    auto start_time = chrono::high_resolution_clock::now();

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
    outfile << "Mole Fraction of CO,Hydrocarbon production,Coverage fraction of CO,Coverage fraction of H,Coverage fraction of C,Coverage fraction of empty sites,Average hydrocarbon production\n";

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

        // Main Loop

        // for (int equilibrium_approached = 0; equilibrium_approached < 2; equilibrium_approached++)
        {
            // int iter_steps = (equilibrium_approached ? STEPS_AFTER_EQUILIBRIUM : STEPS_BEFORE_EQUILIBRIUM);

            // When the iteration with iter_num = 0 ends, equilibrium is supposed to be achieved
            // All required values are calculated in the iteration when iter_num = 1
            // for (int iter_num = 0; iter_num < iter_steps; iter_num++)
            while(!(latest_hydrocarbon_productions.size()==EQUILIBRIUM_VERIFICATION_STEPS && abs(latest_hydrocarbon_productions_sum - (latest_hydrocarbon_productions.front() + latest_hydrocarbon_productions.back()) * EQUILIBRIUM_VERIFICATION_STEPS / 2) <= EQUILIBRIUM_VERIFICATION_THRESHOLD))
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

#define count_trials n_trials
                if (info_messages && count_trials % PRINT_STEP == 0)
                {
                    cerr << "MCS: " << count_trials / LATTICE_POINTS
                         << ", Trials: " << n_trials
                         << ", Successful Trials: " << n_succesful_trials
                         << ", Unsuccessful Trials: " << n_trials - n_succesful_trials;
                        //  << '\n';
                    cerr<<' '<<latest_hydrocarbon_productions.size()<<' '<<EQUILIBRIUM_VERIFICATION_THRESHOLD<<' '<<latest_hydrocarbon_productions_sum<<' '<<latest_hydrocarbon_productions.front()<<' '<<latest_hydrocarbon_productions.back()<<' '<<(latest_hydrocarbon_productions.front() + latest_hydrocarbon_productions.back()) * EQUILIBRIUM_VERIFICATION_STEPS / 2<<' '<<abs(latest_hydrocarbon_productions_sum - (latest_hydrocarbon_productions.front() + latest_hydrocarbon_productions.back()) * EQUILIBRIUM_VERIFICATION_STEPS / 2)<<'\n';
                }

                if (visualize && count_trials % VISUALIZE_STEP == 0)
                {
                    // Visualize current state
                    cout << mole_fraction_of_CO<< '\n'
                         << double(count_trials) / LATTICE_POINTS << '\n'
                         << lattice << '\n'
                         << hydrocarbon_production << '\n';
                }
#undef count_trials

            }

            // if (!equilibrium_approached)
            //     base_HC_production = hydrocarbon_production;

            // if (equilibrium_approached)
            // {
                outfile << mole_fraction_of_CO << ',' << latest_hydrocarbon_productions.back() - latest_hydrocarbon_productions.front() << ',';// hydrocarbon_production - base_HC_production << ',';
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
                outfile << f_CO << ',' << f_H << ',' << f_C << ',' << f_E << ',' << double(latest_hydrocarbon_productions.back() - latest_hydrocarbon_productions.front()) / EQUILIBRIUM_VERIFICATION_STEPS << ",\n";
            // }
        }
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