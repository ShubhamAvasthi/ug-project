// Recommended Compilation Command:
// g++ -O3 -march=native

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <random>
#include <unordered_set>
#include <vector>
using namespace std;

const float BOX_SIZE = 20.0;
const float PARTICLE_DIAMETER = 1.0;
const float POTENTIAL_CUTOFF = 1.1;
// const float MIN_DENSITY = 0.1;
const float MAX_MOVEMENT_MAGNITUDE = PARTICLE_DIAMETER / 2;
// r = Rc + 2*n*d
const float VERLET_RADIUS = 4.0;
const float MAX_VERLET_SAFE_MOVEMENT = (VERLET_RADIUS - POTENTIAL_CUTOFF) / 2;
// const int MIN_VERLET_REUSES = (VERLET_RADIUS - POTENTIAL_CUTOFF) / (2 * MAX_MOVEMENT_MAGNITUDE);
const int COORDINATE_PRECISION = 5;
const int PRINT_STEP = 30;
const int VISUALIZE_STEP = 30;

const float pi = 3.14159265358979323846;

// // Custom ceil function to return constexpr int value
// // This in turn will help to get n_points during compilation time, without compromising code readability and functionality.
// // It allows us to still express n_points as an expression and not hard code the value.
// // n_points, if evaluated during compilation time allows us to initialize a number of arrays during compile time as well.
// // The above point allows us to replace "std::vector"s with faster but static sized "std::array"s
// constexpr int32_t custom_ceil(float num)
// {
//     return (static_cast<float>(static_cast<int32_t>(num)) == num)
//                ? static_cast<int32_t>(num)
//                : static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
// }

const int MAX_POINTS = 5000;
bool visualize = true;
bool info_messages = true;
bool essential_info_messages = true;
// n*(1/6)*pi*PARTICLE_DIAMETER^3 > MIN_DENSITY*BOX_SIZE^3
int n_points; //custom_ceil((6 / pi) * MIN_DENSITY * (BOX_SIZE / PARTICLE_DIAMETER) * (BOX_SIZE / PARTICLE_DIAMETER) * (BOX_SIZE / PARTICLE_DIAMETER));
array<array<float, 3>, MAX_POINTS> points;
array<int, MAX_POINTS> dsu_par;              // Parent array for DSU, note that the parent will also be the color of a DSU tree
array<int, MAX_POINTS> dsu_size;             // Sizes of DSU trees
array<vector<int>, MAX_POINTS> dsu_children; // Children under a leader, used for cluster movement, contains only adjacent children
vector<int> cluster_leader_indices;
array<vector<int>, MAX_POINTS> neighbour;
array<array<float, 3>, MAX_POINTS> previous_verlet_points;
int n_clusters = n_points;

// A useful overload
ostream &operator<<(ostream &out, array<float, 3> &point)
{
    for (auto &x : point)
        out << x << ' ';
    return out;
}

// Calculates Euclidean distance between two points
float dist(array<float, 3> p1, array<float, 3> p2)
{
    array<float, 3> dist;
    for(int i=0;i<3;i++)
        dist[i] = abs(p1[i] - p2[i]);
    return pow(dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2],0.5);
}

// Calculates minimum Euclidean distance amongst two points
// considering all mirrored points in the nearby boxes
float min_dist(array<float, 3> p1, array<float, 3> p2)
{
    array<float, 3> dist;
    for(int i=0;i<3;i++)
    {
        dist[i] = abs(p1[i] - p2[i]);
        dist[i] = min(dist[i], BOX_SIZE - dist[i]);
    }
    return pow(dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2],0.5);
}

// Calculates root of the DSU tree with path compresion
int dsu_root(int x)
{
    if (dsu_par[x] != x)
    {
        dsu_par[x] = dsu_root(dsu_par[x]);
        // dsu_children[dsu_par[x]].push_back(x);
    }
    return dsu_par[x];
}

// Performs union by size of two elements and therefore the DSU trees
void dsu_union(int x, int y)
{
    int r1 = dsu_root(x), r2 = dsu_root(y);
    if (r1 == r2)
        return;

    // Ensure size(r1) > size(r2)
    if (dsu_size[r2] > dsu_size[r1])
        swap(r1, r2);

    dsu_size[r1] += dsu_size[r2];
    dsu_par[r2] = r1;
    dsu_children[r1].push_back(r2);
    cluster_leader_indices.erase(find(cluster_leader_indices.begin(), cluster_leader_indices.end(), r2));
    n_clusters--;
}

// Returns all children (not just adjacent) with the passed node as the leader
void all_children(vector<int> &children, int x, int topmost_parent = -1)
{
    children.push_back(x);
    if (topmost_parent == -1)
        topmost_parent = x;
    dsu_par[x] = topmost_parent;
    for (int c : dsu_children[x])
        all_children(children, c, topmost_parent);
}

int main(int argc, char *argv[])
{
    ios::sync_with_stdio(0);

    // TODO: Use getopt or something similar and give option for single run
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

    ofstream outfile ("radii_of_gyration.csv");
    outfile<<"dist";
    for(int i=0;i<=5*int(BOX_SIZE);i++)
        outfile<<','<<0.1*i;
    outfile<<"\n\nn_points\n";

    if (visualize)
    {
        cout << fixed;
        cout.precision(COORDINATE_PRECISION);
        cout<<BOX_SIZE<<'\n';
    }

    for (n_points = 100; n_points <= 5000; n_points+=100)
    {
        n_clusters = n_points;

        // Initialize the random generator with a time based seed (namely the current time)
        default_random_engine generator(
            chrono::system_clock::now().time_since_epoch().count());
        uniform_real_distribution<float> distribution(0.0, BOX_SIZE);

        if (essential_info_messages)
            cerr << "n_points: " << n_points << '\n';

        // Initialize dsu_size array
        dsu_size.fill(1);

        // Initialize parents for DSU
        for (int i = 0; i < n_points; i++)
            dsu_par[i] = i;

        // Initialize children for DSU
        for (int i = 0; i < n_points; i++)
            dsu_children[i].clear();

        // Initialize cluster leaders
        cluster_leader_indices.clear();
        for (int i = 0; i < n_points; i++)
            cluster_leader_indices.push_back(i);

        // Random Initialization
        for (int i = 0; i < n_points; i++)
        {
            array<float, 3> point;
            for (float &x : point)
                x = distribution(generator);
            bool ok = true;
            for (int j = 0; j < i; j++)
            {
                // Reject configuration if there is an overlap
                if (min_dist(point, points[j]) < PARTICLE_DIAMETER)
                {
                    ok = false;
                    break;
                }
            }
            if (!ok)
            {
                i--;
                continue;
            }
            points[i] = point;
        }

        //Calculation of Initial Energy (in units of epsilon) and Initialization of DSU
        int energy = 0;
        for (int i = 0; i < n_points; i++)
            for (int j = i + 1; j < n_points; j++)
                if (min_dist(points[i], points[j]) < POTENTIAL_CUTOFF)
                {
                    energy--;
                    dsu_union(i, j);
                }

        // Variables for Main Loop
        int n_trials = 0, n_updates = 0, n_rejections = 0, n_overlaps = 0;
        uniform_real_distribution<float> movement_distribution(-MAX_MOVEMENT_MAGNITUDE, MAX_MOVEMENT_MAGNITUDE);

        bool verlet_list_updated = false;

        // Initialize Previous Verlet Points
        for (int i = 0; i < n_points; i++)
            previous_verlet_points[i] = points[i];

        // Create Verlet List
        for (int i = 0; i < n_points; i++)
            for (int j = 0; j < n_points; j++)
            {
                if (i == j)
                    continue;
                if (min_dist(points[i], points[j]) < VERLET_RADIUS)
                    neighbour[i].push_back(j);
            }

        if(visualize)
        {
            // Visualize initial state
            cout << n_points << '\n'
                    << n_trials << ' ' << energy << '\n';
            for (int i = 0; i < n_points; i++)
                cout << points[i] << '\n';
            for (int i = 0; i < n_points; i++)
                cout << dsu_root(i) << ' ';
            cout << '\n';
        }

        // Main Loop
        while (true)
        {
            // Choose a cluster randomly
            int leader_index = cluster_leader_indices[uniform_int_distribution<int>(0, n_clusters - 1)(generator)];

            // Choose a random movement
            array<float, 3> movement;
            for (float &x : movement)
                x = movement_distribution(generator);

            vector<int> cluster_indices;
            all_children(cluster_indices, leader_index);

            // This vector has to be sorted to ignore its elements for overlap
            sort(cluster_indices.begin(), cluster_indices.end());

            bool overlapped = false;
            //Apply movement to the chosen cluster
            //and calculate its new position
            vector<array<float, 3>> new_pos(cluster_indices.size());
            for (int i = 0; i < cluster_indices.size(); i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    new_pos[i][j] = points[cluster_indices[i]][j] + movement[j];
                    if (new_pos[i][j] < 0)
                        new_pos[i][j] += BOX_SIZE;
                    else if (new_pos[i][j] > BOX_SIZE)
                        new_pos[i][j] -= BOX_SIZE;
                }

                // Check for any overlaps
                bool ok = true;
                int next_ignore_index = 0;
                for (int j : neighbour[cluster_indices[i]])
                {
                    // Ensure j is not one of cluster_indices
                    while (j > cluster_indices[next_ignore_index] && next_ignore_index != cluster_indices.size() - 1)
                        next_ignore_index++;
                    if (j == cluster_indices[next_ignore_index])
                    {
                        if (next_ignore_index != cluster_indices.size() - 1)
                            next_ignore_index++;
                        continue;
                    }

                    // Reject configuration if there is an overlap
                    if (min_dist(new_pos[i], points[j]) < PARTICLE_DIAMETER)
                    {
                        ok = false;
                        break;
                    }
                }
                if (!ok)
                {
                    overlapped = true;
                    break;
                }
            }

            if (overlapped)
            {
                n_overlaps++;
                continue;
            }

            // Update trial count
            n_trials++;

            // Calculate new energy
            int old_energy_contrib = 0, new_energy_contrib = 0;
            for (int i = 0; i < cluster_indices.size(); i++)
            {
                int next_ignore_index = 0;
                for (int j : neighbour[cluster_indices[i]])
                {
                    // Ensure j is not one of cluster_indices
                    while (j > cluster_indices[next_ignore_index] && next_ignore_index != cluster_indices.size() - 1)
                        next_ignore_index++;
                    if (j == cluster_indices[next_ignore_index])
                    {
                        if (next_ignore_index != cluster_indices.size() - 1)
                            next_ignore_index++;
                        continue;
                    }

                    if (min_dist(points[cluster_indices[i]], points[j]) < POTENTIAL_CUTOFF)
                        old_energy_contrib--;
                    if (min_dist(new_pos[i], points[j]) < POTENTIAL_CUTOFF)
                        new_energy_contrib--;
                }
            }

            if (new_energy_contrib > old_energy_contrib)
                // Reject configuration if new energy contribution is
                // greater than the old energy contribution
                n_rejections++;
            else
            {
                // Update energy
                energy += new_energy_contrib - old_energy_contrib;

                //Apply the movement
                for (int i = 0; i < cluster_indices.size(); i++)
                    points[cluster_indices[i]] = new_pos[i];

                bool update_verlet_list = false;
                // Update verlet list if maximum safe movement can be reached in the next step
                for (int i : cluster_indices)
                    if (min_dist(points[i], previous_verlet_points[i]) + MAX_MOVEMENT_MAGNITUDE > MAX_VERLET_SAFE_MOVEMENT)
                        update_verlet_list = true; // Verlet list will be updated

                if (update_verlet_list)
                {
                    verlet_list_updated = true;
                    for (int i = 0; i < n_points; i++)
                    {
                        neighbour[i].clear();
                        for (int j = 0; j < n_points; j++)
                        {
                            if (i == j)
                                continue;
                            if (min_dist(points[i], points[j]) < VERLET_RADIUS)
                                neighbour[i].push_back(j);
                        }
                    }
                    previous_verlet_points = points;
                }
                
                // Update DSU data structure
                for (int i : cluster_indices)
                    for (int j : neighbour[i])
                        if (min_dist(points[i], points[j]) < POTENTIAL_CUTOFF)
                            dsu_union(i, j);

                n_updates++;

                // TODO: Asynchronise Python script
                //cout<<n_points<<'\n'<<float(n_trials)/n_points<<' '<<energy<<'\n'<<points;
            }

            // cout<<n_points<<'\n'<<float(n_trials)/n_points<<' '<<energy<<'\n'<<points;
            if (info_messages && (n_clusters == 1 || n_trials % PRINT_STEP == 0))
            {
                cerr //<< "MCS: " << n_trials / n_points
                    << "Trials: " << n_trials
                    << ", Clusters: " << n_clusters
                    << ", Updates: " << n_updates
                    << ", Rejections: " << n_rejections
                    << ", Overlaps: " << n_overlaps;

                if (verlet_list_updated)
                    cerr << " VLU"; // Verlet List Updated

                cerr << '\n';

                verlet_list_updated = false;
            }

            if (visualize && (n_clusters == 1 || n_trials % VISUALIZE_STEP == 0))
            {
                cout << n_points << '\n'
                     << n_trials << ' ' << energy << '\n';
                for (int i = 0; i < n_points; i++)
                    cout << points[i] << '\n';
                // Print root index for colours
                for (int i = 0; i < n_points; i++)
                    cout << dsu_root(i) << ' ';
                cout << '\n';
            }

            if (n_clusters == 1)
                break;
        }
        if (info_messages)
            cerr << "All particles connected to a single cluster.\n";

        // long double r_square = 0;
        // array<float, 3> sum;
        // for(int i=0;i<3;i++)
        //     sum[i] = 0;
        // for(int i=0;i<n_points;i++)
        //     for(int j=0;j<3;j++)
        //         sum[j] += points[i][j];
        // for(int i=0;i<3;i++)
        //     sum[i] /= n_points;
        // for(int i=0;i<n_points;i++)
        // {
        //     array<float, 3> dist;
        //     for(int j=0;j<3;j++)
        //     dist[j] = abs(points[i][j] - sum[j]);
        //     r_square += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
        // }
        // r_square/=n_points;
        // outfile<<n_points<<','<<r_square<<'\n';
        outfile<<n_points;
        array<float,3> centre = {BOX_SIZE/2,BOX_SIZE/2,BOX_SIZE/2};
        for(float cur_dist=0;cur_dist<=BOX_SIZE/2 + 0.001;cur_dist+=0.1)
        {
            int n=0;
            for(int i=0;i<n_points;i++)
                n+=(dist(points[i],centre)<cur_dist);
            outfile<<','<<n;
        }
        outfile<<'\n';
    }
}