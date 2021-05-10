// Recommended Compilation Command:
// g++ -O3 -march=native
// Note to self: Check -g in IORun

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>
using namespace std;

const float BOX_SIZE = 10;
const float PARTICLE_DIAMETER = 1.0;
const float POTENTIAL_CUTOFF = 1.1;
const float MIN_DENSITY = 0.1;
const float MAX_MOVEMENT_MAGNITUDE = PARTICLE_DIAMETER / 2;
// r = Rc + 2*n*d
const float VERLET_RADIUS = 3.0;
const float MAX_VERLET_SAFE_MOVEMENT = (VERLET_RADIUS - POTENTIAL_CUTOFF) / 2;
// const int MIN_VERLET_REUSES = (VERLET_RADIUS - POTENTIAL_CUTOFF) / (2 * MAX_MOVEMENT_MAGNITUDE);
const int COORDINATE_PRECISION = 5;
const int PRINT_MCS_STEP = 1;
const int VISUALIZE_MCS_STEP = 10;

const float pi = 3.14159265358979323846;

// Custom ceil function to return constexpr int value
// This in turn will help to get n_points during compilation time, without compromising code readability and functionality.
// It allows us to still express n_points as an expression and not hard code the value.
// n_points, if evaluated during compilation time allows us to initialize a number of arrays during compile time as well.
// The above point allows us to replace "std::vector"s with faster but static sized "std::array"s
constexpr int32_t custom_ceil(float num)
{
    return (static_cast<float>(static_cast<int32_t>(num)) == num)
               ? static_cast<int32_t>(num)
               : static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
}

// n*(1/6)*pi*PARTICLE_DIAMETER^3 > MIN_DENSITY*BOX_SIZE^3
const int n_points = custom_ceil((6 / pi) * MIN_DENSITY * (BOX_SIZE / PARTICLE_DIAMETER) * (BOX_SIZE / PARTICLE_DIAMETER) * (BOX_SIZE / PARTICLE_DIAMETER));
array<array<float, 3>, n_points> points;
array<int, n_points> dsu_par;      // Parent array for DSU, note that the parent will also be the color of a DSU tree
array<int, n_points> dsu_size;     // Sizes of DSU trees
array<vector<int>, n_points> dsu_children; // Children under a leader, used for cluster movement
unordered_set<int> cluster_set;
int n_clusters = n_points;

// Some useful overloads
ostream &operator<<(ostream &out, array<float, 3> &point)
{
    for (auto &x : point)
        out << x << ' ';
    return out;
}

ostream &operator<<(ostream &out, array<array<float, 3>, n_points> &points)
{
    for (auto &x : points)
        out << x << '\n';
    return out;
}

// Utility Functions

// Calculates minimum Euclidean distance amongst two points
// considering all mirrored points in the nearby boxes
float min_dist(array<float, 3> p1, array<float, 3> p2)
{
    array<float, 3> dist;
    dist[0] = abs(p1[0] - p2[0]);
    dist[0] = min(dist[0], BOX_SIZE - dist[0]);
    dist[1] = abs(p1[1] - p2[1]);
    dist[1] = min(dist[1], BOX_SIZE - dist[1]);
    dist[2] = abs(p1[2] - p2[2]);
    dist[2] = min(dist[2], BOX_SIZE - dist[2]);
    return pow(dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2], 1/3);
}

// Calculates root of the DSU tree with path compresion
int dsu_root(int x)
{                        
    if (dsu_par[x] != x)
    {
        dsu_par[x] = dsu_root(dsu_par[x]);
        dsu_children[dsu_par[x]].push_back(x);
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
        swap(r1,r2);

    dsu_size[r1] += dsu_size[r2];
    dsu_par[r2] = r1;
    dsu_children[r1].push_back(r2);
    cluster_set.erase(r2);
    n_clusters--;
}

int main()
{
    ios::sync_with_stdio(0);
    //cerr.setstate(ios::failbit);
    //cout.setstate(ios::failbit);
    cout << fixed;
    cout.precision(COORDINATE_PRECISION);

    // Initialize the random generator with a time based seed (namely the current time)
    default_random_engine generator(
        chrono::system_clock::now().time_since_epoch().count());
    uniform_real_distribution<float> distribution(0.0, 10.0);

    cerr << "n_points: " << n_points << '\n';

    // Initialize dsu_size array
    dsu_size.fill(1);

    // Initialize parents for DSU
    for (int i = 0; i < n_points; i++)
        dsu_par[i] = i;

    // Initialize cluster set
    for (int i = 0; i < n_points; i++)
        cluster_set.insert(i);

    // Initialize children of each cluster
    for (int i = 0; i < n_points; i++)
        dsu_children[i].push_back(i);

    // Random Initialization
    for (int i = 0; i < n_points; i++)
    {
        array<float, 3> point;
        for (float &x : point)
            x = distribution(generator);
        // point[2] = 0;
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
    {
        for (int j = i + 1; j < n_points; j++)
            if (min_dist(points[i], points[j]) < POTENTIAL_CUTOFF)
                energy--, dsu_union(i, j);
    }

    // Variables for Main Loop
    int n_trials = 0, n_updates = 0, n_rejections = 0, n_overlaps = 0;
    uniform_int_distribution<int> index_distribution(0, n_points - 1);
    uniform_real_distribution<float> movement_distribution(-MAX_MOVEMENT_MAGNITUDE, MAX_MOVEMENT_MAGNITUDE);

    array<vector<int>, n_points> neighbour;
    array<array<float, 3>, n_points> previous_verlet_points = points;

    bool verlet_list_updated = 0;

    // Create Verlet List
    for (int i = 0; i < n_points; i++)
        for (int j = 0; j < n_points; j++)
        {
            if (i == j)
                continue;
            if (min_dist(points[i], points[j]) < VERLET_RADIUS)
                neighbour[i].push_back(j);
        }

    // Visualize initial state
    cout << n_points << '\n'
         << float(n_trials) / n_points << ' ' << energy << '\n'
         << points;
    for (int i = 0; i < n_points; i++)
        cout << dsu_root(i) << ' ';
    cout << '\n';

    // Main Loop
    while (true)
    {
        // Choose a particle randomly
        int index = index_distribution(generator);

        // Choose a random movement
        array<float, 3> movement;
        for (float &x : movement)
            x = movement_distribution(generator);
        // movement[2] = 0;

        //Apply movement to the chosen particle
        //and calculate its new position
        array<float, 3> new_pos;
        for (int i = 0; i < 3; i++)
        {
            new_pos[i] = points[index][i] + movement[i];
            if (new_pos[i] < 0)
                new_pos[i] += BOX_SIZE;
            else if (new_pos[i] > BOX_SIZE)
                new_pos[i] -= BOX_SIZE;
        }

        // Check for any overlaps
        bool ok = true;
        for (int i = 0; i < n_points; i++)
        {
            if (i == index)
                continue;
            // Reject configuration if there is an overlap
            if (min_dist(new_pos, points[i]) < PARTICLE_DIAMETER)
            {
                ok = false;
                break;
            }
        }
        if (!ok)
        {
            n_overlaps++;
            continue;
        }

        // Update trial count
        n_trials++;

        // Calculate new energy
        int old_energy_contrib = 0, new_energy_contrib = 0;
        for (int i : neighbour[index])
        {
            if (min_dist(points[index], points[i]) < POTENTIAL_CUTOFF)
                old_energy_contrib--;
            if (min_dist(new_pos, points[i]) < POTENTIAL_CUTOFF)
                new_energy_contrib--;
        }

        if (new_energy_contrib > old_energy_contrib)
            // Reject configuration if new energy contribution is greater than
            // or equal to the old energy contribution
            n_rejections++;
        else
        {
            // Update energy
            energy += new_energy_contrib - old_energy_contrib;

            //Apply the movement
            points[index] = new_pos;

            // Update verlet list if maximum safe movement can be reached in the next step
            if (min_dist(points[index], previous_verlet_points[index]) + MAX_MOVEMENT_MAGNITUDE > MAX_VERLET_SAFE_MOVEMENT)
            {
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
                verlet_list_updated = 1;
            }

            // Update DSU data structure
            int old_energy_contrib = 0, new_energy_contrib = 0;
            for (int i : neighbour[index])
                if (min_dist(points[index], points[i]) < POTENTIAL_CUTOFF)
                    dsu_union(index, i);

            n_updates++;

            // TODO: Asynchronise Python script
            //cout<<n_points<<'\n'<<float(n_trials)/n_points<<' '<<energy<<'\n'<<points;
        }

        // cout<<n_points<<'\n'<<float(n_trials)/n_points<<' '<<energy<<'\n'<<points;
        if (n_clusters == 1 || n_trials % (n_points * PRINT_MCS_STEP) == 0)
        {
            cerr << "MCS: " << n_trials / n_points
                 << ", Trials: " << n_trials
                 << ", Clusters: " << n_clusters
                 << ", Updates: " << n_updates
                 << ", Rejections: " << n_rejections
                 << ", Overlaps: " << n_overlaps;

            if (verlet_list_updated)
                cerr << " VLU"; // Verlet List Updated

            cerr << '\n';

            verlet_list_updated = 0;
        }

        if (n_clusters == 1 || n_trials % (n_points * VISUALIZE_MCS_STEP) == 0)
        {
            cout << n_points << '\n'
                 << float(n_trials) / n_points << ' ' << energy << '\n'
                 << points;
            // Print root index for colours
            for (int i = 0; i < n_points; i++)
                cout << dsu_root(i) << ' ';
            cout << '\n';
        }

        if (n_clusters == 1)
            break;
    }
    cerr << "All particles connected to a single cluster.\n";
}