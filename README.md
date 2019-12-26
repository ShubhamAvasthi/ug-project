# UG Project
This is my UG Project on Kinetic Monte Carlo Simulation of the reaction of Syngas to produce hydrocarbons.

The presentations of the project I gave in my VI<sup>th</sup> and VII<sup>th</sup> semesters can be found [here](https://drive.google.com/open?id=187FT8hZuwe3KrRsSElPZhlCjB6SiKi99).

## Code Files:
- **simulation.cpp:** The main program that performs all simulations, written in C++.
- **visualizer.py:** A python script to visualize the lattice with the adsorbed species and hydrocarbon production in real time while the simulation is happening. 
- **visualize_results.ipynb:** A Jupyter notebook to visualize the production, distribution and structures of the hydrocarbon molecules formed across all simulations.

## Output Files:
- **hydrocarbon_productions.csv** : Contains the average hydrocarbon productions and coverage fractions of C, H, CO and vacant sites across all simulations.
- **hydrocarbon_sizes_distribution.txt** : Contains the size distribution of product hydrocarbon chains across all simulations. 
- **product_smiles.txt** : Contains the structures of the molecules formed across all simulations as Simplified Molecular Input Line Entry System (SMILES) strings.

## Recommended Compilation Command:
```posh
g++ -O3 -march=native simulation.cpp -o simulation
```

## The Algorithm (Current)
1. Start with a hexagonal lattice of size LATTICE_SIZE * LATTICE_SIZE for some constant LATTICE_SIZE.
1. Start the simulation for a given mole fraction of CO. The sum of the mole fraction of CO and the mole fraction of H<sub>2</sub> is assumed to be equal to 1 throughout the execution of the program.
1. Generate a random floating point number per iteration of the simulation. If this number is greater than the mole fraction of CO, choose H<sub>2</sub> as the reactant for this iteration, else choose CO.
1. If the chosen reactant is CO, generate two random integers in the range [0,LATTICE_SIZE-1], which denote the row and column of a site on the lattice. If this site is empty, a CO molecule adsorbs onto the site. Process this site (explained below).
1. If the chosen reactant is H<sub>2</sub>, generate two random integers in the range [0,LATTICE_SIZE-1], which denote the row and column of one of the two sites on the lattice, where this H<sub>2</sub> molecule will adsorb as individual H atoms. Generate another random integer in the range [0,1], which denotes the direction of the second site relative to the first site (right, bottom-left or bottom-right). If both these sites are empty, the H<sub>2</sub> molecule dissociates to two H atoms which individually adsorb to these sites. Make a list of all the CO sites neighbouring one or both of the just adsorbed H-atoms and shuffle it. Process each of these sites in a random order.
1. Stop the simulation when MAX_ITERATIONS_PER_SIMULATION trials have been completed or at least EQUILIBRIUM_VERIFICATION_STEPS iterations have passed and the following condition is true (This condition greatly decreases the program runtime with almost no loss in accuracy of the results for good enough constants EQUILIBRIUM_VERIFICATION_STEPS and EQUILIBRIUM_VERIFICATION_THRESHOLD):

    ```c++
    // Assuming that the current iteration is iteration i
    abs(sum_of_previous_EQUILIBRIUM_VERIFICATION_STEPS_cumulative_hydrocarbon_productions - (cumulative_hydrocarbon_production[i-EQUILIBRIUM_VERIFICATION_STEPS+1] + cumulative_hydrocarbon_production[i-    EQUILIBRIUM_VERIFICATION_STEPS+1]) * EQUILIBRIUM_VERIFICATION_STEPS / 2) <= EQUILIBRIUM_VERIFICATION_THRESHOLD)
    ```

1. Repeat all the above steps for different mole fractions and plot a graph between the desired properties and the mole fractions of CO.


## Lattice Operations Queues

Throughout the program, three queues are maintained. The purpose of the queues is listed below:
1. **First queue:** This queue consists of the sites, which are to be checked for the reaction of CO and H<sub>2</sub>.
1. **Second queue:** This queue consists of the sites, which are to be checked for the reaction of a CH<sub>x</sub> unit with H<sub>2</sub>.
1. **Third queue:** This queue consists of the sites, which are to be checked for the adsorption of a possible product.


## Some Important Subroutines

### Process site
1. Add the site to the first queue.
1. Process all queues.

### Process all queues
1. Initiate reaction with CO if possible for a site in the first queue. If not possible, proceed to the next step else push the site for which it was possible to the next queue and process all queues.
1. Initiate reaction with CH<sub>x</sub> if possible (0&leq;x&leq;3) for a site in the second queue. If not possible, proceed to the next step else push the site for which it was possible to the next queue and process all queues.
1. Initiate desorption on CH<sub>x</sub> if possible (1&leq;x&leq;3) for a site in the first queue. If not possible, return else process all queues.

### Initiate reaction with CO if possible
1. For a given site, if this site is not occupied by a CO molecule, return.
1. If this site is occupied by a CO molecule, make a list of its neighbouring sites, which are occupied by a hydrogen atom, irrespective of whether the hydrogen atom is already bonded to another C atom. (**Assumption:** Priority of formation of H<sub>2</sub>O over retaining the pre-existing C-H bonds).
1. If the size of the list is greater than or equal to 2, randomly choose two of the sites from this list.
1. The two hydrogen atoms occupying the chosen sites react with the CO molecule in consideration and form H<sub>2</sub>O, which desorbs leaving behind two empty sites, which were previously occupied by the H atoms in consideration (possibly also breaking their bond with another C one or more of them was bonded to previously) and a site occupied by a C atom, which was previously occupied by a CO molecule.
1. For each H that paticipated in the reaction, if it was previously bonded with a C site, add the C site to the first queue (although it should technically be added to the second queue, adding it to the first queue is still correct, we prefer adding to the first queue for a reason which is out of the scope of this document).

### Initiate reaction with CH<sub>x</sub> if possible (0&leq;x&leq;3)
1. For a given site, if this site is not occupied by a C atom, return.
1. If this site is occupied by a C atom, make a list of its neighbouring sites, which are occupied by a hydrogen atom.
1. Shuffle the list.
1. Iterate over the sites in the list.
1. The hydrogen atom occupying the chosen site forms a bond with the C atom in consideration and we add this site to the list of sites of H bonded to the C atom in consideration, effectively increasing x by 1.
1. Continue iterating over the list of H-sites until x becomes equal to 4, i.e. a CH<sub>4</sub> molecule is formed or the list ends.
1. If a CH<sub>4</sub> molecule has been formed, it will desorb leaving behind 5 empty sites and increasing the hydrocarbon production by 1.

### Initiate desorption on CH<sub>x</sub> if possible (1&leq;x&leq;3)
1. For a given site, if this site is not occupied by a C atom, return.
1. With this site as root, search a random branched chain.
1. If a branched chain was found, the branched chain desorps.

### Search branched chains
Finds a random product branched chain containing a given root site.
This is a quite complicated subroutine, so not including its details in the document.


## Assumptions

- Priority of formation of H<sub>2</sub>O over retaining the pre-existing C-H bonds.
