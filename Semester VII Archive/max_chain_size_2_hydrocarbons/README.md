# The Algorithm (Current)
1. Start with a hexagonal lattice of size LATTICE_SIZE * LATTICE_SIZE for some LATTICE_SIZE.
1. Start the simulation for a given mole fraction of CO. The sum of the mole fraction of CO and the mole fraction of H<sub>2</sub> is assumed to be equal to 1.
1. Generate a random floating point number per iteration of the simulation. If this number is greater than the mole fraction of CO, choose H<sub>2</sub> as the reactant for this iteration, else choose CO.
1. If the chosen reactant is CO, generate two random integers in the range [0,LATTICE_SIZE], which denote the row and column of a site on the lattice. If this site is empty, a CO molecule adsorbs onto the site. Now, on this site:
    1. Initiate reaction with CO if possible.
    1. Initiate reaction with CH<sub>x</sub> if possible (0&leq;x&leq;3).
    1. Initiate desorption on CH<sub>x</sub> if possible (1&leq;x&leq;3).
1. If the chosen reactant is H<sub>2</sub>, generate two random integers in the range [0,LATTICE_SIZE], which denote the row and column of one of the two sites on the lattice, where this H<sub>2</sub> molecule will adsorb as individual H atoms. Generate another random integer in the range [0,2], which denotes the direction of the second site relative to the first site (right, bottom-left or bottom-right). If both these sites are empty, the H<sub>2</sub> molecule dissociates to two H atoms which individually adsorb to these sites.
1. Make a list of all the CO sites neighbouring one or both of the just adsorbed H-atoms and shuffle it. Now, on each of these sites:
    1. Initiate reaction with CO if possible.
1. Make a list of all the C sites neighbouring one or both of the just adsorbed H-atoms and shuffle it. Now, on each of these sites:
    1. Initiate reaction with CH<sub>x</sub> if possible (0&leq;x&leq;3).
    1. Initiate desorption on CH<sub>x</sub> if possible (1&leq;x&leq;3).
1. Stop the simulation when at least EQUILIBRIUM_VERIFICATION_STEPS iterations have passed and the following condition is true (This condition greatly decreases the program runtime with almost zero loss in accuracy of the results for say EQUILIBRIUM_VERIFICATION_STEPS = 50 * LATTICE_SIZE * LATTICE_SIZE and EQUILIBRIUM_VERIFICATION_THRESHOLD = 0):

    ```c++
    // Assuming that the current iteration is iteration i
    abs(sum_of_previous_EQUILIBRIUM_VERIFICATION_STEPS_cumulative_hydrocarbon_productions - (cumulative_hydrocarbon_production[i-EQUILIBRIUM_VERIFICATION_STEPS+1] + cumulative_hydrocarbon_production[i-    EQUILIBRIUM_VERIFICATION_STEPS+1]) * EQUILIBRIUM_VERIFICATION_STEPS / 2) <= EQUILIBRIUM_VERIFICATION_THRESHOLD)
    ```

1. Repeat all the above steps for different mole fractions and plot a graph between the desired properties and the mole fractions of CO.


# Some Important Subroutines

## Initiate reaction with CO if possible
1. For a given site, if this site is not occupied by a CO molecule, do nothing.
1. If this site is occupied by a CO molecule, make a list of its neighbouring sites, which are occupied by a hydrogen atom, irrespective of whether the hydrogen atom is already bonded to another C atom. (**Assumption:** Priority of formation of H<sub>2</sub>O over retaining the pre-existing C-H bonds)
1. If the size of the list is greater than or equal to 2, randomly choose two of the sites from this list.
1. The two hydrogen atoms occupying the chosen sites react with the CO molecule in consideration and form H<sub>2</sub>O, which desorbs leaving behind two empty sites, which were previously occupied by the H atoms in consideration (possibly also breaking their bond with another C one or more of them was bonded to previously) and a site occupied by a C atom, which was previously occupied by a CO molecule.

## Initiate reaction with CH<sub>x</sub> if possible (0&leq;x&leq;3)
1. For a given site, if this site is not occupied by a C atom, do nothing.
1. If this site is occupied by a C atom, make a list of its neighbouring sites, which are occupied by a hydrogen atom.
1. Shuffle the list.
1. Iterate over the sites in the list.
1. The hydrogen atom occupying the chosen site forms a bond with the C atom in consideration and we add this site to the list of sites of H bonded to the C atom in consideration, effectively increasing x by 1.
1. Continue iterating over the list of H-sites until x becomes equal to 4, i.e. a CH<sub>4</sub> molecule is formed or the list ends.
1. If a CH<sub>4</sub> molecule has been formed, it will desorb leaving behind 5 empty sites and increasing the hydrocarbon production by 1.

## Initiate desorption on CH<sub>x</sub> if possible (1&leq;x&leq;3)
1. For a given site, if this site is not occupied by a C atom, do nothing.
1. If this site is occupied by a C atom, make a list of its neighbouring sites, which are occupied by a C atom, which is bonded to the same number of H atoms as the C atom in consideration.
1. Randomly choose one of the sites from this list.
1. The chosen CH<sub>x</sub> molecule reacts with the CH<sub>x</sub> molecule in consideration forming C<sub>2</sub>H<sub>2x</sub>, which immediately desorbs leaving behind 2x+2 empty sites.
1. Increase the hydrocarbon production by 1.


# Assumptions

- Priority of formation of H<sub>2</sub>O over retaining the pre-existing C-H bonds.
