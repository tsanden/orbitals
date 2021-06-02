# Orbitals

## What is it?

Orbitals is a python script that calculates the probability densities for atomic orbitals of hydrogen given the quantum numbers. It also allows for dot density and contour plotting of the atomic orbitals of hydrogen.

## Installation

Import the _orbitals.py_ script into your notebook or script:

    import orbitals as orb

The script contains two functions, _prob_orbital()_ and _plot_orbital()_, which returns probabilities and plotting respectively. Note that the functions accepts an array of quantum numbers (_n_, _l_, _m_), where the same position in each quantum number array corresponds to the same orbital. Check out the _orbitals.ipynb_ notebook for examples.