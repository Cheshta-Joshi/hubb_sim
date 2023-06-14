# Quantum Simulation of Lattice Models 
In this project, we study lattice models like **Tight Binding Model**, **Hubbard Model**. We find the basis states of these models, generate their Hamiltonian and study their eigenvalues and eigenstates. 

## Tight Binding Model 

The Tight Binding Model describes the electronic structure of solids by considering the onsite energies and hopping energy of electrons. 

$$ H_{\text{sys}} = \sum_{i\sigma} \epsilon_i c_{i \sigma}^{\dagger} c_{i \sigma} + t \sum_{i \sigma} (c_{i\sigma}^\dagger c_{(i+1) \sigma} +h.c. ) $$

$$ = \sum_{k \sigma} (\epsilon + te^{-2 \pi i k/N}) c^\dagger_{k \sigma} c_{k \sigma} $$

Here, $c^\dagger_{l \sigma}$ is the creation operator which creates an electron with spin $\sigma$ at site l 

The tight binding Hamiltonian gives a Block diagonal matrix, with eigenvalues $\epsilon + 2t\cos k$, where k ranges from 1 to N. 

### Plot of eigenvalues of a tight binding model 

![alt text](https://github.com/[Cheshta-Joshi]/[hubb_sim]/blob/[github.com/Cheshta-Joshi/hubb_sim/images]/eigenvalue plot TB.jpg?raw=true)
