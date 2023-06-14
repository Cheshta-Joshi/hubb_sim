# Quantum Simulation of Lattice Models 
In this project, we study lattice models like **Tight Binding Model**, **Hubbard Model**. We find the basis states of these models, generate their Hamiltonian and study their eigenvalues and eigenstates. 

## 1. Tight Binding Model 

The Tight Binding Model describes the electronic structure of solids by considering the onsite energies and hopping energy of electrons. 

$$ H_{\text{sys}} = \sum_{i\sigma} \epsilon_i c_{i \sigma}^{\dagger} c_{i \sigma} + t \sum_{i \sigma} (c_{i\sigma}^\dagger c_{(i+1) \sigma} +h.c. ) $$

$$ = \sum_{k \sigma} (\epsilon + te^{-2 \pi i k/N}) c^\dagger_{k \sigma} c_{k \sigma} $$

Here, $c^\dagger_{l \sigma}$ is the creation operator which creates an electron with spin $\sigma$ at site l 

The tight binding Hamiltonian gives a Block diagonal matrix, with eigenvalues $\epsilon + 2t\cos k$, where k ranges from 1 to N. 

### 1.1 Plot of eigenvalues of a tight binding model 

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/eigenvalue%20plot%20TB.png" alt="Plot of eigenvalues of TB model" width="400">
</p>
  
### 1.2 Expectation of Occupation Number at Different Sites

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/eig_occ_n4e1TB.png" alt="Example Image" width="400"> <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ_n4e1_TB.png" alt="Example Image" width="400"> 
</p> 

## 2. One Orbital Single Spin Hubbard Model 
N site lattice model, with only one allowed spin. Describes the electronic structure of a solid by considering the onsite, hopping and nearest neighbour interaction term. 
$$H = \sum_i \epsilon _{i} c^\dagger _{i} c _{i} + t \sum _{i} (c^\dagger _{i} c _{i+1} + h.c.) + U \sum _{i} n _{i} n _{i+1}$$

### 2.1 Fock Dimension with increasing number of lattice sites

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/fock_dim_vs_N.png" alt="Plot of eigenvalues of TB model" width="400">
</p>

### 2.2 Subspace Dimension With Number of Electrons Occupied for Different N 

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/sub_dim_vs_n_elec.png" alt="Plot of eigenvalues of TB model" width="400">
</p>

### 2.3 Time Taken to Generate Basis States for Half-Population Subspace for Different N

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/basis_time_vs_lattice.png" width="400">
</p>


