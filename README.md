# Quantum Simulation of Lattice Models 
In this project, we study lattice models like **Tight Binding Model**, **Hubbard Model**. We find the basis states of these models, generate their Hamiltonian and study their eigenvalues and eigenstates. 

## 1. One-Orbital, Spinless Tight Binding Model 

The Tight Binding Model describes the electronic structure of solids by considering the onsite energies and hopping energy of electrons. 

$$ H_{\text{sys}} = \sum_{i\sigma} \epsilon_i c_{i \sigma}^{\dagger} c_{i \sigma} + t \sum_{i \sigma} (c_{i\sigma}^\dagger c_{(i+1) \sigma} +h.c. ) $$

$$ = \sum_{k \sigma} (\epsilon + 2 t \cos{\frac{2 \pi k}{N}}) c^\dagger_{k \sigma} c_{k \sigma} $$

Here, $c^\dagger_{l \sigma}$ is the creation operator which creates an electron with spin $\sigma$ at site l 

The tight binding Hamiltonian gives a Block diagonal matrix, with eigenvalues $\epsilon + 2t\cos k$, where k ranges from 1 to N. 

### 1.1 Eigenvalues of a tight binding model 

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/eigenvalue%20plot%20TB.png" alt="Plot of eigenvalues of TB model" width="400">
</p>
  
### 1.2 Occupation of Eigenstates

#### 1.2.1 Constant Onsite energies 

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/e_const.png" alt="Example Image" width="300">   &nbsp; &nbsp; &nbsp; &nbsp;      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ_const.png" alt="Example Image" width="300"> 
</p> 

#### 1.2.2 Constant Onsite energies but peaks at the middle point

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/e_peak.png" alt="Example Image" width="300">   &nbsp; &nbsp; &nbsp; &nbsp;      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ_peak.png" alt="Example Image" width="300"> 
</p> 

#### 1.2.3 Periodic Onsite energies 

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/e_periodic.png" alt="Example Image" width="300">   &nbsp; &nbsp; &nbsp; &nbsp;      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ_periodic.png" alt="Example Image" width="300"> 
</p> 

#### 1.2.4 Increasing Onsite energies 

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/e_inc.png" alt="Example Image" width="300">   &nbsp; &nbsp; &nbsp; &nbsp;      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ_inc.png" alt="Example Image" width="300"> 
</p> 

### 1.3 Density of States 

$$ D(\omega )= \sum _k \delta (\omega - \epsilon _k)$$
<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos.png" alt="Plot of eigenvalues of TB model" width="400">
</p>

### 1.4 Jordan-Wigner Transformation  

So far, we were representing and solving the tight binding hamiltonian in terms of the fermionic creation $c^{\dagger} _j$ and annihilation $c_j$  operators. To be able to implement this on a quantum computer, we need to convert it into a spin problem. Jordan-Wigner transformation maps these operators to spin operators that we can eventually implement on quantum computers in the forms of quantum gates. We will use the following mapping strategy to do the same :   

$$c^{\dagger} _j = \frac{X_j-iY_j}{2}$$   

$$c_j = \frac{X_j+iY_j}{2}$$  

To preserve the anticommutation relations ($\lbrace \hat{c}^\dagger _\alpha, \hat{c}^\dagger _\beta \rbrace=0$$) of the fermionic operators, Z gates are added to add phase to the spin operators.  

$$c^{\dagger} _j = Z^{\otimes (j-1)} \otimes \Big( \frac{X-iY}{2} \Big) _j \otimes I^{\otimes (N-j)} $$  

$$c _j = Z^{\otimes (j-1)} \otimes \Big( \frac{X+iY}{2} \Big) _j \otimes I^{\otimes (N-j)} $$  

On using this transformation on our fermionic Hamiltonian, we get the following spin Hamiltonian :  

$$H_{sys} = \sum _{k=1}^{N} \textbf{I}^{\otimes (k-1)} \otimes \epsilon _k \big( \sigma ^- . \sigma ^+ \big)_k \otimes \textbf{I}^{\otimes(N-k)} + \sum _{k=1}^{N} \textbf{I}^{\otimes (k-1)} \otimes  \big( \sigma ^- \big) _k \otimes \big( \sigma^+ \big) _{k+1} \otimes \textbf{I}^{\otimes (N-k-1)} - \sum _{k=1}^{N} \textbf{I}^{\otimes (k-1)} \otimes  \big( \sigma ^+ \big) _k \otimes \big( \sigma^- \big) _{k+1} \otimes \textbf{I}^{\otimes (N-k-1)}$$  


where, $\sigma ^\pm = \frac{X \pm i Y}{2}$  


As this expression is expressed in terms of the Pauli operators, we can implement this as a gate in a quantum circuit.  

## 2. One Orbital Spinless Hubbard Model 
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

### 2.4 Occupation and Density of States  

We simulate the Hamiltonian of the 1-orital spinless Hubbard model, calculate its eigenvalues, eigenstates and find some properties related to them. Expectation of occupation number of a site with respect to an eigenstate tells the probability of the eigenstate occupying that site. While, the Density of States tells the distribution of these eigenstates around certain points. We take a Hubbard model with 12 lattice sites and study it for different electron occupation. Here, we study these two properties for the ground state and for an eigenstate with medium eigenvalue, eigvec[N//2]  
Note that the onsite energy is taken to be constant , e= [1,1,1,...,1] , hopping constant, t is 1 and interaction coefficient, U is also 1.  

#### 2.4.1 N=12, r=1  

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ12_1.png" alt="Example Image" width="300">   &nbsp; &nbsp; &nbsp; &nbsp;      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos12_1.png" width="300"> 
</p>   

#### 2.4.2 N=12, r=11  

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ12_11.png" alt="Example Image" width="300">   &nbsp; &nbsp; &nbsp; &nbsp;      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos12_11.png" width="300"> 
</p>   

#### 2.4.3 N=12, r=3  

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ12_3.png" alt="Example Image" width="300">   &nbsp; &nbsp; &nbsp; &nbsp;      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos12_3.png" width="300"> 
</p>   

#### 2.4.4 N=12, r=9  

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ12_9.png" alt="Example Image" width="300">   &nbsp; &nbsp; &nbsp; &nbsp;      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos12_9.png" width="300"> 
</p>  

#### 2.4.5 N=12, r=6  

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ12_6.png" alt="Example Image" width="300">   &nbsp; &nbsp; &nbsp; &nbsp;      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos12_6.png" width="300"> 
</p> 
