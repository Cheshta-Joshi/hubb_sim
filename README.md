# Quantum Simulation of Lattice Models 
In this project, we study lattice models like **Tight Binding Model**, **Hubbard Model**. We generate matrices of these models, solve them using NumPy and discuss some of the properties of the eigenvalues and eigenvectors. We realize the need of quantum simulation of these models and propose a quantum method for solving them. We implement it on a quantum circuit in different ways and compare the solutions with the NumPy results. 

## 1. One-Orbital, Spinless Tight Binding Model 

The Tight Binding Model describes the electronic structure of solids by considering the onsite energies and hopping energy of electrons. 

$$ H_{\text{sys}} = \sum_{i} \epsilon_i c_{i }^{\dagger} c_{i } + t \sum_{i } (c_{i}^\dagger c_{(i+1) } +h.c. ) $$   

Here, $c^\dagger_{l}$ is the creation operator which creates an electron at site l. On exact, diagonalization using discrete fourier transform, we find out that the eigenvalues follow the cos function.   

$H_{\text{sys}} = \sum_{k } (\epsilon + 2 t \cos{\frac{2 \pi k}{N}}) c^\dagger_{k } c_{k }$  


<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/eigenvalue%20plot%20TB.png" alt="Plot of eigenvalues of TB model" width="300">
</p>  


### 1.1 Eigenstate Occupation    

We create the matrices for one electron occupation Tight Binding Model and solve it using NumPy. Using the eigenstates obtained from solving the mkatrix, we find the occupation of that eigenstate at all the lattice sites and plot the probability of the state occupying the sites. In the following figures we show how different onsite energies result in different occupation of the ground state (vec 0) and the eigenstate with maximum energy (vec 1). 

<center>
<table>
  <tr>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/e_peak.png" alt="Image 1" width="200" />
      <br />
      e peak
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/e_periodic.png" alt="Image 2" width="200" />
      <br />
      e periodic
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/e_inc.png" alt="Image 3" width="200" />
      <br />
      e increasing
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ_peak.png" alt="Image 4" width="200" />
      <br />
      occupation peak
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ_periodic.png" alt="Image 5" width="200" />
      <br />
      occupation periodic
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ_inc.png" alt="Image 6" width="200" />
      <br />
      occupation increasing
    </td>
  </tr>
</table>
</center>

### 1.2 Density of States   

It is used to understand the energy distribution in a system.  We find this using a delta or a lorentzian function. Once we have the eigenvalues of a system, we can find the distribution of these eigenvalues across an energy range. We can use the Lorentzian function along with the delta function to plot the density. The Lorentzian function broadens the delta function peaks corresponding to individual energy level giving a continuous density plot.   

The following figure gives the Density of states for N=100, number of centers = 150, gamma = 0.3, t=1 and e=[0]*N

$$ D(\omega )= \sum _k \delta (\omega - \epsilon _k)$$
<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos.png" alt="DOS" width="250">
</p>

## 2. One-Orbital, Spinless Hubbard Model   

N site lattice model, with only one allowed spin. Describes the electronic structure of a solid by considering the onsite, hopping and nearest neighbour interaction term.  

$$H = \sum_i \epsilon _{i} c^\dagger _{i} c _{i} + t \sum _{i} (c^\dagger _{i} c _{i+1} + h.c.) + U \sum _{i} n _{i} n _{i+1}$$  
We now write a generalised matrix for a Hubbard model. This can now handle any electron occupation. It takes the description of the lattice model adn the number of electorns as input and returns the Hamiltonian matrix of that subspace. For an N site lattice model and r occupied electron, the dimension of the Hamiltonian becomes $^NC_r$. The sum of all such dimensions is $2^N$ which is the dimension of the whole fock space. One can divide the fock-space problem into these subspace problems and solve each of them individually to extract more information about the fock-space.  

### 2.1 Eigenstate Occupation for different electron occupation  

In the following figures, we solve a 12 site Hubbard Model with $\epsilon _k = 1$, t=1 and U=1. We solve it for different electron occupation and find the eiegenstate occupation.  

<center>
<table>
  <tr>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ12_1.png" alt="Image 1" width="200" />
      <br />
      1 electron
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ12_3.png" alt="Image 2" width="200" />
      <br />
      1/4th occupation
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ12_6.png" alt="Image 3" width="200" />
      <br />
      1/2 occupation
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ12_11.png" alt="Image 4" width="200" />
      <br />
      1 Hole
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/occ12_9.png" alt="Image 5" width="200" />
      <br />
      3/4th occupation
    </td>
  </tr>
</table>
</center>


### 2.2 Density of States for different electron occupation  

In the following figures, we solve a 12 site Hubbard Model with $\epsilon _k = 1$, t=1 and U=1. We solve it for different electron occupation and find the density of states.  

<center>
<table>
  <tr>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos12_1.png" alt="Image 1" width="200" />
      <br />
      1 electron
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos12_3.png" alt="Image 2" width="200" />
      <br />
      1/4th occupation
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos12_6.png" alt="Image 3" width="200" />
      <br />
      1/2 occupation
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos12_11.png" alt="Image 4" width="200" />
      <br />
      1 Hole
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dos12_9.png" alt="Image 5" width="200" />
      <br />
      3/4th occupation
    </td>
  </tr>
</table>
</center>

We conclude from above data that a lattice problem with r-electron is the same as a problem with r holes. Their eigenvalues, eigenstates, their site occupation and the density of states have similar trends. 

## 3. Why Quantum?  
We can solve a lattice model effectively by breaking down the problem into its subspces and solving them using NumPy. However, as the number of sites increases, the dimension of subspaces increase and the time taken to generate their basis and then generate the Hamiltonian becomes exponentially larger.    

<center>
<table>
  <tr>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/fock_dim_vs_N.png" alt="Image 1" width="200" />
      <br />
      Fock dim vs N
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/sub_dim_vs_n_elec.png" alt="Image 2" width="200" />
      <br />
      Sub dim vs N
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/basis_time_vs_lattice.png" alt="Image 3" width="200" />
      <br />
      Basis time vs N
    </td>
  </tr>
</table>
</center>

In the above figures we can see how the fock space dimesnion increases exponentially, how the subspace dimesnion is maximum for half-filled electron and how that number increases appreciably with increasing lattice size. We can also see how the time taken to generate the half-filled basis increases exponentially with increasing lattice size. 

# 4. Quantum Simulation  

We will first discuss ways in which we can implement the lattice model hamiltonian on a quantum circuit. There are different transformations that we can do to change a fermionic operator to a spin operator. We disscuss one such transformation.

## 4.1 Jordan-Wigner Transformation  

So far, we were representing and solving the tight binding hamiltonian in terms of the fermionic creation $c^{\dagger} _j$ and annihilation $c_j$  operators. To be able to implement this on a quantum computer, we need to convert it into a spin problem. Jordan-Wigner transformation maps these operators to spin operators that we can eventually implement on quantum computers in the forms of quantum gates. We will use the following mapping strategy to do the same :   

$$c^{\dagger} _j = \frac{X_j-iY_j}{2}$$   

$$c_j = \frac{X_j+iY_j}{2}$$  

To preserve the anticommutation relations ($\lbrace \hat{c}^\dagger _\alpha, \hat{c}^\dagger _\beta \rbrace=0$$) of the fermionic operators, Z gates are added to add phase to the spin operators.  

$$c^{\dagger} _j = Z^{\otimes (j-1)} \otimes \Big( \frac{X-iY}{2} \Big) _j \otimes I^{\otimes (N-j)} $$  

$$c _j = Z^{\otimes (j-1)} \otimes \Big( \frac{X+iY}{2} \Big) _j \otimes I^{\otimes (N-j)} $$  

On using this transformation on our fermionic Hamiltonian, we get the following spin Hamiltonian :  

$$H_{sys} = \sum _{k=1}^{N} \textbf{I}^{\otimes (k-1)} \otimes \epsilon _k \big( \sigma ^- . \sigma ^+ \big)_k \otimes \textbf{I}^{\otimes(N-k)} + t \sum _{k=1}^{N} \textbf{I}^{\otimes (k-1)} \otimes  \big( \sigma ^- Z \big) _k \otimes \big( \sigma^+ \big) _{k+1} \otimes \textbf{I}^{\otimes (N-k-1)} + $$  

$$ t \sum _{k=1}^{N} \textbf{I}^{\otimes (k-1)} \otimes  \big( Z \sigma ^+ \big) _k \otimes \big( \sigma^- \big) _{k+1} \otimes \textbf{I}^{\otimes (N-k-1)} + U \sum _{k=1}^{N} \textbf{I}^{\otimes (k-1)} \otimes   \big( \sigma ^- . \sigma ^+ \big) _k \otimes \big( \sigma ^- . \sigma ^+ \big) _{k+1}\otimes \textbf{I}^{\otimes (N-k-1)} $$

where, $\sigma ^\pm = \frac{X \pm i Y}{2}$  

$$ H_{sys} = \sum _{k=1} ^{N} \epsilon _k \big( \sigma ^- . \sigma ^+ \big) _k + t \sum _{k=1} ^{N} \big( \sigma ^- _k \sigma ^+ _{k+1} - \sigma ^+ _k \sigma ^- _{k+1} \big) $$  

$$ H _{JW}= \sum _k \epsilon _k \Big( \frac{\textbf{I}-Z}{2}\Big)_k + t \sum _k \frac{X _k X _{k+1} + Y _k Y _{k+1}}{2} U \sum _{k=1}^N \Big( \frac{\textbf{I}-Z}{2}\Big) _k \Big( \frac{\textbf{I}-Z}{2}\Big) _{k+1} $$

As this expression is expressed in terms of the Pauli operators, we can implement this as a gate in a quantum circuit.  

## 4.2 Quantum Imaginary Time Evolution 

One can use imaginary time evolution to find the ground state and energy of a Hamiltonian. A state $| \psi \rangle$ of the system can be written as : $| \psi \rangle = \sum _i c_i |b _i\rangle$, where $|b_i\rangle$ are the basis states of the system and $c_i$ are complex amplitudes. The idea behind this algorithm is to apply the time evolution operator on this state and evolve it in imaginary time, $\tau=-it$, such that , $U | \psi \rangle =e^{-H\tau} |\psi \rangle$.   

$e^{-H \tau} | \psi \rangle = e^{-E_0 \tau} [c_0 |b_0 \rangle + e^{-(E_1-E_0) \tau} c_1 |b_1 \rangle +....+ e^{-(E_{n-1}-E_0) \tau} c_{n-1} |b_{n-1} \rangle] $  

One can see that for $\tau \xrightarrow{} \infty$, all the other coefficients tend to zero and we are left with the ground state $|b_0 \rangle$. We can also evolve the state for a number of discrete small time steps and end up with a similar result. In this following sections we will perform this imaginary time evolution using different approaches and evaluate ground state energy through them. 

### 4.2.1 Variational Quantum Imaginary Time Evolution  

We use a variational quantum algorithm (VQA) to implement imaginary time evolution. VQAs are an important class of quantum algorithms which are hardware agnostic and utilize both quantum and classical method to maximize efficiency and minimize cost. These algorithms use quantum circuits to evaluate a parameterized function which is optimized by updating the parameters using a classical optimizer and is fed back to the circuit for another iteration. The basic layout of a variational algorithm involves : Choosing a quantum system and constructing its Hamiltonian, Deciding a parameterized trial function, Choosing an appropriate optimizer, Choosing an efficient measurement technique and Reducing Cost and Mitigating errors. 

We make use of the  **VarQITE (Qiskit)** time evolution algorithm and input an evolution problem defined using the lattice Hailtonian. The varitional algorithm takes a parameterized ansatz as input and finds a set of parameters which minimizes the expectation value of the Hamiltonian. It is therefore important to take an ansatz which spans all the wavefunctions and can be trained efficiently. In the following sections we try out two different ansatz, compare it to the benchmarked SciPy and classical results.  

#### 4.2.1.a VarQITE with SU(2)  

When enough information about the interested eigenstate is not known, a heuristic ansatz that is arbitrarily parameterized and entangled, can be used, such as Efficient SU(2) ansatz. Here is a circuit of this ansatz repeated once. 

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/su2ansatz.png" alt="DOS" width="250">
</p>  

One can use the python package **qiskit-symb**, to see the parameterized state obtained from this circuit. One can also observe the amplitudes of different states in this statevector and modify their ansatz in order to achieve specific amplitudes for specific states.  

#### SciPyImaginarySolver  

We compare the efficiency of SU(2) VarQITE method with the SciPy method. This solves the application of the time evolution operator for imaginary time on the quantum state at all the discrete steps. The following figures show the comparison of performance of SciPyImaginaryEvolver and VarQITE with SU(2) ansatz. 
   
<center>
<table>
  <tr>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/su2_U0_N5.png" alt="Image 1" width="200" />
      <br />
      Estimated eigenvalues during time evolution
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/su2_vs_scipy_accuracy.png" alt="Image 2" width="200" />
      <br />
      Accuracy Comparison
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/su2_vs_scipy_time.png" alt="Image 3" width="200" />
      <br />
      Time Comparison
    </td>
  </tr>
</table>
</center>  


#### 4.2.1.b VarQITE with Dicke States' ansatz  

We use the circuit to create Dicke states, which are superposition of states with equal weights, and parameterize it to form our ansatz. 

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dicke_param.png" alt="DOS" width="250">
</p> 

We can see from the comparison below that the number of parameters required for it less than number of parameters required for SU(2) repeated thrice.  

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/param_comp.png" alt="DOS" width="250">
</p>   

We use this ansatz for solving different subspaces of our Hamiltonian and we can see that it gives good results not only for the ground states but for lowest eiegnavlues of all the subspaces. The following figures are for different N for U=0 and U=4 

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/all_occ_u0.png" alt="u0" width="850">
</p>   

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/all_occ_u4.png" alt="u4" width="850">
</p>   

We can also study the eigenvalues for different U  

<p align="center">
<img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/all_occ_diff_U.png" alt="eig_diff_U" width="850">
</p>   

We also test if we get the same density of states using the numpy calculated eigenvalues and using the iegenvalues obtained from VarQITE.  

<center>
<table>
  <tr>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dicke_dos_true.png" alt="Image 1" width="350" />
      <br />
      DOS NumPy
    </td>
    <td align="center">
      <img src="https://github.com/Cheshta-Joshi/hubb_sim/blob/main/images/dicke_dos_varqite.png" alt="Image 2" width="350" />
      <br />
      DOS VarQITE
    </td>
  </tr>
</table>
</center>

