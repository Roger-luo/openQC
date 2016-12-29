# Multiple Algorithm Cooling

## Input state
$$
|\psi\rangle = \sum_{i=1}^{2^n}\sqrt{p_i}|e_i\rangle
$$

where, $|e_i\rangle$ is the eigenvector of the system. $p_i$ is the probability for each eigen state. $n$ is the number of qubits.

Define the cooling operation $\mathbf{C}$


$$
\begin{aligned}
\mathbf{C}|\psi\rangle &= \frac{1}{2}(1-i e^{i\gamma}U)|\psi\rangle|0\rangle+\frac{1}{2}(1+ie^{i\gamma}U)|\psi\rangle|1\rangle\\
 & = \sum_{k=1}^{2^n} \sqrt{p_k(\frac{1}{2}-a_k)} |e_k\rangle|0\rangle + \sqrt{p_k(\frac{1}{2}+a_k)}|e_k\rangle|1\rangle
\end{aligned}
$$

where $U = \exp{(-iH_st)}$, $H_s$ is the Hamiltonian of the system, $a_i = \frac{1}{2}\sin{(E_i t-\gamma)}$, $E_i$ is the eigen value for $|e_i\rangle$, then for $k$ steps cooling we have

$$
\begin{aligned}
C^{\otimes k}|\psi\rangle &= \frac{1}{2^k}\sum_{j = 0}^{k}C_k^j (1-ie^{i\gamma}U)^j(1+ie^{i\gamma}U)^{k-j}|\psi\rangle|\eta_j\rangle\\
&= \sum_{j=0}^{k}\sum_{i=1}^{2^n}\sqrt{p_i C_k^j (\frac{1}{2}-a_i)^j(\frac{1}{2}+a_i)^{k-j}}|e_i\rangle|\eta_j\rangle
\end{aligned}
$$

Then after $k$ steps cooling if the result of an projective measurement on ancilla qubits has $j$ zeros, the probability is

$$
\sum_{i=1}^{2^n} C_k^j p_i (\frac{1}{2}-a_i)^j(\frac{1}{2}+a_i)^{k-j}
$$

with an output state
$$
\frac{1}{2^k}C_k^j (1-ie^{i\gamma}U)^j(1+ie^{i\gamma}U)^{k-j}|\psi\rangle
$$