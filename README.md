# HomogenizedEnergyTheory
Higher-order continuum theory based on Homogenized Energy Density(HED)
# High-Order Deformation Theory Based on Homogenized Strain Energy Density

A finite element implementation of a novel high-order deformation theory for solving inhomogeneous elastic deformation problems with size effects. This theory introduces **homogenized strain energy density** to capture microscopic deformation mechanisms using only **one physical scale parameter** — the size of the Representative Volume Element (RVE).

---

## 1. Limitations of Classical Continuum Mechanics

Classical continuum mechanics assumes that materials can be modeled as infinitely divisible continua composed of material points. Physical quantities (stress, strain, etc.) are defined at the centroid of the Representative Volume Element (RVE) and are assumed to be uniformly distributed within it.

### Key Limitations:

| Issue | Classical Theory Treatment | Real Material Behavior |
|-------|---------------------------|------------------------|
| **Homogenization Error** | Neglects higher-order spatial distribution of fields within RVE | Physical fields exhibit nonlinear spatial distribution in finite-sized RVEs |
| **Size Effects** | Cannot predict size-dependent mechanical response | Micro-scale specimens show enhanced stiffness (hardening) or reduced modulus (softening) |
| **Stress Singularity** | Predicts unphysical singularities at crack tips | Real materials exhibit finite stresses due to microstructure |
| **Nonlocal Interactions** | Absent between neighboring material points | Microscopic deformation gradients create nonlocal coupling |

&gt; **Critical Insight**: When the spatial distribution of physical fields within an RVE is nonlinear, the **averaged quantity** deviates from the **representative quantity** defined at the RVE centroid. This discrepancy is the **homogenization error** [Zhang & Wang, 2022].

---

## 2. Homogenized Strain Energy Density

### Core Idea

Instead of assuming strain energy density is defined at a single point, we compute the **averaged strain energy density** over a finite-sized RVE using Taylor series expansion:

$$\overline{U} = U + \frac{h^2}{24}\Delta U + \frac{h^4}{1920}\Delta^2 U + \cdots$$

where:
- $h$ = RVE size (the **only scale parameter**)
- $\Delta$ = Laplacian operator
- Higher-order terms capture the effect of inhomogeneous deformation

### Strain Energy Decomposition

The homogenized strain energy density naturally decomposes into:

$$\overline{U}(\boldsymbol{\varepsilon}) = U_c(\boldsymbol{\varepsilon}) + U_h(\boldsymbol{\varepsilon})$$

| Component | Expression | Physical Meaning |
|-----------|-----------|----------------|
| **Conventional** $U_c$ | $\frac{\lambda}{2}\varepsilon_{ii}^2 + \mu\varepsilon_{ij}\varepsilon_{ij}$ | Classical elastic strain energy |
| **Higher-Order** $U_h$ | Terms with $h^2$ and $h^4$ | Correction due to strain gradients |

### Key Features

- **Single Scale Parameter**: Only $h$ (RVE size) is needed, regardless of how many gradient terms are included
- **Hardening & Softening**: Cross-multiplication terms ($\varepsilon_{ij}\varepsilon_{ij,kk}$) can be positive or negative, enabling description of both behaviors
- **Physical Interpretation**: $h$ relates to material microstructure (grain size, polymer radius of gyration, etc.)

---

## 3. Solving Inhomogeneous Elastic Deformation Problems

### Augmented Lagrangian Method (ALM)

Direct solution of the governing equations involves **4th-order differential operators**, requiring $C^1$-continuous elements. To overcome this, we employ the **Augmented Lagrangian Method** to reformulate the problem:

**Primary Variables:**
- Displacement field $\mathbf{u} \in \mathbb{R}^3$
- Strain field $\boldsymbol{\varepsilon} \in \mathbb{R}^6$ (auxiliary variable)
- Lagrange multiplier $\boldsymbol{\lambda} \in \mathbb{R}^6$

**Constraint:** $\boldsymbol{\varepsilon} = \nabla_s \mathbf{u}$ (strain-displacement relation enforced via ALM)

### Total Potential Energy Functional

$$\mathcal{J}_{ALM}(\mathbf{u}, \boldsymbol{\varepsilon}, \boldsymbol{\lambda}) = \int_\Omega \overline{U}(\boldsymbol{\varepsilon}) \, dV - \mathcal{W}_{ext} + \int_\Omega \boldsymbol{\lambda}\cdot(\nabla_s\mathbf{u} - \boldsymbol{\varepsilon})\, dV + \frac{r}{2}\int_\Omega |\nabla_s\mathbf{u} - \boldsymbol{\varepsilon}|^2 \, dV$$

where $r$ is the penalty parameter.

### Finite Element Implementation

The weak form yields a **fully coupled system**:

$$
\begin{bmatrix}
\mathbf{K}_{uu} & \mathbf{K}_{u\varepsilon} & \mathbf{K}_{u\lambda} \\
\mathbf{K}_{\varepsilon u} & \mathbf{K}_{\varepsilon\varepsilon} & \mathbf{K}_{\varepsilon\lambda} \\
\mathbf{K}_{\lambda u} & \mathbf{K}_{\lambda\varepsilon} & \mathbf{0}
\end{bmatrix}
\begin{bmatrix}
\mathbf{u}_N \\ \boldsymbol{\varepsilon}_N \\ \boldsymbol{\lambda}_N
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{b} \\ \mathbf{0} \\ \mathbf{0}
\end{bmatrix}
$$

**Implementation Features:**
- Quadratic elements ($C^0$ continuity) sufficient — no need for complex $C^1$ elements
- 15 DOFs per node: 3 (displacement) + 6 (strain) + 6 (multiplier)
- Solved using UMFPACK direct solver
- Built on [MFEM](https://mfem.org/) finite element library

---

## 4. Solving Size Effect Problems

The theory successfully predicts size-dependent mechanical behavior in three canonical problems:

### 4.1 Micro-Cantilever Beam Bending (Hardening Effect)

| Thickness $H$ | Classical Prediction | Present Theory | Experiment |
|--------------|---------------------|----------------|------------|
| 16 μm | 112.8 GPa | ~115 GPa | ~115 GPa |
| 2 μm | 112.8 GPa | ~220 GPa | ~210 GPa |

- **Mechanism**: Higher-order strain energy density $U_h &gt; 0$ dominates as beam thickness decreases
- **Result**: Nominal Young's modulus increases with decreasing thickness

### 4.2 Perforated Plate Under Tension (Strain Concentration)

- Classical theory: Strain concentration factor $K = 3$ (constant, size-independent)
- Present theory: $K$ decreases as hole radius $a$ decreases
- **Mechanism**: Local hardening near hole edge redistributes strain field

### 4.3 Micro-Indentation (Softening Effect)

| Indentation Depth | Classical Theory | Present Theory | Observation |
|------------------|------------------|----------------|-------------|
| Shallow | $E^*$ constant | $E^* \approx E_{true}$ | Matches |
| Deep | $E^*$ constant | $E^* &lt; E_{true}$ | Softening |

- **Mechanism**: Negative cross-terms in $U_h$ dominate in highly deformed zone beneath indenter
- **Result**: Apparent elastic modulus decreases with increasing indentation depth

### Competition Between Conventional and Higher-Order Terms

$$\delta W_c + \delta W_h = \delta W_{ext}$$

- If $\delta W_h &gt; 0$: Conventional term decreases → **Hardening**
- If $\delta W_h &lt; 0$: Total strain energy reduced → **Softening**

---

## Upcoming: Extension to Quasi-Brittle Fracture

We are currently extending this framework to **quasi-brittle fracture problems** by:

- Incorporating **homogenized damage energy release rate**
- Coupling with **gradient-enhanced damage models**
- Addressing **stress singularity regularization** at crack tips using the natural length scale $h$

Stay tuned for updates on damage mechanics implementation and fracture simulations.

---

## References

1. Cao, Y., Zhang, C., & Wang, B. (2024). A New High-order Deformation Theory and Solution Procedure Based on Homogenized Strain Energy Density. *International Journal of Engineering Science*, 195, 103990.
2. Zhang, C., & Wang, B. (2022). Influence of nonlinear spatial distribution of stress and strain on solving problems of solid mechanics. *Applied Mathematics and Mechanics*, 43(9).
3. Fortin, M., & Glowinski, R. (1983). *Augmented Lagrangian Methods: Applications to the Numerical Solution of Boundary-Value Problems*. North-Holland.

---

## Citation

If you use this code in your research, please cite:

```bibtex
@article{cao2024high,
  title={A New High-order Deformation Theory and Solution Procedure Based on Homogenized Strain Energy Density},
  author={Cao, Yuheng and Zhang, Chunyu and Wang, Biao},
  journal={International Journal of Engineering Science},
  volume={195},
  pages={103990},
  year={2024},
  publisher={Elsevier}
}