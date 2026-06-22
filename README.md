# Vibrations and Shape-evolution Control of Hyperelastic Beams Driven by Dynamic Growth

This repository provides the supplementary documents, source code, and validation profiles for the paper:
**Vibrations and Shape-evolution Control of Hyperelastic Beams Driven by Dynamic Growth**.

---

## 📌 Overview

* This paper establishes a finite strain beam theory incorporating **dynamic growth** effects.
* The forward analysis explicitly reveals the nonlinear vibration characteristics of the beam under the influence of base motions and spatiotemporal growth fields. Comparisons with 3D nonlinear finite element simulations validate the proposed model for capturing large-deformation dynamic responses. Nevertheless, the present model is applicable within the framework of plane-strain beams and diagonal through-thickness linear growth, where the accuracy depends on the truncated orders in $\varepsilon$ and $h_0$. Furthermore, an analytical inverse framework is established, enabling high-precision configuration maintenance and continuous time-varying path planning.

---

## 📂 Supplementary Documents

All files are hosted in this repository. You can find:

* **`Movie 1`**: Visualizes the dynamic responses and shape-programming processes. 
    * *Tip*: To download or view the video, click the file in the repository, then click **"View raw"**.
* **`Source Code`**: The complete, runnable MATLAB implementation of the high-order finite difference algorithm and adaptive time integration scheme used to numerically solve the nonlinear dynamic beam equations.

---

## 🛠️ Theory & Numerical Implementation

### 1. Forward Vibration Analysis
Within the framework of nonlinear elasticity, we derive the asymptotic equations with $O(h^2)$ accuracy that seamlessly couples bending stiffness, inertial corrections, and dynamic growth driving forces. Using the provided MATLAB script, the forward dynamic responses under base motions (accelerated/periodic) and prescribed spatio-temporal growth fields can be predicted.
![Theory](https://github.com/JL-WEI-zszj/Dynamic-growth-beam/blob/main/theory.png)

### 2. Inverse Problem & Configuration Control
To achieve precise shape control, we build an analytical inverse framework. For any prescribed target shape evolution path, the explicit required growth functions are derived analytically. 
![Results](https://github.com/JL-WEI-zszj/Dynamic-growth-beam/blob/main/Results.png)

Our results demonstrate that this inverse strategy can:
1.  **High-precision configuration maintenance**.
2.  **Continuous time-varying path planning**.

---

## 💻 3D Finite Element Verification

The theoretical and numerical solutions from our beam model are validated against **3D nonlinear finite element simulations** performed in COMSOL Multiphysics. 

* **Model Setup**: Incompressible Neo-Hookean constitutive model , discretized using second-order tetrahedral elements with plane strain constraints.
* **Results**: The asymptotic model shows excellent quantitative agreement with the 3D FEM results , significantly improving computational efficiency.

## 📝 Citation

If you find this repository or the paper helpful for your research, please cite our work:

```bibtex
