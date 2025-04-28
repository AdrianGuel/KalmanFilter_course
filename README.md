# Kalman Filter Course

This repository contains MATLAB scripts for learning and applying **Kalman Filtering** techniques, including classical Kalman Filters, Extended Kalman Filters (EKF), and Unscented Kalman Filters (UKF), with applications to mechanical systems and stochastic processes.

---

## 📁 Repository Structure and Contents

| Folder | Script | Description |
|:------|:------|:-------------|
| **`MHE/`** | `protomhe.m` | Prototype for Moving Horizon Estimation (MHE). |
| **`Section1/`** | `Gaussian_plot.m` | Visualizes Gaussian distributions. |
|  | `KramersSS.m` | Computes steady-state of Kramers equation. |
|  | `Kramers_euler_maruyama.m` | Simulation of Kramers dynamics with Euler-Maruyama method. |
|  | `Mult_Gaussian_animation.m` | Animation of multivariate Gaussian evolution. |
|  | `brownian_motion.m` | Simulation of basic Brownian motion. |
|  | `chol_LDL.m` | Cholesky decomposition via LDL factorization. |
| **`Session_2/`** | `MSD_Euler.m` | Simulates a mass-spring-damper system using Euler method. |
|  | `MSD_RK4vsRK1.m` | Compares Runge-Kutta 4th order with Euler integration. |
|  | `discrete4order_example.m` | Example of discrete-time 4th order system simulation. |
| **`session_3/`** | `CL_Cart_KF.m` | Closed-loop cart simulation with Kalman filter. |
|  | `Cart_Kalman_Filter.m` | Classical Kalman Filter applied to a cart system. |
|  | `Cart_generative_model.m` | Generative model for cart dynamics simulation. |
| **`Session_4/`** | `CA_KAlmanFilter_Q.m` | Kalman filter with process noise tuning. |
|  | `CA_model.m` | Constant-acceleration model setup. |
|  | `nonliearMSD_Kalman.m` | Kalman filtering on a nonlinear MSD system. |
|  | `nonlinearMSD.m` | Nonlinear MSD system dynamics. |
| **`session_5/`** | `EKF_example.m` | Example of Extended Kalman Filter (EKF) application. |
|  | `MSD_dualSPE.m` | Dual-state estimation using EKF. |
| **`session_6/`** | `Example_UKF.m` | Unscented Kalman Filter (UKF) example simulation. |
|  | `adaptivecontrol_EKF_cart.m` | Adaptive control using EKF on a cart model. |
| **`session_7/`** | `MSD_SRUKF.m` | Square Root UKF applied to a mass-spring-damper system. |
|  | `MSD_UKF.m` | Regular UKF applied to mass-spring-damper system. |
| **`session_8/`** | `nonlinearMSD_SRUKF.m` | Square Root UKF applied to a nonlinear MSD model. |

---

## 📚 Topics Covered

- Basics of **Gaussian processes** and **Brownian motion**
- **Euler-Maruyama simulations** for stochastic differential equations
- **Kalman Filter (KF)** for linear systems
- **Extended Kalman Filter (EKF)** for nonlinear systems
- **Unscented Kalman Filter (UKF)** and **Square Root UKF (SRUKF)**
- Applications to mechanical systems like mass-spring-damper models and cart dynamics
- Introduction to **Moving Horizon Estimation (MHE)**

---

## 🛠 Requirements

- **MATLAB** R2020b or newer
- Recommended Toolboxes:
  - Control System Toolbox
  - Statistics and Machine Learning Toolbox (for advanced filtering tasks)

---

## 🚀 How to Use

1. **Clone the repository**:
   ```bash
   git clone https://github.com/YOUR_USERNAME/KalmanFilter_course.git
   ```
2. **Add to MATLAB path**:
   ```matlab
   addpath(genpath('path_to_cloned_folder'))
   ```
3. **Run examples**, for instance:
   ```matlab
   cd Section1
   run('Gaussian_plot.m')
   ```

---

## ✨ Highlights

- Full pipeline from basic stochastic simulation to advanced Kalman filtering.
- Practical examples on mechanical system state estimation.
- Includes Moving Horizon Estimation (prototype).

---

## 📜 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details (if available).

---

## 👨‍💻 Author

Developed by Adrian Guel.
