# Parameter-estimation-of-damage-using-PINNs
PINNs leverages the governing physical laws that dictate the behavior of time-dependent systems to act as a regularization agent. Parametric PINNs are used to identify
and characterize damage parameters within the structural domain. Ultrasonic imaging with a two-dimensional wave equation for isotropic material is used for simplicity. The stationary ultrasound excitation source is modelled with a five-cycle Hanning window sinusoidal burst which is applied at the center of the plate. The damage parameters are modelled using Gaussian function with critical parameters such as location, size and extent of damage which are to be detected. For simplicity,
only the damage location is parameterized, while damage size and magnitude is kept fixed. The underlying physics governing the behavior of the damage model is
encoded into the loss function in the form of damage equation, complemented by the inclusion of boundary and initial loss terms. The parametric PINN approach is summarized as follows:

![PPINN_working_schematic](https://github.com/user-attachments/assets/9c990501-633a-4125-ba93-a8fcb324acfc)

The repository currently contains forward PINN and parametric PINN for one location parameter (1 pair = 2 damage locations along x and y directions). Traditional PINN implementation faces difficulty to converge, and due to the very low order of displacement field of magnitude, 1e-26 causes exploding loss value. This makes the model optimize for a long duration of about 2-3 days. Hence, two alternative scaling approach is used to achieve convergence. One is the scaling the governing PDE and other is scaling the inputs and outputs of the neural network keeping the PDE same.

The scaled PDE is obtained as below :

![scaled_PDE](https://github.com/user-attachments/assets/370e234e-aa0a-40a9-ae1f-35667ce5a6bc)

The scaling of inputs and outputs of the network is summarized as follows:
![schematic_for_scaled_PDE](https://github.com/user-attachments/assets/7f2bb76b-61f9-448a-ba87-a51a1a700258)

Nevertheless, the PINN needed additional sampling points near source region to initiate the source accurately. The results obtained are attached in this repository.

