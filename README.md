Model for the impedance of a disk electrode according to Newman's models, published in the papers with the following DOIs

* 10.1149/1.2423795
* 10.1149/1.2407464
* 10.1149/1.2424003

## Laplace equation and boundary conditions

![image](https://github.com/user-attachments/assets/4525d6bf-c374-4d53-9bc3-66099560d6d4)

The first boundary condition describes that the capacitive current equals the current due to electric field:

$$
i=C \frac{\partial\left(V-\Phi_0\right)}{\partial t}=-\left.\kappa \frac{\partial \Phi}{\partial z}\right|_{z=0, r \leq r_0}
$$

(later: can add faradaic current too). Potential on the electrode: $V=V_0 e^{j \omega t}$. Potential outside double laver: $\Phi_0$. Newman only considers the diffusion layer (net 0 free charge). Newman considers the double layer as an ideal capacitor.

$$
\begin{aligned}
& \text { BC2: }\left.\kappa \frac{\partial \Phi}{\partial z}\right|_{z=0, r>r_0}=0 \\
& \text { BC3: } \Phi\left(z^2+r^2 \rightarrow \infty\right)=0
\end{aligned}
$$


## Coordinate transformation

![image](https://github.com/user-attachments/assets/031604e6-b638-432b-84cc-5d13d5a6a19b)

New coordinates: ellipsoidal coordinates; $\xi$ : distance from the disk; $\eta$ : angular coordinate.

$$
\begin{aligned}
& z=r_0 \xi \eta \\
& r=r_0 \sqrt{\left(1+\xi^2\right)\left(1-\eta^2\right)}
\end{aligned}
$$

The Laplace equation (Poisson equation for zero net charge in the electrolyte) can then be written as

$$
\frac{\partial}{\partial \xi}\left[\left(1+\xi^2\right) \frac{\partial \Phi}{\partial \xi}\right]+\frac{\partial}{\partial \eta}\left[\left(1-\eta^2\right) \frac{\partial \Phi}{\partial \eta}\right]=0
$$

Quote from Newman: "to obtain a solution by the method of separation of variables we set

$$
\Phi=P(\eta) Q(\xi)
$$

The differential equations for $P$ and $Q$ are

$$
\begin{aligned}
\frac{d}{d \eta}\left[\left(1-\eta^2\right) \frac{d P}{d \eta}\right]+n P & =0 \\
& \frac{d}{d \xi}\left[\left(1+\xi^2\right) \frac{d Q}{d \xi}\right]-n Q=0
\end{aligned}
$$

(unquote). The solutions to these equations are Legendre polynomials. 

We define the potential in the solution as $\Phi=V_0 e^{j \omega t} U(r, z)$. If $\Phi$ satisfies Laplace's equation, then $U$ also satisfies it, because all $r$ and $z$-dependence is in $U$. Thus we obtain

$$
U=\sum_{n=0}^{\infty} B_n P_{2 n}(\eta) M_{2 n}(\xi)
$$

## Solving for the coefficients $B_n$

Using the expression for $\Phi$ in $\mathrm{BC}_1$ we get

$$
j \Omega\left(1-\sum_{n=0}^{\infty} B_n P_{2 n}(\eta)\right)=-\frac{1}{\eta} \sum_{n=0}^{\infty} B_n P_{2 n}(\eta) M_{2 n}^{\prime}(0)
$$

$M_{2 n}^{\prime}(0)$ was given in doi: 10.1149/1.2424003. Multiply the above equation by $\eta P_{2 n}(\eta)$ and integrate from 0 to 1 :

$$
B_m=-\frac{4 m+1}{M_{2 m}^{\prime}(0)} j \Omega\left(\int_0^1 \eta P_{2 m}(\eta) d \eta-\sum_{n=0}^{\infty} B_n \int_0^1 \eta P_{2 m}(\eta) P_{2 n}(\eta) d \eta\right)
$$

Thus we obtain an infinite set of equations: 

$$\left(\frac{1}{j \Omega \mathbf{D}}+\mathbf{M}\right) \vec{B}=\vec{a}$$

Here,

$$D_{mm} = -\frac{4 m+1}{M_{2 m}^{\prime}(0)}$$

(a diagonal matrix),

$$a_m = \int_0^1 \eta P_{2 m}(\eta) d \eta$$

and

$$M_{mn} = \int_0^1 \eta P_{2 m}(\eta) P_{2 n}(\eta) d \eta$$

Can be approximated by taking a finite number of terms, because $B_n$ and the integrals become smaller and smaller for increasing $n$. Then solve with numpy. The impedance only depends on $B_0$. See Newman's papers -- they say that the total current is given by

$$
\begin{aligned}
I=\int_0^{r_0} i 2 \pi r d r=-2 \pi r_0{ }^2 \kappa & \left.\int_0^1 \eta \frac{\partial \Phi}{\partial z}\right|_{z=0} d \eta \\
& =-2 \pi r_0 \kappa V_o e^{j \omega t} B_0 M_{o^{\prime}}(0)
\end{aligned}
$$

which, with $M_0{ }^{\prime}(0)=-2 / \pi$, leads to 

$$
Z=V / I=1 / 4 r_0 \times B_0.
$$
