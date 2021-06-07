# Dynamics-Simulations
This is a collection of problem sets and my simulated results from my Olin 2021 Dynamics Class.

---

## Problem Set 1

### 1. Golf Ball Simulation (Air Resistance)

#### Setup
For this problem, a golf ball, with diameter D and mass m, is launched from point p with an initial velocity vector V<sub>0</sub>=v<sub>x<sub>0</sub></sub>i + v<sub>y<sub>0</sub></sub>j + v<sub>z<sub>0</sub></sub>k. The wind is blowing with air velocity V<sub>w</sub>=W<sub>x</sub>i+W<sub>y</sub>j+W<sub>z</sub>k resulting in a drag force, F<sub>D</sub>, acting on te ball. The drag force acts in the direction opposite its velocity and it proportional to the square of the velocity.

![Golf Ball Problem Setup](/Images/setup/golf.png)


#### Results
For a golf ball with the following characteristics and initial conditions
- Mass: 0.0459 kg
- Diameter: 4.27 cm
- Initial velocity: 130 mph at 45 degrees with respect to ground

and following drag properties
- Air density: 1.29 kg/m<sup>3</sup>
- Drag coefficient: 0.25

![golf simulation graph](/Images/results/golf-results.png)


### 2. Two Mass Hopper Simulation

#### Setup
For this problem, two masses, m<sub>upper</sub> and m<sub>lower</sub>, are connected by a linear spring with an unstreched length, L<sub>0</sub>, and a spring constant, K. The upper mass is given an initial downward displacement and is then released, resulting in the hopper jumping.

![Hopper Problem Setup](/Images/setup/hopper.png)

#### Results

---

## Problem Set 2

### 1. Marble in a Bowl Simulation

#### Setup
For this problem, a small marble is released at the rim of a hemispherical bowl of radius R = 200 mm with an initial velocity of 100 mm/s in the circumferential direction. The marble slides on the inside surface of the bowl without friction. Gravity acts in the downward vertical direction.

![Marble Problem Setup](/Images/setup/marble.png)

#### Results


### 2. Rollercoaster Simulation

#### Setup
For this problem, a roller coaster cart is released along a roller coaster track.  the track is defined by the function, y=h(1-B\*xL)cos<sup>2</sup>(3πx/(2L)) where 0<B<1 controls the height of the second peak of the track. One design constraint is that the cart must always travel over the second peak. 
The cart has air resistance. The magnitude of the force is F<sub>D</sub>=1/2ρC<sub>D</sub>V<sup>2</sup> where ρ is the air dentisty and C<sub>D</sub> is the drag coefficient, A is the frontal area of the car (which has mass, m), and V is the magnitude of the velocity of the car. The direction of the drag force is opposite to that of the car's velocity.

![Rollercoaster Problem Setup](/Images/setup/rollercoaster.png)

#### Results

---

## Problem Set 3

### 1. 3D Spring Pendulum

#### Setup
For this problem, a mass, m, is suspended on a massless linear spring to form a spherical pendulum. The unstreched lentgh of the spring is L<sub>0</sub>. The spring can only chane lenth in the radial direction.

![Spring Pendulum Problem Setup](/Images/setup/spring-pendulum.png)

#### Results

### 2. Planar Pendulum

#### Setup
For this problem, a plan pendulum which consists of mass, m, is suspended by a massless rod of lenth, l. Gravity acts in the downward vertical direction. The position of the mass can be described by the angle,θ, or by the cartesian coordinates (x,y).

![Planar Pendulum Problem Setup](/Images/setup/planar-pendulum.png)

#### Results

---

## Problem Set 4

### Double Pendulum

#### Seup
For this problem, a planar double pendulum consists of masses, m<sub>1</sub> and m<sub>2</sub>, suspended by massless rods of length l<sub>1</sub> and l<sub>2</sub> The position of the masses can be described by the angles, θ<sub>1</sub> and θ<sub>2</sub>, with respect to vertical. The pendulum is supported at point O with a frictionless pivot. Gravity acts in the downward vertical directions. 

![Double Pendulum Problem Setup](/Images/setup/double-pendulum.png)

#### Results

---

## Final Project

### Double Pendulum on a Cart
For this problem, a planar double pendulum consisting of masses, m<sub>2</sub> and m<sub>3</sub>, are suspended by massless rods of length l<sub>1</sub> and l<sub>2</sub> with a frictionless pivot upon a frictionless lateral cart with mass, m<sub>1</sub>. The position of the cart can be described by the distance x with respect to the origin, and the positions of the pendulums can be described by the angles, θ<sub>1</sub> and θ<sub>2</sub>, with respect to vertical. Gravity acts in the downward vertical directions. Furthermore, the cart and pendulum masses are subject to a drag force F<sub>D</sub>=C<sub>D</sub>V<sup>2</sup> where C<sub>D</sub> is the drag coefficient and the force is acting in the oposite direction of velocity. 

<img src="/Images/setup/cart-pendulum.png" width="450">

### Results

#### Drag Coeficient = 0

#### Drag Coeficient = 0.05

#### Drag Coeficient = 0.5
