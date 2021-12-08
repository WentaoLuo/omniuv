# Omni*UV*: An Omnipotent VLBI Simulation Toolkit
<img src="image/logo_flat.png" alt = "Logo" style= "display: block; margin-left:0; margin-right:auto; width:100%;"/>

---

## **Introduction**

We develop Omni*UV*, so as to fulfill the requirement of simulation for both ground and space VLBI observations. 

---

### **Functionality**
The main functionalities of the toolkit are:

- Trajectory calculation;
- Baseline *uv* calculation, by taking the availability of each station into account; 
- Visibility simulation for the given *uv* distribution, source structure and system noise;
- Image and beam reconstruction.

Omni*UV* supports two types of methods for visibility and imaging calculation (taking img -> vis for example):

- **FFT**: point sources are first assigned to the nearby grids in the image plane via certain kind of assignment function, then FFTed to the *uv* plane. Visibilities are reconstructed at given *uv* position via interpolation (assignment) function. 
    - Fast, therefore popular;
    - Suitable for small field;
    - *w* term cannot be taken into account;
    - Artifacts are introduced during gridding.
- **Direct summation**: direct calculation via Eq. 3 and 6 in the paper.
    - Slow, therefore requiring GPU speedup;
    - Suitable for wide field;
    - w* term is naturally supported;
    - No gridding artifact.

---

### **Station**
Omni*UV* supports the following types of ground and space stations:

- Earth fixed (ground);
- Earth orbit;
- Moon fixed;
- Moon orbit;
- Earth-Sun Lagrange 1/2 points; 
- Earth-Moon lagrange 1/2 points.

<img src="image/dataflow.png" alt = "dataflow" style= "display: block; margin-left:0; margin-right:auto; width:100%;"/>

The trajectories of stations are first calculated in ther own coordinate system, and then tranformed to the **Celestial Reference System (CRS)**. After that the station $uvw$ is calculated. The baseline $uvw$ is formed accordingly.

---

### **Source**
A Source in Omni*UV* actually refers to a **phase center**.

Omni*UV* supports $uvw$ calculations for multiple sources (phase centers). For each source, the input R.A. ($\alpha$) and decl. ($\delta$) are first converted to a unit vector in CRS: 

$\vec{s}=(x, y, z)$,

with:

$x=\cos\delta\cos\alpha$

$y=\cos\delta\sin\alpha$

$z=\sin\delta$

The defination of the $uvw$ system follows that in the standard textbook (Thompson et al. 2001):

Given North Pole direction $\vec{n}$ and source direction $\vec{s}$, $w$ is the same direction as $\vec{s}$:

$\vec{w}=\vec{s}$,

 $\vec{u}$ is the direction perpendicularto the plane defined by $\vec{n}$ and $\vec{w}$:

$\vec{u}$ = $\vec{n}\times \vec{w}$,

$\vec{v}$ is defined accordingly:

$\vec{v} = \vec{w}\times\vec{u}$.

The station $uvw$ is the projection of position $\vec{p}$ in the system defined above: 

$(u, v, w) = \vec{p}\cdot(\vec{u}, \vec{v}, \vec{w})$.

---

### **Image**
Visibility simulation requires an image as input. An image is binded to a source, and is described as a list of pixels (point sources). For each pixel, one has to provide the corresponding coordinate and flux. See `run_example.py` for a detailed example.

---

## **Prerequisites**

- Linux system
- Any version of Python 3.6, 3.7, 3.8, 3.9.

---

## **Installation**
 
- The toolkit is provided as precompiled package (`pyc` files). Please select the appropriate Python version and unpack. If you are not sure your Python version, please run the following command:

    `python --version`

- If you are in mainland China, please first run

    `sh pip_mirror.sh` 
  
    to select the pip host provided by aliyun.

- Install dependences:

    `pip install -r requirements.txt`

    Then the toolkit is ready.
