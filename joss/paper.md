---
title: 'EchemFEM: A Firedrake-based Python package for electrochemical transport'
tags:
  - Python
  - Firedrake
  - Finite Element Method
  - electrochemistry
authors:
  - name: Thomas Roy
    orcid: 0000-0000-0000-0000
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Author Without ORCID
    affiliation: 2
affiliations:
 - name: Lawrence Livermore National Laboratory, CA, USA
   index: 1
 - name: Institution Name, Country
   index: 2
date: 12 December 2023
bibliography: paper.bib

---

# Summary
<!---  high-level functionality and purpose of the software for a diverse, non-specialist audience --->
<!--- motivation --->
The shift from fossil fuels towards renewable energy brings about a substantial increase in clean but intermittent electricity.
Thankfully, diverse electrochemical technologies, including energy storage and electrochemical manufacturing, can harness this surplus energy that would otherwise go to waste.
Managing the growing prevalence of renewable energy underscores the importance of developing and scaling up these technologies.
Likewise, the electrification of transport creates an increasing need for energy-dense electrochemical energy storage devices.
Naturally, simulation tools are required to assist in the design of performant and industrial-scale electrochemical devices.

<!--- How modeling is used --->
Modeling and simulation are used extensively to describe the physics of the electrochemical and transport mechanisms in electrochemical devices.
These devices have various applications, such as batteries and supercapacitors, which are used for energy storage.
Flow batteries can have a similar function, but with a flowing electrolyte.
Electrolyzers can be used to transform carbon dioxide into useful products or create hydrogen and oxygen from water.
On the other hand, proton-exchange membrane fuel cells reverse this process to harness electricity from hydrogen.
While these devices vary wildly in purpose, the governing equations used to describe them are very similar.
The transport of charged chemical species in a fluid is often modeled using the Nernst-Planck equation,
which includes the usual advection and diffusion transport as well as *electromigration*, where charged species are transported by an electric field.

<!--- EchemFEM --->
EchemFEM provides a high-level user interface for a finite element implementation of the Nernst-Planck equation.
The user is simply required to provide physical parameters as well as functions describing the chemical reactions.
Then, the desired transport physics are selected using keyword arguments.
Ionic charge can either be modeled using the Poisson equation or the electroneutrality approximation.
The simulated devices can have resolved electrolyte-electrode interfaces or homogenized porous electrodes, in which case electron conduction is also modeled.
Additionally, finite size effects are available to allow the use of models such as Generalized Modified Poisson-Nernst-Planck (GMPNP) [@wang2013simulations].

<!--- Firedrake --->
EchemFEM is based on Firedrake [@FiredrakeUserManual], an open-source finite element package,
enabling straightforward implementation of the governing equations in Python.
Firedake has access to scalable, customizable, solvers through its interface with PETSc [@petsc-user-ref; @petsc-web-page], allowing for parallelization and scalability on computing clusters.
This balance between usability and scalability permits a seamless transition from prototyping to large-scale simulation.
EchemFEM leverages Firedrake's capabilities while further increasing the ease-of-use.
Indeed, since the governing equations are already implemented, little to no knowledge of Firedrake and the finite element method is required to use EchemFEM.

Examples: CO2 reduction. flow reactors

# Statement of need
<!--- section that clearly illustrates the research purpose of the software and places it in the context of related work --->
<!--- Research Need --->
How simulation can help electrochemistry.
What questions are we trying to answer here?
Scaling up technologies: scaling up simulations.
Higher dimensions: many codes are just 1D.
Fluid flow and Architected electrodes: higher-dimensional effects do matter.

<!--- Other codes --->
At the moment, the most common codes used for electrochemistry simulations are commercial software.
COMSOL Multiphysics<sup>&reg;</sup>, with its detailed electrochemistry module, is the most popular, while Simcenter<sup>&trade;</sup> STAR-CCM+<sup>&trade;</sup> is also used commonly for flowing systems.
These programs provide simple graphical user interfaces (GUI), which allow users to quickly set up new simulations.
Additionally, other physics modules such as fluid dynamics are available and can usually be coupled with the electrochemistry simulation.
However, there are several drawbacks to using such commercial software.
For instance, license fees can be prohibitively expensive, therefore limiting collaboration.
Furthermore, the closed nature of the source code limits the flexibility of the software.
Indeed, it is not possible to implement new discretization schemes and preconditioning approaches that may be required for numerical stability or scalability, respectively.
Finally, since everything needs to be set up through the GUI, scripting and coupling to other software are difficult tasks.

There is a growing number of open-source software for electrochemistry, especially Python-based packages [@zheng2023python], many of which are specialized for specific applications, notably batteries.
One such package, PyBaMM [@sulzer2021python], is a battery modelling code with a flexible implementation, allowing for new models and numerical methods to be tested.

<!--- Why echemfem --->
EchemFEM provides a general framework for simulating electrochemical transport: it is not specific to an application.
Since it is based on Firedrake, any additional physics that can be implemented in a finite element framework can be coupled to EchemFEM.
In one of the demos, the Navier-Stokes equations are solved in a reactor with an irregular surface, providing a velocity field for the transport equations.

In some cases, for example for fast flows, stabilization schemes that are not offered in other software may be required.
For continuous Galerkin (CG) elements, a streamline-upwind Petrov-Galerkin (SUPG) method for the Nernst-Planck equation is provided.
For discontinuous Galerking (DG), a custom upwind scheme for the Nernst-Planck equation is used [@roy2023scalable].
In both cases, the upwinding considers the combined advection-migration ``velocity''.

As opposed to commercial software, custom scalable solvers are available in Firedrake.
A plethora of solver options are available through simple PETSc keywords and custom operators for preconditioning can be defined using Firedrake [@Mitusch2019].
In @roy2023scalable, scalable block preconditioners were developed for the electroneutral Nernst-Planck equations with DG and implemented in EchemFEM.

Combining EchemFEM with other Python packages is rather simple.
In @govindarajan2023coupling, multi-scale simulations for CO<sub>2;</sub> reduction in flow reactors are performed by coupling a microkinetics model from CatMAP [@catmap] with the GMPNP transport model from EchemFEM.
The simulations are two orders of magnitude faster than a previous implementation where the transport is done in COMSOL Multiphysics<sup>&reg;</sup>, due to the tedious interface between the commercial software and CatMAP.

Another desirable feature of Firedrake is its automatic adjoint capabilities, facilitating the straightforward solution of PDE-constrained optimization problems [@Mitusch2019], which have already been used for electrochemical applications [@roy2022topology].
Optimization problems using EchemFEM are currently being investigated.

# Acknowledgements

This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory (LLNL) under Contract DE-AC52-07NA27344, and was partially supported by a Cooperative Research and Development Agreement (CRADA) between LLNL and TotalEnergies American Services, Inc. (affiliate of TotalEnergies SE) under agreement number TC02307 and Laboratory Directed Research and Development (LDRD) funding under projects 19-ERD-035 and 22-SI-006.

# References
