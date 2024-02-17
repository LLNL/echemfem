---
title: 'EchemFEM: A Firedrake-based Python package for electrochemical transport'
tags:
  - Python
  - Firedrake
  - Finite Element Method
  - electrochemistry
authors:
  - name: Thomas Roy
    orcid: 0000-0002-4286-4507
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Julian Andrej
    orcid: 0000-0001-7661-4840
    affiliation: 1
  - name: Aymeric Antimes
    affiliation: 2
  - name: Victor A. Beck
    orcid: 0000-0002-0625-9545
    affiliation: 1
  - name: Victoria Ehlinger
    orcid: 0000-0001-7333-1271
    affiliation: 1
  - name: Florian Euzenat
    affiliation: 2
  - name: Nitish Govindarajan
    orcid: 0000-0003-3227-5183
    affiliation: 1
  - name: Jack Guo
    orcid: 0000-0003-4090-9289
    affiliation: 1
  - name: Tiras Y. Lin
    orcid: 0000-0002-3377-9933
    affiliation: 1
  - name: Thomas Moore
    orcid: 0000-0003-0802-5547
    affiliation: 3
affiliations:
 - name: Lawrence Livermore National Laboratory, CA, USA
   index: 1
 - name: Institution Name, Country
   index: 2
 - name: Queensland University of Technology, Australia
   index: 3
date: 8 February 2024
bibliography: paper.bib

---

# Summary
<!---  high-level functionality and purpose of the software for a diverse, non-specialist audience --->
<!--- motivation --->
The transition from fossil fuels to renewable energy has brought about a rapid increase in the availability of clean electricity.
However, electricity generated from sources such as wind and solar are limited to intermittent operation due to daily and seasonal variation.
One solution is to utilize electrochemical devices in energy storage and electrochemical manufacturing applications, where they can harness surplus energy and decarbonize chemical industries traditionally reliant on petrochemical feedstocks.
Managing the growing prevalence of renewable energy underscores the importance of developing and scaling up these technologies.
Likewise, the electrification of transport creates an increasing need for energy-dense electrochemical energy storage devices such as batteries and supercapacitors.
Naturally, simulation tools are required to assist in the design of efficient and industrial-scale electrochemical devices.


<!--- How modeling is used --->
Modeling and simulation are used extensively to describe the physics of the electrochemical and transport mechanisms in electrochemical devices.
These devices have many applications, from miniaturized lithium-ion batteries for medical devices up to industrial-scale hydrogen fuel cells for backup power generation.
Energy storage devices include batteries and supercapacitors, as well as flow batteries which utilize a flowing electrolyte instead of a stationary liquid or polymer electrolyte.
Electrolyzers are devices, which use electrical energy to perform electrochemical reactions.
Some current industrial applications for electrolysis include the color-alkali process for the production of chlorine gas and the Hall-Héroult process for aluminum production.
Active areas of research include the development of electrolyzers that transform carbon dioxide into useful chemicals and electrolyzers that create hydrogen from water.
In the reverse process, fuel cells use fuels such as hydrogen to generate electricity.
While electrochemical devices span many scales and industries, the governing equations and underlying physical phenomena remain similar.

The transport of charged chemical species in a fluid is often modeled using the Nernst-Planck equation,
which includes the usual advection and diffusion transport as well as *electromigration*, where charged species are transported by an electric field.

<!--- EchemFEM --->
EchemFEM provides a high-level user interface for a finite element implementation of the Nernst-Planck equation.
The user is simply required to provide physical parameters as well as functions describing the chemical reactions.
Then, the desired transport physics are selected using keyword arguments.
Ionic charge can be modeled using either the Poisson equation or the electroneutrality approximation.
The simulated devices can have resolved electrolyte-electrode interfaces or homogenized porous electrodes, in which case electron conduction is also modeled.
Additionally, finite size effects are available, which includes models such as Generalized Modified Poisson-Nernst-Planck (GMPNP) [@wang2013simulations].
Lastly, a fluid flow solver for the incompressible Navier-Stokes and Navier-Stokes-Brinkman equations is provided.

<!--- Firedrake --->
EchemFEM is based on Firedrake [@FiredrakeUserManual], an open-source finite element package,
enabling straightforward implementation of the governing equations in Python.
Firedake has access to scalable, customizable, solvers through its interface with PETSc [@petsc-user-ref; @petsc-web-page], allowing for parallelization and scalability on computing clusters.
This balance between usability and scalability permits a seamless transition from prototyping to large-scale simulation.
EchemFEM leverages Firedrake's capabilities while further increasing the ease-of-use.
Indeed, since the governing equations are already implemented, little to no knowledge of Firedrake and the finite element method is required to use EchemFEM.

The repository includes several examples of electrochemical devices such as flow reactors, flow batteries, and CO2 electrolyzers.

# Statement of need
<!--- section that clearly illustrates the research purpose of the software and places it in the context of related work --->
<!--- Research Need --->
Electrochemical phenomena are highly complex, making characterization of electrochemical devices through in-operando experiments challenging.
Simulation is an important tool for predicting the performance of electrochemical devices, as well as designing them.
As technologies get scaled up from the laboratory scale to industrial scale, experiments become less tractable and therefore simulation increasingly important.
Naturally, the scalability of simulators is crucial.
Furthermore, many existing models and codes are just one dimensional.
To capture the effects of fluid flow and non-monolithic, architected electrodes, higher-dimensional effects do matter.
For three-dimensional systems, iterative methods and appropriate preconditioners are required to maintain scalability.

<!--- Other codes --->
Currently, commercial software are the most commonly used codes for electrochemistry simulations.
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
In one of the demos, the incompressible Navier-Stokes equations are solved in a reactor with an irregular surface, providing a velocity field for the transport equations.
Similarly, in a flow battery example, the Navier-Stokes-Brinkman equations are solved.

In some cases, for example for fast flows, stabilization schemes that are not offered in other software may be required.
For continuous Galerkin (CG) elements, a streamline-upwind Petrov-Galerkin (SUPG) method for the Nernst-Planck equation is provided.
For discontinuous Galerkin (DG), a custom upwind scheme for the Nernst-Planck equation is used [@roy2023scalable].
In both cases, the upwinding considers the combined advection-migration ``velocity''.

As opposed to commercial software, custom scalable solvers are available in Firedrake.
A plethora of solver options are available through simple PETSc keywords and custom operators for preconditioning can be defined using Firedrake [@Mitusch2019].
In @roy2023scalable, scalable block preconditioners were developed for the electroneutral Nernst-Planck equations with DG and implemented in EchemFEM.

Combining EchemFEM with other Python packages is rather simple.
In @govindarajan2023coupling, multi-scale simulations for CO2 reduction in flow reactors are performed by coupling a microkinetics model from CatMAP [@catmap] with the GMPNP transport model from EchemFEM.
The simulations are two orders of magnitude faster than a previous implementation where the transport is done in COMSOL Multiphysics<sup>&reg;</sup>, due to the tedious interface between the commercial software and CatMAP.

Firedrake's automatic adjoint capabilities facilitate the straightforward solution of PDE-constrained optimization problems [@Mitusch2019], already employed in electrochemical applications [@roy2022topology; @batista2023design].
We are currently investigating optimization problems using EchemFEM.

# Acknowledgements

This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory (LLNL) under Contract DE-AC52-07NA27344, and was partially supported by a Cooperative Research and Development Agreement (CRADA) between LLNL and TotalEnergies American Services, Inc. (affiliate of TotalEnergies SE) under agreement number TC02307 and Laboratory Directed Research and Development (LDRD) funding under projects 19-ERD-035 and 22-SI-006.

# References
