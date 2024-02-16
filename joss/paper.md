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
The shift from fossil fuels towards renewable energy has brought about a rapid increase in the availability of clean electricity. 
However, electricity generated from sources such as wind and solar are limited to intermittent operation due to daily and seasonal variation. 
One solution is to utilize electrochemical devices in energy storage and electrochemical manufacturing applications, where they can harness surplus energy that would otherwise go to waste and decarbonize chemical industries traditionally reliant on petrochemical feedstocks. 
Managing the growing prevalence of renewable energy underscores the importance of developing and scaling up these technologies.
Likewise, the electrification of transport creates an increasing need for efficient electrochemical energy storage devices such as batteries and supercapacitors.
Naturally, simulation tools are required to assist in the design of efficient and industrial-scale electrochemical devices.

<!--- How modeling is used --->
Modeling and simulation are used extensively to describe the physics of the electrochemical and transport mechanisms in electrochemical devices.
These devices have various applications, from miniaturized lithium-ion batteries for medical devices up to industrial-scale hydrogen fuel cells which act as backup power sources. Energy storage devices include batteries and supercapacitors, as well as flow batteries which utilize a flowing electrolyte.
Electrolyzers are devices which use electrical energy to perform electrochemical reactions. Some current industrial applications for electrolysis include the color-alkali process for the production of chlorine gas and the Hall-Héroult process for aluminum production. Active areas of research include the development of electrolyzers that can be used to transform carbon dioxide into useful chemicals and to create hydrogen from water.
In the reverse process, fuel cells use fuels such as hydrogen to generate electricity.

While these devices vary wildly in purpose, the governing equations used to describe them are very similar.
Transport of charged chemical species in a fluid is often modeled using the Nernst-Planck equation,
which includes the usual advection and diffusion transport as well as *electromigration*, where charged species are transported by an electric field.

<!--- EchemFEM --->
EchemFEM provides a high-level user interface for a finite element implementation of the Nernst-Planck equation.
The user is simply required to provide physical parameters as well as functions describing for chemical reactions.
All transport options are then selected using keyword arguments.
Ionic charge can either be modeled using the Poisson equation or the electroneutrality approximation.
The simulated devices can be resolved electrolyte-electrode interfaces, or have homogenized porous electrodes, in which case electron conduction is also modeled.
Finite size effects are also available to allow the use of models such as Generalized Modified Poisson-Nernst-Planck (GMPNP).

<!--- Firedrake --->
EchemFEM is based on Firedrake [@FiredrakeUserManual], an open-source finite element package,
enabling straightforward implementation of the governing equations in Python.
Firedake has access to scalable, customizable, solvers through its interface with PETSc [@petsc-user-ref,petsc-web-page], allowing for parallelization and scalability on computing clusters.
This balance between usability and scalability permits a seamless transition from prototyping to large-scale simulation.
EchemFEM leverages Firedrake's capabilities while further increasing the ease-of-use.
Indeed, since the governing equations are already implemented, little to no knowledge of Firedrake and the finite element method is required to use EchemFEM



Examples

# Statement of need
<!--- section that clearly illustrates the research purpose of the software and places it in the context of related work --->

[@roy2023scalable]

<!--- Commercial software --->

# Acknowledgements

This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory (LLNL) under Contract DE-AC52-07NA27344, and was partially supported by a Cooperative Research and Development Agreement (CRADA) between LLNL and TotalEnergies American Services, Inc. (affiliate of TotalEnergies SE) under agreement number TC02307 and Laboratory Directed Research and Development (LDRD) funding under projects 19-ERD-035 and 22-SI-006.

# References
