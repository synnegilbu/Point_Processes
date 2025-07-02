## Project Description

This repository contains code developed for my master’s thesis on modeling species distributions in morphological trait space using spatial point processes. The focus is on implementing a scalable framework for fitting Log-Gaussian Cox Processes (LGCPs) via the SPDE-INLA approach in R.

The project provides tools for mesh construction in arbitrary dimensions, simulation of latent Gaussian fields and LGCP realizations, and continuous-space inference using the INLA framework. The method is particularly suited for high-dimensional trait spaces where traditional spatial tools are limited.

The accompanying bird trait dataset (from the AVONET database) is used as a motivating example, but the framework is general and supports both real and synthetic data in 2D and higher.

This work addresses the lack of high-dimensional LGCP infrastructure in R by developing and validating a complete end-to-end pipeline for trait-space point process modeling.

For a more detailed explanation of the methods and usage, see the `explanation.md` file included in the repository.
