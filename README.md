
## Loss Models Computer Labs, 2024

by Katrien Antonio and Jonas Crevecoeur and Jens Robben

Course materials for the *Loss Models Computer Labs* in December
2024.

📆 December, 2024 <br> ⏲ From 11 am to 1 pm <br> 📌 online

Course materials will be posted in the days before the lectures.

## Overview

<p text-align="justify">

This series of lectures introduces the *essential concepts of building
insurance loss models with R*. You will gain insights in the foundations
of handling **loss data**, including useful *data wrangling and
visualization steps*. You will cover a variety of *discrete and
continuous loss distributions*, and techniques to build more flexible
distributions from standard distributions (by mixing and splicing). You
will learn how to fit these models to actual data and inspect their
goodness-of-fit. Then, you will use the fitted model to estimate risk
measures. You will acquire insights in the foundations of these analytic
methods, learn how to set-up the model building process, and focus on
building a good understanding of the resulting model output and
predictions.

</p>

<p align="justify">

Leaving this workshop, you should have a firm grasp of the working
principles of a variety of loss models for frequency and severity data
and be able to explore their use in practical settings. Moreover, you
should have acquired the fundamental insights to explore some other
methods on your own.

</p>

## Schedule and Course Material

The detailed schedule is subject to small changes.

|     Session      | Duration      | Description | Lecture material                                                                                                                                                                                                                        |
| :--------------: | ------------- | ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|      Day 0       | your own pace | Prework     | sheets prework in [html](https://katrienantonio.github.io/loss-models-2020/sheets/prework_day_0.html) and in [pdf](https://katrienantonio.github.io/loss-models-2020/sheets/prework_day_0.pdf)                                          |
| December 14 & 16 |               |             | sheets in [html](https://katrienantonio.github.io/loss-models-2020/sheets/loss_models_2020_computer_labs.html) and in [pdf](https://katrienantonio.github.io/workshop-loss-reserv-fraud-2020/sheets/loss_models_2020_computer_labs.pdf) |

##### Loss modelling analytics

Topics include:

  - data sets used in the course: MTPL and SecuraRe losses
  - data handling and visualization tools with {ggplot} and {dplyr}
  - building frequency models: Poisson, Negative Binomial, ZI and
    Hurdle, maximum likelihood estimation and goodness-of-fit
  - building severity models: simple to complex parametric
    distributions, splicing to construct a global fit, mixing, estimate
    a risk measure.

Download lecture sheets in
[html](https://katrienantonio.github.io/loss-models-2020/sheets/loss_models_2020_computer_labs.html)
and
[pdf](https://katrienantonio.github.io/loss-models-2020/sheets/loss_models_2020_computer_labs.pdf).

## After the lectures: want to more know?

An excellent collection of tutorials, books, workshops, events is
available via

<http://www.actuarialdatascience.org>

## Prework

<p align="justify">

The workshop requires a basic understanding of R. We prepared some
prework sheets that you can take a look at before the workshop
([html](https://katrienantonio.github.io/workshop-loss-reserv-fraud-2020/sheets/prework_day_0.html)
or
[pdf](https://katrienantonio.github.io/workshop-loss-reserv-fraud-2020/sheets/prework_day_0.pdf)).
Yet another good starting level is the material covered in the
[1-Basic](https://github.com/katrienantonio/workshop-R/tree/master/1%20-%20Basic%20R)
folder of the [Programming in
R](https://github.com/katrienantonio/workshop-R) workshop taught at
Nationale Nederlanden in June 2019.

</p>

Being familiar with statistical or machine learning methods is *not*
required. The workshop gradually builds up these concepts, with an
emphasis on hands-on demonstrations and exercises.

## Software Requirements

You have two options to join the coding exercises covered during the
workshop. Either you join the RStudio cloud workspace dedicated to the
workshop, and then you’ll run R in the cloud, from your browser. Or you
use your local installation of R and RStudio.

### R Studio Cloud

If you prefer not to work with a local installation of R (and the
necessary packages), you can join our workspace on R Studio Cloud. This
should be a very accessible set-up for working with R in the cloud for
the less experienced user\!

Here are the steps you should take (before the workshop):

[no longer valid, please download the material from GitHub, eg using the green button on top of the repo]

  - visit the above link
  - log in by creating an account for RStudio or by using your Google or
    GitHub login credentials
  - join the space
  - at the top of your screen you see ‘Projects’, click ‘Projects’
  - with the ‘copy’ button (on the right) you can make your own version
    of the ‘December 14 & 16’ project; in this copy you can work on the
    exercises, add comments etc.
  - you should now be able to visit the project and see the ‘scripts’
    and ‘data’ folders on the right. Open and run the
    ‘installation-instructions.R’ script from the scripts folder, to
    see if everything works fine.

We will have everything set up for you in the correct way. You only have
to login and open your copy of the project\!

### Local installation

Please bring a laptop with a recent version of R and RStudio installed.
Make sure you can connect your laptop to the internet (or download the
course material one day before the start of the workshop). You will
need:

  - R (at least 3.5.2 <https://cloud.r-project.org/bin/windows/base/>)
  - RStudio
    (<https://www.rstudio.com/products/rstudio/download/#download>)

Please run the following script in your R session to install the
required packages

``` r
packages <- c("tidyverse", "here", "gridExtra", "grid", "rstudioapi", "MASS", "actuar", "statmod", "ReIns", "pscl")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

if(sum(!(packages %in% installed.packages()[, "Package"]))) {
  stop(paste('The following required packages are not installed:\n', 
             paste(packages[which(!(packages %in% installed.packages()[, "Package"]))], collapse = ', ')));
} else {
  message("Everything is set up correctly. You are ready to go.")
}
```

## Instructors

<img src="image/Katrien.jpg" width="110"/>

<p align="justify">

[Katrien Antonio](https://katrienantonio.github.io/) is professor in
insurance data science at KU Leuven and professor in insurance data science at
University of Amsterdam. She teaches courses on data science for
insurance, life and non-life insurance mathematics and loss models.
Research-wise Katrien puts focus on pricing, reserving and fraud
analytics, as well as mortality dynamics.

</p>

<p align="justify">

*Jonas Crevecoeur* is a post-doctoral researcher in biostatistics at KU
Leuven. He recently obtained his PhD within the insurance research group
at KU Leuven and holds the degrees of MSc in Mathematics, MSc in
Insurance Studies and MSc in Financial and Actuarial Engineering (KU
Leuven). Before starting the PhD program he worked as an intern with QBE
Re (Belgium office) where he studied multiline products and copulas.
Jonas was a PhD fellow of the Research Foundation - Flanders (FWO, PhD
fellowship fundamental research).

</p>

Happy learning\!

-----
