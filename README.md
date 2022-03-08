# README

<!-- Please add a brief introduction to explain what the project is about    -->

Here I show how to make a regional association plot of some data. We will also
create a MIAMI style plot. And lastly we will check out the co-localization of 
the genetic signal to answer the question whether the signals in two different 
GWAS overlap.

## Where do I start?

You can load this project in RStudio by opening the file called 'racer_coloc_tutorial.Rproj'.

## Project structure

<!--  You can add rows to this table, using "|" to separate columns.         -->
File                       | Description                      | Usage         
-------------------------- | -------------------------------- | --------------
README.md                  | Description of project           | Human editable
LICENSE                    | User permissions                 | Read only     
racer_coloc_tutorial.Rproj | Project file                     | Loads project 
.worcs                     | WORCS metadata YAML              | Read only     
renv.lock                  | Reproducible R environment       | Read only     
report/report.rmd          | Source code for paper            | Human editable
report/references.bib      | BibTex references for manuscript | Human editable
images                     | Generic images directory         | Human editable
scripts                    | Generic scripts directory        | Human editable
scripts/run_me.R           | Script to run analyses           | Human editable
scripts/prepare_data.R     | Script to process raw data       | Human editable
scripts/functions.R        | Some generic functions           | Human editable
scripts/create_worcs.R     | To create WORCS project          | Human editable

<!--  You can consider adding the following to this file:                    -->
<!--  * A citation reference for your project                                -->
<!--  * Contact information for questions/comments                           -->
<!--  * How people can offer to contribute to the project                    -->
<!--  * A contributor code of conduct, https://www.contributor-covenant.org/ -->

# Reproducibility

This project uses the Workflow for Open Reproducible Code in Science (WORCS) to
ensure transparency and reproducibility. The workflow is designed to meet the
principles of Open Science throughout a research project. 

To learn how WORCS helps researchers meet the TOP-guidelines and FAIR principles,
read the preprint at https://osf.io/zcvbs/

## WORCS: Advice for authors

* To get started with `worcs`, see the [setup vignette](https://cjvanlissa.github.io/worcs/articles/setup.html)
* For detailed information about the steps of the WORCS workflow, see the [workflow vignette](https://cjvanlissa.github.io/worcs/articles/workflow.html)

## WORCS: Advice for readers

Please refer to the vignette on [reproducing a WORCS project]() for step by step advice.
<!-- If your project deviates from the steps outlined in the vignette on     -->
<!-- reproducing a WORCS project, please provide your own advice for         -->
<!-- readers here.                                                           -->


# Acknowledgements

Dr. Sander W. van der Laan is funded through grants from the Netherlands CardioVascular Research Initiative of the Netherlands Heart Foundation (CVON 2011/B019 and CVON 2017-20: Generating the best evidence-based pharmaceutical targets for atherosclerosis [GENIUS I&II]). We are thankful for the support of the ERA-CVD program ‘druggable-MI-targets’ (grant number: 01KL1802), the EU H2020 TO_AITION (grant number: 848146), and the Leducq Fondation ‘PlaqOmics’.

The framework was based on the [`WORCS` package](https://osf.io/zcvbs/).

<a href='https://osf.io/zcvbs/'><img src='images/worcs_icon.png' align="center" height="75" /></a> 

#### Changes log
    
    _Version:_      v1.0.0</br>
    _Last update:_  2022-03-08</br>
    _Written by:_   Sander W. van der Laan (s.w.vanderlaan-2[at]umcutrecht.nl).
        
    * v1.0.0 Initial version. 

--------------

#### Creative Commons BY-NC-ND 4.0
Copyright (c) 1979-2022 Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com.

<sup>This is a human-readable summary of (and not a substitute for) the [license](LICENSE).
You are free to share, copy and redistribute the material in any medium or format. 
The licensor cannot revoke these freedoms as long as you follow the license terms.</sup>

<sup>Under the following terms: </sup>

<sup>- Attribution — You must give appropriate credit, provide a link to the license, 
and indicate if changes were made. You may do so in any reasonable manner, but 
not in any way that suggests the licensor endorses you or your use. </sup>
<sup>- NonCommercial — You may not use the material for commercial purposes. </sup>
<sup>- NoDerivatives — If you remix, transform, or build upon the material, you may 
not distribute the modified material. </sup>
<sup>- No additional restrictions — You may not apply legal terms or technological 
measures that legally restrict others from doing anything the license permits.</sup>

<sup>Notices: </sup></br>
<sup>You do not have to comply with the license for elements of the material in the 
public domain or where your use is permitted by an applicable exception or limitation.
No warranties are given. The license may not give you all of the permissions 
necessary for your intended use. For example, other rights such as publicity, 
privacy, or moral rights may limit how you use the material.</sup>